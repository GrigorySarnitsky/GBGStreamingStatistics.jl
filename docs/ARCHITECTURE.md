# GBGStreamingStatistics — Architecture and Rationale

This document explains **why** `GBGStreamingStatistics` is built the way it is, in the style of its
sibling `GBGPeriodicFields` (PF). It is **self-contained**: every symbol — a math symbol or a `code`
identifier — is introduced in plain language at its first use, and it reads top to bottom. All
mathematics is rendered LaTeX in standard probability notation. The naming rationale is in §9;
deliberate absences are stated as invariants in §6; open decisions are flagged **[OPEN]**.

---

## 0. Provenance

The driving requirement is **conditional statistics of the highest quality with no storage of the
underlying state and no recompute-from-raw reruns** — which forces the central architectural choice
of §2.

The package is designed to be shared across at least three consumers — **turbulence** (Lagrangian
tracer statistics: conditional averages, transition densities, Kramers–Moyal coefficients), **cloud
microphysics** (droplet growth conditioned on size and supersaturation), and **molecular dynamics**
(the van Hove and intermediate scattering functions *are* two-time conditional densities) — so it is
domain-neutral from the first commit.

---

## 1. Purpose and scope

A **backend-agnostic library for online nonparametric statistics**: it estimates marginal
densities, conditional expectations, conditional and two-time transition densities, and standard
functionals of these, from a **stream of samples that is never stored**. Each sample is observed
once, folded into a fixed-size summary, and discarded.

The relationship to `GBGPeriodicFields` is exact and deliberate:

> **PF is all mathematics, no physics. This package is all statistics, no physics.**

PF provides field operations a solver consumes; it carries no concept of time, no particle, no
tracer (those are the consumer's). This package provides *estimators* a solver consumes; it
carries no concept of velocity, droplet, trajectory, or turbulence — and, one level deeper, **no
probabilistic model**: the data is not assumed stationary, Markov, diffusive, or parametric. Those
are consumer hypotheses the package's diagnostics *test*, never premises its estimators *require*.

**Boundary — samples in, statistics out.** The package never sees a grid, an FFT, a tracer, or a
droplet. It consumes a stream of values (scalars / vectors, optionally lagged pairs tagged by a lag
``\tau``) and returns estimators. The *sampling and transport layer* — interpolating off an Eulerian
field, reading particle states, walking parcel trajectories — stays in each consumer. The instant a
grid or tracer concept crosses this line, reuse across the three domains dies.

**Execution target**, inherited from PF: a **single node**, **shared-memory** only (CPU
multithreading and a single GPU); distributed / multi-GPU is out of scope. Development is
**GPU-first** but the code is **backend-agnostic** (CPU `Array`, CUDA, ROCm) through
KernelAbstractions — §5.

---

## 2. The central principle: the architecture follows established statistics

Every abstraction is the code image of a standard object of probability and mathematical
statistics, and design decisions are settled by **estimation theory**, not implementation
convenience — the direct analog of PF's "architecture follows established mathematics." As in PF,
the recurring lesson is that an *absence* is a precise statement (§6).

### The spine: empirical measure → summary → post-hoc-smoothed estimator

Let ``X_1, X_2, \dots`` be the observed stream. The object that carries all information is the
**empirical measure**
```math
\hat\mu_n \;=\; \frac1n \sum_{i=1}^{n} \delta_{X_i},
```
where ``\delta_x`` is the unit point mass at ``x`` and ``n`` the number of samples seen. We do **not**
store the ``n`` point masses. We store a **finite-dimensional summary** ``T(\hat\mu_n)`` that is

1. **mergeable** — ``T`` is an *additive monoid* (in fact a *cone* — also scalable by `scale!`; §10)
   over the sample stream: the unnormalized summary of
   a batch is a running sum, so two batches' summaries combine by `merge!` and
   accumulate-then-merge equals accumulate-together (to round-off). This is just additivity of
   measures, and it is what makes multi-device reduction and time-chunk stitching free (§5, §6);
2. **rich enough** that every reported quantity is recoverable from it *post hoc*, with no second
   pass over the data.

Every reported quantity is then an **estimator** — a functional
```math
\hat\theta \;=\; S_h\!\big(T(\hat\mu_n)\big)
```
of that summary, where ``S_h`` is a **smoothing** operator (a choice of kernel, bandwidth ``h``, and any
secondary smoothing). The load-bearing decision:

> **The summary ``T`` is accumulated online and is estimand-agnostic; the smoothing ``S_h`` is applied
> entirely post hoc.** Accumulation knows nothing of bandwidth, kernel, or which response channel
> will be asked for. This separation is the architectural invariant of the package — the analog of
> PF's representation/operator split — and it is exactly what delivers the user's requirement:
> *no state storage* (only ``T`` survives a step) and *no reruns* (every ``h``, every channel, every
> secondary smoothing is computed from the one stored ``T``).

### The algebra of measures — a structural lens

The spine has a compact algebraic reading, offered as a **structural lens, not a second motivation**:
the package is and remains **mathematical statistics grounded in practice** — its targets are
consistency, the ``n^{-2/5}`` convergence rates, the analytic oracles, honest bandwidth selection. The
lens only names *why* the spine takes the shape it does (established structure read back onto the code,
as in PF), and settles nothing the statistics does not already settle.

Read the summaries as living in the **cone of (unnormalized) measures**: a commutative monoid under
`merge!`, closed under nonnegative scaling by `scale!` (the cone of §10). The binning map ``T`` is then
**linear**, ``T(\mu + \nu) = T\mu + T\nu`` — a cone homomorphism — and accumulation is the monoid
homomorphism from the multiset of samples (the free commutative monoid; `observe!` folds a batch,
`merge!` unions two summaries). *This is why* `merge!` is elementwise addition and accumulate-then-merge
equals accumulate-together: not a convenience, but linearity of ``T``.

The estimands sit one level up. Every **reported** quantity — the normalized density ``f_X``, the
conditional mean ``m(x) = \mathbb{E}[Y \mid X = x]``, the conditional density ``f_{Y\mid X}`` — is
**homogeneous of degree zero in the measure**: rescale ``\mu \mapsto \lambda\mu`` and it is unchanged, so
it lives not on the cone but on its **projectivization**, the probability measures. The
`scale!`-invariance §10 verifies (a `conditional_mean` unmoved by reweighting) *is* this degree-zero
homogeneity. And because the readout ``S_h`` is degree zero while the cone operations are linear, ``S_h``
is **not** a cone homomorphism — the estimand-agnostic-accumulation / post-hoc-smoothing split restated:
pooling must happen on the cone (linear, mergeable), and the degree-zero readout must stay post hoc. That
is *why* "no reruns" forces the smoothing post hoc.

The lens also fixes the shape of the estimators. Every estimand is a degree-zero readout of one
**linear** object — the kernel-weighted local-polynomial **jet-smoother** of the binned measure,
returning the local Taylor coefficients ``(\beta_0, \beta_1, \dots)`` at the query point — consumed in
one of three ways: a **ratio** of two smoothings (conditioning, the disintegration of the joint —
`conditional_density`), a **link-capped intercept** (the density's local-likelihood ``e^{\beta_0}``), or
a **higher jet component** (the derivative functionals `drift`, `diffusion`, and the response-aware local
bandwidth ``h^\star(x)``). The package's one nonlinearity — the density's log link, and normalization —
lives in the readout, exactly where the cone homomorphism already breaks. The shared primitive that
realizes this in code is `LocalPolynomial.local_fit` (§11).

### Why a fine histogram is the canonical ``T`` (and the honest caveat)

The canonical summary is a **finely binned empirical measure** — a histogram, the *fine-grained
measure* that post-hoc smoothing (a choice of bandwidth ``h``) later coarse-grains into each reported
estimate — optionally carrying **per-bin moments** (below). For nonparametric estimation there is, strictly, *no* finite
sufficient statistic; the histogram is **lossy**. The justification is quantitative and sets the
bin count: with bin width ``b`` and smoothing bandwidth ``h``, the binning error is ``O(b^2)`` while the
irreducible smoothing bias is ``O(h^2)``, so choosing ``b \ll h`` (in practice ``\gtrsim 10^3`` bins over
``\pm 5\sigma``, so ``b \approx \sigma/100 \ll h``) makes the summary **approximately sufficient**: it
loses nothing at the target accuracy (the WARPing / binned-KDE / average-shifted-histogram result,
Hall & Wand 1996; Scott 2015). We therefore say "approximately-sufficient summary," never
"sufficient statistic" — the latter would be the category overclaim PF's correctness policy exists
to catch.

### Estimand families (the operation taxonomy)

Mirroring PF's operation families (vertical algebra / differential / reductions / transfer), the
package's families are organized by **estimand** — the population quantity targeted. Each is
detailed in §4:

- **Marginal density** ``f_X`` — from a 1-D histogram.
- **Conditional expectation** ``m(x) = \mathbb{E}[Y \mid X = x]`` and conditional variance — from
  per-bin power-sums.
- **Conditional / transition density** ``f_{Y\mid X}``, ``p_\tau(y \mid x)`` — from a 2-D (joint)
  histogram, per lag.
- **Conditional increment moments and their functionals** — Kramers–Moyal coefficients,
  drift/diffusion, Pawula / Q diagnostics — domain-neutral functionals of the above; the SDE
  *reading* stays consumer-side.
- **Bandwidth selection** (cross-cutting) — global and the **response-aware local** ``h^\star(x)`` (both built).
- **[OPEN]** a streaming **projection road** (basis projection à la Stochastic Force Inference /
  ULI) as an alternative accumulator.

### The histogram is the summary; the ECF is deferred (likely never needed)

The empirical measure could also be carried in a **Fourier** representation — the **empirical
characteristic function** ``\hat\varphi_n(s) = \frac1n\sum_i e^{\mathrm{i} s X_i}`` — but the histogram
is the **sole** summary. Tail estimation (the heavy-tailed acceleration / ``(\Delta X)^6`` channels) is
robust only in the histogram, where the ECF inverse rings; the ECF would earn its place only for a
strictly near-Gaussian marginal, which the consumers do not require. It is therefore **deferred and,
per the user, most likely never built.**

**This does not forgo the FFT** — and the distinction is the crux. Storing the ECF as a *sufficient
statistic* (an **accumulation** choice) is a different thing from using the FFT to *smooth the
histogram* (a **post-hoc computation** choice). The former is dropped; the latter is kept as a major
performance lever (§4.5): smoothing a binned density or regression is a convolution of the bin counts
with the kernel, which the FFT does in ``O(B\log B)`` for ``B`` bins and — decisively — **amortizes over
bandwidth** (transform the counts once; each bandwidth ``h`` is then a multiply by the kernel's
transform). So FFT-*smoothing* is recommended; the ECF *summary* is not.

### Naming conventions

Following PF's mathematics-first naming:

- **Estimands / estimated quantities** take the **quantity's noun**: `density`, `conditional_mean`,
  `conditional_variance`, `conditional_density`, `transition_density`, `drift`, `diffusion`.
- **Operations on the summary** that are actions take the **verb**: `observe!` (fold a sample
  batch in), `merge!` (combine summaries), `estimate`/`smooth`.
- **In-place** forms append `!` uniformly (`observe!`, `merge!`).
- **Standard named statistical objects** keep their standard name (`histogram`,
  `characteristic_function`, the eponymous `nadaraya_watson`, `local_linear`).

### Documentation policy — **estimand first, estimator second**

This is the user's explicit requirement and the exact mirror of PF's *"physical-space meaning
first — even for spectral-representation operations."* In PF, because a field is a real object and
spectral is only a representation, a docstring must state what the operation does to the real object
in physical space — with the formula — before the spectral form it actually computes. Here, because
**every operation estimates a population quantity from data**:

> A docstring (and every doc-page entry) for an estimator must **state the estimand — the
> population quantity, in standard probability notation, with its defining formula — *before* the
> estimator** that approximates it from the accumulated summary.

`conditional_mean` is the intended reference implementation. It must lead with the estimand
```math
m(x) \;=\; \mathbb{E}[Y \mid X = x] \;=\; \int y\, f_{Y\mid X}(y \mid x)\,\mathrm{d}y
```
— named as the ``L^2`` projection of ``Y`` onto ``\sigma(X)``, the regression function — and *only then*
give the local-linear estimator ``\hat m_h(x)`` computed from the per-bin power-sums. State the
estimand whenever a clean population definition exists; when an object is intrinsically
finite-sample (a bandwidth selector, a `merge!`), say so plainly rather than forcing a population
form — the analog of PF allowing "this is the symbol of ``-\nabla^2``, not a kernel" when no clean
physical form exists.

As in PF, the symbol is **bound to the code on first mention**: `density` ``\equiv f_X``,
`conditional_mean` ``\equiv m(x)``, `transition_density` ``\equiv p_\tau(y\mid x)``, `drift`
``\equiv D^{(1)}``, `diffusion` ``\equiv D^{(2)}``.

**Correctness is its own check** (PF's separate-from-formatting rule). Verify the estimator is
actually **consistent for the stated estimand** under the stated assumptions, and that **every
normalization is exact**: the ``1/n``, the ``1/(nh)``, ``\int K = 1``, the kernel moments
``\mu_2(K) = \int u^2 K(u)\,\mathrm{d}u`` and roughness ``R(K) = \int K^2``, the bias/variance
constants. A stray ``1/h`` or a kernel-moment factor is the easy slip here — the statistics analog of
PF's stray ``1/N`` or Parseval factor.

---

## 3. The accumulated summary — the universal lower half

These types are the load-bearing, domain-agnostic core; they are pure numerics over a stream and
transfer verbatim to every consumer.

- **`Histogram1D`** — a binned empirical measure of a conditioning variable ``X``: linear (WARPing)
  bins over a fixed range (default ``\pm 5\sigma``, ``\gtrsim 10^3`` bins), each bin carrying its
  **count** and **per-bin power-sums about the bin center** ``S_p = \sum_{i\in\text{bin}}(X_i -
  c_j)^p`` for ``p = 1 \dots P`` (default ``P = 4``; raise ``P`` for the heavy-tailed channels — a high
  marginal central moment such as ``M_6`` needs ``P = 6``, and its *sampling variance* reaches the 12th
  moment; the *conditional* increment moments of the Pawula ``Q`` instead need the `Histogram2D`
  order-raise, §11). Because the
  offset ``X_i - c_j`` is small within a bin (``\lesssim b/2 \approx \sigma/200``), the local
  power-sums are **well-conditioned**: they are accumulated by a plain atomic scatter-add — single-pass
  and numerically stable, with **no compensated summation needed**. The global central moments are
  recovered post hoc by the **shift formula** (the binscatter / WARPing local-moment construction;
  Scott 2015; Cattaneo et al. 2024), *exactly* (the power-sums use the actual ``X_i``, not the center)
  and stably for the even high-order moments the consumers need.
- **Bins are exactly uniform — by design.** The binning map is a single multiply-and-floor,
  `i = floor((x - lower) * scale) + 1` with one `scale = nbins / (upper - lower)` — ``O(1)``, branch-free, and
  allocation-free, the GPU atomic-scatter hot path (§5). This *requires* equal-width bins; arbitrary
  edges would force a per-sample search (``O(\log \texttt{nbins})``) and break the kernel. So the constructor takes
  the explicit `(lower, upper, nbins)`, which *states* the uniform contract, rather than a vector or
  range of edges (which would merely *happen* to be uniform — and floating-point rounding need not even
  guarantee that). Non-uniform binning is out of scope.
- **`Histogram2D`** — a joint binned measure of ``(X, Y)``, one per lag ``\tau``, for conditional and
  transition densities; per-cell counts and bivariate power-sums ``S_{rs}`` about the cell center.
- **`HistogramND` [OPEN]** — the ``n``-D generalization, deferred until a consumer needs conditioning
  on more than one variable (the PF rule against speculative structure).
- **per-bin power-sums are the *one* summary scheme** above, in a univariate form (``S_p`` of ``X``,
  for the marginal moments and density) and a bivariate form (``S_{rs} = \sum_{i\in\text{bin}}(X_i -
  c^x)^{r}(Y_i - c^y)^{s}`` about the cell center, the local sufficient statistic for local-polynomial
  regression; Cattaneo et al. 2024). All bandwidths, kernels, and moments are recovered from them
  post hoc.
- **the mergeable-monoid (cone) interface** — `merge!` on every summary is **elementwise addition** (same
  binning ⇒ same bin centers ⇒ counts and power-sums add); the weighted `merge!(into, from; weight)` and
  `scale!` extend it to the cone (§10). Tested: accumulate-then-merge ``=`` accumulate-together, and the
  cone laws.
- **the backend-agnostic scatter-add primitive** *is* `observe!` / `observe_pairs!` — one
  KernelAbstractions + Atomix kernel (one work-item per sample) that atomically scatters the count and
  the power-sums into the preallocated summary, the direct analog of PF's `_bin_by_shell` /
  `anterpolate!`. Allocation-free up to the kernel launch; verified identical on CPU and the AMD GPU
  under `allowscalar(false)` (§5). Storing power-sums *about the bin center* is precisely what keeps
  this a plain atomic add — no privatization, no workspace.

**Numeric precision — accumulate wider than the data.** Samples often arrive in **Float32** (the
GPU-first, low-memory consumers), but the summaries do **not** accumulate in Float32:

- **Counts are integers** (`UInt64` merged). Float32 represents integers exactly only up to
  ``2^{24} = 16{,}777{,}216`` — beyond which adding one sample is silently lost — yet a single ``256^3``
  snapshot is *exactly* ``2^{24}`` points, and a run accumulates ``10^{10}``+ into the peak (modal) bins.
  Float32 counting is disqualified outright.
- **Power-sum accumulators default to Float64.** The local sums are well-conditioned, but the *global* moment
  recovery (the shift formula combines them with large bin-center offsets) and the high-order channels
  still want full precision; Float32 would discard the ~7 digits the Float32 *input* carries, through
  the reduction over ``n``. The power-sums accumulate in Float64 by default, regardless of input type.

Accumulator precision is a **type parameter** (PF-style element-type genericity) but its **default is
Float64**, counts integer — *accumulate-wider-than-the-data*, the deliberate opposite of PF's
preserve-the-field's-eltype (a **reduction** should out-precision its inputs; a **transform** should
not). It costs almost nothing: the histogram is **small** (a few ``\times10^3`` bins ``\times`` a few
Float64 power-sums — kilobytes per channel, versus a gigabyte field), and on the GPU the Float64
arithmetic is hidden behind the memory- and atomic-bound accumulation (§5), not FP64 throughput.

**Two-time samples / lag handling — settled: the package bins pairs; the consumer produces them.**
For a two-time quantity like ``m_\tau(u) = \langle Y_{t+\tau}\mid X_t = u\rangle`` the conditioning
value ``X(t)`` and the response ``Y(t+\tau)`` must be paired **per source** — two separate marginal
histograms cannot recover the joint (for the OU oracle both marginals are the *same* stationary
Gaussian; the entire ``\tau``-dependence lives in the pair). The one histogram entry is made at the
response step, into the bin of ``X(t)``. The package surface for this is just:

- the **memoryless** `observe_pairs!(stat, x_past, y_now; lag)` — bins a batch of pre-formed,
  lag-labelled pairs into the per-lag summary (§4.2 / §4.3); and
- the **binning map** `bin_index(stat, x)` — exposed so a caller can pre-quantize.

Everything that *produces* the pairs is the **consumer's**, being sample-source policy a
domain-neutral package must not hold: the **snapshot buffer** of past ``X(t)`` (held until the response
arrives — the irreducible ``O(M)``-per-lag state for ``M`` sources, optionally stored as compact
`UInt16` bin indices via `bin_index` when memory is tight); which **sources** to sample; which
**lags** to form (e.g. log-spaced); and the **origin cadence / buffer depth** (a memory-vs-variance
budget). Pairing is by source identity — for a stably-ordered set just element-wise array alignment
(no id tag), an explicit id only if the ordering can change (re-meshing, sorting, birth/death).

What the package *keeps* is the per-lag **histograms** — the accumulated power-sums, the small
*output* (``L \times`` bins ``\times`` moments, a few MB), never the ``O(ML)`` snapshot buffer (the
expensive object, consumer-side). For the **nonstationary** case the running origin ``t`` and lag
``\tau`` are explicit labels on the summary set (``p_{t,\tau}``), never assumed away. (If the
ring-and-pair boilerplate proves worth sharing, it can be a **non-domain-aware utility offered
alongside** the summary — never something the package's vocabulary depends on.)

---

## 4. Estimand families (estimators — the exportable upper half)

Each entry is written **estimand first** (per §2), demonstrating the convention.

### 4.1 Marginal density ``f_X``

**Estimand.** ``f_X = \mathrm{d}\mu_X/\mathrm{d}\lambda``, the Lebesgue density:
``P(X \in A) = \int_A f_X\,\mathrm{d}\lambda``.
**Estimator** (binned KDE / WARPing): from bin counts ``n_j`` at bin locations ``\ell_j``,
```math
\hat f_h(x) \;=\; \frac{1}{n\,h}\sum_{j} n_j\, K\!\Big(\frac{x - \ell_j}{h}\Big), \qquad \int K = 1,
```
the binned approximation to ``\frac{1}{nh}\sum_i K(\frac{x-X_i}{h})``. The location ``\ell_j`` is the
recovered within-bin **centroid** ``c_j + S_1^{(j)}/n_j`` (linear binning, the default — matching the
bin's first moment, so the discretization error is ``O(\Delta^2)`` in the bin width ``\Delta``) or the
geometric center ``c_j`` (`binning = :simple`, ``O(\Delta)``); the centroid is read from the already-stored
``S_1``, so the higher-order accuracy is free. Bandwidth from §4.5.

**Kernels.** ``K`` is pluggable (any ``\int K = 1``); two ship — `gaussian_kernel` (default) and
`epanechnikov_kernel` (the MSE-optimal ``\tfrac34(1-u^2)`` on ``|u| \le 1``). `bandwidth` (§4.5) is
Gaussian-tuned, so a non-Gaussian kernel wants a matched ``h`` (Epanechnikov's AMISE-optimal ``h`` is
``\approx 2.21\times`` the Gaussian's).

**Tail-aware density.** Fixed-``h`` KDE is poor in sparse tails; both fixes are `density` knobs. The
**adaptive (variable-bandwidth)** estimator `density(…; adaptive=true)` (Abramson's
``h_j \propto \tilde f(c_j)^{-1/2}``, widening the kernel where data is thin) cuts the tail *variance*;
the **local-log-polynomial** estimator `density(…; degree=d)` (fit ``\log f`` by a binned local Poisson
likelihood — degree 1 exact for exponential tails, degree 2 for Gaussians) cuts the tail *bias*,
capturing the decay shape a flat KDE smears. They compose. On a **bounded support** (a positive quantity such as a droplet radius), `density(…; bounds=(a,b))` renormalizes the degree-0 KDE by the kernel mass inside ``[a,b]`` (Jones 1993), removing the edge bias that otherwise halves the estimate at the boundary; or pass an **asymmetric kernel** (`gamma_kernel` for ``[0,\infty)``, `beta_kernel` for ``[0,1]``; Chen 1999, 2000) whose support matches the data and so eliminates the boundary bias by construction (usually sharper than the renormalization). For a shape with a known parametric skeleton, a **parametric start** `density(…; start=…)` (Hjort & Glad 1995) multiplies a moment-fitted model (`normal_start`, `lognormal_start`, `gamma_start`, `exponential_start`, `inverse_gaussian_start`) by a KDE of its *correction* ``f_X/p`` — near-exact where the model fits (the bias is driven by ``(f_X/p)''``, not ``f_X''``) and a tail-friendly bias reduction otherwise; `normalize=true` rescales any estimate to ``\int\hat f = 1`` exactly (the asymmetric kernels and a start are not self-normalizing). Finally, `density_standard_error` returns the pointwise ``\sqrt{R(K)\,f_X/(nh)}`` variability band. Two consumer caveats remain: the fixed range
``[\texttt{lower}, \texttt{upper})`` drops events beyond it (set it wide, ``\approx \pm 10\text{–}15\sigma``,
to keep the tail), and uniform bins under-resolve power-law tails (bin a transform such as ``\log|X|``
and change variables back; the ``O(1)`` GPU scatter does not log-bin) — though **variable bandwidth is
the preferred fix** except for genuine multi-decade power laws (it stays post-hoc, mergeable, and
estimand-agnostic; non-uniform bins would break all three).

### 4.2 Conditional expectation ``m(x)`` and variance ``\sigma^2(x)``

**Estimand.** ``m(x) = \mathbb{E}[Y\mid X=x]`` (the ``L^2(\sigma(X))`` projection; §2 reference) and
``\sigma^2(x) = \mathrm{Var}(Y\mid X=x) = \mathbb{E}[Y^2\mid X=x] - m(x)^2``.
**Estimator** (local-linear; Fan & Gijbels), with ``K_h(u) = h^{-1}K(u/h)`` the bandwidth-scaled
kernel: from the per-bin power-sums form the kernel-weighted sums ``T_\ell(x) = \sum_i K_h(X_i-x)(X_i-x)^\ell``
and ``P_\ell(x) = \sum_i K_h(X_i-x)(X_i-x)^\ell Y_i``, then
```math
\hat m_h(x) \;=\; \frac{T_2 P_0 - T_1 P_1}{T_0 T_2 - T_1^{2}}.
```
Nadaraya–Watson is the local-constant (degree-0) fit ``\hat m_h = P_0/T_0`` (kept for reference, but
local-linear is the default — it removes the design-density bias term ``2 m'(x) f'(x)/f(x)``).
Conditional variance is the bias-reduced **residual-based** estimator (Fan & Yao 1998, the default
`method = :residual`): a second local fit of the squared residuals
``\widehat{\mathbb E}[(Y-\hat m)^2\mid x]``, which removes the curved-``\mathbb E[Y^2\mid x]`` smoothing
bias that the simple plug-in ``\hat{\mathbb E}[Y^2\mid x] - \hat m^2`` (`method = :plugin`) carries. A
**one-step robust local M-estimator** for heavy-tailed channels (e.g. ``(\Delta X)^6``; Fan & Jiang 2000)
is forthcoming.

### 4.3 Conditional and transition density

**Estimand.** ``f_{Y\mid X}(y\mid x) = f_{X,Y}(x,y)/f_X(x)``; the **two-time transition density**
``p_\tau(y\mid x) = f_{X_{t+\tau}\mid X_t}(y\mid x)`` (nonstationary: ``p_{t,\tau}``).
**Estimator** (double-kernel; Rosenblatt 1969; Hyndman, Bashtannyk & Grunwald 1996): from the joint histogram,
```math
\hat f(y\mid x) \;=\; \frac{\sum_i K_{h_1}(X_i-x)\,K_{h_2}(Y_i-y)}{\sum_i K_{h_1}(X_i-x)}.
```

### 4.4 Conditional increment moments and their functionals

**Estimand.** The conditional increment moment ``M_n(x,\tau) = \mathbb{E}[(X_{t+\tau}-X_t)^n \mid X_t = x]``
— a pure statistical object (a conditional expectation of a power of an increment), estimated by
§4.2 applied to the response ``(\Delta X)^n``.
**Standard functionals** (named, but model-neutral): the **Kramers–Moyal coefficient**
``D^{(n)}(x) = \lim_{\tau\to 0}\frac{1}{n!\,\tau}M_n(x,\tau)``, reported as the finite-``\tau`` ratio
``M_n/(n!\tau)``, with the **drift** ``D^{(1)}`` and **diffusion** ``D^{(2)}``. Forthcoming: the
**all-orders Bell-polynomial correction** (Rydin Gorjão et al. 2021) that removes the finite-``\tau`` bias
by inverting the ``\tau``-expansion rather than extrapolating ``\tau\to 0``, and the **Pawula / Q
diagnostic** ``Q(x) = \mu_4(x,\tau)/(3\,\mu_2(x,\tau)^2)`` — the standardized fourth conditional-increment
moment (``\mu_n`` central, about the conditional mean), ``= 1`` for the Gaussian increments of a diffusion
(Pawula's theorem), the excess ``Q - 1`` flagging jump-like dynamics.

> **Model-neutrality is the invariant here.** These functions are defined and estimated without
> assuming the data *is* a diffusion. That the conditional moments *are* the coefficients of a
> Fokker–Planck operator, that Pawula's theorem makes ``D^{(n\ge 3)}`` vanish, that ``D^{(1)}`` "is"
> drift — these are the **consumer's interpretation**, exactly as "energy" is the consumer's reading
> of PF's `inner(u,u)`. The package supplies the estimators and — forthcoming — the **tests of the
> hypothesis** (the Pawula ratio, a ``\chi^2`` Markov / Chapman–Kolmogorov test); the consumer asserts
> the model.

### 4.5 Bandwidth selection (cross-cutting)

**Estimand.** the MSE-optimal (AMISE) bandwidth for the requested channel — the ``h`` minimising ``\int \mathbb{E}[(\hat f_h - f_X)^2]`` for a kernel ``K`` and local-polynomial degree ``d``.

For marginals, `bandwidth` (§4.1) defaults to the **kernel- and degree-aware normal-reference rule** ``h = C(K,d)\,\hat\sigma\,n^{-1/(2p+1)}``: the bias order is ``p = 2`` at degree 0/1 (rate ``n^{-1/5}``) and ``p = 4`` at degree 2/3 (a higher-order ``O(h^4)`` fit, so the wider rate ``n^{-1/9}``), and the kernel enters only through the canonical factor ``C(K,d) = (8\sqrt{\pi}/3)^{1/5}[R(K)/\mu_2(K)^2]^{1/5}`` at degree 0/1 (the fourth-order equivalent-kernel analog at degree 2/3; Loader 1996) — reducing to Silverman ``(4/3)^{1/5}\hat\sigma\,n^{-1/5}`` for the Gaussian at degree 0 and ``\approx 2.34\,\hat\sigma\,n^{-1/5}`` for the Epanechnikov. All kernel and degree dependence is thus one post-hoc factor (locfit's `Wikk` discipline), with ``R(K)`` and ``\mu_2(K)`` computed by quadrature so any kernel works. The data-driven **Sheather–Jones** plug-in (`rule = :sj`, degree 0/1) adds the solve-the-equation bandwidth — the curvature ``\psi_4 = R(f'')`` and its ``\psi_6`` pilot estimated as Hermite-weighted, **mergeable** reductions of the bin counts (Sheather & Jones 1991), the eval kernel a clean ``\delta(K)`` rescale — so it adapts to the true shape (much tighter than the normal reference on a bimodal sample); its **improved** variant (`rule = :isj`, Botev, Grotowski & Kroese 2010) reads the curvature off the **DCT** of the bin counts instead of a normal-reference pilot, staying reliable on strongly multimodal data. The fit-based **AIC** (`rule = :aic`) adds an information criterion of the binned KDE minimised over an ``h``-grid — a mergeable functional of the counts, the eval kernel a ``\delta(K)`` rescale, with a stable interior optimum (Loader's preferred route is fit-based selection). Its leave-one-out **cross-validation** cousins (likelihood-CV, least-squares LSCV) are a deliberate coverage gap: a binned summary's ties drive their objective to ``h \to 0``. The **asymmetric** kernels take their own scale ``b`` (rate ``n^{-2/5}``); for them `bandwidth` offers the ``\hat\sigma\,n^{-2/5}`` rule of thumb (`rule = :normal`) and the **Chen plug-in** (`rule = :plugin`), the asymmetric analog of Silverman's — ``b^\star = [\mathbb{E}[w(X)^{-1/2}]/(2\sqrt{\pi}\,n\int w^2 (f'')^2)]^{2/5}`` with local-scale weight ``w(x) = x`` (gamma) / ``x(1-x)`` (beta), the curvature ``\int w^2 (f'')^2`` a grid quadrature of the pilot ``\hat f''`` (boundary-valid because ``w \to 0`` there, unlike the ``\psi_4`` U-statistic which assumes ``f`` vanishes at the edge). It lowers the ISE ``\approx 1.4\text{–}1.7\times`` on skewed positive data. The load-bearing one for conditional moments is the **response-aware *local* plug-in** `bandwidth(summary::Histogram2D, x; rule = :local)`:
```math
h^\star(x) \;=\; \left[\frac{R(K)\,\sigma^2(x)}{\mu_2(K)^2\,m''(x)^2\,f_X(x)\,n}\right]^{1/5}.
```
It depends on the **response's** conditional variance ``\sigma^2(x)`` **and** regression curvature
``m''(x)``, so **each response channel (``\Delta X, \dots, (\Delta X)^6``, the acceleration ``A``) gets
its own ``h^\star``** — a single shared bandwidth is provably suboptimal. This is why accumulation
must be estimand-agnostic and rich (per-bin moments to high order): ``\sigma^2(x)``, ``m''(x)``, and
``f_X(x)`` are all recoverable from the *one* stored summary, per channel, post hoc.

**Fast smoothing (FFT). [performance]** Because the data is binned on a regular grid, every smoother
here — the binned KDE (§4.1) and the local-polynomial sums ``T_\ell, P_\ell`` (§4.2) — is a
**convolution** of bin arrays with the kernel. The FFT computes it in ``O(B\log B)`` (Silverman 1982;
Wand 1994; Fan & Marron 1994 for the binned local-poly form), and the decisive gain is
**amortization over bandwidth**: transform the bin counts once, then each ``h`` — and the bandwidth
*grid* the local ``h^\star(x)`` sweep needs — is a multiply by the kernel's transform
(``e^{-h^2\omega^2/2}`` for the Gaussian) and an inverse FFT, no re-convolution. On the GPU the
transforms are tiny (``B\sim10^3``) and **batched** across all channels / lags / ``t``-windows in one
plan. Lead with direct compact-kernel convolution (simple, exact); switch to batched FFT when the
sweep is large; a recursive Gaussian (Deriche) is the ``O(B)``-per-channel alternative for the
fixed-bandwidth Gaussian. **Priority:** this smoothing is post-hoc on a small grid and is dwarfed —
by ``\sim10^3`` in total work — by the per-sample accumulation over ``N^3\times``steps (§5), so the
dominant performance lever is the accumulation hot loop, *not* the smoother. This is the FFT use the
deferred ECF summary (§2) does **not** preclude.

### 4.6 Projection road [OPEN, alternative]

An alternative to bin-and-smooth: estimate drift/diffusion by **projecting onto a function basis**
using running inner-product accumulators (the Gram matrix and projection vector are themselves a
mergeable summary) — Stochastic Force Inference (Frishman & Ronceray 2020) and the underdamped ULI
(Brückner, Ronceray & Broedersz 2020), which fits on ``(X, \dot X)`` data exactly as tracer sampling provides. Shares
the lower-half monoid foundation. Include if a consumer wants parametric-in-basis drift fields;
otherwise defer.

---

## 5. Backend-agnosticism and the hot-loop contract

Inherited wholesale from PF (its invariants 6 and 13):

- **Storage is generic `AbstractArray`**; no scalar indexing of device arrays in package methods;
  allocate via `similar`.
- **Accumulation is allocation-free.** `observe!(summary, batch)` allocates only a **bounded** amount
  (the kernel launch) on a preallocated summary — never an ``O(n)`` or ``O(\text{bins})`` temporary,
  which on a GPU hot loop would escape the host GC and wedge the run. Enforced by an `@allocated` audit
  **on the CPU** (where a materialized array's bytes are visible), the direct analog of PF's
  `test/allocations.jl`.
- **One portable kernel** (KernelAbstractions + Atomix) does the atomic scatter-add binning, on CPU
  and any GPU backend; the vendor stacks (CUDA/AMDGPU) stay test-only dependencies.
- **Mergeable across devices** — `merge!` is the multi-device / multi-chunk reduction (§2), here just
  elementwise addition. The monoid laws are the spec; tested for associativity and accumulate-then-
  merge equality.
- **Contention, not FP64, is the hot-loop enemy.** The conditioning variable is peaked
  (Gaussian-ish), so the modal bins are atomic hot spots. The power-sums-about-the-bin-center choice
  (§3) makes each update a *plain* atomic add — integer-exact for counts, round-off for the Float64
  power-sums — so global atomics handle even heavy contention directly (the proven PF `anterpolate!`
  pattern, verified on the AMD GPU). **Privatization** (per-block local histograms in fast memory,
  merged hierarchically via the monoid) stays available as a future optimization should a channel's
  modal-bin contention dominate, but is not needed for correctness. This per-sample accumulation over
  ``N^3\times``steps is the package's true performance-critical path; the post-hoc smoothing (§4.5) is
  second-order.

Getting backend-portability and mergeability right **from the first commit** is the one piece of
"design for reuse now" that pays — both are expensive to retrofit and are the reason the package
survives contact with three different execution models.

---

## 6. Invariants — deliberate absences are statements

The PF discipline: a missing capability encodes a decision. The candidate invariants:

1. **Accumulation is estimand-agnostic; smoothing is post hoc.** No summary stores a bandwidth,
   kernel, or channel choice. An API that bakes a bandwidth into accumulation is the banned analog
   of PF's spectral product-as-`*`. *(This is the invariant that makes "no reruns" true.)*
2. **Samples in, statistics out — no domain, no model.** Domain-neutral vocabulary
   (sample / variate / response / lag / conditioning value), never particle / tracer / droplet /
   velocity. No stationarity, Markov, diffusion, or parametric assumption is a *premise*; such
   structure is only ever *tested*. (The analog of PF's "no time, no particles.")
3. **Everything is a mergeable monoid — in fact a cone.** Every summary supports `merge!`; the
   unnormalized summary is additive over the stream, and **scalable** (`scale!`), so the summaries form a
   convex cone (equal-weight `merge!` the built special case; the weighted `merge!` and `scale!` the cone
   operations — §10 hyperstatistics / forgetting). No non-mergeable accumulator is admitted.
4. **The hot path is allocation-free and backend-agnostic** (§5).
5. **The summary is rich enough to recover every target post hoc.** A target that would need a
   finer or different summary is a recorded coverage gap, **not** silently approximated — the "no
   reruns" promise made falsifiable.
6. **Accumulate wider than the data.** Counts are integer; moment accumulators default to Float64
   even on Float32 input (§3). A Float32 count saturates (``2^{24}``); a Float32 high-order moment
   loses the data's own precision through the reduction. (The deliberate opposite of PF's
   preserve-the-eltype, for the opposite reason: this is a reduction, not a transform.)

**Deliberate absences (consumer-side, by scope).** Bulk **scalar time-correlation functions**
``C(\tau) = \langle X(t)\,X(t+\tau)\rangle`` are not provided: a `Histogram2D` per lag would store the
whole joint density to extract one cross-moment (``O(T\,n_x n_y)`` state for ``O(T)`` scalars), and a
scalar correlation is cheaply streamed by a dedicated lagged product-sum or FFT. The
`Histogram2D`-per-lag path is for the *conditional / two-time density* objects — ``p_\tau(y \mid x)``,
``D^{(n)}(x)``, the van Hove ``G(r,t)`` — that genuinely need the joint, over a typically *sparse* set
of lags; the consumer owns the scalar-correlation fast path.

---

## 7. Documentation, correctness, and testing conventions

Mirrors PF; the deltas are noted.

- **Estimand-first docstrings** (§2), self-contained pages, **rendered LaTeX in standard notation**,
  every code object **bound to its math symbol** on first mention.
- **A `FORMALIZATION.md` analog.** PF records a typed statement per operation, expansion-checked
  against the code and **test-linked**. Here the typed entry states the **estimand**, the
  **estimator**, and the **consistency claim** — its bias and variance, e.g.
  ``\mathrm{bias}\,\hat m_h(x) \approx \tfrac12 h^2 \mu_2(K)\,m''(x)`` and
  ``\mathrm{Var}\,\hat m_h(x) \approx \frac{R(K)\,\sigma^2(x)}{n\,h\,f_X(x)}`` for local-linear
  regression — every normalization constant pinned, each linked to the test that verifies the
  **convergence rate** numerically.
- **Testing against analytic oracles** — the statistics analog of PF's analytic round-trips:
  estimators are checked against distributions with closed forms (Gaussian marginals; an
  **Ornstein–Uhlenbeck process**, whose drift, diffusion, and transition density
  ``p_\tau(y\mid x)`` are known exactly — the natural end-to-end oracle for §4.3–4.4), plus
  **consistency / rate tests** (error falls at the theoretical ``n^{-2/5}`` as ``n`` grows),
  **mergeability tests**, the **`@allocated` audit**, and **GPU-vs-CPU parity** under
  `allowscalar(false)`.
- **A `check-docs` battery** — adapt PF's `.claude/commands/check-docs.md`: the same matrix of
  (Correct / Self-contained / Notation / Formatted / Current / Cited) × (every surface), with
  *Correct* extended to **consistency and normalization of every estimator against its estimand**.

---

## 8. Proposed repository layout (mirrors PF)

```
GBGStreamingStatistics/
  CLAUDE.md                      # operational summary + invariants (mirror PF's)
  README.md  Project.toml
  src/
    GBGStreamingStatistics.jl    # umbrella
    Accumulators.jl              # foundation: AbstractAccumulator, observe!/merge! monoid interface
    Histograms.jl                # Histogram1D / Histogram2D, WARP binning, power-sums about the bin
                                 #   center, the KA/Atomix scatter kernel, shift-formula recovery
    LocalPolynomial.jl           # the shared local-polynomial jet primitive (local_fit) — §2, §11
    DensityEstimation.jl         # §4.1   f_X
    Regression.jl                # §4.2   m(x), σ²(x)  [estimand-first reference module]
    ConditionalDensity.jl        # §4.3   f_{Y|X}, transition p_τ
    MarkovInference.jl           # §4.4   increment moments, KM/drift/diffusion, Pawula/Q, Markov test
    Bandwidth.jl                 # §4.5   global + response-aware local plug-in
    # CharacteristicFunctions.jl # [OPEN, secondary] ECF summary
    # ProjectionInference.jl     # [OPEN] SFI / ULI road
  docs/   ARCHITECTURE.md  FORMALIZATION.md  TESTING.md  PERFORMANCE.md  src/…  make.jl
  test/   runtests.jl  per-module  oracles.jl(OU/Gaussian)  mergeability.jl  allocations.jl  aqua.jl  jet.jl  gpu/
  benchmark/
  examples/
  .claude/commands/check-docs.md
```

Internal split, made explicit: **`Accumulators` + `Histograms` are the universal lower half**
(every domain shares them verbatim, and they hold the hard, write-once numerics — GPU atomic
high-order moment recursions); **`DensityEstimation` through `Bandwidth` are the estimator upper half**
(turbulence-leaning but genuinely cross-domain — MD's van Hove function is `ConditionalDensity` with
a different sampler). Keeping the halves separable lets a consumer that only wants marginal PDFs
avoid the SDE machinery.

---

## 9. Naming the package

**Settled: `GBGStreamingStatistics`** (confirmed by the user). Operationally precise
(online/streaming is *the* defining constraint), no overclaim, legible to all three consumer
domains, parallel to PF's descriptive compound. The alternatives weighed, for the record:

- **`GBGEmpiricalMeasures`** — the most PF-like choice, naming the *mathematical object* the package
  is about (as "PeriodicFields" names fields), and honest (the empirical measure and its summaries
  *are* the object). Risk: undersells the regression/conditional scope to a casual reader.
- **`GBGSufficientStatistics`** — names the *architectural principle*. Risk: a mild category
  overclaim (there is no finite sufficient statistic in nonparametrics — §2's caveat), which PF's
  correctness policy would itself flag.

`GBGOnlineStatistics` is avoided to prevent confusion with the existing `OnlineStats.jl` (whose
mergeable streaming-moment design is, however, worth mirroring — §10).

---

## 10. Open questions

- **Name — settled:** `GBGStreamingStatistics` (§9).
- **Numeric precision — settled:** integer counts, Float64 moment accumulators by default
  (type-parameterized), Float32 input promoted at accumulation (§3, invariant 6).
- **ECF summary — settled: deferred, likely never built** (§2). FFT-*smoothing* of the histogram is
  kept regardless (§4.5) — a separate, post-hoc computation, not a stored representation.
- **Lag handling — settled:** the package bins pre-formed pairs (`observe_pairs!`) and exposes the
  `bin_index` map; the snapshot buffer and all lag / sampling / cadence policy are consumer-side (§3).
- **[OPEN] Projection road** — include SFI/ULI from the start or defer (§4.6)?
- **[OPEN] `HistogramND`** — wait for a multi-conditioning consumer (§3)?
- **[OPEN] Prior art to mirror** — study `OnlineStats.jl` for its mergeable streaming-moment design
  and `merge` semantics, a comparison still worth making now that `Accumulators` is built (it lacks the
  GPU-resident path and the conditional/SDE estimators, but its composition model is a good CPU-side
  template).
- **[OPEN] Vector-valued conditioning / responses** — reuse PF's `FiberVector`/`FiberTensor` fiber
  algebra, or stay scalar-channel and let the consumer loop components? (Lean scalar-channel: the
  estimators are per-component, and a fiber dependency would re-import structure the package doesn't
  otherwise need.)
- **Weighted merge — settled, built (the cone).** `merge!` is equal-weight addition (the monoid), which
  already makes **ensemble pooling free**: each run / device / time-chunk accumulates its own summary and
  they add (§2). The unnormalised summary is also **scalable** (`scale!`) — multiply every count and
  power-sum by a nonnegative weight and it is the valid summary of a reweighted measure — so it lives in a
  *cone* (closed under nonnegative linear combinations), not merely an additive monoid. That extra freedom
  is exactly **hyperstatistics** (superstatistics; Beck & Cohen 2003): a parameter-weighted superposition
  ``P(x) = \sum_k w_k f_k(x)`` equals ``\sum_k (w_k / N_k)\, T_k`` over per-parameter summaries ``T_k``.
  Built as the weighted `merge!(into, from; weight)`
  (``\texttt{into} \mathrel{+}= \texttt{weight} \cdot \texttt{from}``, the default ``\texttt{weight} = 1``
  the plain monoid add) and `scale!(summary, weight)`, over a **Float-count summary variant**
  (`count = Float64`; a weighted count ``w\,n`` is not an integer); the conditional / transition estimators
  carry the correct local mixing weight ``w_k f_k(x)`` automatically, being ratios of linear functionals of
  the summary. It fits the post-hoc rule: keep one summary per ``\beta`` and choose ``f(\beta)`` after the
  runs, re-weighting without re-accumulating, just as bandwidth is a post-hoc knob. `scale!` also gives
  **exponential forgetting** (`scale!(s, λ)` with ``0 < \lambda < 1`` before the next block) for a
  non-stationary stream. See the §11 decision-log entry.

---

## 11. Decision log

Append-only; mirrors PF's decision log. The settled decisions of §9–§10 are recorded here as the
first entries so later sessions extend rather than re-litigate them (the overlap with the body is
intentional — do not dedupe it).

- **Package name — `GBGStreamingStatistics`** (§9). Operationally precise; avoids `OnlineStats.jl`
  confusion. Alternatives `GBGEmpiricalMeasures` / `GBGSufficientStatistics` weighed and declined.
- **Numeric precision — integer counts, `Float64` moment accumulators by default** (§3, invariant 6),
  type-parameterized; `Float32` input promoted at accumulation.
- **ECF summary — deferred, likely never built** (§2). FFT-*smoothing* of the histogram is kept
  (§4.5) — a post-hoc computation, not a stored representation.
- **Lag handling — the package bins pre-formed pairs** (`observe_pairs!`) and exposes `bin_index`;
  the snapshot buffer and all lag / sampling / cadence policy stay consumer-side (§3).

- **Summary storage — power-sums about the bin center, not Pébay-about-the-running-mean** (§3).
  *Changed during the GPU pass.* The original design accumulated per-bin central moments by the Pébay
  (2008) running-mean recursion with Kahan compensation. A Pébay update is **not atomic**, so the GPU
  hot loop would have needed per-block *privatization* (private histograms merged hierarchically) — a
  pattern the sibling PF uses nowhere. Switching to **power-sums ``S_p = \sum (X_i - c_j)^p`` about the
  fixed bin center** makes `observe!` a plain atomic scatter-add (PF's proven `anterpolate!` pattern):
  allocation-free, one portable kernel identical on CPU and the AMD GPU, `merge!` reduced to
  elementwise addition, and the marginal moments unified with the regression power-sums. The offset
  ``X_i - c_j`` is tiny, so the local sums are well-conditioned (Kahan no longer needed); global central
  moments are recovered post hoc by the shift formula — exact, and stable for the even high-order
  marginal moments the consumers need (odd moments keep weaker relative precision, an accepted trade).
  Verified by the ROCm GPU parity suite.

- **Kramers–Moyal §4.4 — drift and diffusion first, as conditional increment moments** (§4.4). The
  first slice implements the increment moments and the finite-``\tau`` Kramers–Moyal coefficients
  drift ``D^{(1)} = M_1/\tau`` and diffusion ``D^{(2)} = M_2/(2\tau)`` directly from the **order-2
  `Histogram2D`**, via ``M_1(x) = \hat m_h(x) - x`` and
  ``M_2(x) = \widehat{\mathbb{E}}[Y^2\mid x] - 2x\,\hat m_h(x) + x^2`` — reusing the §4.2 local-linear
  primitives (no new accumulator), reported as the plain finite-``\tau`` ratio ``M_n/(n!\,\tau)``. The
  all-orders Bell-polynomial finite-time correction (Rydin Gorjão et al. 2021) and the ``\tau\to0``
  limit stay forthcoming / consumer-side. ``n \ge 3`` — the higher Kramers–Moyal coefficients and
  `pawula_quotient` (whose ``Q = \mu_4/(3\mu_2^2)`` needs ``M_4``) — is a **recorded coverage gap** (it
  needs the joint summary raised to order-``n`` per-cell ``Y`` power-sums, deferred per invariant 5);
  `markov_test` composes the now-available §4.3 transition density.
  OU-oracle-tested: drift recovers the finite-``\tau`` value ``x(e^{-\theta\tau}-1)/\tau`` exactly (no
  smoothing bias — the OU conditional mean is linear — ``\to -\theta x`` as ``\tau \to 0``), diffusion
  ``\to \sigma^2/2`` carrying the local-linear bias of the curved ``\mathbb{E}[Y^2\mid x]``.

- **Conditional / transition density §4.3 — binned double-kernel from the joint counts** (§4.3). The
  Rosenblatt double-kernel estimator (§4.3) is evaluated from the `Histogram2D` **cell counts** (no
  power-sums needed): at the conditioning value ``x`` the joint is reduced to an ``x``-kernel-weighted
  Y-profile ``P_k = \sum_j n_{jk} K_{h_1}(c^x_j - x)`` with ``D = \sum_k P_k``, and then
  ``\hat f(y\mid x) = \sum_k P_k K_{h_2}(c^y_k - y) / (\sqrt{2\pi}\,h_2 D)``. The ``x``-kernel
  normalization cancels in the conditional ratio — only the ``1/(\sqrt{2\pi}\,h_2)`` survives, giving
  ``\int \hat f(y\mid x)\,\mathrm dy = 1`` exactly. `transition_density` is the same computation on a
  lag-``\tau`` joint (the two-time reading). The two bandwidths ``h_1`` (`hx`, conditioning) and
  ``h_2`` (`hy`, response) are supplied; the conditional MISE rate is the slower ``O(n^{-2/3})``.
  OU-oracle-tested against the closed-form Gaussian ``p_\tau(y\mid x)``.

- **Bandwidth §4.5 — the marginal `bandwidth` made kernel- and degree-aware** (§4.5). *Generalised the
  global rule.* The first slice returned only Silverman ``(4/3)^{1/5}\hat\sigma\,n^{-1/5}`` (Gaussian,
  degree 0; rules `:silverman`/`:scott`). A survey of the field (locfit, KernSmooth, ks, statsmodels)
  confirmed the AMISE-optimal ``h`` depends on the kernel only through the canonical factor
  ``\delta(K) = [R(K)/\mu_2(K)^2]^{1/5}`` and on the local-polynomial degree through both the rate
  (``n^{-1/5}`` at degree 0/1, ``n^{-1/9}`` at degree 2) and the equivalent-kernel constant — exactly
  locfit's `Wikk(K, d)`. `bandwidth(summary; rule = :normal, kernel, degree)` now returns
  ``C(K,d)\,\hat\sigma\,n^{-1/(2p+1)}`` (``R(K)``, ``\mu_2(K)`` by quadrature, so any kernel works) and
  is also `density`'s default `h`, so a degree-2 or Epanechnikov call is matched automatically. The
  kernel-as-standard-deviation shortcut (base-R, KDEpy, KernelDensity.jl — kernel-blind, not
  AMISE-exact) was weighed and declined. The data-driven Sheather–Jones plug-in (`rule = :sj`, degree
  0/1) is now added — the ``\psi_4``/``\psi_6`` curvature functionals as **mergeable** Hermite-weighted
  reductions of the bin counts, the solve-the-equation root-find (base-R `bw.SJ` ste), the eval kernel a
  clean ``\delta(K)`` rescale (the pilot stays Gaussian, kernel-independent) — much tighter than the
  normal reference on a bimodal sample. The fit-based AIC follows (next entry).

- **Bandwidth §4.5 — fit-based AIC added; leave-one-out CV is a deliberate absence** (§4.5). The survey
  flagged Loader's *preferred* route as fit-based selection over an ``h``-grid. Implemented as
  `rule = :aic`: minimise the binned-KDE AIC ``-2\sum_j n_j\log\hat f_h(c_j) + 2\,\mathrm{df}(h)`` over
  the grid — a mergeable functional of the counts, the eval kernel a ``\delta(K)`` rescale. But the
  **cross-validation** selectors (likelihood-CV and least-squares LSCV) are **not** built: they are
  leave-one-out, and a *binned* summary has heavy ties (≈10²–10³ samples per bin) whose same-bin term
  drives the CV objective to ``h \to 0`` (verified — LSCV's ``(1/2 + S/n)/(nh\sqrt\pi)`` numerator goes
  negative, so it diverges to ``-\infty`` as ``h \to 0``). base-R `bw.ucv` survives only via coarse
  bins, a constrained search range, and a boundary warning. AIC needs no leave-one-out — both its
  log-likelihood and df-penalty diverge as ``h \to 0`` — so its optimum is a stable interior point. CV
  recorded as a coverage gap (a within-bin-spread correction from the power-sums could rescue it);
  `:sj`, `:isj`, and `:aic` are the data-driven selectors.

- **Bandwidth §4.5 — Improved Sheather–Jones (`rule = :isj`) added** (§4.5). The classic `:sj` plug-in
  estimates the curvature ``\psi_4`` from a *normal-reference* pilot, which can mislead on strongly
  multimodal data (the Marron–Wand claw). Botev, Grotowski & Kroese (2010) instead read the curvature
  norms ``\|f^{(s)}\|^2`` off the **DCT** of the bin counts and solve the AMISE fixed point — a fully
  nonparametric pilot. Added as `rule = :isj`: a mergeable functional of the counts' cosine transform,
  the eval kernel the same ``\delta(K)`` rescale, robust where `:sj`'s pilot misleads.

- **Density §4.1 — parametric start and exact normalization added** (§4.1). A second pass over the
  downloaded reference software (studying source, not docs) found two capabilities worth mirroring: the
  **parametric start** (Hjort & Glad 1995; `kdensity`'s headline feature) and the **exact
  renormalization** (`kdensity`'s `normalized`). The start is the semiparametric
  ``\hat f = p(\cdot;\hat\theta)\cdot\widehat{(f/p)}`` — a moment-fitted model times a KDE of its
  correction — with five moment-fittable starts (`normal_start`, `lognormal_start`, `gamma_start`,
  `exponential_start`, `inverse_gaussian_start`); it stays in-discipline because ``\hat\theta`` is a
  functional of the recovered moments, so the whole estimate is mergeable, and it is a *bias* reduction
  (the leading bias is driven by ``(f/p)''``), valuable for the skewed positive shapes the package targets
  (droplet radii). `normalize=true` divides any estimate by its integral over the support so
  ``\int\hat f = 1`` exactly — the finishing step for the non-self-normalizing asymmetric kernels and
  starts. Both are post-hoc degree-0 options; `start` does not compose with `degree` ``\ge 1``, `adaptive`,
  `bounds`, or the asymmetric kernels (the local-likelihood-with-start generalization is a recorded gap).
  The cross-validated start bandwidth (`KernelEstimator.jl`) and the Hjort–Jones *locally parametric*
  generalization were weighed and deferred. The same source pass also corrected the gamma/beta kernels
  from Chen's naive ``f_1`` to the bias-reduced ``f_2`` (a fix, committed separately).

- **Bandwidth §4.5 — Chen plug-in for the asymmetric kernels; AIC, U-statistic, and gamma-reference
  routes rejected** (§4.5). The asymmetric (gamma / beta) kernels shipped only the ``\hat\sigma\,n^{-2/5}``
  rule of thumb; a data-driven selector was wanted because the optimal scale is ``0.27\text{–}1.6\times`` it
  across densities (``> 2\times`` ISE on log-normal — the droplet regime). Three routes were tried and
  rejected before the plug-in: the **binned AIC** (as for the symmetric kernels) is *monotonic* in ``b``
  — the gamma self-influence ``\sim 1/\sqrt{xb}`` makes the df penalty too weak against the flat binned
  likelihood of dense data, so it oversmooths without bound; the **``\psi_4`` U-statistic** curvature
  (the ``:sj`` machinery) assumes ``f`` vanishes at the boundary (integration by parts) and goes
  *negative* on the half-line where ``f(0) > 0`` (e.g. the exponential); a **gamma-reference** plug-in is
  pathological for shape ``k \le 1`` (a non-integrable ``f''`` boundary singularity, ``b^\star \to 0``).
  The shipped `rule = :plugin` is the Chen plug-in with the curvature ``\int w^2 (f'')^2`` (weight
  ``w = x`` / ``x(1-x)``) estimated by **grid quadrature of the pilot ``\hat f''``** — boundary-valid
  because ``w \to 0`` there, and always real. The ``f''`` pilot is a normal-reference scale
  ``1.3\,\hat\sigma\,n^{-1/9}`` (constant calibrated for robustness, as Silverman's 1.06 is): never
  materially worse than the rule of thumb, ``1.4\text{–}1.7\times`` lower ISE on skewed data. A
  fully-adaptive (solve-the-equation) curvature pilot is a recorded refinement.

- **Hyperstatistics — the cone: `scale!` + weighted `merge!` over a Float-count variant** (§10). *Built.*
  The unnormalised summary is scalable, so the summaries form a convex cone, not merely an additive monoid.
  Realized as `scale!(summary, weight)` (the scalar action — every count and power-sum ``\times`` weight)
  and the weighted `merge!(into, from; weight)`
  (``\texttt{into} \mathrel{+}= \texttt{weight} \cdot \texttt{from}``; ``\texttt{weight} = 1`` recovers the
  plain monoid add), on a **Float-count variant** (`count = Float64`; the default `UInt64` counts cannot
  hold a fractional ``w\,n``, so `scale!` and a weighted merge `ArgumentError` on them). The summary storage
  needed almost no change — `observe!`'s atomic scatter already adds `one(eltype(counts))`, so a Float count routes
  through the atomic-*float* add (CPU and GPU), and `merge!` was already an elementwise broadcast. Estimand:
  the **hyperstatistical** parameter-mixture density
  ``P(x) = \int f(\beta)\,P(x\mid\beta)\,\mathrm{d}\beta \approx \sum_k w_k f_k(x) = \sum_k (w_k/N_k)\,T_k``
  (superstatistics; Beck & Cohen 2003) — keep one summary ``T_k`` per ``\beta_k``, choose the mixing law
  ``f(\beta)`` post hoc, re-weight without re-accumulating; and **exponential forgetting** (`scale!(s, λ)`
  before the next block). Oracle-tested by the χ²-superstatistics → Student-``t`` result (a Gamma-weighted
  Gaussian superposition has the heavy ``t`` tail), a discrete Gaussian mixture, the cone laws, exact
  forgetting, the `Histogram2D` cone, and GPU parity (`test/hyperstatistics.jl`, `FORMALIZATION.md`). One
  caveat recorded: the mixture's total mass is ``\sum_k w_k`` (``\approx 1`` for a normalized law), so its
  *default* bandwidth reads ``n \approx 1`` — pass an explicit ``h`` (or rescale the weights to a
  representative sample size). An effective-sample-size bandwidth for weighted summaries
  (Kish ``(\sum_i w_i)^2 / \sum_i w_i^2``) is a possible refinement, not built.

- **Linear (centroid) binning — the default for `density`** (§4.1). *Built.* The binned KDE places each
  bin's mass at a representative location ``\ell_j``. The geometric center ``\ell_j = c_j`` (simple binning)
  leaves a nonzero within-bin first moment, an ``O(\Delta)`` discretization of the exact KDE; placing the
  mass at the recovered within-bin **centroid** ``\ell_j = c_j + S_1^{(j)}/n_j`` zeroes that first moment
  (``\sum_{i\in j}(X_i - \ell_j) = 0``), so the leading term cancels and the discretization error drops to
  ``O(\Delta^2)`` — the centroid analog of Hall–Wand linear binning. The first power-sum ``S_1`` is *already
  stored*, so the refinement is free, and it is the new default (`binning = :linear`; `binning = :simple`
  keeps the center). One change threads it everywhere: every `density` path reads a `centers` vector, so the
  KDE, the local-log fit, the parametric start, and the asymmetric kernels all inherit the centroid at no
  extra cost. Oracle-tested by the binning-error convergence sub-test (`test/density_convergence.jl`: a
  log–log slope ``\approx 2`` in ``\Delta`` for `:linear` against ``\approx 1`` for `:simple`, the linear
  estimate an order of magnitude tighter at a fine grid). A higher-order correction using the stored ``S_2``
  is conceivable but not built.

- **Estimator architecture — the linear jet-smoother is the primitive, local likelihood a consumer;
  and the measure-algebra lens** (§2, §4). *A framing decision, landed in the code over Stages 1–3.* Two
  settlements taken together. **(1) The lens.** The spine is given an algebraic reading (§2): the
  summaries are the **cone of measures**, ``T`` is the linear cone homomorphism (so `merge!` is
  addition), and every reported estimand is **homogeneous of degree zero in the measure**, descending
  to the projectivization (the probability measures) — exactly the `scale!`-invariance §10 tests, and
  exactly why the post-hoc readout is not a homomorphism. Adopted as a *structural lens only*; the
  motivation stays mathematical statistics grounded in practice (consistency, rates, oracles).
  **(2) The primitive.** Two abstractions were weighed for sharing code between `density` and
  `conditional_mean`: a **local generalized linear model** (local likelihood — `density` at degree
  ``\ge 1`` is a local Poisson/log fit, `conditional_mean` a local Gaussian/identity fit; Loader 1996,
  Fan & Gijbels 1996), versus a **linear local-polynomial jet-smoother** (smooth the binned measure,
  return the local Taylor coefficients ``(\beta_0, \beta_1, \dots)``; conditioning a ratio of two such).
  They are not competitors: a local GLM is fit by **iteratively reweighted least squares**, so *each*
  Newton step is itself a kernel-weighted local-polynomial solve, and the identity-link case
  (regression) converges in **one** step to the plain solve. The GLM is therefore *built from* the
  linear solve, not parallel to it. So the **linear jet-smoother is the irreducible primitive**
  (provisionally `local_fit`), and the GLM is one of three thin consumers: a **link-cap** (the
  density's ``e^{\beta_0}``, the IRLS loop at degree ``\ge 1``), a **ratio** (conditioning —
  `conditional_density`, the disintegration of the joint), and a **jet component** (the derivative
  functionals `drift`, `diffusion`, and ``h^\star(x)``). This keeps the lone nonlinearity on one call
  site and the hot, allocation-free, GPU-portable path linear. *Landed (Stages 1–3):* the shared
  `LocalPolynomial.local_fit` (the jet-smoother, returning ``(\beta_0, \beta_1, \beta_2)``) now backs
  `conditional_mean` (its ``\beta_0``), the local-likelihood `density` (the IRLS inner solve), and the
  increment moments — replacing the formerly duplicated `Regression._local_linear` (degree-1, intercept
  only) and the bespoke Newton step in `DensityEstimation._local_log`. On it: `conditional_mean` /
  `conditional_variance` gained the `degree` (0/1/2) and `kernel` knobs; `conditional_variance` the
  bias-reduced **Fan & Yao** default (`method = :residual`); and the response-aware local bandwidth
  ``h^\star(x)`` = `bandwidth(::Histogram2D; rule = :local)`, which reads the curvature
  ``m''(x) = 2\beta_2`` of a degree-2 fit together with ``\sigma^2(x)`` and ``f_X(x)``. It landed as the
  substrate of that feature, **not** as a standalone refactor (convenience is not a justification, §6). Settled by estimation theory
  (iteratively reweighted least squares is the linear solve iterated; the smoother is linear on the
  cone), not implementation convenience.

Open decisions still flagged **[OPEN]** in §10 (the projection road, `HistogramND`, vector-valued
conditioning, prior-art mirroring) graduate to entries here when they are settled.

---

## References

The load-bearing works the package implements or pins to:

- P. Pébay (2008), [*Formulas for robust, one-pass parallel computation of covariances and
  arbitrary-order statistical moments*](https://doi.org/10.2172/1028931), Sandia SAND2008-6212.
  (Parallel streaming central moments — weighed, then superseded by local power-sums about the bin
  center; §11.)
- D. W. Scott (2015), [*Multivariate Density Estimation*](https://doi.org/10.1002/9781118575574), 2nd
  ed., Wiley; P. Hall & M. P. Wand (1996), [*On the accuracy of binned kernel density
  estimators*](https://doi.org/10.1006/jmva.1996.0009), J. Multivariate Anal. 56. (WARPing /
  binned-KDE sufficiency.)
- M. P. Wand (1994), [*Fast computation of multivariate kernel
  estimators*](https://doi.org/10.1080/10618600.1994.10474656), J. Comput. Graph. Statist. 3(4); J.
  Fan & J. S. Marron (1994), [*Fast implementations of nonparametric curve
  estimators*](https://doi.org/10.1080/10618600.1994.10474629), J. Comput. Graph. Statist. 3(1).
  (FFT-amortized binned smoothing — §4.5.)
- M. D. Cattaneo, R. K. Crump, M. H. Farrell & Y. Feng (2024),
  [*On binscatter*](https://doi.org/10.1257/aer.20221576), Amer. Econ. Rev. 114(5). (Per-bin
  local-polynomial / binscatter.)
- J. Fan & I. Gijbels (1996), [*Local Polynomial Modelling and Its
  Applications*](https://doi.org/10.1201/9780203748725), Chapman & Hall. (Local-linear regression;
  the response-aware local bandwidth.)
- J. Fan & Q. Yao (1998), [*Efficient estimation of conditional variance functions in stochastic
  regression*](https://doi.org/10.1093/biomet/85.3.645), Biometrika 85(3). (Conditional variance.)
- J. Fan & J. Jiang (2000), [*Variable bandwidth and one-step local
  M-estimator*](https://doi.org/10.1007/BF02903849), Sci. China Ser. A 43. (Robust heavy-tailed
  channels.)
- M. Rosenblatt (1969), *Conditional probability density and regression estimators*, in
  [*Multivariate Analysis II*](https://openlibrary.org/works/OL7299345W); R. J. Hyndman, D. M.
  Bashtannyk & G. K. Grunwald (1996), [*Estimating and visualizing conditional
  densities*](https://doi.org/10.1080/10618600.1996.10474715), J. Comput. Graph. Statist. 5(4).
  (Double-kernel conditional density.)
- B. W. Silverman (1986), [*Density Estimation for Statistics and Data
  Analysis*](https://doi.org/10.1201/9781315140919), Chapman & Hall; B. W. Silverman (1982),
  [*Algorithm AS 176: kernel density estimation using the fast Fourier
  transform*](https://doi.org/10.2307/2347084), J. R. Statist. Soc. C 31(1). (The Silverman global
  rule of thumb; FFT-binned KDE.)
- S. P. Sheather & M. C. Jones (1991), [*A reliable data-based bandwidth selection method for kernel
  density estimation*](https://doi.org/10.1111/j.2517-6161.1991.tb01857.x), J. R. Statist. Soc. B
  53(3). (Plug-in bandwidth.)
- Z. I. Botev, J. F. Grotowski & D. P. Kroese (2010), [*Kernel density estimation via
  diffusion*](https://doi.org/10.1214/10-AOS799), Ann. Statist. 38(5). (Improved Sheather–Jones via the
  DCT — `rule = :isj`.)
- I. S. Abramson (1982), [*On bandwidth variation in kernel estimates—a square root
  law*](https://doi.org/10.1214/aos/1176345986), Ann. Statist. 10(4). (Adaptive / variable-bandwidth
  KDE — the sparse-tail estimator.)
- C. R. Loader (1996), [*Local likelihood density estimation*](https://doi.org/10.1214/aos/1032298287),
  Ann. Statist. 24(4). (Local-log-polynomial / local-likelihood density — the tail-shape estimator.)
- M. C. Jones (1993), [*Simple boundary correction for kernel density
  estimation*](https://doi.org/10.1007/BF00147776), Statist. Comput. 3(3). (Renormalization boundary
  correction — `density(…; bounds)`.)
- S. X. Chen (2000), [*Probability density function estimation using gamma
  kernels*](https://doi.org/10.1023/A:1004165218295), Ann. Inst. Statist. Math. 52(3). (Gamma asymmetric
  kernel for ``[0,\infty)`` — `density(…; kernel=gamma_kernel)`.)
- S. X. Chen (1999), [*Beta kernel estimators for density
  functions*](https://doi.org/10.1016/S0167-9473(99)00010-9), Comput. Statist. Data Anal. 31(2). (Beta
  asymmetric kernel for ``[0,1]`` — `density(…; kernel=beta_kernel)`.)
- N. L. Hjort & I. K. Glad (1995), [*Nonparametric density estimation with a parametric
  start*](https://doi.org/10.1214/aos/1176324627), Ann. Statist. 23(3). (Semiparametric density — the
  parametric model times a KDE of its correction; `density(…; start)`.)
- A. Rydin Gorjão et al. (2021), [*Arbitrary-Order Finite-Time Corrections for the Kramers–Moyal
  Operator*](https://doi.org/10.3390/e23050517), Entropy 23(5) — the all-orders Bell-polynomial
  ``\tau``-correction; the [`kramersmoyal`](https://doi.org/10.21105/joss.01693) software lineage.
- A. Frishman & P. Ronceray (2020), [*Learning force fields from stochastic
  trajectories*](https://doi.org/10.1103/PhysRevX.10.021009), Phys. Rev. X 10; D. B. Brückner,
  P. Ronceray & C. P. Broedersz (2020), [*Inferring the dynamics of underdamped stochastic
  systems*](https://doi.org/10.1103/PhysRevLett.125.058103), Phys. Rev. Lett. 125. (Projection road;
  ULI on ``(X,\dot X)``.)
- L. van Hove (1954), [*Correlations in space and time …*](https://doi.org/10.1103/PhysRev.95.249),
  Phys. Rev. 95. (The MD two-time conditional density — the cross-domain anchor.)
- C. Beck & E. G. D. Cohen (2003),
  [*Superstatistics*](https://doi.org/10.1016/S0378-4371%2803%2900019-0), Physica A 322.
  (Hyperstatistics / superstatistics — the parameter-mixture pooling of §10.)

> Every reference above is a published, reader-accessible source — per PF's policy, *never cite an
> inaccessible source*. The private design notes and literature survey behind this architecture are
> deliberately not cited by the package itself.
