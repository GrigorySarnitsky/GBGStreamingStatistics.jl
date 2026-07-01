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

The lens also fixes the shape of the estimators. The universal object is the **summary-as-measure**; every
estimand is a **degree-zero readout** of it, and the readouts fall into **two families, each with one
shared code atom**:

- **Kernel-weighting** — the *local* / smoothing readouts. A query point ``x``, a bandwidth ``h``, and a
  kernel define the per-bin weights ``w_j`` through the single `Kernels._kernel_weight` atom (symmetric
  ``(1/h)K((c_j-x)/h)``, or asymmetric ``K(c_j; x, h)`` — its own normalized density on a bounded support);
  the readouts consume them: a **count sum** (the marginal density), a **Nadaraya–Watson / local-polynomial
  ratio** returning the jet ``(\beta_0, \beta_1, \dots)`` (the `conditional_mean` is its ``\beta_0``; the
  conditional density a ratio of two such smoothings; the response-aware local bandwidth ``h^\star(x)`` a
  higher jet component ``\beta_2``). `local_fit` is the jet readout built on the atom; an asymmetric kernel
  is degree-0 (it has no jet).
- **Power-sum quadrature** — the *global* / moment readouts. A smooth functional ``\mathbb{E}[g(X)]`` is the
  Taylor contraction ``\tfrac1n \sum_j \sum_p g^{(p)}(c_j)/p!\, S_p^{(j)}`` of the stored power-sums — exact
  for polynomials up to the stored order (the `mean` / `var` / `central_moment` recovery), ``O(\Delta^{p+1})``
  for general smooth ``g``. No kernel, no query point.

The package's one nonlinearity — the density's log link, and normalization — lives in the readout, exactly
where the cone homomorphism already breaks. The shared code atoms are `Kernels._kernel_weight` (with
`LocalPolynomial.local_fit` the jet readout on it) and the shift/Taylor contraction of the ``S_p`` (§11);
`local_fit` is the kernel-weighting family's primitive, **not** the universal one — ``\mathbb{E}[g(X)]`` is
a readout with no jet at all.

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
histogram* (a **post-hoc computation** choice). The former is dropped outright; the latter is the
endorsed post-hoc optimization (§4.5), **deferred (recorded, not built)**: smoothing a binned density or
regression is a convolution of the bin counts with the kernel, which the FFT does in ``O(B\log B)`` for
``B`` bins and — decisively — **amortizes over bandwidth** (transform the counts once; each bandwidth
``h`` is then a multiply by the kernel's transform). It is a *pure* optimization — the direct summation
is already exact and dwarfed by the accumulation (§5) — so it waits for a consumer whose ``h^\star(x)``
bandwidth *grid* makes it bite. So FFT-*smoothing* is the endorsed (if deferred) computation; the ECF
*summary* is rejected as a stored representation.

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
For a two-time quantity like ``m_\tau(u) = \mathbb{E}[Y_{t+\tau}\mid X_t = u]`` the conditioning
value ``X(t)`` and the response ``Y(t+\tau)`` must be paired **per source** — two separate marginal
histograms cannot recover the joint (for the OU oracle both marginals are the *same* stationary
Gaussian; the entire ``\tau``-dependence lives in the pair). The one histogram entry is made at the
response step, into the bin of ``X(t)``. The package surface for this is just:

- the **memoryless** `observe_pairs!(stat, x_past, y_now)` — bins a batch of pre-formed pairs into the
  per-lag summary (the lag is the summary's `Histogram2D(…; lag)` label, not an argument here; §4.2 /
  §4.3); and
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

Each estimand is written **estimand first** (per §2), grouped by the **two readout families** of the lens:
**§4.A kernel-weighting** (the local / smoothing readouts through `Kernels._kernel_weight` — the marginal
density §4.1, the conditional mean and variance §4.2, the conditional / transition density §4.3, and
bandwidth selection §4.5) and **§4.B power-sum quadrature** (the global moment readouts off the stored
``S_p`` — `mean` / `var` / `central_moment`, §3, with ``\mathbb{E}[g(X)]`` and the pushforward density
forthcoming). The §4.4 Kramers–Moyal / SDE-reading functionals are **archived** (`archive/`, restore when a
consumer needs them and `Histogram2D` raises its per-cell order); §4.6 is an [OPEN] alternative.

### 4.A Kernel-weighting readouts

The local / smoothing readouts: each a degree-zero readout of the kernel-weighted measure through the
`Kernels._kernel_weight` atom (§2), written estimand-first.

#### 4.1 Marginal density ``f_X``

**Estimand.** ``f_X = \mathrm{d}\mu_X/\mathrm{d}\lambda``, the Lebesgue density:
``P(X \in A) = \int_A f_X\,\mathrm{d}\lambda``.
**Estimator** (binned KDE / WARPing): ``\hat f_h(x) = \frac1{nh}\sum_j n_j K((x-\ell_j)/h)`` (``\int K = 1``)
from counts ``n_j`` at the within-bin **centroid** ``\ell_j = c_j + S_1^{(j)}/n_j`` (`binning = :linear`,
the default — matching the bin's first moment, so the discretization is ``O(\Delta^2)`` in the bin width
``\Delta``; ``\ell_j = c_j`` is `binning = :simple`, ``O(\Delta)``), the centroid free from the stored
``S_1``. Bandwidth from §4.5; kernel pluggable (`gaussian_kernel` default; `epanechnikov_kernel`, the
MSE-optimal ``\tfrac34(1-u^2)``, wants a matched ``h \approx 2.21\times`` the Gaussian's). Tail and shape
options, all post-hoc on the one summary: the **adaptive** variable bandwidth (`adaptive`; Abramson
``h_j \propto \tilde f(c_j)^{-1/2}``, cutting tail *variance*), the **local-log-polynomial** (`degree = d`;
fit ``\log f`` by a binned local-Poisson likelihood — degree 1 exact for exponential tails, 2 for Gaussians
— cutting tail *bias*), a **bounded support** (`bounds`, the Jones-1993 renormalization, or an
**asymmetric** `gamma_kernel` / `beta_kernel`, Chen 1999/2000, boundary-bias-free by construction), a
**parametric start** (`start`; Hjort & Glad 1995, a moment-fitted model times a KDE of its correction
``f_X/p`` — bias-reducing where the model fits), and exact `normalize`; `density_standard_error` is the
pointwise ``\sqrt{R(K)\,f_X/(nh)}`` band. Caveats: the fixed range ``[\texttt{lower}, \texttt{upper})``
drops events beyond it (set it ``\approx \pm 10\text{–}15\sigma`` wide), and uniform bins under-resolve
power-law tails (the variable bandwidth is the preferred fix — it stays post-hoc, mergeable, and
estimand-agnostic, which non-uniform bins would not).

#### 4.2 Conditional expectation ``m(x)`` and variance ``\sigma^2(x)``

**Estimand.** ``m(x) = \mathbb{E}[Y\mid X=x]`` (the ``L^2(\sigma(X))`` projection; §2) and
``\sigma^2(x) = \mathrm{Var}(Y\mid X=x) = \mathbb{E}[Y^2\mid X=x] - m(x)^2``.
**Estimator** (local-linear; Fan & Gijbels): from the joint `Histogram2D`'s kernel-weighted power-sums
``T_\ell = \sum_i K_h(X_i-x)(X_i-x)^\ell`` and ``P_\ell = \sum_i K_h(X_i-x)(X_i-x)^\ell Y_i``,
```math
\hat m_h(x) \;=\; \frac{T_2 P_0 - T_1 P_1}{T_0 T_2 - T_1^{2}}
```
— the ``\beta_0`` of one `LocalPolynomial.local_fit` (§4.A atom); local-linear removes the design-density
bias ``2 m'(x) f'(x)/f(x)`` that the degree-0 Nadaraya–Watson ``P_0/T_0`` carries. The variance is the
bias-reduced **residual-based** Fan & Yao (1998) estimator (`method = :residual`, the default — a second
local fit of the squared residuals ``\widehat{\mathbb E}[(Y-\hat m)^2\mid x]``, removing the
curved-``\mathbb E[Y^2\mid x]`` smoothing bias the `:plugin` ``\hat{\mathbb E}[Y^2\mid x] - \hat m^2``
carries). A one-step robust local M-estimator for heavy-tailed channels (Fan & Jiang 2000) is forthcoming.

**The general reader** `conditional_expectation` ``\equiv \mathbb{E}[g(Y)\mid x]`` for a smooth ``g``
subsumes both: its response is the per-``X``-column **``y``-quadrature** of ``g`` (the order-2
``Y``-power-sums contracted against ``g``'s local coefficients — the §4.B atom), regressed in ``x`` by the
same `local_fit`. So it **composes the two atoms** — kernel-weighting in ``x`` after power-sum quadrature
in ``y`` — with ``m(x) = \mathbb{E}[Y\mid x]`` its ``g = \mathrm{identity}`` case (`conditional_mean` folded
in exactly) and the `:plugin` variance its ``g = y,\,y^2`` combination. It is exact for ``g`` of degree
``\le 2`` in ``y`` and ``O(\Delta_y^2)`` for a general ``g`` (the fixed order-2 joint summary). The
**conditional pushforward** `conditional_pushforward_density` ``\equiv f_{g(Y)\mid X}(z\mid x)`` — the
density of the transformed response — is then one `conditional_expectation` of the composed kernel,
``\hat{\mathbb E}[K_{h_z}(g(Y)-z)\mid x]`` (no ``g^{-1}`` / Jacobian; a symmetric ``K_z`` integrates to 1);
at ``g = \mathrm{identity}``, `degree = 0` it coincides with the double-kernel §4.3 `conditional_density`.

#### 4.3 Conditional and transition density

**Estimand.** ``f_{Y\mid X}(y\mid x) = f_{X,Y}(x,y)/f_X(x)``, and the two-time transition density
``p_\tau(y\mid x) = f_{X_{t+\tau}\mid X_t}(y\mid x)`` (nonstationary: ``p_{t,\tau}``).
**Estimator** (double-kernel; Rosenblatt 1969; Hyndman, Bashtannyk & Grunwald 1996):
```math
\hat f(y\mid x) \;=\; \frac{\sum_i K_{h_1}(X_i-x)\,K_{h_2}(Y_i-y)}{\sum_i K_{h_1}(X_i-x)},
```
with a separate conditioning kernel `kernelx` (in ``x``) and response kernel `kernely` (in ``y``). For a
**bounded / positive response** — a droplet radius as ``Y`` — `kernely` may be **asymmetric** (gamma / beta;
Chen 1999, 2000), boundary-bias-free in ``y``, with `normalize` restoring ``\int\hat f\,\mathrm dy = 1`` —
the conditional analog of the marginal density's asymmetric kernels. The conditioning kernel cancels in the
``f_{X,Y}/f_X`` ratio, so it stays symmetric (a bounded *conditioning* variable needs none — the
local-linear `conditional_mean` is boundary-adaptive, §4.2).

#### 4.5 Bandwidth selection (cross-cutting)

**Estimand.** the AMISE-optimal bandwidth — the ``h`` minimising ``\int \mathbb{E}[(\hat f_h - f_X)^2]``
for a kernel ``K`` and local-polynomial degree ``d``.
**Estimator.** for a marginal, the **kernel- and degree-aware normal-reference rule**
``h = C(K,d)\,\hat\sigma\,n^{-1/(2p+1)}`` (bias order ``p = 2`` at degree 0/1, ``p = 4`` at degree 2/3; the
kernel enters through one canonical factor ``[R(K)/\mu_2(K)^2]^{1/5}`` by quadrature — Silverman
``(4/3)^{1/5}`` the Gaussian/degree-0 case, Epanechnikov ``\approx 2.34``), plus the data-driven
**Sheather–Jones** `:sj` and DCT-based `:isj` (Botev, Grotowski & Kroese 2010) plug-ins, the fit-based
`:aic`, and the asymmetric-kernel Chen `:plugin` (scale ``b \sim n^{-2/5}``) — each a mergeable functional
of the counts (the eval kernel a clean ``\delta(K)`` rescale); leave-one-out CV is a recorded gap (binned
ties drive ``h \to 0``). The load-bearing one for conditional moments is the **response-aware *local*
plug-in** `bandwidth(::Histogram2D, x; rule = :local)`,
```math
h^\star(x) \;=\; \left[\frac{R(K)\,\sigma^2(x)}{\mu_2(K)^2\,m''(x)^2\,f_X(x)\,n}\right]^{1/5},
```
which depends on the response's variance ``\sigma^2(x)`` **and** curvature ``m''(x)`` — so **each response
channel gets its own ``h^\star``** (a single shared bandwidth is provably suboptimal), all three of
``\sigma^2(x)``, ``m''(x)``, ``f_X(x)`` recovered from the *one* summary, post hoc. *(Performance: every
smoother here is a convolution of bin arrays with the kernel, so an FFT computes it in ``O(B\log B)``
amortized over the bandwidth grid the ``h^\star(x)`` sweep needs — Silverman 1982; Wand 1994 — a deferred
optimization (recorded, not built), and in any case dwarfed by the per-sample accumulation hot loop, §5.)*

### 4.B Power-sum quadrature readouts

The global (non-local) readouts: a smooth functional of ``X`` read straight off the stored power-sums — no
kernel, no query point. Within a bin, ``S_p^{(j)} = \sum_{i\in\text{bin}}(X_i - c_j)^p`` is the within-bin
``p``-th moment, so Taylor-expanding any smooth ``g`` about the center,
```math
\mathbb{E}[g(X)] \;=\; \frac1n \sum_j \sum_{p=0}^{\texttt{order}} \frac{g^{(p)}(c_j)}{p!}\, S_p^{(j)}
\;+\; O\!\big(\Delta^{\,\texttt{order}+1}\big),
```
**exact** for polynomials of degree ``\le \texttt{order}`` (the remainder terminates) — this is the
recovered `mean` / `var` / `central_moment` (§3, the shift formula, the ``g = x^p`` case) — and
``O(\Delta^{\,\texttt{order}+1})`` for general smooth ``g``. **Built** in `Expectations`:

- `expectation` ``\equiv \mathbb{E}[g(X)]`` for a general smooth `g`. The per-bin coefficients
  ``g^{(p)}(c_j)/p!`` are read off a **degree-``\texttt{order}`` polynomial interpolant of ``g`` on each
  bin** — no automatic differentiation and no new dependency (uniform bins ⇒ one Chebyshev–Vandermonde
  inverse serves every bin), which reduces *exactly* to the shift formula on the monomials. Both
  `expectation` and the recovered moments route through the one shared contraction atom
  `Expectations._quadrature` (the moments supplying closed-form coefficients, `expectation` interpolation
  coefficients) — the literal realization of the two-family lens.
- `expectation_standard_error` ``= \sqrt{(\mathbb{E}[g^2] - \mathbb{E}[g]^2)/n}`` — the sampling error of
  ``\hat{\mathbb{E}}[g]``, itself a second quadrature readout (``\mathbb{E}[g^2]`` is `expectation` with
  the smooth integrand ``g^2``).
- `pushforward_density` ``\equiv f_{g(X)}``, the density of ``Y = g(X)``. Rather than form ``g^{-1}`` or a
  Jacobian, it is written as the **expectation of the composed kernel test function**,
  ``\hat f_{g(X)}(y) = \tfrac1n\sum_i K_h(g(X_i) - y) = \hat{\mathbb{E}}[K_h(g(X) - y)]`` — so it is one
  more `expectation` call. A **non-monotone** ``g`` is handled automatically (every preimage contributes),
  a flat spot ``g'=0`` shows up as the expected (integrable) caustic, and a symmetric ``K`` makes
  ``\int \hat f\,\mathrm dy = 1`` exactly.

Limits: ``g`` smooth on the bin scale, mass inside the range, and — for the pushforward — the response
bandwidth ``h \gtrsim |g'|\,\Delta`` so ``K_h\circ g`` is resolved across a bin.

### 4.6 Projection road [OPEN, alternative]

An alternative to bin-and-smooth for the conditional-functional / SDE-inference channels (the archived §4.4
road, or general regression fields): estimate by **projecting onto a function basis** with running
inner-product accumulators (the Gram matrix and projection vector themselves a mergeable summary) —
Stochastic Force Inference (Frishman & Ronceray 2020) and the underdamped ULI (Brückner, Ronceray &
Broedersz 2020), fitting on ``(X, \dot X)`` data. Shares the lower-half monoid foundation; deferred until a
consumer wants parametric-in-basis fields.

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
                                 #   center, the KA/Atomix scatter kernel (the pure summary)
    Kernels.jl                   # all the smoothing kernels (symmetric + asymmetric) + _quad quadrature
    Expectations.jl              # §4.B   E[g(X)], moments (mean/var/…), pushforward — power-sum quadrature
    LocalPolynomial.jl           # the shared local-polynomial jet primitive (local_fit) — §2, §11
    MarginalDensity.jl           # §4.1   f_X
    ConditionalExpectations.jl   # §4.2   m(x), σ²(x)  [estimand-first reference module]
    ConditionalDensity.jl        # §4.3   f_{Y|X}, transition p_τ
    Bandwidth.jl                 # §4.5   global + response-aware local plug-in
    # IncrementMoments.jl        # [archived → archive/] §4.4 increment moments, KM/drift/diffusion
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
high-order moment recursions); **`MarginalDensity` through `Bandwidth` are the estimator upper half**
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
- **ECF summary — settled: deferred, likely never built** (§2). FFT-*smoothing* of the histogram is a
  separate, **deferred** post-hoc optimization (§4.5) — a computation choice, not a stored representation.
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
  runs, re-weighting without re-accumulating, just as bandwidth is a post-hoc choice. `scale!` also gives
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
- **ECF summary — deferred, likely never built** (§2). FFT-*smoothing* of the histogram is a separate,
  **deferred** post-hoc optimization (§4.5) — a computation choice, not a stored representation.
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
  ``\hat f(y\mid x) = \sum_k P_k K_{h_2}(c^y_k - y) / (h_2 D)`` for a symmetric `kernel` ``K`` (any
  ``\int K = 1``, default Gaussian). The ``x``-kernel weighting cancels in the conditional ratio — only
  the ``1/h_2`` survives, and ``\int K = 1`` gives ``\int \hat f(y\mid x)\,\mathrm dy = 1`` exactly.
  `transition_density` is the same computation on a lag-``\tau`` joint (the two-time reading). The two
  bandwidths ``h_1`` (`hx`, conditioning) and ``h_2`` (`hy`, response) and the `kernel` are supplied; the
  conditional MISE rate is the slower ``O(n^{-2/3})``.
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
  (`local_fit`), and the GLM is one of three thin consumers: a **link-cap** (the
  density's ``e^{\beta_0}``, the IRLS loop at degree ``\ge 1``), a **ratio** (conditioning —
  `conditional_density`, the disintegration of the joint), and a **jet component** (the derivative
  functionals `drift`, `diffusion`, and ``h^\star(x)``). This keeps the lone nonlinearity on one call
  site and the hot, allocation-free, GPU-portable path linear. *Landed (Stages 1–3):* the shared
  `LocalPolynomial.local_fit` (the jet-smoother, returning ``(\beta_0, \beta_1, \beta_2)``) now backs
  `conditional_mean` (its ``\beta_0``), the local-likelihood `density` (the IRLS inner solve), and the
  increment moments — replacing the formerly duplicated `ConditionalExpectations._local_linear` (degree-1, intercept
  only) and the bespoke Newton step in `MarginalDensity._local_log`. On it: `conditional_mean` /
  `conditional_variance` gained the `degree` (0/1/2) and `kernel` options; `conditional_variance` the
  bias-reduced **Fan & Yao** default (`method = :residual`); and the response-aware local bandwidth
  ``h^\star(x)`` = `bandwidth(::Histogram2D, x; rule = :local)`, which reads the curvature
  ``m''(x) = 2\beta_2`` of a degree-2 fit together with ``\sigma^2(x)`` and ``f_X(x)``. It landed as the
  substrate of that feature, **not** as a standalone refactor (convenience is not a justification, §6). Settled by estimation theory
  (iteratively reweighted least squares is the linear solve iterated; the smoother is linear on the
  cone), not implementation convenience.

- **Kernels — the asymmetric kernels consolidated into the foundation module** (§4.1, §4.3). The
  varying-shape boundary kernels (`gamma_kernel`, `beta_kernel`, the `AsymmetricKernel` type; Chen 1999,
  2000) moved from `MarginalDensity` to `Kernels`, joining the symmetric `gaussian_kernel` /
  `epanechnikov_kernel`. They had stayed in `MarginalDensity` only because, when `Kernels` was first
  extracted, `density` was their sole consumer — a historical artifact, not a design boundary. **A kernel
  is a kernel**: an asymmetric boundary kernel is a smoothing kernel for bounded / positive support, so it
  belongs in the shared foundation, where any density estimator — marginal or conditional — can reach it.
  This gives the `kernel` argument one documented home (no cross-module split) and enabled the now-built
  asymmetric-kernel **conditional density** (§4.3: a bounded *response* — e.g. a droplet radius as ``Y`` —
  takes a boundary-free ``y``-kernel via `conditional_density(…; kernely = gamma_kernel, normalize = true)`,
  the same first-order boundary-bias fix the marginal density already has). It does **not** motivate asymmetric *conditioning* in `conditional_mean`: local-linear regression
  is boundary-adaptive (Fan 1992), already correcting the ``O(h)`` boundary bias to ``O(h^2)`` with a
  symmetric kernel, so an asymmetric conditioning kernel would be a boundary-*variance* refinement only,
  not a correctness fix — a recorded possibility, not built (convenience / parity is not a justification, §6).

- **Kernels / LocalPolynomial — the kernel-weighting atom `_kernel_weight`; the lens corrected to two
  readout families** (§2, §4). The symmetric/asymmetric kernel protocols, previously re-coded inline at four
  sites (`local_fit`, the density KDE, its asymmetric branch, the conditional density), are unified in one
  `Kernels._kernel_weight` (the ``1/h`` folded in, so every readout — count sum, NW ratio, jet — is uniform;
  the fold is ``\beta``-invariant for the local-polynomial solve, so symmetric results are bit-unchanged).
  `local_fit` carries the **one** guard: an asymmetric kernel has no jet, so it is admitted at ``\texttt{degree} = 0``
  only (a clean `ArgumentError`, not the former `MethodError`). This made degree-0 `conditional_mean` accept
  asymmetric kernels uniformly — **building** the degree-0 conditioning-side asymmetric kernel the entry
  above had left "not built". **Honest finding:** the Chen gamma/beta kernels are a *density* boundary tool;
  their asymmetry interacts with the regression design density, so asymmetric-kernel Nadaraya–Watson does
  **not** generally improve the conditional mean (degree-1 local-linear — boundary-adaptive, unbiased for a
  linear mean — is the regression tool). The capability is admitted for **uniformity of the abstraction**,
  not as the recommended regression estimator; the code is uniform, the guidance separate (in
  `conditional_mean`). The episode also **corrected the lens**: `local_fit` is the kernel-weighting family's
  primitive, not *the* universal one — ``\mathbb{E}[g(X)]`` is a power-sum-quadrature readout with no jet,
  the second family. Settled by the measure-algebra lens, not convenience.
- **IncrementMoments §4.4 — archived** (§4). The SDE-reading drift / diffusion / Kramers–Moyal functionals
  (working, OU-tested; Pawula / Markov were stubs) moved to `archive/` to focus the package on its universal
  statistical core — they are thin functionals of `conditional_mean`, the most physics-flavored,
  least-complete, most consumer-side part, and the identity is sharper without them. Restore when a consumer
  needs them and `Histogram2D` raises its per-cell order for ``n \ge 3``.
- **FFT-amortized smoothing — deferred (recorded, not built)** (§4.5). Every reported quantity is already
  computed *exactly* by direct summation over the few ``\times 10^3`` bins, and that cost is dwarfed by the
  per-sample accumulation (§5), so the FFT is a *pure* amortization win — decisive only once a consumer
  sweeps a *grid* of bandwidths (the local ``h^\star(x)`` road). Recorded here and in §4.5, to be built
  when that cost bites. `FFTW` is retained — the ISJ bandwidth (`rule = :isj`) uses its DCT today; the
  unused `AbstractFFTs` direct dependency (its sole intended consumer was this deferred smoothing) was
  dropped, per the rule against speculative structure (§6).
- **Power-sum quadrature readouts §4.B — built** (`expectation`, `expectation_standard_error`,
  `pushforward_density`). The second readout family lands: a smooth functional ``\mathbb{E}[g(X)]`` of the
  binned measure, its sampling standard error, and the pushforward density of ``g(X)``. All three route
  through the one shared contraction atom `Expectations._quadrature`, into which `mean` / `var` /
  `central_moment` were **folded** as the ``g = x^p`` case (guarded by the unchanged exactness tests) — the
  literal realization of the two-family lens. **Decision — the boundary between the two families is that
  they consume *different sufficient statistics*:** kernel-weighting (§4.A) reads counts + bin locations,
  weights by ``K_h``, and localizes/solves (a function of a query point, needs a bandwidth); power-sum
  quadrature (§4.B) contracts the within-bin power-sums ``S_p`` against a test function's local coefficients
  (a scalar, no bandwidth). The summary stores **both** counts and power-sums precisely because there are
  two readout modes. This settles the tempting over-unification: the kernel estimators are *not* trivial
  cases of ``\mathbb{E}[g(X)]`` (they are ratios / local solves of *joint* moments; only the bare degree-0
  fixed-``h`` KDE is formally ``\mathbb{E}[K_h(\cdot - x_0)]``), so a single universal map would collapse to
  the vacuous "every estimator is a functional of the empirical measure" — the mirror of the corrected
  "the jet is *the* primitive" lens. Two atoms, each estimator a readout of one. `g`'s coefficients come
  from a dependency-free local polynomial interpolant; the pushforward is `expectation` of ``K_h\circ g``.

- **ConditionalExpectations §4.2 — general ``\mathbb{E}[g(Y)\mid x]`` + conditional pushforward, built.**
  The general `conditional_expectation` (with `conditional_mean` folded in as ``g = \mathrm{identity}``, and
  `conditional_pushforward_density` ``= f_{g(Y)\mid X}``) makes the conditional side symmetric with the
  marginal `Expectations`. **Refines the two-family boundary:** conditional functionals *compose* the atoms
  — a per-``X``-column power-sum quadrature in ``y`` (§4.B) fed as the response to a kernel-weighted
  regression in ``x`` (§4.A). So "each estimator is a readout of *one* atom" (the §4.B entry above) softens
  to "readouts are built from the atoms, sometimes composing them" — `conditional_mean` always did (its
  trivial ``\sum Y`` response is the ``y``-quadrature of the identity); the general ``g`` just makes the
  inner ``y``-readout non-trivial. The order-2 ``Y``-power-sums cap a general ``g`` at ``O(\Delta_y^2)``.

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
