!!! info "Written by AI, not regularly reviewed"
    This page is written by AI and is not reviewed regularly by a human.

# GBGStreamingStatistics — Architecture and Rationale

This document explains why `GBGStreamingStatistics` is built the way it is, in the style of its
sibling [`GBGPeriodicFields`](https://grigorysarnitsky.github.io/GBGPeriodicFields.jl/) (PF). It is self-contained: every symbol — a math symbol or a `code`
identifier — is introduced in plain language at its first use, and it reads top to bottom. All
mathematics is rendered LaTeX in standard probability notation. The naming rationale is in §9;
deliberate absences are stated as invariants in §6; open decisions are flagged [OPEN].

---

## 0. Provenance

The driving requirement is conditional statistics of the highest quality with no storage of the
underlying state and no recompute-from-raw reruns — which forces the central architectural choice
of §2.

The package is designed to be shared across at least three consumers — turbulence (Lagrangian
tracer statistics: conditional averages, transition densities, Kramers–Moyal coefficients), cloud
microphysics (droplet growth conditioned on size and supersaturation), and molecular dynamics
(the van Hove and intermediate scattering functions *are* two-time conditional densities) — so it is
domain-neutral from the first commit.

---

## 1. Purpose and scope

A backend-agnostic library for online nonparametric statistics: it estimates marginal
densities, conditional expectations, conditional and two-time transition densities, and standard
functionals of these, from a stream of samples that is never stored. Each sample is observed
once, folded into a fixed-size summary, and discarded.

The relationship to `GBGPeriodicFields` is exact and deliberate:

> **PF is all mathematics, no physics. This package is all statistics, no physics.**

PF provides field operations a solver consumes; it carries no concept of time, no particle, no
tracer (those are the consumer's). This package provides *estimators* a solver consumes; it
carries no concept of velocity, droplet, trajectory, or turbulence — and, one level deeper, no
probabilistic model: the data is not assumed stationary, Markov, diffusive, or parametric. Those
are consumer hypotheses the package's diagnostics *test*, never premises its estimators *require*.

**Boundary — samples in, statistics out.** The package never sees a grid, an FFT, a tracer, or a
droplet. It consumes a stream of values (scalars / vectors, optionally lagged pairs tagged by a lag
``\tau``) and returns estimators. The *sampling and transport layer* — interpolating off an Eulerian
field, reading particle states, walking parcel trajectories — stays in each consumer. The instant a
grid or tracer concept crosses this line, reuse across the three domains dies.

**Execution target**, inherited from PF: a single node, shared-memory only (CPU
multithreading and a single GPU); distributed / multi-GPU is out of scope. Development is
GPU-first but the code is backend-agnostic (CPU `Array`, CUDA, ROCm) through
KernelAbstractions — §5.

---

## 2. The central principle: the architecture follows established statistics

Every abstraction is the code image of a standard object of probability and mathematical
statistics, and design decisions are settled by estimation theory, not implementation
convenience — the direct analog of PF's "architecture follows established mathematics." As in PF,
the recurring lesson is that an *absence* is a precise statement (§6).

### The spine: empirical measure → summary → post-hoc-smoothed estimator

Let ``X_1, X_2, \dots`` be the observed stream. The object that carries all information is the
*empirical measure* ``\mathbb{P}_n``, which gives each set ``A`` the share of the samples that fall in it,
```math
\mathbb{P}_n(A) = \frac{\#\{i : X_i \in A\}}{n},
```
where ``n`` is the number of samples seen. We do not
store the ``n`` samples. We store a finite-dimensional summary ``T(n\mathbb{P}_n)`` — a functional of
the unnormalized empirical measure ``n\mathbb{P}_n``, which counts the samples in each set (the additive object; the
``1/n`` normalization is applied post hoc) — that is

1. **mergeable** — ``T`` is an *additive monoid* (in fact a *cone* — also scalable by `scale!`; §10)
   over the sample stream: the unnormalized summary of
   a batch is a running sum, so two batches' summaries combine by `merge!` and
   accumulate-then-merge equals accumulate-together (to round-off). This is just additivity of
   measures, and it is what makes multi-device reduction and time-chunk stitching free (§5, §6);
2. **rich enough** that every reported quantity is recoverable from it *post hoc*, with no second
   pass over the data.

Every reported quantity is then an estimator, whose estimate is a functional
```math
\mathcal{R}_h\bigl(T(n\mathbb{P}_n)\bigr)
```
of that summary, where ``\mathcal{R}_h`` is a smoothing operator (a choice of kernel, bandwidth ``h``, and any
secondary smoothing). The load-bearing decision:

> **The summary ``T`` is accumulated online and is estimand-agnostic; the smoothing ``\mathcal{R}_h`` is applied
> entirely post hoc.** Accumulation knows nothing of bandwidth, kernel, or which response channel —
> which response variable, each with its own summary — will be asked for. This separation is the architectural invariant of the package — the analog of
> PF's representation/operator split — and it is exactly what delivers the user's requirement:
> *no state storage* (only ``T`` survives a step) and *no reruns* (every ``h``, every channel, every
> secondary smoothing is computed from the one stored ``T``).

### The algebra of measures — a structural lens

The spine has a compact algebraic reading, offered as a structural lens, not a second motivation:
the package is and remains mathematical statistics grounded in practice — its targets are
consistency, the ``n^{-2/5}`` convergence rates, the analytic oracles, honest bandwidth selection. The
lens only names *why* the spine takes the shape it does (established structure read back onto the code,
as in PF), and settles nothing the statistics does not already settle.

Read the summaries as living in the *cone of (unnormalized) measures*: a commutative monoid under
`merge!`, closed under nonnegative scaling by `scale!` (the cone of §10). The binning map ``T`` is then
linear, ``T(\mathcal{A} \uplus \mathcal{B}) = T(\mathcal{A}) \oplus T(\mathcal{B})`` for sample multisets ``\mathcal{A}``, ``\mathcal{B}`` — a cone homomorphism — and accumulation is the monoid
homomorphism from the multiset of samples (the free commutative monoid; `observe!` folds a batch,
`merge!` unions two summaries). *This is why* `merge!` is elementwise addition and accumulate-then-merge
equals accumulate-together: not a convenience, but linearity of ``T``.

The estimands sit one level up. Every reported quantity — the normalized density ``f_X``, the
conditional mean ``m(x) = \mathbb{E}[Y \mid X = x]``, the conditional density ``f_{Y\mid X}`` — is
homogeneous of degree zero in the measure: rescale the summary, ``T \mapsto \lambda \cdot T``, and it is unchanged, so
it lives not on the cone but on its projectivization, the probability measures. The
`scale!`-invariance §10 verifies (a `conditional_mean` unmoved by reweighting) *is* this degree-zero
homogeneity. And because the readout ``\mathcal{R}_h`` is degree zero while the cone operations are linear, ``\mathcal{R}_h``
is not a cone homomorphism — the estimand-agnostic-accumulation / post-hoc-smoothing split restated:
pooling must happen on the cone (linear, mergeable), and the degree-zero readout must stay post hoc. That
is *why* "no reruns" forces the smoothing post hoc.

The lens also fixes the shape of the estimators: the universal object is the summary-as-measure, and
every estimand is a degree-zero readout of it that falls into one of two families, each with one
shared code atom — kernel-weighting (the *local* / smoothing readouts, through the
`Kernels._kernel_weight` atom, with `LocalPolynomial.local_fit` the jet readout on it — the jet being the fit's local Taylor
coefficients) and power-sum
quadrature (the *global* / moment readouts, the Taylor contraction of the stored ``S_p``). §4.A and §4.B
exposit the two families. The jet `local_fit` is the kernel-weighting family's primitive, not the universal
one — ``\mathbb{E}[g(X)]`` is a readout with no jet at all — and the package's one nonlinearity (the
density's log link and normalization) lives in the readout, exactly where the cone homomorphism breaks.

### Why a fine histogram is the canonical ``T`` (and the honest caveat)

The canonical summary is a finely binned empirical measure — a histogram, the *fine-grained
measure* that post-hoc smoothing (a choice of bandwidth ``h``) later coarse-grains into each reported
estimate — optionally carrying per-bin moments (below). For nonparametric estimation there is, strictly, *no* finite
sufficient statistic; the histogram is lossy. The justification is quantitative and sets the
bin count: with bin width ``\Delta`` and smoothing bandwidth ``h``, the binning error is ``O(\Delta^2)`` while the
irreducible smoothing bias is ``O(h^2)``, so choosing ``\Delta \ll h`` (in practice ``\gtrsim 10^3`` bins over
``\pm 5\sigma_X``, so ``\Delta \approx \sigma_X/100 \ll h``) makes the summary *approximately sufficient*: it
loses nothing at the target accuracy (the WARPing / binned-KDE / average-shifted-histogram result,
Hall & Wand 1996; Scott 2015). We therefore say "approximately-sufficient summary," never
"sufficient statistic" — the latter would be the category overclaim PF's correctness policy exists
to catch.

### Estimand families (the operation taxonomy)

Mirroring PF's operation families (vertical algebra / differential / reductions / transfer), the
package's families are organized by estimand — the population quantity targeted. Each is
detailed in §4:

- **Marginal density** ``f_X`` — from a 1-D histogram.
- **Conditional expectation** ``m(x) = \mathbb{E}[Y \mid X = x]`` and conditional variance — from
  per-bin power-sums.
- **Conditional / transition density** ``f_{Y\mid X}``, ``p_\tau(y \mid x)`` — from a 2-D (joint)
  histogram, per lag.
- **Conditional increment moments and their functionals** *(archived, §4.4)* — Kramers–Moyal
  coefficients, drift/diffusion, Pawula / Q diagnostics — domain-neutral functionals of the above; the SDE
  *reading* stays consumer-side.
- **Bandwidth selection** (cross-cutting) — global and the response-aware local ``h^\star(x)`` (both built).
- **[OPEN]** a streaming projection road (basis projection à la Stochastic Force Inference (SFI) /
  Underdamped Langevin Inference (ULI)) as an alternative accumulator.

### The histogram is the summary; the ECF is deferred (likely never needed)

The empirical measure could also be carried in a Fourier representation — the empirical
characteristic function, the sample average of the complex exponential of each sample — but the histogram
is the sole summary. Tail estimation (the heavy-tailed acceleration / sixth-power increment channels) is
robust only in the histogram, where the ECF inverse rings; the ECF would earn its place only for a
strictly near-Gaussian marginal, which the consumers do not require. It is therefore deferred and,
per the user, most likely never built.

**This does not forgo the FFT** — and the distinction is the crux. Storing the ECF as a *sufficient
statistic* (an accumulation choice) is a different thing from using the FFT to *smooth the
histogram* (a post-hoc computation choice). The former is dropped outright; the latter is the
endorsed post-hoc optimization (§4.5), deferred (recorded, not built): smoothing a binned density or
regression is a convolution of the bin counts with the kernel, which the FFT does in ``O(B\log B)`` for
``B`` bins and — decisively — amortizes over bandwidth (transform the counts once; each bandwidth
``h`` is then a multiply by the kernel's transform). It is a *pure* optimization — the direct summation
is already exact and dwarfed by the accumulation (§5) — so it waits for a consumer whose ``h^\star(x)``
bandwidth *grid* makes it bite. So FFT-*smoothing* is the endorsed (if deferred) computation; the ECF
*summary* is rejected as a stored representation.

### Naming conventions

Following PF's mathematics-first naming:

- **Estimands / estimated quantities** take the quantity's noun: `density`, `conditional_mean`,
  `conditional_variance`, `conditional_density`, `transition_density`, `expectation`.
- **Operations on the summary** that are actions take the verb: `observe!` (fold a sample
  batch in), `merge!` (combine summaries), `scale!` (reweight the cone).
- **In-place** forms append `!` uniformly (`observe!`, `merge!`, `scale!`).
- **Standard named statistical objects** would keep their standard name (a `histogram`, a
  `characteristic_function`) — a forward policy with no current instance; note the local-polynomial fits
  are one `conditional_mean` with a `degree` option, not separate `nadaraya_watson` / `local_linear`
  bindings.

### Documentation policy — estimand first, estimator second

The mirror of PF's *physical-meaning-first* rule. Because every operation estimates a population
quantity from data:

> A docstring (and every doc-page entry) for an estimator states the estimand — the population
> quantity, in standard probability notation, with its defining formula — *before* the estimator
> that approximates it from the accumulated summary.

When an object is intrinsically finite-sample (a bandwidth selector, a `merge!`), say so plainly rather
than forcing a population form. The full convention — the `conditional_mean` reference example, the
symbol-binding, and the separate correctness/normalization check — lives with the documentation
conventions (§7).

---

## 3. The accumulated summary — the universal lower half

These types are the load-bearing, domain-agnostic core; they are pure numerics over a stream and
transfer verbatim to every consumer.

- **`Histogram1D`** — a binned empirical measure of a conditioning variable ``X``: linear (WARPing)
  bins over a fixed range (e.g. ``\pm 5\sigma_X``, ``\gtrsim 10^3`` bins), each bin carrying its
  count and per-bin power-sums about the bin center ``S_p = \sum_{i\in\text{bin}}(X_i -
  c_j)^p`` for ``p = 1 \dots \texttt{order}`` (default `order = 4`; raise `order` for the heavy-tailed channels — a high
  marginal central moment such as ``\mu_6(X)`` needs `order = 6`, and its *sampling variance* reaches the 12th
  moment; the *conditional* increment moments of the Pawula quotient ``Q`` (the fourth increment
  moment over thrice the squared second, §11) instead need the `Histogram2D`
  order-raise, §11). Because the
  offset ``X_i - c_j`` is small within a bin (``\lesssim \Delta/2 \approx \sigma_X/200``), the local
  power-sums are well-conditioned: they are accumulated by a plain atomic scatter-add — single-pass,
  numerically stable, with no compensated summation needed. The global central moments are
  recovered post hoc by the shift formula (the binscatter / WARPing local-moment construction;
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
  transition densities; per-cell counts and bivariate power-sums ``S_{pp'}`` about the cell center.
- **`HistogramND` [OPEN]** — the multi-dimensional generalization, deferred until a consumer needs conditioning
  on more than one variable (the PF rule against speculative structure).
- **per-bin power-sums are the one summary scheme** above, in a univariate form (``S_p`` of ``X``,
  for the marginal moments and density) and a bivariate form (``S_{pp'} = \sum_{i\in\text{bin}}(X_i -
  c^x)^{p}(Y_i - c^y)^{p'}`` about the cell center, the local sufficient statistic for local-polynomial
  regression; Cattaneo et al. 2024). All bandwidths, kernels, and moments are recovered from them
  post hoc.
- **the mergeable-monoid (cone) interface** — `merge!` on every summary is elementwise addition (same
  binning ⇒ same bin centers ⇒ counts and power-sums add); the weighted `merge!(into, from; weight)` and
  `scale!` extend it to the cone (§10). Tested: accumulate-then-merge ``=`` accumulate-together, and the
  cone laws.
- **the backend-agnostic scatter-add primitive** *is* `observe!` / `observe_pairs!` — one
  KernelAbstractions + Atomix kernel (one work-item per sample) that atomically scatters the count and
  the power-sums into the preallocated summary, the direct analog of PF's `_bin_by_shell` /
  `anterpolate!`. Allocation-free up to the kernel launch; verified identical on CPU and the AMD GPU
  under `allowscalar(false)` (§5). Storing power-sums *about the bin center* is precisely what keeps
  this a plain atomic add — no privatization, no workspace.

**Numeric precision — accumulate wider than the data.** Samples often arrive in Float32 (the
GPU-first, low-memory consumers), but the summaries do not accumulate in Float32:

- **Counts are integers** (`UInt64` merged). Float32 represents integers exactly only up to
  ``2^{24} = 16{,}777{,}216`` — beyond which adding one sample is silently lost — yet a single ``256^3``
  snapshot is *exactly* ``2^{24}`` points, and a run accumulates ``10^{10}``+ into the peak (modal) bins.
  Float32 counting is disqualified outright.
- **Power-sum accumulators default to Float64.** The local sums are well-conditioned, but the *global* moment
  recovery (the shift formula combines them with large bin-center offsets) and the high-order channels
  still want full precision; Float32 would discard the ~7 digits the Float32 *input* carries, through
  the reduction over ``n``. The power-sums accumulate in Float64 by default, regardless of input type.

Accumulator precision is a type parameter (PF-style element-type genericity) but its default is
Float64, counts integer — *accumulate-wider-than-the-data*, the deliberate opposite of PF's
preserve-the-field's-eltype (a reduction should out-precision its inputs; a transform should
not). It costs almost nothing: the histogram is small (a few ``\times10^3`` bins ``\times`` a few
Float64 power-sums — kilobytes per channel, versus a gigabyte field), and on the GPU the Float64
arithmetic is hidden behind the memory- and atomic-bound accumulation (§5), not FP64 throughput.

**Two-time samples / lag handling — settled: the package bins pairs; the consumer produces them.**
For a two-time quantity like ``m_\tau(x) = \mathbb{E}[Y(t+\tau)\mid X(t) = x]`` the conditioning
value ``X(t)`` and the response ``Y(t+\tau)`` must be paired per source — two separate marginal
histograms cannot recover the joint (for the Ornstein–Uhlenbeck (OU) oracle both marginals are the *same* stationary
Gaussian; the entire ``\tau``-dependence lives in the pair). The one histogram entry is made at the
response step, into the bin of ``X(t)``. The package surface for this is just the memoryless
`observe_pairs!(stat, x_past, y_now)` — it bins a batch of pre-formed pairs into the per-lag summary (the
lag is the summary's `Histogram2D(…; lag)` label, not an argument here; §4.2 / §4.3). (A pre-quantization
path — a public binning map plus an index-taking `observe_pairs!` that folds a compact buffer of x-bin
indices directly, lossless because the conditioning readers use only the x-cell and never its sub-bin
offset — is deferred; the `bin_index` helper is archived until it lands.)

Everything that *produces* the pairs is the consumer's, being sample-source policy a
domain-neutral package must not hold: the snapshot buffer of past ``X(t)`` (held until the response
arrives — the irreducible per-lag state, linear in the number of sources, optionally stored as compact
`UInt16` bin indices when memory is tight); which sources to sample; which
lags to form (e.g. log-spaced); and the origin cadence / buffer depth (a memory-vs-variance
budget). Pairing is by source identity — for a stably-ordered set just element-wise array alignment
(no id tag), an explicit id only if the ordering can change (re-meshing, sorting, birth/death).

What the package *keeps* is the per-lag histograms — the accumulated power-sums, the small
*output* (lags ``\times`` bins ``\times`` moments, a few MB), never the snapshot buffer of every source at every lag (the
expensive object, consumer-side). For the nonstationary case the running origin ``t`` and lag
``\tau`` are explicit labels on the summary set (``p_{t,\tau}``), never assumed away. (If the
ring-and-pair boilerplate proves worth sharing, it can be a non-domain-aware utility offered
alongside the summary — never something the package's vocabulary depends on.)

---

## 4. Estimand families (estimators — the exportable upper half)

Each estimand is written estimand first (per §2), grouped by the two readout families of the lens:
§4.A kernel-weighting (the local / smoothing readouts through `Kernels._kernel_weight` — the marginal
density §4.1, the conditional mean and variance §4.2, the conditional / transition density §4.3, and
bandwidth selection §4.5) and §4.B power-sum quadrature (the global moment readouts off the stored
``S_p`` — `mean` / `var` / `central_moment`, and `expectation` ``\mathbb{E}[g(X)]`` / `pushforward_density`,
all built). The §4.4 Kramers–Moyal / SDE-reading functionals are archived (in git history at commit `37006b5`, restore when a
consumer needs them and `Histogram2D` raises its per-cell order); §4.6 is an [OPEN] alternative.

### 4.A Kernel-weighting readouts

The local / smoothing readouts: each a degree-zero readout of the kernel-weighted measure through the
`Kernels._kernel_weight` atom (§2), written estimand-first.

#### 4.1 Marginal density ``f_X``

**Estimand.** ``f_X = \mathrm{d}\mathbb{P}_X/\mathrm{d}x``, the Lebesgue density of the law ``\mathbb{P}_X`` of ``X``:
``\mathbb{P}_X(A) = \int_A f_X(x)\,\mathrm{d}x``.

**Estimator** (binned KDE / WARPing): ``\hat f_h(x) = \frac{1}{nh}\sum_j n_j K((x-\tilde c_j)/h)`` (``\int K = 1``)
from counts ``n_j`` at the within-bin centroid ``\tilde c_j = c_j + S_1^{(j)}/n_j`` (`binning = :linear`,
the default — matching the bin's first moment, so the discretization is ``O(\Delta^2)`` in the bin width
``\Delta``; ``\tilde c_j = c_j`` is `binning = :simple`, ``O(\Delta)``), the centroid free from the stored
``S_1``. Bandwidth from §4.5; the kernel is pluggable. The tail and shape refinements are all post-hoc on the one
summary — the adaptive variable bandwidth (cutting tail *variance*), the local-log-polynomial
`degree` (cutting tail *bias*), a bounded support (`support`, or an asymmetric boundary kernel), a
parametric start, and exact `normalize` — with `density_standard_error` the pointwise variability band;
`FORMALIZATION.md` and the API page carry each one's estimand and constants. The axes are conceptually
independent (each cures a different defect of the plain KDE), but not every pair combines: the table after
this paragraph records which pairs do, which are refused as redundant, and which are refused only because
the combination is not built.
Caveats: the fixed range ``[\texttt{lower}, \texttt{upper})``
drops events beyond it (cover the data — ``\approx \pm 5\sigma_X`` for light tails, wider for heavy ones,
with proportionally more bins to keep ``\Delta \ll h``), and uniform bins under-resolve
power-law tails (the variable bandwidth is the preferred fix — it stays post-hoc, mergeable, and
estimand-agnostic, which non-uniform bins would not).

The `density` option pairs, as the guards decide them (audited 2026-10-15). "Redundant" means two cures
for the same defect, where the second has no theory behind it; "not built" means a definable estimator,
published or a direct composition, that the package has not implemented; every refusal is an
`ArgumentError` naming the pair.

| Pair | Status | Reason |
|:--|:--|:--|
| `degree` × `adaptive` | combine, at every degree | independent corrections: bias order and tail variance |
| `support` × `degree = 0` (with or without `adaptive`) | combine | the Jones renormalization of the KDE |
| `support` × `degree >= 1` | refused: redundant | the local-log fit is already boundary-adaptive, so the renormalization would be a second correction |
| `support` × asymmetric kernel | refused: redundant | the kernel's mass inside the support is exactly 1, so the renormalization would divide by 1 |
| asymmetric kernel × `degree >= 1` | refused: redundant, theory open | two boundary cures at once; the local-likelihood bias expansion rests on translation-kernel moments, which a varying-shape kernel lacks |
| `start` × `degree >= 1` | refused: not built | the locally parametric estimator of Hjort & Jones (1996), the named first candidate |
| `start` × `adaptive` | refused: not built | Abramson's variable bandwidth on the correction's KDE, a direct composition |
| `start` × asymmetric kernel | refused: not built | a gamma/beta-kernel KDE of the correction, natural for a positive-support start |
| `start` × `support` | refused: not built | definable; the correction ``f_X/f_\theta`` (density over start) near an edge needs study before a claim |
| asymmetric kernel × `adaptive` | refused: theory open | whether the square-root law cancels the leading bias of a varying-shape kernel is unproved |
| `normalize`, `binning` × anything | combine | orthogonal to the estimator choice |

#### 4.2 Conditional expectation ``m(x)`` and variance ``\sigma^2(x)``

**Estimand.** ``m(x) = \mathbb{E}[Y\mid X=x]`` (the ``L^2(\mathbb{P}_X)`` projection; §2) and
``\sigma^2(x) = \operatorname{Var}(Y\mid X=x) = \mathbb{E}[Y^2\mid X=x] - m(x)^2``.

**Estimator** (local-linear; Fan & Gijbels 1996): from the joint `Histogram2D`'s kernel-weighted power-sums
``U_p = \sum_i K_h(X_i-x)(X_i-x)^p`` and ``V_p = \sum_i K_h(X_i-x)(X_i-x)^p Y_i``,
```math
\hat m_h(x) = \frac{U_2 V_0 - U_1 V_1}{U_0 U_2 - U_1^{2}}
```
— the ``\beta_0`` of one `LocalPolynomial.local_fit` (§4.A atom); local-linear removes the design-density
bias ``2 m'(x) f_X'(x)/f_X(x)`` that the degree-0 Nadaraya–Watson ``V_0/U_0`` carries. The variance is the
bias-reduced residual-based Fan & Yao (1998) estimator (`method = :residual`, the default — a second
local fit of the squared residuals ``\widehat{\mathbb{E}}[(Y-\hat m)^2\mid X = x]``, removing the
curved-``\mathbb{E}[Y^2\mid X = x]`` smoothing bias the `:plugin` ``\widehat{\mathbb{E}}[Y^2\mid X = x] - \hat m^2``
carries). A one-step robust local M-estimator for heavy-tailed channels (Fan & Jiang 2000) is forthcoming.

The general reader `conditional_expectation` (``\mathbb{E}[g(Y)\mid X = x]``) subsumes both — a
per-``x``-column ``y``-quadrature of ``g`` (the §4.B atom) regressed in ``x`` by the same `local_fit`, so it
joins the two families — with `conditional_mean` its ``g = \mathrm{identity}`` case and the conditional
pushforward density `conditional_pushforward_density` its transformed-response reader (one
`conditional_expectation` of a composed kernel, no ``g^{-1}``). Each carries a leverage-based `conditional_expectation_standard_error`,
``\sqrt{\ell(x) \operatorname{Var}(g(Y)\mid X = x)}``, the conditional
analog of the marginal `expectation_standard_error`. `FORMALIZATION.md` pins the constants and rates.

#### 4.3 Conditional and transition density

**Estimand.** ``f_{Y\mid X}(y\mid x) = f_{X,Y}(x,y)/f_X(x)``, and the two-time transition density
``p_\tau(y\mid x) = f_{Y(t+\tau)\mid X(t)}(y\mid x)`` (nonstationary: ``p_{t,\tau}``).

**Estimator** (double-kernel; Rosenblatt 1969; Hyndman, Bashtannyk & Grunwald 1996):
```math
\hat f_{Y\mid X}(y\mid x) = \frac{\sum_i K_{h_x}(X_i-x) K_{h_y}(Y_i-y)}{\sum_i K_{h_x}(X_i-x)},
```
with a separate conditioning kernel `kernelx` (``K_x``, in ``x``) and response kernel `kernely` (``K_y``, in ``y``). For a
bounded / positive response — a droplet radius as ``Y`` — `kernely` may be asymmetric (gamma / beta;
Chen 1999, 2000), boundary-bias-free in ``y``, with `normalize` restoring ``\int\hat f_{Y\mid X}\,\mathrm dy = 1`` —
the conditional analog of the marginal density's asymmetric kernels. The conditioning kernel cancels in the
``f_{X,Y}/f_X`` ratio, so it stays symmetric (a bounded *conditioning* variable needs none — the
local-linear `conditional_mean` is boundary-adaptive, §4.2).

#### 4.5 Bandwidth selection (cross-cutting)

**Estimand.** The estimand is the AMISE-optimal bandwidth (AMISE = asymptotic mean integrated squared error), the ``h_{\mathrm{AMISE}}`` asymptotically minimising ``\operatorname{MISE}(h) = \int \mathbb{E}[(\hat f_h(x) - f_X(x))^2]\,\mathrm{d}x``
for a kernel ``K`` and local-polynomial degree ``d``.

**Estimator.** For a marginal the package uses the kernel- and degree-aware normal-reference rule
``h = C(K,d) \hat\sigma_X n^{-1/(2q+1)}`` (bias order ``q = 2`` at degree 0/1, ``q = 4`` at degree 2/3; the
kernel enters through one canonical factor ``[R(K)/\mu_2(K)^2]^{1/5}`` by quadrature — Silverman
``(4/3)^{1/5}`` the Gaussian/degree-0 case, Epanechnikov ``\approx 2.34``), plus the data-driven
Sheather–Jones `:sj` and DCT-based `:isj` (Botev, Grotowski & Kroese 2010) plug-ins, the fit-based
`:aic`, and the asymmetric-kernel Chen `:plugin` (scale ``b \sim n^{-2/5}``) — each a mergeable functional
of the counts (the eval kernel a clean ``\delta(K)`` rescale); leave-one-out CV is a recorded gap (binned
ties drive ``h \to 0``). The load-bearing one for conditional moments is the response-aware *local*
plug-in `bandwidth(::Histogram2D, x; rule = :local)`,
```math
h^\star(x) = \biggl[\frac{R(K) \sigma^2(x)}{\mu_2(K)^2 m''(x)^2 f_X(x) n}\biggr]^{1/5},
```
which depends on the response's variance ``\sigma^2(x)`` and curvature ``m''(x)`` — so each response
channel gets its own ``h^\star`` (a single shared bandwidth is provably suboptimal), all three of
``\sigma^2(x)``, ``m''(x)``, ``f_X(x)`` recovered from the *one* summary, post hoc. *(Performance: every
smoother here is a convolution of bin arrays with the kernel, so an FFT computes it in ``O(B\log B)``
amortized over the bandwidth grid the ``h^\star(x)`` sweep needs — Silverman 1982; Wand 1994 — a deferred
optimization (recorded, not built), and in any case dwarfed by the per-sample accumulation hot loop, §5.)*

### 4.B Power-sum quadrature readouts

The global (non-local) readouts: a smooth functional of ``X`` read straight off the stored power-sums — no
kernel, no query point. Within a bin, ``S_p^{(j)} = \sum_{i\in\text{bin}}(X_i - c_j)^p`` is the within-bin
``p``-th moment, so Taylor-expanding any smooth ``g`` about the center,
```math
\widehat{\mathbb{E}}[g(X)] = \frac{1}{n} \sum_j \sum_{p=0}^{\texttt{order}} \frac{g^{(p)}(c_j)}{p!}  S_p^{(j)}
+ O\bigl(\Delta^{ \texttt{order}+1}\bigr),
```
exact for polynomials of degree ``\le \texttt{order}`` (the remainder terminates) — this is the
recovered `mean` / `var` / `central_moment` (§3, the shift formula, the ``g = x^p`` case) — and
``O(\Delta^{ \texttt{order}+1})`` for general smooth ``g``. Built in `Expectations`:

- `expectation` (``\mathbb{E}[g(X)]``) for a general smooth `g`. The per-bin coefficients
  ``g^{(p)}(c_j)/p!`` are read off a degree-``\texttt{order}`` polynomial interpolant of ``g`` on each
  bin — no automatic differentiation and no new dependency (uniform bins ⇒ one Chebyshev–Vandermonde
  inverse serves every bin), which reduces *exactly* to the shift formula on the monomials. Both
  `expectation` and the recovered moments route through the one shared contraction atom
  `Expectations._quadrature` (the moments supplying closed-form coefficients, `expectation` interpolation
  coefficients) — the literal realization of the two-family lens.
- `expectation_standard_error`, ``\sqrt{(\mathbb{E}[g^2] - \mathbb{E}[g]^2)/n}`` — the sampling error of
  ``\widehat{\mathbb{E}}[g]``, itself a second quadrature readout (``\mathbb{E}[g^2]`` is `expectation` with
  the smooth integrand ``g^2``).
- `pushforward_density` (``f_{g(X)}``), the density of ``Z = g(X)``. Rather than form ``g^{-1}`` or a
  Jacobian, it is written as the expectation of the composed kernel test function,
  ``\hat f_{g(X)}(z) = \tfrac{1}{n}\sum_i K_h(g(X_i) - z) = \widehat{\mathbb{E}}[K_h(g(X) - z)]`` — so it is one
  more `expectation` call. A non-monotone ``g`` is handled automatically (every preimage contributes),
  a flat spot ``g'=0`` shows up as the expected (integrable) caustic, and a symmetric ``K`` makes
  ``\int \hat f_{g(X)}\,\mathrm dz = 1`` exactly.

Limits: ``g`` smooth on the bin scale, mass inside the range, and — for the pushforward — the response
bandwidth ``h \gtrsim |g'| \Delta`` so ``K_h\circ g`` is resolved across a bin.

### 4.6 Projection road [OPEN, alternative]

An alternative to bin-and-smooth for the conditional-functional / SDE-inference channels (the archived §4.4
road, or general regression fields): estimate by projecting onto a function basis with running
inner-product accumulators (the Gram matrix and projection vector themselves a mergeable summary) —
Stochastic Force Inference (Frishman & Ronceray 2020) and the underdamped ULI (Brückner, Ronceray &
Broedersz 2020), fitting on samples of ``X`` paired with its time derivative. Shares the lower-half monoid foundation; deferred until a
consumer wants parametric-in-basis fields.

---

## 5. Backend-agnosticism and the hot-loop contract

Inherited wholesale from PF (its invariants 6 and 13):

- **Storage is generic `AbstractArray`**; no scalar indexing of device arrays in package methods;
  allocate via `similar`.
- **Accumulation is allocation-free.** `observe!(summary, batch)` allocates only a bounded amount
  (the kernel launch) on a preallocated summary — never an ``O(n)`` or ``O(\text{bins})`` temporary,
  which on a GPU hot loop would escape the host GC and wedge the run. Enforced by an `@allocated` audit
  on the CPU (where a materialized array's bytes are visible), the direct analog of PF's
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
  pattern, verified on the AMD GPU). Privatization (per-block local histograms in fast memory,
  merged hierarchically via the monoid) stays available as a future optimization should a channel's
  modal-bin contention dominate, but is not needed for correctness. This per-sample accumulation over
  every grid point at every step is the package's true performance-critical path; the post-hoc smoothing (§4.5) is
  second-order.

Getting backend-portability and mergeability right from the first commit is the one piece of
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
   unnormalized summary is additive over the stream, and scalable (`scale!`), so the summaries form a
   convex cone (equal-weight `merge!` the built special case; the weighted `merge!` and `scale!` the cone
   operations — §10 hyperstatistics / forgetting). No non-mergeable accumulator is admitted.
4. **The hot path is allocation-free and backend-agnostic** (§5).
5. **The summary is rich enough to recover every target post hoc.** A target that would need a
   finer or different summary is a recorded coverage gap, not silently approximated — the "no
   reruns" promise made falsifiable.
6. **Accumulate wider than the data.** Counts are integer; moment accumulators default to Float64
   even on Float32 input (§3). A Float32 count saturates (``2^{24}``); a Float32 high-order moment
   loses the data's own precision through the reduction. (The deliberate opposite of PF's
   preserve-the-eltype, for the opposite reason: this is a reduction, not a transform.)

**Deliberate absences (consumer-side, by scope).** Bulk scalar time-correlation functions,
``\mathbb{E}[X(t) X(t+\tau)]`` as a function of the lag, are not provided: a `Histogram2D` per lag would store the
whole joint density to extract one cross-moment (a full two-dimensional histogram per lag for one scalar per lag), and a
scalar correlation is cheaply streamed by a dedicated lagged product-sum or FFT. The
`Histogram2D`-per-lag path is for the *conditional / two-time density* objects — ``p_\tau(y \mid x)``,
the Kramers–Moyal coefficients, the van Hove function — that genuinely need the joint, over a typically *sparse* set
of lags; the consumer owns the scalar-correlation fast path.

---

## 7. Documentation, correctness, and testing conventions

Mirrors PF; the deltas are noted.

- **Estimand-first docstrings** (§2) — state the estimand (the population quantity, in standard notation,
  with its formula) before the estimator. `conditional_mean` is the reference: lead with
  ``m(x) = \mathbb{E}[Y \mid X = x] = \int y  f_{Y\mid X}(y\mid x)\,\mathrm{d}y`` (the ``L^2`` projection),
  *then* the local-linear ``\hat m_h``. Self-contained pages, rendered LaTeX in standard notation, every
  code object bound to its math symbol on first mention, in apposition (`density` (``f_X``),
  `conditional_mean` (``m(x)``), `transition_density` (``p_\tau(y\mid x)``)).
- **Every estimator reports its standard error.** An estimator approximates a population quantity from a
  finite stream, so it is incomplete without the sampling error of that approximation: every
  population-quantity estimator ships a companion `*_standard_error` (and the one-call paired
  `*_with_standard_error`) — a mergeable functional of the *same* summary, held to the
  `check_standard_error` contract. The exemptions are exactly estimand-first's: an intrinsically
  finite-sample object with no population estimand (`bandwidth`, a selector; `merge!` / `scale!`).
- **A `FORMALIZATION.md` analog.** PF records a typed statement per operation, expansion-checked
  against the code and test-linked. Here the typed entry states the estimand, the
  estimator, and the consistency claim — its bias and variance, e.g.
  ``\operatorname{bias} \hat m_h(x) \approx \tfrac{1}{2} h^2 \mu_2(K) m''(x)`` and
  ``\operatorname{Var} \hat m_h(x) \approx \frac{R(K) \sigma^2(x)}{n h f_X(x)}`` for local-linear
  regression — every normalization constant pinned, each linked to the test that verifies the
  convergence rate numerically.
- **Testing against analytic oracles** — the statistics analog of PF's analytic round-trips:
  estimators are checked against distributions with closed forms (Gaussian marginals; an
  Ornstein–Uhlenbeck process, whose drift, diffusion, and transition density
  ``p_\tau(y\mid x)`` are known exactly — the natural end-to-end oracle for §4.3–4.4), plus
  consistency / rate tests (error falls at the theoretical ``n^{-2/5}`` as ``n`` grows),
  mergeability tests, the `@allocated` audit, and GPU-vs-CPU parity under
  `allowscalar(false)`.
- **The documentation checks are the shared tooling's** — the rulebook of
  [GBGCodingTools](https://grigorysarnitsky.github.io/GBGCodingTools.jl/), run in every docs build, and its
  documentation review, read by a statistician and a software architect. What the package adds is the
  consistency and normalization of every estimator against its estimand — each estimator consistent
  for its estimand, every normalization exact (the ``1/n``, ``1/(nh)``, ``\int K = 1``, ``\mu_2(K)``,
  ``R(K)``); a stray ``1/h`` or kernel-moment factor is the easy slip, separate from formatting.

---

## 8. Proposed repository layout (mirrors PF)

```
GBGStreamingStatistics/
  CLAUDE.md                      # operational summary + invariants (mirror PF's)
  README.md  Project.toml
  src/
    GBGStreamingStatistics.jl    # umbrella
    Accumulators.jl              # foundation: AbstractAccumulator, the
                                 #   observe!/merge! monoid interface
    Histograms.jl                # Histogram1D / Histogram2D, WARP binning,
                                 #   power-sums about the bin center, the
                                 #   KA/Atomix scatter kernel (the pure summary)
    Kernels.jl                   # all the smoothing kernels (symmetric +
                                 #   asymmetric) + _quad quadrature
    Expectations.jl              # §4.B   𝔼[g(X)], moments (mean/var/…),
                                 #   pushforward — power-sum quadrature
    LocalPolynomial.jl           # the shared local-polynomial jet primitive
                                 #   (local_fit) — §2, §11
    MarginalDensity.jl           # §4.1   f_X
    ConditionalExpectations.jl   # §4.2   m(x), σ²(x)  [estimand-first reference]
    ConditionalDensity.jl        # §4.3   f_{Y|X}, transition p_τ
    Bandwidth.jl                 # §4.5   global + response-aware local plug-in
    # IncrementMoments.jl        # [archived → git 37006b5] §4.4 increment
                                 #   moments, KM/drift/diffusion
    # CharacteristicFunctions.jl # [OPEN, secondary] ECF summary
    # ProjectionInference.jl     # [OPEN] SFI / ULI road
  docs/   ARCHITECTURE.md  FORMALIZATION.md  TESTING.md  PERFORMANCE.md
          src/…  make.jl
  test/   runtests.jl  per-module  oracles.jl(OU/Gaussian)  mergeability.jl
          allocations.jl  aqua.jl  jet.jl  gpu/
  benchmark/
  examples/
  .gbg/config.jl  .claude/skills/   # the GBGCodingTools wiring (config,
                                    #   rulebook, review skills)
```

Internal split, made explicit: `Accumulators` + `Histograms` are the universal lower half
(every domain shares them verbatim, and they hold the hard, write-once numerics — the GPU
atomic scatter-add of counts and per-bin power-sums); `MarginalDensity` through `Bandwidth` are the estimator upper half
(turbulence-leaning but genuinely cross-domain — molecular dynamics' van Hove function is `ConditionalDensity` with
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
  separate, deferred post-hoc optimization (§4.5) — a computation choice, not a stored representation.
- **Lag handling — settled:** the package bins pre-formed pairs (`observe_pairs!`); the pre-quantization /
  index-ingestion path (`bin_index`) is deferred (archived); the snapshot buffer and all lag / sampling /
  cadence policy are consumer-side (§3).
- **HDF5 persistence — settled:** a weak-dependency extension (`write_histogram` / `read_histogram`,
  and the `parent[name] = summary` sugar), not a core dependency — a summary maps to a self-describing HDF5
  *group* (arrays as datasets, scalars as attributes, a type/version tag), reloaded exactly. Handle-based:
  the consumer opens and owns the file, so a summary sits beside other data in one file (many groups); the
  package never opens or closes a file. The summary is the durable artifact of the no-storage design, so
  checkpoint / resume / merge-across-runs are the motivating workflows.
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
- **Weighted merge — settled, built (the cone).** `scale!` and the weighted
  `merge!(into, from; weight)` make the unnormalised summary a convex cone (closed under nonnegative
  linear combinations), not merely an additive monoid — realizing hyperstatistics (superstatistics;
  Beck & Cohen 2003) and exponential forgetting as post-hoc reweightings, over a `count_type = Float64`
  summary variant. See the §11 decision-log entry.

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
  deferred post-hoc optimization (§4.5) — a computation choice, not a stored representation.
- **Lag handling — the package bins pre-formed pairs** (`observe_pairs!`); the pre-quantization /
  index-ingestion path (`bin_index`) is deferred (archived); the snapshot buffer and all lag / sampling /
  cadence policy stay consumer-side (§3).

- **Summary storage — power-sums about the bin center, not Pébay-about-the-running-mean** (§3).
  *Changed during the GPU pass.* The original design accumulated per-bin central moments by the Pébay
  (2008) running-mean recursion with Kahan compensation. A Pébay update is not atomic, so the GPU
  hot loop would have needed per-block *privatization* (private histograms merged hierarchically) — a
  pattern the sibling PF uses nowhere. Switching to power-sums ``S_p = \sum (X_i - c_j)^p`` about the
  fixed bin center makes `observe!` a plain atomic scatter-add (PF's proven `anterpolate!` pattern):
  allocation-free, one portable kernel identical on CPU and the AMD GPU, `merge!` reduced to
  elementwise addition, and the marginal moments unified with the regression power-sums. The offset
  ``X_i - c_j`` is tiny, so the local sums are well-conditioned (Kahan no longer needed); global central
  moments are recovered post hoc by the shift formula — exact, and stable for the even high-order
  marginal moments the consumers need (odd moments keep weaker relative precision, an accepted trade).
  Verified by the ROCm GPU parity suite.

- **Kramers–Moyal §4.4 — drift and diffusion first, as conditional increment moments** (§4.4).
  *[Superseded — these functionals were built and OU-tested, then archived; see the archival entries
  below.]* The first slice implements the increment moments and the finite-``\tau`` Kramers–Moyal coefficients
  drift (``M_1/\tau``) and diffusion (``M_2/(2\tau)``) directly from the order-2
  `Histogram2D`, via ``M_1(x) = \hat m_h(x) - x`` and
  ``M_2(x) = \widehat{\mathbb{E}}[Y^2\mid X = x] - 2x \hat m_h(x) + x^2`` — reusing the §4.2 local-linear
  primitives (no new accumulator), reported as the plain finite-``\tau`` ratio ``M_n/(n! \tau)``. The
  all-orders Bell-polynomial finite-time correction (Rydin Gorjão et al. 2021) and the ``\tau\to0``
  limit stay forthcoming / consumer-side. Orders three and above — the higher Kramers–Moyal coefficients and
  `pawula_quotient` (whose ``Q = M_4/(3M_2^2)`` needs ``M_4``) — are a recorded coverage gap (they
  need the joint summary's per-cell ``y``-power-sums raised to that order, deferred per invariant 5);
  `markov_test` composes the now-available §4.3 transition density.
  OU-oracle-tested: drift recovers the finite-``\tau`` value ``x(e^{-\kappa\tau}-1)/\tau`` exactly (no
  smoothing bias — the OU conditional mean is linear — ``\to -\kappa x`` as ``\tau \to 0``), diffusion
  tends to half the squared noise amplitude, carrying the local-linear bias of the curved ``\mathbb{E}[Y^2\mid X = x]``.

- **Conditional / transition density §4.3 — binned double-kernel from the joint counts** (§4.3). The
  Rosenblatt double-kernel estimator (§4.3) is evaluated from the `Histogram2D` cell counts (no
  power-sums needed): at the conditioning value ``x`` the joint is reduced to an ``x``-kernel-weighted
  y-profile ``P_l = \sum_j n_{jl} K_x\bigl((c^x_j - x)/h_x\bigr)`` with ``D = \sum_l P_l``, and then
  ``\hat f_{Y\mid X}(y\mid x) = \sum_l P_l K_y\bigl((c^y_l - y)/h_y\bigr) / (h_y D)`` for symmetric kernels ``K_x``, ``K_y`` (any
  ``\int K = 1``, default Gaussian). The ``x``-kernel weighting cancels in the conditional ratio — only
  the ``1/h_y`` survives, and ``\int K = 1`` gives ``\int \hat f_{Y\mid X}(y\mid x)\,\mathrm dy = 1`` exactly.
  `transition_density` is the same computation on a lag-``\tau`` joint (the two-time reading). The two
  bandwidths ``h_x`` (`hx`, conditioning) and ``h_y`` (`hy`, response) and the kernels `kernelx`, `kernely` are supplied; the
  conditional MISE rate is the slower ``O(n^{-2/3})``.
  OU-oracle-tested against the closed-form Gaussian ``p_\tau(y\mid x)``.

- **Bandwidth §4.5 — the marginal `bandwidth` made kernel- and degree-aware** (§4.5). *Generalised the
  global rule.* The first slice returned only Silverman ``(4/3)^{1/5}\hat\sigma_X n^{-1/5}`` (Gaussian,
  degree 0; rules `:silverman`/`:scott`). A survey of the field (locfit, KernSmooth, ks, statsmodels)
  confirmed the AMISE-optimal ``h_{\mathrm{AMISE}}`` depends on the kernel only through the canonical factor
  ``\delta(K) = [R(K)/\mu_2(K)^2]^{1/5}`` and on the local-polynomial degree through both the rate
  (``n^{-1/5}`` at degree 0/1, ``n^{-1/9}`` at degree 2) and the equivalent-kernel constant — exactly
  locfit's `Wikk(K, d)`. `bandwidth(summary; rule = :normal, kernel, degree)` now returns
  ``C(K,d) \hat\sigma_X n^{-1/(2q+1)}`` (``R(K)``, ``\mu_2(K)`` by quadrature, so any kernel works) and
  is also `density`'s default `h`, so a degree-2 or Epanechnikov call is matched automatically. The
  kernel-as-standard-deviation shortcut (base-R, KDEpy, KernelDensity.jl — kernel-blind, not
  AMISE-exact) was weighed and declined. The data-driven Sheather–Jones plug-in (`rule = :sj`, degree
  0/1) is now added — the ``\psi_4``/``\psi_6`` curvature functionals as mergeable Hermite-weighted
  reductions of the bin counts, the solve-the-equation root-find (base-R `bw.SJ` ste), the eval kernel a
  clean ``\delta(K)`` rescale (the pilot stays Gaussian, kernel-independent) — much tighter than the
  normal reference on a bimodal sample. The fit-based AIC follows (next entry).

- **Bandwidth §4.5 — fit-based AIC added; leave-one-out CV is a deliberate absence** (§4.5). The survey
  flagged Loader's *preferred* route as fit-based selection over an ``h``-grid. Implemented as
  `rule = :aic`: minimise the binned-KDE AIC ``-2\sum_j n_j\log\hat f_h(c_j) + 2\operatorname{df}(h)`` (``\operatorname{df}`` the effective degrees of freedom, the smoother's self-influence trace) over
  the grid — a mergeable functional of the counts, the eval kernel a ``\delta(K)`` rescale. But the
  cross-validation selectors (likelihood-CV and least-squares LSCV) are not built: they are
  leave-one-out, and a *binned* summary has heavy ties (``\approx 10^2\text{–}10^3`` samples per bin) whose same-bin term
  drives the CV objective to ``h \to 0`` (verified — the coefficient of LSCV's leading ``1/(nh\sqrt\pi)`` term goes
  negative, so it diverges to ``-\infty`` as ``h \to 0``). base-R `bw.ucv` survives only via coarse
  bins, a constrained search range, and a boundary warning. AIC needs no leave-one-out — both its
  log-likelihood and df-penalty diverge as ``h \to 0`` — so its optimum is a stable interior point. CV
  recorded as a coverage gap (a within-bin-spread correction from the power-sums could rescue it);
  `:sj`, `:isj`, and `:aic` are the data-driven selectors.

- **Bandwidth §4.5 — Improved Sheather–Jones (`rule = :isj`) added** (§4.5). The classic `:sj` plug-in
  estimates the curvature ``\psi_4`` from a *normal-reference* pilot, which can mislead on strongly
  multimodal data (the Marron–Wand claw). Botev, Grotowski & Kroese (2010) instead read the curvature
  norms ``\|f_X^{(r)}\|^2`` off the DCT of the bin counts and solve the AMISE fixed point — a fully
  nonparametric pilot. Added as `rule = :isj`: a mergeable functional of the counts' cosine transform,
  the eval kernel the same ``\delta(K)`` rescale, robust where `:sj`'s pilot misleads.

- **Density §4.1 — parametric start and exact normalization added** (§4.1). A second pass over the
  downloaded reference software (studying source, not docs) found two capabilities worth mirroring: the
  parametric start (Hjort & Glad 1995; `kdensity`'s headline feature) and the exact
  renormalization (`kdensity`'s `normalized`). The start is the semiparametric
  ``\hat f_X = f_{\hat\theta}\cdot\widehat{f_X/f_\theta}`` — a moment-fitted model times a KDE of its
  correction — with five moment-fittable starts (`normal_start`, `lognormal_start`, `gamma_start`,
  `exponential_start`, `inverse_gaussian_start`); it stays in-discipline because ``\hat\theta`` is a
  functional of the recovered moments, so the whole estimate is mergeable, and it is a *bias* reduction
  (the leading bias is driven by ``(f_X/f_\theta)''``), valuable for the skewed positive shapes the package targets
  (droplet radii). `normalize=true` divides any estimate by its integral over the support so
  ``\int\hat f_X = 1`` exactly — the finishing step for the non-self-normalizing asymmetric kernels and
  starts. Both are post-hoc degree-0 options; `start` does not compose with ``\texttt{degree} \ge 1``, `adaptive`,
  `support`, or the asymmetric kernels (the local-likelihood-with-start generalization is a recorded gap).
  The cross-validated start bandwidth (`KernelEstimator.jl`) and the Hjort–Jones *locally parametric*
  generalization were weighed and deferred. The same source pass also corrected the gamma/beta kernels
  from Chen's naive estimator to the bias-reduced one (a fix, committed separately).

- **Bandwidth §4.5 — Chen plug-in for the asymmetric kernels; AIC, U-statistic, and gamma-reference
  routes rejected** (§4.5). The asymmetric (gamma / beta) kernels shipped only the ``\hat\sigma_X n^{-2/5}``
  rule of thumb; a data-driven selector was wanted because the optimal scale is ``0.27\text{–}1.6\times`` it
  across densities (``> 2\times`` ISE on log-normal — the droplet regime). Three routes were tried and
  rejected before the plug-in: the binned AIC (as for the symmetric kernels) is *monotonic* in ``b``
  — the gamma self-influence ``\sim 1/\sqrt{xb}`` makes the df penalty too weak against the flat binned
  likelihood of dense data, so it oversmooths without bound; the ``\psi_4`` U-statistic curvature
  (the ``:sj`` machinery) assumes ``f_X`` vanishes at the boundary (integration by parts) and goes
  *negative* on the half-line where ``f_X(0) > 0`` (e.g. the exponential); a gamma-reference plug-in is
  pathological for shape ``\alpha \le 1`` (a non-integrable ``f_X''`` boundary singularity, ``b^\star \to 0``).
  The shipped `rule = :plugin` is the Chen plug-in with the curvature ``\int w(x)^2 f_X''(x)^2\,\mathrm{d}x`` (weight
  ``w = x`` / ``x(1-x)``) estimated by grid quadrature of the pilot ``\hat f_X''`` — boundary-valid
  because ``w \to 0`` there, and always real. The ``f_X''`` pilot is a normal-reference scale
  ``1.3 \hat\sigma_X n^{-1/9}`` (constant calibrated for robustness, as Silverman's 1.06 is): never
  materially worse than the rule of thumb, ``1.4\text{–}1.7\times`` lower ISE on skewed data. A
  fully-adaptive (solve-the-equation) curvature pilot is a recorded refinement.

- **Hyperstatistics — the cone: `scale!` + weighted `merge!` over a Float-count variant** (§10). *Built.*
  The unnormalised summary is scalable, so the summaries form a convex cone, not merely an additive monoid.
  Realized as `scale!(summary, weight)` (the scalar action — every count and power-sum ``\times`` weight)
  and the weighted `merge!(into, from; weight)`
  (``\texttt{into} \mathrel{+}= \texttt{weight} \cdot \texttt{from}``; ``\texttt{weight} = 1`` recovers the
  plain monoid add), on a Float-count variant (`count_type = Float64`; the default `UInt64` counts cannot
  hold a fractional ``\lambda n``, so `scale!` and a weighted merge `ArgumentError` on them). The summary storage
  needed almost no change — `observe!`'s atomic scatter already adds `one(eltype(counts))`, so a Float count routes
  through the atomic-*float* add (CPU and GPU), and `merge!` was already an elementwise broadcast. Estimand:
  the hyperstatistical parameter-mixture density
  ``f_X(x) = \int \pi(\theta) f_\theta(x)\,\mathrm{d}\theta \approx \sum_k \pi_k f_{\theta_k}(x)``, read off the mixed summary
  ``T = \sum_k (\pi_k/n(T_k)) T_k`` with ``n(T_k)`` the total count of ``T_k``
  (superstatistics; Beck & Cohen 2003) — keep one summary ``T_k`` per ``\theta_k``, choose the mixing law
  ``\pi(\theta)`` post hoc, re-weight without re-accumulating; and exponential forgetting (`scale!(s, λ)`
  before the next block). Oracle-tested by the ``\chi^2``-superstatistics → Student-``\mathrm{t}`` result (a Gamma-weighted
  Gaussian superposition has the heavy ``\mathrm{t}`` tail), a discrete Gaussian mixture, the cone laws, exact
  forgetting, the `Histogram2D` cone, and GPU parity (`test/hyperstatistics.jl`, `FORMALIZATION.md`). One
  caveat recorded: the mixture's total mass is ``\sum_k \pi_k`` (``\approx 1`` for a normalized law), so its
  *default* bandwidth reads ``n \approx 1`` — pass an explicit ``h`` (or rescale the weights to a
  representative sample size). An effective-sample-size bandwidth for weighted summaries
  (Kish ``(\sum_i \lambda_i)^2 / \sum_i \lambda_i^2``) is a possible refinement, not built.

- **Linear (centroid) binning — the default for `density`** (§4.1). *Built.* The binned KDE places each
  bin's mass at a representative location ``\tilde c_j``. The geometric center ``\tilde c_j = c_j`` (simple binning)
  leaves a nonzero within-bin first moment, an ``O(\Delta)`` discretization of the exact KDE; placing the
  mass at the recovered within-bin centroid ``\tilde c_j = c_j + S_1^{(j)}/n_j`` zeroes that first moment
  (``\sum_{i\in j}(X_i - \tilde c_j) = 0``), so the leading term cancels and the discretization error drops to
  ``O(\Delta^2)`` — the centroid analog of Hall–Wand linear binning. The first power-sum ``S_1`` is *already
  stored*, so the refinement is free, and it is the new default (`binning = :linear`; `binning = :simple`
  keeps the center). One change threads it everywhere: every `density` path reads a `centers` vector, so the
  KDE, the local-log fit, the parametric start, and the asymmetric kernels all inherit the centroid at no
  extra cost. Oracle-tested by the binning-error convergence sub-test (`test/density_convergence.jl`: a
  log–log slope ``\approx 2`` in ``\Delta`` for `:linear` against ``\approx 1`` for `:simple`, the linear
  estimate an order of magnitude tighter at a fine grid). A higher-order correction using the stored ``S_2``
  is conceivable but not built.

- **Estimator architecture — the linear jet-smoother is the primitive, local likelihood a consumer;
  and the measure-algebra lens** (§2, §4). *A framing decision, landed in the code.* Two
  settlements taken together.

  **(1) The lens.** The spine is given an algebraic reading (§2): the
  summaries are the cone of measures, ``T`` is the linear cone homomorphism (so `merge!` is
  addition), and every reported estimand is homogeneous of degree zero in the measure, descending
  to the projectivization (the probability measures) — exactly the `scale!`-invariance §10 tests, and
  exactly why the post-hoc readout is not a homomorphism. Adopted as a *structural lens only*; the
  motivation stays mathematical statistics grounded in practice (consistency, rates, oracles).

  **(2) The primitive.** Two abstractions were weighed for sharing code between `density` and
  `conditional_mean`: a local generalized linear model (local likelihood — `density` at degree
  ``\ge 1`` is a local Poisson/log fit, `conditional_mean` a local Gaussian/identity fit; Loader 1996,
  Fan & Gijbels 1996), versus a linear local-polynomial jet-smoother (smooth the binned measure,
  return the local Taylor coefficients ``(\beta_0, \beta_1, \dots)``; conditioning a ratio of two such).
  They are not competitors: a local GLM is fit by iteratively reweighted least squares, so *each*
  Newton step is itself a kernel-weighted local-polynomial solve, and the identity-link case
  (regression) converges in one step to the plain solve. The GLM is therefore *built from* the
  linear solve, not parallel to it. So the linear jet-smoother is the irreducible primitive
  (`local_fit`), and the GLM is one of three thin consumers: a link-cap (the
  density's ``e^{\beta_0}``, the IRLS loop at degree ``\ge 1``), a ratio (conditioning —
  `conditional_density`, the disintegration of the joint), and a jet component (the derivative
  functionals `drift`, `diffusion`, and ``h^\star(x)``). This keeps the lone nonlinearity on one call
  site and the hot, allocation-free, GPU-portable path linear. *Landed (Stages 1–3):* the shared
  `LocalPolynomial.local_fit` (the jet-smoother, returning ``(\beta_0, \beta_1, \beta_2)``) now backs
  `conditional_mean` (its ``\beta_0``), the local-likelihood `density` (the IRLS inner solve), and the
  increment moments — replacing the formerly duplicated `ConditionalExpectations._local_linear` (degree-1, intercept
  only) and the bespoke Newton step in `MarginalDensity._local_log`. On it: `conditional_mean` /
  `conditional_variance` gained the `degree` (0/1/2) and `kernel` options; `conditional_variance` the
  bias-reduced Fan & Yao default (`method = :residual`); and the response-aware local bandwidth
  ``h^\star(x)`` = `bandwidth(::Histogram2D, x; rule = :local)`, which reads the curvature
  ``\hat m''(x) = 2\hat\beta_2`` of a degree-2 fit together with ``\sigma^2(x)`` and ``f_X(x)``. It landed as the
  substrate of that feature, not as a standalone refactor (convenience is not a justification, §6). Settled by estimation theory
  (iteratively reweighted least squares is the linear solve iterated; the smoother is linear on the
  cone), not implementation convenience.

- **Kernels — the asymmetric kernels consolidated into the foundation module** (§4.1, §4.3). The
  varying-shape boundary kernels (`gamma_kernel`, `beta_kernel`, the `AsymmetricKernel` type; Chen 1999,
  2000) were moved from `MarginalDensity` to `Kernels`, joining the symmetric `gaussian_kernel` /
  `epanechnikov_kernel`. They had stayed in `MarginalDensity` only because, when `Kernels` was first
  extracted, `density` was their sole consumer — a historical artifact, not a design boundary. A kernel
  is a kernel: an asymmetric boundary kernel is a smoothing kernel for bounded / positive support, so it
  belongs in the shared foundation, where any density estimator — marginal or conditional — can reach it.
  This gives the `kernel` argument one documented home (no cross-module split) and enabled the now-built
  asymmetric-kernel conditional density (§4.3: a bounded *response* — e.g. a droplet radius as ``Y`` —
  takes a boundary-free ``y``-kernel via `conditional_density(…; kernely = gamma_kernel, normalize = true)`,
  the same first-order boundary-bias fix the marginal density already has). It does not motivate asymmetric *conditioning* in `conditional_mean`: local-linear regression
  is boundary-adaptive (Fan & Gijbels 1996), already correcting the ``O(h)`` boundary bias to ``O(h^2)`` with a
  symmetric kernel, so an asymmetric conditioning kernel would be a boundary-*variance* refinement only,
  not a correctness fix — a recorded possibility, not built (convenience / parity is not a justification, §6).

- **Kernels / LocalPolynomial — the kernel-weighting atom `_kernel_weight`; the lens corrected to two
  readout families** (§2, §4). The symmetric/asymmetric kernel protocols, previously re-coded inline at four
  sites (`local_fit`, the density KDE, its asymmetric branch, the conditional density), are unified in one
  `Kernels._kernel_weight` (the ``1/h`` folded in, so every readout — count sum, Nadaraya–Watson ratio, jet — is uniform;
  the fold is ``\beta``-invariant for the local-polynomial solve, so symmetric results are bit-unchanged).
  `local_fit` carries the one guard: an asymmetric kernel has no jet, so it is admitted at ``\texttt{degree} = 0``
  only (a clean `ArgumentError`, not the former `MethodError`). This made degree-0 `conditional_mean` accept
  asymmetric kernels uniformly — building the degree-0 conditioning-side asymmetric kernel the entry
  above had left "not built". Honest finding: the Chen gamma/beta kernels are a *density* boundary tool;
  their asymmetry interacts with the regression design density, so asymmetric-kernel Nadaraya–Watson does
  not generally improve the conditional mean (degree-1 local-linear — boundary-adaptive, unbiased for a
  linear mean — is the regression tool). The capability is admitted for uniformity of the abstraction,
  not as the recommended regression estimator; the code is uniform, the guidance separate (in
  `conditional_mean`). The episode also corrected the lens: `local_fit` is the kernel-weighting family's
  primitive, not *the* universal one — ``\mathbb{E}[g(X)]`` is a power-sum-quadrature readout with no jet,
  the second family. Settled by the measure-algebra lens, not convenience.
- **IncrementMoments §4.4 — archived** (§4). The SDE-reading drift / diffusion / Kramers–Moyal functionals
  (working, OU-tested; Pawula / Markov were stubs) moved to `archive/` to focus the package on its universal
  statistical core — they are thin functionals of `conditional_mean`, the most physics-flavored,
  least-complete, most consumer-side part, and the identity is sharper without them. Restore when a consumer
  needs them and `Histogram2D` raises its per-cell order for the Kramers–Moyal orders three and above.
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
  `central_moment` were folded as the ``g = x^p`` case (guarded by the unchanged exactness tests) — the
  literal realization of the two-family lens. Decision — the boundary between the two families is that
  they consume *different sufficient statistics*: kernel-weighting (§4.A) reads counts + bin locations,
  weights by ``K_h``, and localizes/solves (a function of a query point, needs a bandwidth); power-sum
  quadrature (§4.B) contracts the within-bin power-sums ``S_p`` against a test function's local coefficients
  (a scalar, no bandwidth). The summary stores both counts and power-sums precisely because there are
  two readout modes. This settles the tempting over-unification: the kernel estimators are *not* trivial
  cases of ``\mathbb{E}[g(X)]`` (they are ratios / local solves of *joint* moments; only the bare degree-0
  fixed-``h`` KDE is formally ``\mathbb{E}[K_h(X - x)]``), so a single universal map would collapse to
  the vacuous "every estimator is a functional of the empirical measure" — the mirror of the corrected
  "the jet is *the* primitive" lens. Two atoms, each estimator a readout of one. `g`'s coefficients come
  from a dependency-free local polynomial interpolant; the pushforward is `expectation` of ``K_h\circ g``.

- **ConditionalExpectations §4.2 — general ``\mathbb{E}[g(Y)\mid X = x]`` + conditional pushforward, built.**
  The general `conditional_expectation` (with `conditional_mean` folded in as ``g = \mathrm{identity}``, and
  `conditional_pushforward_density` (``f_{g(Y)\mid X}``)) makes the conditional side symmetric with the
  marginal `Expectations`. It refines the two-family boundary: conditional functionals *compose* the atoms
  — a per-``x``-column power-sum quadrature in ``y`` (§4.B) fed as the response to a kernel-weighted
  regression in ``x`` (§4.A). So "each estimator is a readout of *one* atom" (the §4.B entry above) softens
  to "readouts are built from the atoms, sometimes composing them" — `conditional_mean` always did (its
  trivial ``\sum Y`` response is the ``y``-quadrature of the identity); the general ``g`` just makes the
  inner ``y``-readout non-trivial. The order-2 ``y``-power-sums cap a general ``g`` at ``O(\Delta_y^2)``.
- **The archive leaves the tree** (2026-10-09). `archive/` held only the archived §4.4 module, its test,
  and a README with the restore steps — working code waiting for a consumer, not part of the package. The
  shared tooling's `RULE.GIT.TRACKED` admits only the package's own files, so the folder is removed from
  git rather than kept under a switch: git's history keeps it exactly at commit `37006b5`, restored with
  `git checkout 37006b5 -- archive/` and the README's steps.

- **The argument audit** (2026-10-15). Every exported method's arguments and defaults reviewed under the
  rules the owner adopted from the sibling packages (a keyword shows its real default; `nothing` only for
  a genuine absence; one name per concept and one concept per name; no placeholder replaced on entry).
  Done, each its own commit: the histogram constructors' `accumulator_type` / `count_type`; `density`'s
  `bounds` → `support` (the start's and the kernel's name for the same tuple); the asymmetric kernel's
  `local_scale_weight` (the cone's `weight` is a different quantity) and `pilot_bandwidth`; the
  conditional estimators' `degree` defaulting to `0` under an asymmetric kernel, where `1` was refused;
  and `density`'s `h = nothing`, replaced by its measured equal `bandwidth(summary; kernel, degree)`,
  which moved the package's argument order to kernel → degree → bandwidth everywhere (a default sees
  only what precedes it). Not done, with reasons: the local plug-in's `pilot_bandwidth = nothing` stays a
  computed default (a public pilot function for one default is not worth its surface; to reconsider);
  `h` keeps carrying the asymmetric kernels' scale ``b`` (one smoothing parameter per kernel, the kernel
  owning its parametrization, as `kde(…; bandwidth, kernel)` does elsewhere); and the `density` variants
  stay keywords rather than one estimator value — the axes are independent corrections, and the refused
  pairs are redundancies or undeveloped estimators, now tabulated in §4.1, not a reason to redesign.

Open decisions still flagged [OPEN] in §10 (the projection road, `HistogramND`, vector-valued
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
- M. P. Wand & M. C. Jones (1995), [*Kernel Smoothing*](https://doi.org/10.1201/b14876), Chapman & Hall.
  (The two-term finite-sample KDE variance ``R(K)\hat f_X/(nh) - \hat f_X^2/n`` behind the density,
  conditional-density, and pushforward standard errors.)
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
- M. Rosenblatt (1969), [*Conditional probability density and regression
  estimators*](https://openlibrary.org/works/OL7299345W) (no DOI), in *Multivariate Analysis II*,
  Academic Press; R. J. Hyndman, D. M.
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
  correction — `density(…; support)`.)
- S. X. Chen (2000), [*Probability density function estimation using gamma
  kernels*](https://doi.org/10.1023/A:1004165218295), Ann. Inst. Statist. Math. 52(3). (Gamma asymmetric
  kernel for ``[0,\infty)`` — `density(…; kernel=gamma_kernel)`.)
- S. X. Chen (1999), [*Beta kernel estimators for density
  functions*](https://doi.org/10.1016/S0167-9473(99)00010-9), Comput. Statist. Data Anal. 31(2). (Beta
  asymmetric kernel for ``[0,1]`` — `density(…; kernel=beta_kernel)`.)
- N. L. Hjort & I. K. Glad (1995), [*Nonparametric density estimation with a parametric
  start*](https://doi.org/10.1214/aos/1176324627), Ann. Statist. 23(3). (Semiparametric density — the
  parametric model times a KDE of its correction; `density(…; start)`.)
- N. L. Hjort & M. C. Jones (1996), [*Locally parametric nonparametric density
  estimation*](https://doi.org/10.1214/aos/1032298288), Ann. Statist. 24(4). (The `start` × local-likelihood
  combination the package has not built — §4.1's table.)
- A. Rydin Gorjão et al. (2021), [*Arbitrary-Order Finite-Time Corrections for the Kramers–Moyal
  Operator*](https://doi.org/10.3390/e23050517), Entropy 23(5) — the all-orders Bell-polynomial
  ``\tau``-correction; the [`kramersmoyal`](https://doi.org/10.21105/joss.01693) software lineage.
- A. Frishman & P. Ronceray (2020), [*Learning force fields from stochastic
  trajectories*](https://doi.org/10.1103/PhysRevX.10.021009), Phys. Rev. X 10; D. B. Brückner,
  P. Ronceray & C. P. Broedersz (2020), [*Inferring the dynamics of underdamped stochastic
  systems*](https://doi.org/10.1103/PhysRevLett.125.058103), Phys. Rev. Lett. 125. (Projection road;
  ULI on samples paired with their time derivatives.)
- L. van Hove (1954), [*Correlations in space and time …*](https://doi.org/10.1103/PhysRev.95.249),
  Phys. Rev. 95. (The MD two-time conditional density — the cross-domain anchor.)
- C. Beck & E. G. D. Cohen (2003),
  [*Superstatistics*](https://doi.org/10.1016/S0378-4371%2803%2900019-0), Physica A 322.
  (Hyperstatistics / superstatistics — the parameter-mixture pooling of §10.)

> Every reference above is a published, reader-accessible source — per PF's policy, *never cite an
> inaccessible source*. The private design notes and literature survey behind this architecture are
> deliberately not cited by the package itself.
