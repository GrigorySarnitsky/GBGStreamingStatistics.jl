# CLAUDE.md — GBGStreamingStatistics

Julia package: a backend-agnostic library for **online nonparametric statistics**. It estimates
marginal densities, conditional expectations, conditional and two-time transition densities, and
standard functionals of these (Kramers–Moyal coefficients, drift/diffusion, Pawula diagnostics)
from a **stream of samples that is never stored** — each sample is observed once, folded into a
fixed-size summary, and discarded. **Single node, shared-memory only** (CPU multithreading and a
single GPU; distributed/multi-GPU out of scope). Development is **GPU-first** but the code is
**backend-agnostic** (CPU `Array`, CUDA, ROCm) through KernelAbstractions.

The architectural spine is **empirical measure → mergeable summary → post-hoc-smoothed estimator**:
a fine histogram with per-bin moments is accumulated online and is *estimand-agnostic*; every
reported quantity is a post-hoc smoothing functional of that one summary. This is what delivers the
driving requirement — *no storage of the underlying state, and no recompute-from-raw reruns*.

The relationship to the sibling `GBGPeriodicFields` ("PF") is exact: **PF is all mathematics, no
physics. This package is all statistics, no physics.** A reference clone of PF lives under the
gitignored `inspiration/` (alongside the design notes); this package mirrors its layout and
conventions.

Before any structural or design change, read `docs/ARCHITECTURE.md`. It explains why the
architecture is the way it is; do not re-litigate recorded decisions without the user. **Design
changes must be justified by established estimation theory** — not by convenience. "It's
convenient", "for symmetry", or "for completeness" is never a valid justification (cf. the
deliberate absences of the invariants). The conceptual model and the writing/docstring conventions
are recorded once — in `docs/ARCHITECTURE.md` and this file; the umbrella module
`src/GBGStreamingStatistics.jl` only points there.

> **Status: lower half complete (CPU + GPU); estimators landing.** The summaries `Histogram1D` and
> `Histogram2D` (`observe!` / `observe_pairs!` as one KernelAbstractions + Atomix scatter-add of counts
> and power-sums about the bin center, `merge!` by addition — and `scale!` / weighted `merge!` the cone
> for hyperstatistics / forgetting, allocation-free), the marginal density
> `density` (§4.1), local-linear `conditional_mean` / `conditional_variance` (§4.2), the conditional /
> transition density `conditional_density` / `transition_density` (§4.3), the Kramers–Moyal `drift` /
> `diffusion` (§4.4), and global `bandwidth` selection (§4.5) are implemented and tested — Gaussian,
> bivariate-Gaussian, and **Ornstein–Uhlenbeck** oracles, mergeability, allocation, and **GPU-vs-CPU
> parity verified on the AMD RX 6700 XT** (ROCm). The response-aware local bandwidth and the Pawula /
> Markov diagnostics are declared and forthcoming. See *Current state and roadmap*.

## Modules and wiring

Two halves (`docs/ARCHITECTURE.md` §3, §4, §8). The **universal lower half** is pure numerics over a
stream and transfers verbatim to every consumer; the **estimator upper half** is the exportable
estimands.

- `src/Accumulators.jl` — the foundation (lower half): the `AbstractAccumulator` supertype and the
  `observe!` / `merge!` / `scale!` interface only — the mergeable monoid, extended to the **cone** by the
  weighted `merge!(into, from; weight)` and `scale!` (§10 hyperstatistics / forgetting). The concrete
  summaries and the portable KernelAbstractions + Atomix scatter kernels live in `Histograms`.
- `src/Histograms.jl` — the summaries (lower half): `Histogram1D` / `Histogram2D` — linear (WARPing)
  binned empirical measures over a fixed range, each bin/cell carrying its **count** and **power-sums
  about the bin center** ``S_p`` / ``S_{rs}`` (the local sufficient statistic for moments and
  local-polynomial fitting). Accumulated by one **KernelAbstractions + Atomix atomic scatter-add**
  (`observe!` / `observe_pairs!`, one work-item per sample, CPU and GPU alike); `merge!` is elementwise
  addition (and `scale!` / weighted `merge!` elementwise scaling — the cone, on the `count = Float64`
  Float-count variant); `bin_index` is the binning map (exposed for pre-quantization); global central
  moments are recovered post hoc by the **shift formula**. Owns `KernelAbstractions`, `Atomix`.
- `src/DensityEstimation.jl` — §4.1 `density` ``\equiv f_X`` — one call over the binned KDE, the
  local-log-polynomial (`degree`), the Abramson variable bandwidth (`adaptive`), boundary
  renormalization (`bounds`, Jones 1993), the **parametric start** (`start`, Hjort & Glad 1995 — the
  moment-fitted `normal_start` / `lognormal_start` / `gamma_start` / `exponential_start` /
  `inverse_gaussian_start`, or a user-supplied `ParametricStart`), exact `normalize`
  (``\int\hat f = 1``), and the **`binning`** location rule (`:linear` centroid default — ``O(\Delta^2)``
  via the stored ``S_1``; `:simple` bin-center, ``O(\Delta)``) — plus the symmetric `gaussian_kernel` / `epanechnikov_kernel` (or any callable
  with ``\int K = 1``), the **asymmetric** `gamma_kernel` / `beta_kernel` for bounded/positive support
  (boundary-bias-free; Chen 1999/2000, or a user-supplied `AsymmetricKernel`), and
  `density_standard_error` (the pointwise variability band
  ``\mathrm{SE}(x) = \sqrt{R(K)\hat f_h(x)/(nh)}``). `AsymmetricKernel` and `ParametricStart` are
  exported so consumers can plug in custom kernels / parametric skeletons.
- `src/Regression.jl` — §4.2 `conditional_mean` ``\equiv m(x)``, `conditional_variance`
  ``\equiv \sigma^2(x)`` (local-linear, Fan & Gijbels; Nadaraya–Watson kept for reference). The
  **estimand-first reference module** — write its docstrings to the documentation policy below.
- `src/ConditionalDensity.jl` — §4.3 `conditional_density` ``\equiv f_{Y\mid X}``,
  `transition_density` ``\equiv p_\tau(y\mid x)`` (double-kernel, Rosenblatt).
- `src/MarkovInference.jl` — §4.4 `increment_moment` ``\equiv M_n(x,\tau)``, `kramers_moyal`
  ``\equiv D^{(n)}``, `drift` ``\equiv D^{(1)}``, `diffusion` ``\equiv D^{(2)}``, `pawula_quotient`
  ``\equiv Q(x)``, `markov_test`. Model-neutral functionals; the SDE *reading* stays consumer-side.
- `src/Bandwidth.jl` — §4.5 `bandwidth`: the **kernel- and degree-aware normal-reference rule** (the
  default, Silverman the Gaussian/degree-0 case; ``h = C(K,d)\,\hat\sigma\,n^{-1/(2p+1)}``), the
  data-driven **Sheather–Jones** plug-ins (`rule = :sj` classic, `:isj` improved-via-DCT) and the fit-based **AIC** (`rule = :aic`; leave-one-out CV a recorded binned-ties gap), the **Chen plug-in** for the asymmetric gamma/beta kernels (`rule = :plugin`, the asymmetric Silverman analog), and the
  **response-aware local plug-in** ``h^\star(x)``; the post-hoc FFT smoothing helpers (the major
  performance lever, amortized over bandwidth). Owns `FFTW`, `AbstractFFTs`.
- Include order: `Accumulators → Histograms → DensityEstimation → Regression → ConditionalDensity →
  MarkovInference → Bandwidth`. Cross-module access only through explicit `using ..Module: name, …`
  — never a bare `using`.
- Deferred, recorded as commented-out includes (the rule against speculative structure):
  `# CharacteristicFunctions.jl` (the ECF summary — `docs/ARCHITECTURE.md` §2, likely never built),
  `# ProjectionInference.jl` (the SFI / ULI projection road — §4.6).

## Core invariants — do not violate

The PF discipline: a missing capability encodes a decision (`docs/ARCHITECTURE.md` §6).

1. **Accumulation is estimand-agnostic; smoothing is post hoc.** No summary stores a bandwidth,
   kernel, or channel choice. An API that bakes a bandwidth into accumulation is banned. *(This is
   the invariant that makes "no reruns" true.)*
2. **Samples in, statistics out — no domain, no model.** Domain-neutral vocabulary (sample /
   variate / response / lag / conditioning value), never particle / tracer / droplet / velocity. No
   stationarity, Markov, diffusion, or parametric assumption is a *premise*; such structure is only
   ever *tested* (the Pawula ratio, the Markov / Chapman–Kolmogorov test).
3. **Everything is a mergeable monoid — in fact a cone.** Every summary supports `merge!`; the
   unnormalized summary is additive over the stream. No non-mergeable accumulator is admitted. The
   unnormalized summary is also **scalable** (`scale!`), so it lives in a convex cone — equal-weight
   `merge!` is the built special case, the weighted `merge!(into, from; weight)` and `scale!` the cone
   operations (§10 hyperstatistics / exponential forgetting, on the `count = Float64` variant). Tested:
   accumulate-then-merge `=` accumulate-together, and the cone laws.
4. **The hot path is allocation-free and backend-agnostic.** `observe!(summary, batch)` allocates only
   a bounded amount (the kernel launch), never an ``O(n)``/``O(bins)`` temporary (enforced by an
   `@allocated` audit on the CPU). Storage is generic `AbstractArray`; no scalar indexing of device
   arrays; one portable KernelAbstractions + Atomix kernel does the atomic scatter-add binning on CPU
   and any GPU backend. Power-sums about the bin center make each update a plain atomic add, so global
   atomics handle modal-bin **contention** directly (privatization stays a future option).
5. **The summary is rich enough to recover every target post hoc.** A target that would need a finer
   or different summary is a *recorded coverage gap*, not silently approximated.
6. **Accumulate wider than the data.** Counts are integer (`UInt64` merged); moment accumulators
   default to `Float64` even on `Float32` input (a Float32 count saturates at ``2^{24}``; a Float32
   high-order moment loses the data's own precision through the reduction). The deliberate opposite
   of PF's preserve-the-eltype — this is a *reduction*, not a *transform*. Accumulator precision is a
   type parameter; the default is Float64.

## Style

Two-space indentation; two blank lines between top-level definitions; section banner comments;
descriptive names — no cryptic or domain abbreviations (spell out range bounds as `lower`/`upper`, not
`lo`/`hi`), though universal idioms are fine (`n`-prefixed counts like `nbins`/`nx`, standard
single-letter math symbols); moderately long lines are fine. When unsure about
Base/stdlib/package behavior, check official docs — don't guess. For a linear solve use
`\`/factorizations/`ldiv!`, never an explicit matrix inverse (the local-polynomial normal equations
are a 2×2/3×3 solve). Operators and extended Base/`Statistics` functions are imported before being
extended (`import Base: merge!`); no type piracy.

**Docstrings vs comments.** Prefer docstrings; reserve plain comments for non-intuitive *why*
(the atomic scatter-add, the shift-formula recovery, the power-sums-about-the-bin-center choice,
"wider-than-data because this is a reduction"). One docstring per operation, on the bare binding:
`function f end` for package-owned generics, the bare imported name for Base/`Statistics` functions.

**Docstrings (and docs) are a code-agreement guard.** Every claim must be true of the *present* code —
a reader (or `/check-docs`) checks the code against the docstring as a bounded consistency pass. A
capability that is merely *enabled* or *planned* must be marked (forthcoming / future / a recorded
gap), never stated as present; that is what keeps the guard sound.

## Documentation policy — estimand first, estimator second

This is the user's explicit requirement and the exact mirror of PF's *physical-meaning-first* rule
(`docs/ARCHITECTURE.md` §2, §7). Because **every operation estimates a population quantity from
data**:

> A docstring (and every doc-page entry) for an estimator must **state the estimand — the
> population quantity, in standard probability notation, with its defining formula — *before* the
> estimator** that approximates it from the accumulated summary.

`conditional_mean` is the intended reference: lead with ``m(x) = \mathbb{E}[Y \mid X = x] = \int y\,
f_{Y\mid X}(y\mid x)\,\mathrm{d}y`` (the ``L^2`` projection), *then* the local-linear estimator. State
the estimand whenever a clean population definition exists; for an intrinsically finite-sample object
(a bandwidth selector, `merge!`) say so plainly rather than forcing a population form.

- **Symbol bound to code on first mention:** `density` ``\equiv f_X``, `conditional_mean`
  ``\equiv m(x)``, `transition_density` ``\equiv p_\tau(y\mid x)``, `drift` ``\equiv D^{(1)}``,
  `diffusion` ``\equiv D^{(2)}``.
- **Math is rendered LaTeX in standard notation** — Documenter KaTeX: inline `` ``…`` `` (backslashes
  escaped `\\` inside a `"""` docstring), display ` ```math `; never Unicode-glyph approximations,
  never math in backticks (backticks are for code identifiers only). **And the converse — a code
  identifier never goes inside a math span**: it renders as a product of italic single letters
  (`xlo` → x·l·o); backtick it, or wrap it in `\texttt{}` if it must appear inside a formula. **A
  formula is a single LaTeX span — never ``math`` split onto an adjacent `code` span**: ``1 \le i \le``
  `nbins` is wrong, write ``1 \le i \le \texttt{nbins}``; the same name is backticked standalone in
  prose but `\texttt{}`-wrapped inside a formula. **And backticks hold only well-formed Julia** — a
  name, type, call, or valid expression — never math notation (an interval, a product ``a \times b``,
  an inequality, a `·`/`−`/`≤` glyph), which is LaTeX. Plain `#` comments stay Unicode.
- **Correctness is its own check, distinct from formatting.** Verify the estimator is **consistent
  for the stated estimand** under the stated assumptions, and that **every normalization is exact**:
  the ``1/n``, ``1/(nh)``, ``\int K = 1``, kernel moments ``\mu_2(K) = \int u^2 K``, roughness
  ``R(K) = \int K^2``, the bias/variance constants. A stray ``1/h`` or kernel-moment factor is the
  easy slip here (the statistics analog of PF's stray ``1/N`` / Parseval factor).
- **Citations:** context-standard author–year; name up to three authors, `First et al.` for four or
  more; **never cite an inaccessible source** (the private notes under `inspiration/` are working
  material, not citable). The canonical reference list lives in `docs/ARCHITECTURE.md` §References;
  docstrings carry the full reference inline, verbatim-consistent with it.

**Doc-edit checklist — run all on every documentation change** (most strictly on `ARCHITECTURE.md`):
**(1) Correct** — every claim/formula/normalization true and precise, the estimator consistent for
its estimand, run as a *positive* per-formula verification, not a read-through; **(2) Self-contained**
— every symbol *bound* to its meaning in plain language at first use (enumerate symbols, cite the
binding clause); **(3) Notation** — one symbol per concept and one concept per symbol; standard
notation; **(4) Formatted (bidirectional)** — math is LaTeX never backticks, every code identifier
backticked never bare, and inline ``math`` stays on **one source line** (a wrapped continuation that
begins with a list/blockquote leader — `- `/`* `/`> `/`1.` — is parsed as a list item and silently
drops the span out of math mode; `docs/make.jl` lints for this and fails the build, since Documenter
renders the broken page without complaint); a **formula is a single LaTeX span**, never ``math`` split
onto an adjacent `code` continuation (a code name inside a formula is `\texttt{}`-wrapped), and
**backticks hold only well-formed Julia** (a name/type/call/expression), never math notation —
intervals, products, inequalities, `·`/`−`/`≤` glyphs; **(5) Current & complete** — no stale *deferred/none-exist* claims, and a
summary that names only a subset of what the code now provides is stale. "Check the docs" means this
**complete battery across every surface** (the four rationale docs, every `src/` docstring, every
`docs/src/` page) — the `.claude/commands/check-docs.md` command (`/check-docs`) encodes it. A
flagged instance is a sample of a class — fix the class, not the line.

## Verification workflow

Before changing tests, read `docs/TESTING.md`; before optimizing/benchmarking, read
`docs/PERFORMANCE.md`. The nontrivial statements live formally in `docs/FORMALIZATION.md` — the
chain is *prose ↔ formalization ↔ test*.

### Definition of done — for a new or changed operation

Treat the feature as unfinished until each is **done and verified** (not "should work"):

1. **Correctness tests** against an *independent* analytic oracle — Gaussian marginals; an
   **Ornstein–Uhlenbeck process** (drift, diffusion, transition density known exactly) for §4.3–4.4
   — plus a **consistency / rate test** (error falls at the theoretical ``n^{-2/5}``), a
   **mergeability test** (accumulate-then-merge `=` accumulate-together), and — for an in-place op —
   an `@allocated` assertion in `test/allocations.jl`. **GPU-vs-CPU parity** under
   `allowscalar(false)` once a device path exists.
2. **The guiding docs** — a `docs/ARCHITECTURE.md` section (+ decision-log entry), a
   `docs/FORMALIZATION.md` typed entry (estimand / estimator / consistency claim with every constant
   pinned, test-linked), `docs/TESTING.md`, `docs/PERFORMANCE.md`, as applicable.
3. **The rest of the docs + this file** — the API page (`docs/src/api/…` and `make.jl`), the
   `docs/src/index.md` / `api.md` lists, and **`CLAUDE.md`** (module list, include order, roadmap).

Mechanical bar:

- Load: `julia --project -e 'using GBGStreamingStatistics'`.
- `julia --project -e 'using Pkg; Pkg.test()'` — keep it green. `test/aqua.jl` (`Aqua.test_all`:
  ambiguities, piracy, stale deps, compat bounds, undefined exports) and `test/jet.jl` (JET static
  analysis) run against the loadable package.
- Docs build: `julia --project=docs docs/make.jl` (`checkdocs = :exports`, doctests) — must be clean.
- Consult `Project.toml` for dependencies and versions; do not assume.

## Current state and roadmap

**Done (implemented + tested, CPU and GPU).**

1. **`Histograms.Histogram1D` + `Accumulators`** — linear (WARP) binning over `[lower, upper)` (out-of-range
   skipped); per-bin counts + **power-sums about the bin center** ``S_p`` to a configurable `order`.
   `observe!` (batch + scalar) and `merge!` (elementwise addition); `bin_index`; global central moments
   recovered by the **shift formula** and exposed publicly as `mean`/`var`/`std`/`central_moment`
   (extending `Statistics`). Tested:
   Gaussian moment recovery (exact vs direct, and vs the true moments), mergeability, bounded
   allocation, Float32→Float64.
2. **`Histograms.Histogram2D` + `observe_pairs!`** — the joint two-time summary: `nx×ny` cells, per-cell
   counts + order-2 bivariate power-sums ``S_{rs}`` about the cell center (the regression sufficient
   statistic). Tested: independent-Gaussian oracle (marginals + covariance), mergeability, allocation.
3. **The portable scatter kernel** — `observe!` / `observe_pairs!` are one KernelAbstractions + Atomix
   atomic scatter-add (one work-item per sample), allocation-free up to the launch, **verified identical
   on CPU and the AMD RX 6700 XT** (the gated ROCm suite, `test/gpu/`, under `allowscalar(false)`,
   including heavy atomic contention).
4. **`DensityEstimation.density`** (§4.1) — one call over the grid of **local-log-polynomial `degree`**
   (0 = KDE; ≥1 = binned local-Poisson, Loader 1996 — degree 1 exact for exponential tails, 2 for
   Gaussians) × **fixed / Abramson-variable bandwidth** (`adaptive`), with a **pluggable kernel**
   (`gaussian_kernel` default, `epanechnikov_kernel`). The estimand-first reference; tested across the
   grid (recovers the Gaussian pdf into the tail, beats fixed-``h`` there; degree-1 beats the KDE in a
   Laplace tail), and across a **representative analytic-PDF suite** (uniform, exponential, Student-``t``
   heavy tail, log-normal, bimodal) with a **convergence-rate test** — the ISE falls at
   ``n^{-2p/(2p+1)}``, degree-2's steeper ``n^{-8/9}`` vs degree-0/1's ``n^{-4/5}``, the Epanechnikov
   reproducing the kernel-independent ``n^{-4/5}``, plus a **binning-error convergence** test (the binned
   estimate → the exact KDE at ``O(\\Delta^2)`` for the default **linear/centroid `binning`** — mass at the
   recovered within-bin centroid ``c_j + S_1^{(j)}/n_j`` — and ``O(\\Delta)`` for `binning = :simple`,
   mass at the bin center) — in `test/density_convergence.jl`. Also the **bounded/positive-support** kernels (`gamma_kernel` `[0,∞)`,
   `beta_kernel` `[0,1]`; Chen's bias-reduced ``f_2``), the **parametric start** (`start`, Hjort & Glad
   1995 — five moment-fitted starts; a bias reduction, near-exact where the model fits, ``>3\times``
   better than the KDE on log-normal data), and exact **`normalize`** (``\int\hat f = 1``); tested in
   `test/densityestimation.jl` (oracle recovery, the start's wide-``h`` bias removal, normalization,
   mergeability, rejected knob combinations).
5. **`Bandwidth.bandwidth`** (§4.5) — the **kernel- and degree-aware normal-reference** rule
   ``h = C(K,d)\\,\\hat\\sigma\\,n^{-1/(2p+1)}`` from the recovered ``\\hat\\sigma``: the kernel enters
   through the canonical factor ``[R(K)/\\mu_2(K)^2]^{1/5}`` (by quadrature), the degree through the
   bias order ``p`` (``n^{-1/5}`` at degree 0/1, ``n^{-1/9}`` at degree 2); reduces to Silverman
   ``(4/3)^{1/5}`` (Gaussian, degree 0), ``\\approx 2.34`` (Epanechnikov). It is also `density`'s
   default `h`. Tested (`test/bandwidth.jl`): Silverman reproduction, the Epanechnikov ``2.214\\times``
   ratio, the degree-2 ``n^{-1/9}`` closed form. The data-driven **Sheather–Jones** plug-in
   (`rule = :sj`, degree 0/1) is also available — ``\\psi_4``/``\\psi_6`` as mergeable Hermite-weighted
   reductions of the bin counts (Sheather & Jones 1991), the solve-the-equation root-find, the eval
   kernel a clean ``\\delta(K)`` rescale; it adapts to shape (much tighter than the normal reference on
   a bimodal sample); its DCT-based **improved** variant (`rule = :isj`, Botev 2010) reads the curvature
   off the bin counts, robust on strongly multimodal data. Tested (both): ``\\approx`` normal on Gaussian, the exact ``\\delta(K)`` ratio,
   mergeability, the degree ``\\ge 2`` gap. The fit-based **AIC** (`rule = :aic`, degree 0/1) over an
   ``h``-grid is also available — a mergeable information criterion of the binned KDE with a stable
   interior optimum; the leave-one-out CV selectors (LCV, LSCV) are a recorded coverage gap (binned ties
   drive them to zero bandwidth). For the **asymmetric gamma/beta kernels**, `bandwidth` adds the
   ``\\hat\\sigma\\,n^{-2/5}`` rule of thumb (`rule = :normal`) and the **Chen plug-in** (`rule = :plugin`,
   the asymmetric Silverman analog ``b^\\star = [\\mathbb E[w(X)^{-1/2}]/(2\\sqrt\\pi\\,n\\int w^2 f''^2)]^{2/5}``,
   weight ``w = x`` / ``x(1-x)``, curvature from a grid-quadrature pilot ``\\hat f''`` — boundary-valid,
   unlike the ``\\psi_4`` U-statistic); ``\\approx 1.4\\text{–}1.7\\times`` lower ISE on skewed positive data.
   Tested in `test/bandwidth.jl`; the AIC, U-statistic, and gamma-reference routes were tried and rejected
   (see the §11 decision log), and a fully-adaptive curvature pilot is a recorded refinement.
6. **`Regression`** (§4.2) — `conditional_mean` ``\equiv m(x)`` and `conditional_variance``
   ``\equiv \sigma^2(x)`` by local-linear regression from `Histogram2D` (the ``T_\ell``/``P_\ell``
   recovered exactly from per-X-column counts and response totals). Tested against the bivariate-
   Gaussian oracle (``m`` recovers ``\rho x`` to ``\lesssim 2\%``; ``\sigma^2`` recovers ``1-\rho^2``
   to a few percent — the simple ``\hat{\mathbb E}[Y^2\mid x]-\hat m^2`` estimator carries the
   local-linear smoothing bias, Fan–Yao refinement forthcoming).
7. **`MarkovInference`** (§4.4 — drift & diffusion) — `increment_moment` ``\equiv M_n(x,\tau)``,
   `kramers_moyal` ``\equiv D^{(n)}``, `drift` ``\equiv D^{(1)}``, `diffusion` ``\equiv D^{(2)}`` for
   ``n \in \{1,2\}`` from the order-2 `Histogram2D` (``M_1 = \hat m_h(x) - x``,
   ``M_2 = \hat{\mathbb E}[Y^2\mid x] - 2x\hat m_h + x^2``), reusing the §4.2 local-linear primitives;
   reported as the finite-``\tau`` ratio ``M_n/(n!\tau)``. Tested against the **Ornstein–Uhlenbeck**
   oracle (exact stationary lagged pairs): `drift` recovers ``x(e^{-\theta\tau}-1)/\tau`` exactly
   (linear OU mean ⇒ no smoothing bias); `diffusion` recovers the closed-form ``D^{(2)}`` to ``8\%``
   (the curved-``\mathbb E[Y^2\mid x]`` local-linear bias); plus the variance decomposition, the
   ``1/(n!\tau)`` scaling, coverage-gap errors (``n\ge3``, static lag), and mergeability.
8. **`ConditionalDensity`** (§4.3) — `conditional_density` ``\equiv f_{Y\mid X}`` and `transition_density`
   ``\equiv p_\tau(y\mid x)`` by the binned double-kernel (Rosenblatt) from the `Histogram2D` cell
   counts: the joint is reduced at ``x`` to an X-kernel-weighted Y-profile, then smoothed in ``y`` with
   the surviving ``1/(\sqrt{2\pi}\,h_2)`` normalization (``\int \hat f\,\mathrm dy = 1`` exactly). Two
   supplied bandwidths `hx`, `hy`. Tested: the bivariate-Gaussian conditional ``N(\rho x, 1-\rho^2)``
   (normalization + pointwise) and the **Ornstein–Uhlenbeck** transition density
   ``N(x e^{-\theta\tau}, s^2(1-e^{-2\theta\tau}))``; scalar/vector forms, `transition_density`
   equivalence, off-support NaN, and mergeability.
9. **The cone — `scale!` + weighted `merge!`** (§10 hyperstatistics / forgetting) — the unnormalized
   summaries are scalable, so they form a convex cone, not merely an additive monoid: `scale!(summary, w)`
   reweights every count and power-sum (the scalar action), and `merge!(into, from; weight)` forms
   ``\text{into} \mathrel{+}= \text{weight}\cdot\text{from}`` (the default `weight = 1` the plain monoid
   add). On a **Float-count variant** (`count = Float64`; the default `UInt64` cannot hold a fractional
   ``w\,n``) for both histograms. Realizes the **hyperstatistical** parameter-mixture density
   ``P(x) = \sum_k w_k f_k(x) = \sum_k (w_k/N_k) T_k`` (superstatistics; Beck & Cohen 2003) — keep one
   summary per ``\beta``, choose the mixing law ``f(\beta)`` *post hoc*, re-weight without re-accumulating —
   and **exponential forgetting** (`scale!(s, λ)` before the next block). Tested (`test/hyperstatistics.jl`):
   the cone laws, the discrete Gaussian-mixture oracle, the **χ²-superstatistics → Student-``t``** Beck–Cohen
   oracle (a Gamma-weighted Gaussian superposition recovers the heavy ``t`` tail), exact exponential
   forgetting, the `Histogram2D` cone with `conditional_mean` scale-invariance, allocation, and GPU parity
   (incl. the atomic-**float** add of Float-count `observe!`).

**Next, in dependency order:**

1. **The response-aware local bandwidth** ``h^\\star(x)`` (§4.5) — estimand-first, OU-oracle-tested,
   with its `FORMALIZATION.md` entry, from the `Histogram2D` joint summary. Plus the Fan–Yao
   conditional-variance refinement and the FFT-amortized smoothing (§4.5).
2. **Scheduled — the conditional convergence-rate tests** (the one characteristic not yet uniform across
   the testing contract; `docs/TESTING.md`, the coverage matrix). `conditional_mean` at the ``n^{-2/5}``
   ``L^2``-projection rate and `conditional_density` at the ``n^{-2/3}`` rate, each mirroring the
   implemented `density` ISE-slope test (`test/density_convergence.jl`) — accumulate the bivariate-Gaussian
   / OU oracle at growing ``n``, fit the log–log error slope, check it against the exponent. `drift` /
   `diffusion` inherit `conditional_mean`'s rate (they are linear functionals of it); the `bandwidth`
   selectors and `density_standard_error` need none (a selector and a variability band, not function
   estimators). This closes the rate column of the matrix.
3. **The rest of `MarkovInference` §4.4** — `kramers_moyal` for ``n \ge 3``, the **Pawula `Q`**
   diagnostic (``M_4``, ``M_6``) and the all-orders Bell-polynomial finite-time correction, and
   `markov_test` (now unblocked — it composes the §4.3 transition density). The higher moments need the
   `Histogram2D` order-raise (Deferred, below).
4. **Density-estimation methods from the prior-art survey** (`inspiration/density-software/ALGORITHMS.md`),
   in priority order — all mergeable functionals of the summary: the **POT–GPD heavy-tail splice**
   (domain-critical for intermittent tails; a binned GPD-MLE on the tail bin counts that extrapolates
   beyond the histogram's `[lower, upper)` range, plus an optional mergeable top-k max-heap for exact
   Hill/Pickands tail indices); the **P-spline penalized-Poisson** density (Eilers–Marx, counts-native
   global smoothing); **max-entropy from the moments** (uses the power-sums); the **self-consistent**
   Fourier-filter KDE (Bernacchia–Pigolotti) and **orthogonal series** (Hermite/DCT, seeding §4.6). The
   asymmetric gamma/beta kernels are **done** (Done 4). Deliberate skips recorded there: shape-constrained
   MLE (shape-as-premise, invariant 2), taut-string (a mode-count diagnostic), DEFT, kNN, neural/DP-mixture.

Deferred (recorded, not built): the ECF summary (`docs/ARCHITECTURE.md` §2 — likely never); the
projection / SFI–ULI road (§4.6); `HistogramND` (§3, wait for a multi-conditioning consumer);
vector-valued conditioning (lean scalar-channel — let the consumer loop components); the **bulk scalar time-correlation function**
(consumer-side — a `Histogram2D` per lag is the wrong cost for one cross-moment; §6); raising
`Histogram2D` beyond its order-2 per-cell power-sums (the higher Kramers–Moyal / Pawula moments (``M_4``–``M_6``) of *Next* 3 will be
its first consumer); the **leave-one-out cross-validation** bandwidth selectors
(likelihood-CV, least-squares LSCV — a binned summary's ties drive their objective to zero bandwidth; the
built `:aic` / `:sj` are the data-driven selectors instead, and a within-bin-spread correction from the
power-sums could rescue CV later).

## Working agreements

- The API is not frozen; do not build backward-compatibility shims.
- **Commit ⇒ push, on `master`.** This is a *private* repo synced across machines through GitHub, so
  work is committed **directly on `master`** (do not branch) and **pushed in the same step**; the
  docs site deploys on push. End every commit message with `Co-Authored-By: Claude Opus 4.8
  <noreply@anthropic.com>`; never use `[skip ci]`.
- **CI tests do not run on every push — manual + release tags** (`CI.yml` is `workflow_dispatch`
  plus `push: tags: '*'`): per-commit GitHub Actions test runs are kept off to avoid CI cost, so run
  the full test matrix **locally** before every push. The matrix runs on GitHub only when triggered
  by hand (Actions tab) or on a **release tag** — the from-scratch, clean-environment resolve a local
  run cannot give. This is the shared policy of the sibling GBG packages (GBGIsotropicTurbulence,
  GBGCloudPhysics). **The docs workflow builds + deploys on every push to `master`.**
