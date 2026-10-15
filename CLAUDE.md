<!-- gbg:managed begin v2.0.0-rc.69 — rewritten by `gbg sync`; do not edit -->
@.gbg/policy/00-index.md

The rulebook above is the whole of it: the rules, the terms they use, and what GBGStreamingStatistics has
switched off with its reason. Nothing restates a rule outside it — a rule quoted in two files is a
rule that will be reworded in one.

Unasked, run `Pkg.test()` and nothing else; before every push the docs build runs too, through the
push hook, and a failing build stops the push. What a rule finds is reported, never fixed silently:
documentation may be fixed without asking, and anything else waits for an answer.

The readers of a review are derived from the rules and the tool's catalogue of disciplines. A
reader or a discipline this package needs and the tool lacks is requested of the tool, as a request
file the owner has approved, never written here.

A run an earlier session left open is continued before anything else, and the chat says so:
`gbg run next` says what to launch or that the stage is ready to collect, and the run goes on to
its close. It waits for the owner only when `gbg run next` says the tree moved since the launch.
<!-- gbg:managed end -->

# CLAUDE.md — GBGStreamingStatistics

Julia package: a backend-agnostic library for **online nonparametric statistics**. It estimates
marginal densities, conditional expectations, conditional and two-time transition densities, and
functionals of these (expectations, moments, and pushforward densities)
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
deliberate absences of the invariants). The conceptual model is recorded in `docs/ARCHITECTURE.md`; the
writing and docstring conventions are the rulebook's (`.gbg/policy/`, imported at the head of this file)
plus the package's own, below. The umbrella module `src/GBGStreamingStatistics.jl` only points there.

> **Status: lower half complete (CPU + GPU); the two readout families landed.** The summaries `Histogram1D`
> and `Histogram2D` (`observe!` / `observe_pairs!` as one KernelAbstractions + Atomix scatter-add of counts
> and power-sums about the bin center, `merge!` by addition — and `scale!` / weighted `merge!` the cone
> for hyperstatistics / forgetting, allocation-free), the marginal density `density` (§4.1), the shared
> local-polynomial jet primitive `LocalPolynomial.local_fit` and the local-polynomial `conditional_mean` /
> Fan–Yao `conditional_variance` / general `conditional_expectation` (§4.2), the conditional / transition
> density `conditional_density` / `transition_density` (§4.3), the **power-sum quadrature** readouts
> `expectation` / `pushforward_density` and the recovered moments (§4.B), and global **and response-aware
> local** (``h^\star(x)``) `bandwidth` selection (§4.5) — each with its **standard error** and an
> estimate-plus-error paired reader — are implemented and tested (Gaussian, bivariate-Gaussian, and
> **Ornstein–Uhlenbeck** oracles, mergeability, allocation, and **GPU-vs-CPU parity verified on the AMD RX
> 6700 XT**, ROCm). The Kramers–Moyal `drift` / `diffusion`, Pawula, and Markov diagnostics (§4.4) are
> **archived** (retrievable) / forthcoming. See *Current state and roadmap*.

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
  addition (and `scale!` / weighted `merge!` elementwise scaling — the cone, on the `count_type = Float64`
  Float-count variant); `bin_index` (the binning map) is **archived**, the pre-quantization / index-ingestion path deferred. The **pure summary**
  (accumulate / `merge!` / `scale!` / bin) — the shift-formula moment recovery (`mean` / `var` / `std` /
  `central_moment`) and the general quadrature moved to `Expectations` (§4.B); the power-sums it stores are
  their sufficient statistic. Owns `KernelAbstractions`, `Atomix`.
- `src/Kernels.jl` — **all the smoothing kernels** and the kernel-moment quadrature (lower half): the
  **symmetric** `gaussian_kernel` / `epanechnikov_kernel` (any callable with ``\int K = 1``) with the
  internal `_quad` (``R(K)``, ``\mu_2(K)``, … by quadrature, so any kernel works), and the **asymmetric**
  varying-shape boundary kernels `gamma_kernel` (``[0,\infty)``) / `beta_kernel` (``[0,1]``) plus the
  `AsymmetricKernel` type for bounded/positive support (boundary-bias-free; Chen 1999/2000, or a
  user-supplied one). The cross-cutting weighting primitives — shared by `density`, the local-polynomial
  fits, the double-kernel `conditional_density`, and the `bandwidth` rules. Owns `SpecialFunctions`.
- `src/LocalPolynomial.jl` — the shared **local-polynomial jet primitive** (lower half): `local_fit`, the
  kernel-weighted local-polynomial design solve over a binned measure, returning the jet
  ``(\beta_0, \dots, \beta_d)`` (value ``\beta_0``, slope ``\beta_1``, curvature ``2\beta_2``). The single
  *linear* core on which the local smoothers are built — `conditional_mean` is its ``\beta_0``, the
  local-likelihood `density` its IRLS inner solve (`docs/ARCHITECTURE.md` §11). Internal (not exported);
  reached by explicit `using ..LocalPolynomial: local_fit`. Owns `LinearAlgebra`.
- `src/MarginalDensity.jl` — §4.1 `density` ``\equiv f_X`` — one call over the binned KDE, the
  local-log-polynomial (`degree`), the Abramson variable bandwidth (`adaptive`), boundary
  renormalization (`support`, Jones 1993), the **parametric start** (`start`, Hjort & Glad 1995 — the
  moment-fitted `normal_start` / `lognormal_start` / `gamma_start` / `exponential_start` /
  `inverse_gaussian_start`, or a user-supplied `ParametricStart`), exact `normalize`
  (``\int\hat f = 1``), and the **`binning`** location rule (`:linear` centroid default — ``O(\Delta^2)``
  via the stored ``S_1``; `:simple` bin-center, ``O(\Delta)``) — accepting any `kernel` from `Kernels`
  (the symmetric `gaussian_kernel` / `epanechnikov_kernel`, or the asymmetric `gamma_kernel` /
  `beta_kernel` for bounded/positive support), and `density_standard_error` (the pointwise variability
  band ``\mathrm{SE}(x) = \sqrt{R(K)\hat f_h(x)/(nh) - \hat f_h(x)^2/n}``, the two-term KDE variance — it
  **tracks the estimator's options**: the local Abramson ``h(x)`` for `adaptive`, the boundary-restricted
  ``R_a(x)/a(x)`` for `support`, `start` leading-order-neutral; degree-0/1 symmetric-kernel, the asymmetric
  variance a recorded gap), plus `density_with_standard_error` returning the estimate and its band from one
  option specification (a `(value, standard_error)` NamedTuple, drift-proof). `ParametricStart` is exported so consumers can
  plug in custom parametric skeletons (the `AsymmetricKernel` constructor for custom boundary kernels
  lives in `Kernels`).
- `src/ConditionalExpectations.jl` — §4.2 the conditional expectation `conditional_expectation`
  ``\equiv \mathbb{E}[g(Y)\mid x]`` of a smooth response functional (the general reader), with
  `conditional_mean` ``\equiv m(x)`` and `conditional_variance` ``\equiv \sigma^2(x)`` as its
  ``g = \text{identity}`` / ``y^2`` cases, and the conditional pushforward density
  `conditional_pushforward_density` ``\equiv f_{g(Y)\mid X}`` (the density of a transformed response, with
  its `conditional_pushforward_density_standard_error`), plus
  the sampling `conditional_expectation_standard_error` (``= \sqrt{\ell(x)\,\mathrm{Var}(g(Y)\mid x)}``, the
  local-regression leverage ``\ell(x)`` from `LocalPolynomial._leverage` — the reference SE mechanism for
  the kernel family), and the paired `conditional_expectation_with_standard_error` (the
  `(value, standard_error)` NamedTuple from one fit). All
  by local-polynomial regression **in x** (Fan & Gijbels, over `LocalPolynomial.local_fit`) of a
  per-X-column **y-quadrature** response (the `Expectations` coefficient primitive over the order-2
  Y-power-sums) — **composing both atoms** (kernel-weighting in x ∘ quadrature in y); with the `degree`
  (0 Nadaraya–Watson / 1 local-linear default / 2 local-quadratic) and `kernel` options. `conditional_mean`
  is folded in exactly as `conditional_expectation(identity, …)`; `conditional_variance` defaults to the
  bias-reduced **Fan & Yao** residual estimator (`method = :residual`; the plug-in is `:plugin`, ``=
  \mathbb{E}[Y^2\mid x]-\mathbb{E}[Y\mid x]^2`` via the general reader), with the residual-based
  `conditional_variance_standard_error` ``= \sqrt{\ell(x)\,(\mu_4-\sigma^4)}`` (``\mu_4`` the order-2
  quadrature of the smooth ``(Y-\hat m)^4``). The general `g` is `O(Δ_y^2)`
  (order-2 Y-power-sums), exact for `g` of degree ≤ 2. The **estimand-first reference module** — write its
  docstrings to the documentation policy below.
- `src/ConditionalDensity.jl` — §4.3 `conditional_density` ``\equiv f_{Y\mid X}``,
  `transition_density` ``\equiv p_\tau(y\mid x)`` (binned double-kernel, Rosenblatt) — a separate
  `kernelx` / `kernely`, with an **asymmetric** `kernely` (`gamma_kernel` / `beta_kernel`) for a
  bounded/positive response and `normalize` for exact ``\int\hat f\,\mathrm dy = 1`` — each with its
  pointwise `conditional_density_standard_error` (the two-term double-kernel band
  ``\sqrt{R(K_1)R(K_2)\hat f/(h_1 h_2 D) - \hat f^2/n}``, ``D = n\hat f_X(x)`` the estimator's own
  denominator; symmetric-kernel, the asymmetric-`kernely` variance a recorded gap) and the paired
  `conditional_density_with_standard_error`, `transition_density` aliased.
- `src/Bandwidth.jl` — §4.5 `bandwidth`: for the marginal (`Histogram1D`) the **kernel- and degree-aware
  normal-reference rule** (the default, Silverman the Gaussian/degree-0 case;
  ``h = C(K,d)\,\hat\sigma\,n^{-1/(2p+1)}``), the data-driven **Sheather–Jones** plug-ins (`rule = :sj`
  classic, `:isj` improved-via-DCT) and the fit-based **AIC** (`rule = :aic`; leave-one-out CV a recorded
  binned-ties gap), and the **Chen plug-in** for the asymmetric gamma/beta kernels (`rule = :plugin`, the
  asymmetric Silverman analog); for the joint (`Histogram2D`) the **response-aware local plug-in**
  ``h^\star(x)`` (`rule = :local`, composing the Fan–Yao ``\sigma^2(x)``, the degree-2 curvature
  ``m''(x)``, and ``f_X(x)``). The post-hoc FFT smoothing helpers (the major performance lever, amortized
  over bandwidth) are **deferred** (recorded, not built — a pure optimization over the exact direct
  summation). Owns `FFTW` (the ISJ `rule = :isj` DCT).
- `src/Expectations.jl` — §4.B the **power-sum quadrature** readouts (the second of the two readout
  families), and the recovered marginal moments they generalize: `mean` / `var` / `std` / `central_moment`
  (the ``g = x^p`` case, **moved here from `Histograms`**), `expectation` ``\equiv \mathbb{E}[g(X)]`` (a
  general smooth functional of `X`), its sampling `expectation_standard_error` (with the paired
  `expectation_with_standard_error`), and `pushforward_density`
  ``\equiv f_{g(X)}`` (the density of a transformed variable ``g(X)``, with its
  `pushforward_density_standard_error`). All read off the stored power-sums
  through the one shared contraction atom `_quadrature` (**owned here**) — the moments supplying closed-form
  binomial coefficients, the general `expectation` a dependency-free local-polynomial interpolant of `g`
  (exact on monomials); the pushforward is `expectation` of the composed kernel ``K_h \circ g`` (no
  ``g^{-1}``, non-monotone `g` handled). `MarginalDensity` / `Bandwidth` reach its `_sample_std` etc. by
  explicit `using ..Expectations`; reuses `LinearAlgebra`.
- Include order: `Accumulators → Histograms → Kernels → Expectations → LocalPolynomial → MarginalDensity →
  ConditionalExpectations → ConditionalDensity → Bandwidth`. Cross-module access only through explicit `using ..Module: name, …`
  — never a bare `using`.
- `ext/GBGStreamingStatisticsHDF5Ext.jl` — **HDF5 persistence** (a weak-dependency package extension, so
  `HDF5` stays out of the core deps): `write_histogram` / `read_histogram` (+ the `parent[name] = summary`
  sugar) map a summary to a self-describing HDF5 **group** (arrays → datasets, scalars → attributes, a
  type/version tag). The core declares the bare generics in the umbrella; the methods (on `HDF5.File` /
  `HDF5.Group`) live in the extension. Handle-based only — the consumer opens and owns the file.
- Deferred / archived, recorded as commented-out includes (the rule against speculative structure):
  `# IncrementMoments.jl` (the §4.4 SDE-reading `drift` / `diffusion` / Kramers–Moyal functionals —
  **archived** — kept in git history at commit `37006b5` (`archive/`), working and OU-tested; restore with
  `git checkout 37006b5 -- archive/` and its README when a consumer needs them and `Histogram2D`
  raises its order for ``n \ge 3``),
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
   ever *tested*, never assumed (e.g. by a Pawula ratio or Markov / Chapman–Kolmogorov test — the archived §4.4 diagnostics).
3. **Everything is a mergeable monoid — in fact a cone.** Every summary supports `merge!`; the
   unnormalized summary is additive over the stream. No non-mergeable accumulator is admitted. The
   unnormalized summary is also **scalable** (`scale!`), so it lives in a convex cone — equal-weight
   `merge!` is the built special case, the weighted `merge!(into, from; weight)` and `scale!` the cone
   operations (§10 hyperstatistics / exponential forgetting, on the `count_type = Float64` variant). Tested:
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

**Argument order — one convention, everywhere (internal and exported).** Every function orders its
arguments **object → query point(s) → smoothing controls → modifiers**: (1) the data operated on (the
`summary`, or the binned arrays — `counts` / `centers` / `locations` / `mass` / `response` / `P` / `D` /
`n` — in their natural order); (2) the **query point(s)** `x` (then `y` / the response point) where the
estimand is evaluated; (3) the smoothing controls — a **primary algorithm-selector** first where the function has one (`rule`,
which picks the whole estimator), then the fixed sub-order **kernel → degree → bandwidth** (`kernel` /
`kernelx` / `kernely` / `kernelz`, then `degree`, then `h` / `hx` / `hy` / `hz` / `pilot_bandwidth`): the
bandwidth comes last because it is a function of the kernel and the degree — `density`'s default is
`h = bandwidth(summary; kernel, degree)`, and Julia evaluates keyword defaults left to right, so the
inputs of a default precede it; (4) the **secondary** modifiers and refinements (`method`, `adaptive`,
`support`, `start`, `normalize`, `binning`, `maxorder`, IRLS scratch, …) in a stable order — a
`method`-style refinement (it tweaks *how*, not *which* estimator) trails, unlike a primary `rule`-style
selector, which leads. For an exported function, (1)–(2) are positional and (3)–(4) keyword — which is the
rulebook's `RULE.API.POSITIONAL` (GBGCodingTools rc.69): a required argument is positional or keyword as the
package chooses, and every argument with a default is a keyword unless its type selects the method. The
same order holds for the internal helpers' positional lists (`local_fit(locations, mass, response, x,
kernel, degree, h)`), so there is one order to remember. Read it as *what → where → how* (the kernel
is the weight, the degree the polynomial order, the bandwidth the scale they set). A docstring's argument
bullets then follow **signature order**.

**Docstrings vs comments.** Prefer docstrings; reserve plain comments for non-intuitive *why*
(the atomic scatter-add, the shift-formula recovery, the power-sums-about-the-bin-center choice,
"wider-than-data because this is a reduction").

**Marking what is not yet built.** A capability that is merely *enabled* or *planned* is marked —
forthcoming / future / a recorded gap — never stated as present.

## Documentation policy — estimand first, estimator second

Every operation estimates a population quantity from data, so the rulebook's
`RULE.DOC.LEADING_IDENTITY` binds this package with `leading_identity` set to the **estimand**
(`.gbg/config.jl`) — the statistics counterpart of PF's *physical meaning first* (`docs/ARCHITECTURE.md`
§2, §7): the population quantity, in standard probability notation, before the estimator.
`conditional_mean` is the reference: lead with
``m(x) = \mathbb{E}[Y \mid X = x] = \int y f_{Y\mid X}(y\mid x)\,\mathrm{d}y`` (the ``L^2``
projection), *then* the local-linear estimator. The intrinsically finite-sample objects — a bandwidth
selector, `merge!`, `scale!` — say so rather than force a population form.

- **Estimator-variant blocks vs Arguments — by number of variant axes.** When an estimator has a
  **single** variant axis whose branches are genuinely different estimators — a primary `rule` selector
  (`bandwidth`) or a `method` refinement (`conditional_variance`) — give each branch its own
  `**Estimator** (rule = …)` / `**Estimator** (method = …)` block stating that branch's formula. When the
  estimator varies along **several independent axes at once** — `density`'s `degree` (KDE vs
  local-likelihood) × `start` (semiparametric) × `adaptive` × `support` — there is no single axis to break
  out: give **one** unified `**Estimator.**` paragraph that names the branches, and put each axis's formula
  in its own Arguments bullet. Block-vs-Arguments tracks *how many axes*, not primary-vs-secondary — a
  secondary `method` can still earn blocks. *(Recorded so the docstring-structure audit does not re-flag
  the `density`-vs-`bandwidth` asymmetry: it is justified by single- vs multi-axis, not a defect.)*
- **Symbol bound to code on first mention, in apposition** (`RULE.TEX.ONE_SPAN`): `density` (``f_X``),
  `conditional_mean` (``m(x)``), `transition_density` (``p_\tau(y\mid x)``), `expectation`
  (``\mathbb{E}[g(X)]``). Where the formula says more than the name, the name goes inside it as
  `\texttt{…}`.
- **Math:** never Unicode-glyph approximations on a rendered surface; plain `#` comments stay Unicode.
- **Signatures and wide formulas:** a long docstring signature is split with a hanging indent aligned under
  the first argument (PF's `anterpolate!` convention); a wide display formula may also break with `split`.
- **Removing a thin space:** replace it with a plain space, never delete it (`\sigma\,x` deleted glues
  into the invalid `\sigmax`).
- **LaTeX forms that are correct — do not "fix":** `\,\mathrm{d}x`; `\quad`/`\qquad`; `\mathrm{Re}`
  vs `\operatorname{Re}` before a delimiter; `\sim\!` where it is neither before a sized delimiter nor
  doubled. When unsure, consult a standard LaTeX reference and KaTeX's support table, never taste.
- **Register:** the rulebook's "option", never "knob" (`RULE.DOC.REGISTER`) reaches here also the `#`
  comments, the four rationale docs and the `@testset` names; the bandwidth is a post-hoc *choice*.
- **No internal cross-references on the user-facing surface.** Beyond the rulebook's ban on section
  numbers in docstrings: a thrown error message is user-facing too, and neither carries a repo path
  (`docs/ARCHITECTURE.md`, `inspiration/…`); name the published page in prose ("the Architecture page").
  The `§4.x` estimator map stays in `#` comments, `CLAUDE.md` and the design docs.
- **Cross-referencing** follows the rulebook (`RULE.DOC.LINKS`, `RULE.DOC.XLINK`,
  `RULE.DOC.LINK_UNIFORMITY`); the package adds: never link the docstring's own binding, its own argument
  names, `Symbol`/option tokens (`:linear`) or struct fields; the **at-home test** — a sentence equally at
  home in the target's own docstring is cut and linked, but the function's own estimand formula and
  option-incompatibility clauses always stay (estimand-first outranks DRY); the link forms — docstring to
  symbol fully qualified, `` [`name`](@ref GBGStreamingStatistics.Module.name) ``; docstring to page a
  quoted-header prose link, `[the X](@ref "Header")`, never `@ref` on a bare module name; page to symbol
  bare; page to page relative, `[X](X.md)`; never a link inside a formula. A mistyped `@ref` fails the
  build, so a link cannot rot; a link adds depth one hop away, never a fact the reader needs.
- **Standard-error readers delegate their arguments** — "Arguments: as for [`f`](@ref …)", naming any
  argument of `f` they lack ("without `normalize`") — the form the rulebook sanctions since GBGCodingTools
  rc.59 (the owner, 2026-07-18 and 2026-10-04). The parent's link sits on that sentence, its one link in the
  docstring.
- **Correctness — the package's checklists.** Every estimator is consistent for its stated estimand under
  the stated assumptions. The normalizations to pin: ``1/n``, ``1/(nh)``, ``\int K = 1``,
  ``\mu_2(K) = \int u^2 K``, ``R(K) = \int K^2``, the bias and variance constants — a stray ``1/h`` or
  kernel-moment factor is the easy slip (the statistics analog of PF's stray ``1/N``). The category
  confusions to guard against: a **density** vs a **measure** or **CDF** (``f_X`` / ``\mu`` / ``F_X``); an
  **estimator** vs its **estimand** vs a single **estimate** (``\hat m(x)`` / ``m(x)`` / the number);
  ``\mathbb{E}[Y \mid X]`` (a random variable) vs ``\mathbb{E}[Y \mid X = x]`` (a number); a
  **conditional** density vs a **joint** or **marginal**; the overloaded "**kernel**" (the smoothing ``K``
  vs a transition kernel ``p_\tau``); a **pushforward density** vs a **pushforward measure**; and an
  implicit **convergence mode** (a.s. / in probability / ``L^2`` / in distribution) under a consistency or
  rate claim.
- **Citations** follow the rulebook: author–year in the prose, the full entry under `# References`,
  identical to `docs/ARCHITECTURE.md` §References. The private notes under `inspiration/` are working
  material, never cited.

**Doc-edit checklist — the package's own part.** The rulebook's battery runs in every docs build; this
is what it cannot run. On every documentation change, most strictly on `ARCHITECTURE.md`: (1) check each
estimator is consistent for its estimand; (2) enumerate the arguments *from the signature*, not the prose
— that is how `density`'s query point `x`, bandwidth `h` and `kernel` once slipped past a symbol-only
sweep; (3) every code identifier is backticked, never bare; (4) a summary that names only a subset of what
the code now provides is stale — the surfaces are the four rationale docs, every `src/` docstring and every
`docs/src/` page.

## Verification workflow

Before changing tests, read `docs/TESTING.md`; before optimizing/benchmarking, read
`docs/PERFORMANCE.md`. The nontrivial statements live formally in `docs/FORMALIZATION.md` — the
chain is *prose ↔ formalization ↔ test*.

### Definition of done — for a new or changed operation

Treat the feature as unfinished until each is **done and verified**:

1. **Correctness tests** against an *independent* analytic oracle — Gaussian marginals; an
   **Ornstein–Uhlenbeck process** (drift, diffusion, transition density known exactly) for §4.3–4.4
   — plus a **consistency / rate test** (error falls at the theoretical ``n^{-2/5}``), a
   **mergeability test** (accumulate-then-merge `=` accumulate-together, for the estimate **and** any
   standard error — every query-point estimator goes through `check_point_estimator`, `merged` a **two-or-more**-way
   merge), and — for an in-place op — an `@allocated` assertion in `test/allocations.jl`. **Every
   population-quantity estimator ships a standard error** — the `*_standard_error` reader and the one-call
   paired `*_with_standard_error` (a selector like `bandwidth` and the finite-sample `merge!` / `scale!` are
   the only exemptions, as for estimand-first). Each standard error carries the `check_standard_error`
   contract: positivity, scalar `=` vector, the
   ``n^{-1/2}`` **scaling** (merging ``K`` copies shrinks it by ``\sqrt K`` exactly), the ``h^{-1/2}``
   bandwidth scaling for a kernel SE (`check_se_bandwidth_scaling`, asymptotic — flat density, small ``h``),
   and a **calibration** — the predicted SE matches the empirical sampling spread of the estimate across
   independent datasets (Monte-Carlo). **GPU-vs-CPU parity** under `allowscalar(false)` once a device path exists.
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
   `observe!` (batch + scalar) and `merge!` (elementwise addition); global central moments
   recovered by the **shift formula** and exposed publicly as `mean`/`var`/`std`/`central_moment`
   (extending `Statistics`; these now live in `Expectations` §4.B as the ``g = x^p`` quadrature case). Tested:
   Gaussian moment recovery (exact vs direct, and vs the true moments), mergeability, bounded
   allocation, Float32→Float64.
2. **`Histograms.Histogram2D` + `observe_pairs!`** — the joint two-time summary: `nx×ny` cells, per-cell
   counts + order-2 bivariate power-sums ``S_{rs}`` about the cell center (the regression sufficient
   statistic). Tested: independent-Gaussian oracle (marginals + covariance), mergeability, allocation.
3. **The portable scatter kernel** — `observe!` / `observe_pairs!` are one KernelAbstractions + Atomix
   atomic scatter-add (one work-item per sample), allocation-free up to the launch, **verified identical
   on CPU and the AMD RX 6700 XT** (the gated ROCm suite, `test/gpu/`, under `allowscalar(false)`,
   including heavy atomic contention).
4. **`MarginalDensity.density`** (§4.1) — one call over the grid of **local-log-polynomial `degree`**
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
   `test/marginaldensity.jl` (oracle recovery, the start's wide-``h`` bias removal, normalization,
   mergeability, rejected option combinations).
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
6. **`ConditionalExpectations`** (§4.2) — `conditional_mean` ``\equiv m(x)`` and `conditional_variance`
   ``\equiv \sigma^2(x)`` by local-polynomial regression from `Histogram2D` over the shared
   `LocalPolynomial.local_fit` (the ``T_\ell``/``P_\ell`` from per-X-column counts and response totals),
   with the `degree` (0/1/2) and `kernel` options; `conditional_variance` defaults to the bias-reduced
   **Fan & Yao** residual estimator (`method = :residual`; the plug-in is `:plugin`), with the residual-based
   `conditional_variance_standard_error` ``= \sqrt{\ell(x)\,(\mu_4-\sigma^4)}`` (the SE of the squared-residual
   conditional mean; ``\mu_4`` from the order-2 quadrature, MC-calibrated on Gaussian and heavy-tailed
   responses). Tested against the
   bivariate-Gaussian oracle (``m`` recovers ``\rho x`` to ``\lesssim 2\%``; ``\sigma^2`` recovers
   ``1-\rho^2``) and a heteroscedastic curved-mean oracle (Fan–Yao recovers ``\sigma^2(x)`` and beats the
   plug-in in the curved tail). The **general reader** `conditional_expectation`
   ``\equiv \mathbb{E}[g(Y)\mid x]`` and the **conditional pushforward** `conditional_pushforward_density`
   ``\equiv f_{g(Y)\mid X}`` land too — `conditional_mean` folded in as the ``g = \text{identity}`` case,
   both **composing the two atoms** (a per-X-column y-quadrature response regressed in x). The sampling
   `conditional_expectation_standard_error` ``= \sqrt{\ell(x)\,\mathrm{Var}(g(Y)\mid x)}`` lands with them —
   the local-regression **leverage** ``\ell(x) = \sum_j N_j w_j^2`` a new `LocalPolynomial._leverage` readout
   (the design's plain **and squared-kernel** moments) — the reference SE mechanism for the whole kernel
   family; the marginal `expectation_standard_error`'s ``1/n`` leverage becomes the local one. Tested against
   the bivariate-Gaussian ``\mathbb{E}[\cos Y\mid x]`` and lognormal-pushforward oracles, the exact fold, the
   degree-0 identity coincidence with `conditional_density`, and a **Monte-Carlo** SE validation (predicted
   SE = empirical sampling spread across independent datasets).
7. **`IncrementMoments` — archived** (§4.4). The conditional increment moments and Kramers–Moyal
   functionals (`increment_moment`, `kramers_moyal`, `drift`, `diffusion`; `pawula_quotient` / `markov_test`
   were stubs) were built and OU-tested, then **archived** (in git history at commit `37006b5`, `archive/`; out of
   the tree since 2026-10-09) to focus the core on the
   summary-as-measure and its two readout families. Restore when a consumer needs them and `Histogram2D`
   raises its per-cell order for ``n \ge 3``.
8. **`ConditionalDensity`** (§4.3) — `conditional_density` ``\equiv f_{Y\mid X}`` and `transition_density`
   ``\equiv p_\tau(y\mid x)`` by the binned double-kernel (Rosenblatt) from the `Histogram2D` cell
   counts: the joint is reduced at ``x`` to a `kernelx`-weighted Y-profile, then smoothed in ``y`` by
   `kernely` (the surviving ``1/h_2`` gives ``\int \hat f\,\mathrm dy = 1`` exactly for a symmetric
   `kernely`). Supplied bandwidths `hx`, `hy`; for a **bounded/positive response** (a droplet radius as
   ``Y``) an **asymmetric** `kernely` (`gamma_kernel` / `beta_kernel`) is boundary-bias-free in ``y``,
   with `normalize = true` for an exact integral. Each carries its pointwise **standard error**
   `conditional_density_standard_error` / `conditional_density_with_standard_error` (the two-term
   double-kernel variance ``R(K_1)R(K_2)\hat f/(n h_1 h_2 f_X) - \hat f^2/n``, Hyndman–Bashtannyk–Grunwald
   1996 — ``n\hat f_X`` free from the estimator's denominator; `transition_density` aliased), MC-calibrated
   against the bivariate-Gaussian spread. Tested: the
   bivariate-Gaussian conditional ``N(\rho x, 1-\rho^2)``
   (normalization + pointwise) and the **Ornstein–Uhlenbeck** transition density
   ``N(x e^{-\theta\tau}, s^2(1-e^{-2\theta\tau}))``; scalar/vector forms, `transition_density`
   equivalence, off-support NaN, mergeability, and the **bounded-response** oracle (``Y\sim\mathrm{Exp}(1)``:
   the gamma `kernely` recovers the boundary the Gaussian leaks, `normalize` ``\Rightarrow`` ``\int\hat f = 1``).
9. **The cone — `scale!` + weighted `merge!`** (§10 hyperstatistics / forgetting) — the unnormalized
   summaries are scalable, so they form a convex cone, not merely an additive monoid: `scale!(summary, w)`
   reweights every count and power-sum (the scalar action), and `merge!(into, from; weight)` forms
   ``\text{into} \mathrel{+}= \text{weight}\cdot\text{from}`` (the default `weight = 1` the plain monoid
   add). On a **Float-count variant** (`count_type = Float64`; the default `UInt64` cannot hold a fractional
   ``w\,n``) for both histograms. Realizes the **hyperstatistical** parameter-mixture density
   ``P(x) = \sum_k w_k f_k(x) = \sum_k (w_k/N_k) T_k`` (superstatistics; Beck & Cohen 2003) — keep one
   summary per ``\beta``, choose the mixing law ``f(\beta)`` *post hoc*, re-weight without re-accumulating —
   and **exponential forgetting** (`scale!(s, λ)` before the next block). Tested (`test/hyperstatistics.jl`):
   the cone laws, the discrete Gaussian-mixture oracle, the **χ²-superstatistics → Student-``t``** Beck–Cohen
   oracle (a Gamma-weighted Gaussian superposition recovers the heavy ``t`` tail), exact exponential
   forgetting, the `Histogram2D` cone with `conditional_mean` scale-invariance, allocation, and GPU parity
   (incl. the atomic-**float** add of Float-count `observe!`).
10. **`LocalPolynomial.local_fit`** (the shared jet primitive) — the kernel-weighted local-polynomial
    design solve over a binned measure, returning the jet ``(\beta_0, \dots, \beta_d)`` (value, slope,
    ``2\beta_2`` the curvature) by the Hankel normal equations. The single *linear* core:
    `conditional_mean` is its ``\beta_0``, the local-likelihood `density` its IRLS inner solve, and it gave
    `ConditionalExpectations` the `degree` / `kernel` options and the degree-2 curvature ``m''`` (the §11 architecture
    decision, landed). Tested (`test/localpolynomial.jl`): exact jet recovery on a polynomial mean (value
    / slope / curvature), the degree-0 Nadaraya–Watson ratio, degenerate ``\to`` NaN, kernel-normalization
    cancellation.
11. **The response-aware local bandwidth ``h^\star(x)``** (§4.5) — `bandwidth(::Histogram2D, x;
    rule = :local)` ``= [R(K)\sigma^2(x)/(\mu_2(K)^2 m''(x)^2 f_X(x) n)]^{1/5}``, composing the Fan–Yao
    ``\sigma^2(x)``, the degree-2 curvature ``m''(x) = 2\beta_2``, and the X-marginal ``f_X(x)`` from the
    one joint summary at a curvature-rate pilot ``\hat\sigma_X n^{-1/9}``. Tested (`test/localbandwidth.jl`):
    the quadratic-mean Gaussian oracle (closed-form ``h^\star``) recovered near the centre, NaN off
    support, scalar ``=`` vector, mergeability.
12. **The conditional convergence-rate tests + foundation coverage** — `conditional_mean` at the
    ``n^{-2/5}`` ``L^2``-projection rate (ISE ``\\propto n^{-4/5}``; curved ``m(x) = \\sin x``, so the
    local-linear fit carries a genuine ``O(h^2)`` bias) and `conditional_density` at the slower 2-D
    ``n^{-2/3}`` rate (`test/conditional_convergence.jl`) — closing the rate column of the coverage matrix;
    `drift` / `diffusion` / `conditional_variance` inherit `conditional_mean`'s rate (affine and
    second-stage functionals of the same fit), so they need no separate slope test. Plus the direct
    kernel-constant oracle (`test/kernels.jl` — ``\\int K = 1``, ``\\mu_2(K)``, ``R(K)`` for the Gaussian /
    Epanechnikov kernels, by an independent quadrature) and the two remaining parametric starts
    (`exponential_start`, `inverse_gaussian_start`, under matching-distribution oracles in
    `test/marginaldensity.jl`).
13. **`Expectations`** (§4.B) — the **power-sum quadrature** readout family, the second of the two: the
    expectation `expectation` ``\\equiv \\mathbb{E}[g(X)]`` of a general smooth functional, its sampling
    `expectation_standard_error` (``\\sqrt{(\\mathbb{E}[g^2] - \\mathbb{E}[g]^2)/n}``), and
    `pushforward_density` ``\\equiv f_{g(X)}`` — the density of ``Y = g(X)``, written as
    ``\\mathbb{E}[K_h(g(X) - y)]`` so no ``g^{-1}`` / Jacobian is formed and a **non-monotone** ``g`` (and
    ``g'=0`` caustic) is handled. All route through the one shared `_quadrature` atom (**owned in
    `Expectations`**), alongside `mean` / `var` / `std` / `central_moment` — **moved here from `Histograms`**,
    **folded** as the ``g = x^p`` case (the literal two-family lens;
    the boundary is that the families read *different sufficient statistics* — power-sums vs
    counts+locations, settled against the tempting universal-map over-unification). A general ``g``'s
    coefficients come from a **dependency-free** degree-`order` local interpolant (exact on monomials).
    Tested (`test/expectations.jl`): the Gaussian ``\\mathbb{E}[\\cos X] = e^{-1/2}`` / MGF oracles,
    polynomial exactness (``= \\tfrac1n\\sum_i g(X_i)`` to ``10^{-10}``), the lognormal and non-monotone
    ``\\chi^2_1`` pushforwards, the SE reduction to ``s/\\sqrt n``, and mergeability.
14. **Standard errors — the paired readers, the two-term band, and the testing contract.** The
    estimate-plus-error `*_with_standard_error` pairing across all three families
    (`density_with_standard_error`, `expectation_with_standard_error`,
    `conditional_expectation_with_standard_error`) — each the `(value, standard_error)` NamedTuple from one
    shared computation, so the estimate and its error describe the same estimator by construction.
    `density_standard_error` is now the **two-term** Wand–Jones KDE variance
    ``\\sqrt{R(K)\\hat f/(nh) - \\hat f^2/n}`` and **tracks the estimator's options** (`adaptive` → the local
    Abramson ``h(x)``, `support` → the boundary-restricted ``R_a/a^2``, `start` leading-order-neutral;
    degree-0/1 symmetric-kernel). The uniform **standard-error contract** (`check_standard_error`,
    `test/testcontracts.jl`) runs on every SE — positivity, scalar ``=`` vector, the exact ``n^{-1/2}``
    copy-merge scaling, the asymptotic ``h^{-1/2}`` bandwidth scaling (`check_se_bandwidth_scaling`, the
    kernel SEs), and a Monte-Carlo **calibration** (predicted SE ``=`` the empirical sampling spread) — and
    each SE also goes through `check_point_estimator` (SE mergeability, a **two-or-more-way** merge). A
    complementary **coverage** check (across every family) asks the ``\\hat\\theta \\pm 1.96\\,\\mathrm{SE}``
    band to bracket the **oracle truth** (not just the estimate's own biased mean) over an interior grid —
    aggregate and lenient, since a variance band cannot cover the ``O(h^2)`` bias, so it runs in the
    variance-dominated regime (a deliberately small ``n`` for the 2-D `conditional_density` /
    `conditional_pushforward_density`, whose bias dominates the tight large-``n`` SE). Plus
    the systematic **invalid-combination audit**: every `throw(ArgumentError(…))` guard carries a
    `@test_throws` (the valid option space is covered per-option and by the interacting pairs, deliberately
    not exhaustively). Tested in `test/marginaldensity.jl`, `test/expectations.jl`,
    `test/conditionalexpectations.jl`.
15. **HDF5 persistence — the summary is a durable, portable artifact.** A weak-dependency package extension
    (`ext/GBGStreamingStatisticsHDF5Ext.jl`, loaded when `HDF5` is present — `HDF5` is **not** a core dep):
    `write_histogram(parent, name, summary)` / `read_histogram(parent, name)` (+ the `parent[name] =
    summary` setindex! sugar) map a `Histogram1D` / `Histogram2D` to a self-describing HDF5 **group** (the
    arrays → `counts` / `powersums` datasets, the scalars → attributes, a `type` / `version` tag).
    **Handle-based only** — the consumer opens and owns the file, so summaries land beside other data in one
    file, many groups. Device arrays move to host on write; the count element type (the `count_type = Float64`
    cone variant) is preserved; the reload is **bit-exact** and still merges. Tested in `test/hdf5.jl`
    (1-D / cone / 2-D exact round-trip, the setindex! sugar, functional equality after reload,
    reload-preserves-merge, self-describing attributes, the error paths).

**Next, in dependency order:**

The *every estimator ships a standard error* requirement (the Definition of done and
`docs/ARCHITECTURE.md` §7) is now **complete** — every population-quantity estimator carries its
`*_standard_error` and paired `*_with_standard_error`, each under the `check_standard_error` contract. The
last one, `conditional_variance_standard_error` (the residual approximation ``\\sqrt{\\ell(x)\\,(\\mu_4 -
\\sigma^4)}``, with ``\\mu_4`` read off the order-2 quadrature of the smooth ``(Y-\\hat m)^4`` — no
4th-moment storage), is Done 6.

1. **Density-estimation methods from the prior-art survey** (`inspiration/density-software/ALGORITHMS.md`),
   in priority order — all mergeable functionals of the summary: the **POT–GPD heavy-tail splice**
   (domain-critical for intermittent tails; a binned GPD-MLE on the tail bin counts that extrapolates
   beyond the histogram's `[lower, upper)` range, plus an optional mergeable top-k max-heap for exact
   Hill/Pickands tail indices); the **P-spline penalized-Poisson** density (Eilers–Marx, counts-native
   global smoothing); **max-entropy from the moments** (uses the power-sums); the **self-consistent**
   Fourier-filter KDE (Bernacchia–Pigolotti) and **orthogonal series** (Hermite/DCT, seeding §4.6). The
   asymmetric gamma/beta kernels are **done** (Done 4). Deliberate skips recorded there: shape-constrained
   MLE (shape-as-premise, invariant 2), taut-string (a mode-count diagnostic), DEFT, kNN, neural/DP-mixture.

**Audit / accuracy passes (recorded for a return — reviews of the existing surface, not new estimators):**

- **Pushforward coverage — every axis of every estimator.** Enumerate the (estimator × {conditioning ``X``,
  response ``Y``}) grid and settle which pushforwards exist, which are sensible-but-missing, and which are
  deliberate absences. Have: `pushforward_density` (``f_{g(X)}``) and `conditional_pushforward_density`
  (``f_{g(Y)\\mid X}`` — which on a lagged summary already yields the transition pushforward
  ``p_\\tau(g(y)\\mid x)``; verify it is exposed / aliased). Open cases to adjudicate: a pushforward on the
  **conditioning** axis (``f_{Y\\mid g(X)}``, likely ill-posed for non-monotone ``g`` — a candidate
  deliberate absence); a joint / both-axes pushforward (no joint density readout today); a
  `transition_density` pushforward alias. Output: the filled grid + recorded absences with reasons.
- **Option parity — every option, both axes.** Tabulate each shared option (`degree`, `kernel`, `adaptive`,
  `support`, `start`, `normalize`, `method`, `binning`) against every estimator, and for the joint estimators
  against the **``x`` (conditioning) vs ``y`` (response) axis** separately; mark each absence *deliberate*
  (theory) or a *gap to close*. Known asymmetries to settle: `adaptive` / `support` / `start` exist for
  `density` but not the conditional densities; `conditional_density` is a fixed degree-0 double-kernel (no
  local-polynomial `degree` on either axis); the asymmetric kernels reach `density` (`kernel`),
  `conditional_density` (`kernely`), and degree-0 `conditional_mean` (`kernel`) — is the ``x`` axis at
  parity? Output: the parity matrix + a short list of justified additions.
- **Conditional accuracy from the unused per-cell moments.** The conditional readers (`conditional_mean` /
  `conditional_variance` / `conditional_expectation` / `conditional_density`) consume only the cell counts,
  the **``Y``-power-sums** ``S_{01}`` / ``S_{02}``, and the **``X`` cell centers** — never the ``X``-direction
  per-cell moments ``S_{10}`` (within-``x``-bin centroid), ``S_{20}`` (within-``x``-bin variance), ``S_{11}``
  (within-cell ``xy``-covariance) that `observe_pairs!` already accumulates. Investigate using them to raise
  accuracy: an **``X``-centroid correction** to the local-polynomial design (place each column at
  ``c^x_j + S_{10}^{(j)}/N_j`` instead of the bin center — the 2-D analog of `density`'s linear/centroid
  `binning`, ``O(\\Delta_x^2)`` vs ``O(\\Delta_x)``); an ``S_{11}``-based within-cell bilinear correction to
  the regression. Weigh the gain against the cost, and settle the standing question of whether the
  currently-dead ``X``-moments earn their storage (else drop them). Connects to the deferred `Histogram2D`
  order-raise.

Deferred (recorded, not built): the **FFT-amortized post-hoc smoothing** (§4.5 — the ``O(B\\log B)``
convolution of the bin counts with the kernel, amortized over a bandwidth *grid*; a pure optimization over
the exact direct summation, itself dwarfed by accumulation — build it when an ``h^\\star(x)`` grid sweep
makes it bite; `FFTW` stays for the ISJ DCT, the unused `AbstractFFTs` dep was dropped); the ECF summary
(`docs/ARCHITECTURE.md` §2 — likely never); the
projection / SFI–ULI road (§4.6); `HistogramND` (§3, wait for a multi-conditioning consumer);
vector-valued conditioning (lean scalar-channel — let the consumer loop components); the **bulk scalar time-correlation function**
(consumer-side — a `Histogram2D` per lag is the wrong cost for one cross-moment; §6); raising
`Histogram2D` beyond its order-2 per-cell power-sums (the higher Kramers–Moyal / Pawula moments (``M_4``–``M_6``) — the archived §4.4 functionals — will be
its first consumer); the **leave-one-out cross-validation** bandwidth selectors
(likelihood-CV, least-squares LSCV — a binned summary's ties drive their objective to zero bandwidth; the
built `:aic` / `:sj` are the data-driven selectors instead, and a within-bin-spread correction from the
power-sums could rescue CV later).

## Working agreements

- The git and CI conventions are the rulebook's (`.gbg/policy/git-conventions.md`): commit and push
  directly on `master`, the versionless trailer, CI tests by hand and on release tags with the full matrix
  run locally before every push, and the docs built and deployed on every push. The repository is private
  and synced across the two machines through GitHub.
