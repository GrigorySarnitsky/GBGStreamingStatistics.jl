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
deliberate absences of the invariants). The conceptual model and the writing/docstring conventions
are recorded once — in `docs/ARCHITECTURE.md` and this file; the umbrella module
`src/GBGStreamingStatistics.jl` only points there.

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
  addition (and `scale!` / weighted `merge!` elementwise scaling — the cone, on the `count = Float64`
  Float-count variant); `bin_index` is the binning map (exposed for pre-quantization). The **pure summary**
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
  renormalization (`bounds`, Jones 1993), the **parametric start** (`start`, Hjort & Glad 1995 — the
  moment-fitted `normal_start` / `lognormal_start` / `gamma_start` / `exponential_start` /
  `inverse_gaussian_start`, or a user-supplied `ParametricStart`), exact `normalize`
  (``\int\hat f = 1``), and the **`binning`** location rule (`:linear` centroid default — ``O(\Delta^2)``
  via the stored ``S_1``; `:simple` bin-center, ``O(\Delta)``) — accepting any `kernel` from `Kernels`
  (the symmetric `gaussian_kernel` / `epanechnikov_kernel`, or the asymmetric `gamma_kernel` /
  `beta_kernel` for bounded/positive support), and `density_standard_error` (the pointwise variability
  band ``\mathrm{SE}(x) = \sqrt{R(K)\hat f_h(x)/(nh) - \hat f_h(x)^2/n}``, the two-term KDE variance — it
  **tracks the estimator's options**: the local Abramson ``h(x)`` for `adaptive`, the boundary-restricted
  ``R_a(x)/a(x)`` for `bounds`, `start` leading-order-neutral; degree-0/1 symmetric-kernel, the asymmetric
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
- Deferred / archived, recorded as commented-out includes (the rule against speculative structure):
  `# IncrementMoments.jl` (the §4.4 SDE-reading `drift` / `diffusion` / Kramers–Moyal functionals —
  **archived** to `archive/`, working and OU-tested; restore when a consumer needs them and `Histogram2D`
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

**Argument order — one convention, everywhere (internal and exported).** Every function orders its
arguments **object → query point(s) → smoothing controls → modifiers**: (1) the data operated on (the
`summary`, or the binned arrays — `counts` / `centers` / `locations` / `mass` / `response` / `P` / `D` /
`n` — in their natural order); (2) the **query point(s)** `x` (then `y` / the response point) where the
estimand is evaluated; (3) the smoothing controls — a **primary algorithm-selector** first where the function has one (`rule`,
which picks the whole estimator), then the fixed sub-order **bandwidth → kernel → degree** (`h` / `hx` /
`hy` / `b`, then `kernel` / `kernelx` / `kernely`, then `degree`); (4) the **secondary** modifiers and
refinements (`method`, `adaptive`, `bounds`, `start`, `normalize`, `binning`, `maxorder`, IRLS scratch, …)
in a stable order — a `method`-style refinement (it tweaks *how*, not *which* estimator) trails, unlike a
primary `rule`-style selector, which leads. For an exported function, (1)–(2) are positional and (3)–(4) keyword — the
keyword block follows the same bandwidth → kernel → degree → modifiers order. Read it as *what → where →
how* (the bandwidth is the primary scale, the kernel the weight, the degree the polynomial order). A
docstring's argument bullets then follow **signature order** — the natural consequence, and what
`/check-docs` checks.

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
  ``\equiv m(x)``, `transition_density` ``\equiv p_\tau(y\mid x)``, `expectation`
  ``\equiv \mathbb{E}[g(X)]``.
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
- **No horizontal scroll — hand-wrap what the renderer won't.** A rendered **code block** (a docstring
  signature) and a **display-math block** do not auto-wrap, so a long one overflows the page column with
  a horizontal scrollbar. **Signatures:** split across lines with a hanging indent aligned under the
  first argument (the `inspiration/GBGPeriodicFields.jl` convention — its `anterpolate!`), each rendered
  line within the content column ≈ 84 characters (Documenter strips the 4-space code-block indent, so
  ≈ 88 source characters). **Display formulas:** break a wide one with a multi-line LaTeX environment —
  `aligned` (align on `=`), `gathered`, or `split` — never one line packing several equations behind
  `\qquad`/`\quad`. Inline math reflows with the prose, so it is exempt.
- **Neutral reference register — never "knob".** A tunable argument is an **option** (or, where it reads
  better, a *keyword argument*, a *setting*, or a post-hoc *choice* — e.g. the bandwidth), never a "knob":
  the reference docs keep a strict, neutral tone, not a colloquial one. Applies to every authored surface —
  docstrings, the four rationale docs, the doc pages, `#` comments, and `@testset` names.
- **No internal cross-references on the user-facing surface.** A rendered **docstring** or a **thrown
  error message** is user-facing: it must not cite a design-doc **section number** (`§4.5`) or a **repo
  path** (`docs/ARCHITECTURE.md`, `inspiration/…`) — those are dangling, unnavigable, and they renumber.
  Keep it self-contained; where a rationale pointer earns its place, name the published page in prose
  ("the Architecture page"), never a `§`-number or path. The `§4.x` estimator map stays in its developer
  home — `#` comments, `CLAUDE.md`, and the design docs (which cite each other by section freely).
- **Cross-referencing — discriminating-minimalist; a link earns its place.** Live Documenter `@ref`
  hyperlinks are the *encouraged* way to point across the docs (the rule above bans only *dead* `§`-numbers
  and repo paths; a mistyped `@ref` fails the build, so it cannot rot). But a link is navigation, not
  decoration — **backtick every code identifier; promote it to an `@ref` only where a reader would
  genuinely jump and there is one clear destination**: (1) link a given target **at most once per
  docstring** (its first genuine mention — this kills duplicate targets scattered across bullets); (2)
  **never** link the docstring's own binding, its own argument names, `Symbol`/option tokens (`:linear`),
  struct fields, or `Base` names; (3) a **family of built-ins** (the kernels, the starts) is **one**
  catalogue destination, not five — name every member in plain backticks (so `?f` at the REPL still shows
  the menu), `@ref` the **default** member, carry **one** catalogue link (the owning page or a same-page
  section, via the quoted-header form `[the Kernels page](@ref "Kernels")`), and `@ref` the **extension
  type** once with a *supply-your-own* clause (“… or a custom `AsymmetricKernel`”) — 2–3 links per family
  bullet, never five, never one; (4) the **at-home test** for describe-vs-link — if a sentence would be
  equally at home in the *target's* own docstring, cut it and link; the function's **own** estimand formula
  and option-incompatibility clauses always stay (estimand-first outranks DRY); (5) **link forms** —
  docstring→symbol fully qualified `` [`name`](@ref GBGStreamingStatistics.Module.name) ``,
  docstring→page/section a quoted-header prose link `[the X](@ref "Header")` (never `@ref` a bare module
  name for a page — it dangles), page→symbol bare `` [`name`](@ref) ``, page→page relative `[X](X.md)`;
  never a link inside a formula. The docstring must read **complete at the REPL** where links are dead — a
  link adds depth one hop away, never a fact the reader needs.
- **Correctness is its own check, distinct from formatting.** Verify the estimator is **consistent
  for the stated estimand** under the stated assumptions, and that **every normalization is exact**:
  the ``1/n``, ``1/(nh)``, ``\int K = 1``, kernel moments ``\mu_2(K) = \int u^2 K``, roughness
  ``R(K) = \int K^2``, the bias/variance constants. A stray ``1/h`` or kernel-moment factor is the
  easy slip here (the statistics analog of PF's stray ``1/N`` / Parseval factor).
- **Citations:** context-standard author–year; name up to three authors, `First et al.` for four or
  more; **never cite an inaccessible source** (the private notes under `inspiration/` are working
  material, not citable). The canonical reference list lives in `docs/ARCHITECTURE.md` §References;
  docstrings carry the full reference inline, verbatim-consistent with it. **A reference *block* of two
  or more citations is a bullet list — one citation per `-` line, under a `References:` label and a blank
  line — never a semicolon-joined paragraph (which is unreadable in a rendered docstring); a lone
  citation stays inline as `Reference: …`.**

**Doc-edit checklist — run all on every documentation change** (most strictly on `ARCHITECTURE.md`):
**(1) Correct** — every claim/formula/normalization true and precise, the estimator consistent for
its estimand, run as a *positive* per-formula verification, not a read-through; **(2) Self-contained**
— every symbol *bound* to its meaning in plain language at first use (enumerate symbols, cite the
binding clause), **and every signature argument documented** — positional and keyword alike, enumerated
*from the signature* (not the prose) so an "obvious" `x` or `h` cannot be skipped: each says what it is
(and for an option, how to set it), and each as its **own entry in the arguments list** — never relegated
to surrounding prose while its siblings are bulleted; an argument that surfaces only inside a formula, in a
trailing paragraph, or not at all, is the defect — this is how `density`'s query point `x` and bandwidth
`h` slipped past a symbol-only sweep, and its `kernel` sat in a closing paragraph instead of a bullet;
**(3) Notation** — one symbol per concept and one concept per symbol; standard
notation; **(4) Formatted (bidirectional)** — math is LaTeX never backticks, every code identifier
backticked never bare, and inline ``math`` stays on **one source line** (a wrapped continuation that
begins with a list/blockquote leader — `- `/`* `/`> `/`1.` — is parsed as a list item and silently
drops the span out of math mode; `docs/make.jl` lints for this and fails the build, since Documenter
renders the broken page without complaint); a **formula is a single LaTeX span**, never ``math`` split
onto an adjacent `code` continuation (a code name inside a formula is `\texttt{}`-wrapped), and
**backticks hold only well-formed Julia** (a name/type/call/expression), never math notation —
intervals, products, inequalities, `·`/`−`/`≤` glyphs; and **nothing horizontally scrolls** — a long
docstring signature is split with a hanging indent under the first argument, a wide display formula
broken with `aligned`/`gathered`/`split` (neither a rendered code block nor display math auto-wraps);
**(5) Current & complete** — no stale *deferred/none-exist* claims, and a
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
   `observe!` (batch + scalar) and `merge!` (elementwise addition); `bin_index`; global central moments
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
   were stubs) were built and OU-tested, then **moved to `archive/`** to focus the core on the
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
   add). On a **Float-count variant** (`count = Float64`; the default `UInt64` cannot hold a fractional
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
    Abramson ``h(x)``, `bounds` → the boundary-restricted ``R_a/a^2``, `start` leading-order-neutral;
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
