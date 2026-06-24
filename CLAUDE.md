# CLAUDE.md — GBGStreamingStatistics

Julia package: a backend-agnostic library for **online nonparametric statistics**. It estimates
marginal densities, conditional expectations, conditional and two-time transition densities, and
standard functionals of these (Kramers–Moyal coefficients, drift/diffusion, Pawula diagnostics)
from a **stream of samples that is never stored** — each sample is observed once, folded into a
fixed-size summary, and discarded. **Single node, shared-memory only** (CPU multithreading and a
single GPU; distributed/multi-GPU out of scope). Development is **GPU-first** but the code is
**backend-agnostic** (CPU `Array`, CUDA, ROCm) through KernelAbstractions.

The architectural spine is **empirical measure → mergeable summary → post-hoc-smoothed estimator**:
a fine histogram with per-bin moments is accumulated online and is *response-agnostic*; every
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
> and power-sums about the bin center, `merge!` by addition, allocation-free), the marginal density
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
  `observe!` / `merge!` mergeable-monoid interface only. The concrete summaries and the portable
  KernelAbstractions + Atomix scatter kernels live in `Histograms`.
- `src/Histograms.jl` — the summaries (lower half): `Histogram1D` / `Histogram2D` — linear (WARPing)
  binned empirical measures over a fixed range, each bin/cell carrying its **count** and **power-sums
  about the bin center** ``S_p`` / ``S_{rs}`` (the local sufficient statistic for moments and
  local-polynomial fitting). Accumulated by one **KernelAbstractions + Atomix atomic scatter-add**
  (`observe!` / `observe_pairs!`, one work-item per sample, CPU and GPU alike); `merge!` is elementwise
  addition; `bin_index` is the binning map (exposed for pre-quantization); global central moments are
  recovered post hoc by the **shift formula**. Owns `KernelAbstractions`, `Atomix`.
- `src/DensityEstimation.jl` — §4.1 `density` ``\equiv f_X`` (binned KDE / WARPing).
- `src/Regression.jl` — §4.2 `conditional_mean` ``\equiv m(x)``, `conditional_variance`
  ``\equiv \sigma^2(x)`` (local-linear, Fan & Gijbels; Nadaraya–Watson kept for reference). The
  **estimand-first reference module** — write its docstrings to the documentation policy below.
- `src/ConditionalDensity.jl` — §4.3 `conditional_density` ``\equiv f_{Y\mid X}``,
  `transition_density` ``\equiv p_\tau(y\mid x)`` (double-kernel, Rosenblatt).
- `src/MarkovInference.jl` — §4.4 `increment_moment` ``\equiv M_n(x,\tau)``, `kramers_moyal`
  ``\equiv D^{(n)}``, `drift` ``\equiv D^{(1)}``, `diffusion` ``\equiv D^{(2)}``, `pawula_quotient`
  ``\equiv Q(x)``, `markov_test`. Model-neutral functionals; the SDE *reading* stays consumer-side.
- `src/Bandwidth.jl` — §4.5 `bandwidth`: global rules (Silverman; Sheather–Jones) and the
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

1. **Accumulation is response-agnostic; smoothing is post hoc.** No summary stores a bandwidth,
   kernel, or channel choice. An API that bakes a bandwidth into accumulation is banned. *(This is
   the invariant that makes "no reruns" true.)*
2. **Samples in, statistics out — no domain, no model.** Domain-neutral vocabulary (sample /
   variate / response / lag / conditioning value), never particle / tracer / droplet / velocity. No
   stationarity, Markov, diffusion, or parametric assumption is a *premise*; such structure is only
   ever *tested* (the Pawula ratio, the Markov / Chapman–Kolmogorov test).
3. **Everything is a mergeable monoid.** Every summary supports `merge!`; the unnormalized summary
   is additive over the stream. No non-mergeable accumulator is admitted. Tested: accumulate-then-
   merge `=` accumulate-together.
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
descriptive names without abbreviations; moderately long lines are fine. When unsure about
Base/stdlib/package behavior, check official docs — don't guess. For a linear solve use
`\`/factorizations/`ldiv!`, never an explicit matrix inverse (the local-polynomial normal equations
are a 2×2/3×3 solve). Operators and extended Base/`Statistics` functions are imported before being
extended (`import Base: merge!`); no type piracy.

**Docstrings vs comments.** Prefer docstrings; reserve plain comments for non-intuitive *why*
(the atomic scatter-add, the shift-formula recovery, the power-sums-about-the-bin-center choice,
"wider-than-data because this is a reduction"). One docstring per operation, on the bare binding:
`function f end` for package-owned generics, the bare imported name for Base/`Statistics` functions.

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
  never math in backticks (backticks are for code identifiers only). Plain `#` comments stay Unicode.
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
renders the broken page without complaint); **(5) Current & complete** — no stale *deferred/none-exist* claims, and a
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
   recovered by the **shift formula** (`_sample_mean`/`_variance`/`_std`/`_central_moment`). Tested:
   Gaussian moment recovery (exact vs direct, and vs the true moments), mergeability, bounded
   allocation, Float32→Float64.
2. **`Histograms.Histogram2D` + `observe_pairs!`** — the joint two-time summary: `nx×ny` cells, per-cell
   counts + order-2 bivariate power-sums ``S_{rs}`` about the cell center (the regression sufficient
   statistic). Tested: independent-Gaussian oracle (marginals + covariance), mergeability, allocation.
3. **The portable scatter kernel** — `observe!` / `observe_pairs!` are one KernelAbstractions + Atomix
   atomic scatter-add (one work-item per sample), allocation-free up to the launch, **verified identical
   on CPU and the AMD RX 6700 XT** (the gated ROCm suite, `test/gpu/`, under `allowscalar(false)`,
   including heavy atomic contention).
4. **`DensityEstimation.density`** (§4.1) — binned KDE (Gaussian kernel), the estimand-first reference;
   recovers the Gaussian pdf to ``< 10^{-3}`` and integrates to 1.
5. **`Bandwidth.bandwidth`** (§4.5) — global Silverman/Scott plug-in from the recovered ``\\hat\\sigma``.
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

**Next, in dependency order:**

1. **The response-aware local bandwidth** ``h^\\star(x)`` (§4.5) — estimand-first, OU-oracle-tested,
   with its `FORMALIZATION.md` entry, from the `Histogram2D` joint summary. Plus the Fan–Yao
   conditional-variance refinement, the ``n^{-2/5}`` / ``n^{-2/3}`` rate tests, and the FFT-amortized
   smoothing (§4.5).
2. **The rest of `MarkovInference` §4.4** — `kramers_moyal` for ``n \ge 3``, the **Pawula `Q`**
   diagnostic (``M_4``, ``M_6``) and the all-orders Bell-polynomial finite-time correction, and
   `markov_test` (now unblocked — it composes the §4.3 transition density). The higher moments need the
   `Histogram2D` order-raise (Deferred, below).

Deferred (recorded, not built): the ECF summary (`docs/ARCHITECTURE.md` §2 — likely never); the
projection / SFI–ULI road (§4.6); `HistogramND` (§3, wait for a multi-conditioning consumer);
vector-valued conditioning (lean scalar-channel — let the consumer loop components); raising
`Histogram2D` beyond its order-2 per-cell power-sums (the Pawula ``M_4``/``M_6`` of *Next* 2 will be
its first consumer); the **weighted merge** for ensemble / hyperstatistical pooling with
parameter-dependent weights (a Float-count summary variant — equal-weight `merge!` is the built
special case; `docs/ARCHITECTURE.md` §10).

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
