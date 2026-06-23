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

> **Status: scaffold.** The repository layout, conventions, and the full public API surface (as
> estimand-first docstrings on bare bindings and struct skeletons) are in place and load green; the
> estimator numerics are **declared, not yet implemented**. See *Current state and roadmap*.

## Modules and wiring

Two halves (`docs/ARCHITECTURE.md` §3, §4, §8). The **universal lower half** is pure numerics over a
stream and transfers verbatim to every consumer; the **estimator upper half** is the exportable
estimands.

- `src/Accumulators.jl` — the foundation (lower half): the `AbstractAccumulator` supertype, the
  `merge!` mergeable-monoid interface, the **Pébay (2008) central-moment recursion with
  Kahan–Neumaier compensated summation** (per bin), the backend-agnostic **KernelAbstractions +
  Atomix scatter-add primitive**, and `observe!` (fold a sample batch in). Allocation-free on
  preallocated carriers; runs on CPU and any GPU backend. Owns `KernelAbstractions`, `Atomix`.
- `src/Histograms.jl` — the carriers (lower half): `Histogram1D` / `Histogram2D` — linear (WARPing)
  binned empirical measures over a fixed range, each bin carrying its **count**, **per-bin
  standardized central moments to order ≈ 12**, and **per-bin power-sums** ``S_{rs}`` (the local
  sufficient statistic for local-polynomial fitting). `bin_index` (the binning map, exposed for
  pre-quantization), `observe_pairs!` (bin pre-formed lag-labelled pairs), `merge!`.
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
  `# CharacteristicFunctions.jl` (the ECF carrier — `docs/ARCHITECTURE.md` §2, likely never built),
  `# ProjectionInference.jl` (the SFI / ULI projection road — §4.6).

## Core invariants — do not violate

The PF discipline: a missing capability encodes a decision (`docs/ARCHITECTURE.md` §6).

1. **Accumulation is response-agnostic; smoothing is post hoc.** No carrier stores a bandwidth,
   kernel, or channel choice. An API that bakes a bandwidth into accumulation is banned. *(This is
   the invariant that makes "no reruns" true.)*
2. **Samples in, statistics out — no domain, no model.** Domain-neutral vocabulary (sample /
   variate / response / lag / conditioning value), never particle / tracer / droplet / velocity. No
   stationarity, Markov, diffusion, or parametric assumption is a *premise*; such structure is only
   ever *tested* (the Pawula ratio, the Markov / Chapman–Kolmogorov test).
3. **Everything is a mergeable monoid.** Every carrier supports `merge!`; the unnormalized summary
   is additive over the stream. No non-mergeable accumulator is admitted. Tested: accumulate-then-
   merge `=` accumulate-together.
4. **The hot path is allocation-free and backend-agnostic.** `observe!(summary, batch)` allocates
   nothing on a preallocated summary (enforced by an `@allocated` audit on the CPU). Storage is
   generic `AbstractArray`; no scalar indexing of device arrays; one portable KernelAbstractions +
   Atomix kernel does the atomic scatter-add binning on CPU and any GPU backend. The hot-loop enemy
   is **contention, not FP64** — privatize (per-block local histograms, merged hierarchically).
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
(the Pébay update, the Kahan compensation, the privatization choice, "wider-than-data because this
is a reduction"). One docstring per operation, on the bare binding: `function f end` for
package-owned generics, the bare imported name for Base/`Statistics` functions.

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
backticked never bare; **(5) Current & complete** — no stale *deferred/none-exist* claims, and a
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

**Scaffold stage.** Everything below is *declared and documented* (estimand-first docstrings on bare
bindings, struct skeletons with their fields) and loads green; **no estimator numerics are
implemented yet**. The test files are passing placeholders to be filled as each module lands.

Next steps, in dependency order (the lower half is the hard, write-once numerics that everything
depends on — build and lock it first):

1. **`Accumulators` + `Histograms`** (the universal lower half) — the Pébay+Kahan per-bin moment
   recursion, the linear/WARP binning, the per-bin power-sums, the KernelAbstractions + Atomix
   scatter-add primitive, `observe!` / `observe_pairs!` / `bin_index` / `merge!`. Allocation-free CPU
   path; mergeability + `@allocated` + GPU-parity tests.
2. **`DensityEstimation`** (§4.1) — `density` as the first end-to-end estimand-first reference,
   verified against a Gaussian oracle.
3. **`Regression`** (§4.2), **`ConditionalDensity`** (§4.3), **`MarkovInference`** (§4.4),
   **`Bandwidth`** (§4.5) — each estimand-first, OU-oracle-tested, with its `FORMALIZATION.md` entry.

Deferred (recorded, not built): the ECF carrier (`docs/ARCHITECTURE.md` §2 — likely never); the
projection / SFI–ULI road (§4.6); `HistogramND` (§3, wait for a multi-conditioning consumer);
vector-valued conditioning (lean scalar-channel — let the consumer loop components). Prior art to
mirror before writing `Accumulators`: `OnlineStats.jl`'s mergeable streaming-moment design.

## Working agreements

- The API is not frozen; do not build backward-compatibility shims.
- **Commit ⇒ push, on `master`.** This is a *private* repo synced across machines through GitHub, so
  work is committed **directly on `master`** (do not branch) and **pushed in the same step**; the
  docs site deploys on push. End every commit message with `Co-Authored-By: Claude Opus 4.8
  <noreply@anthropic.com>`; never use `[skip ci]`.
- **CI tests never run on GitHub — manual-only** (`workflow_dispatch`, no push or tag trigger):
  GitHub Actions test runs are kept off to avoid CI cost, so run the full test matrix **locally**
  before every push. **Only the docs workflow runs on GitHub** (it builds + deploys on push to
  `master`). Trigger `CI.yml` by hand from the Actions tab if a clean-environment resolve is ever
  wanted.
