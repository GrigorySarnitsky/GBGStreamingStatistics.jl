# Testing

How the tests are built and the per-operation coverage contract. The statistics analog of PF's
analytic round-trips: estimators are checked against distributions with closed forms.

## Layers

- **Analytic oracles.** Estimators are verified against distributions whose target quantities are known
  exactly; the closed forms live once in `test/testoracles.jl` (the `TestOracles` module — `gaussian_pdf`,
  `exponential_pdf`, `lognormal_pdf`, `gamma_pdf`, `beta_pdf`, `student_t_pdf`, `uniform_pdf`, the
  Ornstein–Uhlenbeck `ou_drift` / `ou_diffusion` / `ou_transition`, and the `trapz` / `ise` helpers),
  `include`d by every file that needs one. The marginal density `density` is checked across a
  **representative analytic-PDF suite** — Gaussian, uniform (bounded support), exponential, Student-``t``
  ``\nu = 3`` (heavy tail), log-normal (skew), and a bimodal mixture (`test/density_convergence.jl`) — and
  across its full knob grid (KDE / local-log `degree`, fixed / Abramson `adaptive`, the `bounds` boundary
  correction, the asymmetric `gamma_kernel` / `beta_kernel`, the parametric `start`, and `normalize`) in
  `test/densityestimation.jl`. The conditional / transition estimators use **Gaussian** marginals
  (``f_X``, conditional moments) and an **Ornstein–Uhlenbeck** process, whose drift, diffusion, and
  transition density ``p_\tau(y \mid x)`` are closed-form — the end-to-end oracle for the conditional /
  transition densities (§4.3) and the Kramers–Moyal functionals (§4.4).
- **Bandwidth selectors.** Each `bandwidth` rule is checked against its analytic target
  (`test/bandwidth.jl`): `:normal` reproduces Silverman ``(4/3)^{1/5}`` and the Epanechnikov
  ``2.214\times`` ratio; `:sj` / `:isj` ``\approx`` the normal reference on Gaussian data and tighten on a
  bimodal mixture; `:aic` finds a stable interior optimum; and the asymmetric-kernel Chen `:plugin` lowers
  the ISE on skewed (log-normal / gamma) data and stays neutral on the near-exponential boundary case.
- **Consistency / rate tests.** The estimation error falls at its theoretical rate as ``n`` grows. For
  the **density** this is implemented (`test/density_convergence.jl`): the ISE falls at ``n^{-2p/(2p+1)}``
  — ``n^{-4/5}`` for degree 0/1, the steeper ``n^{-8/9}`` for degree 2 (the higher-order ``n^{-1/9}``
  bandwidth), the Abramson variable bandwidth sharing the degree-0 order. The conditional-mean
  ``n^{-2/5}`` and conditional-density ``n^{-2/3}`` rate tests remain forthcoming.
- **Mergeability tests.** For every summary, *accumulate-then-merge ``=`` accumulate-together* (to
  round-off), plus the incompatible-summary rejection — the monoid laws of invariant 3
  (`test/mergeability.jl` for `Histogram1D`, `test/histograms2d.jl` for `Histogram2D`). Every post-hoc
  estimator is then separately re-checked to be a deterministic functional of the merged summary (the
  `density` start, the `bandwidth` selectors, the conditional moments, the conditional / transition
  densities, drift / diffusion).
- **Allocation audit.** `test/allocations.jl` (`Histogram1D`) and `test/histograms2d.jl` (`Histogram2D`)
  assert via `@allocated` that `observe!` / `observe_pairs!` / `merge!` allocate only a **bounded** amount
  on a preallocated summary — the KernelAbstractions launch, never an ``O(n)`` or ``O(\text{bins})``
  temporary — the GPU-hot-loop contract of invariant 4.
- **GPU-vs-CPU parity.** The gated suite (`test/gpu/`, run on a machine with a device) exercises
  `observe!` / `merge!` / `observe_pairs!` on device arrays under `allowscalar(false)` and compares against
  the CPU result — counts exactly, power-sums to round-off, including heavy atomic contention. Verified on
  the AMD RX 6700 XT (ROCm); the vendor stacks (CUDA / AMDGPU) stay test-only dependencies.
- **Package-quality checks.** `test/aqua.jl` (`Aqua.test_all`: ambiguities, piracy, stale deps, compat
  bounds, undefined exports) and `test/jet.jl` (JET static analysis of every public operation — including
  the asymmetric-kernel, parametric-start, `normalize`, and bandwidth-selector paths, whose function-typed
  fields are the type-stability risk) run against the loaded package.

## Coverage contract

For a new or changed operation, the *Definition of done* in `CLAUDE.md` is the spec: an independent
oracle correctness test, a rate test, a mergeability test, the `@allocated` assertion for in-place ops,
GPU parity, and a `FORMALIZATION.md` entry with its test link. The matrix tracks what each operation has.

## Coverage matrix

| Operation | Oracle | Rate | Merge | Alloc | GPU |
|---|---|---|---|---|---|
| `Histogram1D` (`observe!` / `merge!`) | Gaussian moments ``M_2``–``M_6`` | — | yes | yes | yes |
| `Histogram2D` (`observe_pairs!` / `merge!`) | indep-Gaussian marginals + covariance | — | yes | yes | yes |
| `density` (all knobs) | 6-PDF suite + Gaussian | yes, ``n^{-2p/(2p+1)}`` | yes | n/a | via summary |
| `density_standard_error` | exact SE + band coverage | — | n/a | n/a | via summary |
| `bandwidth` (`:normal` / `:sj` / `:isj` / `:aic`) | Silverman / ``\delta(K)`` / bimodal | — | yes | n/a | via summary |
| `bandwidth` (`:plugin`, asymmetric) | ISE gain on log-normal / gamma / beta | — | yes | n/a | via summary |
| `conditional_mean` / `conditional_variance` | bivariate-Gaussian | forthcoming ``n^{-2/5}`` | yes | n/a | via summary |
| `conditional_density` / `transition_density` | bivariate-Gaussian + OU | forthcoming ``n^{-2/3}`` | yes | n/a | via summary |
| `drift` / `diffusion` / `increment_moment` | Ornstein–Uhlenbeck | — | yes | n/a | via summary |

`n/a` is not applicable (a rate test for a summary; an allocation audit for a read-only post-hoc
estimator); **forthcoming** is a recorded gap. "via summary" means the estimator is a host-side functional
of the summary, so its device coverage is the summary's GPU-parity row.

## Layout

`test/runtests.jl` includes each module's tests in its own wrapper module so top-level helpers and
constants do not collide; every file is also runnable on its own (`julia --project=. test/<file>.jl`).
The per-module files are `accumulators.jl`, `histograms.jl`, `histograms2d.jl`, `densityestimation.jl`,
`density_convergence.jl`, `bandwidth.jl`, `regression.jl`, `conditionaldensity.jl`, `markovinference.jl`;
the cross-cutting files are `oracles.jl` (the Gaussian / OU end-to-end oracles), `mergeability.jl`,
`allocations.jl`, and `testoracles.jl` (the shared `TestOracles` analytic library). `aqua.jl` and `jet.jl`
are the package-quality checks; `gpu/` is the gated device suite. The test environment is
`test/Project.toml` (Aqua, JET, Random, SpecialFunctions, Statistics, Test).

!!! note "Status"
    The summary suites (`Histogram1D`, `Histogram2D`), `density` (with the asymmetric kernels, the
    parametric `start`, and `normalize`), `bandwidth` (the normal-reference, Sheather–Jones `:sj` / `:isj`,
    `:aic`, and the asymmetric Chen `:plugin`), the §4.2 regression, the §4.3
    `conditional_density` / `transition_density`, and the §4.4 `drift` / `diffusion` / `increment_moment`
    are real — the Gaussian and analytic-PDF oracles, the bivariate-Gaussian and **Ornstein–Uhlenbeck**
    oracles, mergeability, the bounded `@allocated` audit, and **CPU-vs-AMD-GPU parity** (`test/gpu/`,
    ROCm); `Aqua.test_all` and JET run for real throughout. The response-aware local bandwidth and the
    Pawula / Markov diagnostics remain placeholders until those land; the conditional-mean ``n^{-2/5}`` and
    conditional-density ``n^{-2/3}`` rate tests are forthcoming (the **density** rate test is implemented —
    `test/density_convergence.jl`).
