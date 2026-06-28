# Testing

How the tests are built and the per-operation coverage contract. The statistics analog of PF's
analytic round-trips: estimators are checked against distributions with closed forms.

## Layers

- **Analytic oracles.** Estimators are verified against distributions whose target quantities are known
  exactly; the closed forms live once in `test/testoracles.jl` (the `TestOracles` module — `gaussian_pdf`,
  `exponential_pdf`, `lognormal_pdf`, `gamma_pdf`, `beta_pdf`, `student_t_pdf`, `uniform_pdf`, the
  Marron–Wand `claw_pdf`, the Ornstein–Uhlenbeck `ou_stationary_std` / `ou_drift` / `ou_diffusion` /
  `ou_transition`, and the `trapz` / `ise` helpers), `include`d by every file that needs one. The marginal density `density` is checked across a
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
  the ISE on skewed (log-normal / gamma) data and stays neutral on the near-exponential boundary case. On
  the **Marron–Wand claw** (the field's standard hard-multimodal benchmark — a broad Gaussian with five
  sharp spikes) the data-driven `:sj` / `:isj` resolve the spikes the global-``\hat\sigma`` normal
  reference oversmooths into one bump (``\approx`` an order of magnitude lower ISE).
- **Structural invariants (property-based, oracle-free).** Beyond the mergeability monoid laws, `density`
  satisfies **affine equivariance** exactly (`test/densityestimation.jl`): for ``Y = aX + b`` with the
  range and bandwidth scaled alike, ``\hat f_Y(ax+b) = \hat f_X(x)/a`` to floating point — a guard against
  a stray normalization / ``\sigma``-scaling factor a single-scale oracle would miss, and a guarantee of
  unit-consistency (the estimate cannot depend on ``\mu``m vs m).
- **Consistency / rate tests.** The estimation error falls at its theoretical rate as ``n`` grows. For
  the **density** this is implemented (`test/density_convergence.jl`): the ISE falls at ``n^{-2p/(2p+1)}``
  — ``n^{-4/5}`` for degree 0/1, the steeper ``n^{-8/9}`` for degree 2 (the higher-order ``n^{-1/9}``
  bandwidth), the Abramson variable bandwidth sharing the degree-0 order, and the **Epanechnikov**
  reproducing the (kernel-independent) ``n^{-4/5}``. The rate is measured on the **Gaussian**, whose fast
  asymptotics make it visible by ``n = 10^6``; a non-Gaussian does not reach the asymptotic regime there
  (its slope muddies to ``\approx -0.66`` and degree 2 fails to separate from degree 0), a recorded
  omission rather than a missing test. A companion **binning-error** test confirms the binned estimate
  converges to the exact un-binned KDE as the grid refines, at ``O(\Delta)`` in the bin width ``\Delta``
  (simple binning; linear binning would be ``O(\Delta^2)``, a refinement latent in the stored ``S_1``).
  The conditional-mean ``n^{-2/5}`` and conditional-density ``n^{-2/3}`` rate tests are scheduled.
- **Mergeability tests.** For every summary, *accumulate-then-merge ``=`` accumulate-together* (to
  round-off), plus the incompatible-summary rejection — the monoid laws of invariant 3
  (`test/mergeability.jl` for `Histogram1D`, `test/histograms2d.jl` for `Histogram2D`). Every post-hoc
  estimator is then re-checked to be a deterministic functional of the merged summary — uniformly, through
  the `check_point_estimator` battery (see *The estimator contract*) for the query-point estimators, and a
  dedicated check for the `bandwidth` selectors and the `density` start.
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

## The estimator contract

Every estimator is held to the **same uniform set** of characteristics, in two tiers plus a few
class-specific variations. The *Definition of done* in `CLAUDE.md` is the parent spec; this section is the
operational form. The universal tier is enforced **by construction** — `check_point_estimator`
(`test/testcontracts.jl`) runs it from one call per estimator, so a new estimator cannot silently omit a
check (the failure mode that once left `conditional_mean` without a mergeability test).

**Universal** — every query-point estimator (a functional of the summary at a query point) satisfies:

1. **Analytic-oracle correctness** — recovers a closed-form target (the oracle layer above).
2. **Scalar ``=`` vector** — the vectorized call equals the per-point scalar calls.
3. **Functional of the merged summary** — `f`(merge) ``=`` `f`(accumulate-together) — mergeability.
4. **Off-support sentinel** — `NaN` where the estimand is undefined off the data's support (the conditional
   estimators); a marginal density is defined everywhere and is exempt.
5. **Argument-domain errors** — a non-positive bandwidth, and any rejected knob combination, throws.

Characteristics 2–4 are the `check_point_estimator` battery; 1 and 5 are per-estimator (the oracle, the
`@test_throws`).

**Estimand-specific** — layered on where the estimand demands it: **normalization** ``\int\hat f = 1`` (the
densities), the **Kramers–Moyal identities** (``M_1 = \mathbb{E}[Y \mid x] - x``, the variance
decomposition, the ``1/(n!\tau)`` scaling — `drift` / `diffusion` / `increment_moment`), the
**`transition_density` alias**, and the **convergence rate** (the matrix).

**Class variations.** The **summaries** (`Histogram1D` / `Histogram2D`) are accumulators, not query-point
estimators: their contract is the moment-recovery oracle, mergeability, the `@allocated` audit, and GPU
parity. The **selector** `bandwidth` has no query point: oracle (per rule), mergeability, and argument
errors — no scalar/vector or off-support.

## Coverage matrix

The query-point estimators — the universal contract (oracle / scalar==vector / merge / off-support /
argument errors) plus the convergence rate:

| Estimator | Oracle | Scalar==vec | Merge | Off-support | Arg err | Rate |
|---|---|---|---|---|---|---|
| `density` | 6-PDF suite + Gaussian | yes | yes | n/a, defined everywhere | yes | yes, ``n^{-2p/(2p+1)}`` |
| `density_standard_error` | exact SE + band coverage | yes | yes | n/a | yes | n/a, variability band |
| `conditional_mean` / `conditional_variance` | bivariate-Gaussian | yes | yes | `NaN` | yes | **scheduled** ``n^{-2/5}`` |
| `conditional_density` / `transition_density` | biv-Gaussian + OU | yes | yes | `NaN` off-support ``x`` | yes | **scheduled** ``n^{-2/3}`` |
| `drift` / `diffusion` / `increment_moment` | Ornstein–Uhlenbeck | yes | yes | `NaN` | yes | **scheduled**, ``= conditional_mean`` |

The summaries and the selector carry a different contract:

| Operation | Oracle | Merge | Alloc | GPU |
|---|---|---|---|---|
| `Histogram1D` (`observe!` / `merge!`) | Gaussian moments ``M_2``–``M_6`` | yes | yes | yes |
| `Histogram2D` (`observe_pairs!` / `merge!`) | indep-Gaussian marginals + covariance | yes | yes | yes |
| `bandwidth` (`:normal` / `:sj` / `:isj` / `:aic` / `:plugin`) | Silverman / ``\delta(K)`` / bimodal / ISE gain | yes | n/a | via summary |

`n/a` is not applicable; **scheduled** is a recorded gap with a planned test (the roadmap in `CLAUDE.md`).
The **rate is the one characteristic not yet uniform** — `density` has it, the conditional estimators have
it scheduled, and `drift` / `diffusion` inherit `conditional_mean`'s rate (they are linear functionals of
it). "via summary" means a host-side functional of the summary, so the device coverage is the summary's
row.

**Tolerances are seed-robust.** The ISE-comparison thresholds (the `:plugin` gains, the parametric-start
gains, the convergence-rate slopes) are not cherry-picked to a fixed seed: they are audited across many
seeds, and the threshold sits above the worst case with margin. The audit is the test-suite analog of the
Monte-Carlo–averaged MISE the source papers report — its goal is a non-fragile regression guard, not an
MISE estimate. (It is why the `:plugin`-on-`Beta(2,2)` test asserts only *competitiveness*: the rule of
thumb is already near-optimal on that smooth, symmetric case, so the gain there is real but marginal.)

## Layout

`test/runtests.jl` includes each module's tests in its own wrapper module so top-level helpers and
constants do not collide; every file is also runnable on its own (`julia --project=. test/<file>.jl`).
The per-module files are `accumulators.jl`, `histograms.jl`, `histograms2d.jl`, `densityestimation.jl`,
`density_convergence.jl`, `bandwidth.jl`, `regression.jl`, `conditionaldensity.jl`, `markovinference.jl`;
the cross-cutting files are `oracles.jl` (the Gaussian / OU end-to-end oracles), `mergeability.jl`,
`allocations.jl`, `testoracles.jl` (the shared `TestOracles` analytic library — the closed forms), and
`testcontracts.jl` (the `TestContracts` universal-contract battery — `check_point_estimator`). `aqua.jl`
and `jet.jl` are the package-quality checks; `gpu/` is the gated device suite. The test environment is
`test/Project.toml` (Aqua, JET, Random, SpecialFunctions, Statistics, Test, plus the `LoweredCodeUtils`
compat pin JET needs).

!!! note "Status"
    The summary suites (`Histogram1D`, `Histogram2D`), `density` (with the asymmetric kernels, the
    parametric `start`, and `normalize`), `bandwidth` (the normal-reference, Sheather–Jones `:sj` / `:isj`,
    `:aic`, and the asymmetric Chen `:plugin`), the §4.2 regression, the §4.3
    `conditional_density` / `transition_density`, and the §4.4 `drift` / `diffusion` / `increment_moment`
    are real — the Gaussian and analytic-PDF oracles, the bivariate-Gaussian and **Ornstein–Uhlenbeck**
    oracles, mergeability, the uniform `check_point_estimator` contract, `density`'s exact affine
    equivariance, the Marron–Wand claw stress, the bounded `@allocated` audit, and **CPU-vs-AMD-GPU
    parity** (`test/gpu/`, ROCm); `Aqua.test_all` and JET run for real throughout. The response-aware local
    bandwidth and the Pawula / Markov diagnostics remain placeholders until those land; the conditional-mean
    ``n^{-2/5}`` and conditional-density ``n^{-2/3}`` rate tests are **scheduled** (`CLAUDE.md` roadmap; the
    **density** rate test is implemented — `test/density_convergence.jl`).
