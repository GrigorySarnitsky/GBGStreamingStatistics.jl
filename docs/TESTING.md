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
  across its full option grid (KDE / local-log `degree`, fixed / Abramson `adaptive`, the `bounds` boundary
  correction, the asymmetric `gamma_kernel` / `beta_kernel`, the parametric `start`, and `normalize`) in
  `test/marginaldensity.jl`. The conditional / transition estimators use **Gaussian** marginals
  (``f_X``, conditional moments) and an **Ornstein–Uhlenbeck** process, whose transition density
  ``p_\tau(y \mid x)`` is closed-form — the end-to-end oracle for the conditional / transition densities
  (§4.3).
- **Bandwidth selectors.** Each `bandwidth` rule is checked against its analytic target
  (`test/bandwidth.jl`): `:normal` reproduces Silverman ``(4/3)^{1/5}`` and the Epanechnikov
  ``2.214\times`` ratio; `:sj` / `:isj` ``\approx`` the normal reference on Gaussian data and tighten on a
  bimodal mixture; `:aic` finds a stable interior optimum; and the asymmetric-kernel Chen `:plugin` lowers
  the ISE on skewed (log-normal / gamma) data and stays neutral on the near-exponential boundary case. On
  the **Marron–Wand claw** (the field's standard hard-multimodal benchmark — a broad Gaussian with five
  sharp spikes) the data-driven `:sj` / `:isj` resolve the spikes the global-``\hat\sigma`` normal
  reference oversmooths into one bump (``\approx`` an order of magnitude lower ISE).
- **Structural invariants (property-based, oracle-free).** Beyond the mergeability monoid laws, `density`
  satisfies **affine equivariance** exactly (`test/marginaldensity.jl`): for ``Y = aX + b`` with the
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
  converges to the exact un-binned KDE as the grid refines — at ``O(\Delta^2)`` in the bin width ``\Delta``
  for the default linear/centroid `binning` (mass at the recovered within-bin centroid, via the stored
  ``S_1``) and at ``O(\Delta)`` for `binning = :simple` (mass at the bin center), the two slopes checked
  against ``2`` and ``1``.
  The conditional-mean ``n^{-2/5}`` and conditional-density ``n^{-2/3}`` rates are likewise implemented
  (`test/conditional_convergence.jl`): a curved-mean (``m(x) = \sin x``) local-linear fit at
  ``h \propto n^{-1/5}`` gives ISE ``\propto n^{-4/5}`` (slope ``\approx -0.80``), and the
  bivariate-Gaussian double-kernel at ``h \propto n^{-1/6}`` gives the slower **2-D** rate ISE
  ``\propto n^{-2/3}`` (slope ``\approx -2/3``, measured in ``y`` at ``x_0 = 0``, the n-range starting
  higher than the 1-D test since the 2-D asymptotic regime is reached later — below ``\sim 10^5`` the
  measured slope runs shallow, the conditional analog of the density test's non-Gaussian muddying).
- **Mergeability tests.** For every summary, *accumulate-then-merge ``=`` accumulate-together* (to
  round-off), plus the incompatible-summary rejection — the monoid laws of invariant 3
  (`test/mergeability.jl` for `Histogram1D`, `test/histograms2d.jl` for `Histogram2D`). The **cone laws**
  that extend the monoid (`scale!` homogeneity and distributivity, the weighted `merge!` as
  scale-then-merge, ``\texttt{weight} = 1`` ``=`` plain `merge!`, the Float-count guard) are tested
  alongside the hyperstatistics oracles in `test/hyperstatistics.jl`. Every post-hoc
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
  bounds, undefined exports) and `test/jet.jl` (JET static analysis of every public estimator path — including
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
5. **Argument-domain errors** — a non-positive bandwidth, and **every rejected option combination**, throws.
   The guards are audited for completeness: each `throw(ArgumentError(…))` — the bad `degree` / `binning` /
   `method` / `rule` symbol, the mutually-exclusive `start` / `adaptive` / `bounds` / asymmetric-kernel
   combinations, the `degree ≥ 2` and asymmetric-kernel standard-error gaps — carries a `@test_throws`.

Characteristics 2–4 are the `check_point_estimator` battery; 1 and 5 are per-estimator (the oracle, the
`@test_throws`).

**Not swept: the valid option space.** Every *option* is tested against an oracle (vary one, hold the rest
at default), and the genuinely-**interacting** pairs are tested together (`adaptive`+`bounds` rewrite the
same variance; `degree`×`kernel`) — but the full Cartesian product of valid option combinations is
deliberately **not** exhausted. The options are largely orthogonal (different code paths), the exotic
combinations have no analytic oracle (only a weak "finite" check), and ``N``-way coverage is a combinatorial
explosion with diminishing returns past the pairwise level. The systematic-and-exhaustive effort goes into
the *invalid* combinations (finite, high-value — characteristic 5), not the valid ones.

**Estimand-specific** — layered on where the estimand demands it: **normalization** ``\int\hat f = 1`` (the
densities), the **`transition_density` alias**, and the **convergence rate** (the matrix).

**Class variations.** The **summaries** (`Histogram1D` / `Histogram2D`) are accumulators, not query-point
estimators: their contract is the moment-recovery oracle, mergeability, the `@allocated` audit, and GPU
parity. The **marginal selector** `bandwidth(::Histogram1D)` has no query point: oracle (per rule),
mergeability, and argument errors — no scalar/vector or off-support. The **local selector**
`bandwidth(::Histogram2D, x)` (``h^\star``) *is* a query-point estimator, so it carries the universal
contract (scalar ``=`` vector, merged-summary functional, off-support `NaN`; `test/localbandwidth.jl`).

**Standard errors** carry their own contract, `check_standard_error` (`test/testcontracts.jl`), on top of
the universal one: every `*_standard_error` is checked for **positivity**, **scalar ``=`` vector**, and the
**``n^{-1/2}`` scaling** — merging ``K`` copies of a summary (``K\times`` the samples, the same shape)
shrinks the standard error by *exactly* ``\sqrt K`` (``\mathrm{SE}(K\cdot s) = \mathrm{SE}(s)/\sqrt K``),
the deterministic form of ``\mathrm{SE}\propto n^{-1/2}`` and a direct test of merged-histogram error
behaviour. For the *kernel* standard errors (`density`, `conditional_expectation`) it also checks the
**``h^{-1/2}`` bandwidth scaling** (`check_se_bandwidth_scaling`, *asymptotic* — on a flat density at small
``h``, ``\mathrm{SE}(2h) \approx \mathrm{SE}(h)/\sqrt2``, since changing ``h`` also moves ``\hat f`` and the
``-\hat f^2/n`` term is ``h``-independent). It is *also* run through `check_point_estimator` (so the error, like the estimate, is a
deterministic functional of the merged summary — `SE`(merge) ``=`` `SE`(accumulate-together)), and carries a
**calibration** check per family: the predicted standard error matches the empirical standard deviation of
the estimate across independent datasets (Monte-Carlo — `density` and `conditional_mean` at the density
peak, `expectation` of ``\cos``). A `*_with_standard_error` additionally asserts its
`(value, standard_error)` equals the separate estimate / standard-error calls exactly.

**Merged histograms — two or more.** Mergeability is checked at the query-point-functional level
(`check_point_estimator`, with `merged` a **3-way** merge of disjoint thirds for `density` — associativity
beyond two) *and* at the summary level (`test/mergeability.jl`, `test/histograms2d.jl` — the monoid / cone
laws), for both the estimates and their standard errors.

## Coverage matrix

The query-point estimators — the universal contract (oracle / scalar==vector / merge / off-support /
argument errors) plus the convergence rate:

| Estimator | Oracle | Scalar==vec | Merge | Off-support | Arg err | Rate |
|---|---|---|---|---|---|---|
| `density` | 6-PDF suite + Gaussian | yes | yes | n/a, defined everywhere | yes | yes, ``n^{-2p/(2p+1)}`` |
| `density_standard_error` / `density_with_standard_error` | two-term variance vs Monte-Carlo spread + band coverage; tracks `adaptive` / `bounds` / `start`; the paired estimate-plus-error consistency | yes | yes | n/a | yes | n/a, variability band |
| `conditional_mean` / `conditional_variance` | bivariate-Gaussian (+ heteroscedastic curved-mean for Fan–Yao `conditional_variance`) | yes | yes | `NaN` | yes | yes, ``n^{-4/5}`` ISE |
| `conditional_expectation` / `conditional_pushforward_density` | biv-Gaussian ``\mathbb{E}[\cos Y\mid x]`` + lognormal pushforward; the `conditional_mean` fold + degree-0 `conditional_density` coincidence | yes | yes | `NaN` | yes | inherits `conditional_mean` ``n^{-4/5}`` |
| `conditional_expectation_standard_error` / `conditional_expectation_with_standard_error` | leverage × noise-scale = Monte-Carlo sampling spread (biv-Gaussian) + leading-order closed form; the paired value+error matches the separate calls | yes | yes | `NaN` | yes | n/a, a standard error |
| `conditional_density` / `transition_density` | biv-Gaussian + OU + bounded-response (asym. `kernely`) | yes | yes | `NaN` off-support ``x`` | yes | yes, ``n^{-2/3}`` |
| `conditional_density_standard_error` / `conditional_density_with_standard_error` | two-term closed form (biv-Gaussian) = Monte-Carlo spread; ``n^{-1/2}`` + ``h_2^{-1/2}`` scaling; paired consistency + `transition_density` aliases | yes | yes | `NaN` | yes | n/a, a standard error |
| `bandwidth(::Histogram2D)` ``h^\star(x)`` | quadratic-mean Gaussian (closed-form) | yes | yes | `NaN` | yes | n/a, a selector |
| `expectation` / `pushforward_density` / `expectation_standard_error` / `expectation_with_standard_error` | Gaussian MGF/trig + polynomial exactness; lognormal & non-monotone ``\chi^2_1`` pushforward; SE ``= s/\sqrt n``, scaling / mergeability / Monte-Carlo, and the paired value+error matches the separate calls | yes (`pushforward`) | yes | n/a, mass in range | yes | n/a, exact for degree ``\le`` order else ``O(\Delta^{\texttt{order}+1})`` |

The summaries and the selector carry a different contract:

| Operation | Oracle | Merge | Alloc | GPU |
|---|---|---|---|---|
| `Histogram1D` (`observe!` / `merge!`) | Gaussian moments ``M_2``–``M_6`` | yes | yes | yes |
| `Histogram2D` (`observe_pairs!` / `merge!`) | indep-Gaussian marginals + covariance | yes | yes | yes |
| cone (`scale!` / weighted `merge!`) | χ²-superstatistics → Student-``t``, Gaussian mixture | cone laws | yes | yes |
| `bandwidth(::Histogram1D)` (`:normal` / `:sj` / `:isj` / `:aic` / `:plugin`) | Silverman / ``\delta(K)`` / bimodal / ISE gain | yes | n/a | via summary |
| `LocalPolynomial.local_fit` (the jet primitive) | exact jet (value / slope / curvature) on a polynomial mean | n/a | n/a | via summary |

`n/a` is not applicable. The **rate column is now uniform** — `density` (`test/density_convergence.jl`),
`conditional_mean` (``n^{-4/5}`` ISE), and `conditional_density` (the slower 2-D ``n^{-2/3}``) each carry a
log–log ISE-slope test (the conditional pair in `test/conditional_convergence.jl`); `conditional_variance`
inherits `conditional_mean`'s rate (a second-stage local fit), so it needs no separate test. The power-sum
quadrature readouts (`expectation`, `pushforward_density`) carry a **polynomial-exactness** test (they equal
the raw empirical mean ``\tfrac1n\sum_i g(X_i)`` to ``10^{-10}`` for ``g`` of degree ``\le \texttt{order}``)
and a binning-error bound in place of an asymptotic slope. "via summary" means a host-side functional of the
summary, so the device coverage is the summary's row.

**Tolerances are seed-robust.** The ISE-comparison thresholds (the `:plugin` gains, the parametric-start
gains, the convergence-rate slopes) are not cherry-picked to a fixed seed: they are audited across many
seeds, and the threshold sits above the worst case with margin. The audit is the test-suite analog of the
Monte-Carlo–averaged MISE the source papers report — its goal is a non-fragile regression guard, not an
MISE estimate. (It is why the `:plugin`-on-`Beta(2,2)` test asserts only *competitiveness*: the rule of
thumb is already near-optimal on that smooth, symmetric case, so the gain there is real but marginal.)

## Layout

`test/runtests.jl` includes each module's tests in its own wrapper module so top-level helpers and
constants do not collide; every file is also runnable on its own (`julia --project=. test/<file>.jl`).
The per-module files are `accumulators.jl`, `histograms.jl`, `histograms2d.jl`, `kernels.jl`,
`localpolynomial.jl`, `marginaldensity.jl`, `density_convergence.jl`, `bandwidth.jl`, `localbandwidth.jl`,
`conditionalexpectations.jl`, `conditionaldensity.jl`, `conditional_convergence.jl`, `expectations.jl`;
the cross-cutting files are `oracles.jl` (the Gaussian / OU end-to-end oracles), `mergeability.jl`,
`hyperstatistics.jl` (the cone laws and the superstatistics / forgetting oracles), `allocations.jl`,
`testoracles.jl` (the shared `TestOracles` analytic library — the closed forms), and
`testcontracts.jl` (the `TestContracts` universal-contract battery — `check_point_estimator`). `aqua.jl`
and `jet.jl` are the package-quality checks; `gpu/` is the gated device suite. The test environment is
`test/Project.toml` (Aqua, JET, Random, SpecialFunctions, Statistics, Test, plus the `LoweredCodeUtils`
compat pin JET needs).

!!! note "Status"
    The summary suites (`Histogram1D`, `Histogram2D`), `density` (with the asymmetric kernels, the
    parametric `start`, and `normalize`), `bandwidth` (the normal-reference, Sheather–Jones `:sj` / `:isj`,
    `:aic`, the asymmetric Chen `:plugin`, and the joint-summary response-aware local `:local`
    ``h^\star(x)``), the shared `LocalPolynomial.local_fit` jet primitive, the §4.2 regression (with the
    `degree` / `kernel` options and the Fan–Yao `conditional_variance`), the §4.3
    `conditional_density` / `transition_density`
    are real — the Gaussian and analytic-PDF oracles, the bivariate-Gaussian, **heteroscedastic
    curved-mean**, quadratic-mean, and **Ornstein–Uhlenbeck** oracles, mergeability, the uniform
    `check_point_estimator` contract, `density`'s exact affine equivariance, the Marron–Wand claw stress,
    the bounded `@allocated` audit, and **CPU-vs-AMD-GPU parity** (`test/gpu/`, ROCm); `Aqua.test_all` and
    JET run for real throughout. The Kramers–Moyal / Pawula / Markov diagnostics are archived (see `archive/`), pending a consumer. The
    conditional-mean ``n^{-2/5}`` and conditional-density ``n^{-2/3}`` rate tests are now implemented
    (`test/conditional_convergence.jl`), alongside the **density** rate test (`test/density_convergence.jl`)
    and the direct kernel-constant oracle (`test/kernels.jl`) — closing the rate column of the matrix.
