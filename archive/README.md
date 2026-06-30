# archive/

Working code removed from the active package but kept retrievable. Not compiled, not tested, not
documented — inert until deliberately restored.

## `IncrementMoments.jl` + `incrementmoments_tests.jl` (the §4.4 SDE-reading functionals)

(The test file is named `incrementmoments_tests.jl` here, not `incrementmoments.jl`, to avoid a
case-insensitive-filesystem collision with `IncrementMoments.jl`; on restore it becomes
`test/incrementmoments.jl` again.)

The conditional increment moments and their Kramers–Moyal functionals — `increment_moment`,
`kramers_moyal`, `drift`, `diffusion`, and the declared-but-unimplemented `pawula_quotient` /
`markov_test` — archived to focus the package on its universal statistical core (the summary-as-measure
and its two readout families: kernel-weighting and power-sum quadrature).

These are thin functionals of `conditional_mean` (e.g. `drift = (m̂(x) − x)/τ`), the most physics-flavored
and least-complete part (the Pawula/Markov diagnostics were stubs; `kramers_moyal` was capped at `n ≤ 2`,
the higher moments needing a `Histogram2D` order-raise). This is the last working version — `drift` /
`diffusion` verified against the Ornstein–Uhlenbeck oracle.

**To restore** (when a consumer needs them, and the `Histogram2D` order-raise is in place for `n ≥ 3`):
move both files back to `src/` and `test/`, then re-wire `src/GBGStreamingStatistics.jl` (the `include`,
the `using .IncrementMoments: …`, and the six exports), `test/runtests.jl`, `test/jet.jl`,
`test/oracles.jl` (the OU `drift`/`diffusion` assertions — `ou_drift` / `ou_diffusion` still live in
`test/testoracles.jl`), `docs/make.jl`, `docs/src/index.md`, and a fresh `docs/src/api/IncrementMoments.md`.
