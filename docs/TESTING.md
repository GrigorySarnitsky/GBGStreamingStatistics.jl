# Testing

How the tests are built and the per-operation coverage contract. The statistics analog of PF's
analytic round-trips: estimators are checked against distributions with closed forms.

## Layers

- **Analytic oracles.** Estimators are verified against distributions whose target quantities are
  known exactly: **Gaussian** marginals (``f_X``, conditional moments) and an **Ornstein–Uhlenbeck**
  process, whose drift, diffusion, and transition density ``p_\tau(y\mid x)`` are closed-form — the
  natural end-to-end oracle for the conditional/transition densities (§4.3) and the Kramers–Moyal
  functionals (§4.4).
- **Consistency / rate tests.** The estimation error falls at its theoretical rate as ``n`` grows
  (``n^{-2/5}`` for a density or conditional mean; the slower ``n^{-2/3}`` for a conditional
  density), checked over a sequence of growing ``n``.
- **Mergeability tests.** For every carrier, *accumulate-then-merge `=` accumulate-together* (to
  round-off), and `merge!` is associative — the monoid laws of invariant 3.
- **Allocation audit.** `test/allocations.jl` asserts via `@allocated` (on the CPU, where a
  materialized array's bytes are visible) that `observe!` / `observe_pairs!` allocate nothing on a
  preallocated summary — the GPU-hot-loop contract of invariant 4.
- **GPU-vs-CPU parity.** Once a device path exists, the public surface is exercised on device arrays
  under `allowscalar(false)` and compared against the CPU result; the vendor stacks (CUDA/AMDGPU)
  stay test-only dependencies.
- **Package-quality checks.** `test/aqua.jl` (`Aqua.test_all`: ambiguities, piracy, stale deps,
  compat bounds, undefined exports) and `test/jet.jl` (JET static analysis) run against the loaded
  package.

## Coverage contract

For a new or changed operation, the *Definition of done* in `CLAUDE.md` is the spec: an independent
oracle correctness test, a rate test, a mergeability test, the `@allocated` assertion for in-place
ops, GPU parity, and a `FORMALIZATION.md` entry with its test link.

## Layout

`test/runtests.jl` includes each module's tests in its own wrapper module so top-level helpers and
constants do not collide. The cross-cutting files are `test/oracles.jl` (the Gaussian/OU oracles and
the rate tests), `test/mergeability.jl`, and `test/allocations.jl`.

!!! note "Status"
    The `Histogram1D` foundation, `density`, and `bandwidth` suites are real — the Gaussian moment
    oracle, the Gaussian density oracle, mergeability (accumulate-then-merge == together), and the
    `@allocated == 0` audit; `Aqua.test_all` and JET run for real throughout. The OU / conditional and
    the ``n^{-2/5}`` rate tests remain placeholders until those estimators land.
