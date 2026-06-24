# GPU tests (opt-in, gated)

These verify that the *one portable kernel* behind `observe!` / `merge!` / `observe_pairs!` runs on a
real GPU and produces results identical to the CPU — counts exactly (integer atomics are
order-independent), power-sums to round-off (float atomic-adds accumulate in a nondeterministic order)
— under `allowscalar(false)`, which makes any stray scalar index of a device array throw.

They are **not** part of `Pkg.test()` (which has no GPU); they run only on a machine with a device.

## ROCm / AMDGPU

```sh
julia --project=test/gpu/rocm -e 'using Pkg; Pkg.instantiate()'   # one-time
julia --project=test/gpu/rocm test/gpu/rocm/runtests.jl
```

The shared, backend-parameterized suite lives in `test/gpu/suite.jl`; the entry script
`rocm/runtests.jl` defines `to_device`, `device_functional`, and `GPU_BACKEND`, and calls
`AMDGPU.allowscalar(false)` before including it. A CUDA backend would add a parallel
`test/gpu/cuda/` with the same three definitions. The sub-environment `Manifest.toml` is gitignored
and resolves fresh per machine.
