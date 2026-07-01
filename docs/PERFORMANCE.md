# Performance

The performance model, best practices, and the benchmark coverage contract.

## The model — one critical path

The package has a single performance-critical path: the **per-sample accumulation hot loop**
([`observe!`](@ref GBGStreamingStatistics.Accumulators.observe!)), which runs over every sample of the stream — for the turbulence consumer,
``N^3 \times`` steps. Everything else is second-order:

- **Post-hoc smoothing is dwarfed by accumulation.** Every smoother here — the binned KDE (§4.1) and
  the local-polynomial sums (§4.2) — is a convolution of the small binned summary (a few ``\times 10^3``
  bins) with the kernel, computed once per request by **direct summation** over the bins (``O(B)`` per
  query point for a compact-reach kernel, ``O(B^2)`` for the distance-histogram curvature functionals of
  `:sj`). Even so it is dwarfed by the accumulation over every sample of the stream — a few thousand bins
  versus ``N^3`` samples. The **FFT amortization** — transform the bin counts once in ``O(B \log B)``, then
  each ``h`` (and the bandwidth grid a local ``h^\star(x)`` sweep needs) is a multiply by the kernel's
  transform and an inverse FFT — is **deferred (recorded, not built)**: a pure optimization over the exact
  direct summation, the major lever only once a bandwidth *grid* is swept (`FFTW` is already in place,
  used today for the ISJ DCT).

## Hot-loop best practices

- **Contention, not FP64, is the enemy.** The conditioning variable is peaked (Gaussian-ish), so the
  modal bins are atomic hot spots. Storing power-sums *about the bin center* makes each update a plain
  atomic add (counts integer-exact, power-sums `Float64` to round-off), so global atomics handle even
  heavy contention directly — the proven PF `anterpolate!` pattern, verified on the AMD GPU.
  **Privatization** (per-block local histograms merged via the monoid) stays available as a future
  optimization if a channel's modal-bin contention dominates.
- **Wider-than-data precision is free.** The histogram is small (kilobytes per channel, versus a
  gigabyte field), so `Float64` moment accumulators cost almost nothing, and on the GPU the `Float64`
  arithmetic is hidden behind the memory- and atomic-bound accumulation, not FP64 throughput
  (invariant 6).
- **Allocation-free.** `observe!` allocates only a bounded amount (the KernelAbstractions launch) on a
  preallocated summary (invariant 4) — never an ``O(n)`` temporary that would escape the host GC and
  wedge a long GPU run.
- **The cone operations are off the hot path.** `merge!`, the weighted `merge!`, and `scale!` are
  ``O(B)`` in-place elementwise broadcasts over the ``B`` bins — per combination, not per sample — so the
  hyperstatistical mixture (one weighted `merge!` per parameter value) and exponential forgetting (one
  `scale!` per block) cost a bin-sweep each, negligible beside the per-sample accumulation, and
  allocation-free like `merge!`. The `count = Float64` variant routes `observe!` through the GPU
  atomic-*float* add — the same memory-/atomic-bound profile as the integer-count path.

## Benchmark coverage contract

For a perf-relevant change, add a benchmark in the matching `benchmark/` class and run it before and
after against a saved baseline; confirm on the GPU, since platform optima diverge.

!!! note "Status"
    The benchmark suite is not yet built. The accumulation hot loop (`observe!`) is implemented and runs
    on both CPU and the GPU (the same KernelAbstractions scatter kernel, parity-verified on the AMD RX
    6700 XT) — it is the first thing to benchmark.
