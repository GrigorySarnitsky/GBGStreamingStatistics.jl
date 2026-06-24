# Performance

The performance model, best practices, and the benchmark coverage contract.

## The model — one critical path

The package has a single performance-critical path: the **per-sample accumulation hot loop**
([`observe!`](@ref GBGStreamingStatistics.Accumulators.observe!)), which runs over every sample of the stream — for the turbulence consumer,
``N^3 \times`` steps. Everything else is second-order:

- **Post-hoc smoothing is dwarfed by accumulation.** Every smoother here — the binned KDE (§4.1) and
  the local-polynomial sums (§4.2) — is a convolution of the small binned summary (a few ``\times
  10^3`` bins) with the kernel, computed once per request. The FFT does it in ``O(B \log B)`` for
  ``B`` bins and **amortizes over bandwidth**: transform the bin counts once, then each ``h`` (and
  the bandwidth grid the local ``h^\star(x)`` sweep needs) is a multiply by the kernel's transform
  and an inverse FFT. This total work is smaller than the accumulation by roughly ``10^3`` — so the
  dominant lever is the accumulation loop, **not** the smoother.

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

## Benchmark coverage contract

For a perf-relevant change, add a benchmark in the matching `benchmark/` class and run it before and
after against a saved baseline; confirm on the GPU, since platform optima diverge.

!!! note "Status"
    The benchmark suite is not yet built. The accumulation hot loop (`observe!`) is now implemented
    and is the first thing to benchmark, alongside the forthcoming GPU-resident kernel.
