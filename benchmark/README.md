# Benchmarks

The performance suite. The single performance-critical path is the **per-sample accumulation hot
loop** (`observe!`); post-hoc smoothing is second-order (`docs/PERFORMANCE.md`).

!!! note
    Scaffold stage — not yet populated. The first benchmark class to add is the accumulation hot loop
    (the scatter-add binning + the Pébay moment recursion), measured on CPU and GPU, once
    `Accumulators` lands.
