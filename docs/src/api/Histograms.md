# Histograms

The summary types: [`Histogram1D`](@ref) and [`Histogram2D`](@ref) — finely binned empirical measures that fold the stream into a fixed-size form — plus [`bin_index`](@ref) and [`observe_pairs!`](@ref) for two-time statistics.

## Summaries

```@docs
Histogram1D
Histogram2D
```

## Accumulation

`observe!` / `merge!` are the universal [`Accumulators`](Accumulators.md) interface; `observe_pairs!`
folds two-time pairs into a `Histogram2D`, and `bin_index` exposes the binning map.

```@docs
observe_pairs!
bin_index
```

The recovered marginal moments (`mean` / `var` / `std` / `central_moment`) read off these power-sums live
in [`Expectations`](Expectations.md), as the ``g = x^p`` case of the power-sum quadrature.
