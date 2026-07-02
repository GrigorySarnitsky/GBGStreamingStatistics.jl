# Histograms

The summary types: [`Histogram1D`](@ref) and [`Histogram2D`](@ref) — finely binned empirical measures that fold the stream into a fixed-size form — plus [`observe_pairs!`](@ref) for two-time statistics.

## Summaries

```@docs
Histogram1D
Histogram2D
```

## Accumulation

`observe!` / `merge!` are the universal [`Accumulators`](Accumulators.md) interface; `observe_pairs!`
folds two-time pairs into a `Histogram2D`.

```@docs
observe_pairs!
```
