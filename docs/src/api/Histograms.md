# Histograms

The summary types: `Histogram1D` and `Histogram2D` — finely binned empirical measures that fold the stream into a fixed-size form — plus `bin_index` and `observe_pairs!` for two-time statistics.

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

## Recovered moments

The global marginal moments, recovered post hoc from the power-sums by the shift formula (exact),
extending `Statistics`.

```@docs
mean(::Histogram1D)
var(::Histogram1D)
std(::Histogram1D)
central_moment
```
