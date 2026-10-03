# Histograms

The summary types: [`Histogram1D`](@ref) and [`Histogram2D`](@ref) — finely binned empirical measures that fold the stream into a fixed-size form — plus [`observe_pairs!`](@ref) for two-time statistics.

## Summaries

```@docs
Histogram1D
Histogram2D
```

## Accumulation

The verbs `observe!` / `merge!` are the universal [`Accumulators`](Accumulators.md) interface; `observe_pairs!`
folds two-time pairs into a `Histogram2D`.

```@docs
observe_pairs!
```

## Persistence (HDF5)

With the `HDF5` package loaded, a summary reads from and writes to an HDF5 **group** of a file you open and
own — so it sits beside your other datasets, one file with many groups. This is provided by a package
extension, which keeps `HDF5` out of the core dependency tree. Write with [`write_histogram`](@ref) (or the
equivalent `parent[name] = summary` sugar); read with [`read_histogram`](@ref). The reload is exact and
still merges.

```@docs
write_histogram
read_histogram
```
