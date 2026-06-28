# Accumulators

The accumulation interface: the `AbstractAccumulator` supertype and the two operations every summary supports — `observe!` (fold a batch of samples into the summary, on CPU or GPU) and `merge!` (combine two summaries by addition). The concrete summary types live in the [`Histograms`](Histograms.md) module.

## The accumulator type

```@docs
AbstractAccumulator
```

## Operations

```@docs
observe!
merge!
```
