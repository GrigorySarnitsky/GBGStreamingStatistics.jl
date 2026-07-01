# Accumulators

The accumulation interface: the `AbstractAccumulator` supertype and the operations every summary supports — `observe!` (fold a batch of samples into the summary, on CPU or GPU), `merge!` (combine two summaries by addition, or — weighted — the affine combination), and `scale!` (reweight a summary by a nonnegative scalar). With `merge!` and `scale!` the summaries form a convex **cone**, the algebra behind hyperstatistics and exponential forgetting. The concrete summary types live in the [`Histograms`](Histograms.md) module.

## The accumulator type

```@docs
AbstractAccumulator
```

## Operations

`observe!` and the plain `merge!` are the additive monoid; the weighted `merge!(into, from; weight)` and `scale!` extend it to the cone (the hyperstatistical parameter-mixture density and exponential forgetting).

```@docs
observe!
merge!
scale!
```
