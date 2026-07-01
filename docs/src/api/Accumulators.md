# Accumulators

The accumulation interface: the [`AbstractAccumulator`](@ref) supertype and the operations every summary supports — [`observe!`](@ref) (fold a batch of samples into the summary, on CPU or GPU), [`merge!`](@ref) (combine two summaries by addition, or — weighted — the nonnegative **conic** combination), and [`scale!`](@ref) (reweight a summary by a nonnegative scalar). The concrete summary types live in the [`Histograms`](Histograms.md) module.

## The accumulator type

```@docs
AbstractAccumulator
```

## Operations

The additive monoid (`observe!`, `merge!`) and its conic extension (`scale!` and weighted `merge!`) — the algebra behind hyperstatistics and exponential forgetting.

```@docs
observe!
merge!
scale!
```
