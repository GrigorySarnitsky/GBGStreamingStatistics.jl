# GBGStreamingStatistics.jl

This package is a backend-agnostic library for **online nonparametric statistics**. It estimates marginal densities, conditional expectations, conditional and two-time transition densities, and standard functionals of these (Kramers–Moyal coefficients, drift and diffusion, Pawula diagnostics), from a **stream of samples that is never stored** — each sample is observed once, folded into a fixed-size summary, and discarded.

The library is domain-neutral and model-neutral: it consumes a stream of values (optionally lagged pairs) and returns estimators, with no concept of velocity, droplet, trajectory, or turbulence, and no assumption that the data is stationary, Markov, diffusive, or parametric. The summary is accumulated online and response-agnostic; all smoothing (kernel, bandwidth, channel) is applied post hoc — so there is no storage of the underlying state and no recompute-from-raw reruns. The implementation is designed for both CPU and GPU array backends through KernelAbstractions, on a single shared-memory node.

It is the statistics sibling of [`GBGPeriodicFields.jl`](https://github.com/GrigorySarnitsky/GBGPeriodicFields.jl): that package is all mathematics and no physics; this one is all statistics and no physics.

The `GBG` part in the name `GBGStreamingStatistics.jl` comes from the city of Göteborg (Gothenburg, Sweden), where the primary development took place.

## Installation

This is an unregistered package, so it is not listed in the `General` Julia registry. To install it, run:

```julia
]add https://github.com/GrigorySarnitsky/GBGStreamingStatistics.jl
```

or

```julia
using Pkg

Pkg.add(url="https://github.com/GrigorySarnitsky/GBGStreamingStatistics.jl")
```

Then you can use it as usual:

```julia
using GBGStreamingStatistics
```

## Documentation

Documentation with the full API reference is available at [https://grigorysarnitsky.github.io/GBGStreamingStatistics.jl/](https://grigorysarnitsky.github.io/GBGStreamingStatistics.jl/).
