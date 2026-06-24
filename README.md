# GBGStreamingStatistics.jl

This package is a backend-agnostic library for **online nonparametric statistics**. It estimates marginal densities, conditional expectations, conditional and two-time transition densities, and standard functionals of these (Kramers–Moyal coefficients, drift and diffusion, Pawula diagnostics), from a **stream of samples that is never stored** — each sample is observed once, folded into a fixed-size summary, and discarded.

The summary (a fine histogram with per-bin moments) is accumulated online and is *response-agnostic*: every reported quantity is a post-hoc smoothing functional of that one summary, so there is no storage of the underlying state and no recompute-from-raw reruns. It is a **mergeable monoid** (combine across threads, devices, or time-chunks by addition) and runs on both **CPU and GPU** array backends, on a single shared-memory node. The library is **domain-neutral and model-neutral**: it consumes a stream of values (optionally lagged pairs) and returns estimators.

**When it fits.** Reach for it when the data is a stream and storing or re-reading it is the bottleneck — inside a time-marching solver (fold each step's samples into the summary and discard, with no trajectory storage and no second pass), when one accumulation pass must serve many post-hoc views (any bandwidth, channel, conditioning value, or lag), for parallel or chunked accumulation (mergeable summaries combine across threads, GPUs, and time windows), and for memory-tight or GPU-resident pipelines. It is used, for example, in turbulence and molecular dynamics — but nothing in the library is specific to a domain.

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
