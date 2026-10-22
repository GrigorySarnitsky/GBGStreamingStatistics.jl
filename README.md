# GBGStreamingStatistics.jl

This package is a backend-agnostic library for **online nonparametric statistics**. It estimates marginal densities, conditional means, variances and expectations, conditional and two-time transition densities, and expectations and pushforward densities of a transformed variable — each with its standard error, and with data-driven bandwidth selection — from a **stream of samples that is never stored** — each sample is observed once, folded into a fixed-size summary, and discarded.

The summary (a fine histogram with per-bin moments) is accumulated online and is *estimand-agnostic*: every reported quantity is a post-hoc smoothing functional of that one summary, so there is no storage of the underlying state and no recompute-from-raw reruns. It is a **mergeable monoid**, in fact a cone (combine across threads, devices, or time-chunks by addition, or weight and superpose them), persists to HDF5, and runs on both **CPU and GPU** array backends, on a single shared-memory node. The library is **domain-neutral and model-neutral**: it consumes a stream of values (optionally lagged pairs) and returns estimators.

**When it fits.** Reach for it when the data is a stream and storing or re-reading it is the bottleneck — inside a time-marching solver (fold each step's samples into the summary and discard, with no trajectory storage and no second pass), when one accumulation pass must serve many post-hoc views (every bandwidth, estimand, and conditioning value), for parallel or chunked accumulation (mergeable summaries combine across threads, GPUs, and time windows), and for memory-tight or GPU-resident pipelines. It suits, for example, turbulence and molecular dynamics — but nothing in the library is specific to a domain.

The `GBG` part in the name `GBGStreamingStatistics.jl` comes from the city of Göteborg (Gothenburg, Sweden), where the primary development took place.

## Installation

Install it by name from the GBG packages' own registry,
[GBGRegistry](https://github.com/GrigorySarnitsky/GBGRegistry), added after `General`:

```julia
using Pkg
Pkg.Registry.add("General")
Pkg.Registry.add(url = "https://github.com/GrigorySarnitsky/GBGRegistry")
Pkg.add("GBGStreamingStatistics")
```

or at the package prompt:

```
pkg> registry add General https://github.com/GrigorySarnitsky/GBGRegistry
pkg> add GBGStreamingStatistics
```

Then you can use it as usual:

```julia
using GBGStreamingStatistics
```

## Documentation

Documentation with the full API reference is available at [https://grigorysarnitsky.github.io/GBGStreamingStatistics.jl/](https://grigorysarnitsky.github.io/GBGStreamingStatistics.jl/).
