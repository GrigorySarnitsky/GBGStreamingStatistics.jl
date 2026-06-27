# Bandwidth

Bandwidth selection: the kernel- and degree-aware normal-reference rule (the default; Silverman the
Gaussian/degree-0 case), the data-driven Sheather–Jones (`rule = :sj` classic, `:isj` improved-via-DCT) and fit-based AIC (`rule = :aic`)
selectors (leave-one-out CV a recorded binned-ties gap), and the response-aware local plug-in
``h^\star(x)`` (§4.5).

```@autodocs
Modules = [GBGStreamingStatistics.Bandwidth]
```
