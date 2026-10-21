# Changelog

## 0.1.0 — 2026-10-22

The first release. A stream of samples is folded, once and never stored, into a fixed-size, mergeable
summary — the binned empirical measure `Histogram1D` or the joint `Histogram2D`, accumulated by one
atomic scatter on the CPU or a GPU — and every estimate is read off that summary afterwards: the
marginal density, the conditional mean, variance and expectation, the conditional and transition
densities, expectations and pushforward densities, and the bandwidth selectors, each estimator with its
standard error. Summaries merge, scale and persist to HDF5.
