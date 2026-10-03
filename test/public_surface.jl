# The package's public surface, as every completeness guard reads it.
#
# The three declared harness contracts (RULE.TEST.ALLOC_ZERO, ORACLE_INDEPENDENT, METHOD_PARITY) each
# need the same thing — the set of operations a user can call — and each derives it from OWNERSHIP
# through GBGCodingTools' `TestHarness.package_operations`: a function is on the surface when a package
# module owns at least one of its methods. That reaches the `Statistics` generics this package extends
# (`mean` / `var` / `std`) exactly as it reaches a generic the package declares itself, which a scan of
# `names(GBGStreamingStatistics)` would not. The module list is the one input the derivation needs.

const PACKAGE_MODULES = [
  GBGStreamingStatistics.Accumulators,
  GBGStreamingStatistics.Histograms,
  GBGStreamingStatistics.Kernels,
  GBGStreamingStatistics.Expectations,
  GBGStreamingStatistics.LocalPolynomial,
  GBGStreamingStatistics.MarginalDensity,
  GBGStreamingStatistics.ConditionalExpectations,
  GBGStreamingStatistics.ConditionalDensity,
  GBGStreamingStatistics.Bandwidth,
]
