# Test suite entry point for GBGStreamingStatistics.
#
# Each module's tests live in their own file, included into its own wrapper module so that top-level
# test helpers and constants do not collide across files. Every file is also runnable on its own:
#   julia --project=. test/<module>.jl

module AccumulatorsTests;       include("accumulators.jl");       end
module HistogramsTests;         include("histograms.jl");         end
module Histogram2DTests;        include("histograms2d.jl");       end
module LocalPolynomialTests;    include("localpolynomial.jl");    end
module DensityEstimationTests;  include("densityestimation.jl");  end
module DensityConvergenceTests; include("density_convergence.jl"); end
module RegressionTests;         include("regression.jl");         end
module ConditionalDensityTests; include("conditionaldensity.jl"); end
module MarkovInferenceTests;    include("markovinference.jl");    end
module BandwidthTests;          include("bandwidth.jl");          end
module OraclesTests;            include("oracles.jl");            end
module MergeabilityTests;       include("mergeability.jl");       end
module HyperstatisticsTests;    include("hyperstatistics.jl");    end
module AllocationsTests;        include("allocations.jl");        end
module AquaTests;               include("aqua.jl");               end
module JETTests;                include("jet.jl");                end
