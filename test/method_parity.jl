# Method parity — the method set is frozen, so a change cannot silently add or drop a method.
#
# The package DECLARES method parity to GBGCodingTools (RULE.TEST.METHOD_PARITY): the engine is the
# tool's `TestHarness.method_parity`, and the frozen numbers live in `test/method_parity_baseline.toml`,
# each moved number refrozen there with its reason, in the same commit as the change that moved it.
#
# Two tables, counting different denominators (see the baseline's header): the public surface per
# operation, and the derived per-module count over the package's whole function surface, internals
# included — ownership (`m.module === mod`) over every function the package touches, so a method defined
# under a qualified name cannot hide. The engine also asserts `detect_ambiguities` over the package is
# empty. Extension methods (the HDF5 weak dependency) live in their own module and are not counted.

using Test
using GBGStreamingStatistics
using Statistics
using TOML
import GBGCodingTools
include("public_surface.jl")

const METHOD_BASELINE = TOML.parsefile(joinpath(@__DIR__, "method_parity_baseline.toml"))

# The external modules whose generics the package extends: `Base.merge!`, and `Statistics.mean` / `var` /
# `std`.
GBGCodingTools.TestHarness.method_parity(GBGStreamingStatistics, PACKAGE_MODULES, METHOD_BASELINE;
                                         externals = (Base, Statistics))
