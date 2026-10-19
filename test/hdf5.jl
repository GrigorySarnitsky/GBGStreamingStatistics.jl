# HDF5 persistence (the weak-dependency extension): summaries write to / read from HDF5 *groups* in a file
# the consumer owns, exactly (fields, element types, shapes), and the reload still merges. Loading HDF5
# here triggers the extension `GBGStreamingStatisticsHDF5Ext`.
using Test
using GBGStreamingStatistics
using HDF5
using Random

@testset "HDF5 persistence (extension)" begin
  mktempdir() do dir
    path = joinpath(dir, "summaries.h5")

    Random.seed!(1)
    h1 = Histogram1D(-8.0, 8.0, 200; order = 4); observe!(h1, randn(50_000))
    hc = Histogram1D(-8.0, 8.0, 200; order = 3, count_type = Float64)   # the Float-count cone variant
    observe!(hc, randn(30_000)); scale!(hc, 0.5)
    xs = randn(40_000)
    h2 = Histogram2D(-6.0, 6.0, 60, -6.0, 6.0, 60; lag = 4)
    observe_pairs!(h2, xs, 0.6 .* xs .+ 0.8 .* randn(40_000))

    # several summaries into ONE file, in groups, beside a plain dataset — the consumer owns the file
    h5open(path, "w") do f
      f["note"] = "hello"                  # the user's own data
      group = create_group(f, "stats")
      write_histogram(group, "marg", h1)       # verb form, into a group they made
      group["cone"] = hc                       # the setindex! sugar
      write_histogram(f, "joint", h2)      # directly under the file root
    end

    @testset "exact round-trip: fields, element types, shapes" begin
      h5open(path, "r") do f
        r1 = read_histogram(f, "stats/marg")        # by path
        rc = read_histogram(f["stats"], "cone")     # by group handle
        r2 = read_histogram(f, "joint")

        @test r1 isa Histogram1D
        @test r1.lower == h1.lower && r1.upper == h1.upper
        @test eltype(r1.counts) == UInt64
        @test r1.counts == h1.counts
        @test r1.powersums == h1.powersums

        @test eltype(rc.counts) == Float64          # the cone variant's count eltype is preserved
        @test rc.counts == hc.counts
        @test rc.powersums == hc.powersums
        @test size(rc.powersums, 1) == 3            # order preserved

        @test r2 isa Histogram2D
        @test r2.n_x == h2.n_x && r2.n_y == h2.n_y && r2.lag == 4
        @test r2.lower_x == h2.lower_x && r2.upper_y == h2.upper_y
        @test r2.counts == h2.counts
        @test r2.powersums == h2.powersums
      end
    end

    @testset "functional round-trip: estimators agree after reload" begin
      r1 = h5open(path, "r") do f; read_histogram(f, "stats/marg"); end
      @test mean(r1) == mean(h1)
      @test density(r1, 0.3; h = 0.4) == density(h1, 0.3; h = 0.4)
    end

    @testset "reload is bit-exact and still mergeable" begin
      Random.seed!(2); d = randn(80_000)
      summary_a = Histogram1D(-8.0, 8.0, 200); observe!(summary_a, view(d, 1:40_000))
      summary_b = Histogram1D(-8.0, 8.0, 200); observe!(summary_b, view(d, 40_001:80_000))
      p = joinpath(dir, "half.h5")
      h5open(p, "w") do f; write_histogram(f, "a", summary_a); end
      ra = h5open(p, "r") do f; read_histogram(f, "a"); end
      @test ra.counts == summary_a.counts && ra.powersums == summary_a.powersums   # reload reproduces it bit-for-bit
      # merging into the reload matches merging into the original (reload preserves the mergeable structure;
      # `==` here, not `≈`, because ra is bit-identical to a, so the two merges share operands and order)
      ma = merge!(ra, summary_b); mo = merge!(deepcopy(summary_a), summary_b)
      @test ma.counts == mo.counts && ma.powersums == mo.powersums
    end

    @testset "self-describing: attributes readable as plain HDF5" begin
      h5open(path, "r") do f
        group = f["stats/marg"]
        @test attrs(group)["type"] == "GBGStreamingStatistics.Histogram1D"
        @test attrs(group)["version"] == 1
        @test attrs(group)["lower"] == -8.0
        @test Set(keys(group)) == Set(["counts", "powersums"])
      end
    end

    @testset "error paths" begin
      h5open(path, "r") do f
        @test_throws ArgumentError read_histogram(f, "note")     # a dataset, not a group
        @test_throws ArgumentError read_histogram(f, "absent")   # no such object
        @test_throws ArgumentError read_histogram(f, "stats")    # a group, but no histogram type tag
      end
    end
  end
end
