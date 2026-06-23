using Pkg
Pkg.activate(@__DIR__)
using Documenter
using GBGStreamingStatistics

# Render the design-rationale docs into the site without moving them: ARCHITECTURE.md,
# FORMALIZATION.md, TESTING.md, and PERFORMANCE.md stay the single sources (referenced by CLAUDE.md),
# copied into src/ at build time. The copies (src/architecture.md, …) are gitignored.
function copy_rationale_pages()
  for name in ("ARCHITECTURE.md", "FORMALIZATION.md", "TESTING.md", "PERFORMANCE.md")
    cp(joinpath(@__DIR__, name), joinpath(@__DIR__, "src", lowercase(name)); force = true)
  end
end

copy_rationale_pages()

# Make `using GBGStreamingStatistics` available to any `jldoctest` examples; `doctest = true` runs and
# verifies them on every build so they cannot drift from the code. Add them incrementally.
DocMeta.setdocmeta!(GBGStreamingStatistics, :DocTestSetup, :(using GBGStreamingStatistics); recursive = true)

makedocs(
  sitename  = "GBGStreamingStatistics.jl",
  modules   = [GBGStreamingStatistics,
               GBGStreamingStatistics.Accumulators,
               GBGStreamingStatistics.Histograms,
               GBGStreamingStatistics.DensityEstimation,
               GBGStreamingStatistics.Regression,
               GBGStreamingStatistics.ConditionalDensity,
               GBGStreamingStatistics.MarkovInference,
               GBGStreamingStatistics.Bandwidth],
  format    = Documenter.HTML(
    canonical     = "https://grigorysarnitsky.github.io/GBGStreamingStatistics.jl/stable",
    collapselevel = 1,
  ),
  checkdocs = :exports,
  doctest   = true,
  pages     = [
    "Home" => "index.md",
    "Getting started" => "tutorial.md",
    "API"  => Any[
      "Overview"               => "api.md",
      "GBGStreamingStatistics" => "api/GBGStreamingStatistics.md",
      "Accumulators"           => "api/Accumulators.md",
      "Histograms"             => "api/Histograms.md",
      "DensityEstimation"      => "api/DensityEstimation.md",
      "Regression"             => "api/Regression.md",
      "ConditionalDensity"     => "api/ConditionalDensity.md",
      "MarkovInference"        => "api/MarkovInference.md",
      "Bandwidth"              => "api/Bandwidth.md",
    ],
    "Design & Internals" => Any[
      "Architecture"  => "architecture.md",
      "Formalization" => "formalization.md",
      "Testing"       => "testing.md",
      "Performance"   => "performance.md",
    ],
  ]
)

deploydocs(
  repo      = "github.com/GrigorySarnitsky/GBGStreamingStatistics.jl",
  devbranch = "master",
)
