using Pkg
Pkg.activate(@__DIR__)
using Documenter
using GBGCodingTools
using GBGStreamingStatistics
const ROOT = normpath(joinpath(@__DIR__, ".."))
const TOOLING = GBGCodingTools.Render.pages(GBGCodingTools.ConfigIO.read_config(ROOT), ROOT)


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

GBGCodingTools.makedocs(ROOT;
  sitename  = "GBGStreamingStatistics.jl",
  modules   = [GBGStreamingStatistics,
               GBGStreamingStatistics.Accumulators,
               GBGStreamingStatistics.Histograms,
               GBGStreamingStatistics.Kernels,
               GBGStreamingStatistics.Expectations,
               GBGStreamingStatistics.LocalPolynomial,
               GBGStreamingStatistics.MarginalDensity,
               GBGStreamingStatistics.ConditionalExpectations,
               GBGStreamingStatistics.ConditionalDensity,
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
      "Overview"                => "api.md",
      "GBGStreamingStatistics"  => "api/GBGStreamingStatistics.md",
      "Accumulators"            => "api/Accumulators.md",
      "Histograms"              => "api/Histograms.md",
      "Kernels"                 => "api/Kernels.md",
      "Expectations"            => "api/Expectations.md",
      "MarginalDensity"         => "api/MarginalDensity.md",
      "ConditionalExpectations" => "api/ConditionalExpectations.md",
      "ConditionalDensity"      => "api/ConditionalDensity.md",
      "Bandwidth"               => "api/Bandwidth.md",
    ],
    "Design & Internals" => Any[
      "Notation"      => "notation.md",
      "Architecture"  => "architecture.md",
      "Formalization" => "formalization.md",
      "Testing"       => "testing.md",
      "Performance"   => "performance.md",
      "Internals"     => "internals.md",
      "Documentation tooling" => only(TOOLING),
    ],
  ]
)

deploydocs(
  repo      = "github.com/GrigorySarnitsky/GBGStreamingStatistics.jl",
  devbranch = "master",
)
