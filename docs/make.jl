using Pkg
Pkg.activate(@__DIR__)
using Documenter
using GBGStreamingStatistics


# Doc lint: an inline-math span ``…`` must not be split across a source line by a CommonMark
# list/blockquote leader. A wrapped continuation beginning with "- "/"* "/"+ "/"> "/"1." is parsed as a
# new list item, dropping the math out of math mode — and Documenter's build does NOT catch this (the
# page renders, just wrongly). This guard does; keep inline math on one source line.
function lint_inline_math()
  leader = r"^\s*([-*+>]|\d+[.)])\s"
  twobacktick = r"(?<!`)``(?!`)"
  sources = String[]
  for f in readdir(joinpath(@__DIR__, "..", "src"); join = true)
    endswith(f, ".jl") && push!(sources, f)
  end
  for name in ("ARCHITECTURE.md", "FORMALIZATION.md", "TESTING.md", "PERFORMANCE.md")
    push!(sources, joinpath(@__DIR__, name))
  end
  generated = Set(["architecture.md", "formalization.md", "testing.md", "performance.md"])
  for (root, _, files) in walkdir(joinpath(@__DIR__, "src")), f in files
    endswith(f, ".md") && !(f in generated) && push!(sources, joinpath(root, f))
  end
  problems = String[]
  for file in sources
    isjl = endswith(file, ".jl")
    indoc = !isjl                 # .md renders entirely; .jl renders only inside """ … """ docstrings
    infence = false
    inspan = false                # inside an unclosed ``…`` inline-math span
    for (n, line) in enumerate(eachline(file))
      if isjl && isodd(count("\"\"\"", line))
        indoc = !indoc; inspan = false; continue
      end
      indoc || continue
      if startswith(strip(line), "```")
        infence = !infence; inspan = false; continue
      end
      infence && continue
      if inspan && occursin(leader, line)
        push!(problems, string(relpath(file, joinpath(@__DIR__, "..")), ":", n,
                               "  inline ``math`` broken across a line by a list/blockquote leader"))
      end
      isodd(count(twobacktick, line)) && (inspan = !inspan)
    end
  end
  isempty(problems) || error("doc lint failed — keep inline math on one line:\n  " * join(problems, "\n  "))
  @info "doc lint: inline math OK across $(length(sources)) source files"
end

lint_inline_math()


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
               GBGStreamingStatistics.Kernels,
               GBGStreamingStatistics.LocalPolynomial,
               GBGStreamingStatistics.MarginalDensity,
               GBGStreamingStatistics.Regression,
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
      "Overview"               => "api.md",
      "GBGStreamingStatistics" => "api/GBGStreamingStatistics.md",
      "Accumulators"           => "api/Accumulators.md",
      "Histograms"             => "api/Histograms.md",
      "Kernels"                => "api/Kernels.md",
      "MarginalDensity"        => "api/MarginalDensity.md",
      "Regression"             => "api/Regression.md",
      "ConditionalDensity"     => "api/ConditionalDensity.md",
      "Bandwidth"              => "api/Bandwidth.md",
    ],
    "Design & Internals" => Any[
      "Architecture"  => "architecture.md",
      "Formalization" => "formalization.md",
      "Testing"       => "testing.md",
      "Performance"   => "performance.md",
      "Internals"     => "internals.md",
    ],
  ]
)

deploydocs(
  repo      = "github.com/GrigorySarnitsky/GBGStreamingStatistics.jl",
  devbranch = "master",
)
