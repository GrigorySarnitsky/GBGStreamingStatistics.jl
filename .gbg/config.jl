# .gbg/config.jl — this package's settings, as data.
#
# Read by `ConfigIO`, which parses this file and evaluates nothing: only literals and the constructors
# of `Model` are accepted. Every field not written here takes its default; write a line only where the
# default is wrong for this package.
Config(
    package = "GBGStreamingStatistics",
    tool = v"2.0.0-rc.57",
    docs_build = ["julia", "--project=docs", "docs/make.jl"],
    test = ["julia", "--project", "-e", "using Pkg; Pkg.test()"],
    # The one tooling page, written by `Render.pages` during the docs build. No rule reads a generated page.
    generated = ["docs/src/gbg.md"],
    # The docs build copies each rationale document into `docs/src/` as a lowercase sibling, and the
    # copy is what the site renders. A fix belongs in the source; the change guard refuses one landing
    # in a copy, which the next build would overwrite.
    copies = [
        "docs/src/architecture.md" => "docs/ARCHITECTURE.md",
        "docs/src/formalization.md" => "docs/FORMALIZATION.md",
        "docs/src/testing.md" => "docs/TESTING.md",
        "docs/src/performance.md" => "docs/PERFORMANCE.md",
    ],
)
