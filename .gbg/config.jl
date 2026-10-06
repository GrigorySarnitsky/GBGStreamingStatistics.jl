# .gbg/config.jl — this package's settings, as data.
#
# Read by `ConfigIO`, which parses this file and evaluates nothing: only literals and the constructors
# of `Model` are accepted. Every field not written here takes its default; write a line only where the
# default is wrong for this package.
Config(
    package = "GBGStreamingStatistics",
    tool = v"2.0.0-rc.59",
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
    # The one address this package publishes: its author line and every commit's author field.
    leak_emails = ["grigory@sarnitsky.se"],
    # The fields whose experts read this package's documentation. A statistician reads every estimand,
    # estimator, standard error and rate claim — the whole content; an architect reads the two-half
    # layering, the invariants and the decision log. No physics and no Fourier scheme (invariant 2), so
    # the turbulence, microphysics and spectral disciplines are deliberately absent.
    disciplines = [:statistical_estimation, :software_architecture],
    # The leading identity of every docstring is the estimand: this package is all statistics and no
    # physics, so its counterpart of "the physical meaning first" is the population quantity, in standard
    # probability notation, stated before the estimator (CLAUDE.md, "estimand first, estimator second").
    settings = [
        Setting(rule = Symbol("RULE.DOC.LEADING_IDENTITY"), name = :leading_identity, value = "estimand"),
    ],
    switches = [
        Switch(target = Symbol("RULE.DOC.COORD_FREE"),
               reason = "No structured space to state a map intrinsically on: every variate is a scalar \
                         on an interval, and the package carries no vector or tensor notation at all \
                         (invariant 2; the \\mathbf and \\partial_i house guards are deliberately absent \
                         here). The rule's differential-geometry reader would have nothing to judge, \
                         and without the switch every review reports that nobody judges it. The owner, \
                         2026-10-03."),
    ],
)
