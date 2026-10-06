Status: open — approved by the owner on 2026-10-06; sent to the GBGCodingTools session the same day.

# To GBGCodingTools: reword `RULE.API.POSITIONAL` to the dispatch criterion Julia practises

The rule today: "In an operation, an argument that can change the result is positional, defaults included;
an argument that cannot change the result — a performance hint, a builder's cost/accuracy tunable — is a
keyword." SS's whole exported surface takes its smoothing controls and modifiers as keywords
(`density(summary, x; h, kernel, degree, adaptive, bounds, start, normalize, binning)`), and would need a
`Switch` on every estimator. Before asking for one, the owner had the criterion researched. No recognized
source states it, and the sources that address the question say otherwise:

- **The Julia manual** gives no positional-vs-keyword rule. Its stated reasons for keywords are "a large
  number of arguments, or a large number of behaviors", readability, and "any subset … in any order"
  (Functions chapter); its own example, `plot(x, y, width=2)`, makes a result-changing option a keyword.
  The one hard fact is that keywords do not take part in dispatch (Methods chapter).
- **Base and the stdlib** make result-changing options keywords throughout: `sort(v; by, rev)`,
  `round(x; digits)`, `var(x; corrected)`, `quantile(v, p; alpha, beta)`, `isapprox(x, y; atol, rtol)`. What
  they make positional is a strategy type that is dispatched on: `round(x, RoundingMode)`,
  `cholesky(A, RowMaximum())`, `qr(A, ColumnNorm())` (the `Val`/Bool forms were deprecated for these).
- **SciML's style guide**: "Arguments that do not have defaults should preferably be made into positional
  arguments", and "prefer instances" for the algorithm, `solve(prob, Tsit5())`.
- **The ecosystem**: `kde(data; bandwidth, kernel)` (KernelDensity.jl), `loess(xs, ys; span, degree)`,
  `fit(Histogram, v; nbins, closed)` (StatsBase), `solve(prob, alg; abstol, reltol)` (DifferentialEquations
  — tolerances change the result and are keywords; only the algorithm is positional), GLM's distribution and
  link positional as dispatch types, `weights` keyword.
- **Elsewhere**: PEP 570 ("use keyword-only when names have meaning"), scikit-learn (every hyperparameter a
  keyword with a default), the Swift and Kotlin guidelines, and the "boolean trap" all choose by call-site
  readability, not by whether the argument changes the result.

The defensible core of the rule is the dispatch argument (Rackauckas & Nie, arXiv:1807.06430, on the
DifferentialEquations common interface): an argument whose type selects the method or algorithm is
positional, so that it dispatches and third parties can extend it. Against making tunables positional:
optional positional arguments cannot be skipped (setting `degree` forces `h` and `kernel`), a call of
several unlabelled values is misread and misordered, and n optional positionals generate n + 1 methods,
each a candidate for dispatch ambiguity.

Asked: reword the rule to that core — for instance, "An argument whose type selects the method or algorithm
is positional, so it takes part in dispatch; options and tunables, whether or not they change the result,
are keywords. An exception is a `Switch` on this rule, with its reason." The final wording is the tool's
and the owner's. SS's surface then complies without a switch, and the four `RULE.API.POSITIONAL` switches
GBGPeriodicFields carries may become unnecessary (that is PF's to judge).

The research and its sources (URLs, verbatim quotes, and what could not be verified) can be sent on request.
