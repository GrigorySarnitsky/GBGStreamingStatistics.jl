# Documentation Review — 2026-10-15-documentation-review-01

!!! info "Written by AI, not regularly reviewed"
    This page is written by AI and is not reviewed regularly by a human.

Run `2026-10-15-documentation-review-01` of `GBGStreamingStatistics`, by GBGCodingTools dev, from commit `eb8df5b6`.

It ran to the end. It launched 34 agent(s), with the bound readers, and made 1 of 1 fix pass(es). 48 line(s) went to the work queue, under the heading the end-of-turn hook does not count.

| Stage | Pass | Agents | Cache write | Cache read | Input | Output |
|---|---|---|---|---|---|---|
| find | 0 | 16 | 1162790 | 5901125 | 2304 | 36268 |
| verify | 1 | 18 | 1034154 | 4360252 | 1956 | 13226 |

Tokens, summed from each agent's transcript by `gbg run cost`.

## Not judged — 1

No reader of this package was given what these need. Nothing below says anything about them.

- RULE.DOC.ONE_NOTATION's `notation_page` names docs/src/notation.md, which is not a file here

## Findings — 100

### 2026-10-15-documentation-review-01-001  Three paired `*_with_standard_error` readers claim the estimate and its error come from one quadrature pass / one local-polynomial fit, but the code makes two passes (one for g, one for g^2) and, for the conditional reader, two column reductions and two fits plus a leverage solve.

Rule: `RULE.DOC.VERIFIED`.  
Reported by `full`, decided by an agent.

`expectation_with_standard_error` calls `expectation(summary, g)` and then `expectation(summary, x -> g(x)^2)`, each a full `_quadrature` sweep (only E[g] is shared, as the second sentence correctly says). `pushforward_density_with_standard_error` delegates to it, so it is also two passes per y. `conditional_expectation_with_standard_error` runs `_column_response(summary, g)` and `_column_response(summary, y -> g(y)^2)` (two y-quadrature responses, not 'the one'), then `_ce_value_and_se` does `local_fit` on sumG, `_leverage`, and `local_fit` on sumG2. By contrast `density_with_standard_error` ('the estimate is computed once and reused') and `conditional_density_with_standard_error` ('one x-profile reduction') state their sharing truthfully.

**Fix.** State the structure the code has: 'from two quadrature passes (E[g] and E[g^2]), the estimate computed once and shared' / 'from the two y-quadrature responses g and g^2, the fit's design shared', or restructure the code to one pass and keep the sentence.

- `src/Expectations.jl`:258  (expectation_with_standard_error)
- `src/Expectations.jl`:369  (pushforward_density_with_standard_error)
- `src/ConditionalExpectations.jl`:208  (conditional_expectation_with_standard_error)
- `src/ConditionalExpectations.jl`:213  (conditional_expectation_with_standard_error)
- pass 1: verifier `FIXED`  (commit `bb56595c`)

### 2026-10-15-documentation-review-01-002  The conditional pushforward at g = identity, degree 0 is said to 'coincide' with `conditional_density` (and its standard error with `conditional_density_standard_error`), but the two differ at O(Delta_y^2) in the value and only agree asymptotically in the standard error; the package's own test asserts agreement to rtol 0.02 / 0.03.

Rule: `RULE.DOC.VERIFIED`.  
Reported by `full`, decided by an agent.

`conditional_pushforward_density` goes through `_column_response`, which contracts the degree-2 interpolant of K_hz(g(y)-z) against the cell count, S_01 and S_02 (the within-cell moments enter), whereas `conditional_density._density_at` uses the cell-center value K_hy(c^y_k - y) alone. For the SE, `_cpf_value_and_se` uses the exact leverage sum Σ_j N_j w_j^2 (with w_j = K_h/D at degree 0) while `_cd_value_and_se` uses R(K_1)/(h_1 D), equal only when Σ_j N_j K_h^2 ≈ n f_X R(K_1)/h_1. The test `test/conditionalexpectations.jl` 'coincidence: pushforward(identity, degree = 0) ≈ conditional_density' notes 'they agree up to the within-Y-cell correction ... to O(Δ_y²)' and tests `≈ ... rtol = 0.02` / `rtol = 0.03`; the code comment above `_cpf_value_and_se` even says 'exactly'.

**Fix.** Write 'agrees with `conditional_density` to O(Δ_y²) (the within-cell correction)' and 'reduces asymptotically to `conditional_density_standard_error`'; drop 'exactly' from the `_cpf_value_and_se` comment.

- `src/ConditionalExpectations.jl`:523  (conditional_pushforward_density)
- `src/ConditionalExpectations.jl`:610  (conditional_pushforward_density_standard_error)
- pass 1: verifier `FIXED`  (commit `bb56595c`)

### 2026-10-15-documentation-review-01-003  The notation is not one-to-one and the package has no notation page to settle it: the bin width is `b` in some places and `Δ` in others while `b` is also the asymmetric-kernel scale; `ℓ` is at once the bin location ℓ_j, the moment index in T_ℓ / P_ℓ (both in one formula of `local_fit`) and the leverage ℓ(x); `w` is the design weight w_j, the effective regression weight w_j(x), the AMISE local-scale weight w(x) and the hyperstatistical mixing weight w_k.

Rule: `RULE.DOC.ONE_NOTATION`.  
Reported by `full`, decided by an agent.

`Histogram1D` and ARCHITECTURE §2/§4.1 call the bin width ``b`` (``b \ll h``, ``O(b^2)``); `density`, `expectation`, `conditional_expectation` and ARCHITECTURE §4.1/§4.B/§11 call it ``\Delta``; `density`, `bandwidth`, `AsymmetricKernel`, `gamma_kernel` and ARCHITECTURE §4.5 use ``b`` for the gamma/beta kernel scale, so ARCHITECTURE alone binds `b` to two concepts. `local_fit` writes ``T_\ell = \sum_j K_h(\ell_j - x) w_j (\ell_j - x)^\ell`` with ℓ as both index and location, and `conditional_expectation_standard_error` then uses ``\ell(x)`` for the leverage. `RULE.DOC.ONE_NOTATION`'s `notation_page` is unset, so it names `docs/src/notation.md`, which does not exist.

**Fix.** Pick one symbol per concept (e.g. Δ for the bin width everywhere, b reserved for the asymmetric scale; ℓ_j → c̃_j or x_j for locations and k for the moment index; distinct letters for the four weights) and sweep every surface; record the choice on a notation page.

- `src/Histograms.jl`:24  (Histogram1D)
- `src/MarginalDensity.jl`:486  (density)
- `src/Expectations.jl`:190  (expectation)
- `docs/ARCHITECTURE.md`:139  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / Why a fine histogram is the canonical ``T`` (and the honest caveat))
- `docs/ARCHITECTURE.md`:329  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.1 Marginal density ``f_X``)
- `docs/ARCHITECTURE.md`:412  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `src/MarginalDensity.jl`:435  (density)
- `src/Bandwidth.jl`:366  (bandwidth)
- `src/LocalPolynomial.jl`:75  (local_fit)
- `src/ConditionalExpectations.jl`:164  (conditional_expectation_standard_error)
- `src/Kernels.jl`:85  (AsymmetricKernel)
- `src/Accumulators.jl`:76  (merge!)
- `docs/ARCHITECTURE.md`:817  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:103  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The algebra of measures — a structural lens)
- `docs/ARCHITECTURE.md`:563  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- `docs/ARCHITECTURE.md`:775  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `src/Bandwidth.jl`:367  (bandwidth)
- `src/Kernels.jl`:89  (AsymmetricKernel)
- `src/Kernels.jl`:146  (gamma_kernel)
- `src/Kernels.jl`:175  (beta_kernel)
- `src/MarginalDensity.jl`:451  (density)

### 2026-10-15-documentation-review-01-004  In the two struct docstrings the `# Fields` section is not last: it precedes `# Arguments`, a code example and the reference, whereas the home order puts `# References` and then `# Fields` after all prose.

Rule: `RULE.DOC.HOME_ORDER`.  
Reported by `full`, decided by an agent.

`AsymmetricKernel`'s docstring runs identity → `# Fields` → `# Arguments` → julia example → `# References`; `ParametricStart`'s runs identity → `# Fields` → `# Arguments` → jldoctest → 'Reference: N. L. Hjort & I. K. Glad (1995)'. RULE.DOC.HOME_ORDER: 'then `# References`, then `# Fields` last'.

**Fix.** Move each `# Fields` section to the end of its docstring, after the reference.

- `src/Kernels.jl`:69  (AsymmetricKernel)
- `src/MarginalDensity.jl`:48  (ParametricStart)
- `src/Histograms.jl`:52  (Histogram1D)
- `src/Histograms.jl`:227  (Histogram2D)

### 2026-10-15-documentation-review-01-005  The umbrella module docstring names its operation families in a prose run with no link to any of them (only the two summary types are linked) and points the reader to 'the API reference' instead of one linked one-liner per family.

Rule: `RULE.DOC.MODULE_MAP`.  
Reported by `full`, decided by an agent.

`GBGStreamingStatistics`'s docstring: 'it provides — post hoc, from that one summary — marginal densities, conditional expectations, conditional and two-time transition densities, bandwidth selection, the power-sum quadrature readouts (moments, expectations, and pushforward densities), and standard-error bands for the point estimators. See the API reference for the full, module-by-module surface.' None of `density`, `conditional_expectation`, `conditional_density`, `transition_density`, `bandwidth`, `expectation`, `pushforward_density`, `observe!`, `merge!` is linked; the submodules carry no docstrings at all, so this is the only module map a REPL reader gets.

**Fix.** Replace the prose run with a list of one linked one-liner per family (summaries and accumulation; marginal density; conditional expectations; conditional / transition density; power-sum readouts; bandwidth; HDF5 persistence), each `@ref`-linked to its home, no formula.

- `src/GBGStreamingStatistics.jl`:24  (GBGStreamingStatistics)
- `src/GBGStreamingStatistics.jl`:27  (GBGStreamingStatistics)
- pass 1: verifier `FIXED`  (commit `1705cf9d`)

### 2026-10-15-documentation-review-01-006  The parenthetical attribution of the local-polynomial estimator to Fan & Gijbels carries no year, so it is neither author–year nor an eponymous method name.

Rule: `RULE.CIT.STYLE`.  
Reported by `full`, decided by an agent.

`conditional_mean`: '**Estimator** (local-polynomial; Fan & Gijbels).' and ARCHITECTURE §4.2: '**Estimator** (local-linear; Fan & Gijbels):' — both are citations of the 1996 book (the same docstring later writes '(Fan & Gijbels 1996)' and its reference line gives the year), not a method named after them.

**Fix.** Write 'Fan & Gijbels 1996' at both sites.

- `src/ConditionalExpectations.jl`:256  (conditional_mean)
- `docs/ARCHITECTURE.md`:368  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.2 Conditional expectation ``m(x)`` and variance ``\sigma^2(x)``)
- pass 1: verifier `FIXED`  (commit `bb56595c`)

### 2026-10-15-documentation-review-01-007  `Histogram1D`'s docstring says the power-sums are accumulated in `accumulator`, a name that is neither a constructor argument nor a field; the keyword is `accumulator_type`.

Rule: `RULE.DOC.VERIFIED`.  
Reported by `full`, decided by an agent.

The signature line reads `Histogram1D(lower, upper, nbins; order = 4, accumulator_type = Float64, count_type = UInt64)` and the Arguments bullet names `accumulator_type`; the struct has fields lower, upper, scale, counts, powersums. The prose 'power-sums are accumulated in `accumulator` (`Float64` by default)' names something that does not exist.

**Fix.** Write 'accumulated in `accumulator_type` (`Float64` by default)'.

- `src/Histograms.jl`:35  (Histogram1D)
- pass 1: verifier `FIXED`  (commit `1705cf9d`)

### 2026-10-15-documentation-review-01-008  The internal kernel-moment quadrature `_quad` takes two defaulted arguments (`L`, `npts`) positionally, though neither's type selects a method.

Rule: `RULE.API.POSITIONAL`.  
Reported by `full`, decided by an agent.

`function _quad(g, L::Float64 = 16.0, npts::Int = 32_769)` in `Kernels`; every caller passes `g` alone. The rule makes any defaulted argument a keyword unless its type selects among methods; the package's own CLAUDE.md reading applies it to exported functions, so this may be out of the rule's intended reach — reported for the owner to decide.

**Fix.** `_quad(g; L = 16.0, npts = 32_769)`, or leave as is if the owner reads the rule as exported-only.

!!! warning "The code may be what is wrong"
    The reader believed the documentation and the code disagree, and that the code is the wrong half. That is a question for the owner, not a documentation fix.

- `src/Kernels.jl`:43  (_quad)

### 2026-10-15-documentation-review-01-009  Rendered code blocks carry lines of 93 to 114 characters, well past the 86-to-89-character width of the content column, so each block renders with a horizontal scrollbar; the overflow is almost always a trailing # comment.

Rule: `RULE.TEX.NO_HSCROLL`.  
Reported by `rendered`, decided by an agent.

Measured on the built pages: the article is 52rem wide with 1rem padding (800px); a pre block adds 1.5rem padding per side and a 2px border, a docstring section 0.75rem per side; at 0.875em JuliaMono (0.6em advance) that is ~89 characters at page level and ~86 inside a docstring, and pre is overflow-x:auto with white-space:pre. Lines: Expectations 96 and 95; Histograms 95; Kernels 93; MarginalDensity 108, 98 and 104; Architecture section 8 layout block 100, 100, 104, 107, 98, 102, 114 and 97; tutorial 97, 100 and 94. Every display formula on the fourteen pages fits.

**Fix.** Move the trailing comment onto its own line above the statement, or shorten it; in the Architecture layout block break the long per-file comments onto continuation lines and split the test/ line into two.

- `src/Expectations.jl`:203  (expectation)
- `src/Expectations.jl`:270  (expectation_with_standard_error)
- `src/Histograms.jl`:243  (Histogram2D)
- `src/Kernels.jl`:99  (AsymmetricKernel)
- `src/MarginalDensity.jl`:498  (density)
- `src/MarginalDensity.jl`:68  (ParametricStart)
- `src/MarginalDensity.jl`:69  (ParametricStart)
- `docs/ARCHITECTURE.md`:583  (GBGStreamingStatistics — Architecture and Rationale / 8. Proposed repository layout (mirrors PF))
- `docs/ARCHITECTURE.md`:584  (GBGStreamingStatistics — Architecture and Rationale / 8. Proposed repository layout (mirrors PF))
- `docs/ARCHITECTURE.md`:586  (GBGStreamingStatistics — Architecture and Rationale / 8. Proposed repository layout (mirrors PF))
- `docs/ARCHITECTURE.md`:587  (GBGStreamingStatistics — Architecture and Rationale / 8. Proposed repository layout (mirrors PF))
- `docs/ARCHITECTURE.md`:588  (GBGStreamingStatistics — Architecture and Rationale / 8. Proposed repository layout (mirrors PF))
- `docs/ARCHITECTURE.md`:593  (GBGStreamingStatistics — Architecture and Rationale / 8. Proposed repository layout (mirrors PF))
- `docs/ARCHITECTURE.md`:597  (GBGStreamingStatistics — Architecture and Rationale / 8. Proposed repository layout (mirrors PF))
- `docs/ARCHITECTURE.md`:600  (GBGStreamingStatistics — Architecture and Rationale / 8. Proposed repository layout (mirrors PF))
- `docs/src/tutorial.md`:23  (Getting started / The workflow)
- `docs/src/tutorial.md`:31  (Getting started / The workflow)
- `docs/src/tutorial.md`:65  (Getting started / The workflow)

### 2026-10-15-documentation-review-01-010  Code lines of 87 to 92 characters sit at or a few pixels past the column width inside a docstring block (about 86 characters) or the Architecture page (about 89), so they render with a slight scrollbar on the default theme; marginal, and font-fallback dependent.

Rule: `RULE.TEX.NO_HSCROLL`.  
Reported by `rendered`, decided by an agent.

Same measurement as above; the docstring pre's inner width is ~722px, which at 8.4px per character is 86 characters. ConditionalDensity's three signature continuation lines are 87; Histograms 88, 90 and 92; Kernels 90; Architecture 92 and 88. The tutorial's 87-character line fits at page level and is not listed.

**Fix.** Shorten or drop the trailing comment, or let the aligned comment column start earlier; the three signature lines lose their trailing ≡ comment or move it to the line above.

- `src/ConditionalDensity.jl`:86  (conditional_density)
- `src/ConditionalDensity.jl`:239  (conditional_density_standard_error)
- `src/ConditionalDensity.jl`:178  (transition_density)
- `src/Histograms.jl`:74  (Histogram1D)
- `src/Histograms.jl`:241  (Histogram2D)
- `src/Histograms.jl`:245  (Histogram2D)
- `src/Kernels.jl`:97  (AsymmetricKernel)
- `docs/ARCHITECTURE.md`:585  (GBGStreamingStatistics — Architecture and Rationale / 8. Proposed repository layout (mirrors PF))
- `docs/ARCHITECTURE.md`:590  (GBGStreamingStatistics — Architecture and Rationale / 8. Proposed repository layout (mirrors PF))

### 2026-10-15-documentation-review-01-011  Three display rows pack a second equation on the same line behind \qquad (the formula, then a symbol's definition), the form the rule names; each still fits the column today, so this is the rule's letter rather than a visible scroll.

Rule: `RULE.TEX.NO_HSCROLL`.  
Reported by `rendered`, decided by an agent.

The rule reads 'Never pack several equations on one line behind \quad.' On the built Bandwidth page the row sits inside an aligned block whose other rows hold one equation each, so a '\\' row for \mu_2(K) is the environment's own idiom; the Accumulators and ConditionalDensity rows are single-line displays with a trailing where-clause. The CLAUDE.md note that \quad/\qquad are correct forms is about the spacing command, not about packing. The widest of the three measures roughly 60 visible characters, under the 800px column.

**Fix.** Give the second equation its own row (a '\\' in the aligned block; a gathered or aligned pair for the two single-line displays), or leave as is if the owner reads the where-clause as a sanctioned exception and records that.

- `src/Bandwidth.jl`:336  (bandwidth)
- `src/Accumulators.jl`:76  (Accumulators)
- `src/ConditionalDensity.jl`:250  (conditional_density_standard_error)

### 2026-10-15-documentation-review-01-012  The page's opening sentence has no verb: it is a noun phrase with a participial tail, not a sentence.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api`, decided by an agent.

"The complete public API, grouped by module on the pages that follow." — subject with no finite verb; read aloud it does not close as a sentence.

**Fix.** Give it a verb, e.g. "This is the complete public API, grouped by module on the pages that follow." or "The pages that follow list the complete public API, grouped by module."

- `docs/src/api.md`:3  (API overview)
- pass 1: verifier `FIXED`  (commit `ceea748c`)

### 2026-10-15-documentation-review-01-013  The package noun "summary" (the accumulated object every estimator reads from) is used twice on the page and defined nowhere on it.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api`, decided by an agent.

"before the estimator that approximates it from the accumulated summary" and "holds the summaries and the accumulation primitives" — a reader new to the package does not know what a summary is, what it is accumulated from, or why an estimator approximates from one rather than from the data; neither phrase glosses it and neither links to a definition.

**Fix.** Bind it in a clause at first use, e.g. "the accumulated summary — the fixed-size binned record of the sample stream that every estimator reads instead of the samples themselves", or link the first mention to the page that defines it.

- `docs/src/api.md`:4  (API overview)
- `docs/src/api.md`:8  (API overview)
- pass 1: verifier `FIXED`  (commit `32066856`)

### 2026-10-15-documentation-review-01-014  The range "Expectations through Bandwidth" names an ordering of modules that this page does not show, so the reader cannot tell which modules the upper half contains.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api`, decided by an agent.

"the upper half ([`Expectations`](Expectations/) through [`Bandwidth`](Bandwidth/)) holds the estimands" — the lower half is enumerated in full, the upper half only by its two endpoints; on this page alone there is no list between them to read "through" against.

**Fix.** Enumerate the upper-half modules as the lower half is enumerated, or drop "through" and say "the remaining modules, from Expectations to Bandwidth" with the list.

- `docs/src/api.md`:9  (API overview)
- pass 1: verifier `NOT_FIXED`  (commit `454ccff7`) — the verifier said FIXED, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-015  Label-style fragments stand as sentences without a verb: a dash-joined noun phrase closed by a full stop reads as a heading, not prose.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-Accumulators`, decided by an agent.

"Weighted form — the cone." and "Two uses." have no verb; "Hyperstatistics — the scalar building block ... assembled by weighted `merge!` (...)" and "Like `merge!`, an intrinsically finite-sample operation ...; specified by the cone laws ..." have no subject-verb pair; the Operations lead-in is a noun phrase with a dash and no verb.

**Fix.** Give each a verb: "The weighted form is the cone." / "It has two uses." / "Hyperstatistics uses it as the scalar building block ..." / "Like `merge!`, it is an intrinsically finite-sample operation ..." / "These operations form the additive monoid ... and its conic extension ...".

- `src/Accumulators.jl`:67  (merge!)
- `src/Accumulators.jl`:112  (scale!)
- `src/Accumulators.jl`:114  (scale!)
- `src/Accumulators.jl`:120  (scale!)
- `docs/src/api/Accumulators.md`:13  (Accumulators / Operations)
- pass 1: verifier `FIXED`  (commit `bef5801e`)

### 2026-10-15-documentation-review-01-016  Package terms of art are used before the page defines them, or are never defined or linked on it.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Accumulators`, decided by an agent.

"hyperstatistics" and "exponential forgetting" appear in the Operations lead-in before `merge!` and `scale!` define them further down; "the conditional / transition estimators" and "the local mixing weight" are named nowhere else on the page and carry no link; "order" in "(range, bin count, order)" is the package's power-sum order, not said; "delivers no reruns" assumes the reader knows what a rerun is (the page's own gloss, "no second pass over the data", sits in another docstring); "WARPing" is a named construction left unexpanded.

**Fix.** Define at first use or move the definition up: "hyperstatistics (a post-hoc mixture over parameter-indexed summaries, defined under `merge!`) and exponential forgetting (geometric down-weighting of the past, under `scale!`)"; "(range, bin count, power-sum order)"; "the conditional and transition density estimators of the Histograms readouts" with a link; "no recompute from the raw stream"; "WARPing (weighted averaging of rounded points)".

- `docs/src/api/Accumulators.md`:13  (Accumulators / Operations)
- `src/Accumulators.jl`:82  (merge!)
- `src/Accumulators.jl`:94  (merge!)
- `src/Accumulators.jl`:37  (observe!)
- `src/Accumulators.jl`:49  (observe!)
- pass 1: verifier `FIXED`  (commit `bef5801e`)

### 2026-10-15-documentation-review-01-017  Symbols are used in a docstring that binds them nowhere in itself; they lean on another docstring on the page or are never bound at all.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-Accumulators`, decided by an agent.

`scale!` writes the reweighted measure $\texttt{weight} \cdot \hat\mu$ and the fractional $\texttt{weight} \cdot n$, but $\hat\mu$ and $n$ are bound only in `AbstractAccumulator`'s docstring; `merge!` likewise relies on "the normalized $\hat\mu$" without binding it; in `observe!` the index $p$ of $S_p$ is never said to be the moment order; in `merge!`'s display, $P(x \mid \beta)$ and the spacing $\Delta\beta$ are introduced without a word.

**Fix.** Bind in apposition at first use within each docstring: "the reweighted empirical measure $\texttt{weight} \cdot \hat\mu$ ($\hat\mu$ the empirical measure of the $n$ samples seen)"; "the per-bin power-sums $S_p = \sum (x_i - c_j)^p$ about the bin center $c_j$, $p$ up to the summary's order"; "$P(x \mid \beta)$ the density at fixed parameter, $\Delta\beta$ the spacing of the $\beta_k$".

- `src/Accumulators.jl`:107  (scale!)
- `src/Accumulators.jl`:118  (scale!)
- `src/Accumulators.jl`:61  (merge!)
- `src/Accumulators.jl`:40  (observe!)
- `docs/FORMALIZATION.md`:97  (Formalization / The cone — `scale!`, weighted `merge!`, hyperstatistics (§10))

### 2026-10-15-documentation-review-01-018  A lone equals sign is set as a code span where it means mathematical equality, not Julia (where `=` is assignment).

Rule: `RULE.TEX.CODE_SPAN_JULIA`.  
Reported by `withheld@api-Accumulators`, decided by an agent.

"tested as accumulate-then-merge `=` accumulate-together" backticks a bare `=` standing for "equals".

**Fix.** Write it as prose or mathematics: "tested as accumulate-then-merge equals accumulate-together".

!!! warning "The code may be what is wrong"
    The reader believed the documentation and the code disagree, and that the code is the wrong half. That is a question for the owner, not a documentation fix.

- `src/Accumulators.jl`:88  (merge!)
- pass 1: verifier `NOT_FIXED`  (commit `bef5801e`) — the verifier said FIXED, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-019  A mutator whose signature already carries a `-> into` arrow also states the return in a sentence, so arrow and sentence both appear.

Rule: `RULE.DOC.RETURNS`.  
Reported by `withheld@api-Accumulators`, decided by an agent.

`merge!` has `-> into` in its signature, and its Arguments bullet adds "updated in place and returned".

**Fix.** Drop "and returned" from the bullet (keep "updated in place"), letting the arrow carry the return.

- `src/Accumulators.jl`:93  (merge!)
- pass 1: verifier `FIXED`  (commit `bef5801e`)

### 2026-10-15-documentation-review-01-020  The main input $\hat\sigma$ (and $\hat\sigma_X$) is never bound: nowhere does the page say it is the sample standard deviation, and 'recovered' is the package's verb for reading it off the summary, unexplained.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-Bandwidth`, decided by an agent.

The normal-reference formula, the Arguments bullet and the local pilot all lean on $\hat\sigma$; a first-time reader can guess 'a scale estimate' but is told neither what it estimates nor from what. 'Recovered' is used three times as if it named a known operation.

**Fix.** At first use: 'From the sample standard deviation $\hat\sigma$ read back from the summary's stored moments and the sample size $n$'; in the second docstring bind $\hat\sigma_X$ as the standard deviation of the conditioning variable.

- `src/Bandwidth.jl`:321  (bandwidth)
- `src/Bandwidth.jl`:395  (bandwidth)
- `src/Bandwidth.jl`:523  (bandwidth)
- `docs/ARCHITECTURE.md`:408  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:722  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:728  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:775  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:787  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- pass 1: verifier `NOT_FIXED`  (commit `b4607ee7`) — the verifier said FIXED, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-021  Symbols of the plug-in machinery are used without being bound: the pilot functional $\psi_6$, the derivative order $s$ in $\|f^{(s)}\|^2$, and Chen's 'modified $f_2$ kernels'.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-Bandwidth`, decided by an agent.

$\psi_4$ is bound as $R(f'')$ but $\psi_6$ only as 'its pilot'; $s$ appears once with no statement of what it ranges over; $f_2$ is Chen's own label for a kernel variant and means nothing to a reader of this page alone.

**Fix.** Bind each in apposition: '$\psi_6 = R(f''')$, the pilot functional'; 'the curvature norms $\|f^{(s)}\|^2$, $s = 2, 3, \dots$'; 'the boundary-modified gamma/beta kernels (Chen's $f_2$ variant)'.

- `src/Bandwidth.jl`:346  (bandwidth)
- `src/Bandwidth.jl`:356  (bandwidth)
- `src/Bandwidth.jl`:369  (bandwidth)
- `src/Bandwidth.jl`:347  (bandwidth)
- `docs/ARCHITECTURE.md`:753  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- pass 1: verifier `NOT_FIXED`  (commit `b4607ee7`) — the verifier said FIXED, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-022  'The eval kernel' is a term of the selector machinery used in three Estimator blocks and defined in none; 'Hermite-weighted reductions' and the abbreviation DCT are likewise used before any clause says what they are.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Bandwidth`, decided by an agent.

The page never says that a plug-in selector estimates its curvature with one (pilot, Gaussian) kernel and then rescales for the kernel the density will actually be evaluated with; 'eval kernel' silently assumes that two-kernel picture. DCT is expanded nowhere (AMISE is, so the page knows the convention).

**Fix.** In the `:sj` block, define once: 'the evaluation kernel — the `kernel` the density is smoothed with, as opposed to the Gaussian pilot used for the curvature — then enters as ...'; expand 'discrete cosine transform (DCT)' at its first mention; say 'sums over the bins weighted by Hermite polynomials'.

- `src/Bandwidth.jl`:349  (bandwidth)
- `src/Bandwidth.jl`:347  (bandwidth)
- `src/Bandwidth.jl`:358  (bandwidth)
- `docs/src/api/Bandwidth.md`:5  (Bandwidth)
- `docs/ARCHITECTURE.md`:413  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:725  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:734  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:740  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:755  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `src/Bandwidth.jl`:334  (bandwidth)
- `src/Bandwidth.jl`:335  (bandwidth)
- `src/Bandwidth.jl`:359  (bandwidth)
- `src/Bandwidth.jl`:362  (bandwidth)
- `src/Bandwidth.jl`:401  (bandwidth)
- pass 1: verifier `FIXED`  (commit `b4607ee7`)

### 2026-10-15-documentation-review-01-023  Package vocabulary stands unglossed: 'mergeable', 'reductions of the bin counts', 'the one stored summary, post hoc', and 'a recorded refinement' assume the reader knows the package's summary-and-merge model and its roadmap register.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Bandwidth`, decided by an agent.

A cold reader does not know that summaries are merged across streams, so 'they are mergeable' states a property whose relevance the page never gives; 'a recorded refinement' is the design documents' phrase for planned work and reads as nonsense here.

**Fix.** Gloss once: 'so two summaries accumulated separately give the same selection as one (they merge)'; replace 'a recorded refinement' with 'not yet offered' or 'planned'; 'read back from the stored summary after accumulation' for 'post hoc'.

- `src/Bandwidth.jl`:352  (bandwidth)
- `src/Bandwidth.jl`:343  (bandwidth)
- `src/Bandwidth.jl`:385  (bandwidth)
- pass 1: verifier `FIXED`  (commit `b4607ee7`)

### 2026-10-15-documentation-review-01-024  Two places do not read as sentences aloud: the header's verbless parenthetical 'Silverman the Gaussian/degree-0 case', and the canonical-bandwidth display that ends in a comma and is followed by a capitalised new sentence.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-Bandwidth`, decided by an agent.

'(the default; Silverman the Gaussian/degree-0 case)' has no verb; the aligned display's last line ends ', \end{aligned}' and the prose resumes with 'At `degree` 2/3 the constant is ...', so the comma dangles.

**Fix.** '(the default; Silverman's rule is its Gaussian, degree-0 case)'; end the display with a full stop, or continue it with a lowercase 'and at `degree` 2/3 ...'.

- `docs/src/api/Bandwidth.md`:4  (Bandwidth)
- `src/Bandwidth.jl`:340  (bandwidth)
- pass 1: verifier `FIXED`  (commit `b4607ee7`)

### 2026-10-15-documentation-review-01-025  The joint-summary signature states its result as mathematics inside the code block, a comment `# ≡ h⋆(x)`, while the marginal signature uses `-> h`; the two methods of one function declare their result in two different registers.

Rule: `RULE.TEX.CODE_SPAN_JULIA`.  
Reported by `withheld@api-Bandwidth`, decided by an agent.

A Julia comment is legal code, so this may be intended, but the `≡` and the star are mathematics rendered in a code font, and RULE.DOC.RETURNS expects an arrow or a sentence; the identity sentence already names $h^\star(x)$, so the comment adds nothing.

**Fix.** Drop the comment and, if an arrow is wanted, use the same form as the first method; or write `-> h` on both.

!!! warning "The code may be what is wrong"
    The reader believed the documentation and the code disagree, and that the code is the wrong half. That is a question for the owner, not a documentation fix.

- `src/Bandwidth.jl`:501  (bandwidth)
- pass 1: verifier `FIXED`  (commit `b4607ee7`)

### 2026-10-15-documentation-review-01-026  Three docstrings open with a verbless noun phrase (an em-dash apposition in place of a predicate), and one caveat is a verbless fragment; none reads as a sentence aloud.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-ConditionalDensity`, decided by an agent.

The lead sentence of `conditional_density_with_standard_error`, `transition_density_standard_error` and `transition_density_with_standard_error` is 'The X — identical to / the two-time reading of ...' with no finite verb; 'NaN where no data is local to x.' has no verb either.

**Fix.** Give each a predicate: 'Returns the estimate and its standard error ...', 'The standard error of `transition_density` is `conditional_density_standard_error` on ...', 'The estimate is `NaN` where ...'.

- `src/ConditionalDensity.jl`:301  (conditional_density_with_standard_error)
- `src/ConditionalDensity.jl`:341  (transition_density_standard_error)
- `src/ConditionalDensity.jl`:361  (transition_density_with_standard_error)
- `docs/src/api/ConditionalDensity.md`  (conditional_density_standard_error)  — **the text here has moved and was not found again**
- pass 1: verifier `NOT_FIXED`  (commit `5bc70f91`) — the verifier said FIXED, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-027  'profile' / 'x-profile reduction' is a package term used as if known and never defined on the page.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-ConditionalDensity`, decided by an agent.

'the profile and its total', 'normalize its profile', 'one x-profile reduction', 'share the one profile' — nowhere does the page say a profile is the kernelx-weighted reduction of the joint counts to a function of y at the conditioning value x.

**Fix.** At first use, bind it: 'the profile — the kernelx-weighted sum of the cell counts over the x-columns, a function of y at the fixed x — and its total'.

- `src/ConditionalDensity.jl`:111  (conditional_density)
- `src/ConditionalDensity.jl`:260  (conditional_density_standard_error)
- `src/ConditionalDensity.jl`:306  (conditional_density_with_standard_error)
- pass 1: verifier `FIXED`  (commit `5bc70f91`)

### 2026-10-15-documentation-review-01-028  'the marginal Wand–Jones correction' is named as a known thing but neither defined, linked nor cited on this page.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-ConditionalDensity`, decided by an agent.

The reader is told the -f^2/n term is 'the double-kernel analog of the marginal Wand–Jones correction'; no Wand & Jones entry is in the References and no clause says what the correction is (the exact-variance -f^2/n term of the marginal KDE).

**Fix.** Define in a clause ('the exact finite-sample -f^2/n term of the marginal kernel-density variance, Wand & Jones 1995') and add the reference, or link `density_standard_error` on that sentence instead of naming the correction.

- `src/ConditionalDensity.jl`:253  (conditional_density_standard_error)
- pass 1: verifier `FIXED`  (commit `5bc70f91`)

### 2026-10-15-documentation-review-01-029  'per-lag joint histogram' and 'the summary's lag' presume a lag the page never shows a `Histogram2D` carrying; the reader cannot tell where tau enters.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-ConditionalDensity`, decided by an agent.

The doctest constructs `Histogram2D(-1.0, 1.0, 2, -1.0, 1.0, 2)` with no lag argument, yet `transition_density` speaks of 'the summary's lag tau' and 'the per-lag joint histogram'. Nothing on the page says the lag is fixed by the caller pairing each sample X_t with X_{t+tau} before `observe_pairs!`.

**Fix.** State it once in the Estimator paragraph: 'a `Histogram2D` fed the pairs (X_t, X_{t+tau}) for one fixed lag tau chosen by the caller — the summary itself stores no lag'.

- `docs/src/api/ConditionalDensity.md`  (transition_density)  — **the text here has moved and was not found again**
- `src/ConditionalDensity.jl`:186  (transition_density)
- pass 1: verifier `NOT_FIXED`  (commit `5bc70f91`) — the verifier said FIXED, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-030  The sample pairs X_i, Y_i, the sample size n and the semicolon kernel form K_2(c; y, h_2) appear in formulas without being bound.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-ConditionalDensity`, decided by an agent.

The first display sums over i with X_i, Y_i never introduced as the observed pairs; n first appears in the variance formula as a bare denominator and is never said to be the number of pairs; K_2(c^y_k; y, h_2) is a different kernel notation from K_{h_2}(c^y_k - y) and is not explained.

**Fix.** Before the first display: 'for n observed pairs (X_i, Y_i)'; and gloss the asymmetric form: 'K_2(u; y, h_2), the kernel's density in its data argument u with location y and scale h_2'.

- `src/ConditionalDensity.jl`:99  (conditional_density)
- `src/ConditionalDensity.jl`:248  (conditional_density_standard_error)
- `src/ConditionalDensity.jl`:118  (conditional_density)
- pass 1: verifier `LINE_ONLY`  (commit `5bc70f91`)

### 2026-10-15-documentation-review-01-031  The standard-error docstring uses N_j and c^x_j bare; they are bound only in the parent `conditional_density` docstring, so the docstring does not read on its own (as in the REPL).

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-ConditionalDensity`, decided by an agent.

'the x-window denominator D = sum_j N_j K_{h_1}(c^x_j - x) = n f_X(x)' — N_j (the per-column count) and c^x_j (the column center) are defined nowhere in this docstring.

**Fix.** Re-bind in apposition: 'D = sum_j N_j K_{h_1}(c^x_j - x), the kernel-weighted total of the per-column counts N_j at the column centers c^x_j'.

- `src/ConditionalDensity.jl`:259  (conditional_density_standard_error)
- pass 1: verifier `FIXED`  (commit `5bc70f91`)

### 2026-10-15-documentation-review-01-032  Compressed prose lets a symbol or an ellipsis stand in for words: '$\approx$' used as an adverb, and a rate compared with a density.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-ConditionalDensity`, decided by an agent.

'but only ≈ normalized in y' reads a math symbol as the word 'approximately'; 'The conditional MISE rate is slower than a marginal density' compares a rate to a density (the marginal density's rate is meant).

**Fix.** 'but only approximately normalized in y'; 'The conditional MISE rate, O(n^{-2/3}), is slower than a marginal density's.'

- `src/ConditionalDensity.jl`:119  (conditional_density)
- `src/ConditionalDensity.jl`:120  (conditional_density)

### 2026-10-15-documentation-review-01-033  Package jargon — 'reader' for a function that reads an estimate off the summary, and 'a recorded gap' for a documented unsupported case — is used without a defining clause.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-ConditionalDensity`, decided by an agent.

The module line ends 'the paired estimate-plus-error reader'; the SE docstring says 'the asymmetric-kernely variance is a recorded gap'. A newcomer does not know a reader is a function of the summary or where a gap is recorded.

**Fix.** 'the paired reader `conditional_density_with_standard_error`, returning estimate and error together'; 'the asymmetric-kernely variance is not implemented (a known gap)'.

- `docs/src/api/ConditionalDensity.md`:3  (ConditionalDensity)
- `src/ConditionalDensity.jl`:265  (conditional_density_standard_error)
- pass 1: verifier `FIXED`  (commit `5bc70f91`)

### 2026-10-15-documentation-review-01-034  The signature code blocks end in a comment holding mathematics (LaTeX-style subscripts, a combining hat, a Unicode divider) that renders verbatim, neither as code nor as math.

Rule: `RULE.TEX.CODE_SPAN_JULIA`.  
Reported by `withheld@api-ConditionalDensity`, decided by an agent.

'# ≡ f_{Y|X}', '# ≡ SE[f̂(y∣x)]' and '# ≡ p_τ(y|x)' sit inside the fenced signature; the braces and the combining hat are LaTeX idiom, not Julia, and a reader sees raw source. Possibly a deliberate house apposition — flagged for a decision, not a verdict.

**Fix.** Drop the comment from the signature and keep the estimand apposition in the prose line that already states it (`conditional_density` (f_{Y|X})), or write the comment in plain words: '# the conditional density of Y given X'.

!!! warning "The code may be what is wrong"
    The reader believed the documentation and the code disagree, and that the code is the wrong half. That is a question for the owner, not a documentation fix.

- `src/ConditionalDensity.jl`:86  (conditional_density)
- `src/ConditionalDensity.jl`:239  (conditional_density_standard_error)
- `src/ConditionalDensity.jl`:178  (transition_density)
- pass 1: verifier `FIXED`  (commit `5bc70f91`)

### 2026-10-15-documentation-review-01-035  The term of art 'jet' is used in the `degree` bullets of three docstrings that never define it; the page defines it only inside `conditional_mean`, which comes after `conditional_expectation` on the page.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent.

A cold reader of `conditional_expectation`, `conditional_variance` or `conditional_pushforward_density` meets 'which has no jet' with no clause saying a jet is the vector of local-polynomial coefficients; the gloss '(the jet $(\beta_0, \dots, \beta_{\texttt{degree}})$ — the local-polynomial coefficients at $x$' sits only in `conditional_mean`.

**Fix.** Say what is lost in a clause at each site — e.g. 'which fits no local polynomial (no jet of coefficients) and so admits only `0`' — or link the word to where the jet is defined.

- `src/ConditionalExpectations.jl`:113  (conditional_expectation)
- `src/ConditionalExpectations.jl`:113  (conditional_expectation)
- `src/ConditionalExpectations.jl`:542  (conditional_pushforward_density)
- `src/ConditionalExpectations.jl`:302  (conditional_mean)
- `src/ConditionalExpectations.jl`:381  (conditional_variance)
- pass 1: verifier `FIXED`  (commit `bb56595c`)

### 2026-10-15-documentation-review-01-036  The summary's structural vocabulary — '$y$-quadrature', '$X$-column', 'order-2 $Y$-power-sums', '$Y$-bin scale', 'column sums' — is used as known on this page, which never says what a power-sum is or what a column of the joint summary is.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent.

The page preamble ends on 'a $y$-quadrature response' before any sentence explains a quadrature of what; `conditional_expectation` then speaks of 'the stored order-2 $Y$-power-sums' and 'per-$X$-column' as given. A reader who has not read the Histograms page cannot tell that a power-sum is the per-cell sum of $Y^p$ or that an $X$-column is one $x$-bin's row of cells; `Histogram2D` is linked but the nouns are not.

**Fix.** Bind the nouns once where they first appear: 'the per-cell sums $\sum Y$, $\sum Y^2$ the summary stores (its order-2 $Y$-power-sums)', 'each $x$-bin of the joint summary (an $X$-column)', and 'the $Y$-bin width $\Delta_y$'.

- `docs/src/api/ConditionalExpectations.md`:9  (ConditionalExpectations)
- `src/ConditionalExpectations.jl`:90  (conditional_expectation)
- `src/ConditionalExpectations.jl`:98  (conditional_expectation)
- `src/ConditionalExpectations.jl`:108  (conditional_expectation)
- pass 1: verifier `NOT_FIXED`  (commit `bb56595c`) — the verifier said FIXED, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-037  The symbol $\Delta_y$ is used in three docstrings and bound nowhere on the page.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent.

'$O(\Delta_y^2)$' appears in `conditional_expectation` and `conditional_variance_standard_error`, and '$h_z \gtrsim |g'| \Delta_y$' in `conditional_pushforward_density`, with no clause saying $\Delta_y$ is the $Y$-bin width of the summary.

**Fix.** At first use in each docstring: '$O(\Delta_y^2)$, $\Delta_y$ the $Y$-bin width'.

- `src/ConditionalExpectations.jl`:100  (conditional_expectation)
- `src/ConditionalExpectations.jl`:437  (conditional_variance_standard_error)
- `src/ConditionalExpectations.jl`:526  (conditional_pushforward_density)
- `docs/ARCHITECTURE.md`:928  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)

### 2026-10-15-documentation-review-01-038  Several symbols in the formulas are not bound at first use: $K_h$, $n$, $f_X$, $\mu_2(K)$ and $R(K)$ in `conditional_mean` (with $R(K)$ bound only later, in `conditional_pushforward_density_standard_error`), and $N$ in `conditional_variance`.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent.

The design sums use $K_h(X_i-x)$ without saying $K_h(u) = K(u/h)/h$; the consistency line uses $\mu_2(K)$, $R(K)$, $n$ and $f_X(x)$ unglossed, after the same density was written $f(x)$ a sentence earlier; `conditional_variance`'s residual total ends in '$\hat m^2 N$' with $N$ never named (the per-column count is $N_j$ only in another docstring).

**Fix.** Bind in apposition at first use: '$K_h(u) = K(u/h)/h$ the scaled kernel', '$n$ the sample count', '$f_X$ the marginal density of $X$', '$\mu_2(K) = \int u^2 K$, $R(K) = \int K^2$', and '$N$ the column count'; use one symbol for the design density.

- `src/ConditionalExpectations.jl`:275  (conditional_mean)
- `src/ConditionalExpectations.jl`:258  (conditional_mean)
- `src/ConditionalExpectations.jl`:368  (conditional_variance)

### 2026-10-15-documentation-review-01-039  The package's jargon 'reader' / 'readout' (a function that reads an estimate off the summary) is used without definition.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent.

`conditional_expectation_standard_error` calls its estimator 'The leverage × noise-scale readout' and speaks of 'two `conditional_expectation` readouts' and 'the $g^2$ reader'; nothing on the page says a readout is an evaluation of a function of the stored summary.

**Fix.** Replace with plain words ('the product of the leverage and the noise scale', 'two evaluations of `conditional_expectation`, the second with $g^2$ in place of $g$') or define 'readout' once.

- `src/ConditionalExpectations.jl`:168  (conditional_expectation_standard_error)
- `src/ConditionalExpectations.jl`:171  (conditional_expectation_standard_error)
- pass 1: verifier `NOT_FIXED`  (commit `bb56595c`) — the verifier said FIXED, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-040  The type `AsymmetricKernel` appears in every signature's `degree` default and is never named, linked or explained on the page; 'the Chen kernels' is used as a known name with no citation or definition.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent.

A reader meets 'kernel isa AsymmetricKernel ? 0 : 1' eight times; the prose only ever says 'an asymmetric `kernel`' and links `gamma_kernel` / `beta_kernel`. `conditional_mean` then says 'the Chen kernels are a density boundary tool' — no Chen entry is under its References and the page never ties the name to the gamma/beta kernels.

**Fix.** In the first `kernel` bullet, link the type ('an [`AsymmetricKernel`](@ref …) such as `gamma_kernel` / `beta_kernel`') and write 'the gamma/beta boundary kernels of Chen (1999, 2000)' with the reference added, or drop the eponym.

- `src/ConditionalExpectations.jl`:286  (conditional_mean)
- `src/ConditionalExpectations.jl`:81  (conditional_expectation)
- `src/ConditionalExpectations.jl`:158  (conditional_expectation_standard_error)
- `src/ConditionalExpectations.jl`:206  (conditional_expectation_with_standard_error)
- `src/ConditionalExpectations.jl`:249  (conditional_mean)
- `src/ConditionalExpectations.jl`:423  (conditional_variance_standard_error)
- pass 1: verifier `NOT_FIXED`  (commit `5adff8e5`) — the verifier said FIXED, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-041  The labelled 'Estimand.' / 'Estimator.' openers are verbless noun phrases in six places — a label followed by a description with no finite verb.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent.

'Estimator. A local-polynomial regression in $x$ of a per-$X$-column $y$-quadrature of $g$.' and 'Estimator. The two-term $\sqrt{\dots}$: the local-regression leverage … times … minus …' have no verb; read aloud they are fragments. The openers that carry an equation ('$\mathbb{E}[\dots] = \int\dots$') or a 'since … is' clause do read as sentences and are not counted here.

**Fix.** Give each a verb: 'Estimator. The estimate is a local-polynomial regression in $x$ of …', 'Estimator. The band is the two-term …', 'Estimand. The estimand is the pointwise standard error …'.

- `src/ConditionalExpectations.jl`:90  (conditional_expectation)
- `src/ConditionalExpectations.jl`:168  (conditional_expectation_standard_error)
- `src/ConditionalExpectations.jl`:370  (conditional_variance)
- `src/ConditionalExpectations.jl`:509  (conditional_pushforward_density)
- `src/ConditionalExpectations.jl`:602  (conditional_pushforward_density_standard_error)
- `docs/src/api/ConditionalExpectations.md`  (conditional_pushforward_density_standard_error)  — **the text here has moved and was not found again**
- pass 1: verifier `NOT_FIXED`  (commit `bb56595c`) — the verifier said FIXED, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-042  Verbless fragments stand as sentences in the body prose: the vector-case return notes, the 'Two smoothing directions:' list, 'Symmetric `kernelz` only.', '`NaN` where …', the preamble's closing 'All by …', and a parenthesis of two abutting code spans.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent.

'For a vector `x`, a `NamedTuple` of vectors.' and 'Symmetric `kernelz` only.' have no verb; '(`method` `:residual` default, or `:plugin`)' is not parseable as a phrase. The sibling `conditional_expectation_with_standard_error` already writes the same note with a verb, 'returns a `NamedTuple` of vectors'.

**Fix.** 'For a vector `x` it returns a `NamedTuple` of vectors.'; 'It smooths in two directions: …'; 'Only a symmetric `kernelz` is supported.'; 'It is `NaN` where …'; 'All are computed by …'; '(`method = :residual`, the default, or `:plugin`)'.

- `src/ConditionalExpectations.jl`:475  (conditional_variance_with_standard_error)
- `src/ConditionalExpectations.jl`:648  (conditional_pushforward_density_with_standard_error)
- `src/ConditionalExpectations.jl`:520  (conditional_pushforward_density)
- `src/ConditionalExpectations.jl`:612  (conditional_pushforward_density_standard_error)
- `docs/src/api/ConditionalExpectations.md`  (conditional_expectation_standard_error)  — **the text here has moved and was not found again**
- `src/ConditionalExpectations.jl`:472  (conditional_variance_with_standard_error)
- pass 1: verifier `NOT_FIXED`  (commit `bb56595c`) — the verifier said FIXED, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-043  A code span presents a schematic call as literal Julia: `K_hz` is not a binding the page or the package offers, so the span is pseudo-code standing in for the mathematics $K_{h_z}$ written just above it.

Rule: `RULE.TEX.CODE_SPAN_JULIA`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent.

'i.e. `conditional_expectation(summary, y -> K_hz(g(y) - z), x)`' — a reader who types it gets an undefined `K_hz`; the span mixes a real signature with a mathematical symbol transliterated into an identifier.

**Fix.** Either make it runnable, e.g. `conditional_expectation(summary, y -> kernelz((g(y) - z) / hz) / hz, x; kernel = kernelx, degree, h = hx)`, or keep the statement in mathematics and drop the code span.

!!! warning "The code may be what is wrong"
    The reader believed the documentation and the code disagree, and that the code is the wrong half. That is a question for the owner, not a documentation fix.

- `src/ConditionalExpectations.jl`:517  (conditional_pushforward_density)
- pass 1: verifier `NOT_FIXED`  (commit `bb56595c`) — the verifier said FIXED, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-044  `conditional_variance` and `conditional_variance_standard_error` accept 'scalar or vector' `x` but, unlike every sibling on the page, never state what is returned in the vector case or where the result is `NaN`.

Rule: `RULE.DOC.RETURNS`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent.

`conditional_expectation`, `conditional_mean` and `conditional_pushforward_density` each say 'A vector of query points returns a vector' and name the `NaN` case; `conditional_variance`'s only statement of the result is its `x` bullet, '(scalar or vector)', and the standard error has none. The identity sentence names the estimand $\sigma^2(x)$, so this is a weak finding — a missing 'Returns …' shape, not a missing result.

**Fix.** Add the sibling sentence to both: 'A vector of query points returns a vector; where no data is local to $x$ the estimate is `NaN`.'

- `src/ConditionalExpectations.jl`:378  (conditional_variance)
- `docs/ARCHITECTURE.md`:364  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.2 Conditional expectation ``m(x)`` and variance ``\sigma^2(x)``)
- `docs/ARCHITECTURE.md`:367  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.2 Conditional expectation ``m(x)`` and variance ``\sigma^2(x)``)
- `docs/ARCHITECTURE.md`:417  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:419  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:421  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:556  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- `docs/ARCHITECTURE.md`:707  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:855  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/api/ConditionalExpectations.md`:7  (ConditionalExpectations)
- `docs/src/index.md`:105  (GBGStreamingStatistics.jl / Main components)
- `src/Bandwidth.jl`:388  (bandwidth)
- `src/Bandwidth.jl`:390  (bandwidth)
- `src/Bandwidth.jl`:510  (bandwidth)
- `src/Bandwidth.jl`:513  (bandwidth)
- `src/Bandwidth.jl`:518  (bandwidth)
- `src/Bandwidth.jl`:538  (bandwidth)
- `src/ConditionalExpectations.jl`:277  (conditional_mean)
- `src/ConditionalExpectations.jl`:358  (conditional_variance)
- `src/ConditionalExpectations.jl`:363  (conditional_variance)
- `src/ConditionalExpectations.jl`:364  (conditional_variance)
- `src/ConditionalExpectations.jl`:428  (conditional_variance_standard_error)
- `src/ConditionalExpectations.jl`:429  (conditional_variance_standard_error)
- `src/ConditionalExpectations.jl`:432  (conditional_variance_standard_error)
- `src/MarginalDensity.jl`:90  (normal_start)
- `src/MarginalDensity.jl`:104  (lognormal_start)
- `src/MarginalDensity.jl`:118  (gamma_start)
- `src/MarginalDensity.jl`:145  (inverse_gaussian_start)
- `src/MarginalDensity.jl`:471  (density)
- pass 1: verifier `FIXED`  (commit `bb56595c`)

### 2026-10-15-documentation-review-01-045  The two section-opening paragraphs of the page are noun phrases with no finite verb.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-Expectations`, decided by an agent.

The Expectations intro runs 'The power-sum quadrature readouts — a smooth functional of X read straight off the stored power-sums, with no kernel and no query point.' ('read' is a participle; nothing predicates). The Recovered moments intro likewise: 'The global marginal moments — the g = x^p case — recovered post hoc from the power-sums by the shift formula (exact), extending Statistics.'

**Fix.** Give each a verb: 'The power-sum quadrature readouts read a smooth functional of X straight off the stored power-sums, with no kernel and no query point.' / 'The global marginal moments — the g = x^p case — are recovered post hoc from the power-sums by the shift formula, exactly, extending Statistics.'

- `docs/src/api/Expectations.md`:3  (Expectations)
- `docs/src/api/Expectations.md`:23  (Expectations / Recovered moments)
- pass 1: verifier `NOT_FIXED`  (commit `454ccff7`) — the verifier said FIXED, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-046  The labelled 'Estimand.' / 'Estimator.' paragraphs and the paired readers' identity sentences are verbless: a label followed by a noun phrase or a bare formula with appositions.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-Expectations`, decided by an agent.

'Estimator. $\sqrt{(\widehat{\mathbb{E}}[g^2] - \widehat{\mathbb{E}}[g]^2)/n}$ — both moments the same power-sum quadrature (...), $n$ the sample count.' has no verb at all; 'Estimand. The population mean $\mathbb{E}[X]$.' is a noun phrase; 'Estimator. The sample mean ... — recovered exactly ... host-side ...' predicates nothing; 'The estimate and its standard error from one power-sum quadrature pass — the `NamedTuple` ...' opens expectation_with_standard_error with no verb, and pushforward_density_with_standard_error repeats the form. Read aloud, none is a sentence.

**Fix.** Supply the verb once per paragraph: 'The estimand is the population mean ...', 'The estimator is the square root ... of both moments by the same power-sum quadrature, with $n$ the sample count.', 'Returns the estimate and its standard error from one pass as the `NamedTuple` ...'.

- `src/Expectations.jl`:237  (expectation_standard_error)
- `src/Expectations.jl`:96  (mean)
- `src/Expectations.jl`:98  (mean)
- `src/Expectations.jl`:258  (expectation_with_standard_error)
- `src/Expectations.jl`:367  (pushforward_density_with_standard_error)
- `src/Expectations.jl`:111  (var)
- `src/Expectations.jl`:151  (central_moment)
- `docs/ARCHITECTURE.md`:446  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.B Power-sum quadrature readouts)
- `docs/ARCHITECTURE.md`:447  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.B Power-sum quadrature readouts)
- `src/Expectations.jl`:238  (expectation_standard_error)
- `src/Expectations.jl`:109  (var)
- `src/Expectations.jl`:149  (central_moment)
- pass 1: verifier `NOT_FIXED`  (commit `bb56595c`) — the verifier said LINE_ONLY, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-047  The package's own nouns 'power-sum(s)', 'readouts' and 'contraction' / 'contraction atom' are used from the first line without a definition on the page; the power-sums get only a partial binding later.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Expectations`, decided by an agent.

The page opens 'The power-sum quadrature readouts — a smooth functional of X read straight off the stored power-sums' before any sentence says what a power-sum is; 'All three are the same power-sum contraction the recovered moments use' names a 'contraction' nowhere explained; central_moment ends 'with which it shares the contraction atom', an internal name a cold reader cannot parse. The only binding, in expectation, is 'its per-bin power-sums $S_p^{(j)}$ about the center (with $S_0^{(j)} = n_j$)' — which never states $S_p^{(j)} = \sum_{i \in j} (X_i - c_j)^p$, so a reader must infer the definition from the formula that uses it.

**Fix.** In the intro, bind the term once: 'the power-sums — per bin, the sums $S_p^{(j)} = \sum_{i\in j}(X_i - c_j)^p$ of the samples' offsets from the bin center — stored by `Histogram1D`'; replace 'power-sum contraction' / 'contraction atom' with 'the same weighted sum over the power-sums' or drop the internal name.

- `docs/src/api/Expectations.md`:7  (Expectations)
- `src/Expectations.jl`:179  (expectation)
- `src/Expectations.jl`:156  (central_moment)
- `docs/ARCHITECTURE.md`:429  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.B Power-sum quadrature readouts)
- `docs/ARCHITECTURE.md`:432  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.B Power-sum quadrature readouts)
- `src/Expectations.jl`:180  (expectation)
- `src/Expectations.jl`:185  (expectation)
- pass 1: verifier `FIXED`  (commit `bb56595c`)

### 2026-10-15-documentation-review-01-048  'The shift formula' is named six times as the mechanism of every recovered moment and defined nowhere on the page.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Expectations`, decided by an agent.

First met in expectation as 'this is the shift formula behind `mean` / `var` / `central_moment`', then the Recovered moments intro 'by the shift formula (exact)', mean 'by the shift formula', var/std 'Exact via the shift formula', central_moment 'recovered exactly from the power-sums by the shift formula'. The term is the package's (the binomial re-expansion of moments about the bin center into moments about the sample mean), not standard vocabulary, and no clause states it.

**Fix.** Define it at first use: 'the shift formula — the binomial re-expansion $\sum_i (X_i - \bar X)^p = \sum_j \sum_{q \le p} \binom{p}{q} (c_j - \bar X)^{p-q} S_q^{(j)}$ that moves the per-bin power-sums about the center to the sample mean —' and leave later mentions bare.

- `src/Expectations.jl`:189  (expectation)
- `docs/src/api/Expectations.md`:23  (Expectations / Recovered moments)
- `src/Expectations.jl`:114  (var)
- pass 1: verifier `FIXED`  (commit `bb56595c`)

### 2026-10-15-documentation-review-01-049  'host-side' is GPU jargon a statistician reading this page has no way to parse, and nothing on the page mentions a device.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Expectations`, decided by an agent.

mean says 'recovered exactly (not a smoothed estimate) host-side from the per-bin counts and first power-sums'; neither 'host' nor a GPU or device array appears anywhere else on the page, so the adverb is an unexplained term from another layer.

**Fix.** Drop 'host-side' (it is an implementation detail), or say it plainly: 'on the CPU, from the per-bin counts ...'.

- `src/Expectations.jl`:98  (mean)
- pass 1: verifier `NOT_FIXED`  (commit `bb56595c`) — the verifier said FIXED, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-050  Three symbols carry the main formulas and are bound nowhere on the page: `order` (set as $\texttt{order}$), $n$ inside expectation, and $K_h$.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-Expectations`, decided by an agent.

expectation's display sums 'to $p = \texttt{order}$' and the text says 'exact for a polynomial $g$ of degree $\le \texttt{order}$', 'a degree-$\texttt{order}$ polynomial interpolant', central_moment 'for $2 \le p \le \texttt{order}$' — no sentence says `order` is the highest power-sum the summary stores, nor that it is a property of `Histogram1D`. The same display divides by $n$, bound only two docstrings later ('$n$ the sample count' in expectation_standard_error). pushforward_density writes $K_h(g(X_i) - y)$ and then 'the surviving $1/h$', but never states $K_h(u) = K(u/h)/h$; a reader who takes $K_h$ as $K$ with a parameter cannot see where the $1/h$ comes from.

**Fix.** In expectation: 'with $n = \sum_j n_j$ the sample count and `order` the highest power stored by the summary'; in pushforward_density: 'with $K_h(u) = K(u/h)/h$ the kernel scaled to the bandwidth $h$'.

- `src/Expectations.jl`:188  (expectation)
- `src/Expectations.jl`:152  (central_moment)
- `src/Expectations.jl`:303  (pushforward_density)
- `docs/ARCHITECTURE.md`:327  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.1 Marginal density ``f_X``)
- `docs/ARCHITECTURE.md`:570  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- `docs/ARCHITECTURE.md`:714  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:715  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/api/Kernels.md`:4  (Kernels)
- `docs/src/index.md`:88  (GBGStreamingStatistics.jl / Main components)
- `src/Bandwidth.jl`:536  (bandwidth)
- `src/ConditionalDensity.jl`:94  (conditional_density)
- `src/ConditionalDensity.jl`:114  (conditional_density)
- `src/ConditionalDensity.jl`:130  (conditional_density)
- `src/ConditionalExpectations.jl`:111  (conditional_expectation)
- `src/ConditionalExpectations.jl`:297  (conditional_mean)
- `src/ConditionalExpectations.jl`:379  (conditional_variance)
- `src/ConditionalExpectations.jl`:519  (conditional_pushforward_density)
- `src/ConditionalExpectations.jl`:538  (conditional_pushforward_density)
- `src/Expectations.jl`:314  (pushforward_density)
- `src/Kernels.jl`:4
- `src/Kernels.jl`:19  (gaussian_kernel)
- `src/Kernels.jl`:30  (epanechnikov_kernel)
- `src/MarginalDensity.jl`:424  (density)
- `src/MarginalDensity.jl`:437  (density)
- `src/MarginalDensity.jl`:700  (density_standard_error)

### 2026-10-15-documentation-review-01-051  Two docstrings admit an input for which they never say what comes back: pushforward_density_standard_error with a vector `y`, and central_moment with $p$ outside $2 \le p \le$ `order`.

Rule: `RULE.DOC.RETURNS`.  
Reported by `withheld@api-Expectations`, decided by an agent.

pushforward_density_standard_error delegates its arguments to pushforward_density, whose `y` is 'a scalar, or a vector for the curve', yet its own text speaks only of the pointwise $\mathrm{SE}(y)$; its paired reader does say 'For a vector `y`, a `NamedTuple` of vectors', so the SE reader alone is silent. central_moment states its estimator 'for $2 \le p \le \texttt{order}$' and never says what $p = 1$, $p = 0$ or $p >$ `order` return (an error, `NaN`, or a value).

**Fix.** Add 'For a vector `y`, a vector of standard errors.' to pushforward_density_standard_error, and to central_moment a sentence on the out-of-range $p$ ('throws an `ArgumentError`' or whatever the code does).

- `src/Expectations.jl`:337  (pushforward_density_standard_error)
- `src/Expectations.jl`:152  (central_moment)
- pass 1: verifier `NOT_FIXED`  (commit `1705cf9d`) — the verifier said FIXED, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-052  The page's opening sentence has no verb: a noun phrase, a dash, a second noun phrase, a full stop.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-GBGStreamingStatistics`, decided by an agent.

"The umbrella module — the assembled public API." is the first thing the visitor reads and cannot be read aloud as a sentence; the rule asks for a subject and a verb in every sentence, with compactness coming from few words, not from a dropped verb.

**Fix.** Give it a verb, e.g. "The umbrella module assembles the public API." or "This is the umbrella module: the assembled public API."

- `docs/src/api/GBGStreamingStatistics.md`:3  (GBGStreamingStatistics)
- pass 1: verifier `FIXED`  (commit `454ccff7`)

### 2026-10-15-documentation-review-01-053  One parenthetical stacks three package terms of art — "cone variant", "hyperstatistics", "forgetting" — and a bare code identifier `count_type`, none defined or linked on the page.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-GBGStreamingStatistics`, decided by an agent.

A scientist meeting the package for the first time cannot know that a "cone" is the summary made scalable as well as mergeable, that "forgetting" means exponential down-weighting of old samples, or what "hyperstatistics" denotes (the standard name in the field is superstatistics, and even that is not said); `count_type = Float64` is a keyword of a constructor the page never names, so the symbol is bound nowhere (RULE.DOC.SELF_CONTAINED applies to it as well). The page gives no clause for any of the four and links none of them.

**Fix.** Either drop the parenthetical from the module map (it belongs in the Histograms home) or gloss it in a clause: "a `count_type = Float64` variant of either histogram whose counts can be rescaled, for mixing summaries with post-hoc weights (superstatistics) and for exponentially forgetting old samples", linking the constructor whose keyword it is.

- `src/GBGStreamingStatistics.jl`:23  (GBGStreamingStatistics)
- pass 1: verifier `FIXED`  (commit `1705cf9d`)

### 2026-10-15-documentation-review-01-054  The package's nouns "mergeable" and "power-sum quadrature readouts" are used before any clause says what they mean.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-GBGStreamingStatistics`, decided by an agent.

"Mergeable" names the package's property that two summaries of two streams add to the summary of the joined stream; the page says only "the mergeable summaries" and never that summaries combine at all. "Power-sum" presupposes that the summary stores per-bin sums of powers of the samples, which the page never states (it says only "a fixed-size summary"), and "readout" is the package's word for an estimator computed from the summary; the parenthetical lists examples but defines neither noun.

**Fix.** Bind each at first use in a clause: "the mergeable summaries (two summaries of two streams combine into the summary of both)", and "the readouts computed by quadrature over the stored per-bin power sums (moments, expectations, and pushforward densities)" — or link each to the docstring that defines it.

- `src/GBGStreamingStatistics.jl`:22  (GBGStreamingStatistics)
- `src/GBGStreamingStatistics.jl`:26  (GBGStreamingStatistics)
- pass 1: verifier `FIXED`  (commit `1705cf9d`)

### 2026-10-15-documentation-review-01-055  Symbols bound nowhere on the page: the map T and the measure \hat\mu_n in 'the canonical summary T(\hat\mu_n)', and M_4 / M_6 (presumably the central moments of order 4 and 6, but the page never says so).

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-Histograms`, decided by an agent.

A cold reader meets $T(\hat\mu_n)$ in the identity sentence of Histogram1D with no clause saying \hat\mu_n is the empirical measure of the n samples and T the summarizing map; $M_4$ appears three times (prose and the `order` bullet) with M never bound.

**Fix.** Bind in apposition at first use: 'the summary T(\hat\mu_n) of the empirical measure \hat\mu_n of the n samples seen so far' and 'the central moment M_4 of order 4'; or drop T(\hat\mu_n) from the identity sentence.

- `src/Histograms.jl`:21  (Histogram1D)
- `src/Histograms.jl`:34  (Histogram1D)
- `src/Histograms.jl`:47  (Histogram1D)
- `docs/ARCHITECTURE.md`:73  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The spine: empirical measure → summary → post-hoc-smoothed estimator)
- `docs/ARCHITECTURE.md`:87  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The spine: empirical measure → summary → post-hoc-smoothed estimator)

### 2026-10-15-documentation-review-01-056  The package's algebraic vocabulary — 'the additive monoid', 'the cone', 'the hyperstatistical parameter-mixture', 'exponential forgetting' — is used as if known; the page never says that scaling a summary is what makes it a cone, nor what the mixture or the forgetting are.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Histograms`, decided by an agent.

'(the additive monoid)' is a bare parenthetical after the count type; 'the cone the hyperstatistical parameter-mixture and exponential forgetting are built from' names three things defined nowhere on the page and then 'Float-count / cone variant' reuses 'cone' as settled.

**Fix.** One clause each at first use: 'summaries add (merge!) and, on the Float-count variant, scale — a convex cone — which is what re-weighting a family of summaries post hoc (a parameter mixture) or down-weighting old data (exponential forgetting) needs'; or link the term to its home with @ref.

- `src/Histograms.jl`:35  (Histogram1D)
- `src/Histograms.jl`:38  (Histogram1D)
- `src/Histograms.jl`:213  (Histogram2D)
- pass 1: verifier `FIXED`  (commit `1705cf9d`)

### 2026-10-15-documentation-review-01-057  Two named techniques are invoked without a gloss or a link: 'WARPing' (a parenthetical acronym for the binning scheme) and 'the shift formula' (the moment-recentering identity by which moments are recovered).

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Histograms`, decided by an agent.

'(WARPing)' is dropped after 'equal-width' with no expansion or citation; 'recovered post hoc (by the shift formula)' names a formula the page never states or references, so the reader cannot tell what makes the recovery exact.

**Fix.** Expand or cite at first use ('equal-width bins — the WARP binning of Härdle & Scott') and give the shift formula in one clause ('re-centering the power-sums from the bin center to the sample mean by the binomial expansion') or drop the parentheticals.

- `src/Histograms.jl`:22  (Histogram1D)
- `src/Histograms.jl`:32  (Histogram1D)
- pass 1: verifier `FIXED`  (commit `1705cf9d`)

### 2026-10-15-documentation-review-01-058  The pair-formation vocabulary of the two-time path — 'source', 'snapshot buffer', 'origin cadence', 'the memoryless surface' — is the consumer's and the package's, not the field's, and the page never says what a source or an origin is.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Histograms`, decided by an agent.

'paired per source' (Histogram2D, and the doctest comment) and 'which sources to sample, … the origin cadence' assume the reader knows the stream is made of sources each sampled at time origins; 'the memoryless surface' uses 'surface' for an API with 'memoryless' unexplained.

**Fix.** Bind once: 'each pair is formed from one source — one trajectory of the stream — at a time origin t and the later time t + τ; how often origins are taken (the cadence) is the caller's choice'; replace 'the memoryless surface' with 'the entry point that stores no sample'.

- `src/Histograms.jl`:211  (Histogram2D)
- `src/Histograms.jl`:346  (observe_pairs!)
- `src/Histograms.jl`:344  (observe_pairs!)
- pass 1: verifier `FIXED`  (commit `1705cf9d`)

### 2026-10-15-documentation-review-01-059  A sentence that does not parse: 'For a two-time quantity X(t) and the response Y(t+τ) are paired per source' lacks the comma or subject that would make 'X(t) and the response' the subject of 'are paired'.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-Histograms`, decided by an agent.

Read aloud, 'For a two-time quantity X(t) and the response Y(t+τ) are paired' first attaches X(t) to 'quantity' and then has no subject for 'are paired'.

**Fix.** 'For a two-time quantity, the value X(t) and the later response Y(t+τ) of one source are paired; the lag τ is an explicit label on the summary.'

- `src/Histograms.jl`:211  (Histogram2D)
- `docs/ARCHITECTURE.md`:281  (GBGStreamingStatistics — Architecture and Rationale / 3. The accumulated summary — the universal lower half)
- `docs/src/tutorial.md`:52  (Getting started / The workflow)
- pass 1: verifier `NOT_FIXED`  (commit `1705cf9d`) — the verifier said FIXED, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-060  Verbless sentences: the Histogram2D 'Same … as Histogram1D.' sentence, the twice-repeated 'Provided by a package extension: …' fragment, and the Fields lead-in of both types have no finite verb.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-Histograms`, decided by an agent.

'Same atomic-scatter accumulation, merge!-by-addition, Float64 accumulation, and count_type = Float64 Float-count / cone variant (…) as Histogram1D.' is a noun phrase; 'Provided by a package extension: available only when the HDF5 package is loaded.' has no subject; 'The summary's whole internal state — … — each reachable by dot (…), for serialization, inspection, or building one by hand:' never reaches a verb.

**Fix.** Supply the verb: 'Accumulation, merge!, the Float64 accumulator and the count_type = Float64 cone variant are as for Histogram1D.'; 'A package extension provides this method; it is available only when HDF5 is loaded.'; 'The fields hold the summary's whole internal state … ; each is reachable by dot …'.

- `src/Histograms.jl`:212  (Histogram2D)
- `src/Histograms.jl`:412  (write_histogram)
- `src/Histograms.jl`:412  (write_histogram)
- `src/Histograms.jl`:54  (Histogram1D)
- `src/Histograms.jl`:435  (read_histogram)
- pass 1: verifier `FIXED`  (commit `1705cf9d`)

### 2026-10-15-documentation-review-01-061  The two persistence functions state their result both with a signature arrow and with a 'returns' sentence, and the arrow names a thing ('group', 'Histogram1D or Histogram2D') rather than a type fixed by the signature.

Rule: `RULE.DOC.RETURNS`.  
Reported by `withheld@api-Histograms`, decided by an agent.

write_histogram carries '-> group' and its identity sentence ends 'returning that group'; read_histogram carries '-> Histogram1D or Histogram2D' and its sentence says 'to return a Histogram1D or a Histogram2D' — the rule says arrow and sentence never both appear, and an arrow only where the signature fixes the type.

**Fix.** Keep one form per function: drop the arrows (the return depends on the stored type / on HDF5's group type) and let the sentences state the result, or keep a typed arrow and strip the clause from the sentence.

- `src/Histograms.jl`:399  (write_histogram)
- `src/Histograms.jl`:402  (write_histogram)
- `src/Histograms.jl`:425  (read_histogram)
- `src/Histograms.jl`:429  (read_histogram)
- pass 1: verifier `FIXED`  (commit `1705cf9d`)

### 2026-10-15-documentation-review-01-062  The gamma_kernel, beta_kernel and epanechnikov_kernel docstrings use the symbols b, n, n_j, c_j, \hat\sigma and h without binding them; only the AsymmetricKernel docstring, further down the page, names b the bandwidth and n_j / c_j the bin counts and centers.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-Kernels`, decided by an agent.

In gamma_kernel, $b$ first appears inside $\mathrm{Gamma}(\rho(x), b)$ with no gloss; the estimate $\hat f(x) = \tfrac{1}{n}\sum_j n_j\,\mathrm{Gamma}(c_j; \rho(x), b)$ introduces $n$, $n_j$ and $c_j$ cold; the bandwidth default $\hat\sigma n^{-2/5}$ leaves $\hat\sigma$ unnamed. In epanechnikov_kernel, 'the bins within $h$ of the query point' is the first and only $h$ on the page, unbound. beta_kernel inherits the unbound $b$ in $(x/b, (1-x)/b)$. A reader who opens gamma_kernel directly (the docstring is what `?gamma_kernel` shows) has none of AsymmetricKernel's bindings.

**Fix.** Bind each symbol in apposition at its first use in each docstring: 'the bandwidth $b$', 'the sample size $n$', 'the bin counts $n_j$ at the bin centers $c_j$', 'the sample standard deviation $\hat\sigma$', 'the bandwidth $h$'.

- `src/Kernels.jl`:137  (gamma_kernel)
- `src/Kernels.jl`:142  (gamma_kernel)
- `src/Kernels.jl`:145  (gamma_kernel)
- `src/Kernels.jl`:32  (epanechnikov_kernel)
- `src/Kernels.jl`:172  (beta_kernel)
- `docs/ARCHITECTURE.md`:775  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `src/Kernels.jl`:157  (Kernels)
- `src/Kernels.jl`:173  (beta_kernel)

### 2026-10-15-documentation-review-01-063  The package term `summary` (the data object every estimator takes) is used in two code spans and in the runnable example without being defined or constructed anywhere on the page.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Kernels`, decided by an agent.

The page never says what a 'summary' is — the binned histogram the estimators read — yet `bandwidth(summary; rule = :plugin, kernel = gamma_kernel)` and the example's `density(summary, 1.0; kernel = mygamma)` both take it as their first argument. In the example `summary` is a free name: the block cannot be run as written, and a cold reader does not know what to put there. The module docstring, which mentions bins nowhere, does not supply the binned-measure premise either.

**Fix.** In the module docstring, bind the term once ('every estimator reads a `summary`, the binned histogram of the stream, …' with an @ref to its type), and either construct a summary at the top of the example or state in a comment what `summary` is.

- `src/Kernels.jl`:147  (gamma_kernel)
- `src/Kernels.jl`:100  (AsymmetricKernel)
- pass 1: verifier `NOT_FIXED`  (commit `5d603309`) — the verifier said FIXED, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-064  The options `degree`, `adaptive`, `support` and `h` of `density` are named as 'not applying' or 'to override' without being introduced on this page, and `support` collides with the AsymmetricKernel field of the same name, so the 'support renormalization' reads as the field.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Kernels`, decided by an agent.

gamma_kernel calls the kernel 'an alternative to (and usually sharper than) the `support` renormalization' and says 'The options `degree`, `adaptive`, and `support` do not apply'; the page defines none of these options, and the only `support` it defines is AsymmetricKernel's '(lower, upper) tuple' — which does apply, since the kernel carries one. beta_kernel then defers to 'the option restrictions' of gamma_kernel, compounding the reference. A cold reader cannot tell which `support` is meant or what renormalization is being compared against.

**Fix.** Name each option as 'the `degree` / `adaptive` / `support` options of [`density`](@ref …)' at first mention, and give the renormalization a prose name ('the boundary renormalization of a symmetric kernel') that distinguishes it from the `support` field.

!!! warning "The code may be what is wrong"
    The reader believed the documentation and the code disagree, and that the code is the wrong half. That is a question for the owner, not a documentation fix.

- `src/Kernels.jl`:144  (gamma_kernel)
- `src/Kernels.jl`:148  (gamma_kernel)
- `src/Kernels.jl`:174  (beta_kernel)
- `docs/ARCHITECTURE.md`:103  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The algebra of measures — a structural lens)
- `docs/ARCHITECTURE.md`:412  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:563  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- `docs/ARCHITECTURE.md`:775  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `src/Bandwidth.jl`:366  (bandwidth)
- `src/Bandwidth.jl`:367  (bandwidth)
- `src/Kernels.jl`:89  (AsymmetricKernel)
- `src/Kernels.jl`:146  (gamma_kernel)
- `src/Kernels.jl`:175  (beta_kernel)
- `src/MarginalDensity.jl`:451  (density)
- pass 1: verifier `FIXED`  (commit `5d603309`)

### 2026-10-15-documentation-review-01-065  'post hoc and mergeable' asserts two package properties of the estimate that the page never defines.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Kernels`, decided by an agent.

The sentence closes the `pdf` argument bullet with '— post hoc and mergeable.' Neither term is bound on the page: 'mergeable' (two summaries' estimates combine by adding the summaries) and 'post hoc' (the kernel is chosen after accumulation) are the package's architecture vocabulary, not the field's; a first-time reader sees two adjectives with no referent.

**Fix.** Either drop the clause from the kernel's field bullet (it belongs to the summary's own docstring) or gloss it: 'computed after accumulation from the stored counts, so estimates from merged summaries agree'.

- `src/Kernels.jl`:82  (AsymmetricKernel)
- pass 1: verifier `FIXED`  (commit `5d603309`)

### 2026-10-15-documentation-review-01-066  Every identity sentence on the page — the module's and all five docstrings' — is a noun phrase with no main verb, held together by em-dashes and colons.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-Kernels`, decided by an agent.

'The smoothing kernels shared across the estimators, in two families.' has no verb. 'The standard Gaussian smoothing kernel $K(u) = …$ (…) — the density default.' has none. The Epanechnikov, gamma, beta and AsymmetricKernel openers follow the same pattern: 'The gamma kernel (Chen 2000) for a density on the half-line $[0,\infty)$ — a varying-shape asymmetric kernel: at the query point…' — the only finite verb ('it is') sits in a subordinate clause after the colon, and the main clause never predicates anything. Read aloud, none is a sentence.

**Fix.** Give each opener a verb: 'The smoothing kernels shared across the estimators come in two families.', '`gaussian_kernel` is the standard Gaussian kernel $K(u) = …$', 'The gamma kernel (Chen 2000) estimates a density on the half-line …'.

- `docs/src/api/Kernels.md`:3  (Kernels)
- `src/Kernels.jl`:19  (gaussian_kernel)
- `src/Kernels.jl`:29  (epanechnikov_kernel)
- `src/Kernels.jl`:136  (gamma_kernel)
- `src/Kernels.jl`:170  (beta_kernel)
- `src/Kernels.jl`:63  (AsymmetricKernel)
- `docs/ARCHITECTURE.md`:1014  (GBGStreamingStatistics — Architecture and Rationale / References)
- `docs/src/api/Kernels.md`:6  (Kernels)
- pass 1: verifier `FIXED`  (commit `5d603309`)

### 2026-10-15-documentation-review-01-067  beta_kernel describes its edge correction as 'a $\sqrt{\cdot}$ correction', a placeholder standing in for a formula the page never states.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-Kernels`, decided by an agent.

gamma_kernel spells its edge shape out ('$(x/b)^2/4 + 1$ within $2b$ of the edge'), so the reader expects the same of the beta kernel, but gets '$\sqrt{ \cdot }$' with nothing inside: no symbol is bound, no formula given, and the reference to 'Chen's $f_2$' is the only way to learn what the correction is.

**Fix.** State the edge shape explicitly — the shape parameter Chen's $f_2$ substitutes near each edge — or say in prose which argument is replaced and by what.

- `src/Kernels.jl`:173  (beta_kernel)
- pass 1: verifier `NOT_FIXED`  (commit `5d603309`)

### 2026-10-15-documentation-review-01-068  Identity sentences across the page are verbless noun phrases (the page lead, the section lead, the Estimator. blocks, the paired reader, every parametric start and the type) — none reads as a sentence aloud.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-MarginalDensity`, decided by an agent.

Each opens the docstring as a noun phrase with modifiers and no finite verb: 'A log-normal ParametricStart on (0,∞), parameters by method of moments (…)'; 'A parametric start for the Hjort–Glad semiparametric density (Hjort & Glad 1995).'; 'Estimator. The smoothed estimate … from the accumulated summary — a binned local smoother …'; 'The density estimate and its pointwise density_standard_error, from a single option specification — so …' (only the subordinate clause carries a verb).

**Fix.** Give each a verb: 'normal_start is the Gaussian ParametricStart on …'; 'ParametricStart is a parametric start for …'; 'The estimator is the smoothed estimate …'; 'density_with_standard_error returns the density estimate and its pointwise standard error from a single option specification, so …'.

- `src/MarginalDensity.jl`:104  (lognormal_start)
- `src/MarginalDensity.jl`:36  (ParametricStart)
- `src/MarginalDensity.jl`:416  (density)
- `src/MarginalDensity.jl`:775  (density_with_standard_error)
- `docs/src/api/MarginalDensity.md`:4  (MarginalDensity)
- `docs/src/api/MarginalDensity.md`:19  (MarginalDensity / Parametric starts)
- `docs/ARCHITECTURE.md`:327  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.1 Marginal density ``f_X``)
- `src/Bandwidth.jl`:317  (bandwidth)
- `src/MarginalDensity.jl`:692  (density_standard_error)
- `src/MarginalDensity.jl`:695  (density_standard_error)
- `src/MarginalDensity.jl`:698  (density_standard_error)
- `src/MarginalDensity.jl`:702  (density_standard_error)
- pass 1: verifier `NOT_FIXED`  (commit `0f1053e0`) — the verifier said LINE_ONLY, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-069  A clause after a semicolon has no subject.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-MarginalDensity`, decided by an agent.

'Pass `start = exponential_start`; pairs with the `degree = 1` log-linear fit …' — 'pairs' has no subject; the imperative 'Pass' cannot supply one.

**Fix.** 'Pass `start = exponential_start`; it pairs with the `degree = 1` log-linear fit as the two exponential-tail tools.'

- `src/MarginalDensity.jl`:132  (exponential_start)
- pass 1: verifier `FIXED`  (commit `0f1053e0`)

### 2026-10-15-documentation-review-01-070  Symbols in density's formulas are bound nowhere on the page: $h_j$ and $n$ in the degree-0 formula, $\hat\beta_0$ in the local-log readout, $a$, $b$ and $K_h$ in the support renormalization, and $m$, $s$ in the log-normal moment fit.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-MarginalDensity`, decided by an agent.

The degree-0 KDE is written with a per-bin $h_j$ although `h` is one number and `adaptive` is a later bullet (the phrase 'the sample-point $h_j$ per bin' arrives two bullets on and does not say what $h_j$ is under the default); $\hat\beta_0$ is never named as the fitted intercept; the integral's limits $a$, $b$ are never tied to `(lower, upper)` and $K_h$ is never written out; $m$ and $s$ are never named as the log-normal's log-scale mean and standard deviation.

**Fix.** Write the default formula with $h$ and say 'with $h_j = h$ unless `adaptive`'; bind 'the fitted intercept $\hat\beta_0$'; write $\int_{\text{lower}}^{\text{upper}} K((x-u)/h)\,\mathrm{d}u/h$ or bind $K_h(u) = K(u/h)/h$; name $m$, $s$ as the log-mean and log-standard-deviation; bind $n$ as the total count on first use.

- `src/MarginalDensity.jl`:437  (density)
- `src/MarginalDensity.jl`:439  (density)
- `src/MarginalDensity.jl`:462  (density)
- `src/MarginalDensity.jl`:104  (lognormal_start)

### 2026-10-15-documentation-review-01-071  The symbol $b$ is bound twice with different meanings in one docstring: the asymmetric kernel's scale, then the upper limit of the support integral.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-MarginalDensity`, decided by an agent.

'it takes its own scale $b$ in place of `h`' (kernel bullet, repeated in the h bullet) and two bullets later '$\int_a^b K_h(x-u)\,\mathrm{d}u$', where $b$ is the upper support bound.

**Fix.** Use the tuple's own names for the limits ($\int_{\text{lower}}^{\text{upper}}$) or rename the asymmetric scale.

- `src/MarginalDensity.jl`:435  (density)
- `src/MarginalDensity.jl`:462  (density)

### 2026-10-15-documentation-review-01-072  The estimand calls $X$ 'the conditioning variable' on a page about a marginal density of a one-dimensional summary, where nothing is conditioned on; the phrase assumes the package's joint/conditional context.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-MarginalDensity`, decided by an agent.

'The marginal (Lebesgue) density $f_X = \mathrm{d}\mu_X/\mathrm{d}\lambda$ of the conditioning variable $X$' — a cold reader asks: conditioning what? No conditional estimator appears on this page.

**Fix.** 'of the variate $X$' (or 'of the observed variable $X$').

- `src/MarginalDensity.jl`:413  (density)
- pass 1: verifier `NOT_FIXED`  (commit `0f1053e0`) — the verifier said FIXED, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-073  Package terms of art are used without definition: 'power-sum' ($S_1$) and 'mergeable'.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-MarginalDensity`, decided by an agent.

'from the stored first power-sum $S_1$' — a reader who has not met the summary does not know that $S_1^{(j)}$ is the sum of the samples' offsets from the bin center; 'keeps the whole estimate mergeable' — nothing on the page says what merging is (two summaries added).

**Fix.** Gloss in a clause: 'the stored first power-sum $S_1^{(j)}$, the sum over the bin's samples of their offsets from the bin center'; 'mergeable — two summaries' estimates combine by adding the summaries, since $\hat\theta$ is …'.

- `src/MarginalDensity.jl`:484  (density)
- `src/MarginalDensity.jl`:59  (ParametricStart)
- pass 1: verifier `FIXED`  (commit `0f1053e0`)

### 2026-10-15-documentation-review-01-074  'WARPing' is used in the Estimator sentence before its gloss, which only arrives in the `degree` bullet.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-MarginalDensity`, decided by an agent.

'a KDE / WARPing, or a local-log-polynomial fit' precedes '(weighted averaging of rounded points)' by a whole paragraph.

**Fix.** Move the gloss to the first use: 'a binned KDE (WARPing, weighted averaging of rounded points)'.

- `src/MarginalDensity.jl`:418  (density)
- pass 1: verifier `FIXED`  (commit `0f1053e0`)

### 2026-10-15-documentation-review-01-075  Mathematics is set in code spans: the returned NamedTuple is written with the hatted symbol f̂(x) and SE(x), which are not Julia, and the signature comments carry formulas.

Rule: `RULE.TEX.CODE_SPAN_JULIA`.  
Reported by `withheld@api-MarginalDensity`, decided by an agent.

'`(value = f̂(x), standard_error = SE(x))`' — f̂(x) and SE(x) are mathematics, not code that could be typed; the signature blocks end in '# ≡ f_X' and '# ≡ SE[f̂(x)]', formulas in Unicode glyphs inside code.

**Fix.** 'the `NamedTuple` `(value, standard_error)` with `value` $= \hat f(x)$ and `standard_error` $= \mathrm{SE}(x)$'; drop the formula comments from the signatures (the Estimand sentence already states them).

!!! warning "The code may be what is wrong"
    The reader believed the documentation and the code disagree, and that the code is the wrong half. That is a question for the owner, not a documentation fix.

- `src/MarginalDensity.jl`:777  (density_with_standard_error)
- `src/MarginalDensity.jl`:688  (density_standard_error)
- pass 1: verifier `FIXED`  (commit `0f1053e0`)

### 2026-10-15-documentation-review-01-076  The Kernels decision-log entry opens with a sentence that has no main verb and an unclosed parenthesis: the subject 'The varying-shape boundary kernels (... Chen 1999,' runs straight into 'moved from MarginalDensity to Kernels, joining ...' with the citation cut off and no finite verb.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@architecture`, decided by an agent.

Rendered text reads: 'The varying-shape boundary kernels (`gamma_kernel`, `beta_kernel`, the `AsymmetricKernel` type; Chen 1999, moved from `MarginalDensity` to `Kernels`, joining the symmetric `gaussian_kernel` / `epanechnikov_kernel`.' The parenthesis opened before `gamma_kernel` never closes, the second Chen year is missing, and the sentence cannot be read aloud as a statement.

**Fix.** Close the parenthesis and give the sentence its verb: '... the `AsymmetricKernel` type; Chen 1999, 2000) were moved from `MarginalDensity` to `Kernels`, joining ...'.

- `docs/ARCHITECTURE.md`:861  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- pass 1: verifier `NOT_FIXED`  (commit `542c6ddf`)

### 2026-10-15-documentation-review-01-077  The bandwidth section's 'Estimand.' and 'Estimator.' labels are followed by verbless noun phrases, each starting lowercase after a full stop, so the two 'sentences' never say what the estimand or estimator is.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@architecture`, decided by an agent.

'Estimand. the AMISE-optimal bandwidth (AMISE = asymptotic mean integrated squared error) — the $h$ asymptotically minimising ...' and 'Estimator. for a marginal, the kernel- and degree-aware normal-reference rule $h = ...$, plus the data-driven ...' — long appositive lists with no finite verb; the other sections' labels at least lead into an equation.

**Fix.** Give each a verb: 'Estimand. The estimand is the AMISE-optimal bandwidth ...'; 'Estimator. For a marginal the package uses the kernel- and degree-aware normal-reference rule ...'.

- `docs/ARCHITECTURE.md`:405  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:407  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- pass 1: verifier `NOT_FIXED`  (commit `542c6ddf`) — the verifier said FIXED, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-078  Several mathematical symbols are used before the page binds them, or are never bound: $S_p$ (§2 lens, bound only in §3), $\ell(x)$ in the conditional standard error (never bound as the leverage), $\delta(K)$ (§4.5, bound only in the §11 decision log), $\beta_0$ (§4.2, bound only in §11), $L$ and $N$ in the size estimates, and $D^{(n)}(x)$ (§6, bound only in §11).

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@architecture`, decided by an agent.

§2: 'the Taylor contraction of the stored $S_p$' precedes §3's '$S_p = \sum_{i\in\text{bin}}(X_i - c_j)^p$'. §4.2: '$\sqrt{\ell(x) \mathrm{Var}(g(Y)\mid x)}$' with $\ell(x)$ unstated (and $\ell_j$ already means the bin location in §4.1). §4.5: 'the eval kernel a clean $\delta(K)$ rescale' before §11's '$\delta(K) = [R(K)/\mu_2(K)^2]^{1/5}$'. §4.2: 'the $\beta_0$ of one `LocalPolynomial.local_fit`' before §11's '(\beta_0, \beta_1, \dots)'. §3: '($L \times$ bins $\times$ moments, a few MB)' and §5: 'over $N^3\times$steps' with $L$, $N$ unbound. §6: '$D^{(n)}(x)$' before §11's '$D^{(1)} = M_1/\tau$'.

**Fix.** Bind each at first use in apposition: 'the stored per-bin power-sums $S_p$ (§3)'; 'the local-regression leverage $\ell(x) = \sum_j N_j w_j^2$'; 'the kernel's canonical factor $\delta(K) = [R(K)/\mu_2(K)^2]^{1/5}$'; 'the intercept $\beta_0$ of the local-polynomial fit'; '$L$ lags', '$N^3$ grid points per step'; 'the Kramers–Moyal coefficients $D^{(n)}(x)$'.

- `docs/ARCHITECTURE.md`:128  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The algebra of measures — a structural lens)
- `docs/ARCHITECTURE.md`:384  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.2 Conditional expectation ``m(x)`` and variance ``\sigma^2(x)``)
- `docs/ARCHITECTURE.md`:413  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:373  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.2 Conditional expectation ``m(x)`` and variance ``\sigma^2(x)``)
- `docs/ARCHITECTURE.md`:299  (GBGStreamingStatistics — Architecture and Rationale / 3. The accumulated summary — the universal lower half)
- `docs/ARCHITECTURE.md`:492  (GBGStreamingStatistics — Architecture and Rationale / 5. Backend-agnosticism and the hot-loop contract)
- `docs/ARCHITECTURE.md`:530  (GBGStreamingStatistics — Architecture and Rationale / 6. Invariants — deliberate absences are statements)
- `docs/ARCHITECTURE.md`:385  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.2 Conditional expectation ``m(x)`` and variance ``\sigma^2(x)``)
- `docs/ARCHITECTURE.md`:725  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:734  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:740  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:755  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `src/Bandwidth.jl`:334  (bandwidth)
- `src/Bandwidth.jl`:335  (bandwidth)
- `src/Bandwidth.jl`:349  (bandwidth)
- `src/Bandwidth.jl`:359  (bandwidth)
- `src/Bandwidth.jl`:362  (bandwidth)
- `src/Bandwidth.jl`:401  (bandwidth)
- `docs/ARCHITECTURE.md`:300  (GBGStreamingStatistics — Architecture and Rationale / 3. The accumulated summary — the universal lower half)
- `docs/ARCHITECTURE.md`:493  (GBGStreamingStatistics — Architecture and Rationale / 5. Backend-agnosticism and the hot-loop contract)
- `docs/ARCHITECTURE.md`:156  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / Estimand families (the operation taxonomy))
- `docs/ARCHITECTURE.md`:531  (GBGStreamingStatistics — Architecture and Rationale / 6. Invariants — deliberate absences are statements)

### 2026-10-15-documentation-review-01-079  The package term 'jet' (and 'jet readout', 'jet-smoother') is used in §2 as if known; its meaning — the vector of local Taylor coefficients returned by the local-polynomial fit — is stated only in the §11 decision log.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@architecture`, decided by an agent.

§2: 'with `LocalPolynomial.local_fit` the jet readout on it' and 'The jet `local_fit` is the kernel-weighting family's primitive ... a readout with no jet at all'; the definition 'return the local Taylor coefficients $(\beta_0, \beta_1, \dots)$' appears only in §11.

**Fix.** At the §2 first use, define in a clause: 'the jet — the local Taylor coefficients $(\beta_0, \beta_1, \dots)$ of the kernel-weighted polynomial fit — read out by `LocalPolynomial.local_fit`'.

- `docs/ARCHITECTURE.md`:127  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The algebra of measures — a structural lens)
- `docs/ARCHITECTURE.md`:129  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The algebra of measures — a structural lens)
- pass 1: verifier `FIXED`  (commit `542c6ddf`)

### 2026-10-15-documentation-review-01-080  'Channel' ('response channel', 'heavy-tailed channels', 'kilobytes per channel', 'a channel's modal-bin contention') is a package term the page never defines; a first-time reader cannot tell whether a channel is a response variable, a summary, or a consumer's data stream.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@architecture`, decided by an agent.

§2: 'Accumulation knows nothing of bandwidth, kernel, or which response channel will be asked for' and 'every $h$, every channel, every secondary smoothing'; §3: 'raise $P$ for the heavy-tailed channels', 'kilobytes per channel'; §5: 'should a channel's modal-bin contention dominate' — no clause binds the noun anywhere.

**Fix.** Define once at the §2 first use: 'which response channel — which response variable $Y$, each with its own summary — will be asked for', and use that meaning consistently.

- `docs/ARCHITECTURE.md`:93  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The spine: empirical measure → summary → post-hoc-smoothed estimator)
- `docs/ARCHITECTURE.md`:223  (GBGStreamingStatistics — Architecture and Rationale / 3. The accumulated summary — the universal lower half)
- pass 1: verifier `NOT_FIXED`  (commit `542c6ddf`) — the verifier said FIXED, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-081  Abbreviations are used before or without their expansion: 'OU' in §3 (Ornstein–Uhlenbeck is spelled out only in §7), 'ULI' (never expanded), 'SFI' (never bound to Stochastic Force Inference), 'NW' (§11), 'MD' (§8), and 'the Pawula $Q$' in §3 ($Q$ is defined only in §11).

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@architecture`, decided by an agent.

§3: 'for the OU oracle both marginals are the same stationary Gaussian'; §2: 'à la Stochastic Force Inference / ULI' and §4.6 'the underdamped ULI'; §8/§10: 'SFI / ULI road', 'include SFI/ULI'; §11: 'count sum, NW ratio, jet'; §8: 'MD's van Hove function'; §3: 'the conditional increment moments of the Pawula $Q$' precede §11's '$Q = \mu_4/(3\mu_2^2)$'.

**Fix.** Expand each at its first use: 'the Ornstein–Uhlenbeck (OU) oracle', 'Stochastic Force Inference (SFI) / Underdamped Langevin Inference (ULI)', 'the Nadaraya–Watson (NW) ratio', 'molecular dynamics (MD)', and state the Pawula quotient $Q = \mu_4/(3\mu_2^2)$ where it first appears in §3.

- `docs/ARCHITECTURE.md`:282  (GBGStreamingStatistics — Architecture and Rationale / 3. The accumulated summary — the universal lower half)
- `docs/ARCHITECTURE.md`:162  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / Estimand families (the operation taxonomy))
- `docs/ARCHITECTURE.md`:606  (GBGStreamingStatistics — Architecture and Rationale / 8. Proposed repository layout (mirrors PF))
- `docs/ARCHITECTURE.md`:878  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:225  (GBGStreamingStatistics — Architecture and Rationale / 3. The accumulated summary — the universal lower half)
- pass 1: verifier `LINE_ONLY`  (commit `542c6ddf`)

### 2026-10-15-documentation-review-01-082  The §4.1 option-pair table sets the inequality 'degree ≥ 1' inside a code span in three rows; an inequality is mathematics and the ≥ glyph is not Julia code.

Rule: `RULE.TEX.CODE_SPAN_JULIA`.  
Reported by `withheld@architecture`, decided by an agent.

Rows '`support` × `degree ≥ 1`', 'asymmetric kernel × `degree ≥ 1`' and '`start` × `degree ≥ 1`' each put a ≥ inequality in backticks, while the same page writes the comparable condition as math elsewhere ('$\texttt{degree} \ge 1$').

**Fix.** Set the condition as math with the option name in \texttt: `support` × $\texttt{degree} \ge 1$ (or the Julia-valid `degree >= 1` if a literal is wanted).

- `docs/ARCHITECTURE.md`:354  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.1 Marginal density ``f_X``)
- `docs/ARCHITECTURE.md`:356  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.1 Marginal density ``f_X``)
- `docs/ARCHITECTURE.md`:357  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.1 Marginal density ``f_X``)
- pass 1: verifier `FIXED`  (commit `542c6ddf`)

### 2026-10-15-documentation-review-01-083  'Stages 1–3' is a project-internal label the page never explains; a reader cannot tell what the stages were or where they are described.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@architecture`, decided by an agent.

§11: 'A framing decision, landed in the code over Stages 1–3.' and 'Landed (Stages 1–3): the shared `LocalPolynomial.local_fit` ...' — no sentence on the page says what a Stage is.

**Fix.** Either drop the label ('landed in the code in three steps: ...') or name the three steps in the entry.

- `docs/ARCHITECTURE.md`:827  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- pass 1: verifier `NOT_FIXED`  (commit `542c6ddf`) — the verifier said LINE_ONLY, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-084  Sentences with no finite verb: a module entry runs 'With ... to persist ...' as a sentence, a clause is reduced to the past participle 'Read off ...', and a parenthetical 'Silverman the Gaussian/degree-0 case' drops its verb.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@home`, decided by an agent.

RULE.DOC.GRAMMAR: every sentence has a subject and a verb; concise means few words, not missing verbs. 'With `observe_pairs!` for two-time statistics, and ... `write_histogram` / `read_histogram` to persist a summary to an HDF5 group.' has no main verb; 'Read off the stored power-sums;' has no subject and reads aloud as an imperative; 'Silverman the Gaussian/degree-0 case' lacks 'is' or 'being'.

**Fix.** Supply the verbs: 'It also offers `observe_pairs!` for two-time statistics and, when `HDF5` is loaded (a package extension), `write_histogram` / `read_histogram`, which persist a summary to an HDF5 group.'; 'All are read off the stored power-sums; ...'; '(Silverman's rule is the Gaussian, degree-0 case)'.

- `docs/src/index.md`:84  (GBGStreamingStatistics.jl / Main components)
- `docs/src/index.md`:95  (GBGStreamingStatistics.jl / Main components)
- `docs/src/index.md`:115  (GBGStreamingStatistics.jl / Main components)
- pass 1: verifier `FIXED`  (commit `8a9d122e`)

### 2026-10-15-documentation-review-01-085  Lists that do not read aloud as sentences: the Status sentence places its 'and' before the penultimate item so the final subject dangles; the MarginalDensity entry has an 'and' mid-list and uses 'the local-log-polynomial' as a noun with no head; the ConditionalExpectations entry stacks two parentheticals so the method clause attaches to the wrong item.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@home`, decided by an agent.

RULE.DOC.GRAMMAR: compact prose that still reads as sentences aloud. 'and the conditional mean and variance (...), the conditional and transition densities (...) are implemented' — the conjunction precedes the wrong item; 'the local-log-polynomial (`degree`), and the Abramson variable bandwidth (`adaptive`, for tails), the asymmetric' — a second 'and' is needed and 'local-log-polynomial' needs 'estimator'; '`conditional_pushforward_density` ($f_{g(Y)\mid X}$) (local-polynomial regression in $x$ of a $y$-quadrature response) — each with its standard error' — the second parenthetical describes all four estimators but is appended to the last.

**Fix.** Move the 'and' before the last item of the Status list ('..., the conditional mean and variance (...), and the conditional and transition densities (...) are implemented'); write 'the local-log-polynomial estimator (`degree`)' and move 'and' to the last item; move the method clause to a sentence of its own ('All are local-polynomial regressions in $x$ of a $y$-quadrature response.').

- `docs/src/index.md`:35  (GBGStreamingStatistics.jl)
- `docs/src/index.md`:98  (GBGStreamingStatistics.jl / Main components)
- `docs/src/index.md`:109  (GBGStreamingStatistics.jl / Main components)
- `docs/ARCHITECTURE.md`:922  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- pass 1: verifier `NOT_FIXED`  (commit `8a9d122e`) — the verifier said FIXED, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-086  The package's terms 'hyperstatistics', 'forgetting' and 'Float-count cone variant' are used without being defined on the page.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@home`, decided by an agent.

RULE.DOC.TERMS_DEFINED: a term of art is defined in a clause, or linked to its definition, before its first use on the same page. 'Hyperstatistics' is not the field's word (the literature says 'superstatistics'); 'forgetting' (exponential down-weighting of old blocks) and 'Float-count cone variant' (a summary whose counts are floating-point so they can be scaled) are the package's own names, and the page nowhere says what they are.

**Fix.** Gloss each at first use: '... or weight and superpose for hyperstatistics (a post-hoc mixture over summaries with different parameters, the superstatistics of Beck & Cohen)'; '... with a variant whose counts are floating-point (the cone variant), so summaries can be re-weighted post hoc for hyperstatistics or exponentially forgotten'.

- `docs/src/index.md`:11  (GBGStreamingStatistics.jl)
- `docs/src/index.md`:33  (GBGStreamingStatistics.jl)
- pass 1: verifier `LINE_ONLY`  (commit `8a9d122e`)

### 2026-10-15-documentation-review-01-087  The package's vocabulary 'readout' / 'reader' and 'quadrature' (as in 'power-sum quadrature readouts', 'paired reader', '$y$-quadrature response') is used before any clause says what it means.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@home`, decided by an agent.

RULE.DOC.TERMS_DEFINED. A 'readout' or 'reader' here is a function that evaluates an estimand from the stored summary, and 'quadrature' means a contraction of the stored power-sums with polynomial coefficients; neither is the field's standard sense, and the first use ('plus the power-sum quadrature readouts') precedes the only hint ('Read off the stored power-sums') by a page. 'Paired reader' and '$y$-quadrature response' have no gloss at all.

**Fix.** Bind once at first use: 'the power-sum quadrature readouts (functions that evaluate an expectation by contracting the stored per-bin power-sums with the coefficients of a polynomial fit to $g$)'; replace 'paired reader' with 'a `_with_standard_error` reader returning the estimate and its standard error together'; gloss '$y$-quadrature response' as 'a per-bin expectation of $g(Y)$ read off the $y$-power-sums'.

- `docs/src/index.md`:40  (GBGStreamingStatistics.jl)
- `docs/src/index.md`:92  (GBGStreamingStatistics.jl / Main components)
- `docs/src/index.md`:113  (GBGStreamingStatistics.jl / Main components)
- `docs/src/index.md`:110  (GBGStreamingStatistics.jl / Main components)
- pass 1: verifier `FIXED`  (commit `8a9d122e`)

### 2026-10-15-documentation-review-01-088  'Channel' is listed among the things computed post hoc, but the page never says what a channel is.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@home`, decided by an agent.

RULE.DOC.TERMS_DEFINED. 'every bandwidth, channel, and conditioning value is computed post hoc' — bandwidth and conditioning value are the field's words; 'channel' is the package's (presumably which variate or component is being estimated) and is used only here, undefined.

**Fix.** Replace with the field's word or gloss it: 'every bandwidth, estimand, and conditioning value' or 'every bandwidth, channel (which quantity is read off the summary), and conditioning value'.

- `docs/src/index.md`:21  (GBGStreamingStatistics.jl)
- pass 1: verifier `FIXED`  (commit `8a9d122e`)

### 2026-10-15-documentation-review-01-089  'Consumer' is used in the package's sense (a downstream package or code that calls this library) without definition.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@home`, decided by an agent.

RULE.DOC.TERMS_DEFINED. 'restorable when a consumer needs them' and 'shared verbatim by every consumer' — a scientist on a first read does not know whether a consumer is a user, a program, or a package.

**Fix.** Say it once: 'shared verbatim by every consumer (a package or program built on this one)' and reuse the word thereafter.

- `docs/src/index.md`:44  (GBGStreamingStatistics.jl)
- `docs/src/index.md`:74  (GBGStreamingStatistics.jl / Main components)
- pass 1: verifier `FIXED`  (commit `8a9d122e`)

### 2026-10-15-documentation-review-01-090  Symbols appear unbound: the lag $\tau$ in $p_\tau(y\mid x)$, the response $Y$ and conditioning value $x$ in $\mathbb{E}[g(Y)\mid x]$ (and the pair $(X, Y)$ generally), and the kernel $K$ in $\int K = 1$.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@home`, decided by an agent.

RULE.DOC.SELF_CONTAINED: every symbol is bound in plain language at first use. The page's only bound symbols are $\hat\mu_n$ (the empirical measure), $g$ (the smooth functional), $m(x)$, $\sigma^2(x)$ and $h^\star(x)$ by apposition; $\tau$ is never said to be the time lag between the two members of a pair, $X$ and $Y$ are never said to be the conditioning variate and the response, and $K$ is never said to be the kernel (the sentence names 'smoothing kernels', then writes '$\int K = 1$' of 'any callable').

**Fix.** Bind each at first use: 'two-time transition density `transition_density` ($p_\tau(y\mid x)$, the density of the response $Y$ at lag $\tau$ after the conditioning value $x$)'; 'the conditional expectation ($\mathbb{E}[g(Y)\mid x]$, $Y$ the response and $x$ the conditioning value)'; 'any callable $K$ with $\int K = 1$'.

- `docs/src/index.md`:112  (GBGStreamingStatistics.jl / Main components)
- `docs/src/index.md`:103  (GBGStreamingStatistics.jl / Main components)
- `docs/src/index.md`:88  (GBGStreamingStatistics.jl / Main components)
- `docs/ARCHITECTURE.md`:391  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.3 Conditional and transition density)
- `docs/ARCHITECTURE.md`:545  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- `docs/ARCHITECTURE.md`:562  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- `docs/ARCHITECTURE.md`:719  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/api/ConditionalDensity.md`:3  (ConditionalDensity)
- `docs/src/index.md`:113  (GBGStreamingStatistics.jl / Main components)
- `docs/ARCHITECTURE.md`:380  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.2 Conditional expectation ``m(x)`` and variance ``\sigma^2(x)``)
- `docs/ARCHITECTURE.md`:920  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/index.md`:104  (GBGStreamingStatistics.jl / Main components)
- `docs/ARCHITECTURE.md`:327  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.1 Marginal density ``f_X``)
- `docs/ARCHITECTURE.md`:570  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- `docs/ARCHITECTURE.md`:714  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:715  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/api/Kernels.md`:4  (Kernels)
- `src/Bandwidth.jl`:536  (bandwidth)
- `src/ConditionalDensity.jl`:94  (conditional_density)
- `src/ConditionalDensity.jl`:114  (conditional_density)
- `src/ConditionalDensity.jl`:130  (conditional_density)
- `src/ConditionalExpectations.jl`:111  (conditional_expectation)
- `src/ConditionalExpectations.jl`:297  (conditional_mean)
- `src/ConditionalExpectations.jl`:379  (conditional_variance)
- `src/ConditionalExpectations.jl`:519  (conditional_pushforward_density)
- `src/ConditionalExpectations.jl`:538  (conditional_pushforward_density)
- `src/Expectations.jl`:303  (pushforward_density)
- `src/Expectations.jl`:314  (pushforward_density)
- `src/Kernels.jl`:4
- `src/Kernels.jl`:19  (gaussian_kernel)
- `src/Kernels.jl`:30  (epanechnikov_kernel)
- `src/MarginalDensity.jl`:424  (density)
- `src/MarginalDensity.jl`:437  (density)
- `src/MarginalDensity.jl`:700  (density_standard_error)

### 2026-10-15-documentation-review-01-091  Three sentences on the page have no finite verb: the page's opening line, the section's opening line, and the docstring's identity sentence are all noun phrases carried by a dash or a colon.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@internals`, decided by an agent.

Read aloud, 'Shared building blocks below the public estimator API — documented here for development, review, and extension (...)' never says what the building blocks do (the verb 'are' is missing before 'documented'); 'The single linear primitive on which the local smoothers are built (the formalization: ...)' is a noun phrase with a subordinate clause and no main verb; the docstring opens 'The local-polynomial jet of a binned 1-D measure at the query point `x`: the kernel-weighted least-squares fit of ... returning the coefficient vector' — a noun phrase, a colon, another noun phrase, no finite verb.

**Fix.** Give each a verb: 'These are the shared building blocks below the public estimator API, documented here for ...'; 'This is the single linear primitive on which the local smoothers are built'; '`local_fit` computes the local-polynomial jet of a binned 1-D measure at the query point `x`: it fits ... and returns the coefficient vector ...'.

!!! warning "The code may be what is wrong"
    The reader believed the documentation and the code disagree, and that the code is the wrong half. That is a question for the owner, not a documentation fix.

- `docs/src/internals.md`:3  (Internals)
- `docs/src/internals.md`:16  (Internals / The local-polynomial jet — `LocalPolynomial.local_fit`)
- `src/LocalPolynomial.jl`:66  (local_fit)
- pass 1: verifier `FIXED`  (commit `32066856`)

### 2026-10-15-documentation-review-01-092  The letter \ell is bound to two different things inside one display: the bin location \ell_j and the moment order (the subscript of T_\ell, P_\ell and the exponent), so the formula reads as if the exponent were a location.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@internals`, decided by an agent.

The docstring binds 'bin locations $\ell_j$', then writes 'T_\ell &= \sum_j K_h(\ell_j - x)  w_j  (\ell_j - x)^\ell' and the normal equations '[ T_{i+\ell} ]_{i,\ell=0}^{d}', where \ell is now a running index from 0 to d. A first-time reader has to work out from the summation structure that the \ell in the exponent is not the \ell of \ell_j.

**Fix.** Rename one of the two: write the locations as $c_j$ (or $x_j$) and keep \ell as the moment order, or keep $\ell_j$ and index the moments by $k$ or $p$ ($T_k$, $(\ell_j - x)^k$, $[T_{i+k}]_{i,k=0}^d$). Use the same choice in the `locations` argument bullet.

!!! warning "The code may be what is wrong"
    The reader believed the documentation and the code disagree, and that the code is the wrong half. That is a question for the owner, not a documentation fix.

- `src/LocalPolynomial.jl`:75  (local_fit)
- `src/LocalPolynomial.jl`:80  (local_fit)

### 2026-10-15-documentation-review-01-093  Symbols are used before the page binds them, or never bound: $\beta_0$ appears in the section prose five lines before the docstring defines the jet $\beta$; $\hat m(x)$ is used as the fitted regression function without ever being named as such; $Y$ appears as the summand of the response total with no sentence saying it is the response variable.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@internals`, decided by an agent.

Reading top to bottom: 'The estimator `conditional_mean` is its $\beta_0$' comes first, and $\beta_0$ is bound only inside the docstring below ('returning the coefficient vector (jet) $\beta = (\beta_0, \dots, \beta_d)$'). Then 'so a regression fit gives $\hat m(x) = \beta_0$ and $\hat m''(x) = 2\beta_2$' — $\hat m$ is nowhere said to be the estimate of the regression function $m(x) = \mathbb{E}[Y \mid X = x]$. In the `response` bullet, '(e.g. $\sum Y$ over the bin)' introduces $Y$ bare.

**Fix.** In the section prose, write 'is its fitted value $\beta_0$ (the leading coefficient of the jet defined below)' or defer the symbol to the docstring. In the docstring, bind $\hat m$ at first use: 'a regression fit of a response $Y$ on $X$ gives the regression-function estimate $\hat m(x) = \beta_0$'. In the `response` bullet, 'the sum $\sum Y$ of the response variable $Y$ over the bin'.

!!! warning "The code may be what is wrong"
    The reader believed the documentation and the code disagree, and that the code is the wrong half. That is a question for the owner, not a documentation fix.

- `docs/src/internals.md`:17  (Internals / The local-polynomial jet — `LocalPolynomial.local_fit`)
- `src/LocalPolynomial.jl`:82  (local_fit)
- `src/LocalPolynomial.jl`:96  (local_fit)
- pass 1: verifier `FIXED`  (commit `32066856`)

### 2026-10-15-documentation-review-01-094  Package nouns are used on the page without a defining clause: 'summary' (the `Histogram` summary types), 'recovered' (raw recovered arrays, recovered centroids), and 'atom' (the `Kernels._kernel_weight` atom), the last naming an internal function the page describes nowhere.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@internals`, decided by an agent.

'which take the `Histogram` summary types; `local_fit` works one level lower, on raw recovered arrays' — a cold reader does not know that a 'summary' is the package's binned accumulator nor what the arrays are 'recovered' from; '(geometric centers or recovered centroids)' repeats 'recovered' with the same gap; 'The weight is the shared `Kernels._kernel_weight` atom (its normalization cancels in the solve)' uses 'atom' in the package's own sense (a shared primitive) and leans on a private name the reader cannot look up on this page.

**Fix.** One clause each: 'the `Histogram1D` / `Histogram2D` summaries — the binned accumulators the package keeps instead of the samples'; 'on the arrays read back from such a summary (bin locations, counts, response totals)'; 'the within-bin centroids read back from the stored first power-sum'; 'The weight is the package's one shared kernel-weighting routine, `Kernels._kernel_weight`, so its normalization cancels in the solve' — or drop the private name and keep only 'the kernel weight's normalization cancels in the solve'.

!!! warning "The code may be what is wrong"
    The reader believed the documentation and the code disagree, and that the code is the wrong half. That is a question for the owner, not a documentation fix.

- `docs/src/internals.md`:11  (Internals)
- `src/LocalPolynomial.jl`:94  (local_fit)
- `src/LocalPolynomial.jl`:83  (local_fit)
- pass 1: verifier `FIXED`  (commit `32066856`)

### 2026-10-15-documentation-review-01-095  'Fit the design itself' uses 'the design' in a sense the page never fixes: the bullet tells the reader to pass `mass` as `response` without saying what that fit estimates.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@internals`, decided by an agent.

The `response` bullet ends 'pass `mass` to fit the design itself.' Earlier, 'design' names the weights $w_j$ and the moments $T_\ell$; a first-time reader cannot tell whether 'fitting the design' means a fit of the counts (a local-polynomial smooth of the bin mass, as in the density use) or something about the design matrix.

**Fix.** Say what it yields: 'pass `mass` here to smooth the bin counts themselves (the jet of the binned mass, as the local-likelihood density does), in which case $\beta_0$ is the local count level rather than a response mean.'

!!! warning "The code may be what is wrong"
    The reader believed the documentation and the code disagree, and that the code is the wrong half. That is a question for the owner, not a documentation fix.

- `src/LocalPolynomial.jl`:96  (local_fit)
- pass 1: verifier `FIXED`  (commit `32066856`)

### 2026-10-15-documentation-review-01-096  The package term 'channel' is used twice and never defined on the page.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@tutorial`, decided by an agent.

Step 2 says accumulation 'commits to no channel, bandwidth, or kernel' and step 3 promises estimates 'for any channel'; bandwidth and kernel are the field's words, 'channel' is not — a statistician reading this page alone cannot say what a channel is (a variate? a component? a response?) or what committing to one would mean.

**Fix.** Bind it in a clause at first use (e.g. 'no channel — which variate of the stream is being summarised —') or replace it with the field's word (variate / response).

- `docs/src/tutorial.md`:27  (Getting started / The workflow)
- `docs/src/tutorial.md`:44  (Getting started / The workflow)
- pass 1: verifier `FIXED`  (commit `32066856`)

### 2026-10-15-documentation-review-01-097  Step 1 calls the variate of a one-dimensional marginal summary a 'conditioning variable', a role the page only introduces in step 4.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@tutorial`, decided by an agent.

At step 1 nothing has been conditioned on anything; the reader asks 'conditioning what?' and finds the (X, Y) pairing three steps later. The same sentence calls Histogram1D 'the canonical summary' without saying canonical with respect to what.

**Fix.** In step 1 say 'of a scalar variate X' (or 'of the variate X'), and either drop 'canonical' or say in a clause why this summary is the reference one.

- `docs/src/tutorial.md`:17  (Getting started / The workflow)
- pass 1: verifier `NOT_FIXED`  (commit `32066856`) — the verifier said FIXED, but its quote is not in the text; read as NOT_FIXED

### 2026-10-15-documentation-review-01-098  'snapshot buffer', 'lag policy' and 'consumer' are package terms used without definition.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@tutorial`, decided by an agent.

Step 4 says 'the snapshot buffer and lag policy stay consumer-side'; the page never says that a buffer of past values must be kept to pair X(t) with Y(t+τ), nor what a 'lag policy' decides, and 'consumer' (also in the code comments 'the consumer produces the batches') is never bound to 'the calling code'.

**Fix.** Bind each in a clause: 'the buffer of past snapshots needed to form the pairs, and the choice of lags, are the calling code's (the consumer's) business'.

- `docs/src/tutorial.md`:54  (Getting started / The workflow)
- pass 1: verifier `FIXED`  (commit `32066856`)

### 2026-10-15-documentation-review-01-099  Keyword arguments in the code are bound nowhere in the prose: `lag`, `hx`, `hy`.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@tutorial`, decided by an agent.

The constructor passes 'lag = 1' right after the prose said the lag policy stays consumer-side, so the reader cannot tell what the summary does with a lag or whether 1 is the τ of the comment 'at lag τ'; 'hx = 0.3, hy = 0.2' are left to be guessed as the bandwidths in x and y by analogy with 'h = bw', which the step-3 comment does bind.

**Fix.** A clause or a code comment per keyword: '# lag records τ in the summary's metadata', '# bandwidths in x and y'.

- `docs/src/tutorial.md`:58  (Getting started / The workflow)
- `docs/src/tutorial.md`:65  (Getting started / The workflow)
- pass 1: verifier `FIXED`  (commit `32066856`)

### 2026-10-15-documentation-review-01-100  Step 4 opens with a verbless fragment where steps 1–3 open with imperative sentences.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@tutorial`, decided by an agent.

'Allocate a summary.', 'Observe the stream.', 'Estimate, post hoc.' are sentences; 'Two-time statistics.' has no verb and breaks the list's pattern.

**Fix.** 'Form the two-time statistics.' or 'Pair the stream for two-time statistics.'

- `docs/src/tutorial.md`:51  (Getting started / The workflow)
- pass 1: verifier `FIXED`  (commit `32066856`)

## Locations that could not be found again — 4

The text each of these named has changed since it was read. They are named rather than
pointed at whatever now sits on the old line, and no fix was attempted at any of them.

## Questions

### 2026-10-15-documentation-review-01-090-code

This finding needs a change only you may make. Will you make it?

Binding the lag, the response, the conditioning value and the kernel symbol on the home page adds inline math spans (the lag symbol, the response symbol, the kernel symbol), which the change guard refuses as a mathematics change (clause 6); the owner decides whether these bindings are added.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-090-code "<option>"`

### 2026-10-15-documentation-review-01-092-code

This finding needs a change only you may make. Will you make it?

Renaming the bin location or the moment index in local_fit's display formula and normal equations rewrites mathematics, which the change guard refuses (clause 6); the owner chooses which letter changes.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-092-code "<option>"`

### 2026-10-15-documentation-review-01-017-code

This finding needs a change only you may make. Will you make it?

Binding the empirical measure, the per-bin power-sum with its defining sum, and the mixture density's parameter spacing inside the accumulation docstrings adds inline math spans, which the change guard refuses as a mathematics change (clause 6); the owner decides whether these bindings are added.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-017-code "<option>"`

### 2026-10-15-documentation-review-01-032-code

This finding needs a change only you may make. Will you make it?

Replacing the approximately-equal sign by the word removes an inline math span, which the change guard refuses as a mathematics change (clause 6); the rate sentence's rewording follows once the owner rules.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-032-code "<option>"`

### 2026-10-15-documentation-review-01-037-code

This finding needs a change only you may make. Will you make it?

Binding the response bin width in the three conditional docstrings adds an inline math span, which the change guard refuses as a mathematics change (clause 6); the owner decides whether the binding is added.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-037-code "<option>"`

### 2026-10-15-documentation-review-01-038-code

This finding needs a change only you may make. Will you make it?

Binding the scaled kernel, the sample count, the marginal density, the kernel moments and the column count in conditional_mean and conditional_variance adds inline math spans, which the change guard refuses as a mathematics change (clause 6); the owner decides whether these bindings are added.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-038-code "<option>"`

### 2026-10-15-documentation-review-01-050-code

This finding needs a change only you may make. Will you make it?

Binding the summary order, the sample count and the scaled kernel in expectation and pushforward_density adds inline math spans, which the change guard refuses as a mathematics change (clause 6); the owner decides whether these bindings are added.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-050-code "<option>"`

### 2026-10-15-documentation-review-01-055-code

This finding needs a change only you may make. Will you make it?

Binding the summary map, the empirical measure and the central moments of order 4 and 6 in the Histogram1D docstring adds inline math spans, which the change guard refuses as a mathematics change (clause 6); the owner decides whether these bindings are added or the identity sentence drops the map.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-055-code "<option>"`

### 2026-10-15-documentation-review-01-062-code

This finding needs a change only you may make. Will you make it?

Binding the bandwidth, the sample size, the bin counts and centers, the sample standard deviation and the symmetric bandwidth in the kernel docstrings adds inline math spans, which the change guard refuses as a mathematics change (clause 6); the owner decides whether these bindings are added.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-062-code "<option>"`

### 2026-10-15-documentation-review-01-003-code

This finding needs a change only you may make. Will you make it?

One symbol per concept across every surface (the bin width, the asymmetric scale, the bin location, the moment index, the four weights) rewrites formulas on every page and docstring, which the change guard refuses as a mathematics change (clause 6), and a notation page is a new page the owner decides on.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-003-code "<option>"`

### 2026-10-15-documentation-review-01-004-code

This finding needs a change only you may make. Will you make it?

Three of the four Fields sections now close their docstrings (Histogram1D, Histogram2D, ParametricStart, committed). AsymmetricKernel's Fields section holds a formula that would move past the Arguments section's formulas, which the change guard refuses as reordered mathematics (clause 6); the owner moves that one or rules it stays.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-004-code "<option>"`

### 2026-10-15-documentation-review-01-008-code

This finding needs a change only you may make. Will you make it?

Making the quadrature's two defaulted arguments keywords changes code (the internal _quad signature and its call sites); the owner decides whether the positional rule reaches internals.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-008-code "<option>"`

### 2026-10-15-documentation-review-01-009-code

This finding needs a change only you may make. Will you make it?

The tutorial, the Architecture layout block and the AsymmetricKernel example are rewrapped (committed under other findings). The remaining long lines sit inside jldoctest blocks (expectation, expectation_with_standard_error, Histogram2D, density, ParametricStart), which the change guard refuses to edit as executed blocks (clause 5); the owner rewraps those or rules them acceptable.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-009-code "<option>"`

### 2026-10-15-documentation-review-01-010-code

This finding needs a change only you may make. Will you make it?

The three signature comments and the Architecture layout block are fixed (committed under other findings). The remaining lines sit inside jldoctest blocks (Histogram1D, Histogram2D), which the change guard refuses to edit as executed blocks (clause 5); the owner rewraps those or rules them acceptable.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-010-code "<option>"`

### 2026-10-15-documentation-review-01-011-code

This finding needs a change only you may make. Will you make it?

Giving the second equation of each display its own row rewrites the display formulas (bandwidth, merge!, conditional_density_standard_error), which the change guard refuses as a mathematics change (clause 6); the owner splits them or records the where-clause as the sanctioned form.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-011-code "<option>"`

### 2026-10-15-documentation-review-01-078-code

This finding needs a change only you may make. Will you make it?

Binding the stored power-sums, the leverage, the canonical kernel factor, the fit intercept, the lag and grid counts and the Kramers–Moyal coefficients at their first use on the Architecture page adds inline math spans, which the change guard refuses as a mathematics change (clause 6); the owner decides whether these bindings are added.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-078-code "<option>"`

### 2026-10-15-documentation-review-01-070-code

This finding needs a change only you may make. Will you make it?

Binding the per-bin bandwidth, the sample count, the fitted intercept, the support limits and the scaled kernel in density's formulas adds inline math spans or rewrites the displayed formulas, which the change guard refuses as a mathematics change (clause 6); the log-normal start's parameters are already named in prose (committed). The owner decides whether these bindings are added.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-070-code "<option>"`

### 2026-10-15-documentation-review-01-071-code

This finding needs a change only you may make. Will you make it?

Writing the support integral's limits as the tuple's own names, or renaming the asymmetric scale, rewrites a formula in density's docstring, which the change guard refuses as a mathematics change (clause 6); the owner chooses which symbol changes.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-071-code "<option>"`

