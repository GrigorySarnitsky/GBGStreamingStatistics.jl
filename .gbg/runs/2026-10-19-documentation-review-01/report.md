# Documentation Review — 2026-10-19-documentation-review-01

!!! info "Written by AI, not regularly reviewed"
    This page is written by AI and is not reviewed regularly by a human.

Run `2026-10-19-documentation-review-01` of `GBGStreamingStatistics`, by GBGCodingTools 2.0.0, from commit `33f22de1`.

It ran to the end. It launched 43 agent(s), with the bound readers, and made 1 of 1 fix pass(es). 20 line(s) went to the work queue, under the heading the end-of-turn hook does not count.

| Stage | Pass | Agents | Cache write | Cache read | Input | Output |
|---|---|---|---|---|---|---|
| find | 0 | 17 | 1233902 | 9718317 | 210 | 45853 |
| verify | 1 | 26 | 1681024 | 10514124 | 376 | 23012 |

Tokens, summed from each agent's transcript by `gbg run cost`.

## Findings — 140

### 2026-10-15-documentation-review-01-003  The notation is not one-to-one and the package has no notation page to settle it: the bin width is `b` in some places and `Δ` in others while `b` is also the asymmetric-kernel scale; `ℓ` is at once the bin location ℓ_j, the moment index in T_ℓ / P_ℓ (both in one formula of `local_fit`) and the leverage ℓ(x); `w` is the design weight w_j, the effective regression weight w_j(x), the AMISE local-scale weight w(x) and the hyperstatistical mixing weight w_k.

Rule: `RULE.DOC.ONE_NOTATION`.  
Reported by `full`, decided by an agent; carried over from an earlier run, which left it open.

`Histogram1D` and ARCHITECTURE §2/§4.1 call the bin width ``b`` (``b \ll h``, ``O(b^2)``); `density`, `expectation`, `conditional_expectation` and ARCHITECTURE §4.1/§4.B/§11 call it ``\Delta``; `density`, `bandwidth`, `AsymmetricKernel`, `gamma_kernel` and ARCHITECTURE §4.5 use ``b`` for the gamma/beta kernel scale, so ARCHITECTURE alone binds `b` to two concepts. `local_fit` writes ``T_\ell = \sum_j K_h(\ell_j - x) w_j (\ell_j - x)^\ell`` with ℓ as both index and location, and `conditional_expectation_standard_error` then uses ``\ell(x)`` for the leverage. `RULE.DOC.ONE_NOTATION`'s `notation_page` is unset, so it names `docs/src/notation.md`, which does not exist.

**Fix.** Pick one symbol per concept (e.g. Δ for the bin width everywhere, b reserved for the asymmetric scale; ℓ_j → c̃_j or x_j for locations and k for the moment index; distinct letters for the four weights) and sweep every surface; record the choice on a notation page.

- `src/Histograms.jl`:21  (Histogram1D)
- `src/MarginalDensity.jl`:490  (density)
- `src/Expectations.jl`:194  (expectation)
- `docs/ARCHITECTURE.md`:138  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / Why a fine histogram is the canonical ``T`` (and the honest caveat))
- `docs/ARCHITECTURE.md`:332  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.1 Marginal density ``f_X``)
- `docs/ARCHITECTURE.md`:418  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `src/MarginalDensity.jl`:439  (density)
- `src/Bandwidth.jl`:367  (bandwidth)
- `src/LocalPolynomial.jl`:76  (local_fit)
- `src/ConditionalExpectations.jl`:166  (conditional_expectation_standard_error)
- `src/Kernels.jl`:85  (AsymmetricKernel)
- `src/Accumulators.jl`:76  (merge!)
- `docs/ARCHITECTURE.md`:833  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:103  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The algebra of measures — a structural lens)
- `docs/ARCHITECTURE.md`:569  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- `docs/ARCHITECTURE.md`:790  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `src/Bandwidth.jl`:367  (bandwidth)
- `src/Kernels.jl`:89  (AsymmetricKernel)
- `src/Kernels.jl`:149  (gamma_kernel)
- `src/Kernels.jl`:179  (beta_kernel)
- `src/MarginalDensity.jl`:455  (density)
- `docs/ARCHITECTURE.md`:140  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / Why a fine histogram is the canonical ``T`` (and the honest caveat))
- `src/Bandwidth.jl`:368  (bandwidth)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-004  In the two struct docstrings the `# Fields` section is not last: it precedes `# Arguments`, a code example and the reference, whereas the home order puts `# References` and then `# Fields` after all prose.

Rule: `RULE.DOC.HOME_ORDER`.  
Reported by `full`, decided by an agent; carried over from an earlier run, which left it open.

`AsymmetricKernel`'s docstring runs identity → `# Fields` → `# Arguments` → julia example → `# References`; `ParametricStart`'s runs identity → `# Fields` → `# Arguments` → jldoctest → 'Reference: N. L. Hjort & I. K. Glad (1995)'. RULE.DOC.HOME_ORDER: 'then `# References`, then `# Fields` last'.

**Fix.** Move each `# Fields` section to the end of its docstring, after the reference.

- `src/Kernels.jl`:69  (AsymmetricKernel)
- `src/MarginalDensity.jl`:74  (ParametricStart)
- `src/Histograms.jl`:70  (Histogram1D)
- `src/Histograms.jl`:245  (Histogram2D)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-009  Rendered code blocks carry lines of 93 to 114 characters, well past the 86-to-89-character width of the content column, so each block renders with a horizontal scrollbar; the overflow is almost always a trailing # comment.

Rule: `RULE.TEX.NO_HSCROLL`.  
Reported by `rendered`, decided by an agent; carried over from an earlier run, which left it open.

Measured on the built pages: the article is 52rem wide with 1rem padding (800px); a pre block adds 1.5rem padding per side and a 2px border, a docstring section 0.75rem per side; at 0.875em JuliaMono (0.6em advance) that is ~89 characters at page level and ~86 inside a docstring, and pre is overflow-x:auto with white-space:pre. Lines: Expectations 96 and 95; Histograms 95; Kernels 93; MarginalDensity 108, 98 and 104; Architecture section 8 layout block 100, 100, 104, 107, 98, 102, 114 and 97; tutorial 97, 100 and 94. Every display formula on the fourteen pages fits.

**Fix.** Move the trailing comment onto its own line above the statement, or shorten it; in the Architecture layout block break the long per-file comments onto continuation lines and split the test/ line into two.

- `src/Expectations.jl`:207  (expectation)
- `src/Expectations.jl`:275  (expectation_with_standard_error)
- `src/Histograms.jl`:239  (Histogram2D)
- `src/Kernels.jl`:101  (AsymmetricKernel)
- `src/MarginalDensity.jl`:502  (density)
- `src/MarginalDensity.jl`:63  (ParametricStart)
- `src/MarginalDensity.jl`:64  (ParametricStart)
- `docs/ARCHITECTURE.md`:589  (GBGStreamingStatistics — Architecture and Rationale / 8. Proposed repository layout (mirrors PF))
- `docs/ARCHITECTURE.md`:591  (GBGStreamingStatistics — Architecture and Rationale / 8. Proposed repository layout (mirrors PF))
- `docs/ARCHITECTURE.md`:594  (GBGStreamingStatistics — Architecture and Rationale / 8. Proposed repository layout (mirrors PF))
- `docs/ARCHITECTURE.md`:589  (GBGStreamingStatistics — Architecture and Rationale / 8. Proposed repository layout (mirrors PF))
- `docs/ARCHITECTURE.md`:598  (GBGStreamingStatistics — Architecture and Rationale / 8. Proposed repository layout (mirrors PF))
- `docs/ARCHITECTURE.md`:604  (GBGStreamingStatistics — Architecture and Rationale / 8. Proposed repository layout (mirrors PF))
- `docs/ARCHITECTURE.md`:610  (GBGStreamingStatistics — Architecture and Rationale / 8. Proposed repository layout (mirrors PF))
- `docs/ARCHITECTURE.md`:614  (GBGStreamingStatistics — Architecture and Rationale / 8. Proposed repository layout (mirrors PF))
- `docs/src/tutorial.md`:23  (Getting started / The workflow)
- `docs/src/tutorial.md`:32  (Getting started / The workflow)
- `docs/src/tutorial.md`:65  (Getting started / The workflow)

### 2026-10-15-documentation-review-01-010  Code lines of 87 to 92 characters sit at or a few pixels past the column width inside a docstring block (about 86 characters) or the Architecture page (about 89), so they render with a slight scrollbar on the default theme; marginal, and font-fallback dependent.

Rule: `RULE.TEX.NO_HSCROLL`.  
Reported by `rendered`, decided by an agent; carried over from an earlier run, which left it open.

Same measurement as above; the docstring pre's inner width is ~722px, which at 8.4px per character is 86 characters. ConditionalDensity's three signature continuation lines are 87; Histograms 88, 90 and 92; Kernels 90; Architecture 92 and 88. The tutorial's 87-character line fits at page level and is not listed.

**Fix.** Shorten or drop the trailing comment, or let the aligned comment column start earlier; the three signature lines lose their trailing ≡ comment or move it to the line above.

- `src/ConditionalDensity.jl`:86  (conditional_density)
- `src/ConditionalDensity.jl`:244  (conditional_density_standard_error)
- `src/ConditionalDensity.jl`:180  (transition_density)
- `src/Histograms.jl`:66  (Histogram1D)
- `src/Histograms.jl`:237  (Histogram2D)
- `src/Histograms.jl`:241  (Histogram2D)
- `src/Kernels.jl`:97  (AsymmetricKernel)
- `docs/ARCHITECTURE.md`:592  (GBGStreamingStatistics — Architecture and Rationale / 8. Proposed repository layout (mirrors PF))
- `docs/ARCHITECTURE.md`:601  (GBGStreamingStatistics — Architecture and Rationale / 8. Proposed repository layout (mirrors PF))

### 2026-10-15-documentation-review-01-011  Three display rows pack a second equation on the same line behind \qquad (the formula, then a symbol's definition), the form the rule names; each still fits the column today, so this is the rule's letter rather than a visible scroll.

Rule: `RULE.TEX.NO_HSCROLL`.  
Reported by `rendered`, decided by an agent; carried over from an earlier run, which left it open.

The rule reads 'Never pack several equations on one line behind \quad.' On the built Bandwidth page the row sits inside an aligned block whose other rows hold one equation each, so a '\\' row for \mu_2(K) is the environment's own idiom; the Accumulators and ConditionalDensity rows are single-line displays with a trailing where-clause. The CLAUDE.md note that \quad/\qquad are correct forms is about the spacing command, not about packing. The widest of the three measures roughly 60 visible characters, under the 800px column.

**Fix.** Give the second equation its own row (a '\\' in the aligned block; a gathered or aligned pair for the two single-line displays), or leave as is if the owner reads the where-clause as a sanctioned exception and records that.

- `src/Bandwidth.jl`:337  (bandwidth)
- `src/Accumulators.jl`:76  (merge!)
- `src/ConditionalDensity.jl`:256  (conditional_density_standard_error)

### 2026-10-15-documentation-review-01-014  The range "Expectations through Bandwidth" names an ordering of modules that this page does not show, so the reader cannot tell which modules the upper half contains.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api`, decided by an agent; carried over from an earlier run, which left it open.

"the upper half ([`Expectations`](Expectations/) through [`Bandwidth`](Bandwidth/)) holds the estimands" — the lower half is enumerated in full, the upper half only by its two endpoints; on this page alone there is no list between them to read "through" against.

**Fix.** Enumerate the upper-half modules as the lower half is enumerated, or drop "through" and say "the remaining modules, from Expectations to Bandwidth" with the list.

- `docs/src/api.md`:10  (API overview)

### 2026-10-15-documentation-review-01-017  Symbols are used in a docstring that binds them nowhere in itself; they lean on another docstring on the page or are never bound at all.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-Accumulators`, decided by an agent; carried over from an earlier run, which left it open.

`scale!` writes the reweighted measure $\texttt{weight} \cdot \hat\mu$ and the fractional $\texttt{weight} \cdot n$, but $\hat\mu$ and $n$ are bound only in `AbstractAccumulator`'s docstring; `merge!` likewise relies on "the normalized $\hat\mu$" without binding it; in `observe!` the index $p$ of $S_p$ is never said to be the moment order; in `merge!`'s display, $P(x \mid \beta)$ and the spacing $\Delta\beta$ are introduced without a word.

**Fix.** Bind in apposition at first use within each docstring: "the reweighted empirical measure $\texttt{weight} \cdot \hat\mu$ ($\hat\mu$ the empirical measure of the $n$ samples seen)"; "the per-bin power-sums $S_p = \sum (x_i - c_j)^p$ about the bin center $c_j$, $p$ up to the summary's order"; "$P(x \mid \beta)$ the density at fixed parameter, $\Delta\beta$ the spacing of the $\beta_k$".

- `src/Accumulators.jl`:108  (scale!)
- `src/Accumulators.jl`:122  (scale!)
- `src/Accumulators.jl`:59  (merge!)
- `src/Accumulators.jl`:40  (observe!)
- `docs/FORMALIZATION.md`:96  (Formalization / The cone — `scale!`, weighted `merge!`, hyperstatistics (§10))
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-018  A lone equals sign is set as a code span where it means mathematical equality, not Julia (where `=` is assignment).

Rule: `RULE.TEX.CODE_SPAN_JULIA`.  
Reported by `withheld@api-Accumulators`, decided by an agent; carried over from an earlier run, which left it open.

"tested as accumulate-then-merge `=` accumulate-together" backticks a bare `=` standing for "equals".

**Fix.** Write it as prose or mathematics: "tested as accumulate-then-merge equals accumulate-together".

!!! warning "The code may be what is wrong"
    The reader believed the documentation and the code disagree, and that the code is the wrong half. That is a question for the owner, not a documentation fix.

- `src/Accumulators.jl`:89  (merge!)

### 2026-10-15-documentation-review-01-020  The main input $\hat\sigma$ (and $\hat\sigma_X$) is never bound: nowhere does the page say it is the sample standard deviation, and 'recovered' is the package's verb for reading it off the summary, unexplained.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-Bandwidth`, decided by an agent; carried over from an earlier run, which left it open.

The normal-reference formula, the Arguments bullet and the local pilot all lean on $\hat\sigma$; a first-time reader can guess 'a scale estimate' but is told neither what it estimates nor from what. 'Recovered' is used three times as if it named a known operation.

**Fix.** At first use: 'From the sample standard deviation $\hat\sigma$ read back from the summary's stored moments and the sample size $n$'; in the second docstring bind $\hat\sigma_X$ as the standard deviation of the conditioning variable.

- `src/Bandwidth.jl`:321  (bandwidth)
- `src/Bandwidth.jl`:396  (bandwidth)
- `src/Bandwidth.jl`:525  (bandwidth)
- `docs/ARCHITECTURE.md`:414  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:737  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:737  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:790  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:790  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:743  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:802  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/notation.md`:61  (Notation / Smoothing kernels)
- `docs/src/notation.md`:106  (Notation / The moments of the data)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-021  Symbols of the plug-in machinery are used without being bound: the pilot functional $\psi_6$, the derivative order $s$ in $\|f^{(s)}\|^2$, and Chen's 'modified $f_2$ kernels'.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-Bandwidth`, decided by an agent; carried over from an earlier run, which left it open.

$\psi_4$ is bound as $R(f'')$ but $\psi_6$ only as 'its pilot'; $s$ appears once with no statement of what it ranges over; $f_2$ is Chen's own label for a kernel variant and means nothing to a reader of this page alone.

**Fix.** Bind each in apposition: '$\psi_6 = R(f''')$, the pilot functional'; 'the curvature norms $\|f^{(s)}\|^2$, $s = 2, 3, \dots$'; 'the boundary-modified gamma/beta kernels (Chen's $f_2$ variant)'.

- `src/Bandwidth.jl`:346  (bandwidth)
- `src/Bandwidth.jl`:356  (bandwidth)
- `src/Bandwidth.jl`:368  (bandwidth)
- `src/Bandwidth.jl`:347  (bandwidth)
- `docs/ARCHITECTURE.md`:768  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-026  Three docstrings open with a verbless noun phrase (an em-dash apposition in place of a predicate), and one caveat is a verbless fragment; none reads as a sentence aloud.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-ConditionalDensity`, decided by an agent; carried over from an earlier run, which left it open.

The lead sentence of `conditional_density_with_standard_error`, `transition_density_standard_error` and `transition_density_with_standard_error` is 'The X — identical to / the two-time reading of ...' with no finite verb; 'NaN where no data is local to x.' has no verb either.

**Fix.** Give each a predicate: 'Returns the estimate and its standard error ...', 'The standard error of `transition_density` is `conditional_density_standard_error` on ...', 'The estimate is `NaN` where ...'.

- `src/ConditionalDensity.jl`:305  (conditional_density_with_standard_error)
- `src/ConditionalDensity.jl`:345  (transition_density_standard_error)
- `src/ConditionalDensity.jl`:365  (transition_density_with_standard_error)
- `docs/src/api/ConditionalDensity.md`  (conditional_density_standard_error)  — **the text here has moved and was not found again**
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-029  'per-lag joint histogram' and 'the summary's lag' presume a lag the page never shows a `Histogram2D` carrying; the reader cannot tell where tau enters.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-ConditionalDensity`, decided by an agent; carried over from an earlier run, which left it open.

The doctest constructs `Histogram2D(-1.0, 1.0, 2, -1.0, 1.0, 2)` with no lag argument, yet `transition_density` speaks of 'the summary's lag tau' and 'the per-lag joint histogram'. Nothing on the page says the lag is fixed by the caller pairing each sample X_t with X_{t+tau} before `observe_pairs!`.

**Fix.** State it once in the Estimator paragraph: 'a `Histogram2D` fed the pairs (X_t, X_{t+tau}) for one fixed lag tau chosen by the caller — the summary itself stores no lag'.

- `docs/src/api/ConditionalDensity.md`  (transition_density)  — **the text here has moved and was not found again**
- `src/ConditionalDensity.jl`:189  (transition_density)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-030  The sample pairs X_i, Y_i, the sample size n and the semicolon kernel form K_2(c; y, h_2) appear in formulas without being bound.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-ConditionalDensity`, decided by an agent; carried over from an earlier run, which left it open.

The first display sums over i with X_i, Y_i never introduced as the observed pairs; n first appears in the variance formula as a bare denominator and is never said to be the number of pairs; K_2(c^y_k; y, h_2) is a different kernel notation from K_{h_2}(c^y_k - y) and is not explained.

**Fix.** Before the first display: 'for n observed pairs (X_i, Y_i)'; and gloss the asymmetric form: 'K_2(u; y, h_2), the kernel's density in its data argument u with location y and scale h_2'.

- `src/ConditionalDensity.jl`:100  (conditional_density)
- `src/ConditionalDensity.jl`:253  (conditional_density_standard_error)
- `src/ConditionalDensity.jl`:113  (conditional_density)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-032  Compressed prose lets a symbol or an ellipsis stand in for words: '$\approx$' used as an adverb, and a rate compared with a density.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-ConditionalDensity`, decided by an agent; carried over from an earlier run, which left it open.

'but only ≈ normalized in y' reads a math symbol as the word 'approximately'; 'The conditional MISE rate is slower than a marginal density' compares a rate to a density (the marginal density's rate is meant).

**Fix.** 'but only approximately normalized in y'; 'The conditional MISE rate, O(n^{-2/3}), is slower than a marginal density's.'

- `src/ConditionalDensity.jl`:113  (conditional_density)
- `src/ConditionalDensity.jl`:126  (conditional_density)

### 2026-10-15-documentation-review-01-036  The summary's structural vocabulary — '$y$-quadrature', '$X$-column', 'order-2 $Y$-power-sums', '$Y$-bin scale', 'column sums' — is used as known on this page, which never says what a power-sum is or what a column of the joint summary is.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent; carried over from an earlier run, which left it open.

The page preamble ends on 'a $y$-quadrature response' before any sentence explains a quadrature of what; `conditional_expectation` then speaks of 'the stored order-2 $Y$-power-sums' and 'per-$X$-column' as given. A reader who has not read the Histograms page cannot tell that a power-sum is the per-cell sum of $Y^p$ or that an $X$-column is one $x$-bin's row of cells; `Histogram2D` is linked but the nouns are not.

**Fix.** Bind the nouns once where they first appear: 'the per-cell sums $\sum Y$, $\sum Y^2$ the summary stores (its order-2 $Y$-power-sums)', 'each $x$-bin of the joint summary (an $X$-column)', and 'the $Y$-bin width $\Delta_y$'.

- `docs/src/api/ConditionalExpectations.md`:9  (ConditionalExpectations)
- `src/ConditionalExpectations.jl`:90  (conditional_expectation)
- `src/ConditionalExpectations.jl`:99  (conditional_expectation)
- `src/ConditionalExpectations.jl`:109  (conditional_expectation)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-037  The symbol $\Delta_y$ is used in three docstrings and bound nowhere on the page.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent; carried over from an earlier run, which left it open.

'$O(\Delta_y^2)$' appears in `conditional_expectation` and `conditional_variance_standard_error`, and '$h_z \gtrsim |g'| \Delta_y$' in `conditional_pushforward_density`, with no clause saying $\Delta_y$ is the $Y$-bin width of the summary.

**Fix.** At first use in each docstring: '$O(\Delta_y^2)$, $\Delta_y$ the $Y$-bin width'.

- `src/ConditionalExpectations.jl`:100  (conditional_expectation)
- `src/ConditionalExpectations.jl`:445  (conditional_variance_standard_error)
- `src/ConditionalExpectations.jl`:535  (conditional_pushforward_density)
- `docs/ARCHITECTURE.md`:947  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-038  Several symbols in the formulas are not bound at first use: $K_h$, $n$, $f_X$, $\mu_2(K)$ and $R(K)$ in `conditional_mean` (with $R(K)$ bound only later, in `conditional_pushforward_density_standard_error`), and $N$ in `conditional_variance`.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent; carried over from an earlier run, which left it open.

The design sums use $K_h(X_i-x)$ without saying $K_h(u) = K(u/h)/h$; the consistency line uses $\mu_2(K)$, $R(K)$, $n$ and $f_X(x)$ unglossed, after the same density was written $f(x)$ a sentence earlier; `conditional_variance`'s residual total ends in '$\hat m^2 N$' with $N$ never named (the per-column count is $N_j$ only in another docstring).

**Fix.** Bind in apposition at first use: '$K_h(u) = K(u/h)/h$ the scaled kernel', '$n$ the sample count', '$f_X$ the marginal density of $X$', '$\mu_2(K) = \int u^2 K$, $R(K) = \int K^2$', and '$N$ the column count'; use one symbol for the design density.

- `src/ConditionalExpectations.jl`:281  (conditional_mean)
- `src/ConditionalExpectations.jl`:264  (conditional_mean)
- `src/ConditionalExpectations.jl`:375  (conditional_variance)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-039  The package's jargon 'reader' / 'readout' (a function that reads an estimate off the summary) is used without definition.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent; carried over from an earlier run, which left it open.

`conditional_expectation_standard_error` calls its estimator 'The leverage × noise-scale readout' and speaks of 'two `conditional_expectation` readouts' and 'the $g^2$ reader'; nothing on the page says a readout is an evaluation of a function of the stored summary.

**Fix.** Replace with plain words ('the product of the leverage and the noise scale', 'two evaluations of `conditional_expectation`, the second with $g^2$ in place of $g$') or define 'readout' once.

- `src/ConditionalExpectations.jl`:172  (conditional_expectation_standard_error)
- `src/ConditionalExpectations.jl`:175  (conditional_expectation_standard_error)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-040  The type `AsymmetricKernel` appears in every signature's `degree` default and is never named, linked or explained on the page; 'the Chen kernels' is used as a known name with no citation or definition.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent; carried over from an earlier run, which left it open.

A reader meets 'kernel isa AsymmetricKernel ? 0 : 1' eight times; the prose only ever says 'an asymmetric `kernel`' and links `gamma_kernel` / `beta_kernel`. `conditional_mean` then says 'the Chen kernels are a density boundary tool' — no Chen entry is under its References and the page never ties the name to the gamma/beta kernels.

**Fix.** In the first `kernel` bullet, link the type ('an [`AsymmetricKernel`](@ref …) such as `gamma_kernel` / `beta_kernel`') and write 'the gamma/beta boundary kernels of Chen (1999, 2000)' with the reference added, or drop the eponym.

- `src/ConditionalExpectations.jl`:292  (conditional_mean)
- `src/ConditionalExpectations.jl`:81  (conditional_expectation)
- `src/ConditionalExpectations.jl`:211  (conditional_expectation_with_standard_error)
- `src/ConditionalExpectations.jl`:211  (conditional_expectation_with_standard_error)
- `src/ConditionalExpectations.jl`:255  (conditional_mean)
- `src/ConditionalExpectations.jl`:255  (conditional_mean)
- `src/ConditionalExpectations.jl`:162  (conditional_expectation_standard_error)
- `src/ConditionalExpectations.jl`:431  (conditional_variance_standard_error)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-041  The labelled 'Estimand.' / 'Estimator.' openers are verbless noun phrases in six places — a label followed by a description with no finite verb.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent; carried over from an earlier run, which left it open.

'Estimator. A local-polynomial regression in $x$ of a per-$X$-column $y$-quadrature of $g$.' and 'Estimator. The two-term $\sqrt{\dots}$: the local-regression leverage … times … minus …' have no verb; read aloud they are fragments. The openers that carry an equation ('$\mathbb{E}[\dots] = \int\dots$') or a 'since … is' clause do read as sentences and are not counted here.

**Fix.** Give each a verb: 'Estimator. The estimate is a local-polynomial regression in $x$ of …', 'Estimator. The band is the two-term …', 'Estimand. The estimand is the pointwise standard error …'.

- `src/ConditionalExpectations.jl`:90  (conditional_expectation)
- `src/ConditionalExpectations.jl`:172  (conditional_expectation_standard_error)
- `src/ConditionalExpectations.jl`:377  (conditional_variance)
- `src/ConditionalExpectations.jl`:518  (conditional_pushforward_density)
- `src/ConditionalExpectations.jl`:609  (conditional_pushforward_density_standard_error)
- `docs/src/api/ConditionalExpectations.md`  (conditional_pushforward_density_standard_error)  — **the text here has moved and was not found again**
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-042  Verbless fragments stand as sentences in the body prose: the vector-case return notes, the 'Two smoothing directions:' list, 'Symmetric `kernelz` only.', '`NaN` where …', the preamble's closing 'All by …', and a parenthesis of two abutting code spans.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent; carried over from an earlier run, which left it open.

'For a vector `x`, a `NamedTuple` of vectors.' and 'Symmetric `kernelz` only.' have no verb; '(`method` `:residual` default, or `:plugin`)' is not parseable as a phrase. The sibling `conditional_expectation_with_standard_error` already writes the same note with a verb, 'returns a `NamedTuple` of vectors'.

**Fix.** 'For a vector `x` it returns a `NamedTuple` of vectors.'; 'It smooths in two directions: …'; 'Only a symmetric `kernelz` is supported.'; 'It is `NaN` where …'; 'All are computed by …'; '(`method = :residual`, the default, or `:plugin`)'.

- `src/ConditionalExpectations.jl`:484  (conditional_variance_with_standard_error)
- `src/ConditionalExpectations.jl`:656  (conditional_pushforward_density_with_standard_error)
- `src/ConditionalExpectations.jl`:526  (conditional_pushforward_density)
- `src/ConditionalExpectations.jl`:620  (conditional_pushforward_density_standard_error)
- `docs/src/api/ConditionalExpectations.md`  (conditional_expectation_standard_error)  — **the text here has moved and was not found again**
- `src/ConditionalExpectations.jl`:481  (conditional_variance_with_standard_error)
- `src/Expectations.jl`:375  (pushforward_density_with_standard_error)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-043  A code span presents a schematic call as literal Julia: `K_hz` is not a binding the page or the package offers, so the span is pseudo-code standing in for the mathematics $K_{h_z}$ written just above it.

Rule: `RULE.TEX.CODE_SPAN_JULIA`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent; carried over from an earlier run, which left it open.

'i.e. `conditional_expectation(summary, y -> K_hz(g(y) - z), x)`' — a reader who types it gets an undefined `K_hz`; the span mixes a real signature with a mathematical symbol transliterated into an identifier.

**Fix.** Either make it runnable, e.g. `conditional_expectation(summary, y -> kernelz((g(y) - z) / hz) / hz, x; kernel = kernelx, degree, h = hx)`, or keep the statement in mathematics and drop the code span.

!!! warning "The code may be what is wrong"
    The reader believed the documentation and the code disagree, and that the code is the wrong half. That is a question for the owner, not a documentation fix.

- `src/ConditionalExpectations.jl`:526  (conditional_pushforward_density)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-045  The two section-opening paragraphs of the page are noun phrases with no finite verb.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-Expectations`, decided by an agent; carried over from an earlier run, which left it open.

The Expectations intro runs 'The power-sum quadrature readouts — a smooth functional of X read straight off the stored power-sums, with no kernel and no query point.' ('read' is a participle; nothing predicates). The Recovered moments intro likewise: 'The global marginal moments — the g = x^p case — recovered post hoc from the power-sums by the shift formula (exact), extending Statistics.'

**Fix.** Give each a verb: 'The power-sum quadrature readouts read a smooth functional of X straight off the stored power-sums, with no kernel and no query point.' / 'The global marginal moments — the g = x^p case — are recovered post hoc from the power-sums by the shift formula, exactly, extending Statistics.'

- `docs/src/api/Expectations.md`:3  (Expectations)
- `docs/src/api/Expectations.md`:23  (Expectations / Recovered moments)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-046  The labelled 'Estimand.' / 'Estimator.' paragraphs and the paired readers' identity sentences are verbless: a label followed by a noun phrase or a bare formula with appositions.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-Expectations`, decided by an agent; carried over from an earlier run, which left it open.

'Estimator. $\sqrt{(\widehat{\mathbb{E}}[g^2] - \widehat{\mathbb{E}}[g]^2)/n}$ — both moments the same power-sum quadrature (...), $n$ the sample count.' has no verb at all; 'Estimand. The population mean $\mathbb{E}[X]$.' is a noun phrase; 'Estimator. The sample mean ... — recovered exactly ... host-side ...' predicates nothing; 'The estimate and its standard error from one power-sum quadrature pass — the `NamedTuple` ...' opens expectation_with_standard_error with no verb, and pushforward_density_with_standard_error repeats the form. Read aloud, none is a sentence.

**Fix.** Supply the verb once per paragraph: 'The estimand is the population mean ...', 'The estimator is the square root ... of both moments by the same power-sum quadrature, with $n$ the sample count.', 'Returns the estimate and its standard error from one pass as the `NamedTuple` ...'.

- `src/Expectations.jl`:241  (expectation_standard_error)
- `src/Expectations.jl`:98  (mean)
- `src/Expectations.jl`:100  (mean)
- `src/Expectations.jl`:262  (expectation_with_standard_error)
- `src/Expectations.jl`:372  (pushforward_density_with_standard_error)
- `src/Expectations.jl`:113  (var)
- `src/Expectations.jl`:152  (central_moment)
- `docs/ARCHITECTURE.md`:452  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.B Power-sum quadrature readouts)
- `docs/ARCHITECTURE.md`:452  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.B Power-sum quadrature readouts)
- `src/Expectations.jl`:242  (expectation_standard_error)
- `src/Expectations.jl`:111  (var)
- `src/Expectations.jl`:150  (central_moment)
- `docs/ARCHITECTURE.md`:453  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.B Power-sum quadrature readouts)
- `docs/src/notation.md`:104  (Notation / The moments of the data)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-049  'host-side' is GPU jargon a statistician reading this page has no way to parse, and nothing on the page mentions a device.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Expectations`, decided by an agent; carried over from an earlier run, which left it open.

mean says 'recovered exactly (not a smoothed estimate) host-side from the per-bin counts and first power-sums'; neither 'host' nor a GPU or device array appears anywhere else on the page, so the adverb is an unexplained term from another layer.

**Fix.** Drop 'host-side' (it is an implementation detail), or say it plainly: 'on the CPU, from the per-bin counts ...'.

- `src/Expectations.jl`:100  (mean)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-050  Three symbols carry the main formulas and are bound nowhere on the page: `order` (set as $\texttt{order}$), $n$ inside expectation, and $K_h$.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-Expectations`, decided by an agent; carried over from an earlier run, which left it open.

expectation's display sums 'to $p = \texttt{order}$' and the text says 'exact for a polynomial $g$ of degree $\le \texttt{order}$', 'a degree-$\texttt{order}$ polynomial interpolant', central_moment 'for $2 \le p \le \texttt{order}$' — no sentence says `order` is the highest power-sum the summary stores, nor that it is a property of `Histogram1D`. The same display divides by $n$, bound only two docstrings later ('$n$ the sample count' in expectation_standard_error). pushforward_density writes $K_h(g(X_i) - y)$ and then 'the surviving $1/h$', but never states $K_h(u) = K(u/h)/h$; a reader who takes $K_h$ as $K$ with a parameter cannot see where the $1/h$ comes from.

**Fix.** In expectation: 'with $n = \sum_j n_j$ the sample count and `order` the highest power stored by the summary'; in pushforward_density: 'with $K_h(u) = K(u/h)/h$ the kernel scaled to the bandwidth $h$'.

- `src/Expectations.jl`:191  (expectation)
- `src/Expectations.jl`:152  (central_moment)
- `src/Expectations.jl`:308  (pushforward_density)
- `docs/ARCHITECTURE.md`:330  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.1 Marginal density ``f_X``)
- `docs/ARCHITECTURE.md`:576  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- `docs/ARCHITECTURE.md`:729  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:729  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/api/Kernels.md`:4  (Kernels)
- `docs/src/index.md`:85  (GBGStreamingStatistics.jl / Main components)
- `src/Bandwidth.jl`:538  (bandwidth)
- `src/ConditionalDensity.jl`:95  (conditional_density)
- `src/ConditionalDensity.jl`:95  (conditional_density)
- `src/ConditionalDensity.jl`:95  (conditional_density)
- `src/ConditionalExpectations.jl`:112  (conditional_expectation)
- `src/ConditionalExpectations.jl`:303  (conditional_mean)
- `src/ConditionalExpectations.jl`:387  (conditional_variance)
- `src/ConditionalExpectations.jl`:528  (conditional_pushforward_density)
- `src/ConditionalExpectations.jl`:547  (conditional_pushforward_density)
- `src/Expectations.jl`:319  (pushforward_density)
- `src/Kernels.jl`:19  (gaussian_kernel)
- `src/Kernels.jl`:19  (gaussian_kernel)
- `src/Kernels.jl`:30  (epanechnikov_kernel)
- `src/MarginalDensity.jl`:428  (density)
- `src/MarginalDensity.jl`:441  (density)
- `src/MarginalDensity.jl`:704  (density_standard_error)
- `docs/ARCHITECTURE.md`:730  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/notation.md`:46  (Notation / Smoothing kernels)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-051  Two docstrings admit an input for which they never say what comes back: pushforward_density_standard_error with a vector `y`, and central_moment with $p$ outside $2 \le p \le$ `order`.

Rule: `RULE.DOC.RETURNS`.  
Reported by `withheld@api-Expectations`, decided by an agent; carried over from an earlier run, which left it open.

pushforward_density_standard_error delegates its arguments to pushforward_density, whose `y` is 'a scalar, or a vector for the curve', yet its own text speaks only of the pointwise $\mathrm{SE}(y)$; its paired reader does say 'For a vector `y`, a `NamedTuple` of vectors', so the SE reader alone is silent. central_moment states its estimator 'for $2 \le p \le \texttt{order}$' and never says what $p = 1$, $p = 0$ or $p >$ `order` return (an error, `NaN`, or a value).

**Fix.** Add 'For a vector `y`, a vector of standard errors.' to pushforward_density_standard_error, and to central_moment a sentence on the out-of-range $p$ ('throws an `ArgumentError`' or whatever the code does).

- `src/Expectations.jl`:342  (pushforward_density_standard_error)
- `src/Expectations.jl`:152  (central_moment)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-055  Symbols bound nowhere on the page: the map T and the measure \hat\mu_n in 'the canonical summary T(\hat\mu_n)', and M_4 / M_6 (presumably the central moments of order 4 and 6, but the page never says so).

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-Histograms`, decided by an agent; carried over from an earlier run, which left it open.

A cold reader meets $T(\hat\mu_n)$ in the identity sentence of Histogram1D with no clause saying \hat\mu_n is the empirical measure of the n samples and T the summarizing map; $M_4$ appears three times (prose and the `order` bullet) with M never bound.

**Fix.** Bind in apposition at first use: 'the summary T(\hat\mu_n) of the empirical measure \hat\mu_n of the n samples seen so far' and 'the central moment M_4 of order 4'; or drop T(\hat\mu_n) from the identity sentence.

- `src/Histograms.jl`:21  (Histogram1D)
- `src/Histograms.jl`:32  (Histogram1D)
- `src/Histograms.jl`:49  (Histogram1D)
- `docs/ARCHITECTURE.md`:72  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The spine: empirical measure → summary → post-hoc-smoothed estimator)
- `docs/ARCHITECTURE.md`:87  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The spine: empirical measure → summary → post-hoc-smoothed estimator)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-059  A sentence that does not parse: 'For a two-time quantity X(t) and the response Y(t+τ) are paired per source' lacks the comma or subject that would make 'X(t) and the response' the subject of 'are paired'.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-Histograms`, decided by an agent; carried over from an earlier run, which left it open.

Read aloud, 'For a two-time quantity X(t) and the response Y(t+τ) are paired' first attaches X(t) to 'quantity' and then has no subject for 'are paired'.

**Fix.** 'For a two-time quantity, the value X(t) and the later response Y(t+τ) of one source are paired; the lag τ is an explicit label on the summary.'

- `src/Histograms.jl`:214  (Histogram2D)
- `docs/ARCHITECTURE.md`:282  (GBGStreamingStatistics — Architecture and Rationale / 3. The accumulated summary — the universal lower half)
- `docs/src/tutorial.md`:54  (Getting started / The workflow)
- `docs/ARCHITECTURE.md`:283  (GBGStreamingStatistics — Architecture and Rationale / 3. The accumulated summary — the universal lower half)
- `docs/ARCHITECTURE.md`:395  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.3 Conditional and transition density)
- `docs/src/notation.md`:17  (Notation / The samples)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-062  The gamma_kernel, beta_kernel and epanechnikov_kernel docstrings use the symbols b, n, n_j, c_j, \hat\sigma and h without binding them; only the AsymmetricKernel docstring, further down the page, names b the bandwidth and n_j / c_j the bin counts and centers.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-Kernels`, decided by an agent; carried over from an earlier run, which left it open.

In gamma_kernel, $b$ first appears inside $\mathrm{Gamma}(\rho(x), b)$ with no gloss; the estimate $\hat f(x) = \tfrac{1}{n}\sum_j n_j\,\mathrm{Gamma}(c_j; \rho(x), b)$ introduces $n$, $n_j$ and $c_j$ cold; the bandwidth default $\hat\sigma n^{-2/5}$ leaves $\hat\sigma$ unnamed. In epanechnikov_kernel, 'the bins within $h$ of the query point' is the first and only $h$ on the page, unbound. beta_kernel inherits the unbound $b$ in $(x/b, (1-x)/b)$. A reader who opens gamma_kernel directly (the docstring is what `?gamma_kernel` shows) has none of AsymmetricKernel's bindings.

**Fix.** Bind each symbol in apposition at its first use in each docstring: 'the bandwidth $b$', 'the sample size $n$', 'the bin counts $n_j$ at the bin centers $c_j$', 'the sample standard deviation $\hat\sigma$', 'the bandwidth $h$'.

- `src/Kernels.jl`:139  (gamma_kernel)
- `src/Kernels.jl`:145  (gamma_kernel)
- `src/Kernels.jl`:147  (gamma_kernel)
- `src/Kernels.jl`:32  (epanechnikov_kernel)
- `src/Kernels.jl`:175  (beta_kernel)
- `docs/ARCHITECTURE.md`:790  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `src/Kernels.jl`:160  (Kernels)
- `src/Kernels.jl`:176  (beta_kernel)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-063  The package term `summary` (the data object every estimator takes) is used in two code spans and in the runnable example without being defined or constructed anywhere on the page.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Kernels`, decided by an agent; carried over from an earlier run, which left it open.

The page never says what a 'summary' is — the binned histogram the estimators read — yet `bandwidth(summary; rule = :plugin, kernel = gamma_kernel)` and the example's `density(summary, 1.0; kernel = mygamma)` both take it as their first argument. In the example `summary` is a free name: the block cannot be run as written, and a cold reader does not know what to put there. The module docstring, which mentions bins nowhere, does not supply the binned-measure premise either.

**Fix.** In the module docstring, bind the term once ('every estimator reads a `summary`, the binned histogram of the stream, …' with an @ref to its type), and either construct a summary at the top of the example or state in a comment what `summary` is.

- `src/Kernels.jl`:150  (gamma_kernel)
- `src/Kernels.jl`:103  (AsymmetricKernel)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-067  beta_kernel describes its edge correction as 'a $\sqrt{\cdot}$ correction', a placeholder standing in for a formula the page never states.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-Kernels`, decided by an agent; carried over from an earlier run, which left it open.

gamma_kernel spells its edge shape out ('$(x/b)^2/4 + 1$ within $2b$ of the edge'), so the reader expects the same of the beta kernel, but gets '$\sqrt{ \cdot }$' with nothing inside: no symbol is bound, no formula given, and the reference to 'Chen's $f_2$' is the only way to learn what the correction is.

**Fix.** State the edge shape explicitly — the shape parameter Chen's $f_2$ substitutes near each edge — or say in prose which argument is replaced and by what.

- `src/Kernels.jl`:175  (beta_kernel)
- pass 1: verifier `NOT_FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-068  Identity sentences across the page are verbless noun phrases (the page lead, the section lead, the Estimator. blocks, the paired reader, every parametric start and the type) — none reads as a sentence aloud.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-MarginalDensity`, decided by an agent; carried over from an earlier run, which left it open.

Each opens the docstring as a noun phrase with modifiers and no finite verb: 'A log-normal ParametricStart on (0,∞), parameters by method of moments (…)'; 'A parametric start for the Hjort–Glad semiparametric density (Hjort & Glad 1995).'; 'Estimator. The smoothed estimate … from the accumulated summary — a binned local smoother …'; 'The density estimate and its pointwise density_standard_error, from a single option specification — so …' (only the subordinate clause carries a verb).

**Fix.** Give each a verb: 'normal_start is the Gaussian ParametricStart on …'; 'ParametricStart is a parametric start for …'; 'The estimator is the smoothed estimate …'; 'density_with_standard_error returns the density estimate and its pointwise standard error from a single option specification, so …'.

- `src/MarginalDensity.jl`:104  (lognormal_start)
- `src/MarginalDensity.jl`:36  (ParametricStart)
- `src/MarginalDensity.jl`:419  (density)
- `src/MarginalDensity.jl`:779  (density_with_standard_error)
- `docs/src/api/MarginalDensity.md`:4  (MarginalDensity)
- `docs/src/api/MarginalDensity.md`:18  (MarginalDensity / Parametric starts)
- `docs/ARCHITECTURE.md`:330  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.1 Marginal density ``f_X``)
- `src/Bandwidth.jl`:317  (bandwidth)
- `src/MarginalDensity.jl`:696  (density_standard_error)
- `src/MarginalDensity.jl`:696  (density_standard_error)
- `src/MarginalDensity.jl`:699  (density_standard_error)
- `src/MarginalDensity.jl`:702  (density_standard_error)
- `docs/ARCHITECTURE.md`:410  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `src/MarginalDensity.jl`:706  (density_standard_error)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-070  Symbols in density's formulas are bound nowhere on the page: $h_j$ and $n$ in the degree-0 formula, $\hat\beta_0$ in the local-log readout, $a$, $b$ and $K_h$ in the support renormalization, and $m$, $s$ in the log-normal moment fit.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-MarginalDensity`, decided by an agent; carried over from an earlier run, which left it open.

The degree-0 KDE is written with a per-bin $h_j$ although `h` is one number and `adaptive` is a later bullet (the phrase 'the sample-point $h_j$ per bin' arrives two bullets on and does not say what $h_j$ is under the default); $\hat\beta_0$ is never named as the fitted intercept; the integral's limits $a$, $b$ are never tied to `(lower, upper)` and $K_h$ is never written out; $m$ and $s$ are never named as the log-normal's log-scale mean and standard deviation.

**Fix.** Write the default formula with $h$ and say 'with $h_j = h$ unless `adaptive`'; bind 'the fitted intercept $\hat\beta_0$'; write $\int_{\text{lower}}^{\text{upper}} K((x-u)/h)\,\mathrm{d}u/h$ or bind $K_h(u) = K(u/h)/h$; name $m$, $s$ as the log-mean and log-standard-deviation; bind $n$ as the total count on first use.

- `src/MarginalDensity.jl`:440  (density)
- `src/MarginalDensity.jl`:440  (density)
- `src/MarginalDensity.jl`:466  (density)
- `src/MarginalDensity.jl`:104  (lognormal_start)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-071  The symbol $b$ is bound twice with different meanings in one docstring: the asymmetric kernel's scale, then the upper limit of the support integral.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-MarginalDensity`, decided by an agent; carried over from an earlier run, which left it open.

'it takes its own scale $b$ in place of `h`' (kernel bullet, repeated in the h bullet) and two bullets later '$\int_a^b K_h(x-u)\,\mathrm{d}u$', where $b$ is the upper support bound.

**Fix.** Use the tuple's own names for the limits ($\int_{\text{lower}}^{\text{upper}}$) or rename the asymmetric scale.

- `src/MarginalDensity.jl`:439  (density)
- `src/MarginalDensity.jl`:466  (density)

### 2026-10-15-documentation-review-01-072  The estimand calls $X$ 'the conditioning variable' on a page about a marginal density of a one-dimensional summary, where nothing is conditioned on; the phrase assumes the package's joint/conditional context.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-MarginalDensity`, decided by an agent; carried over from an earlier run, which left it open.

'The marginal (Lebesgue) density $f_X = \mathrm{d}\mu_X/\mathrm{d}\lambda$ of the conditioning variable $X$' — a cold reader asks: conditioning what? No conditional estimator appears on this page.

**Fix.** 'of the variate $X$' (or 'of the observed variable $X$').

- `src/MarginalDensity.jl`:416  (density)

### 2026-10-15-documentation-review-01-076  The Kernels decision-log entry opens with a sentence that has no main verb and an unclosed parenthesis: the subject 'The varying-shape boundary kernels (... Chen 1999,' runs straight into 'moved from MarginalDensity to Kernels, joining ...' with the citation cut off and no finite verb.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@architecture`, decided by an agent; carried over from an earlier run, which left it open.

Rendered text reads: 'The varying-shape boundary kernels (`gamma_kernel`, `beta_kernel`, the `AsymmetricKernel` type; Chen 1999, moved from `MarginalDensity` to `Kernels`, joining the symmetric `gaussian_kernel` / `epanechnikov_kernel`.' The parenthesis opened before `gamma_kernel` never closes, the second Chen year is missing, and the sentence cannot be read aloud as a statement.

**Fix.** Close the parenthesis and give the sentence its verb: '... the `AsymmetricKernel` type; Chen 1999, 2000) were moved from `MarginalDensity` to `Kernels`, joining ...'.

- `docs/ARCHITECTURE.md`:880  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-077  The bandwidth section's 'Estimand.' and 'Estimator.' labels are followed by verbless noun phrases, each starting lowercase after a full stop, so the two 'sentences' never say what the estimand or estimator is.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@architecture`, decided by an agent; carried over from an earlier run, which left it open.

'Estimand. the AMISE-optimal bandwidth (AMISE = asymptotic mean integrated squared error) — the $h$ asymptotically minimising ...' and 'Estimator. for a marginal, the kernel- and degree-aware normal-reference rule $h = ...$, plus the data-driven ...' — long appositive lists with no finite verb; the other sections' labels at least lead into an equation.

**Fix.** Give each a verb: 'Estimand. The estimand is the AMISE-optimal bandwidth ...'; 'Estimator. For a marginal the package uses the kernel- and degree-aware normal-reference rule ...'.

- `docs/ARCHITECTURE.md`:410  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:412  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))

### 2026-10-15-documentation-review-01-078  Several mathematical symbols are used before the page binds them, or are never bound: $S_p$ (§2 lens, bound only in §3), $\ell(x)$ in the conditional standard error (never bound as the leverage), $\delta(K)$ (§4.5, bound only in the §11 decision log), $\beta_0$ (§4.2, bound only in §11), $L$ and $N$ in the size estimates, and $D^{(n)}(x)$ (§6, bound only in §11).

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@architecture`, decided by an agent; carried over from an earlier run, which left it open.

§2: 'the Taylor contraction of the stored $S_p$' precedes §3's '$S_p = \sum_{i\in\text{bin}}(X_i - c_j)^p$'. §4.2: '$\sqrt{\ell(x) \mathrm{Var}(g(Y)\mid x)}$' with $\ell(x)$ unstated (and $\ell_j$ already means the bin location in §4.1). §4.5: 'the eval kernel a clean $\delta(K)$ rescale' before §11's '$\delta(K) = [R(K)/\mu_2(K)^2]^{1/5}$'. §4.2: 'the $\beta_0$ of one `LocalPolynomial.local_fit`' before §11's '(\beta_0, \beta_1, \dots)'. §3: '($L \times$ bins $\times$ moments, a few MB)' and §5: 'over $N^3\times$steps' with $L$, $N$ unbound. §6: '$D^{(n)}(x)$' before §11's '$D^{(1)} = M_1/\tau$'.

**Fix.** Bind each at first use in apposition: 'the stored per-bin power-sums $S_p$ (§3)'; 'the local-regression leverage $\ell(x) = \sum_j N_j w_j^2$'; 'the kernel's canonical factor $\delta(K) = [R(K)/\mu_2(K)^2]^{1/5}$'; 'the intercept $\beta_0$ of the local-polynomial fit'; '$L$ lags', '$N^3$ grid points per step'; 'the Kramers–Moyal coefficients $D^{(n)}(x)$'.

- `docs/ARCHITECTURE.md`:129  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The algebra of measures — a structural lens)
- `docs/ARCHITECTURE.md`:388  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.2 Conditional expectation ``m(x)`` and variance ``\sigma^2(x)``)
- `docs/ARCHITECTURE.md`:419  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:377  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.2 Conditional expectation ``m(x)`` and variance ``\sigma^2(x)``)
- `docs/ARCHITECTURE.md`:301  (GBGStreamingStatistics — Architecture and Rationale / 3. The accumulated summary — the universal lower half)
- `docs/ARCHITECTURE.md`:498  (GBGStreamingStatistics — Architecture and Rationale / 5. Backend-agnosticism and the hot-loop contract)
- `docs/ARCHITECTURE.md`:536  (GBGStreamingStatistics — Architecture and Rationale / 6. Invariants — deliberate absences are statements)
- `docs/ARCHITECTURE.md`:389  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.2 Conditional expectation ``m(x)`` and variance ``\sigma^2(x)``)
- `docs/ARCHITECTURE.md`:740  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:740  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:740  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:755  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `src/Bandwidth.jl`:335  (bandwidth)
- `src/Bandwidth.jl`:335  (bandwidth)
- `src/Bandwidth.jl`:350  (bandwidth)
- `src/Bandwidth.jl`:360  (bandwidth)
- `src/Bandwidth.jl`:363  (bandwidth)
- `src/Bandwidth.jl`:402  (bandwidth)
- `docs/ARCHITECTURE.md`:301  (GBGStreamingStatistics — Architecture and Rationale / 3. The accumulated summary — the universal lower half)
- `docs/ARCHITECTURE.md`:499  (GBGStreamingStatistics — Architecture and Rationale / 5. Backend-agnosticism and the hot-loop contract)
- `docs/ARCHITECTURE.md`:157  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / Estimand families (the operation taxonomy))
- `docs/ARCHITECTURE.md`:537  (GBGStreamingStatistics — Architecture and Rationale / 6. Invariants — deliberate absences are statements)
- `docs/ARCHITECTURE.md`:749  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:770  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/notation.md`:60  (Notation / Smoothing kernels)
- `docs/src/notation.md`:61  (Notation / Smoothing kernels)
- `src/Bandwidth.jl`:336  (bandwidth)
- `docs/src/notation.md`:81  (Notation / Estimands)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-080  'Channel' ('response channel', 'heavy-tailed channels', 'kilobytes per channel', 'a channel's modal-bin contention') is a package term the page never defines; a first-time reader cannot tell whether a channel is a response variable, a summary, or a consumer's data stream.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@architecture`, decided by an agent; carried over from an earlier run, which left it open.

§2: 'Accumulation knows nothing of bandwidth, kernel, or which response channel will be asked for' and 'every $h$, every channel, every secondary smoothing'; §3: 'raise $P$ for the heavy-tailed channels', 'kilobytes per channel'; §5: 'should a channel's modal-bin contention dominate' — no clause binds the noun anywhere.

**Fix.** Define once at the §2 first use: 'which response channel — which response variable $Y$, each with its own summary — will be asked for', and use that meaning consistently.

- `docs/ARCHITECTURE.md`:92  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The spine: empirical measure → summary → post-hoc-smoothed estimator)
- `docs/ARCHITECTURE.md`:222  (GBGStreamingStatistics — Architecture and Rationale / 3. The accumulated summary — the universal lower half)

### 2026-10-15-documentation-review-01-081  Abbreviations are used before or without their expansion: 'OU' in §3 (Ornstein–Uhlenbeck is spelled out only in §7), 'ULI' (never expanded), 'SFI' (never bound to Stochastic Force Inference), 'NW' (§11), 'MD' (§8), and 'the Pawula $Q$' in §3 ($Q$ is defined only in §11).

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@architecture`, decided by an agent; carried over from an earlier run, which left it open.

§3: 'for the OU oracle both marginals are the same stationary Gaussian'; §2: 'à la Stochastic Force Inference / ULI' and §4.6 'the underdamped ULI'; §8/§10: 'SFI / ULI road', 'include SFI/ULI'; §11: 'count sum, NW ratio, jet'; §8: 'MD's van Hove function'; §3: 'the conditional increment moments of the Pawula $Q$' precede §11's '$Q = \mu_4/(3\mu_2^2)$'.

**Fix.** Expand each at its first use: 'the Ornstein–Uhlenbeck (OU) oracle', 'Stochastic Force Inference (SFI) / Underdamped Langevin Inference (ULI)', 'the Nadaraya–Watson (NW) ratio', 'molecular dynamics (MD)', and state the Pawula quotient $Q = \mu_4/(3\mu_2^2)$ where it first appears in §3.

- `docs/ARCHITECTURE.md`:282  (GBGStreamingStatistics — Architecture and Rationale / 3. The accumulated summary — the universal lower half)
- `docs/ARCHITECTURE.md`:162  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / Estimand families (the operation taxonomy))
- `docs/ARCHITECTURE.md`:620  (GBGStreamingStatistics — Architecture and Rationale / 8. Proposed repository layout (mirrors PF))
- `docs/ARCHITECTURE.md`:897  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:222  (GBGStreamingStatistics — Architecture and Rationale / 3. The accumulated summary — the universal lower half)

### 2026-10-15-documentation-review-01-083  'Stages 1–3' is a project-internal label the page never explains; a reader cannot tell what the stages were or where they are described.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@architecture`, decided by an agent; carried over from an earlier run, which left it open.

§11: 'A framing decision, landed in the code over Stages 1–3.' and 'Landed (Stages 1–3): the shared `LocalPolynomial.local_fit` ...' — no sentence on the page says what a Stage is.

**Fix.** Either drop the label ('landed in the code in three steps: ...') or name the three steps in the entry.

- `docs/ARCHITECTURE.md`:843  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-085  Lists that do not read aloud as sentences: the Status sentence places its 'and' before the penultimate item so the final subject dangles; the MarginalDensity entry has an 'and' mid-list and uses 'the local-log-polynomial' as a noun with no head; the ConditionalExpectations entry stacks two parentheticals so the method clause attaches to the wrong item.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@home`, decided by an agent; carried over from an earlier run, which left it open.

RULE.DOC.GRAMMAR: compact prose that still reads as sentences aloud. 'and the conditional mean and variance (...), the conditional and transition densities (...) are implemented' — the conjunction precedes the wrong item; 'the local-log-polynomial (`degree`), and the Abramson variable bandwidth (`adaptive`, for tails), the asymmetric' — a second 'and' is needed and 'local-log-polynomial' needs 'estimator'; '`conditional_pushforward_density` ($f_{g(Y)\mid X}$) (local-polynomial regression in $x$ of a $y$-quadrature response) — each with its standard error' — the second parenthetical describes all four estimators but is appended to the last.

**Fix.** Move the 'and' before the last item of the Status list ('..., the conditional mean and variance (...), and the conditional and transition densities (...) are implemented'); write 'the local-log-polynomial estimator (`degree`)' and move 'and' to the last item; move the method clause to a sentence of its own ('All are local-polynomial regressions in $x$ of a $y$-quadrature response.').

- `docs/src/index.md`:35  (GBGStreamingStatistics.jl)
- `docs/src/index.md`:96  (GBGStreamingStatistics.jl / Main components)
- `docs/src/index.md`:102  (GBGStreamingStatistics.jl / Main components)
- `docs/ARCHITECTURE.md`:941  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/index.md`:106  (GBGStreamingStatistics.jl / Main components)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-086  The package's terms 'hyperstatistics', 'forgetting' and 'Float-count cone variant' are used without being defined on the page.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@home`, decided by an agent; carried over from an earlier run, which left it open.

RULE.DOC.TERMS_DEFINED: a term of art is defined in a clause, or linked to its definition, before its first use on the same page. 'Hyperstatistics' is not the field's word (the literature says 'superstatistics'); 'forgetting' (exponential down-weighting of old blocks) and 'Float-count cone variant' (a summary whose counts are floating-point so they can be scaled) are the package's own names, and the page nowhere says what they are.

**Fix.** Gloss each at first use: '... or weight and superpose for hyperstatistics (a post-hoc mixture over summaries with different parameters, the superstatistics of Beck & Cohen)'; '... with a variant whose counts are floating-point (the cone variant), so summaries can be re-weighted post hoc for hyperstatistics or exponentially forgotten'.

- `docs/src/index.md`:11  (GBGStreamingStatistics.jl)
- `docs/src/index.md`:32  (GBGStreamingStatistics.jl)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-090  Symbols appear unbound: the lag $\tau$ in $p_\tau(y\mid x)$, the response $Y$ and conditioning value $x$ in $\mathbb{E}[g(Y)\mid x]$ (and the pair $(X, Y)$ generally), and the kernel $K$ in $\int K = 1$.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@home`, decided by an agent; carried over from an earlier run, which left it open.

RULE.DOC.SELF_CONTAINED: every symbol is bound in plain language at first use. The page's only bound symbols are $\hat\mu_n$ (the empirical measure), $g$ (the smooth functional), $m(x)$, $\sigma^2(x)$ and $h^\star(x)$ by apposition; $\tau$ is never said to be the time lag between the two members of a pair, $X$ and $Y$ are never said to be the conditioning variate and the response, and $K$ is never said to be the kernel (the sentence names 'smoothing kernels', then writes '$\int K = 1$' of 'any callable').

**Fix.** Bind each at first use: 'two-time transition density `transition_density` ($p_\tau(y\mid x)$, the density of the response $Y$ at lag $\tau$ after the conditioning value $x$)'; 'the conditional expectation ($\mathbb{E}[g(Y)\mid x]$, $Y$ the response and $x$ the conditioning value)'; 'any callable $K$ with $\int K = 1$'.

- `docs/src/index.md`:110  (GBGStreamingStatistics.jl / Main components)
- `docs/src/index.md`:101  (GBGStreamingStatistics.jl / Main components)
- `docs/src/index.md`:85  (GBGStreamingStatistics.jl / Main components)
- `docs/ARCHITECTURE.md`:395  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.3 Conditional and transition density)
- `docs/ARCHITECTURE.md`:551  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- `docs/ARCHITECTURE.md`:568  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- `docs/ARCHITECTURE.md`:734  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/api/ConditionalDensity.md`:3  (ConditionalDensity)
- `docs/src/index.md`:111  (GBGStreamingStatistics.jl / Main components)
- `docs/ARCHITECTURE.md`:384  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.2 Conditional expectation ``m(x)`` and variance ``\sigma^2(x)``)
- `docs/ARCHITECTURE.md`:939  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/index.md`:102  (GBGStreamingStatistics.jl / Main components)
- `docs/ARCHITECTURE.md`:330  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.1 Marginal density ``f_X``)
- `docs/ARCHITECTURE.md`:576  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- `docs/ARCHITECTURE.md`:729  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:729  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/api/Kernels.md`:4  (Kernels)
- `src/Bandwidth.jl`:538  (bandwidth)
- `src/ConditionalDensity.jl`:95  (conditional_density)
- `src/ConditionalDensity.jl`:95  (conditional_density)
- `src/ConditionalDensity.jl`:95  (conditional_density)
- `src/ConditionalExpectations.jl`:112  (conditional_expectation)
- `src/ConditionalExpectations.jl`:303  (conditional_mean)
- `src/ConditionalExpectations.jl`:387  (conditional_variance)
- `src/ConditionalExpectations.jl`:528  (conditional_pushforward_density)
- `src/ConditionalExpectations.jl`:547  (conditional_pushforward_density)
- `src/Expectations.jl`:308  (pushforward_density)
- `src/Expectations.jl`:319  (pushforward_density)
- `src/Kernels.jl`:19  (gaussian_kernel)
- `src/Kernels.jl`:19  (gaussian_kernel)
- `src/Kernels.jl`:30  (epanechnikov_kernel)
- `src/MarginalDensity.jl`:428  (density)
- `src/MarginalDensity.jl`:441  (density)
- `src/MarginalDensity.jl`:704  (density_standard_error)
- `docs/ARCHITECTURE.md`:730  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/notation.md`:46  (Notation / Smoothing kernels)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-092  The letter \ell is bound to two different things inside one display: the bin location \ell_j and the moment order (the subscript of T_\ell, P_\ell and the exponent), so the formula reads as if the exponent were a location.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@internals`, decided by an agent; carried over from an earlier run, which left it open.

The docstring binds 'bin locations $\ell_j$', then writes 'T_\ell &= \sum_j K_h(\ell_j - x)  w_j  (\ell_j - x)^\ell' and the normal equations '[ T_{i+\ell} ]_{i,\ell=0}^{d}', where \ell is now a running index from 0 to d. A first-time reader has to work out from the summation structure that the \ell in the exponent is not the \ell of \ell_j.

**Fix.** Rename one of the two: write the locations as $c_j$ (or $x_j$) and keep \ell as the moment order, or keep $\ell_j$ and index the moments by $k$ or $p$ ($T_k$, $(\ell_j - x)^k$, $[T_{i+k}]_{i,k=0}^d$). Use the same choice in the `locations` argument bullet.

!!! warning "The code may be what is wrong"
    The reader believed the documentation and the code disagree, and that the code is the wrong half. That is a question for the owner, not a documentation fix.

- `src/LocalPolynomial.jl`:76  (local_fit)
- `src/LocalPolynomial.jl`:81  (local_fit)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-15-documentation-review-01-097  Step 1 calls the variate of a one-dimensional marginal summary a 'conditioning variable', a role the page only introduces in step 4.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@tutorial`, decided by an agent; carried over from an earlier run, which left it open.

At step 1 nothing has been conditioned on anything; the reader asks 'conditioning what?' and finds the (X, Y) pairing three steps later. The same sentence calls Histogram1D 'the canonical summary' without saying canonical with respect to what.

**Fix.** In step 1 say 'of a scalar variate X' (or 'of the variate X'), and either drop 'canonical' or say in a clause why this summary is the reference one.

- `docs/src/tutorial.md`:17  (Getting started / The workflow)

### 2026-10-19-documentation-review-01-048  The conditional_mean docstring and ARCHITECTURE section 4.2 define the local-fit design and response moments as per-sample sums U_p = sum_i K_h(X_i - x)(X_i - x)^p and V_p = sum_i K_h(X_i - x)(X_i - x)^p Y_i, said to be 'recovered' from the summary or to be the Histogram2D's 'kernel-weighted power-sums'. The code computes something else: local_fit(c_x, N, sum_G, ...) uses the x-column bin centers, so U_p = sum_j N_j K_h(c^x_j - x)(c^x_j - x)^p. Every sample of column j sits at c^x_j. The x-direction power-sums S_10, S_20 and S_11 are never read. The x side is therefore an O(Delta_x) mass-at-center approximation, not an exact recovery. Only the response totals sum Y are exact.

Rule: `RULE.DOC.VERIFIED`.  
Reported by `full`, decided by an agent.

ConditionalExpectations.conditional_expectation calls _column_response, which returns c_x = [_xbincenter(summary, i_x) ...], and passes it as `locations` to local_fit. LocalPolynomial._design_moments uses locations[j] - x. Nothing in ConditionalExpectations reads powersums rows 1, 3 or 4 (S_10, S_20, S_11). local_fit's own docstring states the binned form correctly, with tilde c_j and n_j.

**Fix.** State the binned form the code computes, U_p = sum_j N_j K_h(c^x_j - x)(c^x_j - x)^p and V_p = sum_j K_h(c^x_j - x)(c^x_j - x)^p (sum_{i in column j} Y_i), and say the x locations are the column centers, an O(Delta_x) discretization. Alternatively, present the per-sample formula as the target and the binned one as what is evaluated, as conditional_density does. Drop 'kernel-weighted power-sums' from ARCHITECTURE.

- `src/ConditionalExpectations.jl`:264  (conditional_mean)
- `docs/ARCHITECTURE.md`:372  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.2 Conditional expectation ``m(x)`` and variance ``\sigma^2(x)``)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-049  Several places say the x-kernel weighting 'cancels' in the f_{X,Y}/f_X ratio and conclude that a bounded conditioning variable needs no asymmetric kernel_x. Only the kernel's normalization (the 1/h_x) cancels. The weighting itself sets the estimate. In x the double-kernel estimator is a degree-0 (Nadaraya-Watson-type) weighting, which per the package's own conditional_mean docstring 'keeps the O(h) edge bias' near a bounded edge. ARCHITECTURE 4.3 justifies the claim by citing the boundary adaptivity of the local-linear conditional_mean, which conditional_density does not have.

Rule: `RULE.DOC.VERIFIED`.  
Reported by `full`, decided by an agent.

ConditionalDensity._y_profile weights each column by _kernel_weight(c^x_j, x, kernel_x, h_x), and _density_at divides by D, so only the scale of kernel_x drops out. The code comment itself says 'so its normalization is irrelevant'. conditional_mean's docstring says: 'The degree = 0 Nadaraya-Watson fit does not - it keeps the O(h) edge bias'.

**Fix.** Say that the kernel's normalization cancels, not its weighting. Replace 'a bounded conditioning variable needs no asymmetric kernel_x' with an accurate caveat: the conditional density is degree-0 in x and keeps an O(h_x) bias at a bounded conditioning edge. Drop the appeal to conditional_mean's local-linear boundary adaptivity.

- `src/ConditionalDensity.jl`:113  (conditional_density)
- `src/ConditionalDensity.jl`:115  (conditional_density)
- `src/ConditionalDensity.jl`:135  (conditional_density)
- `docs/ARCHITECTURE.md`:405  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.3 Conditional and transition density)
- `docs/ARCHITECTURE.md`:729  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `src/ConditionalDensity.jl`:22  (ConditionalDensity)
- `src/ConditionalDensity.jl`:47  (ConditionalDensity)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-050  Several overview surfaces present conditional_variance (sigma^2(x)) as the 'y^2 case' of conditional_expectation, and ARCHITECTURE says the general reader 'subsumes both'. That is a category error. sigma^2(x) = E[Y^2 | X = x] - m(x)^2 is not E[g(Y) | X = x] for g(y) = y^2. The default estimator (method = :residual, Fan-Yao) is a second local fit of squared residuals and is not a conditional_expectation call at all. Only the :plugin method is built from the y^2 case, as a difference of two conditional expectations.

Rule: `RULE.DOC.VERIFIED`.  
Reported by `full`, decided by an agent.

ConditionalExpectations._conditional_variance with :residual fits _residual_response, built from _column_sums; :plugin computes ey2 - m^2. The conditional_expectation docstring states it correctly: 'for g(y) = y^2 it is E[Y^2 | X = x], the plug-in ingredient of conditional_variance'.

**Fix.** Say that the conditional mean is the g = identity case and that E[Y^2 | X = x] (the y^2 case) is the plug-in ingredient of the conditional variance. Do not call sigma^2(x) itself a case of the general reader.

- `docs/src/api/ConditionalExpectations.md`:8  (ConditionalExpectations)
- `docs/src/index.md`:104  (GBGStreamingStatistics.jl / Main components)
- `docs/ARCHITECTURE.md`:384  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.2 Conditional expectation ``m(x)`` and variance ``\sigma^2(x)``)
- `docs/ARCHITECTURE.md`:386  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.2 Conditional expectation ``m(x)`` and variance ``\sigma^2(x)``)
- `docs/ARCHITECTURE.md`:940  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:939  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/index.md`:102  (GBGStreamingStatistics.jl / Main components)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-051  The gamma_kernel docstring says the estimator's variance 'self-adapts, tightening toward the edge'. For Chen's gamma-kernel estimator the variance grows toward the boundary: it is O(n^-1 b^-1/2 x^-1/2) in the interior and O(n^-1 b^-1) within O(b) of 0. The package's own AsymmetricKernel docstring says the variance is proportional to w(x)^{-1/2} with w(x) = x, which diverges as x -> 0. What tightens toward the edge is the kernel's own spread (variance alpha b^2 = x b), not the estimator's variance.

Rule: `RULE.DOC.VERIFIED`.  
Reported by `full`, decided by an agent.

AsymmetricKernel docstring: 'the leading bias is 1/2 b w(x) f_X''(x) and the variance ∝ w(x)^{-1/2} (w(x) = x for the gamma kernel'. Bandwidth's Chen plug-in uses V = (2 sqrt(pi))^-1 E[w(X)^-1/2], so the variance weight grows as w -> 0.

**Fix.** Write 'the kernel narrows toward the edge (its spread is x b)' if that is the intended point, or state the estimator variance correctly: larger near 0, O(1/(n b)) against O(1/(n sqrt(b x))) in the interior.

- `src/Kernels.jl`:147  (gamma_kernel)
- `src/Kernels.jl`:86  (AsymmetricKernel)
- `src/Bandwidth.jl`:371  (bandwidth)
- `src/Kernels.jl`:95  (AsymmetricKernel)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-052  Two interface details stated in docstrings are not true of the methods. (a) observe! says it folds 'a vector, or a single sample' and that 'It runs identically on a CPU Array and any GPU backend'. The single-sample method is a separate host-only path: it scalar-indexes summary.counts and summary.powersums, so on a device summary it fails under allowscalar(false). It also does not use the scatter kernel the paragraph describes. (b) write_histogram lists the HDF5 attributes as `nbins` / `n_x` / `n_y`, but the extension writes the 2-D cell counts as "nx" / "ny" and the 2-D range edges as "xlower", "xupper", "ylower", "yupper". A reader opening the file with another HDF5 tool, as the docstring invites, will not find `n_x`.

Rule: `RULE.DOC.VERIFIED`.  
Reported by `full`, decided by an agent.

src/Histograms.jl, observe!(summary::Histogram1D, batch::Real): 'Scalar fold — a direct host update (scalar-indexes the summary, so host arrays only; batch on GPU)'. ext/GBGStreamingStatisticsHDF5Ext.jl write_histogram(::Histogram2D): a["nx"] = h.n_x; a["ny"] = h.n_y; a["xlower"] = h.lower_x.

**Fix.** (a) Say the single-sample form is a host-side convenience for CPU summaries and that a device summary takes a batch. (b) Name the attributes as the extension writes them: `lower`/`upper`/`nbins`/`order` for 1-D and `xlower`/`xupper`/`nx`/`ylower`/`yupper`/`ny`/`lag` for 2-D. Alternatively, rename the attributes in the extension, which is a schema change and needs the owner.

- `src/Accumulators.jl`:33  (observe!)
- `src/Accumulators.jl`:43  (observe!)
- `src/Histograms.jl`:414  (write_histogram)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-053  The status of the archived Kramers-Moyal module is stated three different ways, and the index's version is false. index.md says the Kramers-Moyal functionals 'the higher coefficients' and 'the Pawula / Markov diagnostics are archived — built and tested'. The archive (commit 37006b5, archive/README.md) records pawula_quotient and markov_test as 'declared-but-unimplemented' stubs and kramers_moyal as capped at n <= 2, so neither the higher coefficients nor the diagnostics were built. The tutorial's Status admonition says the same functionals 'are not currently provided but are planned', which contradicts 'archived — built'. index.md also advertises use in turbulence for 'drift and diffusion', which the package does not currently provide.

Rule: `RULE.DOC.VERIFIED`.  
Reported by `full`, decided by an agent.

git show 37006b5:archive/README.md: 'the declared-but-unimplemented `pawula_quotient` / `markov_test`', 'the Pawula/Markov diagnostics were stubs; `kramers_moyal` was capped at `n ≤ 2`'. ARCHITECTURE 11 (IncrementMoments archived) also says 'Pawula / Markov were stubs'.

**Fix.** Make index.md and the tutorial say the same thing: drift and diffusion (Kramers-Moyal orders 1-2) were built and OU-tested, then archived and are restorable; the higher coefficients and the Pawula/Markov diagnostics were never built (stubs) and wait on the Histogram2D order-raise. Drop 'drift and diffusion' from the list of current uses or mark it as archived.

- `docs/src/index.md`:43  (GBGStreamingStatistics.jl)
- `docs/src/tutorial.md`:12  (Getting started)
- `docs/src/index.md`:28  (GBGStreamingStatistics.jl)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-054  The notation page's bindings are broken both ways. (1) One concept, two symbols. The page says an estimate is the hatted estimand symbol (hat f_X, hat m(x)), yet the density estimate is also written hat f_h(x), in the same docstring as hat f_X(x) in density. The conditional-mean estimate is also written hat m_h(x) beside hat m. In scale! the cone scalar appears as texttt{weight} in one formula and lambda in the next, though the page assigns lambda. (2) One symbol, two concepts. The page reserves alpha(x) for the gamma kernel's shape, but gamma_start writes hat alpha for the gamma start's shape and ARCHITECTURE writes alpha for a gamma-reference density's shape. The doctests name a summary `h`, which is the bandwidth symbol: density's doctest calls `density(h, 0.4; h = 0.2)`, with `h` the histogram and `h =` the bandwidth on one line. The conditional_density doctest names a joint summary `g`, which is the response functional symbol.

Rule: `RULE.DOC.ONE_NOTATION`.  
Reported by `full`, decided by an agent.

docs/src/notation.md: 'A density symbol always carries its subscript, and an estimate its hat: hat f_X, hat m(x)'; 'alpha(x), w(x) | the shape of the gamma kernel'; 'lambda | the scalar of the cone: rescaling, weighted merging'; 'h | a bandwidth'; 'E[g(X)] | the expectation of a smooth functional g'.

**Fix.** Write hat f_X(x) and hat m(x) for the estimates. Where bandwidth dependence must be shown, add hat f_h to the notation page as an explicitly defined alias, or write hat f_X(x; h). Give the gamma start's shape and the gamma-reference shape another letter, or add a row. Use lambda consistently in scale! after first binding it to `weight`. Rename the doctest summaries to `summary` / `joint`, as the other doctests and the tutorial do.

- `src/MarginalDensity.jl`:419  (density)
- `src/MarginalDensity.jl`:696  (density_standard_error)
- `docs/ARCHITECTURE.md`:330  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.1 Marginal density ``f_X``)
- `src/ConditionalExpectations.jl`:269  (conditional_mean)
- `src/ConditionalExpectations.jl`:378  (conditional_variance)
- `src/Bandwidth.jl`:508  (bandwidth)
- `src/Accumulators.jl`:125  (scale!)
- `src/MarginalDensity.jl`:119  (gamma_start)
- `docs/ARCHITECTURE.md`:798  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `src/MarginalDensity.jl`:502  (density)
- `src/MarginalDensity.jl`:728  (density_standard_error)
- `src/MarginalDensity.jl`:728  (density_standard_error)
- `src/Expectations.jl`:205  (expectation)
- `src/Expectations.jl`:273  (expectation_with_standard_error)
- `src/ConditionalDensity.jl`:146  (conditional_density)
- `src/Bandwidth.jl`:317  (bandwidth)
- `src/MarginalDensity.jl`:699  (density_standard_error)
- `src/MarginalDensity.jl`:702  (density_standard_error)
- `src/MarginalDensity.jl`:706  (density_standard_error)
- `docs/ARCHITECTURE.md`:576  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- `docs/ARCHITECTURE.md`:729  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:730  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/api/Kernels.md`:4  (Kernels)
- `docs/src/index.md`:85  (GBGStreamingStatistics.jl / Main components)
- `docs/src/notation.md`:46  (Notation / Smoothing kernels)
- `src/Bandwidth.jl`:538  (bandwidth)
- `src/ConditionalDensity.jl`:95  (conditional_density)
- `src/ConditionalExpectations.jl`:112  (conditional_expectation)
- `src/ConditionalExpectations.jl`:303  (conditional_mean)
- `src/ConditionalExpectations.jl`:387  (conditional_variance)
- `src/ConditionalExpectations.jl`:528  (conditional_pushforward_density)
- `src/ConditionalExpectations.jl`:547  (conditional_pushforward_density)
- `src/Expectations.jl`:308  (pushforward_density)
- `src/Expectations.jl`:319  (pushforward_density)
- `src/Kernels.jl`:19  (gaussian_kernel)
- `src/Kernels.jl`:30  (epanechnikov_kernel)
- `src/MarginalDensity.jl`:428  (density)
- `src/MarginalDensity.jl`:441  (density)
- `src/MarginalDensity.jl`:704  (density_standard_error)
- `src/ConditionalExpectations.jl`:273  (conditional_mean)
- `src/ConditionalExpectations.jl`:281  (conditional_mean)
- `src/ConditionalExpectations.jl`:282  (conditional_mean)
- `src/Accumulators.jl`:116  (scale!)
- `src/Accumulators.jl`:117  (scale!)
- `src/MarginalDensity.jl`:106  (lognormal_start)
- `src/MarginalDensity.jl`:788  (density_with_standard_error)

### 2026-10-19-documentation-review-01-055  Two defining formulas live outside their operation's home docstring. (1) The full defining formula of the joint-summary selector h*(x) = [R(K) sigma^2(x)/(mu_2(K)^2 m''(x)^2 f_X(x) n)]^{1/5} is restated in the docstring of the Histogram1D bandwidth method, as well as in its own bandwidth(::Histogram2D, x) docstring. The pointer should carry only a link and a one-line meaning. (2) The density estimator's asymmetric-kernel branch, hat f_X(x) = (1/n) sum_j n_j K(c_j; x, b), is stated in the gamma_kernel docstring and again in the AsymmetricKernel docstring, but not in density's own docstring. density's formulas cover only the symmetric (n_j/h_j) K((x - tilde c_j)/h_j) form, so the operation's home lacks one of its estimators and two other docstrings each carry a copy.

Rule: `RULE.DOC.ONE_HOME`.  
Reported by `full`, decided by an agent.

Bandwidth.jl bandwidth(::Histogram1D) docstring, the paragraph before # Arguments. Kernels.jl gamma_kernel and AsymmetricKernel (pdf bullet). MarginalDensity.jl density, whose degree bullet holds the only estimator formula and is symmetric.

**Fix.** In bandwidth(::Histogram1D) replace the h*(x) formula with a one-line pointer to the joint-summary method. Move the asymmetric-kernel estimator formula into density's kernel bullet, and leave the kernel docstrings with a link to density.

- `src/Bandwidth.jl`:389  (bandwidth)
- `src/Kernels.jl`:145  (gamma_kernel)
- `src/Kernels.jl`:81  (AsymmetricKernel)

### 2026-10-19-documentation-review-01-056  The umbrella module docstring, the package's map, does not list one linked one-liner per operation family. The accumulation verbs (observe!, observe_pairs!, merge!, scale!) appear only as adjectives ('mergeable', 'rescalable') under the Histogram entry and are never linked. The conditional mean and variance, and the recovered moments (mean/var/std/central_moment), are named in prose without links. conditional_pushforward_density is absent. Two exported families have no entry at all: the smoothing kernels (gaussian_kernel, epanechnikov_kernel, gamma_kernel, beta_kernel, AsymmetricKernel) and the parametric starts (ParametricStart and the five *_start values).

Rule: `RULE.DOC.MODULE_MAP`.  
Reported by `full`, decided by an agent.

src/GBGStreamingStatistics.jl module docstring, against its own export list.

**Fix.** Add linked one-line entries for the accumulation operations, the smoothing kernels and the parametric starts. Link conditional_mean, conditional_variance, central_moment and conditional_pushforward_density where they are named, or name the family by its linked home.

- `src/GBGStreamingStatistics.jl`:22  (GBGStreamingStatistics)
- `src/GBGStreamingStatistics.jl`:29  (GBGStreamingStatistics)
- `src/GBGStreamingStatistics.jl`:34  (GBGStreamingStatistics)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-057  The AsymmetricKernel docstring is out of order. Its `# Fields` section comes second, right after the identity, ahead of `# Arguments`, the usage example and `# References`. The rule puts `# References` after the prose and `# Fields` last. The other type docstrings (Histogram1D, Histogram2D, ParametricStart) follow the order, with Fields last.

Rule: `RULE.DOC.HOME_ORDER`.  
Reported by `full`, decided by an agent.

src/Kernels.jl AsymmetricKernel: identity, then '# Fields', then '# Arguments', then the julia example, then '# References'.

**Fix.** Move the `# Fields` section to the end of the docstring, after `# References`.

- `src/Kernels.jl`:69  (AsymmetricKernel)
- `src/Histograms.jl`:70  (Histogram1D)
- `src/Histograms.jl`:245  (Histogram2D)
- `src/MarginalDensity.jl`:74  (ParametricStart)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-058  On the built API pages, six docstring example code blocks scroll sideways. Their longest lines (91 to 108 characters, almost all because of a trailing comment) are wider than the 784 px code column. The end of each comment, and in places the code, is hidden until the reader scrolls the block.

Rule: `RULE.TEX.NO_HSCROLL`.  
Reported by `rendered`, decided by an agent.

I rendered the built pages in headless Chromium at 1400x900 with KaTeX and the site CSS loaded. The column is 784 px, set in 14 px JuliaMono at about 8.4 px per character plus 24 px padding, so 90 characters is the limit. Measured overflow (scrollWidth > clientWidth) of the code element: api/Expectations, two blocks, 830>784 and 822>784; api/Histograms 822>784; api/Kernels 906>784; api/MarginalDensity 931>784 and 898>784. No other page overflows. Every display formula fits, the widest being the ConditionalDensity variance formula at 0.89 of the column, and none packs equations behind \quad. Tables and the notation page fit too. The tutorial's 93-character line sits in the wider 812 px page column and fits.

**Fix.** Move each trailing comment onto its own line above the input, or shorten it, so that every line is 90 characters or fewer. In a jldoctest, a full-line comment can go in its own `julia>` input or in plain prose above the block. In the AsymmetricKernel example, break the long `#` comment over two lines.

- `src/Expectations.jl`:207  (expectation)
- `src/Expectations.jl`:275  (expectation_with_standard_error)
- `src/Histograms.jl`:239  (Histogram2D)
- `src/Histograms.jl`:241  (Histogram2D)
- `src/Kernels.jl`:94  (AsymmetricKernel)
- `src/MarginalDensity.jl`:502  (density)
- `src/MarginalDensity.jl`:64  (ParametricStart)
- `src/MarginalDensity.jl`:63  (ParametricStart)

### 2026-10-19-documentation-review-01-059  The page uses 'the lower half' and 'the upper half' as package terms of art without ever saying what is being halved or what the split means. A reader seeing only this page cannot tell that it is the package's two-layer architecture (the accumulation layer and the estimation layer), and the page links no definition.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api`, decided by an agent.

Neither 'half' is introduced earlier on the page, and the page links no definition. The only context the page gives is the module lists inside each 'half'.

**Fix.** Define the split in a clause before using it, e.g. 'The package has two layers: a lower half that accumulates the stream into summaries, and an upper half that reads estimates off them.' Otherwise drop 'half' and say 'the accumulation modules' / 'the estimator modules'.

- `docs/src/api.md`:8  (API overview)
- `docs/src/api.md`:10  (API overview)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-060  'the stream' and 'the samples' are definite references with nothing earlier on the page to point back to. The page never says that the package reads a stream of samples once and does not store them. That fact is what 'read in place of the samples' depends on, so a first-time reader meets those words without their premise.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api`, decided by an agent.

This is the first and only use of 'stream' and 'samples' on the page, and both take 'the' as if already introduced.

**Fix.** Bind them in plain language first, e.g. '... from the accumulated summary: a binned record of a stream of samples, each observed once and discarded, which the estimators read in place of the samples.'

- `docs/src/api.md`:5  (API overview)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-061  The upper-half modules are said to 'hold the estimands', which turns the page's own term around. In the previous sentence an estimand is the population quantity and the estimator approximates it. A module holds estimators, so to a reader the sentence misuses the term the page has just bound.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api`, decided by an agent.

The first sentence defines the estimand as 'the population quantity' and distinguishes it from 'the estimator that approximates it'. The last sentence then says code modules hold estimands.

**Fix.** Say the upper half 'holds the estimators' (or 'the estimators of those estimands').

- `docs/src/api.md`:12  (API overview)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-062  Two sentences have no subject or no main verb: the AbstractAccumulator identity is a bare noun phrase, and the scale! requirement opens on a verb with no subject.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-Accumulators`, decided by an agent.

Read aloud, 'Supertype of every mergeable summary ... of an observed stream ..., which gives a set A the share of samples in it' never reaches a main verb; 'Requires a Float-count summary ... — an `ArgumentError` otherwise.' has no subject.

**Fix.** 'Every mergeable summary ... subtypes `AbstractAccumulator`.' / 'It requires a Float-count summary ...; otherwise it throws an `ArgumentError`.'

- `src/Accumulators.jl`:15  (AbstractAccumulator)
- `src/Accumulators.jl`:122  (scale!)
- `docs/ARCHITECTURE.md`:73  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The spine: empirical measure → summary → post-hoc-smoothed estimator)
- `docs/ARCHITECTURE.md`:87  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The spine: empirical measure → summary → post-hoc-smoothed estimator)
- `docs/src/notation.md`:36  (Notation / The summary)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-063  Symbols in merge! are used without a binding: the operator \oplus is never said to be the merge, T is bound on the page as a map of the measure n P_n but is here applied to sample blocks and only called 'accumulation T' after the formula, and f_\theta (the per-parameter density) is never bound at all.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-Accumulators`, decided by an agent.

The display defines f_X as a mixture of f_\theta, and the following sentence binds \theta and \pi(\theta) but not f_\theta; \oplus appears only in the formula and 'to summaries (under `merge!`)' never names it.

**Fix.** Write '... = T(\mathcal{A}) \oplus T(\mathcal{B}), with \oplus the merge' (or write merge! itself), say T maps a sample block to its summary before the formula, and bind f_\theta as 'the density of the samples drawn at parameter \theta'.

- `docs/src/api/Accumulators.md`  (merge!)  — **the text here has moved and was not found again**
- `src/Accumulators.jl`:80  (merge!)
- `docs/ARCHITECTURE.md`:816  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:819  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/notation.md`:96  (Notation / Parameters)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-064  The package's own type-parameter vocabulary is used with no definition or link on this page: 'Float-count' and the keyword `count_type` (a constructor keyword of the histogram types, which this page never shows), and 'the summary's accumulator type'.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Accumulators`, decided by an agent.

A cold reader cannot tell where `count_type` is set or what an 'accumulator type' is as distinct from the count type; neither is linked to the Histogram1D / Histogram2D constructors that own them.

**Fix.** Gloss once: 'a summary whose counts are floating point (built with the keyword `count_type = Float64` of [`Histogram1D`](@ref ...))', and say the accumulator type is the element type of the power-sums.

- `src/Accumulators.jl`:85  (merge!)
- `src/Accumulators.jl`:44  (observe!)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-065  The merge! hyperstatistics paragraph uses package nouns no sentence on the page defines: 'readouts', 'a joint summary', and 'the local mixing weight'.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Accumulators`, decided by an agent.

The page has introduced only one-dimensional summaries of a single stream; a 'joint summary' and its 'readouts', and what a weight 'local' in x would be, assume the module's context.

**Fix.** Name and link the joint summary (Histogram2D, the summary of pairs), and say what the local mixing weight is (the mixing weight reweighted by each component's share of mass near the conditioning value), or cut the parenthesis.

- `src/Accumulators.jl`:83  (merge!)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-066  'Kernel' carries two meanings on adjacent paragraphs with no disambiguation: the smoothing kernel of a density estimator, then the GPU compute kernel that does the binning.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Accumulators`, decided by an agent.

The preceding paragraph says accumulation 'knows nothing ... of the bandwidth or kernel' (the smoothing kernel), and the next sentence says the binning is done by 'one KernelAbstractions + Atomix kernel', which a statistician reads as a smoothing kernel contradicting the line before.

**Fix.** Say 'one KernelAbstractions + Atomix compute kernel (a GPU-style parallel loop)' or 'one parallel scatter-add launch', and 'smoothing kernel' in the earlier sentence.

- `src/Accumulators.jl`:40  (observe!)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-067  Package slogans stand in for definitions: 'estimand-agnostic' is used in AbstractAccumulator before the page explains it (the explanation comes only in observe!), and observe! appeals to 'the architectural invariant that delivers no reruns' as if the reader knew the package's invariants.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Accumulators`, decided by an agent.

On a top-to-bottom read the first 'estimand-agnostic' has no gloss; 'reruns' and 'the architectural invariant' name things the page never introduced.

**Fix.** Define the term at its first use ('estimand-agnostic: it does not depend on which quantity will later be estimated from it') and drop the 'architectural invariant ... no reruns' phrasing for the plain clause that follows it.

- `src/Accumulators.jl`:22  (AbstractAccumulator)
- `src/Accumulators.jl`:37  (observe!)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-068  Several sentences have no main verb: both 'Estimand.' paragraphs are a noun phrase extended by a dash clause, and the module docstring's 'For the marginal (...): ...' and 'For the joint (...): ...' sentences are verbless lists of selectors.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-Bandwidth`, decided by an agent.

'The AMISE-optimal ... smoothing bandwidth for the marginal density ... for the smoothing kernel K and local-polynomial degree d.' has no finite verb ('minimising' sits in a participial modifier). The same holds for 'The response-aware local AMSE-optimal bandwidth h*(x) for the conditional mean ... — the bandwidth that ... minimises ...', where the only verb belongs to the relative clause. 'For the joint (Histogram2D): the response-aware local plug-in h*(x) (rule = :local).' has no verb either.

**Fix.** Give each one a verb, e.g. 'Estimand. The bandwidth estimated is the AMISE-optimal ...'; 'For the joint summary, the selector is the response-aware local plug-in ...'.

- `src/Bandwidth.jl`:315  (bandwidth)
- `src/Bandwidth.jl`:505  (bandwidth)
- `docs/src/api/Bandwidth.md`:4  (Bandwidth)
- `docs/src/api/Bandwidth.md`:7  (Bandwidth)
- `docs/ARCHITECTURE.md`:162  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / Estimand families (the operation taxonomy))
- `docs/ARCHITECTURE.md`:182  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The histogram is the summary; the ECF is deferred (likely never needed))
- `docs/ARCHITECTURE.md`:423  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:429  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:865  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:873  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:918  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/api/Bandwidth.md`:8  (Bandwidth)
- `docs/src/index.md`:39  (GBGStreamingStatistics.jl)
- `docs/src/index.md`:115  (GBGStreamingStatistics.jl / Main components)
- `docs/src/notation.md`:52  (Notation / Smoothing kernels)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-069  Some symbols are used before the page binds them in words. In the first docstring, f_X and \hat f_h appear in the MISE formula but neither is named as the true density of X or the bandwidth-h estimate, and X is bound only later, in the Arguments bullet. In the second docstring, \hat\sigma_X is used in the prose before its binding in the pilot_bandwidth bullet.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-Bandwidth`, decided by an agent.

The Estimand gives only 'the marginal density density of summary', then '\operatorname{MISE}(\hat f_h) = \int \mathbb{E}[(\hat f_h(x) - f_X(x))^2] dx'. The Histogram2D docstring's 'the default pilot is the curvature-rate rule of thumb \hat\sigma_X n^{-1/9}' comes before 'the conditioning variable's sample standard deviation'.

**Fix.** Bind the symbols at first use: 'where f_X is the density of X and \hat f_h its kernel estimate at bandwidth h'; '\hat\sigma_X n^{-1/9}, with \hat\sigma_X the sample standard deviation of X'.

- `src/Bandwidth.jl`:317  (bandwidth)
- `src/Bandwidth.jl`:525  (bandwidth)
- `docs/ARCHITECTURE.md`:802  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-070  The page refers to a constant for degree 2/3 that it never states. The prose promises 'at degree 2/3 the analogous fourth-order equivalent-kernel one', but the displayed C(K, d) = (8 sqrt(pi)/3)^{1/5} delta(K) does not depend on d, so a reader cannot get the degree-2/3 constant from the page.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-Bandwidth`, decided by an agent.

The display binds C(K, d) through d-free expressions only. Nothing on the page gives the degree-2/3 factor, so the reader has to go to Loader (1996).

**Fix.** State the degree-2/3 constant in the display, e.g. as a case of C(K, d), or say plainly that the display is the degree-0/1 constant.

- `src/Bandwidth.jl`:330  (bandwidth)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-071  The page uses the package's own nouns ('summary', 'accumulation', 'the marginal', 'the joint') as if the reader already knew them. The page never says that a summary is the accumulated binned histogram of the sample stream, or what accumulation is.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Bandwidth`, decided by an agent.

The module docstring opens with 'For the marginal (Histogram1D)' and 'For the joint (Histogram2D)' without saying marginal or joint what. The first docstring says 'read back from the summary's stored moments' and 'Every input is read back from the one stored summary after accumulation.'

**Fix.** Define the term in a clause at first use, e.g. 'the summary (the binned histogram accumulated from the sample stream, Histogram1D)', or link it to its definition.

- `src/Bandwidth.jl`:343  (bandwidth)
- `docs/src/api/Bandwidth.md`:3  (Bandwidth)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-072  Two abbreviations are used without being expanded on the page: 'ISE', although MISE and AMISE are expanded, and 'eval kernel', a clipped form of the 'evaluation kernel' named earlier.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Bandwidth`, decided by an agent.

'It lowers the ISE by approx 1.4-1.7x on skewed positive data'. In the AIC block: 'the eval kernel enters only as the canonical delta(K) rescale'.

**Fix.** Write 'integrated squared error (ISE)' and 'evaluation kernel'.

- `src/Bandwidth.jl`:380  (bandwidth)
- `src/Bandwidth.jl`:363  (bandwidth)
- `docs/ARCHITECTURE.md`:419  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:740  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:749  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:755  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:770  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/notation.md`:60  (Notation / Smoothing kernels)
- `docs/src/notation.md`:61  (Notation / Smoothing kernels)
- `src/Bandwidth.jl`:335  (bandwidth)
- `src/Bandwidth.jl`:336  (bandwidth)
- `src/Bandwidth.jl`:350  (bandwidth)
- `src/Bandwidth.jl`:360  (bandwidth)
- `src/Bandwidth.jl`:402  (bandwidth)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-073  The Histogram2D signature's arrow '-> h' names a single bandwidth, but the result's type depends on x: a scalar for a scalar x, a vector for a vector x. A later sentence has to say so. Here the arrow is not fixed by the signature, and the result ends up stated twice.

Rule: `RULE.DOC.RETURNS`.  
Reported by `withheld@api-Bandwidth`, decided by an agent.

The signature reads 'bandwidth(summary::Histogram2D, x; ...) -> h', and the prose adds 'A vector `x` returns a vector; `NaN` marks an $x$ off the data's support.'

**Fix.** Drop the arrow and state the result in one 'Returns ...' sentence (the bandwidth h*(x), a vector for a vector x, NaN off support), or make the arrow cover both shapes.

!!! warning "The code may be what is wrong"
    The reader believed the documentation and the code disagree, and that the code is the wrong half. That is a question for the owner, not a documentation fix.

- `src/Bandwidth.jl`:527  (bandwidth)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-074  The sample size n is used in formulas and rates on this page but never bound to a name: the conditional_density docstring says the pairs' 'count is the sample size' without calling it n, then uses n in the rate, and the standard-error docstring uses n in its variance display and in D = n \hat f_X(x) with no binding at all.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-ConditionalDensity`, decided by an agent.

A cold reader meets $O(n^{-2/3})$ and $\frac{R(K_x)}{n h_x f_X(x)}$ with n defined nowhere; the only nearby lowercase-n symbol is the bin count $n_{jl}$, which invites the wrong reading.

**Fix.** Bind it where the sample size is first mentioned, e.g. 'the sums running over the n observed pairs (X_i, Y_i)', and in the standard-error docstring 'n the number of observed pairs'.

- `src/ConditionalDensity.jl`:96  (conditional_density)
- `src/ConditionalDensity.jl`:126  (conditional_density)
- `src/ConditionalDensity.jl`:254  (conditional_density_standard_error)
- `docs/ARCHITECTURE.md`:733  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-075  Symbols used before, or without, a binding: the lag tau appears in $p_\tau(y\mid x)$ in the page's opening paragraph and is bound only in the transition_density docstring much further down; the standard-error docstring speaks of an $O(h^2)$ bias where only h_x and h_y exist.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-ConditionalDensity`, decided by an agent.

Top to bottom, tau is first read in the intro with no gloss; 'h' alone is bound nowhere on the page.

**Fix.** Gloss tau in the intro ('the density of the state a lag tau later'); write the bias as $O(h_x^2 + h_y^2)$ or name h.

- `docs/src/api/ConditionalDensity.md`:3  (ConditionalDensity)
- `src/ConditionalDensity.jl`:269  (conditional_density_standard_error)
- `docs/ARCHITECTURE.md`:395  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.3 Conditional and transition density)
- `docs/ARCHITECTURE.md`:551  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- `docs/ARCHITECTURE.md`:568  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- `docs/ARCHITECTURE.md`:734  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/index.md`:111  (GBGStreamingStatistics.jl / Main components)
- `docs/ARCHITECTURE.md`:115  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The algebra of measures — a structural lens)
- `docs/ARCHITECTURE.md`:157  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / Estimand families (the operation taxonomy))
- `docs/ARCHITECTURE.md`:394  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.3 Conditional and transition density)
- `docs/ARCHITECTURE.md`:399  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.3 Conditional and transition density)
- `docs/ARCHITECTURE.md`:403  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.3 Conditional and transition density)
- `docs/ARCHITECTURE.md`:548  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- `docs/ARCHITECTURE.md`:728  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:730  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/index.md`:110  (GBGStreamingStatistics.jl / Main components)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-076  Verbless sentences and clauses: the standard-error 'Estimator.' paragraph opens with a noun phrase and no verb ('The two-term sqrt(...) from the conditional density ..., the bandwidths ..., and the x-window denominator D ...'); 'Symmetric kernels only (...)' stands alone as a sentence; and in the conditional_density estimator sentence the clause after the semicolon, 'both $\int K = 1$, default gaussian_kernel, the sums running over ...', is a fragment.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-ConditionalDensity`, decided by an agent.

Read aloud, none of these has a main verb; the first runs some fifty words before its full stop without one.

**Fix.** E.g. 'Estimator. The standard error is estimated by the two-term sqrt(...), built from ...'; 'Only symmetric kernels are supported (...)'; 'both integrate to 1 and default to gaussian_kernel; the sums run over the n observed pairs'.

- `src/ConditionalDensity.jl`:262  (conditional_density_standard_error)
- `src/ConditionalDensity.jl`:270  (conditional_density_standard_error)
- `src/ConditionalDensity.jl`:95  (conditional_density)
- `docs/ARCHITECTURE.md`:330  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.1 Marginal density ``f_X``)
- `docs/ARCHITECTURE.md`:576  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- `docs/ARCHITECTURE.md`:729  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:730  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/api/Kernels.md`:4  (Kernels)
- `docs/src/index.md`:85  (GBGStreamingStatistics.jl / Main components)
- `docs/src/notation.md`:46  (Notation / Smoothing kernels)
- `src/Bandwidth.jl`:538  (bandwidth)
- `src/ConditionalExpectations.jl`:112  (conditional_expectation)
- `src/ConditionalExpectations.jl`:303  (conditional_mean)
- `src/ConditionalExpectations.jl`:387  (conditional_variance)
- `src/ConditionalExpectations.jl`:528  (conditional_pushforward_density)
- `src/ConditionalExpectations.jl`:547  (conditional_pushforward_density)
- `src/Expectations.jl`:308  (pushforward_density)
- `src/Expectations.jl`:319  (pushforward_density)
- `src/Kernels.jl`:19  (gaussian_kernel)
- `src/Kernels.jl`:30  (epanechnikov_kernel)
- `src/MarginalDensity.jl`:428  (density)
- `src/MarginalDensity.jl`:441  (density)
- `src/MarginalDensity.jl`:704  (density_standard_error)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-077  The telegraphic parenthetical '(nonstationary: $p_{t,\tau}$)' has no verb and introduces a new symbol whose meaning, the transition density depending on the start time t too when the process is not stationary, the reader must reconstruct.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-ConditionalDensity`, decided by an agent.

It interrupts the estimand sentence as a label plus a symbol, with no clause saying what changes when the process is nonstationary.

**Fix.** E.g. '(for a nonstationary process it also depends on t, written $p_{t,\tau}$)'.

- `src/ConditionalDensity.jl`:185  (transition_density)
- `docs/ARCHITECTURE.md`:304  (GBGStreamingStatistics — Architecture and Rationale / 3. The accumulated summary — the universal lower half)
- `docs/ARCHITECTURE.md`:395  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.3 Conditional and transition density)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-078  A comparison that compares a rate with a density, followed by a parenthetical rate that, by position, attaches to the marginal density rather than to the conditional estimator.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-ConditionalDensity`, decided by an agent.

'The conditional MISE rate is slower than a marginal density ($O(n^{-2/3})$)' reads as though $O(n^{-2/3})$ were the marginal density's rate; the intended meaning is the reverse.

**Fix.** E.g. 'The conditional MISE converges at $O(n^{-2/3})$, slower than a marginal density's $O(n^{-4/5})$.'

- `src/ConditionalDensity.jl`:126  (conditional_density)
- `docs/ARCHITECTURE.md`:733  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-079  'the summary stores the lag only as a label' refers to a lag label on the summary that this page never introduces: the Histogram2D example constructor takes no lag, and nothing says where such a label lives or how it is set.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-ConditionalDensity`, decided by an agent.

A cold reader is told the summary stores a label but has no way to tell whether one exists, or whether the sentence means the lag is the caller's bookkeeping alone.

**Fix.** Say what is meant, e.g. 'the summary does not record the lag; the caller keeps track of which lag a histogram holds'.

- `src/ConditionalDensity.jl`:191  (transition_density)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-080  Several opening sentences are noun phrases with no main verb: the page's lead sentence, most 'Estimand.' paragraphs, the openers of two paired readers, and conditional_mean's 'Consistency:' line, which is only a formula.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent.

The page's first sentence lists items ('The conditional expectation ... with its sampling standard error ..., the conditional mean ... and the conditional pushforward density ...') and never reaches a verb. Each 'Estimand.' paragraph is an appositive noun phrase ending in a formula. conditional_pushforward_density_standard_error shows the full-sentence form ('The estimand is the pointwise standard error ...'), and the others do not follow it.

**Fix.** Give each one a verb, e.g. 'The module estimates the conditional expectation ...'; 'Estimand. The estimand is the conditional variance ...'; 'Returns the conditional variance and its standard error as the NamedTuple ...'; 'Consistency: the bias is approximately ... and the variance approximately ...'.

- `docs/src/api/ConditionalExpectations.md`:3  (ConditionalExpectations)
- `src/ConditionalExpectations.jl`:83  (conditional_expectation)
- `src/ConditionalExpectations.jl`:165  (conditional_expectation_standard_error)
- `src/ConditionalExpectations.jl`:257  (conditional_mean)
- `src/ConditionalExpectations.jl`:281  (conditional_mean)
- `src/ConditionalExpectations.jl`:364  (conditional_variance)
- `src/ConditionalExpectations.jl`:434  (conditional_variance_standard_error)
- `src/ConditionalExpectations.jl`:478  (conditional_variance_with_standard_error)
- `src/ConditionalExpectations.jl`:514  (conditional_pushforward_density)
- `src/ConditionalExpectations.jl`:653  (conditional_pushforward_density_with_standard_error)
- `docs/src/notation.md`:16  (Notation / The samples)
- `src/ConditionalExpectations.jl`:543  (conditional_pushforward_density)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-081  The term 'jet' first appears in conditional_expectation's degree bullet without a definition. The page defines it only further down, in conditional_mean ('the jet $(\hat\beta_0, \dots)$ — the local-polynomial coefficients at $x$'). conditional_variance and conditional_pushforward_density use it again.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent.

A first-time reader meets 'which fits no polynomial — no jet' before any sentence says what a jet is. In this package 'jet' means the vector of local-polynomial coefficients, which is not the field's usual word for it.

**Fix.** At the first use, gloss it: '(which fits no polynomial, so it has no jet — no local-polynomial coefficients beyond the fitted value)'. Or drop the word where 'fits no polynomial' already says enough.

- `src/ConditionalExpectations.jl`:118  (conditional_expectation)
- `src/ConditionalExpectations.jl`:390  (conditional_variance)
- `src/ConditionalExpectations.jl`:308  (conditional_mean)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-082  The symbol $\Delta_y$ (presumably the response-bin width) appears in three docstrings, and nothing on the page binds it.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent.

'$O(\Delta_y^2)$ for a general smooth $g$' in conditional_expectation; '$O(\Delta_y^2)$ quadrature' in conditional_variance_standard_error; '$h_z \gtrsim |g'| \Delta_y$' in conditional_pushforward_density. None of them says what $\Delta_y$ is.

**Fix.** Bind it at first use: '... and $O(\Delta_y^2)$, with $\Delta_y$ the $y$-bin width of the summary, for a general smooth $g$'.

- `src/ConditionalExpectations.jl`:101  (conditional_expectation)
- `src/ConditionalExpectations.jl`:445  (conditional_variance_standard_error)
- `src/ConditionalExpectations.jl`:535  (conditional_pushforward_density)
- `docs/ARCHITECTURE.md`:947  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-083  conditional_mean's estimator and consistency formulas use $K_h$, $\mu_2(K)$, $R(K)$ and $n$ without binding them. $\mu_2(K)$, $K_h$ and $n$ are bound nowhere on the page. $R(K)$ is bound only later, in the pushforward standard-error docstring. The pushforward docstring likewise uses $K_{h_z}$ without stating the scaling.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent.

'$U_p(x) = \sum_i K_h(X_i-x)(X_i-x)^p$' never says $K_h(u) = K(u/h)/h$. 'Consistency: ... \mu_2(K) ... R(K) \sigma^2(x)/(n h f_X(x))' introduces two kernel constants and the sample size bare.

**Fix.** Add one binding clause, e.g. 'with $K_h(u) = K(u/h)/h$, $\mu_2(K) = \int u^2 K$, $R(K) = \int K^2$ and $n$ the sample count'.

- `src/ConditionalExpectations.jl`:264  (conditional_mean)
- `src/ConditionalExpectations.jl`:282  (conditional_mean)
- `src/ConditionalExpectations.jl`:523  (conditional_pushforward_density)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-084  conditional_pushforward_density names the arguments `kernel_y` and `normalize` as if the reader knew them. Neither is an argument of any function on this page; both belong to conditional_density on another page.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent.

'the native estimator for the untransformed response (with the bounded-response asymmetric `kernel_y` and `normalize`)' — no signature on this page has either name, and the parenthesis does not say whose options they are.

**Fix.** Attribute them: '(which offers an asymmetric response kernel and exact normalization)'. Or link them to conditional_density's arguments explicitly.

- `src/ConditionalExpectations.jl`:534  (conditional_pushforward_density)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-085  The pushforward docstrings use package terms nowhere defined on the page: 'the within-cell correction', 'the two-term' band, and the optics word 'caustic' for $g'=0$.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent.

'agrees with `conditional_density` to second order in the response bin width (the within-cell correction)' does not say what is corrected. 'The band is the two-term $\sqrt{...}$' uses 'two-term' as a known label. 'caustic' is not standard statistics vocabulary.

**Fix.** Gloss each one: 'the difference is the within-cell spread of the response, which the binned double kernel ignores'; 'the two-term variance (the kernel term minus the squared-density term)'; '$g'=0$, where the change-of-variables Jacobian diverges, is integrable'.

- `src/ConditionalExpectations.jl`:533  (conditional_pushforward_density)
- `src/ConditionalExpectations.jl`:527  (conditional_pushforward_density)
- `src/ConditionalExpectations.jl`:615  (conditional_pushforward_density_standard_error)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-086  The page's lead paragraph speaks of a 'per-column' expectation 'read off the stored response power-sums' before the page has said what is stored, or what the columns are. The two-dimensional binned summary is first named, and its columns first called conditioning bins, in conditional_expectation's docstring.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent.

The module paragraph never mentions the Histogram2D summary or its $x$-bins, so 'per-column' and 'the stored response power-sums' assume the package's data model.

**Fix.** Name the summary in the lead: '... of a $y$-quadrature response (the expectation of the functional within each $x$-bin of the joint histogram, read off its stored per-cell sums of the response and its square)'.

- `docs/src/api/ConditionalExpectations.md`:10  (ConditionalExpectations)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-087  The signature comment of conditional_pushforward_density puts LaTeX subscript markup inside a code block. It renders literally with its braces and underscore, so the comment does not read as prose or as code.

Rule: `RULE.TEX.CODE_SPAN_JULIA`.  
Reported by `withheld@api-ConditionalExpectations`, decided by an agent.

'# ≡ f_{g(Y)|X}' shows as raw TeX. The sibling signature comments ('# ≡ m(x)', '# ≡ σ²(x)', '# ≡ 𝔼[g(Y) | X = x]') use plain Unicode, which reads cleanly.

**Fix.** Write the comment in plain Unicode like its siblings, e.g. '# ≡ density of g(Y) given X = x', or drop it.

!!! warning "The code may be what is wrong"
    The reader believed the documentation and the code disagree, and that the code is the wrong half. That is a question for the owner, not a documentation fix.

- `src/ConditionalExpectations.jl`:508  (conditional_pushforward_density)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-088  The symbol order (set as \texttt{order}) bounds the quadrature sum, the exactness degree, the bias order and the admissible central-moment orders, but nothing on the page says what it is (the summary's highest stored power-sum order). A reader cannot tell which polynomials are exact or which p central_moment accepts.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-Expectations`, decided by an agent.

Every docstring and both page introductions use \texttt{order}, and none binds it. The expectation docstring binds the per-bin count n_j, the center c_j and the power-sums S_p^{(j)}, but not the summation limit.

Also reported by `withheld@api-Expectations`: The total sample count n is unbound in the formulas of expectation and pushforward_density; the page binds it only later, in expectation_standard_error. In pushforward_density, the samples X_i and the scaled kernel K_h are also never bound, so the remark about 'the surviving 1/h' points at a 1/h that nothing has introduced. Evidence: The phrase 'with $n$ the sample count' first appears in the expectation_standard_error docstring, after both formulas. Nothing on the page defines K_h(u) = K(u/h)/h. Fix proposed: In each docstring, bind n as the total sample count (n = \sum_j n_j), X_i as the samples, and K_h(u) = K(u/h)/h at its first use.

**Fix.** Bind it at first use in each docstring, for example 'up to the summary's power-sum order \texttt{order} (the highest p stored per bin)'.

- `src/Expectations.jl`:187  (expectation)
- `src/Expectations.jl`:191  (expectation)
- `src/Expectations.jl`:153  (central_moment)
- `src/Expectations.jl`:301  (pushforward_density)
- `src/Expectations.jl`:308  (pushforward_density)
- `docs/ARCHITECTURE.md`:330  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.1 Marginal density ``f_X``)
- `docs/ARCHITECTURE.md`:576  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- `docs/ARCHITECTURE.md`:729  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:730  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/api/Kernels.md`:4  (Kernels)
- `docs/src/index.md`:85  (GBGStreamingStatistics.jl / Main components)
- `docs/src/notation.md`:46  (Notation / Smoothing kernels)
- `src/Bandwidth.jl`:538  (bandwidth)
- `src/ConditionalDensity.jl`:95  (conditional_density)
- `src/ConditionalExpectations.jl`:112  (conditional_expectation)
- `src/ConditionalExpectations.jl`:303  (conditional_mean)
- `src/ConditionalExpectations.jl`:387  (conditional_variance)
- `src/ConditionalExpectations.jl`:528  (conditional_pushforward_density)
- `src/ConditionalExpectations.jl`:547  (conditional_pushforward_density)
- `src/Expectations.jl`:319  (pushforward_density)
- `src/Kernels.jl`:19  (gaussian_kernel)
- `src/Kernels.jl`:30  (epanechnikov_kernel)
- `src/MarginalDensity.jl`:428  (density)
- `src/MarginalDensity.jl`:441  (density)
- `src/MarginalDensity.jl`:704  (density_standard_error)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-089  The page introduction reads a functional 'of $X$' off 'the stored power-sums' without saying what X is or which object stores the power-sums (a Histogram1D summary of samples of X). It also speaks of 'offsets' without saying what they are offsets from (the bin center).

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-Expectations`, decided by an agent.

The bin center as the reference point, and Histogram1D as the summary, appear only later, inside the expectation docstring.

**Fix.** Name the summary and the reference point in the introduction, for example 'off the power-sums a Histogram1D of samples of $X$ stores (per bin, the summed powers of the samples' offsets from the bin center)'.

- `docs/src/api/Expectations.md`:3  (Expectations)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-090  Several package terms appear before the page defines them, or are never defined: 'the recovered moments' and 'the composed kernel' in the introduction (the composed kernel K_h(g(x) - z) is explained only in pushforward_density); 'the shift formula', named bare in the mean and var docstrings and glossed only in passing inside expectation; 'the bandwidth rules', defined nowhere on the page; and 'readouts' / 'reader' used as nouns for functions.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Expectations`, decided by an agent.

On a cold read of this page alone, 'recovered moments', 'composed kernel', 'shift formula' and 'bandwidth rules' are nouns that assume the package's context. The shift formula is glossed only inside a dash-chained parenthetical in another docstring, and nothing on this page defines or links the bandwidth rules.

**Fix.** Define each term in a clause where it first appears, or @ref-link it: 'the composed kernel x -> K_h(g(x) - z)'; 'the shift formula (the binomial re-expansion of each bin's powers from its center to the sample mean)'; link bandwidth to its docstring.

- `docs/src/api/Expectations.md`:8  (Expectations)
- `docs/src/api/Expectations.md`:10  (Expectations)
- `src/Expectations.jl`:100  (mean)
- `src/Expectations.jl`:114  (var)
- `docs/src/api/Expectations.md`:5  (Expectations)
- `docs/ARCHITECTURE.md`:414  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:737  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:743  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:790  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:802  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/notation.md`:61  (Notation / Smoothing kernels)
- `docs/src/notation.md`:106  (Notation / The moments of the data)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-091  Most of the 'Estimand.' and 'Estimator.' run-in paragraphs are noun phrases with no verb. The mean docstring shows that the same content can be written as a full sentence ('The estimand is the population mean').

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-Expectations`, decided by an agent.

Verbless openers appear in expectation, expectation_standard_error, pushforward_density, pushforward_density_standard_error, var, std and central_moment. In mean, 'Estimand. The estimand is the population mean $\mathbb{E}[X]$.' is a full sentence.

**Fix.** Give each one a verb, as mean does: 'The estimand is the population expectation ...', 'The estimator is the recovered variance ...'.

- `src/Expectations.jl`:177  (expectation)
- `src/Expectations.jl`:293  (pushforward_density)
- `src/Expectations.jl`:111  (var)
- `src/Expectations.jl`:113  (var)
- `src/Expectations.jl`:136  (std)
- `docs/ARCHITECTURE.md`:455  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.B Power-sum quadrature readouts)
- `docs/src/notation.md`:16  (Notation / The samples)
- `src/Expectations.jl`:5
- `src/Expectations.jl`:316  (pushforward_density)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-092  Some caveat sentences lack a subject or a verb.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-Expectations`, decided by an agent.

'`NaN` for an empty summary.' (expectation_standard_error) has no verb. 'Needs $g$ smooth ...' (pushforward_density) has no subject. The expectation docstring writes the same ideas in full: 'An empty summary gives `NaN`.' and 'The estimand needs $g$ smooth ...'.

**Fix.** Write 'An empty summary gives `NaN`.' and 'The estimator needs $g$ smooth on the bin scale ...'.

- `src/Expectations.jl`:244  (expectation_standard_error)
- `src/Expectations.jl`:309  (pushforward_density)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-093  In expectation, a parenthetical chains three em-dash clauses, so a reader cannot tell what 'the $g = x^p$ case' or 'the binomial re-expansion' is in apposition to. The parenthetical does not read as a sentence aloud.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-Expectations`, decided by an agent.

The full parenthetical reads '(the remainder terminates — this is the shift formula behind ... — the binomial re-expansion ... — the $g = x^p$ case)'. The dashes nest three appositions with no punctuation to separate their scopes.

**Fix.** Split it, for example: '(the Taylor remainder terminates). For $g = x^p$ this is the shift formula behind `mean`, `var` and `central_moment`, the binomial re-expansion from the bin center to the sample mean.'

- `src/Expectations.jl`:192  (expectation)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-094  Two signature code blocks carry mathematics in a trailing comment. In pushforward_density it is raw LaTeX source (f_{g(X)}), which shows literally with its braces instead of rendering as a formula.

Rule: `RULE.TEX.CODE_SPAN_JULIA`.  
Reported by `withheld@api-Expectations`, decided by an agent.

Both comments sit inside fenced signature blocks, where '# ≡ f_{g(X)}' displays as literal TeX. Each docstring's prose already binds the symbol in rendered math.

**Fix.** Drop the comments from the signatures; each estimand's prose paragraph already binds the symbol.

- `src/Expectations.jl`:175  (expectation)
- `src/Expectations.jl`:290  (pushforward_density)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-095  Verbless clauses stand as sentences: the Histogram1D/Histogram2D entry's second clause has no verb ('the ... variant also rescalable ...'), and the list lead-in is a bare noun phrase with no main verb.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-GBGStreamingStatistics`, decided by an agent.

'the count_type = Float64 variant also rescalable for post-hoc mixtures' lacks 'is'; 'Its operation families, each read post hoc from the one summary:' has no finite verb, unlike a lead-in such as 'It has these operation families, each read post hoc from the one summary:'.

**Fix.** Write 'the `count_type = Float64` variant is also rescalable, for ...' and 'It offers these operation families, each read post hoc from the one summary:'.

- `src/GBGStreamingStatistics.jl`:24  (GBGStreamingStatistics)
- `src/GBGStreamingStatistics.jl`:20  (GBGStreamingStatistics)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-096  The package's own terms for the rescaling operations are used without definition: 'hyperstatistics' (not a standard term in the field; the standard neighbour is superstatistics), 'post-hoc mixtures' (a mixture of what, weighted how) and bare 'forgetting' (exponential forgetting of older samples is not stated).

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-GBGStreamingStatistics`, decided by an agent.

On a cold read the parenthetical '(hyperstatistics)' names a concept nowhere explained on the page, and 'forgetting' does not say what is forgotten or how rescaling achieves it.

**Fix.** Define in a clause, e.g. 'rescalable, so summaries can be reweighted and added into a mixture whose weights are chosen after accumulation (hyperstatistics), or older data down-weighted geometrically (exponential forgetting)'; or link the definitions.

- `src/GBGStreamingStatistics.jl`:24  (GBGStreamingStatistics)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-097  'the stored power sums' presupposes what the summary stores, which the page never says: the summary is only called 'a fixed-size summary' and named Histogram1D/Histogram2D, never said to carry per-bin counts and power sums; 'Quadrature readouts' is likewise package jargon unexplained here.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-GBGStreamingStatistics`, decided by an agent.

The definite article 'the stored power sums' refers to an object the page has not introduced; a reader cannot tell power sums of what, or about which point, nor what a 'readout' is as opposed to the other families.

**Fix.** Say once, where the summaries are introduced, that each bin carries its count and the power sums of its samples, and gloss the family as 'expectations computed by quadrature over those power sums'.

- `src/GBGStreamingStatistics.jl`:33  (GBGStreamingStatistics)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-098  The code identifier count_type is used without being bound: the page does not say it is a constructor option selecting the element type of the bin counts, nor why the default variant cannot be rescaled.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-GBGStreamingStatistics`, decided by an agent.

A reader meets 'the `count_type = Float64` variant' with no prior mention of count_type or of a default integer count type, so the contrast the clause draws is opaque.

**Fix.** Bind it in apposition: 'the variant constructed with floating-point counts (`count_type = Float64`; the default integer counts cannot hold a fractional weight)'.

- `src/GBGStreamingStatistics.jl`:24  (GBGStreamingStatistics)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-099  Framing terms assume the package's context: 'backend-agnostic' (which backends: CPU arrays, GPUs?) and 'the underlying state' (the page speaks only of samples; a statistician cannot tell whether 'state' means the raw samples or a modelled system's state).

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-GBGStreamingStatistics`, decided by an agent.

'backend-agnostic' opens the docstring with no clause saying what a backend is here; 'so the underlying state is never stored or re-read' introduces 'state' with a definite article and no referent on the page.

**Fix.** Gloss 'backend-agnostic (the same code runs on CPU arrays and on a GPU)' and replace 'the underlying state' with 'the samples' or define it.

- `src/GBGStreamingStatistics.jl`:12  (GBGStreamingStatistics)
- `src/GBGStreamingStatistics.jl`:17  (GBGStreamingStatistics)
- `docs/ARCHITECTURE.md`:30  (GBGStreamingStatistics — Architecture and Rationale / 1. Purpose and scope)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-100  The symbol T in the identity sentence of Histogram1D is bound nowhere on the page, and 'the canonical summary' is a term of the package that the page never defines; n (the sample count) is also unbound at this first use.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-Histograms`, decided by an agent.

The reader meets $T(n\mathbb{P}_n)$ as a map applied to the unnormalized empirical measure, but nothing says what T is (the binning map? the summary construction?) or why it is 'canonical'; the rest of the page never uses T again.

**Fix.** Bind T in plain words at first use (e.g. 'the binning map T applied to the unnormalized empirical measure nP_n of the n samples'), or drop the formula and keep the prose.

- `src/Histograms.jl`:21  (Histogram1D)
- `docs/ARCHITECTURE.md`:73  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The spine: empirical measure → summary → post-hoc-smoothed estimator)
- `docs/ARCHITECTURE.md`:87  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The spine: empirical measure → summary → post-hoc-smoothed estimator)
- `docs/src/notation.md`:36  (Notation / The summary)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-101  The page's opening sentence has no verb: it is a colon-headed noun phrase ('The summary types: X and Y — ... — plus Z for two-time statistics.').

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-Histograms`, decided by an agent.

Read aloud, it names things without saying anything about them; the only verb ('fold') belongs to a relative clause inside the dash aside.

**Fix.** Give it a verb, e.g. 'This page documents the summary types X and Y, ..., and observe_pairs!, which folds two-time pairs into a Histogram2D.'

- `docs/src/api/Histograms.md`:3  (Histograms)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-102  Mathematical relations and products are set as code spans instead of mathematics, inconsistently with the page's own math-set shapes ($\texttt{order} \times \texttt{nbins}$, $5 \times (\texttt{n\_x} \cdot \texttt{n\_y})$).

Rule: `RULE.TEX.CODE_SPAN_JULIA`.  
Reported by `withheld@api-Histograms`, decided by an agent.

The range constraint is an inequality, and the counts length is a product; both are statements about quantities, not code to be typed. The powersums shape on the same page is set as math, so the counts length in backticks reads as a different kind of object.

**Fix.** Set them as math with the names in \texttt: $\texttt{upper} > \texttt{lower}$, length $\texttt{n\_x} \cdot \texttt{n\_y}$.

!!! warning "The code may be what is wrong"
    The reader believed the documentation and the code disagree, and that the code is the wrong half. That is a question for the owner, not a documentation fix.

- `src/Histograms.jl`:47  (Histogram1D)
- `src/Histograms.jl`:253  (Histogram2D)
- `src/Histograms.jl`:75  (Histogram1D)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-103  'cell counts' is used for the numbers of cells n_x, n_y, while the same page uses 'counts' (the `counts` field, 'per-cell counts') for the number of samples in each cell; a cold reader cannot tell which is meant.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Histograms`, decided by an agent.

In the Fields list, n_x and n_y are called 'the box and cell counts', two lines above 'counts — per-cell counts'; the merge! requirement 'matching box, cell counts, and lag' would read as requiring equal sample counts per cell, which is wrong.

**Fix.** Say 'the numbers of cells' (or 'the cell grid') for n_x, n_y in both places.

- `src/Histograms.jl`:250  (Histogram2D)
- `src/Histograms.jl`:234  (Histogram2D)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-104  The bandwidth is never bound as a symbol. The asymmetric-kernel docstrings write it as $b$ in every formula but name the keyword `h` ('pass `h` to override'), and never say that b is that bandwidth. The Epanechnikov docstring uses a bare $h$. A cold reader has to guess that $b$ and $h$ are the same quantity.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-Kernels`, decided by an agent.

gamma_kernel uses $\mathrm{Gamma}(\alpha(x), b)$, $x/b$ and $2b$, then says 'pass `h` to override'. beta_kernel uses $(x/b, (1-x)/b)$ and $2b$. epanechnikov_kernel says 'within $h$ of the query point'. Nowhere on the page is b or h said to be the bandwidth.

**Fix.** Bind each symbol at its first use in apposition, e.g. 'the bandwidth $b$ (the keyword `h`)' in gamma_kernel and beta_kernel, and 'the bandwidth $h$' in epanechnikov_kernel.

- `src/Kernels.jl`:140  (gamma_kernel)
- `src/Kernels.jl`:175  (beta_kernel)
- `src/Kernels.jl`:32  (epanechnikov_kernel)
- `src/Kernels.jl`:160  (Kernels)
- `src/Kernels.jl`:176  (beta_kernel)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-105  Several symbols appear before anything on the page binds them. The intro uses $R(K)$ and $\mu_2(K)$ (they are glossed only later, in gaussian_kernel). gamma_kernel's estimator formula uses $n$, $n_j$ and $c_j$, which are bound only in the later AsymmetricKernel docstring. $f_X'$ appears with $f_X$ never introduced. $\hat\sigma_X$ is never bound in either docstring that uses it.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-Kernels`, decided by an agent.

The intro says 'computes $R(K)$ and $\mu_2(K)$ by quadrature' before the gaussian_kernel gloss. gamma_kernel gives '$\hat f_X(x) = \tfrac{1}{n} \sum_j n_j \mathrm{Gamma}(c_j; ...)$' with no gloss of n, n_j or c_j. gamma_kernel and AsymmetricKernel both write '$\hat\sigma_X n^{-2/5}$' and never say it is the sample standard deviation.

**Fix.** Gloss each symbol in apposition at its first use, e.g. 'the roughness $R(K)$ and second moment $\mu_2(K)$', 'the bin counts $n_j$ at centers $c_j$, $n$ the total count', 'the density $f_X$', and 'the sample standard deviation $\hat\sigma_X$'.

- `docs/src/api/Kernels.md`:4  (Kernels)
- `src/Kernels.jl`:145  (gamma_kernel)
- `src/Kernels.jl`:142  (gamma_kernel)
- `src/Kernels.jl`:148  (gamma_kernel)
- `src/Kernels.jl`:89  (AsymmetricKernel)
- `docs/ARCHITECTURE.md`:415  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:423  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:561  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- `docs/ARCHITECTURE.md`:576  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- `docs/ARCHITECTURE.md`:740  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:743  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/notation.md`:60  (Notation / Smoothing kernels)
- `src/Bandwidth.jl`:336  (bandwidth)
- `src/Bandwidth.jl`:337  (bandwidth)
- `src/Bandwidth.jl`:341  (bandwidth)
- `src/Bandwidth.jl`:390  (bandwidth)
- `src/Bandwidth.jl`:402  (bandwidth)
- `src/Bandwidth.jl`:511  (bandwidth)
- `src/Bandwidth.jl`:516  (bandwidth)
- `src/Bandwidth.jl`:523  (bandwidth)
- `src/ConditionalExpectations.jl`:282  (conditional_mean)
- `src/Kernels.jl`:20  (gaussian_kernel)
- `src/Kernels.jl`:30  (epanechnikov_kernel)
- `docs/ARCHITECTURE.md`:790  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-106  The phrase 'the $\sqrt{\cdot}$ correction' in beta_kernel uses a definite article, as if the reader already knew it, but the page never defines it. The clause gives no formula, and 'a quadratic's root in the bandwidth and position' does not say which quadratic.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Kernels`, decided by an agent.

beta_kernel says 'the $\sqrt{ \cdot }$ correction within $2b$ of an edge replaces that shape by a quadratic's root in the bandwidth and position'. gamma_kernel gives its boundary shape explicitly, as $(x/b)^2/4 + 1$, but beta_kernel gives none.

**Fix.** State the boundary shape explicitly, as gamma_kernel does, or describe it without the definite article and link to Chen (1999) for the formula.

- `src/Kernels.jl`:176  (beta_kernel)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-107  'AMISE local-scale weight' is the package's own term, and its first use comes in the Fields list, which renders before the Arguments bullet that defines it (through the leading bias and the variance). On a top-to-bottom read the term arrives undefined.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Kernels`, decided by an agent.

The Fields bullet reads 'the AMISE local-scale weight $w(x)$, a callable, or `nothing` when unset.' The meaning ('the leading bias is $\tfrac{1}{2} b w(x) f_X''(x)$ and the variance $\propto w(x)^{-1/2}$') comes only in the later Arguments bullet.

**Fix.** Define the weight where it is first used, or shorten the Fields bullet to a pointer to the Arguments definition.

- `src/Kernels.jl`:73  (AsymmetricKernel)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-108  gamma_kernel calls itself 'its own degree-0' estimator and says the `degree` and `adaptive` options do not apply. Nothing on this page says what 'degree' (the local-polynomial order) or 'adaptive' means. They are the package's option names, and the link to density is on an earlier sentence, not on this clause.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-Kernels`, decided by an agent.

'The `degree`, `adaptive` and `support` options of `density` do not apply — the kernel is its own degree-0, boundary-free estimator.' Only `support` is glossed earlier ('boundary renormalization').

**Fix.** Gloss the options in a clause (e.g. 'the local-polynomial `degree` and variable-bandwidth `adaptive` options'), or put the density link on this sentence.

- `src/Kernels.jl`:152  (gamma_kernel)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-109  'Left unset' is a dangling participle. The sentence's subject is `rule = :plugin`, but what is left unset is `local_scale_weight`, so read literally the sentence says the rule is left unset.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-Kernels`, decided by an agent.

'Left unset (`nothing`, the default), `rule = :plugin` errors rather than use a wrong weight'

**Fix.** Give the participle its subject, e.g. 'With the weight left unset (`nothing`, the default), `rule = :plugin` errors rather than use a wrong weight'.

- `src/Kernels.jl`:88  (AsymmetricKernel)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-110  Several lead sentences have no finite verb: a labelled noun phrase or a formula stands where a sentence should be.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@api-MarginalDensity`, decided by an agent.

The density_standard_error Estimator paragraph is a formula followed by 'from the recovered ...' with no verb. 'Restricted to degree 0/1 ...' has no subject or verb. density_with_standard_error opens with 'The density estimate and its pointwise ..., from a single option specification', a noun phrase. density's Estimand is also a verbless noun phrase ('The marginal (Lebesgue) density ... of the variate X with law P_X, defined by ...').

**Fix.** Give each one a verb, e.g. 'The estimator is sqrt(...), computed from ...'; 'The band is restricted to degree 0/1 ...'; 'Returns the density estimate and its pointwise standard error from one option specification'; 'The estimand is the marginal (Lebesgue) density f_X ...'.

- `src/MarginalDensity.jl`:702  (density_standard_error)
- `src/MarginalDensity.jl`:720  (density_standard_error)
- `src/MarginalDensity.jl`:779  (density_with_standard_error)
- `src/MarginalDensity.jl`:416  (density)
- `docs/ARCHITECTURE.md`:330  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.1 Marginal density ``f_X``)
- `docs/ARCHITECTURE.md`:410  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `src/Bandwidth.jl`:317  (bandwidth)
- `src/MarginalDensity.jl`:419  (density)
- `src/MarginalDensity.jl`:696  (density_standard_error)
- `src/MarginalDensity.jl`:699  (density_standard_error)
- `src/MarginalDensity.jl`:706  (density_standard_error)
- `docs/ARCHITECTURE.md`:327  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.1 Marginal density ``f_X``)
- `docs/src/notation.md`:35  (Notation / The summary)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-111  In the density docstring, some symbols are used before they are bound, and some are never bound.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@api-MarginalDensity`, decided by an agent.

h_j in the degree-0 formula is bound only later, in the adaptive bullet's parenthetical ('sample-point h_j per bin'), and the fixed-h case h_j = h is never stated. beta_0 in e^{beta_0} is never bound as the intercept of the local fit. sigma_X and n first appear in the h bullet's Silverman formula; sigma_X is bound only later, in the start bullet. c_j is used in the :linear formula but bound as the bin center only in the following :simple clause. f_theta, with an unhatted theta, appears in ParametricStart but is never bound.

Also reported by `withheld@api-MarginalDensity`: The same estimate goes by two symbols, f_h(x) with a hat and f_X(x) with a hat, and the page never says they are one. Evidence: The Estimator paragraph introduces '$\hat f_h(x)$'. The degree, support and start bullets then write '$\hat f_X(x)$'. The standard-error docstring returns to $\hat f_h(x)$, so a cold reader cannot tell whether these are different quantities. Fix proposed: Use one symbol throughout, or bind the second as the same estimate.

**Fix.** Bind each symbol where it first appears: h_j (equal to h when adaptive = false), beta_0 the fitted intercept, n the total count and sigma_X the recovered standard deviation, c_j the bin center, and f_theta (or use f_{theta-hat} consistently).

- `src/MarginalDensity.jl`:440  (density)
- `src/MarginalDensity.jl`:443  (density)
- `src/MarginalDensity.jl`:451  (density)
- `src/MarginalDensity.jl`:488  (density)
- `src/MarginalDensity.jl`:38  (ParametricStart)
- `src/MarginalDensity.jl`:419  (density)
- `docs/ARCHITECTURE.md`:737  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:331  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.1 Marginal density ``f_X``)
- `docs/ARCHITECTURE.md`:831  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:363  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.1 Marginal density ``f_X``)
- `docs/ARCHITECTURE.md`:776  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:780  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:330  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.1 Marginal density ``f_X``)
- `docs/ARCHITECTURE.md`:410  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `src/Bandwidth.jl`:317  (bandwidth)
- `src/MarginalDensity.jl`:696  (density_standard_error)
- `src/MarginalDensity.jl`:699  (density_standard_error)
- `src/MarginalDensity.jl`:702  (density_standard_error)
- `src/MarginalDensity.jl`:706  (density_standard_error)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-112  The package's term 'recovered' is used throughout and never defined on the page; in one place it even describes the computed estimate.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-MarginalDensity`, decided by an agent.

'the recovered $\bar X$, $\hat\sigma_X^2$', 'the recovered sample mean and variance' and 'the recovered moments' all assume the reader knows that moments are read back off the summary's stored power-sums. 'from the recovered $\hat f_h(x)$' applies the word to the density estimate itself, which nothing recovers.

**Fix.** Define it once, e.g. 'recovered, that is, computed exactly from the summary's stored counts and power-sums', and write 'the estimate' for f_h(x).

- `src/MarginalDensity.jl`:474  (density)
- `src/MarginalDensity.jl`:702  (density_standard_error)
- `docs/ARCHITECTURE.md`:330  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.1 Marginal density ``f_X``)
- `docs/ARCHITECTURE.md`:410  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `src/Bandwidth.jl`:317  (bandwidth)
- `src/MarginalDensity.jl`:419  (density)
- `src/MarginalDensity.jl`:696  (density_standard_error)
- `src/MarginalDensity.jl`:699  (density_standard_error)
- `src/MarginalDensity.jl`:706  (density_standard_error)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-113  'Band' is defined as the interval f_h(x) ± 1.96 SE(x), yet elsewhere the word names the standard error itself; the page also uses it before defining it.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@api-MarginalDensity`, decided by an agent.

The page intro calls density_standard_error 'its pointwise variability band', before any definition. The standard-error docstring then defines 'The pointwise variability band is $\hat f_h(x) \pm 1.96 \operatorname{SE}(x)$'. density_with_standard_error then says '`standard_error` its band', which conflates the half-width scale with the interval.

**Fix.** Call the reader's output the pointwise standard error, and keep 'band' for the ±1.96 SE interval.

- `docs/src/api/MarginalDensity.md`:3  (MarginalDensity)
- `src/MarginalDensity.jl`:782  (density_with_standard_error)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-114  The density signature block carries mathematics in a code comment, '# ≡ f_X', so it shows as a raw underscore and a glyph instead of rendered math.

Rule: `RULE.TEX.CODE_SPAN_JULIA`.  
Reported by `withheld@api-MarginalDensity`, decided by an agent.

On the rendered page the signature block shows '# ≡ f_X' literally. f_X is mathematics, and inside code it cannot render as a subscript; the Estimand paragraph directly below already binds f_X.

**Fix.** Drop the comment from the signature, since the Estimand line binds f_X.

- `src/MarginalDensity.jl`:411  (density)

### 2026-10-19-documentation-review-01-115  Symbols are used several sections before the page binds them, so a top-to-bottom read meets them bare.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@architecture`, decided by an agent.

$S_p$ appears in the lens paragraph of section 2 but is defined only in section 3 (Histogram1D). $Q$ (and 'Pawula') appears in the estimand-family list, but the Pawula quotient $Q$ is glossed only in section 3. $\delta(K)$ is used in section 4.5 and defined only in the section 11 bandwidth entry. $h^\star(x)$ appears in the family list and in section 2's FFT paragraph, but its formula comes only in section 4.5.

**Fix.** Bind each one in plain words at its first use (e.g. 'the stored per-bin power-sums $S_p$ (section 3)', 'the Pawula quotient $Q = M_4/(3M_2^2)$', '$\delta(K) = [R(K)/\mu_2(K)^2]^{1/5}$', 'the local bandwidth $h^\star(x)$'), or move the later definition forward.

- `docs/ARCHITECTURE.md`:129  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The algebra of measures — a structural lens)
- `docs/ARCHITECTURE.md`:159  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / Estimand families (the operation taxonomy))
- `docs/ARCHITECTURE.md`:419  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:162  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / Estimand families (the operation taxonomy))
- `docs/ARCHITECTURE.md`:740  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:749  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:755  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:770  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/notation.md`:60  (Notation / Smoothing kernels)
- `docs/src/notation.md`:61  (Notation / Smoothing kernels)
- `src/Bandwidth.jl`:335  (bandwidth)
- `src/Bandwidth.jl`:336  (bandwidth)
- `src/Bandwidth.jl`:350  (bandwidth)
- `src/Bandwidth.jl`:360  (bandwidth)
- `src/Bandwidth.jl`:363  (bandwidth)
- `src/Bandwidth.jl`:402  (bandwidth)
- `docs/ARCHITECTURE.md`:182  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The histogram is the summary; the ECF is deferred (likely never needed))
- `docs/ARCHITECTURE.md`:423  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:429  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:865  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:873  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:918  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/api/Bandwidth.md`:8  (Bandwidth)
- `docs/src/index.md`:39  (GBGStreamingStatistics.jl)
- `docs/src/index.md`:115  (GBGStreamingStatistics.jl / Main components)
- `docs/src/notation.md`:52  (Notation / Smoothing kernels)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-116  Some symbols are never bound anywhere on the page.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@architecture`, decided by an agent.

$\lambda$ in the hyperstatistics entry has no definition, and the surrounding code calls the factor `weight`. $\kappa$ in the Kramers–Moyal entry is the OU mean-reversion rate, but the page never introduces the OU parametrization. $\ell(x)$ is only implied by 'leverage-based'. $C(K,d)$ is never defined, and its relation to the 'canonical factor' that follows is not stated.

**Fix.** Gloss each one in apposition: '$\lambda$ the scale factor (`weight`)', 'the OU rate $\kappa$', 'the local-regression leverage $\ell(x)$', and '$C(K,d)$, the kernel- and degree-dependent constant'.

- `docs/ARCHITECTURE.md`:811  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:720  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:389  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.2 Conditional expectation ``m(x)`` and variance ``\sigma^2(x)``)
- `docs/ARCHITECTURE.md`:413  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:812  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:414  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/src/notation.md`:61  (Notation / Smoothing kernels)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-117  The page refers to a section 4.4 that it does not contain. Its sections run 4.A, 4.1–4.3, 4.5, 4.B and 4.6, so the reader cannot follow these references.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@architecture`, decided by an agent.

There is no 4.4 heading anywhere on the page. The archived increment-moment material survives only in the section 11 log entries.

**Fix.** Point these references at the section 11 archival entry, or add a one-line 4.4 stub saying the section is archived.

- `docs/ARCHITECTURE.md`:159  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / Estimand families (the operation taxonomy))
- `docs/ARCHITECTURE.md`:317  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half))
- `docs/ARCHITECTURE.md`:568  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-118  Process vocabulary from inside the project appears without definition. A visitor cannot tell who 'the user' or 'the owner' is, or what the 'later sessions', 'Stages 1–3' and 'the guards' are.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@architecture`, decided by an agent.

'The user' here means the package's owner, not the reader who uses the package, and on a public page it reads as the reader. 'Later sessions' and 'Stages 1–3' refer to the development process, and 'the guards' to argument-validation code. None of them is introduced.

**Fix.** Say 'the package's author' (or drop the attribution). Replace 'later sessions' with 'later decisions'. Drop 'Stages 1–3'. Write 'as the argument checks decide them'.

- `docs/ARCHITECTURE.md`:172  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The histogram is the summary; the ECF is deferred (likely never needed))
- `docs/ARCHITECTURE.md`:629  (GBGStreamingStatistics — Architecture and Rationale / 9. Naming the package)
- `docs/ARCHITECTURE.md`:954  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:682  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:866  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:348  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.1 Marginal density ``f_X``)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-119  Internals of the sibling package GBGPeriodicFields (PF) are used as explanations without being defined. A reader new to the package has no access to PF's operator split, its 'spectral product-as-*', its private functions or its numbered invariants.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@architecture`, decided by an agent.

Each of these is used to justify a design choice, so the justification is lost to a reader who has not read PF.

**Fix.** Gloss each one in a clause (e.g. 'PF's rule that a spectral product is never spelled `*`'), or cut the analogy where the sentence stands without it.

- `docs/ARCHITECTURE.md`:513  (GBGStreamingStatistics — Architecture and Rationale / 6. Invariants — deliberate absences are statements)
- `docs/ARCHITECTURE.md`:257  (GBGStreamingStatistics — Architecture and Rationale / 3. The accumulated summary — the universal lower half)
- `docs/ARCHITECTURE.md`:478  (GBGStreamingStatistics — Architecture and Rationale / 5. Backend-agnosticism and the hot-loop contract)
- `docs/ARCHITECTURE.md`:94  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The spine: empirical measure → summary → post-hoc-smoothed estimator)
- `docs/ARCHITECTURE.md`:667  (GBGStreamingStatistics — Architecture and Rationale / 10. Open questions)
- pass 1: verifier `LINE_ONLY`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-120  Terms coined by the package are used before the page defines them.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@architecture`, decided by an agent.

'Hyperstatistics' and 'forgetting' first appear in the section 6 invariant, but are explained only in section 10 ('superstatistics; Beck & Cohen 2003') and the section 11 entry. 'Readers' (the conditioning readers) is never defined; the page defines 'readout' instead.

**Fix.** Gloss each one at first use: 'hyperstatistics (pooling summaries under a mixing law over a parameter)', 'exponential forgetting (down-weighting older data)'. Write 'readouts' rather than 'readers'.

- `docs/ARCHITECTURE.md`:521  (GBGStreamingStatistics — Architecture and Rationale / 6. Invariants — deliberate absences are statements)
- `docs/ARCHITECTURE.md`:290  (GBGStreamingStatistics — Architecture and Rationale / 3. The accumulated summary — the universal lower half)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-121  Many prose sentences lack a subject or a verb. These are note-form fragments in the body and in the section 11 decision log.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@architecture`, decided by an agent.

Each quoted string is a whole sentence as it renders, with no finite verb or no subject. Further instances follow the same pattern ('Built.', 'Tested: ...', 'Enforced by an `@allocated` audit ...', 'Done, each its own commit:').

**Fix.** Give each sentence a subject and a verb (e.g. 'This mirrors PF's physical-meaning-first rule.', 'The log is append-only and mirrors PF's.', 'The name is operationally precise ...').

- `docs/ARCHITECTURE.md`:202  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / Documentation policy — estimand first, estimator second)
- `docs/ARCHITECTURE.md`:544  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- `docs/ARCHITECTURE.md`:681  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:629  (GBGStreamingStatistics — Architecture and Rationale / 9. Naming the package)
- `docs/ARCHITECTURE.md`:843  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:258  (GBGStreamingStatistics — Architecture and Rationale / 3. The accumulated summary — the universal lower half)
- `docs/ARCHITECTURE.md`:696  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:736  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:462  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.B Power-sum quadrature readouts)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-122  Source line wraps break sentences in the rendered page. In the Kernels log entry, the parenthesis opened at '(`gamma_kernel`, ...; Chen 1999,' never closes, and the verb 'were moved' starts a separate paragraph. In section 5, 'accumulate-then- merge' carries a stray space.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@architecture`, decided by an agent.

A wrapped source line that starts with '2000)' is parsed by Markdown as an ordered-list marker. That swallows 'Chen 1999, 2000)' and splits the subject from its verb across two blocks. The section 5 hyphen-wrap leaves 'then- merge'.

**Fix.** Rewrap so that no line begins with '2000)' (or escape it as '2000\)'), and join 'accumulate-then-merge'.

- `docs/ARCHITECTURE.md`:880  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:881  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:490  (GBGStreamingStatistics — Architecture and Rationale / 5. Backend-agnosticism and the hot-loop contract)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-123  A Julia symbol is set as mathematics, which is the converse of the code/math boundary this rule draws. `:sj` is a literal `rule` value, but it is written in dollar signs, so it renders as math (colon, italic s j) and not as code.

Rule: `RULE.TEX.CODE_SPAN_JULIA`.  
Reported by `withheld@architecture`, decided by an agent.

Everywhere else the page writes the option as `:sj` in a code span, e.g. 'the data-driven Sheather–Jones `:sj`'.

**Fix.** Write it as the code span `:sj`.

- `docs/ARCHITECTURE.md`:795  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)

### 2026-10-19-documentation-review-01-124  Symbols appear on the page that no sentence binds in plain language: the starred bandwidth h^star(x) (first met in the Status paragraph, again under Bandwidth), the subscript tau of the transition density, the X_i and n of the empirical-measure formula, and the exponent p of the g = x^p case. A first-time reader must guess that tau is the lag, n the sample count, and what makes h^star(x) different from h.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@home`, decided by an agent.

Status names 'the response-aware local bandwidth $h^\star(x)$' with no gloss; ConditionalDensity writes '$p_\tau(y\mid x)$' and never says what tau is; Accumulators writes '#{i : X_i in A}/n' without naming X_i or n; Expectations ends 'the $g = x^p$ case' with p unbound.

**Fix.** Bind each at first use in a clause: 'h^star(x), the bandwidth that minimizes the local error at x', 'p_tau(y|x), the density of the value a lag tau later given x', 'over n samples X_1, ..., X_n', 'the moments of order p (g = x^p)'.

- `docs/src/index.md`:38  (GBGStreamingStatistics.jl)
- `docs/src/index.md`:110  (GBGStreamingStatistics.jl / Main components)
- `docs/src/index.md`:78  (GBGStreamingStatistics.jl / Main components)
- `docs/src/index.md`:94  (GBGStreamingStatistics.jl / Main components)
- `docs/ARCHITECTURE.md`:162  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / Estimand families (the operation taxonomy))
- `docs/ARCHITECTURE.md`:182  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The histogram is the summary; the ECF is deferred (likely never needed))
- `docs/ARCHITECTURE.md`:423  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:429  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:865  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:873  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:918  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/api/Bandwidth.md`:8  (Bandwidth)
- `docs/src/index.md`:39  (GBGStreamingStatistics.jl)
- `docs/src/index.md`:115  (GBGStreamingStatistics.jl / Main components)
- `docs/src/notation.md`:52  (Notation / Smoothing kernels)
- `docs/ARCHITECTURE.md`:395  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.3 Conditional and transition density)
- `docs/ARCHITECTURE.md`:551  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- `docs/ARCHITECTURE.md`:568  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- `docs/ARCHITECTURE.md`:734  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/api/ConditionalDensity.md`:3  (ConditionalDensity)
- `docs/src/index.md`:111  (GBGStreamingStatistics.jl / Main components)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-125  The Status paragraph, and once the overview before it, use the package's own nouns that this page defines only later or never: 'Float-count cone variant' (what a Float-count is, and why the cone needs it), 'forgetting' (exponential forgetting is never named), 'readouts' (explained only in the later Expectations bullet), 'response-aware' (never explained), 'asymmetric bounded-support kernels' (introduced only in the later Kernels bullet), 'AMD-ROCm parity' (parity of the GPU result with the CPU result is never said), and 'the underlying state' (the data are called samples everywhere else; 'state' reads as a physical system's state).

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@home`, decided by an agent.

A cold reader meets 'with a Float-count cone variant, re-weightable post hoc for hyperstatistics / forgetting' with no prior definition of either term, and 'the readouts by quadrature' before the Expectations bullet says what a readout is.

**Fix.** Gloss each in a clause where it first appears, or drop it from Status and let the Main components list carry it: 'a variant with floating-point counts, so a summary can be scaled', 'exponential forgetting (down-weighting older samples)', 'GPU results identical to the CPU's (AMD, ROCm)', 'the samples themselves are never stored'.

- `docs/src/index.md`:33  (GBGStreamingStatistics.jl)
- `docs/src/index.md`:40  (GBGStreamingStatistics.jl)
- `docs/src/index.md`:37  (GBGStreamingStatistics.jl)
- `docs/src/index.md`:9  (GBGStreamingStatistics.jl)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-126  Several Main-components bullets are long, period-terminated prose with no finite verb, so they read as sentences missing their verb rather than as short labels: MarginalDensity ('the marginal density ... — one call over the binned KDE, ... — plus ...'), ConditionalDensity, Kernels and Bandwidth. The Accumulators, Histograms and ConditionalExpectations bullets show the same list can carry verbs ('Each is', 'It also offers', 'All are').

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@home`, decided by an agent.

The MarginalDensity bullet runs over five clauses and two dashes and ends with a period, with no verb in its main clause; 'one call over the binned KDE' also leaves the reader to supply what the call does.

**Fix.** Give each long bullet a verb after its label, e.g. 'MarginalDensity: `density` estimates the marginal density f_X in one call, as the binned KDE, ...'; 'ConditionalDensity: `conditional_density` and `transition_density` estimate ...'; 'Bandwidth: `bandwidth` selects the bandwidth by ...'.

- `docs/src/index.md`:95  (GBGStreamingStatistics.jl / Main components)
- `docs/src/index.md`:109  (GBGStreamingStatistics.jl / Main components)
- `docs/src/index.md`:113  (GBGStreamingStatistics.jl / Main components)
- `docs/src/index.md`:84  (GBGStreamingStatistics.jl / Main components)
- `docs/ARCHITECTURE.md`:115  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The algebra of measures — a structural lens)
- `docs/ARCHITECTURE.md`:157  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / Estimand families (the operation taxonomy))
- `docs/ARCHITECTURE.md`:394  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.3 Conditional and transition density)
- `docs/ARCHITECTURE.md`:399  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.3 Conditional and transition density)
- `docs/ARCHITECTURE.md`:403  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.3 Conditional and transition density)
- `docs/ARCHITECTURE.md`:548  (GBGStreamingStatistics — Architecture and Rationale / 7. Documentation, correctness, and testing conventions)
- `docs/ARCHITECTURE.md`:728  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:730  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/api/ConditionalDensity.md`:3  (ConditionalDensity)
- `docs/src/index.md`:110  (GBGStreamingStatistics.jl / Main components)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-127  The `response` bullet's second clause makes `mass` the subject of 'smooths the counts themselves'. The words saying that `mass` is passed as `response` are missing, so the clause reads as a claim about the `mass` argument and contradicts the bullet before it.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@internals`, decided by an agent.

`mass` is the design weight, not a smoother. The intended reading, that passing the counts as `response` gives a count level, is not on the page.

**Fix.** Supply the missing words, for example: 'passing `mass` itself as `response` smooths the counts, as the local-likelihood density does, giving a count level, not a response mean.'

- `src/LocalPolynomial.jl`:97  (local_fit)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-128  The page's own prose uses package terms and a symbol before anything on the page introduces them. The opening paragraph says 'this primitive' and 'the solve' before any primitive or solve is named. The section prose uses 'jet' and $\hat\beta_0$ ahead of the docstring that defines them, pointing only with 'below'.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@internals`, decided by an agent.

A cold reader meets 'this primitive', 'the solve' and 'the jet's leading coefficient' with no antecedent. 'jet' is first glossed in the docstring as 'the coefficient vector (jet)', and $\hat\beta = (\hat\beta_0, \dots)$ is bound only there.

**Fix.** In the opening paragraph, name the primitive: 'the local-polynomial fit `local_fit`, a kernel-weighted least-squares solve'. In the section prose, gloss the jet before using it: '(the vector of fitted local Taylor coefficients $\hat\beta_0, \dots, \hat\beta_d$)'.

- `docs/src/internals.md`:4  (Internals)
- `docs/src/internals.md`:18  (Internals / The local-polynomial jet — `LocalPolynomial.local_fit`)
- `docs/src/internals.md`:19  (Internals / The local-polynomial jet — `LocalPolynomial.local_fit`)
- `docs/src/notation.md`:72  (Notation / The local fit)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-129  The `locations` bullet relies on 'the stored first power-sum', a term of the package's summary that this page never defines. A reader cannot tell which power-sum is stored or what it is taken about, so cannot tell what a centroid location is.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@internals`, decided by an agent.

The page calls `Histogram1D` / `Histogram2D` only 'the binned record of the stream' and never mentions power-sums anywhere else.

**Fix.** Gloss it in a clause, for example: 'centroids $c_j + S_1^{(j)}/n_j$ recovered from the stored within-bin sum of offsets from the bin center'. Alternatively, link the term to its definition on the Histograms page.

- `src/LocalPolynomial.jl`:94  (local_fit)
- `docs/ARCHITECTURE.md`:330  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.1 Marginal density ``f_X``)
- `docs/ARCHITECTURE.md`:331  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.1 Marginal density ``f_X``)
- `docs/ARCHITECTURE.md`:333  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.1 Marginal density ``f_X``)
- `docs/ARCHITECTURE.md`:829  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:831  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:832  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/notation.md`:28  (Notation / The summary)
- `src/LocalPolynomial.jl`:68  (local_fit)
- `src/LocalPolynomial.jl`:70  (local_fit)
- `src/LocalPolynomial.jl`:76  (local_fit)
- `src/LocalPolynomial.jl`:77  (local_fit)
- `src/MarginalDensity.jl`:57  (ParametricStart)
- `src/MarginalDensity.jl`:420  (density)
- `src/MarginalDensity.jl`:441  (density)
- `src/MarginalDensity.jl`:471  (density)
- `src/MarginalDensity.jl`:488  (density)
- `src/MarginalDensity.jl`:490  (density)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-130  Several symbols in the docstring's defining sentence and formula are first used unbound and are bound only later, in the Arguments list. $K$ and $h$ first appear inside $K_h(v) \propto K(v/h)$, and $\bar y_j$ first appears in 'response totals $n_j \bar y_j$'. The rule asks for a plain-language binding at first use.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@internals`, decided by an agent.

$K$ is bound later as 'a symmetric callable $K$' in the `kernel` bullet. $h$ is bound as 'the smoothing bandwidth' in the `h` bullet, and $\bar y_j$ as 'the bin's mean response' in the `response` bullet. All three come after the display formula that uses them.

**Fix.** Bind each one at first use, for example: 'response totals $n_j \bar y_j$ ($\bar y_j$ the bin's mean response), and the kernel $K$ scaled by the bandwidth $h$, $K_h(v) \propto K(v/h)$'.

- `src/LocalPolynomial.jl`:70  (local_fit)
- `docs/src/notation.md`:34  (Notation / The summary)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-131  The page relies on terms that belong to the package and never defines or links them: 'the shift formula', 'the cone', 'the archived module', the '$y$-quadrature', a 'kernel-weighted $y$-profile', 'a joint summary' (used before anything says what a summary is), and 'The rule' (which rule is meant is never said). A reader new to the package cannot tell what these refer to.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@notation`, decided by an agent.

On a cold read nothing on the page says what is being shifted, which cone is meant or what its scalar acts on, which module was archived, what is quadrature along y, what the y-profile is a profile of, or which rule 'The rule' is. 'Joint summary' appears in the first table, before the page has shown what a summary contains.

**Fix.** Add a defining clause at first use of each term, or link it to its definition. For example: 'the shift formula that recovers central moments from the power-sums about bin centers'; 'the cone of summaries, closed under addition and scaling by a nonnegative scalar'; and either drop 'of the archived module' or say what it refers to.

- `docs/src/notation.md`:5  (Notation)
- `docs/src/notation.md`:38  (Notation / The summary)
- `docs/src/notation.md`:89  (Notation / Estimands)
- `docs/src/notation.md`:19  (Notation / The samples)
- `docs/src/notation.md`:87  (Notation / Estimands)
- `docs/src/notation.md`:12  (Notation / The samples)
- `docs/src/notation.md`:4  (Notation)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-132  Several rows use a symbol that only a later row binds, so a reader going top to bottom meets it unbound: $g$ in the $Z$ row (bound only under Estimands), $\hat\sigma_X$, $q$ and $d$ in the $C(K,d)$ row (bound only under the local fit, the moments and the orders), the subscript $h$ of $\mathcal{R}_h$ (bound only under smoothing kernels, and its role in $\mathcal{R}_h$ is never said), and $S_1^{(j)}$ and $n_j$ in the $\tilde c_j$ row (bound in the rows after it).

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@notation`, decided by an agent.

The page states its own convention that each symbol has one fixed meaning, but it binds these symbols in plain language only after their first use.

**Fix.** Gloss each one at first use (for example, 'a transformed variable $Z = g(X)$ for a smooth function $g$'; 'with $d$ the local-polynomial degree and $q$ the bias order'; '$\mathcal{R}_h$, at bandwidth $h$'), or reorder the rows so the binding row comes first.

- `docs/src/notation.md`:16  (Notation / The samples)
- `docs/src/notation.md`:61  (Notation / Smoothing kernels)
- `docs/src/notation.md`:39  (Notation / The summary)
- `docs/src/notation.md`:28  (Notation / The summary)
- `docs/ARCHITECTURE.md`:455  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.B Power-sum quadrature readouts)
- `src/Expectations.jl`:5
- `src/Expectations.jl`:293  (pushforward_density)
- `src/Expectations.jl`:316  (pushforward_density)
- `src/ConditionalExpectations.jl`:514  (conditional_pushforward_density)
- `src/ConditionalExpectations.jl`:543  (conditional_pushforward_density)
- `docs/ARCHITECTURE.md`:414  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:419  (GBGStreamingStatistics — Architecture and Rationale / 4. Estimand families (estimators — the exportable upper half) / 4.A Kernel-weighting readouts / 4.5 Bandwidth selection (cross-cutting))
- `docs/ARCHITECTURE.md`:740  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:749  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:755  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/ARCHITECTURE.md`:770  (GBGStreamingStatistics — Architecture and Rationale / 11. Decision log)
- `docs/src/notation.md`:60  (Notation / Smoothing kernels)
- `src/Bandwidth.jl`:335  (bandwidth)
- `src/Bandwidth.jl`:336  (bandwidth)
- `src/Bandwidth.jl`:350  (bandwidth)
- `src/Bandwidth.jl`:360  (bandwidth)
- `src/Bandwidth.jl`:363  (bandwidth)
- `src/Bandwidth.jl`:402  (bandwidth)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-133  The page never binds two symbols that are not universally standard. $\oplus$ is the summary's merge operation, the central identity of the $\mathcal{A}$, $\mathcal{B}$ row. $e_0$ appears in the leverage formula.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@notation`, decided by an agent.

No row or sentence says that $\oplus$ is merging, meaning elementwise addition of summaries, or that $e_0$ is the first standard basis vector of $\mathbb{R}^{d+1}$.

**Fix.** Add a gloss: '... $= T(\mathcal{A}) \oplus T(\mathcal{B})$, with $\oplus$ the merge of two summaries' and '$e_0$ the first unit vector'.

- `docs/src/notation.md`:37  (Notation / The summary)
- `docs/src/notation.md`:73  (Notation / The local fit)
- `docs/ARCHITECTURE.md`:109  (GBGStreamingStatistics — Architecture and Rationale / 2. The central principle: the architecture follows established statistics / The algebra of measures — a structural lens)
- pass 1: verifier `LINE_ONLY`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-134  In the $M_n$ row, $n$ is the order of the increment moment ($(X(t+\tau)-x)^n$, $n!$), yet the page binds $n$ as 'the sample count' and says each symbol denotes one concept. A reader who holds the page's own binding reads $M_n$ as a moment indexed by the sample size.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@notation`, decided by an agent.

The samples table binds '$n$ the sample count'. The Orders table assigns orders to $p$, $q$ and $r$. The $M_n$ row nonetheless uses $n$ as an order, and the intro exempts only dummy indices bound inside one formula, whereas $M_n$ is a named quantity.

**Fix.** Index the Kramers-Moyal moments by a letter the page reserves for orders (for example $M_p$, with $M_p/(p!\,\tau)$), or state the exception in the row.

- `docs/src/notation.md`:89  (Notation / Estimands)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-135  The page's opening sentence uses an unhyphenated, unquoted slogan as a compound modifier, so on a first read it parses as 'shows the intended samples in' followed by a stray fragment 'statistics out workflow'.

Rule: `RULE.DOC.GRAMMAR`.  
Reported by `withheld@tutorial`, decided by an agent.

Nothing marks 'samples in, statistics out' as one adjective, and the comma splits the sentence in two.

**Fix.** Write it as a unit, e.g. 'the intended workflow: samples in, statistics out.' or 'the intended samples-in, statistics-out workflow'.

- `docs/src/tutorial.md`:3  (Getting started)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-136  'The consumer' is the package's own term for the calling code, and it is used in step 2 (a code comment) before step 4 defines it as 'the consumer (the calling code)'.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@tutorial`, decided by an agent.

The first use comes in step 2's comment. The gloss comes only in step 4.

**Fix.** Gloss it at first use in step 2 ('the consumer, the calling code, produces the batches'), and leave step 4 bare.

- `docs/src/tutorial.md`:32  (Getting started / The workflow)
- `docs/src/tutorial.md`:53  (Getting started / The workflow)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-137  'The snapshot buffer' appears with the definite article as if the reader knew it, but the page never says what is buffered or why the consumer would hold snapshots.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@tutorial`, decided by an agent.

The only mention is in step 4. It is not defined there or anywhere else on the page.

**Fix.** Say it in plain words, e.g. 'keeping the earlier samples needed to form the pairs, and choosing the lags, is up to the consumer'.

- `docs/src/tutorial.md`:55  (Getting started / The workflow)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-138  'SDE-reading functionals' is the package's own coinage. It is not standard vocabulary, and the Status note gives no gloss for it beyond the parenthetical names.

Rule: `RULE.DOC.TERMS_DEFINED`.  
Reported by `withheld@tutorial`, decided by an agent.

A reader in the field knows Kramers–Moyal coefficients, but 'SDE-reading' does not say what is read from what.

**Fix.** Say 'the Kramers–Moyal functionals that estimate an SDE's drift and diffusion', or drop 'SDE-reading'.

- `docs/src/tutorial.md`:10  (Getting started)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-139  Symbols are used that the page never binds. σ̂_X ('the recovered σ̂_X') is not stated to be the estimated standard deviation of X, and 'recovered' (from what?) is unexplained. Y in the pair (X(t), Y(t+τ)) is never introduced: earlier the page has only the scalar variate X, so the reader cannot tell whether Y is a second variate or X at a later time.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@tutorial`, decided by an agent.

σ̂_X appears only in the bandwidth comment. Y first appears inside the pair formula with no gloss.

**Fix.** Bind each one at first use: 'the standard deviation σ̂_X of X, read off the summary', and 'a response Y (for a transition density, X itself at the later time)'.

- `docs/src/tutorial.md`:49  (Getting started / The workflow)
- `docs/src/tutorial.md`:54  (Getting started / The workflow)
- `docs/src/notation.md`:17  (Notation / The samples)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

### 2026-10-19-documentation-review-01-140  The transition_density call passes two positional query points, 0.0 and 0.0, and the comment p_τ(0 | 0) is symmetric. A reader therefore cannot learn from the page which argument is the conditioning value x and which is the response point y.

Rule: `RULE.DOC.SELF_CONTAINED`.  
Reported by `withheld@tutorial`, decided by an agent.

Both values are 0.0, so the example does not show the order. No prose states it either.

**Fix.** Use distinct values (e.g. transition_density(joint, 0.5, 0.0; ...) with the comment p_τ(0 | 0.5)), or name the order in the comment.

- `docs/src/tutorial.md`:68  (Getting started / The workflow)
- pass 1: verifier `FIXED`  (commit `f49d9a67`)

## Locations that could not be found again — 5

The text each of these named has changed since it was read. They are named rather than
pointed at whatever now sits on the old line, and no fix was attempted at any of them.

## Questions

### 2026-10-15-documentation-review-01-009-code

This finding needs a change only you may make. Will you make it?

Needs doctest edits: the long lines left are trailing comments inside jldoctest blocks (density, ParametricStart, expectation), which a fix pass may not touch; the fix moves each comment onto its own line, output unchanged. The other locations are already wrapped.

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

Needs doctest edits: the remaining 87-92-column lines are trailing comments inside the Histogram1D and Histogram2D jldoctest blocks; moving them onto their own lines changes no output but is a doctest change.

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

Needs math changes: three display formulas put a second equation behind \qquad (bandwidth's R(K) and mu_2(K) row, conditional_density_standard_error's R(K) where-clause, merge!'s pi_k definition); each must be split into aligned rows or the where-clause moved into prose.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-011-code "<option>"`

### 2026-10-15-documentation-review-01-014-code

This finding needs a change only you may make. Will you make it?

Not a defect at this commit: the text it names was already rewritten before this run. not a defect at HEAD: api.md already listed every upper-half module, with no 'through' range. The paragraph was reworded this pass for -059/-060/-061.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
  RULE.ASK.NO_COORDINATES       context: a qualified name in the context: api.md
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-014-code "<option>"`

### 2026-10-15-documentation-review-01-018-code

This finding needs a change only you may make. Will you make it?

Not a defect at this commit: the text it names was already rewritten before this run. not a defect at HEAD: merge! already reads 'accumulate-then-merge equals accumulate-together', with no `=` code span.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-018-code "<option>"`

### 2026-10-15-documentation-review-01-032-code

This finding needs a change only you may make. Will you make it?

Needs a math change: 'only ``\approx`` normalized in ``y``' uses a math span as an adverb; the fix replaces the span with the word 'approximately', removing a math span. The rate sentence of the finding is fixed.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-032-code "<option>"`

### 2026-10-15-documentation-review-01-071-code

This finding needs a change only you may make. Will you make it?

Not a defect at this commit: the text it names was already rewritten before this run. not a defect at HEAD: the support integral already uses \texttt{lower}/\texttt{upper} as its limits, so b means only the asymmetric scale.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-071-code "<option>"`

### 2026-10-15-documentation-review-01-072-code

This finding needs a change only you may make. Will you make it?

Not a defect at this commit: the text it names was already rewritten before this run. not a defect at HEAD: density's estimand already reads 'of the variate X' (only the verb was added, under 110).

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-072-code "<option>"`

### 2026-10-15-documentation-review-01-077-code

This finding needs a change only you may make. Will you make it?

Not a defect at this commit: the text it names was already rewritten before this run. already resolved at HEAD: both 4.5 run-in paragraphs have verbs ('The estimand is', 'the package uses')

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-077-code "<option>"`

### 2026-10-15-documentation-review-01-080-code

This finding needs a change only you may make. Will you make it?

Not a defect at this commit: the text it names was already rewritten before this run. already resolved at HEAD: 'channel' is defined at its first use (L93-94, 'which response variable, each with its own summary')

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-080-code "<option>"`

### 2026-10-15-documentation-review-01-081-code

This finding needs a change only you may make. Will you make it?

Not a defect at this commit: the text it names was already rewritten before this run. already resolved at HEAD: OU, SFI, ULI, Nadaraya-Watson, molecular dynamics and the Pawula quotient are all expanded at first use; the remaining bare 'Q' in the section 2 taxonomy was removed under 115

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-081-code "<option>"`

### 2026-10-15-documentation-review-01-097-code

This finding needs a change only you may make. Will you make it?

Not a defect at this commit: the text it names was already rewritten before this run. not a defect at HEAD: step 1 already says 'of a scalar variate X — the summary every estimator reads', with no 'conditioning variable' or 'canonical'.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-15-documentation-review-01-097-code "<option>"`

### 2026-10-19-documentation-review-01-054-code

This finding needs a change only you may make. Will you make it?

Needs math and doctest changes: \hat f_h and \hat m_h should read \hat f_X and \hat m (or the notation page gains rows for the bandwidth-indexed estimates); the gamma-reference shape \hat\alpha clashes with the table's \alpha(x); merge! and scale! write \texttt{weight} where the table has \lambda; and doctests name a summary h or g. Bindings in words are done.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
  RULE.ASK.NO_HOUSE_TERMS       context: a word of this package's own: hat
```

Answer with: `gbg answer 2026-10-19-documentation-review-01-054-code "<option>"`

### 2026-10-19-documentation-review-01-055-code

This finding needs a change only you may make. Will you make it?

Needs math changes: the restated h*(x) formula in the Histogram1D docstring and the asymmetric estimator formula in gamma_kernel and AsymmetricKernel should leave those docstrings, keeping only a pointer to bandwidth and density, which hold them.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-19-documentation-review-01-055-code "<option>"`

### 2026-10-19-documentation-review-01-058-code

This finding needs a change only you may make. Will you make it?

Needs doctest edits: the remaining long lines are trailing comments inside jldoctest blocks (density, ParametricStart, expectation); moving them onto their own lines changes no output.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-19-documentation-review-01-058-code "<option>"`

### 2026-10-19-documentation-review-01-114-code

This finding needs a change only you may make. Will you make it?

Needs a prose-baseline decision: dropping the '# ≡ f_X' comment from density's signature changes the docstring's first line, which keys its prose-shape row; under a new key its long Arguments paragraph, over the cap since before the cap, is refused. Needs the owner's rebaseline, as for the conditional_pushforward_density signature.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-19-documentation-review-01-114-code "<option>"`

### 2026-10-19-documentation-review-01-123-code

This finding needs a change only you may make. Will you make it?

Needs the owner: setting the math span :sj back as a code span in Architecture section 11 is ruled 'mathematics changed' by the guard (clause 6), although it only respells a Julia symbol set as math.

- **I make the change** (recommended) — the change is yours: the finding goes on your list, and the run leaves it
- **fix it on my ruling** — the fixer makes it as your --direction says; the finding is open again, and a verifier closes it
- **not a defect** — nothing changes, and the finding is closed

This question does not read as the rules ask:

```
  RULE.ASK.ONE_DECISION         decision: more than one sentence: This finding needs a change only you may make. Will you make it?
```

Answer with: `gbg answer 2026-10-19-documentation-review-01-123-code "<option>"`

