---
description: Run the COMPLETE documentation-quality battery — every check across every doc surface, fix all findings, rebuild.
---

Run the **entire** documentation-quality battery. This is the full matrix of (every check) × (every
surface), run to completion in one invocation — not a single check, not a single file. Do not stop after
one check or one file, and do not ask which check to run. Cover every cell, fix every genuine finding in
one pass, rebuild, and report the matrix.

## Surfaces (ALL of them)

1. `docs/ARCHITECTURE.md` — every section §0–§11 and §References.
2. `docs/FORMALIZATION.md`, `docs/TESTING.md`, `docs/PERFORMANCE.md`.
3. Every docstring in `src/*.jl` — `Accumulators`, `Histograms`, `DensityEstimation`, `Regression`,
   `ConditionalDensity`, `MarkovInference`, `Bandwidth`, and the umbrella `GBGStreamingStatistics.jl`
   (Public **and** Internal docstrings; `#` comments are exempt — not rendered).
4. Every page in `docs/src/` — `index.md`, `tutorial.md`, `api.md`, and all `docs/src/api/*.md`
   (hand-written prose; skip `@autodocs` / `@index` / `@contents` blocks).

## Checks (ALL of them, on every surface)

Run each as a POSITIVE procedure — enumerate, then verify each item. A defect that is absent-by-omission
(an unbound symbol; a false-but-plausible formula) does not announce itself on a read-through, so passive
reading is not a check.

1. **Correct** — enumerate every formula, normalization constant (the `1/n`, `1/(nh)`, `∫K=1`, the kernel
   moments `μ₂(K)`, the roughness `R(K)`, the bias/variance constants), consistency claim, complexity
   claim, and quantifier (*exact / all / never*); verify EACH against the implementation in `src/` (read
   it, quote the decisive line) or by derivation. **Extended for this package: every estimator must be
   consistent for its stated estimand** — check the population definition is right and the estimator
   actually converges to it (a stray `1/h` or kernel-moment factor is the easy slip). The verified
   statements live in `docs/FORMALIZATION.md`: check each nontrivial prose claim **against its entry**
   (the chain is *prose ↔ formalization ↔ test*), and flag any nontrivial statement with **no entry**, or
   an entry with **no test link**, as a coverage gap.
2. **Self-contained** — enumerate every symbol, sub/superscript, operator name, and bold/blackboard
   letter; cite the clause that binds each; flag any unbound. **Emit the enumeration as a `symbol |
   binding clause (or UNBOUND)` table** — never report "Self-contained: clean" without it (an unbound
   symbol that reads as obvious — a bandwidth `h`, a stray `\sigma` — is exactly what a hand-waved pass
   misses). A symbol that is really a code identifier is a *Formatted* defect (it should be backticked or
   `\texttt`-wrapped), not a binding. (Universally standard math — `∫`, `ℝ`, `π`,
   `𝔼`, the imaginary unit `i`, plain summation indices — needs no binding; the estimand symbols
   `f_X`, `m(x)`, `p_τ`, `D⁽ⁿ⁾`, `σ²(x)` are bound on first use.)
3. **Estimand first** — every estimator surface states the **population quantity, in standard probability
   notation, with its defining formula, BEFORE** the estimator that approximates it, and binds the code
   name to its math symbol. Flag any docstring/page that leads with the estimator.
4. **Formatted (bidirectional)** — math is LaTeX (`` ``…`` `` inline, ```math``` display), never backticks
   (incl. ASCII-math that carries no glyph); every code identifier / type / module name is backticked,
   never bare; check BOTH directions. **Code identifiers never sit inside a `` ``…`` `` math span** —
   they render as a product of italic single letters (`xlo` → x·l·o); backtick them, or wrap in
   `\texttt{}` if a name must appear inside a formula. Markdown tables top-level, never indented in a list. Inline
   ``math`` stays on ONE source line — a wrapped continuation beginning with a list/blockquote leader
   (`- `/`* `/`> `/`1.`) is parsed as a list item and silently drops the span out of math mode
   (`docs/make.jl` lints for this; Documenter renders the broken page without complaint). **A formula
   is a single LaTeX span** — flag any ``math`` span split onto an adjacent `code` span that continues
   it (e.g. ``1 \le i \le`` `nbins`, which should be ``1 \le i \le \texttt{nbins}``); code names inside a
   formula are `\texttt{}`-wrapped. **Backticks hold only well-formed Julia** — a name, type, call, or
   expression; flag math notation in backticks (a half-open interval like `[lo, hi)`, a product `a × b`,
   an inequality, or a `·`/`−`/`≤` glyph), which belongs in LaTeX. **No horizontal scroll** — a rendered
   code block and a display-math block do not auto-wrap: flag a docstring **signature** line wider than
   ~88 source columns (rendered > ~84; split it with a hanging indent aligned under the first argument,
   the `inspiration/GBGPeriodicFields.jl` `anterpolate!` convention), and a **display formula** that
   packs several equations on one line behind `\qquad`/`\quad` (break it with `aligned`/`gathered`/`split`).
5. **Current & complete** — every claim matches the code's present state; no stale *future / deferred /
   none-exist* claims about what is now implemented; and a docstring or summary that names only a *subset*
   of what the code now provides is stale even if literally true. Check each module/operation/summary
   still names everything added since it was written, against the current `export`s / method set. (At
   scaffold stage, the honest "declared, not yet implemented" notes are correct — flag them only once the
   code lands.)
6. **Cited** — context-standard author–year; ≤3 authors named, `First et al.` for 4+; never cite an
   inaccessible source; inline references verbatim-consistent with the `ARCHITECTURE.md` §References entry.
   **A multi-citation reference block is a bullet list** — one citation per `-` line under a `References:`
   label and a blank line, never a semicolon-joined paragraph (a lone citation stays inline as
   `Reference: …`); flag any `Reference(s):` block that crams two or more citations into one paragraph.

`CLAUDE.md` (the doc-edit checklist and the standing notation/structure rules) holds the authoritative
definitions; this command runs all of them.

## Mechanism

- Fan out subagents across the surfaces in parallel — each runs ALL checks on its surface via the positive
  procedures above, reading `src/` as ground truth, and returns genuine findings with evidence and a
  proposed fix.
- Run the deterministic grep guards too: the table guard `grep -rnE '^[[:space:]]+\|[-:| ]+\|' docs` (must
  be empty); a stray-`$` scan over `docs/` (math must be `` `` ``/```math```, not TeX `$…$`); a
  bare-module-name scan; a glyph-vs-backtick scan; a **code-in-math scan** — enumerate the `` ``…`` ``
  math spans and flag any containing a 2+-letter lowercase run that is not a `\command` nor
  `\texttt`/`\mathrm`/`\text`-wrapped (those render as letter-products, e.g. `[lo, hi)`, `scale`);
  a **backticked-non-Julia scan** — flag backtick spans holding math not code: a single-backtick
  half-open interval (``grep -rnE '`\[[^`]*\)`' docs``) and the glyphs `×`/`·`/`−`/`≤`/`≥`/`…` inside a
  backtick span (drop `#`-comment hits and valid array/tuple literals); a **split-formula scan** — a
  math span ending in `\le`/`\ge`/`<`/`>`/`\dots` immediately followed by a backticked `code` span,
  which must instead be one `\texttt{}`-using LaTeX span; triage the legitimate LaTeX command bodies; a
  **crammed-references scan** — `grep -rnzE 'References?:[^"]*;' src` flags a `Reference(s):` docstring block
  that joins two citations with a `;` (anywhere before its closing `"""`), which must instead be a `- `
  bullet list, one citation per line (`CLAUDE.md` Citations); a **no-horizontal-scroll scan** —
  `awk '/^"""/{d=!d;s=d;next} s&&/^[^ ]/{s=0} s&&/^    [A-Za-z_].*\(/&&length>88{print FILENAME":"NR}' src/*.jl`
  flags a docstring signature wider than 88 source columns (rendered > 84; split with a hanging indent
  under the first argument), and `grep -rnE '\\q?quad.*\\q?quad' src docs` flags a math line queuing two
  equations behind spacers (rebreak with `aligned`/`gathered`/`split` unless already `\\`-separated).
- Triage (drop false positives — universally standard math; kwarg names correctly backticked), then fix
  EVERY genuine finding in one pass. A flagged instance is a sample of a class — fix the class, not the line.
- Rebuild: `julia --project=docs docs/make.jl` (the inline-math lint + doctests + cross-references +
  render) — must be clean.
- Commit + push (no `[skip ci]`; the docs site must redeploy).

## Report

Report the full matrix: which surfaces × checks were run, the findings (class and count) with fixes, and
what was verified clean — including the `FORMALIZATION.md` coverage (every nontrivial operation has a typed
entry, each test-linked). Done only when every cell is covered and the build is clean.
