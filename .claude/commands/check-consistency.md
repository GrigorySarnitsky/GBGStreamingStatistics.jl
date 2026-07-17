# /check-consistency — the fresh-eyed multi-hat RENDERED-docs panel

Fire an expert panel of fresh-eyed agents at the **rendered** documentation and collect every
internal-consistency and expert-instinct defect they find. This is the **judgment-level** counterpart to
`/check-docs`: where `/check-docs` runs a *deterministic* battery (greps, the build lint, the
argument-coverage sweep) against the **source**, `/check-consistency` reads the **rendered site as a
reader would** and hunts what a grep cannot — an estimand defined two ways across pages, a statistical or
category error, a normalization that disagrees between two formulas, a structure/design smell, an
unreachable citation, a cold-read gap.

Each hat is a **fresh agent** that reads ONLY the reader-text manuscript (no `src/`, no `CLAUDE.md`, no
memory — it judges the docs as written, from its own expertise) and banks its findings to a file. The one
exception is the `latex-typesetter` hat, whose substrate deliberately also grants the authored source and
`CLAUDE.md §Style`: TeX idiom defects live in the markup (a `.jl` docstring's doubled backslashes, an
unbraced group, a `\left`), so it reads the rubric first, then the source, and reports at source file:line. The
panel is **collect-only**: it never edits `src/` or `docs/`; you act only when the user says fix. It runs
**sequentially**, so a token-out loses at most the current hat; the hat order is **priority** (most
critical first).

## What it catches — and what it does not

Catches: cross-page contradictions, statistical/estimand/category errors, method disagreements, design and
API smells in the *documented* interface, structure non-uniformity, citation gaps/overclaims, and
self-containedness ("cold read") failures. These need judgment and cross-referencing.

Does **not** reliably catch: mechanical, pattern-matchable idioms — gratuitous `\,`, unbraced
`\tfrac`/`\boldsymbol`, formula-in-backticks, a code identifier set inside a math span, inline math split
across a `code` boundary. A hat that *reads* the TeX may not *flag* these (a panel once read the `\,`
overuse and waved it through). **Those belong to `/check-docs`'s deterministic scans** — do not lean on
the panel for them.

## The panel (11 hats)

Substrate **FULL** = every rendered page incl. the Design & Internals pages; **SITEONLY** = the Design &
Internals pages (Architecture / Formalization / Testing / Performance / Internals) withheld (so the two
no-design hats judge self-containedness and structure of the user-facing surface without the design docs
as a crutch).

| # | Hat | Lens | Substrate |
|---|-----|------|-----------|
| 1 | `consistency-matrix` | the spine — matrix every recurring entity (estimand, symbol, formula, constant, method, docstring shape) and diff it across all pages; the category/type-error pass | FULL |
| 2 | `software-architect` | the *documented* design as a senior Julia dev — signatures, dispatch, in-place/bang convention, decomposition, single-source-of-truth, decision-log soundness | FULL |
| 3 | `statistics-theoretician` | every estimand definition / estimator / consistency claim / constant / rate vs established estimation theory *and* itself | FULL |
| 4 | `streaming-stability-practitioner` | the one-pass / fixed-memory / no-rerun claims, numerical stability of the streaming updates, the mergeable-monoid/cone laws, the allocation-free hot path | FULL |
| 5 | `structure-consistency` | infer the implicit (estimand-first) docstring template, then flag every outlier against it | SITEONLY |
| 6 | `latex-typesetter` | every formula's TeX idiom read as source (operators, upright differentials, gratuitous `\,`, unbraced groups, the `\left`/`\right` ban, no manual spacing around relations, `\mathbb` only on capitals/`k`) — read against `CLAUDE.md §Style`, reported at source file:line | FULL + source + `CLAUDE.md` |
| 7 | `api-user` | build/drive the library from the docs alone — arg coverage, return shapes, exported-vs-documented agreement, estimand-first | FULL |
| 8 | `normalization-convention-auditor` | every normalization factor correct *and* identical across pages; the density-vs-mass-vs-count and estimand-vs-estimator-vs-estimate category distinctions | FULL |
| 9 | `citations-integrity` | consistent author-year, inline↔References matching, no inaccessible source, no overclaim | FULL |
| 10 | `style-editor` | internal uniformity (majority style) **+** the raw built HTML for silent-render failures | FULL + raw HTML |
| 11 | `cold-reader` | a first-year statistics grad student's top-to-bottom cold read — undefined symbols/terms, self-containedness | SITEONLY |

The hats live in `.claude/commands/check-consistency/panel.js` (their personas are the source of truth;
this table is the map). Scale the panel to the ask: this full 11-hat run is the "thorough audit"; for a
quick check, run a subset.

## Run it

Let `$SCRATCH` be this session's scratchpad and `$REPO` the repo root (absolute). The large, ephemeral
manuscripts go in `$SCRATCH`; the findings pile is the repo's `.claude/doc-review/`.

1. **Build the docs fresh** — the panel reads the *rendered* output, so it must be current, and the build
   lint must be clean first:
   ```bash
   julia --project=docs docs/make.jl          # inline-math lint + checkdocs + doctests must pass
   ```
2. **Consolidate** the built HTML into the two reader-text manuscripts (needs `beautifulsoup4`):
   ```bash
   python3 .claude/commands/check-consistency/consolidate.py docs/build \
       "$SCRATCH/manuscript-full.md" "$SCRATCH/manuscript-siteonly.md"
   ```
   Confirm the reported page count matches the site and no page is empty.
3. **Archive** any prior findings so a re-run starts clean (keep the last round for diffing):
   ```bash
   if [ -n "$(ls -A .claude/doc-review 2>/dev/null)" ]; then
     rm -rf .claude/doc-review-prev && mv .claude/doc-review .claude/doc-review-prev
   fi
   mkdir -p .claude/doc-review
   ```
4. **Run the panel** — a background `Workflow` on the companion script, sequential, collect-only. Pass
   **absolute** paths (each hat `Read`s them):
   ```
   Workflow({ scriptPath: ".claude/commands/check-consistency/panel.js",
              args: { full:     "$SCRATCH/manuscript-full.md",
                      siteonly: "$SCRATCH/manuscript-siteonly.md",
                      rawhtml:  "$REPO/docs/build",
                      pile:     "$REPO/.claude/doc-review" } })
   ```
   Each hat banks `.claude/doc-review/NN-<key>.md` and returns a structured findings object; you are
   notified when the workflow completes (~2 M tokens, tens of minutes for the full 11).
5. **Consolidate the findings** into `.claude/doc-review/00-SUMMARY.md` **and an Artifact dashboard**:
   read every `NN-*.md`, **dedup** (cluster cross-hat — a defect found by many hats is a strong signal),
   **rank** by severity × corroboration, and record which batteries came back **clean** and which are
   confirmed **false positives**.

## Discipline (the hard-won rules — follow them)

- **Collect-only.** The panel proposes; the user disposes. Do not touch `src/` or `docs/` during a panel
  run, and do not "quietly fix while I'm here."
- **Verify before you trust — a rendered-text read is not ground truth.** Before reporting any Tier-1/2
  finding, confirm it against the **current `src/`** and the **built HTML** (quote the decisive line).
  Rendered-text reads have produced *false* blockers: an ordered list that renders numbered but reads as
  unnumbered in the consolidation (so a "decision N" cross-reference looks dangling); a grep-escaping slip
  that "found" zero of something that was plentiful. Re-check the crux, do not relay it raw.
- **Judgment, not mechanics.** See "what it does not catch" — route mechanical idioms to `/check-docs`.
- **Don't launder a flagged smell.** Never close a finding by declaring it "by design" / "for symmetry" /
  "it's convenient" (the exact non-justifications the invariants forbid). Fix it, or escalate it to the
  user with the trade-off — the user decides.
- **Record dispositions and false positives** in `00-SUMMARY.md` (a short block at the top): what was
  fixed, escalated, left-with-a-first-principle-reason, and confirmed-not-a-defect. A later run reads this
  so it does not re-flag settled items; a still-open finding stays visible.
- **A flagged instance is a sample of a class** — when you do fix (later, on the user's word), fix the
  class across every surface, not just the quoted line.

## Output

- `.claude/doc-review/NN-<key>.md` — one file per hat (the raw banked findings).
- `.claude/doc-review/00-SUMMARY.md` — the consolidated, deduplicated, ranked report with dispositions,
  the verified-clean batteries, and the confirmed false positives.
- An **Artifact** dashboard (filterable by severity/surface) for the user to scan.
- `.claude/doc-review-prev/` — the previous round, kept for diffing.

`.claude/doc-review*` is untracked working notes (never commit it — stage tracked files with `git add -u`,
never `-A`). This command file and its `check-consistency/` scripts ARE tracked.
