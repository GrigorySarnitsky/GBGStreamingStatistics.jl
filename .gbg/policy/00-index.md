<!-- gbg:managed begin v2.0.0-rc.57 — rewritten by `gbg sync`; do not edit -->
# The rulebook

One text for every package that pins GBGCodingTools, at the version this package pins. Every rule has
a stable id, `RULE.<FAMILY>.<NAME>`, which is never renamed and never reused: a finding, a switch and
a report all name a rule by it.

This directory is WRITTEN, not edited. Each rule is a Julia value carrying its statement, the check
that decides a violation, and the examples that show the check matches the statement; these files are
rendered from those values by `gbg sync`, so a rule cannot be reworded here without being reworded in
the program that enforces it, and the two cannot drift.

A package that a rule does not fit does not edit it. It declares a `Switch` in its own `.gbg/config.jl`,
naming the rule and the reason, for the whole package or for named files; `switches.md` lists what this
package has declared, and a session reads a switched rule as the switch states it. A rule a package
adds for itself is a `PatternRule` in the same file, and is nobody else's.

The files:

- `docs-rubric.md` — What a docstring and an API page contain, in what order and in what shape. (47 rules)
- `latex-idiom.md` — What renders, and how mathematics is written so that it renders the same way everywhere. (18 rules)
- `citations.md` — Author–year in the prose, the full reference where the reader can reach it, a DOI on every title. (6 rules)
- `testing-protocol.md` — What a test is for, where its oracle comes from, and when a thing is done. (12 rules)
- `instruments.md` — What a session runs or sends unasked, and how it reads what a scan prints. (3 rules)
- `git-conventions.md` — What a commit carries, where it goes, and what is never committed. (9 rules)
- `asking.md` — One decision, short, in his words rather than this package's, with every identifier moved out of the sentence. (5 rules)
- `publication.md` — What a repository may not publish when it is made public: no credential, nothing personal, nothing never meant for a reader. (9 rules)
- `glossary.md` — the terms the rules use, bound once.
- `switches.md` — what this package has switched off, and why.

@.gbg/policy/docs-rubric.md
@.gbg/policy/latex-idiom.md
@.gbg/policy/citations.md
@.gbg/policy/testing-protocol.md
@.gbg/policy/instruments.md
@.gbg/policy/git-conventions.md
@.gbg/policy/asking.md
@.gbg/policy/publication.md
@.gbg/policy/glossary.md
@.gbg/policy/switches.md
<!-- gbg:managed end -->
