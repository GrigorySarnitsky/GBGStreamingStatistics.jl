<!-- gbg:managed begin v2.0.0-rc.57 — rewritten by `gbg sync`; do not edit -->
# Citations — `RULE.CIT.*`

Author–year in the prose, the full reference where the reader can reach it, a DOI on every title. 6 rules, written from the rule values themselves,
so a rule reworded here was reworded in the program that checks it.

### RULE.CIT.BLOCK_LIST
References are listed one to a line, never run together with semicolons.

A reader looking for one reference in a semicolon-separated run has to parse the whole paragraph.

Covered by a program.

### RULE.CIT.DOI
Every reference carries a clickable DOI, or an Open Library link for a book, on its title. A work whose entry states it has no DOI — the words "no DOI", in any case — may carry on its title instead an https link to its stable public record.

A DOI a reader cannot click is a string to retype; one written as a bare label is indistinguishable from prose.

Covered by a program.

### RULE.CIT.FULL_INLINE
A docstring cites author–year in its prose and carries the full reference in a `# References` section, identical to the canonical entry in title and DOI.

A full reference mid-sentence takes the reader's eye out of the sentence, and a hand-copied entry drifts from the one the package keeps.

Covered by a program.
This package may set `references`; unset, it is `docs/ARCHITECTURE.md`.

### RULE.CIT.NO_INACCESSIBLE
Never cite an inaccessible source on any user-facing surface: not a gitignored monograph, an accompanying thesis, an appendix or figure the reader cannot open, nor working notes. Whether a link answers is RULE.CIT.RESOLVES's question, not this rule's.

A citation the reader cannot follow is the author telling themselves where they got it; the reader is left with an assertion and nowhere to check it.

Covered by a reader given other documents, because this needs meaning.

### RULE.CIT.RESOLVES
Every DOI or Open Library link on a title answers, and answers with the work the title names.

A link is a promise that the reader can reach the work. A DOI that answers 404, or answers with another paper's record, breaks it and nothing in the text looks different. This reaches the network, so it runs only under `gbg check --network` — on request, never in a build or a test — which is why the reader's rule about inaccessible sources does not ask it.

Covered by a program, which reaches the network and runs only under `gbg check --network`.

### RULE.CIT.STYLE
An in-text citation is author–year. An author-prominent citation puts the year in parentheses. Up to three authors are named in full; four or more are written 'First et al.'. The author count is read from the package's reference list, the file RULE.CIT.FULL_INLINE's `references` names.

Keep enough surnames to disambiguate two works that would otherwise collide. An eponymous method name takes no year unless the work is actually cited.

Covered by a reader given other documents, because this needs meaning.

<!-- gbg:managed end -->
