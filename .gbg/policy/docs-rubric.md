<!-- gbg:managed begin v2.0.0-rc.58 — rewritten by `gbg sync`; do not edit -->
# Documentation and API — `RULE.DOC.*`, `RULE.API.*`

What a docstring and an API page contain, in what order and in what shape. 47 rules, written from the rule values themselves,
so a rule reworded here was reworded in the program that checks it.

### RULE.API.POSITIONAL
In an operation, an argument that can change the result is positional, defaults included; an argument that cannot change the result — a performance hint, a builder's cost/accuracy tunable — is a keyword. An exception is a `Switch` on this rule in the package's own config, with its reason.

In an operation, every argument that affects the result is positional, defaults included; an argument that cannot change the result is a keyword. An exception is a `Switch` on this rule, with its reason.

Covered by a reader given source + other documents, because this needs meaning.

### RULE.DOC.ADMONITION_TITLE
Every !!! block carries a custom descriptive title. Where this rule's `admonition_titles` is `load-bearing only`, a block need carry one only when text elsewhere relies on the block's content; a general aside may keep its default title.

For this package the rule is what `admonition_titles` says: either every `!!!` block carries a custom descriptive title, or a block carrying one load-bearing claim does while a general aside may stay a bare `!!! note`. A caveat admonition is the last substantive block of a docstring, before the fields list.

Covered by a reader given other documents + siblings.
This package may set `admonition_titles`; unset, it is `always`.

### RULE.DOC.AI_BANNER
A rationale page and a rendered Literate example carry the provenance banner, verbatim; a script that is never rendered is not a page and carries none.

A banner reworded in one package is a banner a reader learns to skim past in all of them.

Covered by a program.
This package may set `banner_pages`; unset, it is `String[]`.

### RULE.DOC.ALIASING_BULLET
An aliasing statement specific to one operation — one that names that operation's own arguments — lives in the bullet for its destination argument. The package-wide aliasing guarantee lives once, in the module docstring, and is not repeated per operation.

The package-wide guarantee lives in the module docstring and is not repeated per operation; an operation states it for itself only when it has an operation-specific reason, or when its contract is the contrary one.

Covered by a reader given siblings, because this needs meaning.

### RULE.DOC.API_GROUPING
One docstring has one home on the built site, and a Base function the package extends is public.

Two homes for one docstring is two things to keep true, and the second teaches a reader nothing.

Covered by a program.

### RULE.DOC.ARGCOV
Every argument in a docstring's signature, positional and keyword, is named as code in its body.

A reader passing an argument the docstring never mentions has to read the source to learn what it does.

Covered by a program.

### RULE.DOC.ARGS_FORM
Three or more arguments, or any keyword, are listed as bullets under an `# Arguments` heading, where several may head one bullet together; fewer are bound in running prose, as in "The arguments are `a`, the first, and `b`, the second."; the two are never mixed. A type lists its constructor's arguments so too, even where its `# Fields` section repeats them.

A reader looking for one argument of many needs a list; a reader of a two-argument call needs a sentence.

Covered by a program.

### RULE.DOC.ARG_BULLET
An argument or field bullet names its argument in a plain code span, never a bold one.

Bold marks emphasis; a name set in bold reads as a heading over the paragraph rather than as the argument.

Covered by a program.

### RULE.DOC.BATTERY
'Check the docs' means the complete battery over every surface, never one check or one file.

When asked to check, review or audit the documentation in any words, run the full matrix — every check across every surface — and report it; do not run one check, sweep one file, and stop.

Covered by the session, as an instruction: followed while working, checked by nothing. Nothing checks it: an instruction to the session: followed while working, checked by nothing

### RULE.DOC.BUILD_VERIFY
A doc edit is verified by the full docs build exiting 0, never by a diff, a lint subset or a log grep.

A lint run standalone cannot see a dangling reference, and grepping a log for error text passes a build that died at load. After editing a Literate example, regenerate it and read the output.

Covered by the session, as an instruction: followed while working, checked by nothing. Nothing checks it: an instruction to the session: followed while working, checked by nothing

### RULE.DOC.CLASS_SWEEP
A flagged instance is a sample of a class: sweep every surface for the class and fix all in one pass.

Fixing only the pointed-at line and waiting to be handed the next instance is the failure this rule exists to stop. After a class fix, enumerate the surfaces from the finding's own text, not from memory, and re-measure.

Covered by the session, as an instruction: followed while working, checked by nothing. Nothing checks it: an instruction to the session: followed while working, checked by nothing

### RULE.DOC.COORD_FREE
Maps and pairings are stated intrinsically, without flattening a structured space to a bare vector; matrix and flattening machinery appears only in the package's formalization document — the file this rule's `formalization` names — and nowhere when none is named. Where a basis form is needed to mirror the code, the docstring states the basis and the space it acts on.

Component notation for a genuine fiber algebra is fine; the rule targets flattening a structured space to a bare vector. Where a basis form is needed to mirror the code, state the basis and the space it acts on.

Covered by a differential_geometry reader given other documents, because this needs meaning.
This package may set `formalization`; unset, it is ``.

### RULE.DOC.COVERAGE_BASES
An operation that reduces a field over a measure documents every base the code supports for it — the reduction over a grid and the same reduction over a point set are both stated where the operation is documented.

A measure reduction over a grid and the same reduction over a point set are both stated where the operation is documented.

Covered by a reader given source, because this needs meaning.

### RULE.DOC.CURRENT
No stale claim: a rule cited by id names a rule that exists, and a retired phrasing is not used.

A citation of something renamed or removed reads exactly like one that resolves, which is why it survives.

Covered by a program.

### RULE.DOC.FENCED_CODE
A code block in a page or a rationale document is fenced with three backticks, never indented; a docstring's indented signature is exempt.

An indented block is code only by context, so a program reading the page cannot tell it from an indented paragraph without the parser's whole reading, and a list item left indented after a blank line renders as code nobody wrote.

Covered by a program.

### RULE.DOC.FIELDS
A type with named fields documents them in a `# Fields` section, in bullets: every bullet a real field, and every field at the head of a bullet, alone or with others; a field named only inside another's bullet has none. A sentence listing them, "The fields are …", is no such section.

A field list that has drifted is worse than none: it names what does not exist and hides what does.

Covered by a program.

### RULE.DOC.GRAMMAR
Every sentence has a subject and a verb; concise means few words, not missing verbs.

Concise means few words, not missing verbs — compact prose that still reads as sentences aloud.

Covered by a reader given nothing, because this needs meaning.

### RULE.DOC.HOME_ORDER
A function docstring is ordered: identity, when-to-use, arguments (under `# Arguments` where RULE.DOC.ARGS_FORM asks for the list), one formula, citation, links, caveats, then `# References`, then `# Fields` last. Implementation internals and design rationale are excluded; so is any term that neither the surrounding sentence nor the package's notation page — the file RULE.DOC.ONE_NOTATION's `notation_page` names — defines.

In full: (1) a plain one-sentence identity in the user's terms; (2) when to use it, versus the alternatives; (3) each argument (positional and keyword) and its effect, bound in prose or listed under `# Arguments` as `RULE.DOC.ARGS_FORM` decides; (4) one formula if a clean one exists, never overflowing; (5) an inline citation for a named or published method and none for a classical one; (6) cross-links; (7) the caveats a user needs, as an admonition where load-bearing — the last block of the home's PROSE. Then the sections that are not prose, in this order and no other: the `# References` section of `RULE.CIT.FULL_INLINE`, then `# Fields`, which is last (`RULE.DOC.FIELDS`). Three rules each used to claim the end of a docstring; the order lives here, and they defer to it. Excluded: implementation internals (a `#` comment), unexplained specialist jargon, architectural rationale and the formal apparatus (the rationale documents).

Covered by a reader given other documents.

### RULE.DOC.INTERNAL_API
An API page puts its internal section last, holds no empty block of collected docstrings, and is never indexed.

Internals in the middle of the public surface, and an index on a hand-grouped page, both make a reader work around the page.

Covered by a program.

### RULE.DOC.JULIA_MANUAL
A docstring or page that uses a keyword splat (`kwargs...`) or a `do` block links the Julia manual's section on it.

A reader who has not met the construct learns it where it is used, and the manual is its official home: a paraphrase of it drifts, and a link does not.

Covered by a program.

### RULE.DOC.LEADING_IDENTITY
Item (1) of a docstring is one plain sentence saying what the operation computes — the physical meaning first, or the estimand first, as this rule's `leading_identity` sets for the package. Whether that sentence is true of the code is RULE.DOC.VERIFIED's question, not this rule's.

For this package the identity leads with what `leading_identity` names: an operation computed in a derived representation still states what it does to the real object first — with the explicit relation — and gives the derived form it uses after; an estimator states the population quantity it estimates, in standard probability notation, before the estimator. Where no clean form exists, say so briefly rather than forcing one.

Covered by a reader given other documents, because this needs meaning.
This package may set `leading_identity`; unset, it is `physical meaning`.

### RULE.DOC.LINKS
Cross-referencing is discriminating: a target is linked at most once in one docstring, at its first genuine mention.

A second link to the same target is invisible in the rendered text and tells the reader nothing new.

Covered by a program.

### RULE.DOC.LINK_UNIFORMITY
In one run of names, either every name is linked or none is.

A list where one name links and the next does not reads as though the difference meant something.

Covered by a program.

### RULE.DOC.LITERATE
In the source file of a Literate example, every full-line # comment is prose that will render, any mathematics in such a comment is LaTeX, and a comment inside control flow is written ## so that it does not render. The rule is judged on the source file; what the page made of it is RULE.DOC.LITERATE_RENDERED's question.

A full-line comment inside a control-flow block splits the code fence, so reserve full-line `# ` for narrative prose at statement boundaries. An underscored identifier in that prose is backticked, or Markdown parses it as an italic span.

Covered by a reader given nothing.

### RULE.DOC.LITERATE_RENDERED
The built page of a Literate example shows its prose as prose: no stray `#`, no italics or code block that a comment's underscores or indentation produced, no mathematics split or broken across comment lines, and no comment landing inside a code fence. The rule is judged on the built page.

What a `#` comment becomes is decided by Literate and then by Markdown, and neither is visible in the source: `a_b and c_d` renders as italic, a comment beginning `# 1.` becomes a list item, four spaces of indentation become a code block, and mathematics split across two comment lines is not mathematics. The habits in the source are RULE.DOC.LITERATE's; this rule reads what they produced.

Covered by a reader given rendered.

### RULE.DOC.MODULE_MAP
A module docstring lists one linked one-liner per operation family and no formula.

Its order is identity, then the contents list, then an optional trailing remark.

Covered by a reader given siblings.

### RULE.DOC.NO_DEAD_REFS
A docstring cites no section number and no inaccessible source; it reads complete at the REPL.

At the REPL there are no links: a section number cannot be followed, and it is renumbered besides.

Covered by a program.

### RULE.DOC.ONE_HOME
One concept, one home: one docstring per operation on the bare binding; the ! form owns the formula and is listed first; every other variant is a thin wrapper that @ref-links the home and carries no second formula.

Per-scheme or per-representation specifics live in the scheme's or the representation's own docstring, linked. Across an `f`/`f!` pair the in-place `!` form is the home — it owns the formula and the kernel and is listed first — and the allocating form and any operator alias is a thin wrapper that leads with the one-sentence meaning and `@ref`-links the home, carrying no second copy of the formula.

Covered by a reader given siblings.

### RULE.DOC.ONE_NOTATION
One symbol denotes one concept and one concept is denoted by one symbol, as the package's notation page records — the file this rule's `notation_page` names; a local use that conflicts with the notation page is the one that changes.

Before introducing a symbol, check it is not a synonym for one already in use, and check the section for the glyph and confirm it is free. When an edit introduces a symbol, retire the old one.

Covered by a reader given other documents + siblings, because this needs derivation.
This package may set `notation_page`; unset, it is `docs/src/notation.md`.

### RULE.DOC.REF_ABUT
A cross-reference never abuts a second code span; the whole call is the link text.

The two render as adjacent monospace chips with a seam, as though the name had been cut in half.

Covered by a program.

### RULE.DOC.REGISTER
A tunable is an option, a setting or a keyword argument, never a `knob` or a `tweak`.

The house register is plain and exact; `knob` is neither, and reads as though the value did not matter.

Covered by a program.

### RULE.DOC.RETIRED_CLAIMS
A phrase a package has retired does not reappear in its documentation.

A claim retired once is retired because it was wrong, and the way a wrong claim comes back is that somebody copies the sentence that carried it. A DECISION LOG quotes what it retires, in order to record that it was retired, so this rule fires on a package's own history and is right to: keep the log in a file of its own and switch the rule off for that file, never for the design document the log was buried in, which is where the guard is worth most. The check is only as good as the list, and the list is kept by the party being checked: nothing in a repository knows that a claim was retired, so a claim dropped from the list is invisible here, and a ban whose reason was itself later overturned stays in force until someone notices.

Covered by a program.
This package may set `retired_claims`; unset, it is `GBGCodingTools.Model.RetiredClaim[]`.

### RULE.DOC.RETURNS
A type docstring leads with its constructor. A function docstring whose identity sentence does not name the result states it with a '-> Type' arrow or a 'Returns x.' sentence. Whether the stated result is what the code returns is RULE.DOC.VERIFIED's question, not this rule's.

A mutator whose signature carries no arrow says what it returns in a sentence before any caveat; a signature carries an arrow only where the return type is fixed by the signature; arrow and sentence never both appear.

Covered by a reader given nothing.

### RULE.DOC.SECTION_ORDER
A caveat admonition at a docstring's margin comes after its `# Arguments` section and before its `# Fields` section.

`RULE.DOC.HOME_ORDER` orders a docstring: the arguments, then the caveats that close its prose, then `# References`, then `# Fields` last. These two places are where a caveat most often lands wrongly, and a program can find them; an admonition indented under a list item is that item's.

Covered by a program.

### RULE.DOC.SELF_CONTAINED
A docstring, a page, and the one rationale document this rule's `architecture_doc` names read top to bottom without another file: every symbol is bound in plain language at first use. Every other rationale document is exempt, since it explains one architecture and may assume the rest of itself.

Later uses may be bare; universally standard mathematics may stand unglossed. The exemption is not indulgence: holding a design document to a cold read means re-binding the same vocabulary on every page, which is what the one-home rule forbids.

Covered by a reader given nothing, because this needs reading order.
This package may set `architecture_doc`; unset, it is `docs/ARCHITECTURE.md`.

### RULE.DOC.SENTENCE_OPENER
A sentence opens with a capitalised ordinary word, never a lowercase identifier or a formula; bullets, comments, headings and fenced code are exempt.

At a full stop a reader cannot tell a new sentence beginning with a lowercase name from the old one running on.

Covered by a program.

### RULE.DOC.SHAPE
A paragraph of prose stays inside the house cap, an admonition count stays inside its own, and prose that was here before the cap may shrink but never grow.

Prose grows by accretion: each addition is small, nobody rewrites, and the paragraph nobody can read was never written in one sitting.

Covered by a program.

### RULE.DOC.SIBLING_LINK
The first mention of a sibling package links to that package's site; later mentions are bare.

A reader meeting the name for the first time must be able to reach it, and a fourth link teaches nothing.

Covered by a program.

### RULE.DOC.SIGNATURE_FIDELITY
The signature a docstring displays and the methods the module defines agree, argument for argument.

A displayed signature is what a reader calls; one the method table has moved past is a call that fails.

Covered by a program.

### RULE.DOC.SIGNATURE_MIRROR
A docstring's signature block fits the house column of 88 characters, as rendered, and repeats only as much of the method as a reader needs.

A signature wider than the column wraps in the rendered page and in every terminal that shows it.

Covered by a program.

### RULE.DOC.STEER
Every ! form with an allocating sibling carries the in-place steer sentence, verbatim from the template this rule's `why` gives, with `destination` replaced by the ! form's own destination argument and `f` by the allocating sibling's name.

The sentence is a template, not a link: `destination` is replaced by the argument the signature actually has and `f` by the allocating sibling's own name, written as an `@ref` link, so the link resolves. The template: Use this in-place form when `destination` is already allocated; the allocating `f` returns a fresh result instead.

Covered by a reader given siblings + other documents.

### RULE.DOC.TABLE_TOP_LEVEL
A Markdown table is top-level, never indented inside a list item.

Markdown reads an indented block as code, so an indented table renders as the pipes and dashes themselves.

Covered by a program.

### RULE.DOC.TERMS_DEFINED
A term of art is defined in a clause, or @ref-linked to its definition, before its first use on the same page. The test is a cold read of that page alone, and it is asked of docstrings, of pages, and of the one rationale document this rule's `architecture_doc` names; every other rationale document is exempt.

A docstring can bind every symbol and still be unreadable to someone who has not read the module, because its terms of art are nouns, not symbols; parse it with the module's context withheld, and a noun that silently assumes the package's context is the defect. A design document is read the other way round, by someone who has the architecture in hand, and holding it to a cold read means defining the same nouns on every page.

Covered by a reader given nothing, because this needs meaning.
This package may set `architecture_doc`; unset, it is `docs/ARCHITECTURE.md`.

### RULE.DOC.TWO_TIER
A module docstring is a map: one linked one-line entry per operation family, and no formula. A function docstring is a home: the one place where that operation's identity, arguments and defining formula are stated in full.

The map states what the module is, its scope and place in the architecture, and its invariants — then a linked one-liner per operation family, and no formulas; every formula lives in the operation's own docstring.

Covered by a reader given siblings, because this needs reading order.

### RULE.DOC.VARIANT_PAIR
Across an f/f! pair or a getter/setter pair, the defining mathematics appears once, in the docstring of the member that is the home under RULE.DOC.ONE_HOME; the other member carries no display formula. The pair is one unit for this rule.

Typed twice, the math drifts — one member writing a partial one way and its twin another — and a within-string check cannot see it; a lane that groups each operation's generics and diffs them can.

Covered by a reader given other documents.

### RULE.DOC.VERIFIED
Every formula, normalization, adjoint relation and cross-scheme equivalence stated in the documentation is true of what the code computes, established by expanding the statement against the code, never taken on the docstring's word.

A plausible 'X is just Y' is a false claim until the expansion matches. Run this as a positive verification — enumerate every formula, constant, relation and quantifier and verify each — not as a read-through; a false-but-plausible statement is as invisible as an unbound symbol. Category and type errors are their own class: an object conflated with its representation, an operator with its parameter, a rate with its time-integral, a scalar with the field it reduces; type every object from the code.

Covered by a reader given source, because this needs meaning.

### RULE.DOC.XLINK
A name another module owns, mentioned on its own, is linked at its first mention in a docstring, qualified so the link resolves; a name inside a code expression owes no link.

A reader meeting a name from elsewhere has nowhere to go: the docstring in front of him does not define it.

Covered by a program.

<!-- gbg:managed end -->
