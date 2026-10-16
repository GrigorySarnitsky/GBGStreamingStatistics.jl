<!-- gbg:managed begin v2.0.0-rc.76 — rewritten by `gbg sync`; do not edit -->
# Glossary

The terms the rules use. A rule may use any of these without defining it; a term not here is defined
where it is used, or it is not a term of art.

- **consumer** — a package that pins this one and receives its wiring.
- **rule** — one statement with one check and the examples that show the check matches it. Its id is
  `RULE.<FAMILY>.<NAME>`, stable for the life of the rule.
- **check** — what decides that a rule was broken: a program, or an agent where no program can.
- **clause** — one requirement of a statement. A rule with three clauses owes three breaking examples
  and three complying ones.
- **switch** — a package turning a rule off, wholly or for named files, with its reason. It replaces
  the exemption ledger: there is no expiry and no ceremony, only a reason a reader can weigh.
- **deferral** — a package saying a rule is owed here and not yet met, with what would close it. A
  switch is an answer; a deferral is a question with a condition on it. It fails no build.
- **file set** — which of a package's files a rule runs over: docstrings, code, pages, rationale,
  Literate examples, tests.
- **finding** — one violation at one location, with the rule that found it.
- **location** — where a finding is: the file, the node (a binding, or a heading path), the quote as
  it read then, and the line as a hint only.
- **generated file** — a file written by a Julia function of this package. Editing one by hand is
  drift, which the wiring check reports and `gbg sync` repairs.
- **pin** — the one version of this package a consumer uses, named in its sub-environments.
- **intention** — whether a change or a finding is in documentation or in code. Documentation is safe
  to fix without asking; code is not.
- **operation family** — an operation together with its variants: the allocating and in-place
  forms, a getter and its setter, the methods of one name across representations. A module map
  lists one entry per family; the family's home docstring is the one `RULE.DOC.ONE_HOME` names.
- **instrument** — one review procedure with one question, started only on the owner's request.
- **run** — one execution of one instrument, with its own record, its own report and its own id.
  A run never starts another.
- **stage** — one step of a run, in a fixed order. The session launches each; a stage out of order
  is refused.
- **reader** — an agent given some of the documentation and asked what is wrong with it. Which
  readers a review has follows from the rules, not from a list.
- **bundle** — what a reader is given: one built page alone, the whole documentation with the
  source, or the built site. A rule's tag says which bundle can judge it.
- **unit** — the one page a reader of the first bundle is given. That reader is launched once per
  unit, so that what it cannot see is what the rule needs it not to see.
- **pass** — one round of fixing and verifying inside a run. A documentation review makes one,
  after which what is still open goes to the owner's queue.
- **verifier** — an agent that judges whether the fixes of one place fixed their findings, with
  the text as it now stands and the diff of the pass in front of it.
- **cover** — what decides a rule is broken: a program, a reader, or nothing, the last being an
  instruction the session follows while working.
- **discipline** — a field of expertise a package declares it needs a reader to have. A rule naming
  a discipline the package does not declare is a gap.
- **gap** — a rule no reader of this package can be given what it needs to judge. A run says so at
  every launch and first in its report.
- **settled** — a decision the owner has already given, kept in `.gbg/settled.jl` and given to every
  reader as "do not report this again".
<!-- gbg:managed end -->
