<!-- gbg:managed regressions-note v2.0.0-rc.57 — rewritten by `gbg sync`; do not edit -->
# `test/regressions/` — the failing tests of defects that were fixed

Each file here began as the FAILING TEST of a bug hunt: a self-contained `@testset` asserting the
correct behaviour, written before the fix and observed to fail on the tree that had the defect. It is
moved here unchanged when the fix lands, so the test that proves the fix is the test that proved the
bug.

That "unchanged" is the whole discipline. A test edited to suit the fix has stopped being evidence: it
no longer records what was wrong, only what the fixer decided to assert afterwards. If one needs
changing to pass, either the fix is wrong or the test was — settle that, do not adjust the test.

Each file must also:

- be runnable on its own (`julia --project test/regressions/<file>.jl`), since it was, as a witness;
- assert the correct behaviour against an ORACLE — an identity, a documented contract, a hand
  computation — never against another operation of the same package alone, which proves only that
  two things agree;
- carry a header comment naming the defect, the date, and the commit that fixed it.

`escalated/` holds the failing tests of defects the owner has not yet ruled on. A runner that walks
this directory must not descend into it: a confirmed defect nobody has approved a fix for must not
turn a package red.

Naming: `<operation>-<one-line-symptom>.jl`, as in `weighted_inner-complex-weight-accepted.jl`.
