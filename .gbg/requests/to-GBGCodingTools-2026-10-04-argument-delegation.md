Status: open — approved by the owner on 2026-10-04; sent to the GBGCodingTools session the same day.

# To GBGCodingTools: admit argument delegation to a named sibling in `ARGCOV` and `ARGS_FORM`

SS has fifteen standard-error readers — `density_standard_error`, `density_with_standard_error`,
`expectation_standard_error`, `conditional_density_standard_error`, and their siblings — whose signatures
are the parent estimator's, less at most one option. Since 2026-07-18 (SS commit `30854bb`, an owner
decision) each documented its arguments by one sentence pointing to the parent: "Arguments: as for
[`density`](@ref …)", naming any option it drops ("without `normalize` (the band is for the un-normalized
estimate)"). The argument names then live in one place and are not maintained in fifteen copies.

`RULE.DOC.ARGCOV` and `RULE.DOC.ARGS_FORM` admit no delegation: the wiring reported 30 findings for it, and
SS met them by repeating each list (commit `d5c06e5`). The owner wants the delegation kept, as a form the
rulebook sanctions rather than a deviation.

Asked: a delegation form both rules accept, under conditions a program can check —

1. the docstring says "as for" followed by an `@ref` link to a function of the same package, the
   delegate;
2. the delegate's own docstring documents its arguments in full under the rules as they stand (no chains
   of delegation);
3. every argument of the delegating signature, positional and keyword, is an argument of the delegate's
   signature;
4. every argument of the delegate that the delegating signature lacks is named, as code, in the
   delegating sentence.

Condition 4 keeps the one fact a pointer would otherwise hide — what the reader does not take — and
conditions 2 and 3 keep a pointer from covering an undocumented or drifted argument. The final shape is
the tool's and the owner's; the fleet gains the same option.

When it lands, SS restores the fifteen pointers and moves to the release in one commit, and deletes this
file.
