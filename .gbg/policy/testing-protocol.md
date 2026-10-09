<!-- gbg:managed begin v2.0.0-rc.69 — rewritten by `gbg sync`; do not edit -->
# Testing protocol — `RULE.TEST.*`

What a test is for, where its oracle comes from, and when a thing is done. 12 rules, written from the rule values themselves,
so a rule reworded here was reworded in the program that checks it.

### RULE.TEST.ALLOC_ZERO
In a package that declares the allocation audit, every in-place ! operator is allocation-free and in the audit; a deliberate materializer gets a bounded guard; where a bound covers another nondeterminism a configuration pins the true zero.

The audit runs on the CPU, where a materialized temporary is visible; the completeness of its rows is derived from the export list, never remembered.

Covered by the session, as an instruction: followed while working, checked by nothing. Nothing checks it: an instruction to the session of a package that has declared this contract by calling `TestHarness.allocation_completeness` from its suite, where that suite checks it; a package that has not declared it is not bound, and nothing checks it here

### RULE.TEST.BOTH_DIRECTIONS
Test the check in both directions; a check that cannot fail is not a check; a self-reported field is acceptable only when the missing case is the failing case.

The flattering direction is the one nobody guards: a check that demotes sound answers is stopped being read, and a check nobody reads passes everything.

Covered by the session, as an instruction: followed while working, checked by nothing. Nothing checks it: an instruction to the session: followed while working, checked by nothing

### RULE.TEST.DEFINITION_OF_DONE
A feature is done when its correctness tests, performance test, four guiding docs, remaining docs, benchmarks (built AND run) and examples (run) are each verified.

'Should work' is not done: a step is done when its verification command exits 0 and the exit code was observed. Grep the docs for stale deferred, none-exist and future claims about whatever was just implemented.

Covered by the session, as an instruction: followed while working, checked by nothing. Nothing checks it: an instruction to the session: followed while working, checked by nothing

### RULE.TEST.DEVICE_COST
In a package that declares the device-cost audit, each audited device operation is timed against the memory it must move, in a ratio to a copy measured in the same run, and holds at or below the ratio frozen for it on that device and element type.

A count of reads and writes cannot see computation that costs time; only a timing against the least traffic can. The copy measured beside it keeps the ratio honest under the device's clocks and power limit; the fastest of many runs is the cost, and a frozen bound is refrozen in the commit that changes the operation, with the reason.

Covered by the session, as an instruction: followed while working, checked by nothing. Nothing checks it: an instruction to the session of a package that has declared this contract by calling `TestHarness.device_cost` from its GPU tests, where those tests check it; without a device they record a skip, never a pass, so it cannot run on CI; a package that has not declared it is not bound, and nothing checks it here

### RULE.TEST.EXAMPLES_RUN
Shipped examples are executed after any public-surface change; nothing else runs them.

An example is both a rendered page and a downloadable script; the docs build renders it without running it, and a call that no longer dispatches is found only by running it.

Covered by the session, as an instruction: followed while working, checked by nothing. Nothing checks it: an instruction to the session: followed while working, checked by nothing

### RULE.TEST.GATE_VS_SURVEY
Pkg.test() is a gate; the exhaustive suite is a survey; a specific oracle and every completeness guard go in the gate.

A gate re-checks claims already known and must be green or the change does not ship; a survey searches for combinations nobody enumerated and is read, not passed. A universal oracle over a product goes in the survey, and anything that may legitimately be non-green goes there too, because a suppressed failure reads as a pass.

Covered by the session, as an instruction: followed while working, checked by nothing. Nothing checks it: an instruction to the session: followed while working, checked by nothing

### RULE.TEST.METHOD_PARITY
In a package that declares method parity, the public method set is frozen per generic and per module; any change moves a number.

The per-module count is derived from ownership over the whole function surface, internals included, so a method defined under a qualified name cannot hide.

Covered by the session, as an instruction: followed while working, checked by nothing. Nothing checks it: an instruction to the session of a package that has declared this contract by calling `TestHarness.method_parity` from its suite, where that suite checks it; a package that has not declared it is not bound, and nothing checks it here

### RULE.TEST.MISSING_IS_FAILING
The missing case is the failing case: derived denominators, absence arms, zero count = dead.

An empty witness directory, an empty ledger, an unparsable status — each fails; 'nothing to check' is not 'all clear'.

Covered by the session, as an instruction: followed while working, checked by nothing. Nothing checks it: an instruction to the session: followed while working, checked by nothing

### RULE.TEST.MUTATION
A symmetric fixture cannot see an asymmetric defect; a self-referential assertion cannot see a shared defect; a stated convention must be falsifiable by its fixture.

Inject the defect and require the suite to catch it; a fixture that cannot reach the broken branch proves nothing about it.

Covered by the session, as an instruction: followed while working, checked by nothing. Nothing checks it: an instruction to the session: followed while working, checked by nothing

### RULE.TEST.ORACLE_INDEPENDENT
In a package that declares the oracle registry, every public value-returning operation is classified in it; relational-only coverage is a finding.

The classes are analytic (a closed form derived independently of the package), independent (a direct computation written in the test), property (an invariant that pins the operation without computing it), relational (compared only against other package operations — legitimate as a second layer, never as the only one) and structural (returns metadata or storage). A new operation fails the registry until someone classifies it, which is the moment to notice it has no independent oracle.

Covered by the session, as an instruction: followed while working, checked by nothing. Nothing checks it: an instruction to the session of a package that has declared this contract by calling `TestHarness.oracle_registry` from its suite, where that suite checks it; a package that has not declared it is not bound, and nothing checks it here

### RULE.TEST.PASS_AUDIT
In a package that declares the traffic audit, every audited operation has a frozen per-array traffic row with a named holder for any gap.

A raised count is a fusion regression unless a new holder justifies it; the frozen values are refrozen in the same commit as the change, with the reason.

Covered by the session, as an instruction: followed while working, checked by nothing. Nothing checks it: an instruction to the session of a package that has declared this contract by calling `TestHarness.traffic_completeness` from its suite, where that suite checks it; a package that has not declared it is not bound, and nothing checks it here

### RULE.TEST.SEED_PARITY
In a package that declares seed parity, one seed gives the same draws on every backend it runs: a draw the run uses as drawn is bitwise equal, one it sees through a transform agrees within a stated bound in machine epsilon, and a run from another seed is told apart on every draw.

Runs on two machines are then comparable and a regression is reproducible on either. How a package meets it is its own: drawing on the host and moving the numbers, or one counter-based generator on both, which only the comparison shows to agree.

Covered by the session, as an instruction: followed while working, checked by nothing. Nothing checks it: an instruction to the session of a package that has declared this contract by calling `TestHarness.seed_parity` from its GPU tests, where those tests check it; without a GPU they record a skip, never a pass, so it cannot run on CI; a package that has not declared it is not bound, and nothing checks it here

<!-- gbg:managed end -->
