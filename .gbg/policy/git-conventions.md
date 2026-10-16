<!-- gbg:managed begin v2.0.0-rc.76 — rewritten by `gbg sync`; do not edit -->
# Git conventions — `RULE.GIT.*`

What a commit carries, where it goes, and what is never committed. 11 rules, written from the rule values themselves,
so a rule reworded here was reworded in the program that checks it.

### RULE.GIT.BRANCH
The working branch is the one `branch` names, `master` unless the package sets it; a rename is the package's decision.

A rename is a decision about the package, never bundled with a tooling adoption.

Covered by a program in the end-of-turn hook, which reads the branch and the commits not yet pushed — weaker than a build check, since a hook can be skipped and a failing test cannot.
This package may set `branch`; unset, it is `master`.

### RULE.GIT.CI_MANUAL
CI is manual-only, `workflow_dispatch` and release tags; the matrix runs locally before a commit; docs deploy on push.

Before committing a code change, run the gate locally; a documentation-only change needs a load and the docs build. The docs build is therefore the per-push enforcement point, and the rules that read the built site run inside it.

Covered by a program.

### RULE.GIT.LICENSE
A package carries the MIT licence in one file at its root, `LICENSE` or `LICENSE.md`: a copyright line naming the year and the holder, then the licence's text word for word.

A repository published without a licence can be read but not legally reused. MIT is what Julia, its standard libraries and nearly every package of the General registry carry, and what each package of the fleet depends on, so it adds nothing a user must check (the owner's ruling, 2026-10-07). The text is compared word for word because an edited licence is a different licence.

Covered by a program.

### RULE.GIT.LOCAL_PATH
No committed environment file, a `Project.toml`'s `[sources]` or a `Manifest.toml`, names a path outside the repository.

A path outside the repository exists on the machine that wrote it and nowhere else: GBGPeriodicFields committed docs and test projects whose `[sources]` named this package's Dropbox checkout, and its docs build failed on the runner while every local run passed (2026-09-18).

Covered by a program.

### RULE.GIT.MANIFEST_TRACKED
Sub-environment manifests are gitignored; a tracked root manifest is a declared deviation, a switch naming the file.

A stale sub-environment manifest is resolved locally when precompilation breaks; CI resolves its own.

Covered by a program.

### RULE.GIT.NO_SHIMS
The API is not frozen; no backward-compatibility shims.

Errors found later are fixed one by one; architectural soundness is judged on the full picture first.

Covered by the session, as an instruction: followed while working, checked by nothing. Nothing checks it: an instruction to the session: what counts as a shim is judgement; followed while working, checked by nothing

### RULE.GIT.PUSH
Commit ⇒ push, directly on the configured branch; no feature branches.

The repository is private and synced across machines through GitHub, so work is committed on the default branch and pushed in the same step; GitHub is the cross-machine sync path and the docs site deploys on push. Stage with `git add -u` plus every new file by name, never `-A`.

Covered by a program in the end-of-turn hook, which reads the branch and the commits not yet pushed — weaker than a build check, since a hook can be skipped and a failing test cannot.

### RULE.GIT.RECORDS
A record the tool keeps for the owner — `.gbg/questions.json`, and a closed run under `.gbg/runs/` — is committed, never left untracked or modified; an open run's is committed at its close.

It holds the owner's rulings and the runs' findings. Uncommitted, they exist on one machine only, and a session on the other neither sees them nor knows that a question's id is taken. No verb commits it: `gbg answer` and the review's skill say so, and the session commits.

Covered by a program in the end-of-turn hook, which reads the branch and the commits not yet pushed — weaker than a build check, since a hook can be skipped and a failing test cannot.

### RULE.GIT.REQUESTS
A cross-package request is a tracked file .gbg/requests/to-<PKG>-<date>-<topic>.md in the requester, whose first paragraph begins `Status:` and says when the owner approved sending it, as `approved by the owner on <date>`; a concluded request is deleted.

The requests directory holds what is still awaiting the other package and nothing else.

Covered by a program.

### RULE.GIT.TRACKED
A package's git tracks only the package and the tool's files: `src/`, `ext/`, `test/`, of `docs/` what its build and site take part in (the pages under `docs/src/`, the build's script, `docs/Project.toml`, `docs/prose_shape.csv`, and each document declared in `copies` or in a rationale file set the config writes out), `examples/`, `benchmark/`, `.github/`, `Project.toml`, `README.md`, a `LICENSE`, `CHANGELOG.md`, `.gitignore` and `CLAUDE.md`; `.gbg/config.jl`, `.gbg/policy/`, `.gbg/requests/`, `.gbg/questions.json`, `.gbg/settled.jl`, `.gbg/runs/` and `.gbg/leak-accepted.txt`; `.claude/settings.json`, `.claude/hooks/`, `.claude/skills/` and `.claude/WORK_QUEUE.md`. A draft or a note between sessions lives in `inspiration/`, which git ignores; any other file is declared with its reason.

Git carried notes between machines because it was the easiest way, and a public release would have published them: one package held a design proposal, notes to other sessions and old review rounds. `inspiration/` reaches both machines through Dropbox, so no note needs git, and a file kept on purpose says why.

Covered by a program in the end-of-turn hook, which reads the branch and the commits not yet pushed — weaker than a build check, since a hook can be skipped and a failing test cannot.

### RULE.GIT.TRAILER
Every commit ends with the versionless 'Co-Authored-By: Anthropic Claude <noreply@anthropic.com>' trailer, and no other Co-Authored-By or Claude-Session line; never [skip ci].

The trailer names no model version on purpose (ruled 2026-09-08): the model in use changes often and the versions supersede each other, so a version in the trailer is a fact about the day rather than about the work, and it goes stale in the one place git never lets anyone correct it.

Covered by a program in the end-of-turn hook, which reads the branch and the commits not yet pushed — weaker than a build check, since a hook can be skipped and a failing test cannot.

<!-- gbg:managed end -->
