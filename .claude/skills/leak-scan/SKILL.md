<!-- gbg:managed begin v2.0.0-rc.73 — rewritten by `gbg sync`; do not edit -->
# /leak-scan

Nothing secret or personal is published when a repository is made public; it records only accepted fingerprints.

Run it when the owner asks for it, or as soon as he considers making the repository public — before he
decides, since it says what publishing would cost — and not otherwise; it ends after one pass over every branch, commit and message.

## The steps

1. `gbg leak` in the package, `gbg leak --all` for every package of this machine's fleet, with
   `--github` to read the GitHub side too — issues, pull requests, comments, release notes, workflow
   logs — which needs the network. Where Gitleaks is missing it says so: `gbg install --gitleaks`.
2. Read what it printed, all of it: the count by rule and by file, every finding — its rule, where it
   is, a fingerprint — and what it could not read. A finding never shows the text it found, and the
   session does not go and read a credential, or a word the owner keeps private (`RULE.LEAK.TERM`),
   in: whatever a session reads, it sends to its model.
3. Sort the findings into the three tiers the scan prints, and bring the owner ONE decision per tier,
   with your recommendation — never a list of findings to judge one by one:
   - **Blocking** (`RULE.LEAK.CREDENTIAL`, `RULE.LEAK.NEVER_PUBLISHED`). Not publishable as it stands. A
     credential is revoked first, wherever it was issued; then the file or the secret leaves the
     history. The route to propose: rewrite the history once (`git filter-repo`, on a copy), removing
     every blocking file and every session link, and publish the rewritten history as a NEW public
     repository, keeping this one private as the archive; the commit map `git filter-repo` writes
     lets every commit hash the package's documents cite be updated mechanically. The alternative,
     a new repository started from today's tree, loses the public history.
   - **In the current files**, and on the GitHub side, which turns public with the repository. Cleaned.
     A path, an address or a session link is no secret: open today's file, find it, and fix it — in a
     document a documentation fix; in code it waits for the owner's answer like any code change; on
     GitHub an edit or a deletion there. An email address stays only if the owner approves it in
     `leak_emails`. A private word (`RULE.LEAK.TERM`) you do not read: name its file to the owner.
     A workflow log is the output GitHub keeps of one run, and it deletes the log by itself once
     the repository's retention period has passed (90 days unless its Actions settings say
     otherwise; `gh api repos/<owner>/<repo>/actions/permissions/artifact-and-log-retention`
     shows it). A run the scan lists as one whose log is no longer kept is one GitHub no longer
     serves — past its retention, or deleted — so it has nothing to publish and needs no check
     by the owner. A log is removed only for a finding the scan names in that run's own log, and
     only on the owner's answer; propose deleting the log alone
     (`gh api -X DELETE repos/<owner>/<repo>/actions/runs/<id>/logs`), which keeps the run's
     record, rather than the run (`gh run delete`). A log is how a failed run is diagnosed, and a
     pass that prints `INCOMPLETE`, or whose findings move from one pass to the next, is the
     scan's defect to report, not a reason to delete anything.
   - **Only in the history** (a path, a word, an address, a session link in an old commit or its
     message). A note, not a secret — a session link opens only to the owner's account — and so
     accepted once the current files are clean (the owner's three tiers, 2026-10-02; GBGCodingTools
     accepted its own on 2026-10-03): on the owner's one answer, `gbg leak --accept` records them all
     in `.gbg/leak-accepted.txt`, and you commit that file.
4. The package is ready to publish when a scan, `--github` included, prints no blocking finding, none
   in the current files, and every history note accepted. Say so, with the counts, and with what the
   scan says it could not read (the wiki, release assets and the like), which the owner checks himself;
   the runs whose log is no longer kept you give as one count, since there is nothing in them to check.
5. Say also what publishing needs beyond the scan: the package's test and docs environments fetch
   GBGCodingTools by its URL, so while GBGCodingTools is private nobody else can run them.

A session link enters the history with every commit that carries the harness's `Claude-Session:`
line. `RULE.GIT.TRAILER` forbids it, by the owner's ruling, whatever the harness says, and the push
hook refuses a push that would publish one.

## What it may not do

- It writes no tracked file but `.gbg/leak-accepted.txt`, and that only on `--accept`, which records
  fingerprints alone: every writer of a record refuses text a leak scan printed.
- It never removes a leak, from a file or from the history, on its own: a rewrite is made on the
  owner's answer, and a credential found is revoked, not only removed.
- It never switches a `RULE.LEAK.*` rule off, nor approves an address in `leak_emails`: the owner does.

<!-- gbg:managed end -->
