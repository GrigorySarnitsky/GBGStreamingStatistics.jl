<!-- gbg:managed begin v2.0.0-rc.69 — rewritten by `gbg sync`; do not edit -->
# /leak-scan

Nothing secret or personal is published when a repository is made public; it writes nothing but session output.

Run it when the owner asks for it, and not otherwise; it ends after one pass over every branch, commit and message.

## The steps

1. `gbg leak` in the package, `gbg leak --all` for every package of this machine's fleet, with
   `--github` to read the GitHub side too — issues, pull requests, comments, release notes, workflow
   logs — which needs the network. Where Gitleaks is missing it says so: `gbg install --gitleaks`.
2. Read what it printed, all of it: the count by rule and by file, every finding — its rule, where it
   is, a fingerprint — and what it could not read. A finding never shows the text it found, and the
   session does not go and read that text in: whatever a session reads, it sends to its model.
3. Report to the owner in the session, and stop. What to do about each finding is his.

## What it may not do

- It writes no tracked file: every writer of a record refuses text a leak scan printed.
- It never removes a leak, from a file or from the history: a rewrite is the owner's to make, and a
  credential found is revoked, not only removed.
- It never switches a `RULE.LEAK.*` rule off, nor approves an address in `leak_emails`: the owner does.

<!-- gbg:managed end -->
