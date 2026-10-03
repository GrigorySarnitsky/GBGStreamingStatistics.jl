<!-- gbg:managed begin v2.0.0-rc.57 — rewritten by `gbg sync`; do not edit -->
# Publication — `RULE.LEAK.*`

What a repository may not publish when it is made public: no credential, nothing personal, nothing never meant for a reader. 9 rules, written from the rule values themselves,
so a rule reworded here was reworded in the program that checks it.

### RULE.LEAK.CREDENTIAL
No credential — a token, a key, a password, a connection string with one in it — is in any commit of a package that is to be published, in a commit's metadata, or on its GitHub side.

A repository made public publishes every commit it ever held, on every branch and tag, deleted files included; a later commit that removes a thing does not unpublish it. A credential found must be revoked, not only removed.

Covered by the session, as an instruction: followed while working, checked by nothing. Nothing checks it: an instruction to the session, never to commit one; `gbg leak` decides it with Gitleaks' own credential rules, over every commit, the commit metadata and the GitHub side, which no program here restates

### RULE.LEAK.EMAIL
No email address is published but those the package's config approves in `leak_emails`, GitHub's own and no-reply addresses, the commit trailer's, and those in a domain reserved for documentation.

A repository made public publishes every commit it ever held, on every branch and tag, deleted files included; a later commit that removes a thing does not unpublish it. An address in a commit's author field is published with the commit; only a rewrite of the history removes it, which is the owner's to make.

Covered by a program.

### RULE.LEAK.HOME_PATH
No home-directory path — `/home/<name>/`, `/Users/<name>/`, `C:\Users\<name>` — is published, but a hosted runner's.

A repository made public publishes every commit it ever held, on every branch and tag, deleted files included; a later commit that removes a thing does not unpublish it. A home path names the person and the machine, and a credential scanner allows it on purpose.

Covered by a program.

### RULE.LEAK.MANIFEST_PATH
No `Manifest.toml` or `Project.toml` in any commit names an absolute path.

A repository made public publishes every commit it ever held, on every branch and tag, deleted files included; a later commit that removes a thing does not unpublish it. `RULE.GIT.LOCAL_PATH` holds the tree as it is; a path committed once and removed since is still in the history.

Covered by a program.

### RULE.LEAK.NEVER_PUBLISHED
No file that is never published — an environment file, a private key or certificate, a credentials file, a session transcript, a local settings file — is in any commit.

A repository made public publishes every commit it ever held, on every branch and tag, deleted files included; a later commit that removes a thing does not unpublish it. Such a file is judged by its name and kind, whatever it holds: its holding nothing today says nothing of the commit that held something.

Covered by a program.

### RULE.LEAK.SCRATCH_PATH
No path into an assistant's scratch directory, `/tmp/claude-<uid>/…`, is published.

A repository made public publishes every commit it ever held, on every branch and tag, deleted files included; a later commit that removes a thing does not unpublish it. A scratch path carries the user's id and the name of the session's project, and points at nothing a reader can open.

Covered by a program.

### RULE.LEAK.SESSION_LINK
No link to an assistant's session is published.

A repository made public publishes every commit it ever held, on every branch and tag, deleted files included; a later commit that removes a thing does not unpublish it. A session link opens a conversation that was never meant to be public, with whatever it held.

Covered by a program.

### RULE.LEAK.SYNCED_PATH
No path inside a synced folder — Dropbox, OneDrive, Google Drive, iCloud Drive, Nextcloud — is published.

A repository made public publishes every commit it ever held, on every branch and tag, deleted files included; a later commit that removes a thing does not unpublish it. A synced folder's path tells where the owner keeps his work, and which service holds it.

Covered by a program.

### RULE.LEAK.TERM
No word this machine's `~/.gbg/leak-terms` lists — the machines' names, or anything else the owner keeps private — is published, in any case.

A repository made public publishes every commit it ever held, on every branch and tag, deleted files included; a later commit that removes a thing does not unpublish it. What is private to one owner no generic scanner knows, and a list of it committed with a package would publish it: it is kept per machine, outside every repository.

Covered by a program.

<!-- gbg:managed end -->
