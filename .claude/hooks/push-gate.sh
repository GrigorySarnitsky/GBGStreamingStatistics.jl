#!/usr/bin/env bash
# PUSH HOOK — the docs build runs before every push of this package, and a build that fails stops it.
# Written by GBGCodingTools v2.0.0-rc.83; do not edit.
#
# CONTRACT. Claude Code runs this before every shell command (PreToolUse), with the call as JSON on
# standard input. A command that pushes this package — `git push` in any spelling, or `gbg release`,
# which pushes a tag — first runs the package's docs build, in the package root; exit 2 stops the push
# and hands stderr back as the reason, exit 0 lets it through.
#
# WHY. The owner's ruling of 2026-09-29: every push does the docs build. A session ran only
# `Pkg.test()` before pushing, so a broken cross-reference reached GitHub sixteen times in one package.
#
# WHAT IS READ. Only the command and the directory it runs in, decoded from the JSON: the call's
# description is prose, and matching the whole payload ran a build for any command that merely named a
# push. A quoted string and a heredoc body are data, not commands, unless a shell is told to run them.
# A push from another repository is that repository's to judge, not this one's; and a push while
# tracked files are uncommitted is stopped, since the build reads the tree and the push sends the commit
# (all four found in review, 2026-09-29).
set -uo pipefail
payload=$(cat)

# The string value of a top-level key of the call, its escapes decoded.
field() {
  printf '%s' "$payload" | awk -v key="\"$1\"" '
    { s = s $0 "\n" }
    END {
      i = index(s, key); if (i == 0) exit
      s = substr(s, i + length(key)); sub(/^[ \t\n]*:[ \t\n]*"/, "", s)
      out = ""
      for (k = 1; k <= length(s); k++) {
        c = substr(s, k, 1)
        if (c == "\\") { k++; e = substr(s, k, 1)
          if (e == "n") out = out "\n"; else if (e == "t") out = out "\t"; else out = out e
          continue }
        if (c == "\"") break
        out = out c
      }
      printf "%s", out
    }'
}
cmd=$(field command)
here=$(field cwd)
[ -n "$here" ] || here="${CLAUDE_PROJECT_DIR:-.}"

# Whether the command pushes, and every directory it moves to (`cd`, `git -C`) on the way.
verdict=$(printf '%s\n' "$cmd" | awk '
  tag != "" { t = $0; gsub(/^[ \t]+|[ \t]+$/, "", t); if (t == tag) tag = ""; next }
  {
    raw = $0; line = $0
    if (match(line, /<<-?[ \t]*["\047]?[A-Za-z_][A-Za-z0-9_]*["\047]?/)) {
      h = substr(line, RSTART, RLENGTH); gsub(/^<<-?[ \t]*["\047]?|["\047]?$/, "", h); tag = h
    }
    if (line ~ /(^|[ \t;&|(])(ba|z)?sh[ \t]+-c[ \t]/) gsub(/["\047]/, " ", line)     # the shell runs the quoted text
    else { gsub(/"[^"]*"/, "\"\"", line); gsub(/\047[^\047]*\047/, "\047\047", line) }
    s = raw
    while (match(s, /(^|[ \t;&|(])cd[ \t]+("[^"]*"|\047[^\047]*\047|[^ \t;&|)]+)/)) {
      d = substr(s, RSTART, RLENGTH); sub(/^[ \t;&|(]*cd[ \t]+/, "", d); print "DIR\t" d; s = substr(s, RSTART + RLENGTH)
    }
    s = raw
    while (match(s, /git[ \t]+-C[ \t]+("[^"]*"|\047[^\047]*\047|[^ \t;&|)]+)/)) {
      d = substr(s, RSTART, RLENGTH); sub(/^git[ \t]+-C[ \t]+/, "", d); print "DIR\t" d; s = substr(s, RSTART + RLENGTH)
    }
    if (!pushes && line ~ /(^|[;&|( \t])git([ \t]+(-C[ \t]+[^ \t;&|]+|-c[ \t]+[^ \t;&|]+|--[^ \t;&|]+))*[ \t]+commit([^A-Za-z0-9_-]|$)/) commits = 1
    if (line ~ /(^|[;&|( \t])git([ \t]+(-C[ \t]+[^ \t;&|]+|-c[ \t]+[^ \t;&|]+|--[^ \t;&|]+))*[ \t]+add[ \t]+[^;&|]*[^ \t;&|]/) {
      a = line; sub(/.*git([ \t]+(-C[ \t]+[^ \t;&|]+|-c[ \t]+[^ \t;&|]+|--[^ \t;&|]+))*[ \t]+add[ \t]+/, "", a); sub(/[;&|].*/, "", a)
      gsub(/(^|[ \t])(-u|--update)([ \t]|$)/, " ", a); if (a ~ /[^ \t]/) adds = 1
    }
    if (line ~ /(^|[;&|( \t])git([ \t]+(-C[ \t]+[^ \t;&|]+|-c[ \t]+[^ \t;&|]+|--[^ \t;&|]+))*[ \t]+push([^A-Za-z0-9_-]|$)/ ||
        line ~ /(^|[;&|( \t])gbg(dev)?([ \t]+--dev)?[ \t]+release([ \t;&|]|$)/) pushes = 1
  }
  END { if (pushes) print "PUSH"; if (pushes && commits) print "COMMITS"; if (pushes && commits && adds) print "ADDS" }')
case "$verdict" in *PUSH*) ;; *) exit 0 ;; esac

# Which repository is pushed: the call's directory, moved by each `cd` and `-C` in turn.
root=$(cd "${CLAUDE_PROJECT_DIR:-.}" 2>/dev/null && git rev-parse --show-toplevel 2>/dev/null) || exit 0
where="$here"
while IFS="$(printf '\t')" read -r kind d; do
  [ "$kind" = "DIR" ] || continue
  d=${d#\"}; d=${d%\"}; d=${d#\'}; d=${d%\'}
  case "$d" in
    '$CLAUDE_PROJECT_DIR'|'${CLAUDE_PROJECT_DIR}') d="${CLAUDE_PROJECT_DIR:-.}" ;;
    '~'|'~/'*) d="$HOME${d#\~}" ;;
    '$'*) d="" ;;
  esac
  [ -n "$d" ] || { where=""; break; }
  case "$d" in /*) where="$d" ;; *) where="$where/$d" ;; esac
done <<EOF
$verdict
EOF
if [ -n "$where" ]; then
  target=$(cd "$where" 2>/dev/null && git rev-parse --show-toplevel 2>/dev/null) || target=""
  [ -n "$target" ] && [ "$target" != "$root" ] && exit 0      # another repository's push: not this one's
fi

# A command that commits before it pushes is checked on the tree it is about to commit; any other push
# while tracked files are uncommitted is stopped.
case "$verdict" in *COMMITS*) committing=1 ;; *) committing=0 ;; esac
if [ "$committing" = 0 ] && [ -n "$(git -C "$root" status --porcelain --untracked-files=no 2>/dev/null)" ]; then
  echo "PUSH STOPPED: tracked files are uncommitted. The docs build reads the tree and the push sends the commit; commit first, then push." >&2
  exit 2
fi
cd "$root" || exit 0
# A session link stops the push (the owner, 2026-10-09): `RULE.GIT.TRAILER` forbids the harness's
# `Claude-Session:` line, and the end-of-turn hook that reported it holds a session only once, so a
# package published one with nearly every commit. Read first, since it needs no Julia: every commit of
# every local branch and tag that no remote holds, and the commit this command makes — its `-m`,
# `--trailer` and `-F` file (found in review, 2026-10-09).
links=$(git rev-list --branches --tags --not --remotes 2>/dev/null | while read -r c; do
  git log -1 --format=%B "$c" | grep -q '^Claude-Session:' && git log -1 --format='%h %s' "$c"
done)
if [ "$committing" = 1 ]; then
  made=$(printf '%s\n' "$cmd")
  for f in $(printf '%s\n' "$cmd" | grep -oE '(^|[[:space:]])(-[A-Za-z]*F|--file)[[:space:]=]+[^[:space:];&|]+' | sed -E 's/^[[:space:]]*(-[A-Za-z]*F|--file)[[:space:]=]+//'); do
    f=${f#\"}; f=${f%\"}; f=${f#\'}; f=${f%\'}
    case "$f" in /*) ;; *) f="${where:-$root}/$f" ;; esac
    [ -f "$f" ] && made=$(printf '%s\n%s' "$made" "$(cat "$f")")
  done
  printf '%s\n' "$made" | grep -qE '(^|["='"'"'])[[:space:]]*Claude-Session:' &&
    links=$(printf '%s\n%s' "$links" "the commit this command makes")
fi
links=$(printf '%s\n' "$links" | sed '/^$/d')
if [ -n "$links" ]; then
  {
    echo "PUSH STOPPED: a commit carries a Claude-Session line, a link to an assistant's session, which"
    echo "RULE.GIT.TRAILER forbids by the owner's ruling, whatever the harness asks, and which would be published:"
    printf '%s\n' "$links" | sed 's/^/  /'
    echo "Remove the line from each: 'git commit --amend' for the last commit, or for every commit no remote holds"
    echo "  git filter-branch -f --msg-filter \"sed '/^Claude-Session:/d'\" -- --branches --tags --not --remotes"
    echo "and end every commit with the one trailer 'Co-Authored-By: Anthropic Claude <noreply@anthropic.com>'."
  } >&2
  exit 2
fi
# What git tracks is judged before the build, by `RULE.GIT.TRACKED` under the package's own switches:
# a stray file stops the push (the owner, 2026-10-05), since publishing is when it does harm.
log=$(mktemp)
# a command that stages new files, commits and pushes stages them only after this runs: files git
# neither tracks nor ignores count then too, and only then — `git add -u` stages none (review, 2026-10-05)
case "$verdict" in *ADDS*) adding=1 ;; *) adding=0 ;; esac
if ! strays=$(GBG_UNTRACKED="$adding" julia --project=docs --startup-file=no -e 'using GBGCodingTools; const G = GBGCodingTools; foreach(println, G.Engine.Checks.tracked_strays(pwd(), G.ConfigIO.read_config(pwd()); untracked = get(ENV, "GBG_UNTRACKED", "0") == "1"))' 2>"$log"); then
  {
    echo "PUSH STOPPED: what git tracks could not be checked (RULE.GIT.TRACKED):"
    printf '%s\n' "$strays"; grep -E 'ERROR|Error' "$log"
    echo "The whole output: $log"
  } >&2
  exit 2
fi
if [ -n "$strays" ]; then
  {
    echo "PUSH STOPPED: files the rulebook does not allow in git (RULE.GIT.TRACKED):"
    printf '%s\n' "$strays" | sed 's/^/  /'
    echo "Move a draft or a note to inspiration/ (git ignores it) and git rm --cached it, delete it, or ask"
    echo "the owner for a switch naming it, with its reason, in .gbg/config.jl; then commit, and push in a"
    echo "separate command — this check reads the index as it stands before the command runs."
  } >&2
  rm -f "$log"
  exit 2
fi
if julia --project=docs docs/make.jl >"$log" 2>&1; then
  rm -f "$log"
  exit 0
fi
{
  echo "PUSH STOPPED: the docs build failed, and every push does the docs build (the owner, 2026-09-29)."
  echo "What it reported, every finding and error line:"
  grep -E '^  (RULE\.|WIRING|CONFIG)|ERROR|Error' "$log"
  echo "The whole output: $log"
} >&2
exit 2
