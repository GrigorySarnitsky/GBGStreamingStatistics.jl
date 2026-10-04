#!/usr/bin/env bash
# PUSH HOOK — the docs build runs before every push of this package, and a build that fails stops it.
# Written by GBGCodingTools v2.0.0-rc.58; do not edit.
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
    if (line ~ /(^|[;&|( \t])git([ \t]+(-C[ \t]+[^ \t;&|]+|-c[ \t]+[^ \t;&|]+|--[^ \t;&|]+))*[ \t]+push([^A-Za-z0-9_-]|$)/ ||
        line ~ /(^|[;&|( \t])gbg(dev)?([ \t]+--dev)?[ \t]+release([ \t;&|]|$)/) pushes = 1
  }
  END { if (pushes) print "PUSH"; if (pushes && commits) print "COMMITS" }')
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
log=$(mktemp)
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
