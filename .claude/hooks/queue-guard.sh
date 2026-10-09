#!/usr/bin/env bash
# STOP HOOK — refuses to let the session end while committed work is open.
# Written by GBGCodingTools v2.0.0-rc.69; do not edit.
#
# CONTRACT. `.claude/WORK_QUEUE.md` lists only what is still to do, one line per item, `- [ ]`. A
# finished item is DELETED, in the commit that finishes it: the deletion shows in the diff as plainly
# as a tick would, and git keeps the history. While any item is open, or any done line `- [x]` is
# left, this exits 2, which Claude Code treats as "do not stop", feeding stderr back as the reason. A
# blocked item stays open, its text beginning `BLOCKED:` and the reason, and is not counted, so the block
# stays in sight instead of among the done ones. Ticked lines kept as history grew GBGPeriodicFields'
# queue to 206 lines with nothing open, and nothing removed them (the owner, 2026-10-08).
#
# WHY IT IS A HOOK AND NOT A RULE. Told "continue with everything, don't stop in the middle", a
# session finished three of five queued items and reported status with the other two labelled
# "remaining queue". Stopping to narrate IS stopping. Every earlier defence was a sentence in
# CLAUDE.md, and a rule written down reads as compliance while nothing checks it.
#
# THE LOOP GUARD IS LOAD-BEARING. `stop_hook_active` is true when the session is already being held
# open by this hook; blocking again there would trap it forever. It is the difference between a
# harness and a cage.
#
# AN OPEN RUN HOLDS THE SESSION TOO, and for the reason this hook exists at all: a run is walked
# stage by stage BY THE SESSION — Julia cannot launch an agent — so a session that stops in the middle
# has spent its agents for nothing and left a run that blocks the next one. The queue cannot catch
# that, because a run writes its queue lines at CLOSE, which is the point it has not reached. Found
# when the first real review stopped between stages to report (GBGPeriodicFields, 2026-09-22).
#
# ABANDONING IS NOT THE WAY OUT and is not offered as one: it DROPS the findings. A session that
# cannot go on leaves the run OPEN and says where it stands; the next session continues it, and a
# stage whose agents already answered is not relaunched. Abandoning is the owner's call, for a run he
# does not want continued.
#
# THE UNCOUNTED HEADINGS. Lines under "## From the last review", and under "## Code the owner asked
# for" — his list of code to change, which no session is held to — are NOT counted, up to the next
# heading: a run's leftovers are the owner's list to act on when he chooses, and holding a session to
# them would be a review setting a session to work nobody asked for. A section the owner adds below
# it is his, and counted.
#
# EVERY FILE NEWLY IN GIT IS NAMED, and nothing is stopped for it: the files the index holds that it
# did not at the last turn's end, against a list kept in the ignored `.gbg/scratch/`. Shown to the
# owner whether the turn ends or is held, and to the session beside the reason when it is held. A package's git once carried a
# design proposal, notes to other sessions and old review rounds, and no commit of them was ever
# remarked on (GBGPeriodicFields' request, 2026-10-05); `RULE.GIT.TRACKED` judges them, this only
# makes each one visible.
set -uo pipefail
payload=$(cat)
gbg_root="${CLAUDE_PROJECT_DIR:-.}"
note=""; plain=""; news=0
gbg_json() {
  local s=$1 out="" c code i bs dq
  printf -v bs '\134'; printf -v dq '\042'
  for ((i = 0; i < ${#s}; i++)); do
    c=${s:i:1}
    if [[ $c == "$dq" || $c == "$bs" ]]; then
      out+="$bs$c"
    else
      printf -v code '%d' "'$c"
      (( code >= 0 && code < 32 )) && printf -v c '\\u%04x' "$code"
      out+=$c
    fi
  done
  printf '%s' "$out"
}
if git -C "$gbg_root" rev-parse --git-dir >/dev/null 2>&1; then
  seen="$gbg_root/.gbg/scratch/tracked-files"
  mkdir -p "$gbg_root/.gbg/scratch"
  git -C "$gbg_root" ls-files -z 2>/dev/null | LC_ALL=C sort -z > "$seen.now"
  if [[ -f "$seen" ]]; then
    while IFS= read -r -d '' f; do
      news=$((news + 1))
      if (( news <= 20 )); then
        note+="${note:+, }$(gbg_json "$f")"; plain+="${plain:+, }$f"
      fi
    done < <(LC_ALL=C comm -z -13 "$seen" "$seen.now" 2>/dev/null)
  fi
  mv -f "$seen.now" "$seen"
fi
if (( news > 20 )); then note+=", and $((news - 20)) more"; plain+=", and $((news - 20)) more"; fi
(( news > 0 )) && { note="new in git: $note"; plain="new in git: $plain"; }
gbg_note() { [[ -n "$note" ]] && printf '{"systemMessage": "%s"}\n' "$note"; return 0; }
finish() { gbg_note; exit 0; }
held() { [[ -n "$plain" ]] && printf '%s\n' "$plain" >&2; gbg_note; exit 2; }

if grep -q '"stop_hook_active"[[:space:]]*:[[:space:]]*true' <<<"$payload"; then
  finish
fi

runs="${CLAUDE_PROJECT_DIR:-.}/.gbg/runs"
if [[ -d "$runs" ]]; then
  for d in "$runs"/*/; do
    [[ -f "$d/stage" && ! -f "$d/run.json" ]] || continue
    id=$(basename "$d")
    {
      echo "A RUN IS OPEN — $id, at stage $(cat "$d/stage" 2>/dev/null || echo '?')."
      echo
      echo "A run is walked stage by stage by YOU: the verb writes each stage's arguments and exits,"
      echo "because Julia cannot launch an agent. Stopping in the middle spends its agents for nothing"
      echo "and leaves a run that blocks the next one."
      echo
      echo "  Continue it:  gbg run next   (or collect, or close — whichever it is waiting for)"
      echo "  Or stop here: leave it OPEN and say which run and which stage, so the next session"
      echo "                picks it up. Nothing is lost; a stage whose agents answered is not relaunched."
      echo
      echo "Do NOT abandon it to end the turn: 'gbg run abandon' DROPS its findings. That is the"
      echo "owner's call, for a run he does not want continued."
    } >&2
    held
  done
fi

# THE GIT RULES. RULE.GIT.BRANCH, RULE.GIT.PUSH and RULE.GIT.TRAILER ask about the branch and the
# commits — what no build or test may read, since a fresh clone and CI must answer alike — so they
# are checked HERE, at the end of the turn, on what is not yet pushed. A hook check is weaker than
# a build check: a hook can be skipped and a failing test cannot, and the rulebook says so beside
# each. Every function below is the rule's own lines, spliced from the rule value, so `gbg check`
# and this script are one program and cannot disagree.
    gbg_rule_RULE_GIT_BRANCH() {
      b=$(git rev-parse --abbrev-ref HEAD 2>/dev/null) || return 0
      [ "$b" = "$GBG_BRANCH" ] || echo "on branch '$b'; the working branch is '$GBG_BRANCH'"
      return 0
    }
    gbg_rule_RULE_GIT_PUSH() {
      git rev-parse --verify -q HEAD >/dev/null 2>&1 || return 0
      if ! git rev-parse --abbrev-ref '@{upstream}' >/dev/null 2>&1; then
        n=$(git rev-list --count HEAD 2>/dev/null) || return 0
        [ "$n" = "0" ] || echo "no upstream; $n commit(s) never pushed: commit ⇒ push"
        return 0
      fi
      n=$(git rev-list --count '@{upstream}..HEAD' 2>/dev/null) || return 0
      [ "$n" = "0" ] || echo "$n commit(s) not pushed: commit ⇒ push"
      return 0
    }
    gbg_rule_RULE_GIT_TRAILER() {
      git rev-parse --verify -q HEAD >/dev/null 2>&1 || return 0
      range='@{upstream}..HEAD'
      git rev-parse --abbrev-ref '@{upstream}' >/dev/null 2>&1 || range=HEAD
      for c in $(git rev-list "$range" 2>/dev/null); do
        msg=$(git log -1 --format=%B "$c")
        printf '%s\n' "$msg" | grep -Fxq 'Co-Authored-By: Anthropic Claude <noreply@anthropic.com>' || echo "$(git log -1 --format=%h "$c") does not carry the trailer 'Co-Authored-By: Anthropic Claude <noreply@anthropic.com>'"
        printf '%s\n' "$msg" | grep -Eq '^Co-Authored-By: .*<noreply@anthropic\.com>' && ! printf '%s\n' "$msg" | grep -Eq '^Co-Authored-By: Anthropic Claude <noreply@anthropic\.com>$' && echo "$(git log -1 --format=%h "$c") carries a versioned Co-Authored-By trailer"
        printf '%s\n' "$msg" | grep -Eq '^Claude-Session:' && echo "$(git log -1 --format=%h "$c") carries a Claude-Session line"
        printf '%s\n' "$msg" | grep -Fq '[skip ci]' && echo "$(git log -1 --format=%h "$c") says [skip ci]"
      done
      return 0
    }
    gbg_rule_RULE_GIT_RECORDS() {
      git rev-parse --git-dir >/dev/null 2>&1 || return 0
      seen=""
      while IFS= read -r line; do
        p=${line:3}
        case "$p" in
          .gbg/runs/*) d=${p#.gbg/runs/}; d=".gbg/runs/${d%%/*}"
                       [ -f "$d/stage" ] && [ ! -f "$d/run.json" ] && continue
                       case " $seen " in *" $d "*) continue ;; esac
                       seen="$seen $d"; p="$d/" ;;
        esac
        echo "$p is a record the tool keeps for the owner, not committed: commit it"
      done < <(git status --porcelain --untracked-files=all -- .gbg/questions.json .gbg/runs 2>/dev/null)
      return 0
    }

GBG_BRANCH="master"
git_said=""
for r in RULE.GIT.BRANCH RULE.GIT.PUSH RULE.GIT.TRAILER RULE.GIT.RECORDS; do
  fn="gbg_rule_${r//./_}"
  out=$(cd "${CLAUDE_PROJECT_DIR:-.}" && "$fn" 2>/dev/null) || out=""
  [[ -n "$out" ]] && git_said+="$(printf '%s\n' "$out" | sed "s/^/  $r  /")"$'\n'
done
if [[ -n "$git_said" ]]; then
  {
    echo "THE GIT RULES ARE NOT MET — the turn does not end on a branch or a commit the rulebook refuses:"
    echo
    printf '%s' "$git_said"
    echo
    echo "Push what is committed (commit ⇒ push), on the working branch, with the trailer on every"
    echo "commit — or say why you cannot, so the owner decides. A hook can be skipped; do not."
  } >&2
  held
fi

queue="${CLAUDE_PROJECT_DIR:-.}/.claude/WORK_QUEUE.md"
[[ -f "$queue" ]] || finish

# The leftovers section is the EXACT heading line, up to the next heading of ANY level, and the
# lines printed carry their number in the file (found in review, 2026-09-23: a section of the
# owner's under a `#` heading, or one whose name began with the leftovers', went uncounted).
# An item is any list bullet (`-`, `*`, `+`, `1.`) with a box; a line inside a fenced block is an example, not
# an item; an open item is BLOCKED only when its text BEGINS `BLOCKED:` and gives a reason (found in review).
open_items=$(awk '/^[[:space:]]*(```|~~~)/ {f = !f; next} f {next}
                    { h = $0; sub(/[ 	]+$/, "", h) }
                    h == "## From the last review" || h == "## Code the owner asked for" {skip=1; next}
                    /^#/ {skip=0}
                    !skip && /^[[:space:]]*([-*+]|[0-9]+[.)])[[:space:]]*\[[[:space:]]*\]/ &&
                    !/^[[:space:]]*([-*+]|[0-9]+[.)])[[:space:]]*\[[[:space:]]*\][[:space:]]+BLOCKED:[[:space:]]*[^[:space:]]/ {print NR ":" $0}' "$queue")
# A done line in ANY section, the owner's and the last review's included: the queue keeps no history.
done_items=$(awk '/^[[:space:]]*(```|~~~)/ {f = !f; next} f {next}
                    /^[[:space:]]*([-*+]|[0-9]+[.)])[[:space:]]*\[[[:space:]]*[xX][[:space:]]*\]/ {print NR ":" $0}' "$queue")
[[ -z "$open_items" && -z "$done_items" ]] && finish

{
  if [[ -n "$open_items" ]]; then
    echo "WORK QUEUE IS NOT EMPTY — you committed to finishing these before reporting back:"
    echo
    sed 's/^/  /' <<<"$open_items"
    echo
    echo "Continue working. DELETE each item's line from .claude/WORK_QUEUE.md, in the commit that"
    echo "finishes it, ONLY when it is actually done and verified — an item removed to end the turn is"
    echo "the defect this hook exists to prevent. If an item is genuinely blocked, keep it '- [ ]' and"
    echo "begin its text with 'BLOCKED:' and the reason; such a line is not counted, and stays in sight."
  fi
  if [[ -n "$done_items" ]]; then
    [[ -n "$open_items" ]] && echo
    echo "WORK QUEUE HOLDS DONE ITEMS — it lists only what is still to do. Delete these lines (git keeps"
    echo "the history), in any section, the owner's included:"
    echo
    sed 's/^/  /' <<<"$done_items"
  fi
} >&2
held
