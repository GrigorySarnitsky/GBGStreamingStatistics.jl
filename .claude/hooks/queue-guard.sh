#!/usr/bin/env bash
# STOP HOOK — refuses to let the session end while committed work is open.
# Written by GBGCodingTools v2.0.0-rc.58; do not edit.
#
# CONTRACT. `.claude/WORK_QUEUE.md` holds one line per item, `- [ ]` open, `- [x]` done. While any
# item is open this exits 2, which Claude Code treats as "do not stop", feeding stderr back as the
# reason. Closing an item is an edit that shows in the diff, so the queue cannot be silently emptied.
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
set -uo pipefail
payload=$(cat)
if grep -q '"stop_hook_active"[[:space:]]*:[[:space:]]*true' <<<"$payload"; then
  exit 0
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
    exit 2
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
  exit 2
fi

queue="${CLAUDE_PROJECT_DIR:-.}/.claude/WORK_QUEUE.md"
[[ -f "$queue" ]] || exit 0

# The leftovers section is the EXACT heading line, up to the next heading of ANY level, and the
# lines printed carry their number in the file (found in review, 2026-09-23: a section of the
# owner's under a `#` heading, or one whose name began with the leftovers', went uncounted).
open_items=$(awk '{ h = $0; sub(/[ 	]+$/, "", h) }
                    h == "## From the last review" || h == "## Code the owner asked for" {skip=1; next}
                    /^#/ {skip=0}
                    !skip && /^[[:space:]]*- \[ \]/ {print NR ":" $0}' "$queue")
[[ -z "$open_items" ]] && exit 0

{
  echo "WORK QUEUE IS NOT EMPTY — you committed to finishing these before reporting back:"
  echo
  sed 's/^/  /' <<<"$open_items"
  echo
  echo "Continue working. Close each item in .claude/WORK_QUEUE.md by changing '- [ ]' to '- [x]'"
  echo "ONLY when it is actually done and verified — an item closed to end the turn is the defect"
  echo "this hook exists to prevent. If an item is genuinely blocked, mark it '- [x]' and record the"
  echo "blocker on the same line, so the block is visible in the diff rather than silent."
} >&2
exit 2
