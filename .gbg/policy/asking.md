<!-- gbg:managed begin v2.0.0-rc.76 — rewritten by `gbg sync`; do not edit -->
# Asking the owner — `RULE.ASK.*`

One decision, short, in his words rather than this package's, with every identifier moved out of the sentence. 5 rules, written from the rule values themselves,
so a rule reworded here was reworded in the program that checks it.

### RULE.ASK.NO_COORDINATES
A question's decision and context name no rule, file or piece of code; those go in its items and its exhibit.

A decision carrying an identifier asks the owner to hold a coordinate in his head before he can reach the question, and the question is the part he is being asked about. A name in backticks is exempt: quoting one is how a question says which exact name it means.

Covered by a program.

### RULE.ASK.NO_HOUSE_TERMS
A question is written in the owner's words, not in this package's.

He is being asked about his own work. A question written in the tool's vocabulary is answerable only by someone who has read the tool, which makes the answer the tool's rather than his.

Covered by a program.
This package may set `house_terms`; unset, it is `["fileset", "file set", "ratchet", "anchor", "wiring check", "managed block", "pin", "switch", "rulebook", "instrument", "lane", "hat", "round", "harness", "consumer"]`.

### RULE.ASK.ONE_DECISION
A question's decision is one sentence, and it ends in a question mark.

Two sentences in a decision are two decisions, and the owner answers with one option: whichever he meant, the other is recorded as answered and was not.

Covered by a program.

### RULE.ASK.PLAIN
A question reads plainly, stands on its own, and invents no word for something that already has one.

A question is read once, by someone who was doing something else. Whether it can be answered without opening anything is a matter of meaning, which no program decides.

Covered by an agent, because this needs meaning.

### RULE.ASK.SHORT
A question's decision and its context are short enough to read in one pass; its exhibit is not measured.

What does not fit belongs in the exhibit, which holds the commits, the diffs and the definitions the owner reads to decide. A limit on the exhibit would push the evidence out of the question instead of shortening the question.

Covered by a program.
This package may set `decision_words`; unset, it is `40`.
This package may set `context_words`; unset, it is `150`.

<!-- gbg:managed end -->
