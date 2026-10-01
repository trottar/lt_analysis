# KaonLT workflow checkpoint Fix.1 — synchronize all authoritative policy surfaces

## Purpose

Repair one concrete blocker found by ChatGPT actual-diff review of:

```text
kaonlt_review(20261001-111102).diff
```

The first memory-bracketed-throughput candidate correctly updates CURRENT, MEMORY,
USER, CODEX, MAINTENANCE, LEARNINGS, roadmap status, and the new decision record,
but it leaves several authoritative startup/contract/farm surfaces on the old
recursive policy. It also leaves CODEX/MAINTENANCE requiring `--fail-on-warning`
as the ordinary task-final health gate while simultaneously declaring some
warnings nonblocking.

That contradiction would allow the old stagnation behavior to re-enter the next
session or future Codex contract.

This is a **narrow memory-policy surface repair only**. Do not reopen scientific
architecture, F.4/F.6/E.8 implementation, farm orchestration, or accepted
scientific evidence.

---

## Starting state

Required branch:

```text
test
```

Required committed HEAD:

```text
712ba32b062772d87fe44efa6865346e5c438827
```

The worktree must contain the uncommitted memory-bracketed-throughput candidate
reviewed in `kaonlt_review(20261001-111102).diff`, including exactly these
candidate paths:

```text
docs/memory/CODEX.md
docs/memory/CURRENT.md
docs/memory/LEARNINGS.md
docs/memory/MAINTENANCE.md
docs/memory/MEMORY.md
docs/memory/USER.md
docs/memory/decisions/memory-bracketed-scientific-throughput.md
docs/memory/manifest.json
docs/memory/phases/memory-bracketed-scientific-throughput-task-contract.md
docs/memory/roadmap/STATUS.md
```

plus the temporary review bundle.

Before editing, verify:

```bash
git branch --show-current
git rev-parse HEAD
git status --short --untracked-files=all
```

If the committed HEAD differs, or unrelated versionable changes are present,
STOP.

Do not discard, reset, or rewrite the already reviewed candidate. Apply this fix
on top of it.

---

## Why this repair is required

The previous candidate establishes the intended durable rule:

```text
memory checkpoint
-> substantive science / implementation
-> local deterministic test
-> ChatGPT actual-diff review
-> user commit/push
-> ChatGPT pushed-state synchronization review
-> narrow farm/runtime test when required
-> evidence review
-> memory checkpoint
```

and says only hard failures or **material** active-state/provenance ambiguity
should block progress.

However, authoritative surfaces omitted from the first contract still encode
the old behavior:

1. `docs/memory/AGENTS.md`
   - startup authority still says unresolved warnings block the next gate;
   - workflow still omits ChatGPT pushed-state synchronization.

2. `docs/memory/COMMUNICATION.md`
   - farm readiness still says any unresolved memory warning blocks;
   - execution workflow still omits pushed-state synchronization.

3. `docs/memory/templates/CODEX_CONTRACT.md`
   - future contracts still require `--fail-on-warning` as the ordinary
     task-final gate;
   - still mandates a final-pre-push reconciliation formulation;
   - workflow still omits pushed-state synchronization.

4. `docs/memory/CODEX.md`
   - the candidate prose distinguishes blocking and nonblocking warnings, but
     the unchanged contract-policy paragraph still mandates the strict
     `--fail-on-warning` task-final command.

5. `docs/memory/MAINTENANCE.md`
   - the candidate distinguishes blocking and nonblocking warnings, but the
     unchanged required-maintenance sequence still mandates
     `--fail-on-warning`;
   - it still says the non-strict health command cannot complete a substantial
     task.

These are material workflow contradictions because AGENTS is read first at
startup and the contract template is reused for future implementation tasks.

---

## Allowed modifications

Modify only:

```text
docs/memory/AGENTS.md
docs/memory/CODEX.md
docs/memory/COMMUNICATION.md
docs/memory/MAINTENANCE.md
docs/memory/templates/CODEX_CONTRACT.md
docs/memory/manifest.json
```

Create exactly:

```text
docs/memory/phases/memory-bracketed-scientific-throughput-policy-surface-fix-task-contract.md
```

Do not modify any other existing candidate file unless a deterministic checker
proves one of the six allowlisted files cannot be made coherent without it. If
that occurs, STOP and report the blocker rather than expanding scope.

---

## Frozen files

All scientific/runtime/testing/plotting/production files remain frozen,
including:

```text
src/
testing/
tools/
run_Prod_Analysis.sh
main.py
```

Also freeze the already reviewed candidate content of:

```text
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/USER.md
docs/memory/LEARNINGS.md
docs/memory/roadmap/STATUS.md
docs/memory/decisions/memory-bracketed-scientific-throughput.md
docs/memory/phases/memory-bracketed-scientific-throughput-task-contract.md
```

apart from manifest entries necessarily changing because this repair changes
versioned memory files.

Do not change scientific status, authority, blocker, or NEXT.

---

## Required behavior after repair

### A. AGENTS startup authority

Update `docs/memory/AGENTS.md` so its execution-authority section preserves:

```text
Codex local changes
-> ChatGPT actual-diff review
-> user commit/push
-> ChatGPT pushed-state synchronization review
-> user farm run when required
-> ChatGPT evidence review
```

State explicitly that commit/push and pushed-state review are required
synchronization stages, not independent scientific phases.

Replace the blanket warning blocker with:

- hard integrity failures block;
- warnings block only when they create material ambiguity about current source
  identity, accepted evidence/status, frozen scientific interfaces, active
  scientific ownership/blocker, or exact NEXT;
- nonblocking warnings are recorded and batched at the next checkpoint/milestone.

Preserve all scientific boundaries and startup ordering.

### B. COMMUNICATION farm policy

Update `docs/memory/COMMUNICATION.md` consistently:

- all executable farm steps still must be tracked, reviewed, pushed, and
  pushed-state reviewed;
- hard failures/material active-state or provenance ambiguity block farm
  readiness;
- a nonblocking memory warning alone does not block an otherwise ready farm
  gate;
- nonblocking warning cleanup is batched;
- execution workflow includes ChatGPT pushed-state synchronization.

Do not weaken the tracked-owner requirement or farm safety rules.

### C. CODEX ordinary health semantics

In `docs/memory/CODEX.md`, replace the ordinary task-final strict-warning rule.

The normal task-final memory-health command becomes:

```text
<PYTHON> -B tools/check_memory_health.py --root .
```

The command still returns nonzero for hard health failures.

`--fail-on-warning` remains available and should be required only when:

- an explicit milestone audit contract asks for zero warnings;
- a memory-hardening task specifically owns warning elimination;
- or the active warning has already been classified as materially blocking.

A normal source/science cycle may complete with a recorded **nonblocking**
warning if all hard checks pass and the warning is explicitly scheduled for the
next checkpoint/milestone.

Preserve manifest write/check, byte-count reporting, actual-diff review, user
commit/push, and pushed-state review.

### D. MAINTENANCE ordinary health semantics

Update `docs/memory/MAINTENANCE.md` to match CODEX:

- ordinary task-final health uses the non-strict command;
- `--fail-on-warning` is a milestone/explicit-zero-warning mode, not the
  universal completion gate;
- remove the statement that the non-strict command cannot complete a
  substantial task;
- require warnings to be classified as blocking or nonblocking in the memory
  health report;
- blocking warnings/hard failures must be repaired before the next gate;
- nonblocking warnings may be carried to the next checkpoint/milestone.

Do not change size thresholds or schema integrity requirements.

### E. Future contract template

Update `docs/memory/templates/CODEX_CONTRACT.md` so future tasks cannot recreate
the old recursive policy.

The template must:

- include pushed-state synchronization in the execution workflow;
- use the ordinary non-strict health command by default;
- reserve `--fail-on-warning` for explicit milestone/zero-warning tasks;
- require classification of warnings as blocking/nonblocking;
- require a push-stable substantive NEXT without mandating a separate
  final-pre-push reconciliation task;
- state that a matching pushed candidate with materially accurate CURRENT
  proceeds directly to the substantive gate;
- state that cosmetic/historical wording drift is batched unless it changes
  active meaning.

### F. Checker/tool boundary

Do **not** modify `tools/check_memory_health.py` or
`testing/test_memory_health.py`.

The tool already supports both modes:

```text
check_memory_health.py --root .
check_memory_health.py --root . --fail-on-warning
```

The repair changes **when each existing mode is policy-required**, not the
checker implementation.

---

## Positive checks

Prove all of the following:

1. `docs/memory/AGENTS.md` contains pushed-state synchronization.
2. `docs/memory/COMMUNICATION.md` contains pushed-state synchronization.
3. `docs/memory/templates/CODEX_CONTRACT.md` contains pushed-state
   synchronization.
4. The three surfaces no longer say every unresolved warning blocks the next
   gate.
5. CODEX and MAINTENANCE use non-strict health as the ordinary completion gate.
6. CODEX, MAINTENANCE, and the template retain `--fail-on-warning` only as an
   explicit stricter/milestone mode.
7. No file says a nonblocking warning alone requires a new reconciliation
   phase.
8. CURRENT remains byte-identical to the reviewed candidate and still has one
   NEXT.
9. The decision record remains byte-identical to the reviewed candidate.
10. Scientific/runtime/testing/plotting/production source remains unchanged.

---

## Local validation

Use the repository-authoritative Python from `TOOLS.md`.

Run:

```bash
python -B tools/update_memory_manifest.py --root . --write
python -B tools/update_memory_manifest.py --root . --check
python -B tools/memory_bootstrap.py --root .
python -m unittest testing.test_memory_health
python -B tools/check_memory_health.py --root .
python -B tools/check_memory_health.py --root . --fail-on-warning
git -c core.safecrlf=false diff --check
```

For this one repair checkpoint, both health modes must pass because the current
candidate has zero warnings. Running strict mode here verifies cleanliness; it
does **not** make strict mode the default policy for future ordinary cycles.

Run targeted text checks that fail if any of these stale policy forms remain in
the live authoritative surfaces:

```text
unresolved warnings block the next phase
unresolved memory-health warnings also block farm readiness
Do not hand off to another gate with unresolved warnings
Final-pre-push reconciliation must make
Workflow: Codex local changes -> ChatGPT audit -> user commit/push -> user farm run
```

Do not apply that stale-text scan to historical phase/task-contract records.

---

## Diff audit

The final cumulative worktree may contain only:

### First candidate

```text
docs/memory/CODEX.md
docs/memory/CURRENT.md
docs/memory/LEARNINGS.md
docs/memory/MAINTENANCE.md
docs/memory/MEMORY.md
docs/memory/USER.md
docs/memory/decisions/memory-bracketed-scientific-throughput.md
docs/memory/manifest.json
docs/memory/phases/memory-bracketed-scientific-throughput-task-contract.md
docs/memory/roadmap/STATUS.md
```

### This narrow repair

```text
docs/memory/AGENTS.md
docs/memory/COMMUNICATION.md
docs/memory/templates/CODEX_CONTRACT.md
docs/memory/phases/memory-bracketed-scientific-throughput-policy-surface-fix-task-contract.md
```

`docs/memory/CODEX.md`, `docs/memory/MAINTENANCE.md`, and
`docs/memory/manifest.json` are shared candidate+repair paths.

The only other allowed worktree file is the temporary repository-root review
bundle.

Create one fresh cumulative review bundle in the repository root that contains:

- full `git status --short --untracked-files=all`;
- committed branch/HEAD;
- complete tracked diff;
- complete `git diff --no-index /dev/null ...` for each new contract/decision
  file;
- exact final changed-path inventory;
- all required check outputs.

Use a new timestamped name such as:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

Do not stage merely for review.

---

## Acceptance criteria

PASS only if:

1. committed HEAD remains exactly
   `712ba32b062772d87fe44efa6865346e5c438827`;
2. the first candidate remains intact except the three shared policy/manifest
   files intentionally refined here;
3. AGENTS, CODEX, COMMUNICATION, MAINTENANCE, and the contract template express
   one coherent workflow;
4. pushed-state synchronization remains mandatory;
5. push synchronization does not become an independent scientific phase;
6. hard failures/material ambiguity block;
7. nonblocking warnings are batchable and do not automatically block;
8. non-strict health is the ordinary completion mode;
9. strict warning mode remains available for explicit milestone/zero-warning
   audits;
10. no checker/test/scientific/runtime source changes;
11. CURRENT scientific state/NEXT is unchanged from the first reviewed
    candidate;
12. manifest check passes;
13. bootstrap passes;
14. all 36 existing memory-health tests pass;
15. both ordinary and strict health commands pass on this zero-warning
    candidate;
16. diff check passes;
17. no commit, push, or farm run occurs.

---

## Hard stop

Stop after this narrow policy-surface repair, deterministic checks, manifest
regeneration, and creation of one fresh cumulative review bundle.

Do not:

- change science;
- change checker implementation;
- redesign memory architecture;
- add another policy layer;
- commit or push;
- run the farm;
- begin F.4.Refresh.2 materialization.

Return the new cumulative review bundle for one ChatGPT actual-diff review.
After PASS, the user performs one commit/push for the complete checkpoint,
ChatGPT performs one lightweight pushed-state synchronization review, and the
project proceeds directly to the F.4.Refresh.2 farm gate.
