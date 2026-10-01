# KaonLT Post-Push CURRENT Continuity and Push-Stable NEXT Hardening

## 1. Purpose

Repair one residual repository-memory continuity defect found during pushed-state
review of the already accepted memory-health / farm-readiness hardening commit.

Pushed commit:

```text
3fd4edfcd05fd4cc25fb4e2119fb5f212ee916a1
Harden memory health and farm readiness gates
```

Its parent is:

```text
86590fa655512926f2e4d0c50bf12b57d5198da5
```

The pushed commit itself matches the independently reviewed 19-path candidate.
No rollback or scientific repair is required.

However, pushed `docs/memory/CURRENT.md` still contains pre-push transitional
language:

```text
The candidate is local to committed test base 86590fa...
```

and:

```text
NEXT — user-controlled commit/push of the independently reviewed
memory-health and farm-gate operational-completeness hardening candidate.
```

Those statements became stale as soon as commit `3fd4edf...` was pushed.
Because CURRENT is the sole ordinary active-state authority, this is a live
continuity defect even though structural/strict memory-health checks and the
manifest pass.

This task is memory/workflow hardening only. It must:

1. reconcile CURRENT to the actual pushed state;
2. make the NEXT durable across future user-controlled commit/push transitions;
3. prevent final-pre-push reconciliation from writing a NEXT whose sole action
   is the very commit/push that will immediately consume and stale it.

Do not begin the missing F.4.Refresh.2 execution-driver task in this change.

## 2. Exact starting state

Required branch:

```text
test
```

Required committed HEAD:

```text
3fd4edfcd05fd4cc25fb4e2119fb5f212ee916a1
```

Required subject:

```text
Harden memory health and farm readiness gates
```

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
git log -1 --oneline
```

Hard requirements:

- branch must be `test`;
- committed HEAD must be exactly `3fd4edfcd05fd4cc25fb4e2119fb5f212ee916a1`;
- inspect the actual worktree;
- pre-existing `kaonlt_review*.diff` files may remain untracked and must not be
  deleted, staged, rewritten, or treated as candidate files;
- root `AGENTS.md` and `.codex/` remain local-only/untracked;
- if unrelated tracked changes exist, STOP;
- do not reset, stash, clean, discard, commit, push, or run the farm.

## 3. Mandatory startup

Read in exact order:

1. root `AGENTS.md`;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read only:

- `docs/memory/CODEX.md`
- `docs/memory/MAINTENANCE.md`
- `docs/memory/COMMUNICATION.md`
- `docs/memory/templates/CODEX_CONTRACT.md`
- `docs/memory/phases/memory-health-operational-completeness-hardening.md`
- `docs/memory/investigations/f4-refresh2-validation1-operational-readiness-failure.md`
- `docs/memory/phases/f4-refresh2-validation1-farm-materialization-bundle-profile.md`
- `docs/memory/roadmap/STATUS.md`

Do not reopen scientific architecture.

## 4. Established pushed-state evidence

Preserve these facts:

```text
memory-health / operational-completeness hardening — SOURCE REVIEWED
F.4.Refresh.2.Validation.1                           — SOURCE REVIEWED
F.4.Refresh.2 farm execution gate                   — BLOCKED
```

The pushed hardening commit contains exactly the reviewed 19 paths and no
temporary `kaonlt_review*.diff` files.

The reviewed/pushed hardening established:

- strict `--fail-on-warning` memory-health mode;
- unresolved warnings block new phase/fix/farm gates unless the user approves a
  recorded maintenance exception;
- Codex and ChatGPT both surface health after substantial gates;
- farm commands require `Farm readiness: PASS` and complete tracked ownership;
- missing multi-step orchestration yields `Farm readiness: BLOCKED`;
- no interactive farm-shell reconstruction of missing orchestration.

Do not weaken any of those rules.

## 5. Scientific/runtime freeze

Do not modify:

- any `src/` file;
- `src/main.py`;
- `run_Prod_Analysis.sh`;
- any scientific builder;
- F.4.Refresh.2 materializer/comparator;
- validation collector/profile/wrapper;
- any accepted authority;
- any F.6.3/E.8.4 scientific/runtime source;
- Method B;
- production physics.

Do not implement the missing F.4.Refresh.2 execution owner.

## 6. Allowed files

Edit only:

```text
docs/memory/CURRENT.md
docs/memory/CODEX.md
docs/memory/MAINTENANCE.md
docs/memory/templates/CODEX_CONTRACT.md
docs/memory/manifest.json
```

Create exactly:

```text
docs/memory/phases/post-push-current-continuity-task-contract.md
docs/memory/phases/post-push-current-continuity.md
```

No other path may change.

## 7. Push-stable NEXT rule

Add one concise durable rule to `CODEX.md`, `MAINTENANCE.md`, and the contract
template:

A final-pre-push reconciliation must not make the sole ordinary NEXT be the
commit/push action itself.

Instead, CURRENT's NEXT must be **push-stable**: it must remain accurate both
before and immediately after the user-controlled commit/push.

Preferred form:

```text
NEXT — after user-controlled commit/push and pushed-state review of this
candidate, <the actual next repository/scientific gate>.
```

After push, the condition is satisfied and the substantive next action remains
correct; CURRENT does not become stale merely because the push occurred.

A pushed-state review must explicitly check CURRENT for consumed pre-push
language. If CURRENT says work is still local/unpushed or tells the user to
perform a push that is already complete, memory continuity is `BLOCKED` until
reconciled. Do not begin a new scientific/farm gate from stale CURRENT.

Do not turn CURRENT into a commit ledger. This rule concerns only live-state
truth and NEXT durability.

## 8. CURRENT reconciliation

Update `docs/memory/CURRENT.md` minimally.

### 8.1 Current work item

Replace the stale statement that the hardening candidate is local to base
`86590fa...`.

Record instead that:

```text
memory-health / operational-completeness hardening — SOURCE REVIEWED
pushed source — 3fd4edfcd05fd4cc25fb4e2119fb5f212ee916a1
pushed-state source/provenance review — PASS
```

Do not claim runtime validation.

### 8.2 Preserve blocker

Keep:

```text
F.4.Refresh.2 farm execution gate — BLOCKED
```

until the missing complete tracked execution owner is separately implemented,
locally checked, independently source reviewed, pushed, and pushed-state
reviewed.

### 8.3 Exact push-stable NEXT

Set exactly one ordinary NEXT:

```text
NEXT — after user-controlled commit/push and pushed-state review of this continuity reconciliation, audit and contract the missing tracked F.4.Refresh.2 materialize -> verify -> package execution owner.
```

This task does not perform that later audit/implementation.

Keep CURRENT <= 7168 bytes.

## 9. Phase record

Create:

```text
docs/memory/phases/post-push-current-continuity.md
```

During implementation status is:

```text
ACTIVE
```

Record:

- pushed hardening commit/source identity;
- source push itself passed;
- strict checker/manifest evidence passed but did not detect the external
  push-consumed NEXT;
- the live-state mismatch is therefore a continuity blocker rather than a
  scientific/runtime failure;
- push-stable NEXT policy prevents recurrence.

Do not mark SOURCE REVIEWED; independent ChatGPT review owns that transition.

## 10. Required deterministic checks

Use the established `<PYTHON>`.

Run:

```bash
<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/update_memory_manifest.py --root . --check

<PYTHON> -B tools/check_memory_health.py --root .
<PYTHON> -B tools/check_memory_health.py --root . --fail-on-warning

<PYTHON> -B tools/memory_bootstrap.py --root . --json

wc -c \
  docs/memory/CURRENT.md \
  docs/memory/MEMORY.md \
  docs/memory/handoffs/CURRENT_HANDOFF.md

git -c core.safecrlf=false diff --check
```

Acceptance requires:

- manifest PASS;
- non-strict health PASS;
- strict health PASS;
- zero warnings;
- CURRENT <= 7168 bytes;
- bootstrap PASS;
- diff check PASS.

No unit suite is required because no checker implementation changes. If Codex
runs tests, report them accurately.

## 11. Health report

Return the exact task-final block:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
manifest check: PASS | FAIL
live-state continuity: PASS | BLOCKED
```

Acceptance requires both memory health and live-state continuity PASS.

## 12. Diff audit and review bundle

Before stopping:

```bash
git status --short
git -c core.safecrlf=false diff --stat
git -c core.safecrlf=false diff --no-ext-diff
```

For new files:

```bash
git diff --no-index /dev/null <path> || true
```

Create one fresh:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

containing:

1. branch and committed HEAD;
2. status;
3. complete cumulative diff/stat;
4. no-index diffs for both new files;
5. exact validation outputs;
6. exact health/continuity report;
7. final candidate path inventory.

Keep the review bundle untracked.

Expected candidate paths:

```text
docs/memory/CODEX.md
docs/memory/CURRENT.md
docs/memory/MAINTENANCE.md
docs/memory/manifest.json
docs/memory/templates/CODEX_CONTRACT.md
docs/memory/phases/post-push-current-continuity-task-contract.md
docs/memory/phases/post-push-current-continuity.md
```

## 13. Acceptance criteria

The candidate passes only if:

1. committed base is exactly `3fd4edfcd05fd4cc25fb4e2119fb5f212ee916a1`;
2. no scientific/runtime source changes;
3. pushed hardening remains SOURCE REVIEWED;
4. F.4.Refresh.2 farm execution remains BLOCKED;
5. CURRENT no longer says the hardening work is local/unpushed;
6. CURRENT no longer asks for the already-completed hardening push;
7. CURRENT has exactly one push-stable NEXT;
8. push-stable NEXT policy is documented in CODEX/MAINTENANCE/template;
9. strict health passes with zero warnings;
10. live-state continuity reports PASS;
11. CURRENT <= 7168 bytes;
12. manifest/bootstrap/diff checks pass;
13. missing execution owner is not implemented;
14. no farm command/run occurs;
15. one fresh review bundle is produced.

## 14. Hard stop

After memory reconciliation, checks, health/continuity report, diff audit, and
fresh review bundle:

**STOP.**

Do not commit or push.

Do not start the F.4.Refresh.2 execution-owner implementation.

Do not provide or run a farm command.

Return the fresh review bundle for independent ChatGPT actual-diff review.
