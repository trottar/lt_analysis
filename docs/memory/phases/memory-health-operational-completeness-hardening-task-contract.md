# KaonLT Memory Health and Farm-Gate Operational-Completeness Hardening

## 1. Purpose

Stop all further scientific/farm progression and repair two workflow failures
before any F.4.Refresh.2 farm validation continues:

1. **Operational-completeness failure:** F.4.Refresh.2.Validation.1 was allowed
   to reach pushed/source-reviewed profile status and ChatGPT then began
   composing an interactive farm sequence even though the repository had no
   tracked, reviewed owner for the complete multi-step
   `materialize -> verify -> package` operation.
2. **Memory-health failure:** repeated `CURRENT.md` soft-size warnings were
   reported but treated as acceptable instead of triggering the focused
   maintenance required by `docs/memory/MAINTENANCE.md`.

This is a repository-memory/workflow/tooling hardening task only.

It must **not** implement the missing F.4.Refresh.2 execution driver yet.
It must **not** run or prepare a farm command.
It must **not** modify scientific/runtime/production source.

The purpose of this task is to make the failure durable, make health reporting
mandatory and machine-enforceable, compact CURRENT back to a healthy size, and
make future farm-readiness handoff explicitly dependent on a complete reviewed
execution path rather than merely a reviewed profile.

---

## 2. Exact starting state

Required branch:

```text
test
```

Required committed starting HEAD:

```text
86590fa655512926f2e4d0c50bf12b57d5198da5
```

Required subject:

```text
F4 Refresh.2 Validation.1: add farm materialization bundle profile
```

Before editing, run:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
git log -1 --oneline
```

Hard requirements:

- branch must be `test`;
- committed HEAD must be exactly
  `86590fa655512926f2e4d0c50bf12b57d5198da5`;
- inspect the actual worktree first;
- root `AGENTS.md` and `.codex/` remain local-only/untracked;
- pre-existing temporary `kaonlt_review*.diff` files may remain untracked and
  must not be deleted, rewritten, staged, or treated as candidate files;
- if unrelated tracked changes exist, STOP and report them;
- do not reset, stash, clean, discard, commit, push, or run the farm.

---

## 3. Mandatory startup

Read in this exact order:

1. root `AGENTS.md`;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read only the task-relevant records/source:

- `docs/memory/AGENTS.md`
- `docs/memory/CODEX.md`
- `docs/memory/COMMUNICATION.md`
- `docs/memory/MAINTENANCE.md`
- `docs/memory/LEARNINGS.md`
- `docs/memory/TOOLS.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/templates/CODEX_CONTRACT.md`
- `docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization.md`
- `docs/memory/phases/f4-refresh2-validation1-farm-materialization-bundle-profile.md`
- `docs/memory/decisions/farm-validation-bundle-procedure.md`
- `tools/check_memory_health.py`
- `testing/test_memory_health.py`

Inspect the current pushed source for:

- `testing/materialize_method_a_current_baseline_authority.py`
- `testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json`
- `testing/package_pion_hgcer_validation_bundle.tcsh`

only to establish the operational gap. Do not edit them.

---

## 4. Established facts that must be preserved

### 4.1 Pushed source identity

Current remote/pushed `test` source is:

```text
86590fa655512926f2e4d0c50bf12b57d5198da5
```

with parent:

```text
141a3d04f9e5d07be21dba14e0e63212c3990bf1
```

### 4.2 F.4.Refresh.2 source status

Preserve:

```text
F.4.Refresh.1              — CLOSED / RUNTIME VALIDATED
F.4.Refresh.1.Fix.1        — CLOSED / RUNTIME VALIDATED
                              detached comparison gate only

F.4.Refresh.2              — SOURCE REVIEWED
F.4.Refresh.2.Fix.1        — SOURCE REVIEWED
F.4.Refresh.2.Validation.1 — SOURCE REVIEWED

F.6.3                      — SOURCE REVIEWED
E.8.4                      — SOURCE REVIEWED
E.8                        — ACTIVE
final E.8                  — BLOCKED
F.6.4                      — BLOCKED
lifecycle hook             — BLOCKED / DEFERRED
```

Do not downgrade the source-reviewed materializer or profile.

### 4.3 No F.4.Refresh.2 materialization occurred

The attempted farm invocation failed before Python could open:

```text
testing/materialize_method_a_current_baseline_authority.py
```

in the stale ordinary farm checkout.

Therefore no F.4.Refresh.2 materializer execution, candidate output, bundle ZIP,
accepted-authority update, F.6.3/E.8.4 runtime acceptance, or Method-A promotion
resulted from that attempt.

### 4.4 The materializer itself is not missing from pushed source

The pushed repository contains:

```text
testing/materialize_method_a_current_baseline_authority.py
```

with the reviewed blob inherited from the F.4.Refresh.2 materializer commit.

The actual repository-level gap is that no tracked, reviewed F.4.Refresh.2
operation owns the complete later farm sequence:

```text
validated inputs
  -> candidate materialization
  -> output verification
  -> validation-bundle collection
  -> returned review ZIP
```

The existing generic package wrapper is intentionally bundle-only and cannot
run the materializer. This task must record that distinction exactly.

### 4.5 Memory-health state

At pushed HEAD `86590fa...`, `docs/memory/CURRENT.md` is 16,338 bytes.

Existing memory policy defines:

```text
CURRENT soft warning = 8 KiB
CURRENT hard failure = 16 KiB
```

so CURRENT is only 46 bytes below the 16-KiB hard threshold.

`docs/memory/MAINTENANCE.md` already says a soft-limit warning calls for focused
consolidation. Repeated local checks reported the warning, but the workflow
continued. This task must treat that as a process failure, not as an acceptable
long-term warning.

---

## 5. Scientific/runtime freeze

Do not modify:

- any `src/` file;
- `src/main.py`;
- `run_Prod_Analysis.sh`;
- any F.1/F.2/F.3/F.4 scientific builder;
- `testing/materialize_method_a_current_baseline_authority.py`;
- `testing/test_materialize_method_a_current_baseline_authority.py`;
- `testing/compare_method_a_current_baseline_authority.py`;
- `testing/collect_pion_hgcer_validation_bundle.py`;
- `testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json`;
- `testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py`;
- `testing/package_pion_hgcer_validation_bundle.tcsh`;
- accepted F.2/F.3/F.4 authority constants/artifacts;
- F.5/F.6.3/E.8.4 scientific or runtime source;
- Method B;
- production physics.

Do not implement the missing execution driver in this task.

---

## 6. Allowed files

Existing files that may be edited:

```text
docs/memory/AGENTS.md
docs/memory/CODEX.md
docs/memory/COMMUNICATION.md
docs/memory/CURRENT.md
docs/memory/LEARNINGS.md
docs/memory/MAINTENANCE.md
docs/memory/MEMORY.md
docs/memory/TOOLS.md
docs/memory/USER.md
docs/memory/manifest.json
docs/memory/roadmap/STATUS.md
docs/memory/templates/CODEX_CONTRACT.md
docs/memory/phases/f4-refresh2-validation1-farm-materialization-bundle-profile.md
tools/check_memory_health.py
testing/test_memory_health.py
```

Create exactly:

```text
docs/memory/investigations/f4-refresh2-validation1-operational-readiness-failure.md
docs/memory/phases/memory-health-operational-completeness-hardening-task-contract.md
docs/memory/phases/memory-health-operational-completeness-hardening.md
```

No other path may change.

---

## 7. Required durable failure record

Create:

```text
docs/memory/investigations/f4-refresh2-validation1-operational-readiness-failure.md
```

It must record, concisely and factually:

1. Validation.1 correctly source-reviewed a bundle profile, not a complete
   execution driver.
2. ChatGPT incorrectly treated profile/pushed-state readiness as sufficient to
   begin interactive farm orchestration.
3. The pushed materializer existed; the stale ordinary farm checkout was an
   observed condition but **not the repository-design root cause**.
4. The root workflow defect was failure to audit the full farm execution path
   before authorizing a farm step.
5. The generic package wrapper is bundle-only and cannot own materialization.
6. No candidate materialization or accepted-authority mutation occurred.
7. Repeated CURRENT soft-size warnings were allowed to persist even though
   MAINTENANCE requires focused consolidation.
8. Future prevention is owned by the operational-completeness gate and strict
   health-report gate introduced by this task.

Do not write emotional or blame-oriented prose. Preserve the exact technical
failure and corrective rule.

---

## 8. Operational-completeness gate

Reinforce tracked memory so that **before any farm command is provided**,
ChatGPT/Codex must establish a complete farm-operation ownership chain.

The durable rule must require an explicit audit of:

```text
input authority
  -> producer/materializer/analyzer/renderer, if any
  -> verification/checker, if any
  -> collector/packager
  -> invocation owner
  -> expected returned artifact
```

### 8.1 Required rule

A farm gate is not ready merely because a profile or collector exists.

Every executable step required by the requested operation must already be:

- represented by tracked source;
- locally deterministic where possible;
- independently source reviewed;
- pushed;
- pushed-state reviewed.

If the requested farm operation requires multi-step orchestration and no
tracked/reviewed driver owns that orchestration, the farm gate is:

```text
BLOCKED
```

and work returns to the source-changing workflow.

Do **not** assemble the missing orchestration interactively in the farm shell.

### 8.2 Direct-command exception

Preserve simple direct farm operations when one reviewed repository CLI already
owns the complete requested action and no ad-hoc staging/orchestration is
required.

Do not force wrapper proliferation for genuinely single-command operations.

### 8.3 Required pre-farm communication

Before ChatGPT provides a farm command, the response must state:

```text
Farm readiness: PASS
```

and name the tracked source path(s) that own the complete requested operation.

If those source owners do not exist or are not pushed/source-reviewed, state:

```text
Farm readiness: BLOCKED
```

and provide no farm command.

Add this ownership rule to the appropriate specialized records without
duplicating full procedures:

- `CODEX.md` — source/task completeness and contract requirement;
- `COMMUNICATION.md` — pre-farm handoff requirement;
- `USER.md` — stable collaboration boundary;
- `MEMORY.md` — concise durable cross-phase rule;
- `LEARNINGS.md` — generalized lesson;
- tracked `AGENTS.md` — concise governing pointer/rule;
- `templates/CODEX_CONTRACT.md` — mandatory contract section.

Update the F.4.Refresh.2.Validation.1 phase and roadmap only enough to state that
the profile remains SOURCE REVIEWED but the farm execution gate is BLOCKED
until the missing tracked execution owner is implemented and reviewed.

---

## 9. Mandatory memory-health update gate

The absence of explicit health reporting must also be repaired.

### 9.1 Health update required after every substantial gate

After every substantial implementation, source review reconciliation, closure,
or pushed-state handoff, Codex and ChatGPT must surface a concise health update.

The update must contain:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
manifest check: PASS | FAIL
```

Codex reports the deterministic local result.
ChatGPT independently verifies the available evidence and must not silently omit
the health state from a substantial-gate conclusion.

### 9.2 Warnings are workflow-stopping

A memory-health warning is not a scientific/runtime failure, but it is a
workflow stop for beginning a new phase/fix/farm gate.

No new gate may proceed until the warning is resolved or an explicit,
user-approved maintenance exception is recorded.

The current repeated practice of saying the CURRENT soft warning is merely
"nonfatal" and then proceeding is forbidden.

### 9.3 CURRENT compaction now

Compact `docs/memory/CURRENT.md` substantially.

Requirements:

- preserve schema 3;
- preserve the exact required section order;
- preserve exactly one ordinary NEXT;
- preserve the active E.8 objective;
- preserve only concise active/verified state needed to resume;
- replace historical prose with links to the existing phase/evidence records;
- preserve the current F.4.Refresh.1 closure and F.4.Refresh.2 source-reviewed
  statuses;
- record the new operational-readiness blocker and memory-hardening task;
- do not lose frozen F.6.2 identity/boundary facts needed for safe resume;
- do not duplicate the roadmap or phase history.

Target:

```text
CURRENT.md <= 7 KiB
```

This leaves real margin below the existing 8-KiB soft warning.

---

## 10. Machine-enforced strict health mode

Modify:

```text
tools/check_memory_health.py
testing/test_memory_health.py
```

narrowly.

### 10.1 CLI

Add an explicit CLI option:

```text
--fail-on-warning
```

Behavior:

- existing structural/semantic errors remain failures;
- without the flag, preserve current warning-report behavior for compatibility;
- with `--fail-on-warning`, any health warning causes nonzero CLI exit;
- warnings must still be printed/reported distinctly from structural errors.

Do not change the existing 8/16-KiB thresholds in this task.

### 10.2 Tests

Add deterministic coverage that proves:

1. a healthy fixture passes with `--fail-on-warning`;
2. a CURRENT fixture above the soft limit but below hard limit still reports a
   warning;
3. the same fixture returns nonzero under `--fail-on-warning`;
4. the default compatibility mode retains its existing warning-only behavior;
5. a hard-limit failure remains a failure with or without strict mode.

Use the existing test structure; no external dependencies.

### 10.3 Canonical task-final health command

Update `TOOLS.md`, `MAINTENANCE.md`, `CODEX.md`, and the contract template so
the task-final gate uses:

```text
<PYTHON> -B tools/check_memory_health.py --root . --fail-on-warning
```

The non-strict form may remain documented for diagnostic inspection, but it is
not sufficient for completing a substantial task.

---

## 11. Contract-template reinforcement

Update:

```text
docs/memory/templates/CODEX_CONTRACT.md
```

so every future source-changing contract must explicitly contain:

### Operational readiness

- whether the task leads to a farm gate;
- the complete intended execution chain;
- the tracked source owner of every executable step;
- whether a driver/wrapper is required;
- the condition that missing orchestration blocks farm handoff.

### Memory health

- strict final health command;
- byte-size report for CURRENT/MEMORY/CURRENT_HANDOFF;
- manifest regeneration/check when memory changes;
- prohibition on handing off with unresolved health warnings.

Update `CODEX.md` consistently.

---

## 12. Current status and exact NEXT during this task

During implementation record:

```text
Memory health / operational-completeness hardening — ACTIVE
```

Preserve F.4.Refresh.2.Validation.1 as `SOURCE REVIEWED`, but record:

```text
F.4.Refresh.2 farm execution gate — BLOCKED
```

because the complete tracked execution owner is missing.

No farm command is authorized.

Set exactly one CURRENT NEXT:

```text
NEXT — independent ChatGPT actual-diff review of the memory-health and farm-gate operational-completeness hardening candidate.
```

Do not advance the hardening task to SOURCE REVIEWED; ChatGPT owns that review.

---

## 13. Required local validation

Discover one working `<PYTHON>` and use it consistently.

Run:

```bash
<PYTHON> -B -m py_compile \
  tools/check_memory_health.py \
  testing/test_memory_health.py

<PYTHON> -B -m unittest testing.test_memory_health -v

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

- strict health exits 0;
- no health warnings remain;
- CURRENT is <= 7 KiB;
- manifest check passes;
- memory bootstrap exits 0;
- memory tests pass with exact count/skips reported;
- diff check passes.

If strict health reports any warning, STOP. Do not call it nonfatal.

---

## 14. Diff audit

Before stopping:

```bash
git status --short
git -c core.safecrlf=false diff --stat
git -c core.safecrlf=false diff --no-ext-diff
```

Audit every intended new file with:

```bash
git diff --no-index /dev/null <path> || true
```

`|| true` is permitted only for `git diff --no-index` status 1 after displaying
the diff.

No source file outside the allowlist may appear.

---

## 15. Required review bundle

Create one fresh repository-root bundle:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

It must contain byte-faithfully:

1. branch and committed HEAD;
2. `git status --short`;
3. complete cumulative tracked diff/stat;
4. complete no-index diffs for all new files;
5. exact local validation commands/results;
6. exact health update:
   - strict/non-strict health result;
   - CURRENT/MEMORY/HANDOFF byte counts;
   - warning list;
   - manifest result;
7. final changed-path inventory.

The review bundle remains temporary and untracked.

---

## 16. Acceptance criteria

This hardening candidate is acceptable for ChatGPT review only if:

1. committed base remains exactly `86590fa655512926f2e4d0c50bf12b57d5198da5`;
2. no scientific/runtime/production source changes;
3. the operational-readiness failure is durably recorded;
4. future farm handoff requires a complete tracked/pushed/reviewed execution
   chain, not merely a profile;
5. multi-step unowned orchestration is explicitly farm-BLOCKING;
6. ChatGPT/Codex health updates are mandatory after substantial gates;
7. unresolved health warnings block progression;
8. strict health mode is implemented and tested;
9. canonical task-final health command uses `--fail-on-warning`;
10. CURRENT is compacted to <= 7 KiB with one NEXT and valid schema 3;
11. strict health produces zero warnings and exits 0;
12. manifest/bootstrap/memory tests/diff check pass;
13. F.4.Refresh.2.Validation.1 remains SOURCE REVIEWED;
14. F.4.Refresh.2 farm execution remains BLOCKED;
15. no missing execution driver is implemented in this task;
16. no farm command/run occurs;
17. one fresh cumulative review bundle is produced.

---

## 17. Hard stop

After memory/tooling hardening, deterministic checks, manifest regeneration,
health report, diff audit, and creation of the fresh review bundle:

**STOP.**

Do not commit or push.

Do not implement the F.4.Refresh.2 execution driver.

Do not provide or run a farm command.

Do not run the materializer.

Do not run the collector.

Do not create a farm ZIP.

Do not modify accepted authority or scientific/runtime source.

Return the fresh review bundle for independent ChatGPT actual-diff review.
