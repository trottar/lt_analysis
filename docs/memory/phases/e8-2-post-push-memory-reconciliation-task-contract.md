# KaonLT E.8.2 Post-Push Memory Reconciliation — Task Contract

## 1. Objective

Perform one **memory-only post-push reconciliation** after the accepted E.8.2 source/test/memory set was pushed to `test`.

This task does **not** change analysis source, tests, physics, runtime behavior, profiles, collectors, wrappers, farm scripts, or validation artifacts.

The sole purpose is to make the durable repository memory agree with the already-pushed and independently reviewed E.8.2 state:

- E.8 remains `ACTIVE`.
- E.8.2 is `SOURCE REVIEWED`.
- E.8.3 is `NEXT`.
- F.6.3 is `BLOCKED` pending source-reviewed E.8.3.
- E.8.4 remains `BLOCKED` pending F.6.3.
- final E.8 remains `BLOCKED` pending E.8.4 and its later runtime/visual gate.
- F.6.4 remains `BLOCKED` pending completed F.6.3/E.8.4 production-impact evidence.
- the remaining E.8.1 canonical-five expansion remains `DEFERRED` by user decision.
- no E.8.2 farm/runtime acceptance is created by this documentation task.

After reconciliation, the exact next action must be the E.8.3 detached/non-production Method-A artifact/source-path audit and contract preparation. Do **not** implement E.8.3 in this task.

---

## 2. Exact Starting Identity

Required branch:

```text
test
```

Required starting HEAD:

```text
91bb7809d27d84d7709a6600dbb0dc9ab514a458
```

Commit subject:

```text
Add E8.2 baseline full-analysis stage audit
```

Parent:

```text
5b2ba7f92dc51b303257dbf209df687415364602
```

Before editing, run:

```bash
git status --short --branch
git branch --show-current
git rev-parse HEAD
git show --stat --oneline --decorate HEAD
```

Hard stop if:

- the branch is not `test`;
- HEAD is not exactly `91bb7809d27d84d7709a6600dbb0dc9ab514a458`;
- the worktree contains unrelated changes beyond this user-placed task contract.

Do not reset, stash, clean, discard, or rewrite unrelated user work.

---

## 3. Required Startup Reading

Read the normal startup sequence first, in this exact order:

1. root `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read only the task-relevant records:

- `docs/memory/CODEX.md`
- `docs/memory/MAINTENANCE.md`
- `docs/memory/TOOLS.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/phase-f6-method-a-production-promotion.md`
- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/phases/e8-2-baseline-stage-audit.md`
- this task contract

Do not broaden into unrelated phase history.

---

## 4. Established Evidence to Preserve

The pushed `test` source identity is already:

```text
91bb7809d27d84d7709a6600dbb0dc9ab514a458
```

The pushed-state review established that the E.8.2 implementation/source/test/memory set matches the independently reviewed candidate.

E.8.2 is therefore:

```text
SOURCE REVIEWED
```

only.

This does **not** establish:

- ROOT/PyROOT validation;
- full `main.py` integration;
- Jefferson Lab farm execution;
- procedure-PDF runtime rendering;
- E.8.2 runtime acceptance;
- Method-A production promotion;
- F.6.3 implementation;
- final E.8 closure.

The accepted scientific boundaries remain frozen:

- baseline production branch unchanged;
- active profile remains `no_empirical_residual`;
- legacy empirical residual Fit 1/Fit 2 remain dormant/excluded;
- Method A remains detached/non-production;
- Method B remains diagnostic/cross-check only and numerically excluded from Method-A application;
- frozen F.6.2 JSON and fingerprints remain unchanged;
- no canonical-child independent renormalization is permitted.

---

## 5. Problem to Repair

Several durable records still describe the pre-push/pre-source-review dependency state.

The stale forms include:

- E.8.2 labeled `NEXT`;
- E.8.3 labeled `BLOCKED` pending E.8.2 source review;
- F.6.3 labeled as blocked on both E.8.2 and E.8.3;
- E.8.2 phase/CURRENT next-action text still saying the user must commit/push the already-pushed E.8.2 set.

These statements are now stale because E.8.2 is already pushed and `SOURCE REVIEWED`.

This task must reconcile those records without changing any scientific decision.

---

## 6. Allowed Files

Only the following existing files may be edited:

```text
docs/memory/CURRENT.md
docs/memory/roadmap/STATUS.md
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/decisions/e8-full-analysis-procedure-roadmap.md
docs/memory/phases/e8-2-baseline-stage-audit.md
docs/memory/manifest.json
```

This new contract file is also allowed:

```text
docs/memory/phases/e8-2-post-push-memory-reconciliation-task-contract.md
```

No other file may change.

---

## 7. Frozen Files and Interfaces

Everything outside the allowlist is frozen.

In particular, do not edit:

```text
src/**
testing/**
config/**
scripts/**
*.sh
```

Do not edit:

- analysis code;
- plotting code;
- background profiles;
- proton subtraction;
- pion subtraction;
- HGCer Method A;
- HGCer Method B;
- SIMC;
- yield extraction;
- efficiencies;
- acceptance;
- cross sections;
- collectors;
- validation profiles;
- bundle wrappers;
- farm launchers;
- E.8.1/F.6.2 evidence records.

Do not modify any scientific artifact or fingerprint.

---

## 8. Required Memory Changes

### 8.1 `docs/memory/CURRENT.md`

Preserve E.8 as `ACTIVE` and E.8.2 as `SOURCE REVIEWED`.

Update the stale post-review action so that the exact next action is no longer:

```text
user-controlled commit/push ... -> pushed-state review -> E.8.3
```

because that push/review has already occurred.

Replace it with a concise exact NEXT equivalent to:

```text
NEXT — audit the accepted F.4/F.5/F.6.1/F.6.2 Method-A artifact/source path needed by E.8.3, then write the standalone E.8.3 implementation contract.
```

The audit is source/artifact-path analysis only. No E.8.3 implementation occurs in this task.

Record the pushed E.8.2 source identity where appropriate:

```text
91bb7809d27d84d7709a6600dbb0dc9ab514a458
```

Do not imply this commit has farm/runtime acceptance.

Keep:

- E.8.3 = `NEXT`;
- F.6.3 = `BLOCKED` pending source-reviewed E.8.3;
- E.8.4/final E.8/F.6.4 downstream blockers unchanged;
- E.8.1 canonical-five expansion = `DEFERRED`.

Do not expand CURRENT unnecessarily; it is the concise active-state authority.

### 8.2 `docs/memory/roadmap/STATUS.md`

Change only dependency/status wording required by the accepted pushed state:

```text
E.8.2 — SOURCE REVIEWED
E.8.3 — NEXT
F.6.3 — BLOCKED pending source-reviewed E.8.3
```

Preserve the approved sequence:

```text
E.8.2 -> E.8.3 -> F.6.3 -> E.8.4 -> final E.8 -> F.6.4
```

Do not change scientific ownership or farm-validation semantics.

### 8.3 `docs/memory/phases/phase-f6-method-a-production-promotion.md`

Reconcile the same status/dependency transition:

- E.8.2: `SOURCE REVIEWED`, not `NEXT`;
- E.8.3: `NEXT`, not blocked on E.8.2;
- F.6.3: remains `BLOCKED`, now only pending source-reviewed E.8.3.

Preserve all existing F.1–F.6.2 scientific/runtime closures and the F.6.3/F.6.4 ownership boundaries.

Where the record says the “active successor” is the full `E.8.2 -> E.8.3 ...` forward sequence, update only enough wording to make clear that E.8.2 is complete at source-review level and E.8.3 is the current successor.

Do not rewrite the roadmap architecture.

### 8.4 `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`

Reconcile status labels only:

- E.8.2 = `SOURCE REVIEWED`;
- E.8.3 = `NEXT`;
- F.6.3 = `BLOCKED` pending source-reviewed E.8.3.

Do not alter the governing E.8 scientific decision, page/content requirements, or milestone cadence.

### 8.5 `docs/memory/phases/e8-2-baseline-stage-audit.md`

Preserve its `SOURCE REVIEWED` status and implementation chronology.

Replace the stale `## Next` text that still requests the user-controlled E.8.2 commit/push.

The updated next section must state that:

- the reviewed E.8.2 set was pushed at `91bb7809d27d84d7709a6600dbb0dc9ab514a458`;
- pushed-state review passed;
- no E.8.2 farm/runtime validation is claimed;
- E.8.3 artifact/source-path audit is now next;
- F.6.3 remains blocked pending source-reviewed E.8.3.

Do not rewrite the implementation history.

### 8.6 `docs/memory/manifest.json`

Regenerate after all versionable memory changes, including this task contract.

Do not hand-edit semantic content into the manifest.

---

## 9. Scientific Ownership

This task has **zero scientific ownership**.

It must not:

- recompute or reinterpret Method A;
- change `C_j`;
- change `w0_j`;
- construct `B_pi^A`;
- alter pion templates;
- alter production histograms/yields;
- alter E.8.2 rendering/source;
- alter accepted F.6.2 evidence;
- alter Method-B treatment;
- decide production promotion.

It only reconciles durable status/dependency documentation after the accepted push.

---

## 10. Before / After Behavior

### Before

Repository memory contains a mixed state:

- CURRENT already recognizes E.8.2 as `SOURCE REVIEWED` and E.8.3 as `NEXT`;
- CURRENT and the E.8.2 phase record still contain a stale “commit/push E.8.2” next action;
- roadmap/Phase-F/decision records still label E.8.2 `NEXT` and E.8.3 blocked on it.

### After

All durable records agree on:

```text
E.8      — ACTIVE
E.8.2    — SOURCE REVIEWED
E.8.3    — NEXT
F.6.3    — BLOCKED pending source-reviewed E.8.3
E.8.4    — BLOCKED pending F.6.3
final E.8 — BLOCKED pending E.8.4 + later runtime/visual gate
F.6.4    — BLOCKED pending completed F.6.3/E.8.4 evidence
```

The exact next action in CURRENT is the E.8.3 accepted-artifact/source-path audit and later contract preparation.

No source/runtime behavior changes.

---

## 11. Positive Checks

Verify by direct inspection that:

1. all five status-bearing records agree on E.8.2/E.8.3/F.6.3;
2. no record still asks for the already-completed E.8.2 commit/push;
3. pushed source identity `91bb7809d27d84d7709a6600dbb0dc9ab514a458` is recorded only as source/push provenance, not runtime evidence;
4. E.8.3 remains detached/non-production;
5. F.6.3 remains the sole owner of constructing the parallel Method-A full-analysis branch;
6. Method B remains numerically excluded;
7. empirical residual Fit 1/Fit 2 remain dormant/excluded;
8. E.8.1 canonical-five expansion remains `DEFERRED`;
9. no F.1–F.6.2 runtime closure is downgraded.

---

## 12. Negative Checks

Search the edited memory records for stale or forbidden state.

At minimum inspect for occurrences equivalent to:

```text
E.8.2 — NEXT
E.8.3 — BLOCKED pending E.8.2
user-controlled commit/push of the reviewed E.8.2
F.6.3 — BLOCKED pending source-reviewed E.8.2 and E.8.3
```

Any remaining active-state occurrence must be justified as historical quotation/context rather than current status.

Do not “fix” historical task-contract text that intentionally records earlier state.

---

## 13. Local Validation

Discover a working interpreter as instructed by `docs/memory/TOOLS.md`; do not assume a workstation-specific Python executable.

Run the memory checks in the repository root.

Regenerate the manifest:

```bash
<PYTHON> -B tools/update_memory_manifest.py --root . --write
```

Then run:

```bash
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/check_memory_health.py --root .
<PYTHON> -B tools/memory_bootstrap.py --root . --json
git diff --check
```

Also run scoped diff inspection:

```bash
git status --short
git diff --name-only
git diff --stat
git diff -- \
  docs/memory/CURRENT.md \
  docs/memory/roadmap/STATUS.md \
  docs/memory/phases/phase-f6-method-a-production-promotion.md \
  docs/memory/decisions/e8-full-analysis-procedure-roadmap.md \
  docs/memory/phases/e8-2-baseline-stage-audit.md \
  docs/memory/manifest.json
```

There are no analysis/unit tests required because this task changes memory only.

Do not run ROOT, PyROOT, `main.py`, or farm analysis.

---

## 14. Farm Boundary

No Jefferson Lab farm execution is authorized or required.

This task must not claim:

```text
CLOSED / RUNTIME VALIDATED
```

for E.8.2 or E.8.3.

No procedure PDF or runtime artifact is produced or accepted here.

---

## 15. Diff Audit / Review Bundle

This is a small memory-only diff, but produce a complete review bundle in the repository root so ChatGPT can inspect the actual proposal.

Use:

```bash
git diff -- \
  docs/memory/CURRENT.md \
  docs/memory/roadmap/STATUS.md \
  docs/memory/phases/phase-f6-method-a-production-promotion.md \
  docs/memory/decisions/e8-full-analysis-procedure-roadmap.md \
  docs/memory/phases/e8-2-baseline-stage-audit.md \
  docs/memory/manifest.json \
  > kaonlt_review.diff

git diff --no-index -- /dev/null \
  docs/memory/phases/e8-2-post-push-memory-reconciliation-task-contract.md \
  >> kaonlt_review.diff || true
```

`|| true` is allowed only for `git diff --no-index` returning status 1 because it displayed a difference.

Do not use `|| true` to suppress any validation failure.

Before stopping, confirm:

```bash
git status --short
git diff --check
```

Do not commit or push.

---

## 16. Acceptance Criteria

The task passes local implementation review only if:

- branch and starting HEAD are exact;
- only allowlisted files changed;
- the task contract itself is present;
- E.8.2 is consistently `SOURCE REVIEWED`;
- E.8.3 is consistently `NEXT`;
- F.6.3 is consistently blocked only on source-reviewed E.8.3;
- CURRENT no longer requests the completed E.8.2 push/review step;
- E.8.2 phase chronology no longer requests the completed push;
- pushed identity is not mistaken for runtime validation;
- all scientific/runtime ownership remains unchanged;
- manifest is regenerated and passes its check;
- memory-health/bootstrap checks pass;
- `git diff --check` passes;
- the complete `kaonlt_review.diff` is produced for independent ChatGPT review;
- Codex does not commit, push, or run the farm.

After Codex stops, ChatGPT reviews the actual diff. Only after that review passes may the user commit/push this memory-only reconciliation.

---

## 17. Hard Stop

Stop and report without making broader changes if:

- starting branch/HEAD is wrong;
- unrelated dirty worktree changes exist;
- satisfying the reconciliation would require editing analysis source or tests;
- a memory-health failure indicates a broader migration/design problem;
- the status transition cannot be made without changing an accepted scientific dependency;
- any direct evidence contradicts E.8.2 = `SOURCE REVIEWED`.

Do not reinterpret the architecture, reopen closed phases, or begin E.8.3 implementation inside this task.
