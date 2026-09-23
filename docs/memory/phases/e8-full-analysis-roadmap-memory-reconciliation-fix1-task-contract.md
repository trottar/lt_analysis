# KaonLT E.8 Roadmap Memory Reconciliation — Narrow Repair Task Contract

## Objective

Repair the current uncommitted E.8 roadmap memory reconciliation without changing its approved new roadmap.

The first reconciliation diff is largely correct, but it over-compressed
`docs/memory/phases/phase-f6-method-a-production-promotion.md` and deleted durable
closed-phase scientific/validation constraints that the original contract explicitly
required to retain. It also rewrote the lifecycle-hook status in `CURRENT.md` as prose
instead of the exact project status label.

This repair is memory/documentation-only.

Do **not** revert the new E.8 roadmap, the new Fix.6 runtime evidence, or the new
dependency order.

---

## Exact starting source/worktree state

Committed branch/HEAD must remain:

- branch: `test`
- HEAD: `0ec29d4e1bb345eb37e8cca35b8b7e5cbe1b4d5b`

The worktree is expected to contain the already-reviewed first-pass reconciliation
changes from:

`docs/memory/phases/e8-full-analysis-roadmap-memory-reconciliation-task-contract.md`

plus this newly placed narrow repair contract.

Do not pull, reset, stash, clean, switch branches, commit, push, or run the farm.

Before editing, show:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
```

Hard stop if the branch or HEAD differs, or if unrelated changes exist outside the
current E.8 memory-reconciliation set.

---

## Required startup/read order

Read:

1. `docs/memory/AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read:

- `docs/memory/phases/e8-full-analysis-roadmap-memory-reconciliation-task-contract.md`
- `docs/memory/phases/phase-f6-method-a-production-promotion.md`
- the committed pre-edit version of
  `docs/memory/phases/phase-f6-method-a-production-promotion.md`
  at `0ec29d4e1bb345eb37e8cca35b8b7e5cbe1b4d5b`
- `docs/memory/CURRENT.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/manifest.json`

The committed `0ec29d4...` Phase-F record is the source for the durable details that
must be restored. Do not improvise new scientific rules.

---

## Problem found by independent diff review

The first-pass reconciliation correctly records the new E.8 program, but its rewrite
of `phase-f6-method-a-production-promotion.md` removed durable accepted information
from the detailed Phase-F record, including closed F.6.1/F.6.2 scientific and
validation boundaries.

The original reconciliation contract required:

> retain every frozen scientific boundary

The detailed Phase-F record must therefore preserve the accepted closed-phase content
while updating only the forward dependency/roadmap portions.

Examples of durable information that must remain represented in the detailed Phase-F
record include:

### F.6.1

- accepted runtime/bundle and reviewed-science identities;
- five settings / 15 parents / 135 cells closure scope;
- low-response population is a reference shape, not absolute leakage probability;
- explicit shape-normalization semantics;
- accepted Method-A model variables versus independent validation variables;
- validation variables do not become Method-A inputs;
- detached/shadow-only ownership;
- exact application identity matching and failure on missing/extra/duplicate records;
- frozen producer-radians to application-degrees phi contract;
- aggregate reweighting reproduces accepted F.5 closure;
- factors remain transient; no new authoritative event-factor table is persisted;
- F.6.1.Validation.1 bundle/profile provenance and fail-closed source gate.

### F.6.2

- baseline remains the guide for broad physical missing-mass structure;
- Method A is an acceptance-dependent refinement, not a replacement;
- low response remains a reference shape;
- pion leakage in kaon-bearing regions is permitted;
- accepted diagnostics remain detached from production;
- no composite score;
- no frozen numerical case thresholds;
- no Method-B numerical dependency;
- no claim of final-yield uncertainty reduction;
- F.6.2.Fix.5 presentation provenance remains distinct from scientific provenance.

### F.6.4 / uncertainty boundary

Retain the durable requirement that any later promotion contract must define, at
minimum:

- supported kinematics;
- enable/disable behavior;
- fail-closed authority;
- setting-wide atomicity;
- provenance;
- uncertainty treatment;
- downstream validation;
- regression against the unchanged baseline.

Also retain:

- this roadmap defines no Method-A systematic uncertainty;
- `|Method-A result - baseline result|` is a correction effect, not automatically an
  uncertainty;
- no automatic production promotion follows from detached validation or presentation.

These are examples of the deleted durable content; restore the complete warranted
closed-phase detail from the committed `0ec29d4...` version rather than only these
sentences.

---

## Required repair

### 1. `docs/memory/phases/phase-f6-method-a-production-promotion.md`

Reconcile additively.

Preserve/restore the detailed accepted F.6.1, F.6.1.Validation.1, F.6.2,
F.6.2.Fix.5, F.6.4, uncertainty, and forbidden-shortcut content from the committed
`0ec29d4...` version.

At the same time, retain the already approved new forward roadmap from the current
worktree:

- accepted narrow E.8.1 Fix.5/Fix.6 runtime evidence;
- E.8.2 baseline full-analysis audit;
- E.8.3 detached Method-A reweighting audit;
- F.6.3 parallel full procedure plus Method A;
- E.8.4 baseline-versus-Method-A production-impact audit;
- final E.8 closure;
- F.6.4 explicit production-promotion decision;
- milestone farm cadence.

The final detailed Phase-F file should read as a cumulative durable record, not a
compressed replacement.

Do not duplicate contradictory old dependency text such as “F.6.3 BLOCKED pending
E.8”. Replace that obsolete forward dependency with the new sequence while preserving
closed-phase technical detail.

### 2. `docs/memory/CURRENT.md`

Normalize the lifecycle status to the exact allowed wording:

`Lifecycle-hook dispatch remains BLOCKED / DEFERRED.`

Do not otherwise redesign CURRENT.

### 3. `docs/memory/manifest.json`

Regenerate after the repair and after adding this contract.

### 4. Track this repair contract

Add:

`docs/memory/phases/e8-full-analysis-roadmap-memory-reconciliation-fix1-task-contract.md`

to the memory manifest.

---

## Frozen files

Do not modify any file outside `docs/memory/`.

Within `docs/memory/`, do not change the already-correct first-pass files unless
manifest regeneration requires their unchanged entries.

Expected content edits are limited to:

- `docs/memory/CURRENT.md`
- `docs/memory/phases/phase-f6-method-a-production-promotion.md`
- `docs/memory/manifest.json`
- this task contract

All other first-pass reconciliation files must remain content-identical to the current
worktree.

No `src/`, `testing/`, launcher, profile, collector, configuration, JSON analysis
artifact, PDF, or runtime source changes.

---

## Scientific boundaries

Preserve exactly:

- accepted baseline production remains authoritative until F.6.4;
- Method A uses exactly the accepted F.4 parent-preserving correction;
- no independent `(t,phi)` child renormalization;
- Method B remains diagnostic/cross-check only and never changes pion weights;
- E.8 presentation never recomputes fits, factors, corrections, normalizations, or
  yields;
- F.6.3 alone constructs the actual parallel Method-A branch;
- accepted F.6.2 JSON/fingerprints remain frozen;
- source review does not establish farm/runtime validation.

Do not alter the statuses of F.1 through F.6.2.

---

## Local validation

Run the established repository memory tooling for:

- manifest regeneration;
- manifest check;
- memory integrity/health;
- memory bootstrap;
- `git diff --check`.

Also verify:

```bash
git diff --name-only
git status --short
```

Every changed/untracked path must be under `docs/memory/`, except an optional temporary
root-level `kaonlt_review.diff` used only for review.

No ROOT/PyROOT or farm execution.

---

## Review bundle

Produce a refreshed root-level:

`kaonlt_review.diff`

containing:

1. the complete tracked diff from committed HEAD;
2. complete `git diff --no-index -- /dev/null ...` content for every intended
   untracked memory file.

Do not add `kaonlt_review.diff` to Git.

---

## Acceptance criteria

PASS only if:

1. HEAD remains exactly `0ec29d4e1bb345eb37e8cca35b8b7e5cbe1b4d5b`;
2. the approved E.8.2 -> E.8.3 -> F.6.3 -> E.8.4 -> final E.8 -> F.6.4 roadmap remains
   unchanged;
3. the direct `w0 -> w0*C` and baseline-vs-reweighted pion-background plot requirements
   remain unchanged;
4. the complete baseline random/dummy/proton/pion/Fit1/Fit2/final-MM/yield audit remains
   unchanged;
5. the deleted durable F.6.1/F.6.2/F.6.4/uncertainty boundaries are restored in the
   detailed Phase-F record;
6. lifecycle status uses exactly `BLOCKED / DEFERRED`;
7. no analysis/runtime/test/profile source changes;
8. manifest and memory checks pass.

---

## Hard stop

Stop and report if satisfying this repair would require any scientific redesign,
production-source change, or work outside the allowed memory scope.

Codex must not commit, push, or run the farm.
