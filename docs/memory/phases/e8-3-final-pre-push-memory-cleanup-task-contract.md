# KaonLT E.8.3 Final Pre-Push Memory Cleanup — Task Contract

## Objective

Perform one final **memory-only** cleanup before committing/pushing the already
independently source-reviewed E.8.3 + Fix.1 cumulative candidate.

Two stale current-state statements remain after the E.8.3 source-review
reconciliation:

1. `docs/memory/phases/e8-2-baseline-stage-audit.md` still says the next action is
   the E.8.3 artifact/source-path audit and that F.6.3 is blocked pending
   source-reviewed E.8.3.
2. `docs/memory/phases/phase-f6-method-a-production-promotion.md` still says
   E.8.3 is the active successor, even though E.8.3 is now SOURCE REVIEWED and
   F.6.3 is NEXT.

Repair only those stale statements. Do not change source, tests, physics,
runtime behavior, E.8.3 implementation, accepted artifacts, or validation
products.

After this task, all current-state records must agree:

```text
E.8      — ACTIVE
E.8.2    — SOURCE REVIEWED
E.8.3    — SOURCE REVIEWED
F.6.3    — NEXT
E.8.4    — BLOCKED pending F.6.3
final E.8 — BLOCKED pending E.8.4 + later runtime/visual gate
F.6.4    — BLOCKED pending completed F.6.3/E.8.4 evidence
```

The exact immediate action remains:

```text
user-controlled commit/push of the reviewed cumulative set
-> ChatGPT pushed-state review
-> F.6.3 source/runtime-path audit and standalone implementation contract
```

No farm run and no F.6.3 implementation are authorized here.

---

## Exact starting identity

Branch:

```text
test
```

Committed HEAD must still be:

```text
91bb7809d27d84d7709a6600dbb0dc9ab514a458
```

The worktree is intentionally dirty with the accepted cumulative
post-E.8.2-reconciliation + E.8.3 + Fix.1 + source-review reconciliation.

Do not reset, stash, clean, discard, or rewrite any existing cumulative work.

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
git diff --name-only
```

Hard stop on a different branch/HEAD or unrelated dirty work.

---

## Required startup reading

Read the normal five-file startup sequence:

1. root `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read:

- `docs/memory/roadmap/STATUS.md`
- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/phases/phase-f6-method-a-production-promotion.md`
- `docs/memory/phases/e8-2-baseline-stage-audit.md`
- `docs/memory/phases/e8-3-detached-method-a-reweighting-audit.md`
- `docs/memory/phases/e8-3-source-review-reconciliation-task-contract.md`
- this contract

---

## Allowed files

Only these existing files may be edited:

```text
docs/memory/phases/e8-2-baseline-stage-audit.md
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/manifest.json
```

Add this contract:

```text
docs/memory/phases/e8-3-final-pre-push-memory-cleanup-task-contract.md
```

Everything else is frozen.

In particular, do not edit:

```text
src/**
testing/**
docs/memory/CURRENT.md
docs/memory/roadmap/STATUS.md
docs/memory/decisions/e8-full-analysis-procedure-roadmap.md
docs/memory/phases/e8-3-detached-method-a-reweighting-audit.md
```

Those records are already correct and reviewed.

---

## Required edits

### 1. `docs/memory/phases/e8-2-baseline-stage-audit.md`

Preserve the full E.8.2 implementation/source-review history.

Replace the stale current `## Next` content that says:

```text
NEXT — audit the accepted F.4/F.5/F.6.1/F.6.2 Method-A artifact/source path
needed by E.8.3 ...
...
F.6.3 remains BLOCKED pending source-reviewed E.8.3
```

with a concise historical handoff consistent with the current state:

- E.8.2 remains SOURCE REVIEWED only;
- E.8.3 has since completed independent source review;
- F.6.3 is now the dependency NEXT;
- the immediate repository action is the user-controlled commit/push of the
  reviewed cumulative candidate, then ChatGPT pushed-state review;
- no E.8.2/E.8.3 farm/runtime acceptance is implied.

Do not rewrite earlier E.8.2 chronology.

### 2. `docs/memory/phases/phase-f6-method-a-production-promotion.md`

Preserve:

```text
E.8.2 — SOURCE REVIEWED
E.8.3 — SOURCE REVIEWED
F.6.3 — NEXT
```

Replace the stale sentence equivalent to:

```text
E.8.2 is SOURCE REVIEWED; its active successor is E.8.3 ...
```

with wording equivalent to:

```text
E.8.2 and E.8.3 are SOURCE REVIEWED; F.6.3 is the current successor in the
approved E.8.2 -> E.8.3 -> F.6.3 -> E.8.4 sequence.
```

Preserve all scientific ownership:

- F.6.3 alone owns the future production branch;
- only `w0_j -> w0_j * C_j` changes in that branch;
- Method B remains diagnostic only;
- empirical residual Fit 1/Fit 2 remain dormant;
- no F.1-F.6.2 closure is downgraded.

### 3. `docs/memory/manifest.json`

Regenerate after versioned memory changes.

---

## Source/test immutability

Before editing, snapshot the complete cumulative implementation/test patches:

```bash
git diff -- src/cuts/full_background_subtraction_plots.py src/cuts/rand_sub.py \
  > /tmp/e8_3_final_cleanup_source_before.diff

git diff -- testing/test_e8_3_detached_method_a_reweighting_audit.py \
  > /tmp/e8_3_final_cleanup_test_before.diff
```

After memory edits, repeat and require:

```bash
git diff -- src/cuts/full_background_subtraction_plots.py src/cuts/rand_sub.py \
  > /tmp/e8_3_final_cleanup_source_after.diff

git diff -- testing/test_e8_3_detached_method_a_reweighting_audit.py \
  > /tmp/e8_3_final_cleanup_test_after.diff

cmp /tmp/e8_3_final_cleanup_source_before.diff \
    /tmp/e8_3_final_cleanup_source_after.diff

cmp /tmp/e8_3_final_cleanup_test_before.diff \
    /tmp/e8_3_final_cleanup_test_after.diff
```

Both `cmp` commands must exit 0.

---

## Validation

Run:

```bash
<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/check_memory_health.py --root .
<PYTHON> -B tools/memory_bootstrap.py --root . --json
<PYTHON> -B -m unittest testing.test_memory_health -v
git -c core.safecrlf=false diff --check
```

No analysis/unit suite rerun is required because source/test bytes must remain
unchanged.

---

## Complete cumulative review bundle

Refresh root:

```text
kaonlt_review.diff
```

It must remain a complete cumulative diff from committed HEAD
`91bb7809d27d84d7709a6600dbb0dc9ab514a458`.

Include every intended new/untracked file, including this final cleanup
contract.

Do not stage or commit `kaonlt_review.diff`.

---

## Acceptance criteria

Pass only if:

- committed HEAD remains exact;
- only the three memory/manifest files plus this new contract change in this
  task;
- source/test cumulative patch snapshots are byte-identical before/after;
- no current-state record still says E.8.3 is NEXT/ACTIVE/BLOCKED;
- no current-state record still says F.6.3 is blocked pending E.8.3;
- no current-state record still says E.8.3 is the active successor;
- E.8.2 phase `## Next` no longer points to an already-completed E.8.3 audit;
- CURRENT/STATUS/E.8 roadmap/E.8.3 phase record remain untouched;
- manifest/memory checks pass;
- `git diff --check` passes;
- Codex does not commit, push, update remote refs, run the farm, or begin F.6.3.

---

## Hard stop

Stop rather than broaden scope if:

- any source/test bytes change;
- another current-state inconsistency is found outside the allowlist;
- a scientific ownership change appears necessary;
- committed HEAD differs;
- unrelated dirty work appears.

Do not begin F.6.3 in this task.
