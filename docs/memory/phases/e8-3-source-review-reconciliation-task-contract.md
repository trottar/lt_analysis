# KaonLT E.8.3 — Source-Review Reconciliation

## 1. Objective

Record the independent ChatGPT actual-diff review PASS for the current cumulative
E.8.3 detached Method-A reweighting-audit candidate, including Fix.1, **before**
the user commits or pushes it.

This is a **memory-only reconciliation**. It must not change any source, tests,
physics, runtime behavior, presentation implementation, accepted artifacts, or
validation products.

After this reconciliation the durable state must be:

```text
E.8      — ACTIVE
E.8.2    — SOURCE REVIEWED
E.8.3    — SOURCE REVIEWED
F.6.3    — NEXT
E.8.4    — BLOCKED pending F.6.3
final E.8 — BLOCKED pending E.8.4 and its later runtime/visual gate
F.6.4    — BLOCKED pending completed F.6.3/E.8.4 production-impact evidence
```

The exact immediate action in `CURRENT.md` is still the user-controlled
commit/push of this reviewed cumulative set, followed by ChatGPT pushed-state
review. Only after pushed-state review passes does substantive F.6.3 work begin.

No farm run is authorized by this task.

---

## 2. Starting identity and intentional dirty worktree

Committed branch:

```text
test
```

Committed HEAD must remain exactly:

```text
91bb7809d27d84d7709a6600dbb0dc9ab514a458
```

The worktree is intentionally dirty with the cumulative accepted:

1. post-E.8.2 memory reconciliation;
2. E.8.3 implementation;
3. E.8.3 Fix.1 repair and tests;
4. this newly placed reconciliation contract.

Do not reset, stash, clean, checkout away, discard, or rewrite the cumulative
candidate.

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
git diff --name-only
```

Hard stop if the branch/committed HEAD differs or if unrelated changes are
present.

Do not pull, commit, push, update remote refs, or run the Jefferson Lab farm.

---

## 3. Required startup reading

Read the normal five-file startup sequence in order:

1. root `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read only:

- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/phases/phase-f6-method-a-production-promotion.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/e8-3-detached-method-a-reweighting-audit-task-contract.md`
- `docs/memory/phases/e8-3-detached-method-a-reweighting-audit-fix1-task-contract.md`
- `docs/memory/phases/e8-3-detached-method-a-reweighting-audit.md`
- this reconciliation contract

Do not broaden into unrelated phase history.

---

## 4. Independent review result to record

Independent ChatGPT review of the complete cumulative
`kaonlt_review(20260924-064150).diff` passed at the source/diff level.

The reviewed implementation establishes:

- E.8.3 remains detached and presentation-only.
- It consumes only byte-pinned accepted F.4/F.5/F.6.1/F.6.2 authorities.
- F.4/F.5/F.6.1 reader failures are fail-closed.
- F.6.2 remains consumed through the existing frozen E.8 authority.
- F.6.1 persisted baseline/Method-A/delta missing-mass arrays pass through
  without scientific reconstruction.
- F.5 persisted baseline/Method-A/delta `3 x 9` child aggregates pass through
  without child renormalization.
- Fix.1 preserves the persisted F.5 `event_counts` matrix and labels `EMPTY`
  only from persisted `event_count == 0`.
- Fix.1 renders the defined ratio with points only (`AP`), so omitted undefined
  denominator bins are not bridged/interpolated by a line.
- accepted F.5 parent closure is displayed without inventing a new relative
  metric.
- the E.8.3 route remains in the existing E.8.2 pair-safe PDF/manifest
  lifecycle.
- `src/cuts/rand_sub.py` retained the already-reviewed E.8.3 integration and
  Fix.1 did not alter it.
- Method B is not numerical.
- empirical residual Fit 1/Fit 2 remain absent.
- no F.6.3 production branch is constructed.
- baseline production, public yield/error calculation, cuts, binning,
  normalizations, pion/proton/random treatment, SIMC, efficiencies, acceptance,
  L/T separation, and cross sections remain outside E.8.3 ownership.

Codex-reported deterministic checks are recorded in the E.8.3 phase record.
They were **NOT RUN by ChatGPT**.

This source-review PASS does **not** establish:

- ROOT/PyROOT behavior;
- full `main.py` runtime integration;
- Jefferson Lab farm execution;
- procedure-PDF farm rendering;
- runtime acceptance;
- Method-A production promotion;
- F.6.3 implementation;
- final E.8 closure.

---

## 5. Allowed files

Only these existing files may be edited:

```text
docs/memory/CURRENT.md
docs/memory/decisions/e8-full-analysis-procedure-roadmap.md
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/roadmap/STATUS.md
docs/memory/phases/e8-3-detached-method-a-reweighting-audit.md
docs/memory/manifest.json
```

Add this contract:

```text
docs/memory/phases/e8-3-source-review-reconciliation-task-contract.md
```

No other file may change.

---

## 6. Frozen implementation

Everything outside the memory allowlist is frozen byte-for-byte for this task.

In particular, do not edit:

```text
src/**
testing/**
config/**
scripts/**
*.sh
```

Do not alter the already reviewed E.8.3 implementation or tests.

Before memory edits, save exact source/test patch snapshots:

```bash
git diff -- src/cuts/full_background_subtraction_plots.py src/cuts/rand_sub.py \
  > /tmp/e8_3_source_review_source_before.diff

git diff -- testing/test_e8_3_detached_method_a_reweighting_audit.py \
  > /tmp/e8_3_source_review_test_before.diff
```

After all memory edits, regenerate the same snapshots and compare with `cmp`.
Any difference is a hard stop.

---

## 7. Required status reconciliation

### 7.1 `docs/memory/CURRENT.md`

Keep E.8 `ACTIVE`.

Change E.8.3 from `ACTIVE` to:

```text
SOURCE REVIEWED
```

State narrowly that independent ChatGPT actual-diff review passed for the
cumulative E.8.3 + Fix.1 candidate.

Retain the source/runtime distinction:

- Codex-reported checks were `NOT RUN by ChatGPT`;
- no ROOT/PyROOT/full-analysis/farm/runtime acceptance is claimed.

Change F.6.3 from blocked-by-E.8.3 to dependency-ready:

```text
F.6.3 — NEXT
```

but make the **exact immediate Next Action**:

```text
NEXT — user-controlled commit/push of the reviewed cumulative
post-E.8.2-reconciliation + E.8.3/Fix.1 source/test/memory set
-> ChatGPT pushed-state review
-> F.6.3 source/runtime-path audit and standalone implementation contract.
```

Do not begin F.6.3 in this task.

Preserve:

- E.8.4 `BLOCKED` pending F.6.3;
- final E.8 `BLOCKED` pending E.8.4 and later runtime/visual gate;
- F.6.4 `BLOCKED` pending full production-impact evidence;
- deferred E.8.1 canonical-five expansion remains `DEFERRED`.

### 7.2 `docs/memory/roadmap/STATUS.md`

Record:

```text
E.8.2 — SOURCE REVIEWED
E.8.3 — SOURCE REVIEWED
F.6.3 — NEXT
E.8.4 — BLOCKED pending F.6.3
```

Do not change scientific ownership.

### 7.3 `docs/memory/phases/phase-f6-method-a-production-promotion.md`

Update only the E.8.3/F.6.3 dependency state:

- E.8.3 = `SOURCE REVIEWED`;
- F.6.3 = `NEXT`.

Preserve all F.1–F.6.2 closures and all production-promotion boundaries.

F.6.3 remains the sole owner of the future production change:

```text
w0_j -> w0_j * C_j
```

inside the parallel pion-template branch.

### 7.4 `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`

Update only current status/dependency wording:

- E.8.3 = `SOURCE REVIEWED`;
- F.6.3 = `NEXT`.

Preserve the approved order:

```text
E.8.2 -> E.8.3 -> F.6.3 -> E.8.4 -> final E.8 -> F.6.4
```

Do not alter the governing scientific/page requirements.

### 7.5 `docs/memory/phases/e8-3-detached-method-a-reweighting-audit.md`

Change phase status to:

```text
SOURCE REVIEWED
```

Add a concise independent-review closure section recording:

- complete cumulative diff reviewed;
- Fix.1 gap-safe ratio accepted;
- persisted empty-child identity accepted;
- focused contract coverage accepted;
- no source/test changes after this review;
- Codex checks `NOT RUN by ChatGPT`;
- no farm/runtime acceptance.

Update `## Next` to the same immediate commit/push -> pushed-state review ->
F.6.3 audit sequence.

Do not rewrite implementation chronology or Codex-run test evidence.

### 7.6 `docs/memory/manifest.json`

Regenerate after all versioned memory changes.

---

## 8. Validation

Use the repository-appropriate interpreter already established by the current
worktree.

Run:

```bash
<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/check_memory_health.py --root .
<PYTHON> -B tools/memory_bootstrap.py --root . --json
<PYTHON> -B -m unittest testing.test_memory_health -v
git -c core.safecrlf=false diff --check
```

Then verify implementation/test immutability:

```bash
git diff -- src/cuts/full_background_subtraction_plots.py src/cuts/rand_sub.py \
  > /tmp/e8_3_source_review_source_after.diff

git diff -- testing/test_e8_3_detached_method_a_reweighting_audit.py \
  > /tmp/e8_3_source_review_test_after.diff

cmp /tmp/e8_3_source_review_source_before.diff \
    /tmp/e8_3_source_review_source_after.diff

cmp /tmp/e8_3_source_review_test_before.diff \
    /tmp/e8_3_source_review_test_after.diff
```

Both `cmp` commands must return exit 0.

No analysis tests need rerunning because this task is memory-only and source/test
bytes must remain unchanged.

---

## 9. Complete cumulative review bundle

Refresh root:

```text
kaonlt_review.diff
```

It must remain a **complete cumulative diff from committed HEAD
`91bb7809d27d84d7709a6600dbb0dc9ab514a458`**.

Include:

1. complete tracked `git diff`;
2. `git diff --no-index /dev/null ...` sections for every intended untracked/new
   file, including:
   - post-E.8.2 reconciliation contract;
   - original E.8.3 contract;
   - E.8.3 Fix.1 contract;
   - this E.8.3 source-review reconciliation contract;
   - E.8.3 phase record;
   - focused E.8.3 test file.

Do not stage or commit `kaonlt_review.diff`.

---

## 10. Acceptance criteria

This reconciliation is ready for final ChatGPT pre-push verification only if:

- committed HEAD remains exact;
- only memory allowlist files changed during this task;
- source/test patch snapshots are byte-identical before/after;
- E.8.3 is consistently `SOURCE REVIEWED`;
- F.6.3 is consistently `NEXT`;
- no farm/runtime claim is introduced;
- Method A remains non-production;
- Method B remains diagnostic-only;
- empirical residual Fit 1/Fit 2 remain dormant/excluded;
- downstream blockers remain correct;
- manifest/memory checks pass;
- `git diff --check` passes;
- complete cumulative review bundle is refreshed;
- Codex does not commit, push, update remote refs, or run the farm.

---

## 11. Hard stop

Stop and report rather than broadening scope if:

- any source/test byte changes;
- a status transition would require changing scientific ownership;
- a memory check exposes a broader unrelated inconsistency;
- starting committed HEAD differs;
- unrelated dirty work appears;
- any evidence contradicts the independent E.8.3 source-review PASS.

Do not begin F.6.3 inside this reconciliation task.
