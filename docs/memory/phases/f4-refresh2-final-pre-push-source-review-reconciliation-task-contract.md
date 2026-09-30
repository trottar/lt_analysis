# KaonLT F.4.Refresh.2 — Final Pre-Push Source-Review Reconciliation

## 1. Purpose

Perform the required memory/status reconciliation after independent ChatGPT
actual-diff/source-runtime-path review **PASSED** the repaired cumulative
F.4.Refresh.2 current-baseline Method-A candidate materializer.

Reviewed cumulative bundle:

```text
kaonlt_review(20260930-111446).diff
```

Required committed base:

```text
08f8c5be84ab7a54278e8c623eae5d9d3c43d938
```

Independent review established that:

- remote `test` remains exactly
  `08f8c5be84ab7a54278e8c623eae5d9d3c43d938`;
- the cumulative candidate remains scoped to the approved detached materializer,
  focused tests, warranted F.4.Refresh.1 evidence/status updates, and memory;
- no `src/` scientific/runtime source changed;
- the original F.4.Refresh.2 builder/scientific architecture remains intact;
- F.4.Refresh.2.Fix.1 changes only:
  - manifest labeling of provenance-bound F.2/F.3/F.4 stage fingerprints;
  - overwrite completion-marker safety;
  - focused regression tests;
  - the F.4.Refresh.2 phase note;
  - manifest bookkeeping;
  - the Fix.1 contract;
- the reviewed materializer uses exact canonical-five F.1 input identity,
  comparison SHA identity, accepted F.2/F.3/F.4 identity, public
  F.2 -> F.3 -> F.4 builders, deterministic candidate F.2/F.3 byte identities,
  exact provenance-excluded F.2/F.3 scientific equality, the diagnostic
  candidate-F.3 zero-head override, and exact F.4.Refresh.1 F.4
  comparison/summary reproduction;
- candidate artifacts retain distinct non-authoritative filenames and are
  written with the existing public F.2/F.3/F.4 writers;
- Fix.1 correctly distinguishes:
  - `representation_fingerprint` for F.2;
  - `map_fingerprint` for F.3;
  - `correction_fingerprint` for F.4;
  from the explicit `scientific_gate`;
- under overwrite, a prior complete manifest is moved out of its final path
  only after all temporary artifacts and the new manifest have passed their
  pre-publication gates;
- if overwrite publication fails after candidate replacement begins, the final
  manifest remains absent rather than falsely certifying a mixed set;
- if failure occurs before candidate publication begins, restoration of the old
  manifest is supported;
- accepted authority constants, historical accepted artifacts, production
  physics, F.5, F.6.3, E.8.4, and Method A promotion remain untouched.

Codex-reported deterministic checks in the reviewed bundle:

```text
py_compile materializer/test                        PASS
materializer focused suite                          13 tests OK, 0 skips
F.4.Refresh.1 comparator suite                       8 tests OK, 0 skips
F.2 representation suite                            14 tests OK, 0 skips
F.3 map suite                                       10 tests OK, 0 skips
F.4 parent-preserving correction suite               8 tests OK, 0 skips
F.6.3 parallel full-procedure suite                 21 tests OK, 0 skips
memory-health suite                                 35 tests OK, 0 skips
manifest write/check                                PASS
memory health                                       exit 0, CURRENT soft warning
memory bootstrap                                    exit 0
git diff --check                                    PASS
```

Independent ChatGPT review also `py_compile`d the two new Python files
successfully.

The unit-test suites above were **NOT RUN by ChatGPT**.

This task is memory/status reconciliation only. It must not alter the
independently reviewed materializer or tests.

---

## 2. Exact starting state

Required branch:

```text
test
```

Required committed HEAD:

```text
08f8c5be84ab7a54278e8c623eae5d9d3c43d938
```

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
```

Expected cumulative candidate paths before adding this reconciliation contract:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/f4-refresh1-current-baseline-authority-comparison.md
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/roadmap/STATUS.md
docs/memory/evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md
docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization-task-contract.md
docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization.md
docs/memory/phases/f4-refresh2-fix1-materialization-manifest-integrity-task-contract.md
testing/materialize_method_a_current_baseline_authority.py
testing/test_materialize_method_a_current_baseline_authority.py
```

Temporary `kaonlt_review*.diff` files may remain untracked. They must not be
staged, removed, rewritten, or included in the candidate.

Root `AGENTS.md` and `.codex/` remain local-only/untracked.

If committed HEAD differs or unrelated tracked changes exist, **STOP**.

Do not reset, stash, clean, discard, commit, push, run the farm, materialize
farm artifacts, mutate accepted authority, or run the full analysis.

---

## 3. Mandatory repository-memory startup

Read in this exact order:

1. root `AGENTS.md`;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read only task-relevant records:

- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/phase-f6-method-a-production-promotion.md`
- `docs/memory/phases/f4-refresh1-current-baseline-authority-comparison.md`
- `docs/memory/evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md`
- `docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization.md`
- `docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization-task-contract.md`
- `docs/memory/phases/f4-refresh2-fix1-materialization-manifest-integrity-task-contract.md`

Do not reopen adjacent architecture or scientific phases.

---

## 4. Files that must remain byte-identical

Do not edit:

```text
testing/materialize_method_a_current_baseline_authority.py
testing/test_materialize_method_a_current_baseline_authority.py
testing/compare_method_a_current_baseline_authority.py
testing/test_compare_method_a_current_baseline_authority.py
docs/memory/evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md
docs/memory/phases/f4-refresh1-current-baseline-authority-comparison.md
docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization-task-contract.md
docs/memory/phases/f4-refresh2-fix1-materialization-manifest-integrity-task-contract.md
```

Also freeze every scientific/runtime source, including all `src/` files,
`src/main.py`, and `run_Prod_Analysis.sh`.

If the reviewed materializer or tests must change, **STOP**. That requires a new
repair rather than reconciliation.

---

## 5. Allowed reconciliation changes

Only:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization.md
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/roadmap/STATUS.md
```

Create exactly one new reconciliation contract:

```text
docs/memory/phases/f4-refresh2-final-pre-push-source-review-reconciliation-task-contract.md
```

No other file may change.

---

## 6. Required status reconciliation

### 6.1 F.4.Refresh.2

Advance:

```text
F.4.Refresh.2 — ACTIVE
```

to:

```text
F.4.Refresh.2 — SOURCE REVIEWED
```

Record that independent ChatGPT review of:

```text
kaonlt_review(20260930-111446).diff
```

passed.

Record the reviewed boundary:

- detached pure-Python candidate materializer only;
- current canonical-five F.1 inputs are explicit and fail-closed;
- reviewed comparison JSON raw SHA is explicit and fail-closed;
- accepted F.2/F.3/F.4 raw identities are verified;
- candidate F.2/F.3 are rebuilt with public builders;
- candidate F.2/F.3 deterministic writer-byte SHA values reproduce the
  reviewed comparator;
- F.2/F.3 scientific equality is exact and provenance-excluded;
- candidate F.4 is rebuilt with the public F.4 builder using only the explicit
  in-memory zero-head candidate-F.3 diagnostic authority record;
- the complete reviewed F.4 comparison/summary is reproduced;
- candidate outputs are distinct and non-authoritative;
- public writers own F.2/F.3/F.4 output bytes;
- manifest stage identities are correctly labeled as
  representation/map/correction fingerprints;
- the scientific-equality result remains separate in `scientific_gate`;
- overwrite completion-marker behavior is fail-closed;
- no accepted authority changes;
- no scientific/runtime source changes;
- no production application;
- no Method-A promotion.

### 6.2 F.4.Refresh.2.Fix.1

Advance:

```text
F.4.Refresh.2.Fix.1 — ACTIVE
```

to:

```text
F.4.Refresh.2.Fix.1 — SOURCE REVIEWED
```

Record that Fix.1 owns only:

- correct manifest labeling of provenance-bound stage identities;
- separation of those identities from exact scientific equality;
- overwrite completion-marker safety;
- associated focused regressions.

### 6.3 Correct the local-check count

The pre-Fix.1 phase text still says:

```text
materializer 11
```

Update that durable local-check summary to the reviewed current result:

```text
materializer 13
```

with `0 skips`.

Do not alter historical test counts for unchanged suites.

### 6.4 Preserve all other statuses

Keep:

- F.1 through F.6.2, including historical accepted F.4/F.5 and F.6.2.Fix.5 —
  **CLOSED / RUNTIME VALIDATED**
- F.4.Refresh.1 — **CLOSED / RUNTIME VALIDATED** only for its detached
  comparator gate
- F.4.Refresh.1.Fix.1 — **CLOSED / RUNTIME VALIDATED** only as part of that
  detached comparator gate
- E.8.4.Fix.4 — **CLOSED / RUNTIME VALIDATED** only for the narrow
  cache-semantics gate
- E.8 — **ACTIVE**
- E.8.2 — **SOURCE REVIEWED**
- E.8.3 — **SOURCE REVIEWED**
- F.6.3 — **SOURCE REVIEWED**
- E.8.4 — **SOURCE REVIEWED**
- final E.8 — **BLOCKED**
- F.6.4 — **BLOCKED**
- lifecycle-hook dispatch — **BLOCKED / DEFERRED**

The fresh Left/lowe F.6.3/E.8.4 runtime path remains blocked pending the
current-baseline authority lineage. F.4.Refresh.2 source review does not accept
a refreshed F.4 authority.

---

## 7. Validation boundary

Record accurately:

- Codex-reported materializer suite: 13 tests OK, 0 skips;
- comparator: 8 tests OK, 0 skips;
- F.2: 14 tests OK, 0 skips;
- F.3: 10 tests OK, 0 skips;
- F.4: 8 tests OK, 0 skips;
- F.6.3: 21 tests OK, 0 skips;
- memory: 35 tests OK, 0 skips;
- manifest/check/bootstrap/diff checks passed as recorded;
- ChatGPT independently `py_compile`d the two new Python files;
- unit tests were **NOT RUN by ChatGPT**.

Do not claim:

- farm materializer execution;
- ROOT/PyROOT validation;
- full `main.py` validation;
- candidate-artifact runtime validation;
- refreshed F.4 authority acceptance;
- F.5 refresh;
- F.6.3 runtime acceptance;
- E.8.4 runtime acceptance;
- production acceptance;
- Method-A promotion.

---

## 8. CURRENT.md exact NEXT

Replace the pre-review NEXT with exactly:

```text
NEXT — user-controlled commit/push of the independently reviewed F.4.Refresh.2 current-baseline Method-A candidate materializer.
```

There must be exactly one ordinary NEXT.

This does not authorize Codex to commit or push.

---

## 9. Manifest and local reconciliation checks

Regenerate and verify:

```bash
python -B tools/update_memory_manifest.py --root . --write
python -B tools/update_memory_manifest.py --root . --check
python -B tools/check_memory_health.py --root .
python -B tools/memory_bootstrap.py --root . --json
python -B -m unittest testing.test_memory_health -v
git -c core.safecrlf=false diff --check
```

The existing `CURRENT.md` soft-size warning is acceptable only if memory-health
exits 0.

No source/materializer/comparator/F.2/F.3/F.4/F.6.3 test rerun is required for
this memory-only reconciliation. If Codex reruns deterministic tests, report
them accurately without upgrading the evidence boundary.

---

## 10. Required cumulative diff audit

Before stopping:

```bash
git status --short
git -c core.safecrlf=false diff --stat
git -c core.safecrlf=false diff --no-ext-diff
```

The cumulative candidate must contain exactly:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/f4-refresh1-current-baseline-authority-comparison.md
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/roadmap/STATUS.md
docs/memory/evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md
docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization-task-contract.md
docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization.md
docs/memory/phases/f4-refresh2-fix1-materialization-manifest-integrity-task-contract.md
docs/memory/phases/f4-refresh2-final-pre-push-source-review-reconciliation-task-contract.md
testing/materialize_method_a_current_baseline_authority.py
testing/test_materialize_method_a_current_baseline_authority.py
```

Temporary `kaonlt_review*.diff` files may exist but must remain untracked and
outside the candidate inventory.

No `git add -A`.

---

## 11. Required final pre-push review bundle

Create one fresh repository-root bundle:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

It must contain byte-faithfully:

1. branch;
2. committed HEAD;
3. `git status --short`;
4. cumulative `git diff --stat`;
5. complete cumulative tracked `git diff --no-ext-diff`;
6. complete `git diff --no-index /dev/null ...` for every intended untracked
   candidate file, including the new final-reconciliation contract;
7. exact reconciliation-check commands/results;
8. final cumulative changed-path inventory.

Do not stage the review bundle.

---

## 12. Acceptance criteria

The final reconciliation is acceptable only if:

1. committed HEAD remains
   `08f8c5be84ab7a54278e8c623eae5d9d3c43d938`;
2. materializer/test bytes remain identical to
   `kaonlt_review(20260930-111446).diff`;
3. direct F.4.Refresh.1 comparator evidence remains unchanged;
4. both prior F.4.Refresh.2 contracts remain unchanged;
5. every scientific/runtime source remains unchanged;
6. F.4.Refresh.2 becomes **SOURCE REVIEWED**;
7. F.4.Refresh.2.Fix.1 becomes **SOURCE REVIEWED**;
8. the phase local-check count is corrected from 11 to 13 materializer tests;
9. F.4.Refresh.1/Fix.1 retain only their narrow
   **CLOSED / RUNTIME VALIDATED** comparator closure;
10. historical F.1-F.6.2 closures remain unchanged;
11. F.6.3/E.8.4 remain **SOURCE REVIEWED** and runtime blocked;
12. tests are described with the correct ChatGPT/Codex evidence boundary;
13. no farm materialization result is invented;
14. CURRENT has exactly the required single NEXT;
15. manifest/memory checks pass with at most the known nonfatal CURRENT warning;
16. no source, farm, authority, commit, push, full-analysis, or production
    action occurs;
17. a fresh byte-faithful cumulative review bundle is created.

---

## 13. Hard stop

After memory/status reconciliation, manifest regeneration, deterministic memory
checks, cumulative diff audit, and creation of the fresh review bundle:

**STOP.**

Do not commit, push, create a farm validation profile, run the materializer on
the farm, package candidate artifacts, update accepted F.2/F.3/F.4/F.5/F.6.3
authority, rerun the full analysis, or promote Method A.

Return the fresh bundle for independent ChatGPT final pre-push reconciliation
review.
