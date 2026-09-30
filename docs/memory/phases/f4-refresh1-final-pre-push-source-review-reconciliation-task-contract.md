# KaonLT F.4.Refresh.1 — Final Pre-Push Source-Review Reconciliation

## 1. Purpose

Perform the memory/status reconciliation required after independent ChatGPT
actual-diff/source-runtime-path review **PASSED** the repaired cumulative
F.4.Refresh.1 current-baseline Method-A authority comparator candidate.

Reviewed cumulative bundle:

```text
kaonlt_review(20260930-090154).diff
```

Required committed base:

```text
a8bd4dc4e0990dd45cfade0079dc1ce40bdd36ce
```

Independent review established that:

- remote `test` is still exactly `a8bd4dc4e0990dd45cfade0079dc1ce40bdd36ce`;
- the cumulative candidate remains scoped to the approved F.4.Refresh.1
  diagnostic/test plus warranted memory/history;
- all scientific/runtime source remains unchanged;
- the original reviewed F.4.Refresh.1 task contract, direct evidence record,
  CURRENT diff, E.8.4.Fix.4 phase diff, Phase-F roadmap diff, and STATUS diff
  are byte-identical to the prior `kaonlt_review(20260930-000747).diff`;
- Fix.1 changes only:
  - list-recursive F.4 diagnostic deltas;
  - list traversal in `_numeric_deltas()`;
  - realistic list-shaped focused tests;
  - the F.4.Refresh.1 phase note;
  - manifest bookkeeping;
  - the new Fix.1 task contract;
- the Fix.1 task contract in the cumulative bundle is byte-identical to the
  supplied contract, SHA-256:
  `b7aa93644477ad3c16f9576e787080f96ba88f19ee4c9b8c45a598e2ec31b9ce`;
- real `source_diagnostics` and `canonical_phi_diagnostics` list structure is
  now recursively represented without silent truncation;
- numerical leaves expose accepted/candidate/absolute/relative values;
- unequal list lengths emit explicit missing-side entries;
- the public F.2 -> F.3 -> F.4 builder chain, provenance exclusions,
  candidate-F.3 in-memory authority override, whole-payload scientific
  matching, first-changed-stage precedence, CLI, schema, and write safety remain
  unchanged.

Codex-reported Fix.1 checks passed:

```text
py_compile comparator/test                         PASS
focused comparator suite                           8 tests OK
F.4 parent-preserving correction suite             8 tests OK
F.6.3 parallel full-procedure suite                21 tests OK
manifest write/check                               PASS
memory health                                      exit 0, CURRENT soft warning
memory bootstrap                                   exit 0
memory-health suite                                35 tests OK
git diff --check                                   PASS
```

The original candidate had additionally reported F.2 (14 tests OK) and F.3
(10 tests OK); those source sections were not changed by Fix.1.

All test results are **Codex-reported and NOT RUN by ChatGPT**.

This task is memory/status reconciliation only. It must not alter the
independently reviewed comparator/test candidate.

---

## 2. Exact starting state

Required branch:

```text
test
```

Required committed HEAD:

```text
a8bd4dc4e0990dd45cfade0079dc1ce40bdd36ce
```

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
```

Expected cumulative candidate paths:

```text
docs/memory/CURRENT.md
docs/memory/evidence/e8-4-fix4-left-lowe-runtime-closure-and-f4-baseline-divergence.md
docs/memory/manifest.json
docs/memory/phases/e8-4-fix4-alignment-cache-semantics.md
docs/memory/phases/f4-refresh1-current-baseline-authority-comparison-task-contract.md
docs/memory/phases/f4-refresh1-current-baseline-authority-comparison.md
docs/memory/phases/f4-refresh1-fix1-nested-f4-delta-task-contract.md
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/roadmap/STATUS.md
testing/compare_method_a_current_baseline_authority.py
testing/test_compare_method_a_current_baseline_authority.py
```

Temporary repository-root review bundles may remain untracked, including:

```text
kaonlt_review(20260930-000747).diff
kaonlt_review(20260930-090154).diff
```

They must not be staged, removed, rewritten, or treated as candidate files.

Root `AGENTS.md` and `.codex/` remain local-only/untracked.

If committed HEAD differs or unrelated tracked changes exist, **STOP**.

Do not reset, stash, clean, discard, commit, push, run the farm comparator,
rerun the full analysis, or mutate accepted artifacts.

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
- `docs/memory/phases/f4-refresh1-current-baseline-authority-comparison-task-contract.md`
- `docs/memory/phases/f4-refresh1-fix1-nested-f4-delta-task-contract.md`
- `docs/memory/evidence/e8-4-fix4-left-lowe-runtime-closure-and-f4-baseline-divergence.md`

Do not reopen adjacent architecture or scientific phases.

---

## 4. Files that must remain byte-identical

Do not edit:

```text
testing/compare_method_a_current_baseline_authority.py
testing/test_compare_method_a_current_baseline_authority.py
docs/memory/evidence/e8-4-fix4-left-lowe-runtime-closure-and-f4-baseline-divergence.md
docs/memory/phases/e8-4-fix4-alignment-cache-semantics.md
docs/memory/phases/f4-refresh1-current-baseline-authority-comparison-task-contract.md
docs/memory/phases/f4-refresh1-fix1-nested-f4-delta-task-contract.md
```

Also freeze every scientific/runtime source, including:

```text
src/cuts/pion_hgcer_method_a_acceptance_representation.py
src/cuts/pion_hgcer_method_a_acceptance_map.py
src/cuts/pion_hgcer_method_a_parent_preserving_correction.py
src/cuts/pion_hgcer_method_a_tphi_propagation.py
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
src/cuts/pion_component_fits.py
src/cuts/pion_hgcer_method_a_acceptance_contract.py
src/cuts/pion_hgcer_refinement_method_a.py
src/cuts/full_background_subtraction_plots.py
src/cuts/rand_sub.py
src/main.py
run_Prod_Analysis.sh
```

If any reviewed comparator/test byte must change, **STOP**. That requires a new
repair rather than reconciliation.

---

## 5. Allowed reconciliation changes

Only:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/f4-refresh1-current-baseline-authority-comparison.md
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/roadmap/STATUS.md
```

Create exactly one new reconciliation contract:

```text
docs/memory/phases/f4-refresh1-final-pre-push-source-review-reconciliation-task-contract.md
```

No other file may change.

---

## 6. Required status reconciliation

### 6.1 F.4.Refresh.1

Advance:

```text
F.4.Refresh.1 — ACTIVE
```

to:

```text
F.4.Refresh.1 — SOURCE REVIEWED
```

Record that independent ChatGPT review of:

```text
kaonlt_review(20260930-090154).diff
```

passed.

Record the reviewed boundary:

- diagnostic-only comparator;
- exact canonical-five current F.1 inventory;
- raw-byte SHA-256 provenance;
- existing public F.2/F.3/F.4 builders only;
- deterministic candidate F.2/F.3 serialized hashes;
- candidate F.3 authority override in-memory only;
- provenance-excluded exact scientific comparisons;
- recursive F.4 parent/child diagnostic deltas;
- no accepted authority mutation;
- no production/scientific source mutation;
- no Method-A promotion.

### 6.2 Fix.1

Record Fix.1 as source-reviewed within the F.4.Refresh.1 record:

```text
F.4.Refresh.1.Fix.1 — SOURCE REVIEWED
```

The Fix.1 boundary is diagnostic-report completeness only:

- list-recursive `_delta()`;
- explicit unmatched-list-side reporting;
- list traversal in `_numeric_deltas()`;
- realistic focused fixtures;
- no scientific-comparison or builder-chain change.

### 6.3 Preserve current project status

Keep:

- F.1 through F.6.2, including F.6.2.Fix.5 —
  **CLOSED / RUNTIME VALIDATED**
- E.8.4.Fix.4 persisted-alignment semantic-version repair —
  **CLOSED / RUNTIME VALIDATED** only for the narrow cache-semantics gate
- E.8 — **ACTIVE**
- E.8.2 — **SOURCE REVIEWED**
- E.8.3 — **SOURCE REVIEWED**
- F.6.3 — **SOURCE REVIEWED**
- E.8.4 — **SOURCE REVIEWED**
- final E.8 — **BLOCKED**
- F.6.4 — **BLOCKED**
- lifecycle-hook dispatch — **BLOCKED / DEFERRED**

The fresh Left/lowe F.6.3/E.8.4 runtime gate remains **BLOCKED** by:

```text
f6_3_f4_shared_reproduction_failed:f3_fingerprint_input_content_mismatch
```

until the comparator is pushed/pushed-state-reviewed and then run on the
approved farm artifacts.

Do not claim an F.2/F.3/F.4 current-baseline outcome yet. The comparator has not
run against farm artifacts.

---

## 7. Validation boundary

Record accurately:

- Codex-reported Fix.1 checks passed;
- comparator suite: 8 tests OK;
- F.4 suite: 8 tests OK;
- F.6.3 suite: 21 tests OK;
- memory suite: 35 tests OK;
- no skips in required Fix.1 suites;
- original unchanged F.2/F.3 candidate suites had reported 14 and 10 tests OK;
- all were **NOT RUN by ChatGPT**.

Do not claim:

- ROOT/PyROOT validation;
- farm comparator execution;
- full `main.py` validation;
- procedure-PDF acceptance;
- F.6.3 runtime acceptance;
- E.8.4 runtime acceptance;
- production acceptance;
- Method-A promotion;
- current-baseline F.4 authority acceptance.

---

## 8. CURRENT.md exact NEXT

Replace the pre-review NEXT with exactly:

```text
NEXT — user-controlled commit/push of the independently reviewed F.4.Refresh.1 current-baseline Method-A authority comparator.
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

The existing CURRENT soft-size warning is acceptable only if memory-health exits
0.

No source/comparator rerun is required for this memory-only reconciliation. If
Codex reruns deterministic tests, report them accurately without upgrading the
evidence boundary.

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
docs/memory/evidence/e8-4-fix4-left-lowe-runtime-closure-and-f4-baseline-divergence.md
docs/memory/manifest.json
docs/memory/phases/e8-4-fix4-alignment-cache-semantics.md
docs/memory/phases/f4-refresh1-current-baseline-authority-comparison-task-contract.md
docs/memory/phases/f4-refresh1-current-baseline-authority-comparison.md
docs/memory/phases/f4-refresh1-fix1-nested-f4-delta-task-contract.md
docs/memory/phases/f4-refresh1-final-pre-push-source-review-reconciliation-task-contract.md
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/roadmap/STATUS.md
testing/compare_method_a_current_baseline_authority.py
testing/test_compare_method_a_current_baseline_authority.py
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
7. exact reconciliation check commands/results;
8. final cumulative changed-path inventory.

Do not stage the review bundle.

---

## 12. Acceptance criteria

The final reconciliation is acceptable only if:

1. committed HEAD remains
   `a8bd4dc4e0990dd45cfade0079dc1ce40bdd36ce`;
2. comparator/test bytes remain identical to
   `kaonlt_review(20260930-090154).diff`;
3. direct evidence, Fix.4 phase, original task contract, and Fix.1 contract
   remain unchanged;
4. scientific/runtime source remains unchanged;
5. F.4.Refresh.1 becomes **SOURCE REVIEWED**;
6. F.4.Refresh.1.Fix.1 is recorded **SOURCE REVIEWED**;
7. Fix.4 retains only its narrow **CLOSED / RUNTIME VALIDATED** closure;
8. F.6.3/E.8.4 remain **SOURCE REVIEWED** and runtime blocked;
9. historical F.1-F.6.2 closures remain unchanged;
10. tests are described as Codex-reported and `NOT RUN by ChatGPT`;
11. no farm comparator result is invented;
12. CURRENT has exactly the required single NEXT;
13. manifest/memory checks pass with at most the known nonfatal CURRENT warning;
14. no source, farm, authority, commit, push, or production action occurs;
15. a fresh byte-faithful cumulative review bundle is created.

---

## 13. Hard stop

After memory/status reconciliation, manifest regeneration, deterministic memory
checks, cumulative diff audit, and creation of the fresh review bundle:

**STOP.**

Do not commit, push, run the farm comparator, rerun the full analysis, update
accepted F.2/F.3/F.4 authorities, or promote Method A.

Return the fresh bundle for independent ChatGPT final pre-push reconciliation
review.
