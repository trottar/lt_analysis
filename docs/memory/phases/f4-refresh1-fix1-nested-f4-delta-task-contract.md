# KaonLT F.4.Refresh.1.Fix.1 — Nested F.4 Diagnostic Delta Repair

## 1. Purpose

Repair one narrow source-review defect in the local F.4.Refresh.1
current-baseline Method-A authority comparator candidate.

Independent ChatGPT actual-diff/source-runtime-path review of:

```text
kaonlt_review(20260930-000747).diff
```

found that the candidate otherwise follows the approved F.4.Refresh.1
architecture, but its detailed F.4 comparison does not satisfy the contract for
the real nested shapes of:

```text
source_diagnostics
canonical_phi_diagnostics
```

The F.4.Refresh.1 contract requires, for every parent, accepted value, candidate
value, absolute difference, and relative difference where defined for those
diagnostics as well as the other listed parent metrics.

The comparator currently includes both fields in `PARENT_METRICS`, but
`_delta()` recurses only through dictionaries. Real F.4 source persists
`source_diagnostics` and `canonical_phi_diagnostics` as lists of dictionaries.
Therefore the current candidate treats each complete list as one atomic
equality object and does not emit nested numerical absolute/relative
differences for the signed source sums or canonical-phi child sums.

The focused synthetic test currently models both diagnostics as dictionaries,
so it does not exercise the real F.4 list structure.

This is a diagnostic-output completeness defect only. It does not invalidate:

- the F.2/F.3/F.4 public-builder orchestration;
- the canonical-five F.1 inventory;
- raw SHA provenance;
- deterministic F.2/F.3 candidate hashes;
- the in-memory candidate-F.3 authority override;
- F.2/F.3 scientific projections;
- whole-F.4 scientific matching;
- first-changed-stage ordering;
- output write safety;
- any scientific/runtime source;
- the narrow E.8.4.Fix.4 runtime closure.

Repair only the recursive delta/reporting layer and its focused tests.

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

The worktree must contain the cumulative F.4.Refresh.1 candidate represented by:

```text
kaonlt_review(20260930-000747).diff
```

Expected candidate paths before this repair:

```text
docs/memory/CURRENT.md
docs/memory/evidence/e8-4-fix4-left-lowe-runtime-closure-and-f4-baseline-divergence.md
docs/memory/manifest.json
docs/memory/phases/e8-4-fix4-alignment-cache-semantics.md
docs/memory/phases/f4-refresh1-current-baseline-authority-comparison-task-contract.md
docs/memory/phases/f4-refresh1-current-baseline-authority-comparison.md
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/roadmap/STATUS.md
testing/compare_method_a_current_baseline_authority.py
testing/test_compare_method_a_current_baseline_authority.py
```

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
```

Hard requirements:

- committed HEAD remains exactly `a8bd4dc4e0990dd45cfade0079dc1ce40bdd36ce`;
- do not reset, stash, clean, discard, commit, push, or run the farm;
- temporary `kaonlt_review*.diff` files remain temporary and untracked;
- root `AGENTS.md` and `.codex/` remain local-only/untracked;
- if the cumulative worktree differs materially from the reviewed candidate,
  STOP and report the blocker.

---

## 3. Mandatory repository-memory startup

Read in order:

1. root `AGENTS.md`;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read only task-relevant records/source:

- `docs/memory/phases/f4-refresh1-current-baseline-authority-comparison-task-contract.md`
- `docs/memory/phases/f4-refresh1-current-baseline-authority-comparison.md`
- `src/cuts/pion_hgcer_method_a_parent_preserving_correction.py`
- `testing/compare_method_a_current_baseline_authority.py`
- `testing/test_compare_method_a_current_baseline_authority.py`

---

## 4. Independent-review finding to repair

The real accepted/current F.4 parent payload uses:

```python
source_diagnostics = [
    {
        "source_label": ...,
        "event_count": ...,
        "baseline_signed_sum": ...,
        "adjusted_signed_sum": ...,
        "signed_delta": ...,
    },
    ...
]

canonical_phi_diagnostics = [
    {
        "phi_index": ...,
        "phi_low": ...,
        "phi_high": ...,
        "event_count": ...,
        "baseline_signed_sum": ...,
        "adjusted_signed_sum": ...,
        "signed_delta": ...,
    },
    ...
]
```

The reviewed comparator currently behaves schematically as:

```python
if numeric:
    return numeric_delta
if dict:
    recurse
return {
    "accepted": accepted,
    "candidate": candidate,
    "match": accepted == candidate,
}
```

Thus lists are atomic.

Required behavior:

- recursively preserve list order;
- recursively compare every common list element;
- emit numerical accepted/candidate/absolute/relative values at nested numerical
  leaves;
- preserve non-numeric accepted/candidate/match information;
- represent list-length or missing-element differences explicitly rather than
  truncating them;
- retain exact deterministic JSON-safe output.

The existing real F.4 lists are deterministically ordered by source identity and
canonical-phi identity, so index-preserving recursion is appropriate. Do not add
sorting or semantic re-keying that changes their source-defined ordering.

---

## 5. Allowed substantive edits

Only:

```text
testing/compare_method_a_current_baseline_authority.py
testing/test_compare_method_a_current_baseline_authority.py
```

Allowed memory edit:

```text
docs/memory/phases/f4-refresh1-current-baseline-authority-comparison.md
docs/memory/manifest.json
```

Create exactly one new repair contract:

```text
docs/memory/phases/f4-refresh1-fix1-nested-f4-delta-task-contract.md
```

No other file may change during this repair.

In particular, keep byte-identical:

```text
docs/memory/CURRENT.md
docs/memory/evidence/e8-4-fix4-left-lowe-runtime-closure-and-f4-baseline-divergence.md
docs/memory/phases/e8-4-fix4-alignment-cache-semantics.md
docs/memory/phases/f4-refresh1-current-baseline-authority-comparison-task-contract.md
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/roadmap/STATUS.md
```

and every scientific/runtime source file.

---

## 6. Required code repair

### 6.1 `_delta()`

Extend `_delta()` so lists recurse deterministically.

For equal-length lists, every output element must be the recursive delta for the
same accepted/candidate index.

For unequal-length lists, do not use `zip()` truncation. Preserve every index
through the maximum length and explicitly mark a missing accepted or candidate
side for unmatched elements.

A valid design is a list whose entries are recursive deltas for common indexes
and explicit missing-side records for unmatched indexes. Equivalent
deterministic structure is acceptable.

Do not alter numeric delta semantics:

```text
absolute_difference = abs(candidate - accepted)
relative_difference = absolute_difference / abs(accepted)
```

with `relative_difference = null` when the accepted value is zero.

### 6.2 `_numeric_deltas()`

Make `_numeric_deltas()` recurse through both dictionaries and lists.

This preserves correct maxima behavior if a future required maximum is drawn
from a nested list structure. It must not change the existing maxima definitions
or add new scientific thresholds.

### 6.3 No other comparator behavior change

Freeze:

- `ALIASES`;
- F.2/F.3/F.4 provenance exclusion sets;
- required scientific top-level fields;
- F.1 input hashing;
- public builder calls;
- candidate F.2/F.3 serialization and hashes;
- candidate F.3 authority override;
- whole-payload scientific matching;
- first mismatch behavior;
- parent inventory;
- `PARENT_METRICS`;
- maxima names and definitions;
- first-changed-stage precedence;
- CLI;
- overwrite behavior;
- schema version;
- exit behavior.

Do not introduce new acceptance thresholds or promotion logic.

---

## 7. Required focused tests

Update the parent fixture so:

```text
source_diagnostics
canonical_phi_diagnostics
```

use realistic lists of dictionaries rather than dictionaries.

Add explicit regression assertions that:

1. a changed `source_diagnostics[*].baseline_signed_sum` emits:
   - accepted value;
   - candidate value;
   - absolute difference;
   - relative difference;

2. a changed `canonical_phi_diagnostics[*].adjusted_signed_sum` emits the same
   four numerical fields;

3. unchanged identity fields remain represented without losing structure;

4. unequal list lengths do not truncate and explicitly identify the missing
   side;

5. `_numeric_deltas()` reaches numerical delta leaves nested inside lists.

Preserve all existing focused orchestration/write-safety tests.

---

## 8. Local validation

Run:

```bash
python -B -m py_compile \
  testing/compare_method_a_current_baseline_authority.py \
  testing/test_compare_method_a_current_baseline_authority.py

python -B -m unittest \
  testing.test_compare_method_a_current_baseline_authority -v

python -B -m unittest \
  testing.test_pion_hgcer_method_a_parent_preserving_correction -v

python -B -m unittest \
  testing.test_f6_3_parallel_full_procedure_method_a -v

python -B tools/update_memory_manifest.py --root . --write
python -B tools/update_memory_manifest.py --root . --check
python -B tools/check_memory_health.py --root .
python -B tools/memory_bootstrap.py --root . --json
python -B -m unittest testing.test_memory_health -v

git -c core.safecrlf=false diff --check
```

The known CURRENT soft-size warning is acceptable only with exit code 0.

No ROOT/PyROOT, farm comparator, full analysis, or production run.

---

## 9. Required memory update

Append to:

```text
docs/memory/phases/f4-refresh1-current-baseline-authority-comparison.md
```

that independent ChatGPT review of:

```text
kaonlt_review(20260930-000747).diff
```

found one narrow repair:

```text
nested F.4 source_diagnostics/canonical_phi_diagnostics were emitted as atomic
list equality objects rather than recursively reporting numerical deltas
```

Record that Fix.1 addresses only diagnostic report completeness.

Keep F.4.Refresh.1:

```text
ACTIVE
```

pending independent review of the repaired cumulative candidate.

Do not change the existing exact CURRENT NEXT.

Regenerate:

```text
docs/memory/manifest.json
```

---

## 10. Diff audit

Before stopping:

```bash
git status --short
git -c core.safecrlf=false diff --stat
git -c core.safecrlf=false diff --no-ext-diff
```

Audit every intended new/untracked file with `git diff --no-index /dev/null ... || true`.

The repaired cumulative candidate must contain the original F.4.Refresh.1
candidate plus exactly this repair contract and the allowed comparator/test/phase
repair.

No unrelated path may change.

---

## 11. Fresh cumulative review bundle

Create:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

at repository root.

It must contain byte-faithfully:

1. branch;
2. committed HEAD;
3. `git status --short`;
4. complete cumulative `git diff --stat`;
5. complete cumulative tracked `git diff --no-ext-diff`;
6. `git diff --no-index /dev/null ...` for every intended untracked candidate
   file, including:
   - original F.4.Refresh.1 task contract;
   - original F.4.Refresh.1 phase record;
   - direct Fix.4/F.4 divergence evidence;
   - comparator;
   - focused test;
   - this Fix.1 contract;
7. exact validation commands/results;
8. final changed-path inventory.

Do not stage the review bundle.

---

## 12. Acceptance criteria

The repair passes only if:

1. committed HEAD remains `a8bd4dc4e0990dd45cfade0079dc1ce40bdd36ce`;
2. scientific/runtime source remains unchanged;
3. original F.4.Refresh.1 architecture remains unchanged;
4. `_delta()` recursively represents list contents;
5. no unequal-length list is silently truncated;
6. source-diagnostic numerical leaves expose accepted/candidate/absolute/relative
   values;
7. canonical-phi numerical leaves expose the same;
8. `_numeric_deltas()` traverses lists;
9. realistic list-shaped tests cover the repair;
10. all original comparator focused tests remain passing;
11. existing F.4 and F.6.3 regressions pass;
12. F.4.Refresh.1 remains `ACTIVE`;
13. E.8.4.Fix.4 remains `CLOSED / RUNTIME VALIDATED` only for its narrow cache
    semantics gate;
14. F.6.3 and E.8.4 remain `SOURCE REVIEWED` and runtime blocked;
15. no farm run, authority refresh, production change, commit, or push occurs;
16. manifest and memory checks pass;
17. fresh cumulative review bundle is produced.

---

## 13. Hard stop

After the narrow repair, tests, manifest regeneration, phase-record update, diff
audit, and fresh cumulative review bundle:

**STOP.**

Do not commit, push, run the farm comparator, rerun the full analysis, update
accepted F.2/F.3/F.4 authority, or promote Method A.

NEXT remains independent ChatGPT actual-diff/source-runtime-path review of the
repaired F.4.Refresh.1 cumulative candidate.
