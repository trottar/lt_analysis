# KaonLT E.8.4.Fix.2 — Malformed-Sidecar Fail-Closed Repair

## Purpose

Repair one narrow source-review defect remaining after E.8.4.Fix.1.

The refreshed cumulative candidate in:

```text
kaonlt_review(20260928-184452).diff
```

correctly repairs the previous three findings:

1. semantic E.8.2 epsilon `low/high` is mapped explicitly to F.6.3
   `lowe/highe`;
2. wide E.8.2 MM geometry is kept distinct from the narrow F.6.3 analysis-MM
   geometry;
3. the common pion input remains visibly drawn with `B_pi_0` and `B_pi_A`.

Do not reopen those repairs.

This Fix.2 addresses only the remaining fail-closed contract gap: a malformed
but nominally `available=True` F.6.3 sidecar can currently raise an uncaught
Python exception inside the E.8.4 consumer instead of becoming an explicit
E.8.4-unavailable payload.

This is source development only. Do not run the Jefferson Lab farm.

---

# 1. Exact committed base and dirty candidate

Required branch:

```text
test
```

Required committed HEAD:

```text
d84f9610427cf38ccaaa8d6d8b14a8a995ecbac9
Add F6.3 parallel Method-A full procedure
```

Continue from the existing uncommitted cumulative E.8.4 + Fix.1 candidate.

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git log -1 --oneline
git status --short
```

The existing candidate paths must remain the same allowlisted E.8.4 paths, plus
this new contract:

```text
docs/memory/phases/e8-4-fix2-malformed-sidecar-fail-closed-task-contract.md
```

The old temporary review bundles may remain untracked/read-only. Do not stage
them.

Hard stop on any unrelated tracked modification or unexpected untracked work.
Do not reset, stash, clean, commit, push, or overwrite user work.

---

# 2. Mandatory startup reading

Read in the repository-standard order:

1. local `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read only the relevant E.8/F.6 records:

- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/phase-f6-method-a-production-promotion.md`
- `docs/memory/phases/e8-4-production-impact-audit-task-contract.md`
- `docs/memory/phases/e8-4-fix1-source-review-repair-task-contract.md`
- `docs/memory/phases/e8-4-production-impact-audit.md`
- this Fix.2 contract

Inspect the actual current source and focused test before editing.

---

# 3. Independent-review finding

The original E.8.4 contract explicitly requires:

```text
If the F.6.3 source is missing, unavailable, malformed, wrong-setting,
stale in geometry, or violates its flags, return an explicit unavailable
E.8.4 payload.
```

It also requires an unavailable optional F.6.3/E.8.4 branch not to abort or
invalidate the public baseline calculation.

The current Fix.1 candidate catches only `_E8PayloadError`, but performs raw
container coercions such as:

```python
source_window = tuple(source.get("lambda_integration_window") or ())
source_children = tuple(source.get("children") or ())
```

Therefore malformed available sidecars such as:

```python
{"lambda_integration_window": 1, ...}
```

or:

```python
{"children": 1, ...}
```

raise `TypeError: 'int' object is not iterable` rather than returning the
explicit unavailable E.8.4 payload.

This matters because
`finalize_full_background_subtraction_e8_2(...)` calls the E.8.4 builder
directly after the baseline E.8.2 presentation succeeds. The optional E.8.4
consumer must fail closed rather than letting malformed optional-sidecar
structure escape into the full-analysis runtime.

---

# 4. Required repair

Keep the existing E.8.4 algorithm and renderer unchanged except for the narrow
source-container validation needed here.

Before converting the F.6.3 Lambda window or child inventory to tuples:

1. require the Lambda-window object to be a non-string sequence;
2. require the child-inventory object to be a non-string sequence;
3. map malformed Lambda-window container structure to the existing literal
   unavailable reason:

```text
e8_4_lambda_window_invalid
```

4. map malformed child-inventory container structure to the existing literal
   unavailable reason:

```text
e8_4_child_inventory_invalid
```

Use existing `_E8PayloadError` fail-closed flow.

Do not add a blanket `except Exception:` around the full builder merely to hide
programming errors. Validate the external sidecar container fields explicitly.

Do not change:

- F.6.3 producer source;
- setting-token mapping;
- analysis-MM geometry logic;
- histogram detachment;
- baseline `Y0`/error identity;
- Method-A scalar calculations;
- renderer page content/order;
- pair-safe finalization;
- Method B;
- dormant Fit 1/Fit 2;
- public production outputs.

---

# 5. Allowed source/test scope

Source:

```text
src/cuts/full_background_subtraction_plots.py
```

Focused test:

```text
testing/test_e8_4_production_impact_audit.py
```

Warranted memory/history only:

```text
docs/memory/CURRENT.md
docs/memory/decisions/e8-full-analysis-procedure-roadmap.md
docs/memory/manifest.json
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/phases/e8-4-production-impact-audit.md
docs/memory/phases/e8-4-production-impact-audit-task-contract.md
docs/memory/phases/e8-4-fix1-source-review-repair-task-contract.md
docs/memory/phases/e8-4-fix2-malformed-sidecar-fail-closed-task-contract.md
docs/memory/roadmap/STATUS.md
```

No other path may change.

---

# 6. Explicitly frozen files

Remain byte-unchanged:

```text
src/main.py
src/binning/calculate_yield.py
src/cuts/rand_sub.py
src/cuts/pion_component_subtraction.py
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
src/utility/background_config.py
```

Do not change accepted F-stage calculators, Method-B code, collectors, profiles,
wrappers, launchers, yield/cross-section code, or production subtraction.

---

# 7. Focused deterministic coverage

Extend the E.8.4 focused suite with direct malformed-sidecar coverage.

At minimum prove:

1. `lambda_integration_window = 1` returns:
   ```text
   available=False
   reason=e8_4_lambda_window_invalid
   ```
   and does not raise.
2. `lambda_integration_window = {"low": 1.08, "high": 1.18}` is rejected as a
   malformed container, even though it is iterable.
3. `children = 1` returns:
   ```text
   available=False
   reason=e8_4_child_inventory_invalid
   ```
   and does not raise.
4. `children = {"unexpected": "mapping"}` is rejected as a malformed inventory.
5. The existing finalizer, supplied a malformed-but-schema-valid F.6.3 sidecar,
   can complete the baseline procedure transaction successfully with:
   ```text
   status=available
   e8_4_available=False
   e8_4_reason=<literal malformed-sidecar reason>
   ```
   when the E.8.4 unavailable page renders successfully.
6. Public/baseline objects are unchanged.
7. All Fix.1 setting-token, MM-geometry, common-input, rollback/recovery,
   zero-Y0, and static-ownership tests remain passing.

Do not create fake fallback Method-A output for malformed sources.

---

# 8. Required regression checks

Run:

```bash
python -B -m py_compile \
  src/cuts/full_background_subtraction_plots.py \
  testing/test_e8_4_production_impact_audit.py

python -B -m unittest testing.test_e8_4_production_impact_audit -v
python -B -m unittest testing.test_e8_2_baseline_stage_audit -v
python -B -m unittest testing.test_e8_3_detached_method_a_reweighting_audit -v
python -B -m unittest testing.test_f6_3_parallel_full_procedure_method_a -v
python -B -m unittest testing.test_full_background_subtraction_plots -v
python -B -m unittest testing.test_pion_hgcer_phase_e_runtime_contract -v
python -B -m unittest testing.test_pion_hgcer_phase_f_runtime_contract -v
```

Then:

```bash
python -B tools/update_memory_manifest.py --root . --write
python -B tools/update_memory_manifest.py --root . --check
python -B tools/check_memory_health.py --root .
python -B tools/memory_bootstrap.py --root . --json
python -B -m unittest testing.test_memory_health -v

git -c core.safecrlf=false diff --check
```

Record the actual interpreter/version.

These local checks are not ROOT/PyROOT/full-analysis/farm validation.

---

# 9. Memory/status update

Append a narrow:

```text
E.8.4.Fix.2 — malformed-sidecar fail-closed repair
```

section to:

```text
docs/memory/phases/e8-4-production-impact-audit.md
```

Record:

- the uncaught malformed-container review finding;
- explicit Lambda-window/children sequence validation;
- finalizer-level optional-sidecar isolation coverage;
- deterministic checks.

Keep status:

```text
E.8.4.Fix.2 — ACTIVE
```

pending independent ChatGPT review of the refreshed cumulative diff.

Do not self-promote to `SOURCE REVIEWED`,
`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`, or
`CLOSED / RUNTIME VALIDATED`.

Final E.8 and F.6.4 remain `BLOCKED`.

Regenerate `docs/memory/manifest.json`.

---

# 10. Diff audit and review bundle

Before stopping:

```bash
git status --short
git diff --stat
git diff -- src/cuts/full_background_subtraction_plots.py
git diff -- testing/test_e8_4_production_impact_audit.py
git diff -- docs/memory
git -c core.safecrlf=false diff --check
```

Confirm all frozen files remain byte-unchanged.

Create a new complete timestamped review bundle containing the cumulative diff
from committed HEAD `d84f961...` plus complete no-index sections for all intended
untracked/new files.

Do not stage merely to create the bundle.

---

# 11. Acceptance criteria

Fix.2 is locally complete only if:

1. committed base remains `d84f961...`;
2. only allowlisted E.8.4 paths changed;
3. malformed Lambda-window containers fail closed without Python exceptions;
4. malformed child-inventory containers fail closed without Python exceptions;
5. malformed optional E.8.4 source cannot abort an otherwise successful baseline
   finalization transaction;
6. no fallback/reconstruction of Method A occurs;
7. all Fix.1 repairs remain intact;
8. no frozen producer/runtime file changes;
9. required deterministic/regression/memory checks pass;
10. `git diff --check` passes;
11. memory remains `ACTIVE` pending independent source review;
12. no commit, push, farm run, final-E.8 closure, F.6.4 work, or Method-A
    promotion occurs.

---

# 12. Hard stop

Stop after the narrow Fix.2 implementation, deterministic checks, warranted
memory update, manifest regeneration, diff audit, and refreshed review bundle.

Do not commit.
Do not push.
Do not run the farm.
Do not begin final E.8 or F.6.4.

The exact NEXT is:

```text
independent ChatGPT review of the refreshed complete cumulative E.8.4 + Fix.1 + Fix.2 diff/source/runtime path
```
