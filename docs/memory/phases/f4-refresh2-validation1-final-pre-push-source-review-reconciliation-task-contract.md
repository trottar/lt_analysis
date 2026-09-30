# KaonLT F.4.Refresh.2.Validation.1 — Final Pre-Push Source-Review Reconciliation

## 1. Purpose

Perform the required memory/status reconciliation after independent ChatGPT
actual-diff/source-provenance review **PASSED** the cumulative
F.4.Refresh.2.Validation.1 farm materialization bundle/profile candidate.

Reviewed cumulative bundle:

```text
kaonlt_review(20260930-152818).diff
```

Required committed base:

```text
141a3d04f9e5d07be21dba14e0e63212c3990bf1
```

Independent review established that:

- remote `test` remains exactly
  `141a3d04f9e5d07be21dba14e0e63212c3990bf1`;
- the candidate changes only the approved generic bundle profile, focused
  profile test, task/phase records, and warranted durable memory;
- no materializer, comparator, collector, scientific, runtime, or production
  source changed;
- the new profile uses the existing
  `pion_hgcer_validation_bundle_profile/v4` schema and
  `generic_artifacts` collection mode;
- the top-level setting inventory remains exactly the canonical five:
  Left/lowe, Left/highe, Center/lowe, Center/highe, Right/highe;
- the profile declares no setting-scoped artifacts;
- the profile declares exactly five required global JSON artifacts:
  the byte-faithful F.4.Refresh.1 comparison-input copy, candidate F.2,
  candidate F.3, candidate F.4, and the F.4.Refresh.2 materialization manifest;
- source provenance pins the independently reviewed pushed materializer source
  `141a3d04f9e5d07be21dba14e0e63212c3990bf1`;
- only the new profile and focused profile test are allowed committed analysis
  files after that required source, with `docs/memory/` remaining the only
  allowed non-analysis prefix;
- focused synthetic coverage verifies exact declarations, complete atomic
  packaging, no undeclared payload inclusion, missing/invalid required JSON
  behavior, source-ancestry failure, unexpected materializer/scientific
  committed-file failure, allowed profile/test plus memory provenance, and
  output-ZIP no-overwrite behavior;
- the existing generic collector remains unchanged;
- the existing F.4.Refresh.2 materializer remains unchanged;
- no farm materialization, bundle ZIP, accepted-authority update, F.5/F.6.3/
  E.8.4 authority re-pin, full-analysis run, or Method-A promotion occurred.

The exact independently reviewed new-file SHA-256 values are:

```text
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
8490f5192d099853e1ecd37be1ff561601ee3e39526e975c65b2a9deb3714227

testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py
b43bd29bb655208ae79bb303755f0d5e17a48d838b3f5680f5d5e8e0a07623a7

docs/memory/phases/f4-refresh2-validation1-farm-materialization-bundle-profile-task-contract.md
907be106d4b274e0a8254726736e00552be4adf0c4bbf42ca4c4e22c1b290986

docs/memory/phases/f4-refresh2-validation1-farm-materialization-bundle-profile.md
a0ffb3c5aad0bb19f1ec9b8959c2466a905d57f3fefce8cc67e839337fdab1ac
```

The task-contract SHA above is byte-identical to the contract supplied by
ChatGPT.

Codex-reported local checks in the reviewed bundle:

```text
py_compile focused profile test                         PASS
focused F.4.Refresh.2 profile suite                     4 tests OK, 0 skips
generic validation-bundle collector suite              28 tests OK, 0 skips
F.4.Refresh.2 materializer suite                       13 tests OK, 0 skips
memory-health suite                                    35 tests OK, 0 skips
manifest write/check                                   PASS
memory health                                          exit 0, CURRENT soft warning
memory bootstrap                                       exit 0
git diff --check                                       PASS
ROOT/PyROOT and farm materialization/collection        NOT RUN
```

These unit suites were **NOT RUN by ChatGPT**.

This task is memory/status reconciliation only. It must not alter the
independently reviewed profile or test.

---

## 2. Exact starting state

Required branch:

```text
test
```

Required committed HEAD:

```text
141a3d04f9e5d07be21dba14e0e63212c3990bf1
```

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
git log -1 --oneline
```

Expected cumulative candidate before adding this reconciliation contract:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization.md
docs/memory/phases/f4-refresh2-validation1-farm-materialization-bundle-profile-task-contract.md
docs/memory/phases/f4-refresh2-validation1-farm-materialization-bundle-profile.md
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/roadmap/STATUS.md
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py
```

Temporary `kaonlt_review*.diff` files may remain untracked. They must remain
outside the candidate and must not be removed, rewritten, or staged.

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
- `docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization.md`
- `docs/memory/phases/f4-refresh2-validation1-farm-materialization-bundle-profile.md`
- `docs/memory/phases/f4-refresh2-validation1-farm-materialization-bundle-profile-task-contract.md`

Do not reopen adjacent architecture or scientific phases.

---

## 4. Files that must remain byte-identical

Do not edit:

```text
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py
docs/memory/phases/f4-refresh2-validation1-farm-materialization-bundle-profile-task-contract.md
testing/materialize_method_a_current_baseline_authority.py
testing/test_materialize_method_a_current_baseline_authority.py
testing/compare_method_a_current_baseline_authority.py
testing/test_compare_method_a_current_baseline_authority.py
testing/collect_pion_hgcer_validation_bundle.py
testing/test_collect_pion_hgcer_validation_bundle.py
```

Also freeze every scientific/runtime source, including all `src/` files,
`src/main.py`, and `run_Prod_Analysis.sh`.

The independently reviewed profile/test must retain their exact SHA-256 values
listed in Section 1.

If any reviewed substantive file must change, **STOP**. That requires a new
repair rather than reconciliation.

---

## 5. Allowed reconciliation changes

Only:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization.md
docs/memory/phases/f4-refresh2-validation1-farm-materialization-bundle-profile.md
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/roadmap/STATUS.md
```

Create exactly one new reconciliation contract:

```text
docs/memory/phases/f4-refresh2-validation1-final-pre-push-source-review-reconciliation-task-contract.md
```

No other file may change.

---

## 6. Required status reconciliation

### 6.1 F.4.Refresh.2.Validation.1

Advance:

```text
F.4.Refresh.2.Validation.1 — ACTIVE
```

to:

```text
F.4.Refresh.2.Validation.1 — SOURCE REVIEWED
```

Record that independent ChatGPT actual-diff/source-provenance review of:

```text
kaonlt_review(20260930-152818).diff
```

passed.

Record the exact reviewed source boundary:

- generic v4 collector reused unchanged;
- profile-only declaration of the five deterministic F.4.Refresh.2 outputs;
- no setting-scoped artifacts;
- canonical-five top-level setting inventory retained only because the existing
  collector schema requires it;
- required materializer source commit is
  `141a3d04f9e5d07be21dba14e0e63212c3990bf1`;
- only profile/test are allowed committed analysis files after that source,
  with `docs/memory/` the sole allowed non-analysis prefix;
- missing or malformed required artifacts fail bundle completeness;
- unexpected committed materializer/scientific/runtime changes fail provenance;
- no undeclared extra source artifact is silently packaged;
- existing output ZIP is not overwritten;
- `complete=true` will mean only package/source-check completeness, not
  scientific authority acceptance.

### 6.2 Preserve F.4.Refresh.2 status

Keep:

```text
F.4.Refresh.2       — SOURCE REVIEWED
F.4.Refresh.2.Fix.1 — SOURCE REVIEWED
```

Do not upgrade them to runtime validated. The materializer still has not run on
farm artifacts.

### 6.3 Preserve all upstream/downstream statuses

Keep:

- F.4.Refresh.1 and F.4.Refresh.1.Fix.1 —
  **CLOSED / RUNTIME VALIDATED**, detached comparison gate only;
- historical F.1 through F.6.2 —
  **CLOSED / RUNTIME VALIDATED**;
- E.8.4.Fix.4 —
  **CLOSED / RUNTIME VALIDATED** for its narrow cache-semantics gate only;
- E.8 —
  **ACTIVE**;
- E.8.2 —
  **SOURCE REVIEWED**;
- E.8.3 —
  **SOURCE REVIEWED**;
- F.6.3 —
  **SOURCE REVIEWED**;
- E.8.4 —
  **SOURCE REVIEWED**;
- final E.8 —
  **BLOCKED**;
- F.6.4 —
  **BLOCKED**;
- lifecycle-hook dispatch —
  **BLOCKED / DEFERRED**.

The fresh Left/lowe F.6.3/E.8.4 runtime path remains blocked pending the
current-baseline authority lineage. Source review of this bundle profile does
not accept refreshed F.4 authority.

---

## 7. Validation boundary

Record accurately:

```text
Codex-reported focused profile suite     4 tests OK, 0 skips
Codex-reported collector suite          28 tests OK, 0 skips
Codex-reported materializer suite       13 tests OK, 0 skips
Codex-reported memory suite             35 tests OK, 0 skips
manifest/check/bootstrap/diff checks     passed as recorded
ChatGPT unit tests                       NOT RUN
farm materialization                    NOT RUN
farm bundle collection                  NOT RUN
ROOT/PyROOT                              NOT RUN
full main.py                             NOT RUN
```

Do not claim:

- candidate-artifact runtime validation;
- refreshed F.4 authority acceptance;
- F.5 authority refresh;
- F.6.3 runtime acceptance;
- E.8.4 runtime acceptance;
- production acceptance;
- Method-A promotion.

---

## 8. CURRENT.md exact NEXT

Replace the pre-review NEXT with exactly:

```text
NEXT — user-controlled commit/push of the independently reviewed F.4.Refresh.2.Validation.1 farm materialization bundle/profile candidate.
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

If `python` is not the interpreter already established by the implementation
task, reuse that established interpreter instead.

The existing `CURRENT.md` soft-size warning is acceptable only if
memory-health exits 0.

No profile/collector/materializer suite rerun is required for this memory-only
reconciliation. If Codex reruns deterministic tests, report them accurately
without changing the evidence boundary.

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
docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization.md
docs/memory/phases/f4-refresh2-validation1-farm-materialization-bundle-profile-task-contract.md
docs/memory/phases/f4-refresh2-validation1-farm-materialization-bundle-profile.md
docs/memory/phases/f4-refresh2-validation1-final-pre-push-source-review-reconciliation-task-contract.md
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/roadmap/STATUS.md
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py
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
   `141a3d04f9e5d07be21dba14e0e63212c3990bf1`;
2. profile/test bytes remain identical to
   `kaonlt_review(20260930-152818).diff`;
3. the original Validation.1 task contract remains identical;
4. every materializer/comparator/collector/scientific/runtime source remains
   unchanged;
5. F.4.Refresh.2.Validation.1 becomes **SOURCE REVIEWED**;
6. F.4.Refresh.2/Fix.1 remain **SOURCE REVIEWED**;
7. historical F.1-F.6.2 closures remain unchanged;
8. F.6.3/E.8.4 remain **SOURCE REVIEWED** and runtime blocked;
9. tests are described with the correct ChatGPT/Codex evidence boundary;
10. no farm materialization result or bundle is invented;
11. CURRENT has exactly the required single NEXT;
12. manifest/memory checks pass with at most the known nonfatal CURRENT warning;
13. no source, farm, authority, commit, push, full-analysis, or production
    action occurs;
14. one fresh byte-faithful cumulative review bundle is created.

---

## 13. Hard stop

After memory/status reconciliation, manifest regeneration, deterministic memory
checks, cumulative diff audit, and creation of the fresh review bundle:

**STOP.**

Do not commit or push.

Do not run the farm materializer.

Do not run the farm collector.

Do not create the farm ZIP.

Do not update accepted F.2/F.3/F.4/F.5/F.6.3 authority.

Do not rerun the full analysis.

Do not promote Method A.

Return the fresh bundle for independent ChatGPT final pre-push reconciliation
review.
