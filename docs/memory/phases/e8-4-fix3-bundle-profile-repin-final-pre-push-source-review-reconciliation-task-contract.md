# KaonLT E.8.4.Fix.3 — Bundle-Profile Re-pin Final Pre-Push Source-Review Reconciliation

## Purpose

Perform the final **memory/status reconciliation only** after independent ChatGPT review of the complete E.8.4.Fix.3 validation-bundle profile re-pin candidate returned **PASS**.

This is the established final pre-push reconciliation gate.

The reviewed candidate remains based on:

```text
test
29d7b7f9635db899939efeb3508e941e994e8928
E8.4 Fix.3: repair pion alignment determinism
```

The independently reviewed profile re-pin changes only:

```text
testing/pion_hgcer_validation_bundle_profile_e8_1.json
testing/test_pion_hgcer_validation_bundle_profile_e8_1.py
```

plus warranted durable-memory records from the already reviewed candidate.

This task must not alter those substantive profile/test changes. It only records the independent PASS in durable memory and advances the exact NEXT to the user-controlled Git handoff.

---

## 1. Exact starting state

Required branch:

```text
test
```

Required committed HEAD:

```text
29d7b7f9635db899939efeb3508e941e994e8928
```

The worktree is expected to contain the already-reviewed uncommitted E.8.4.Fix.3 profile-repin candidate.

Before editing, verify:

```bash
git branch --show-current
git rev-parse HEAD
git log -1 --oneline
git status --short
```

Hard stop if committed HEAD differs.

Do not reset, stash, clean, commit, push, package, rerender, or run the farm.

Do not discard or rewrite the already-reviewed substantive profile/test changes.

---

## 2. Mandatory startup reading

Read in this order:

1. root `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read only the directly relevant records:

- `docs/memory/CODEX.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/e8-4-fix3-pion-alignment-determinism.md`
- `docs/memory/phases/e8-4-fix3-bundle-profile-repin-task-contract.md`
- `docs/memory/phases/e8-4-fix3-bundle-profile-repin.md`
- this reconciliation contract
- `testing/pion_hgcer_validation_bundle_profile_e8_1.json`
- `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`

Do not reopen scientific architecture.

---

## 3. Independent review result to record

Independent ChatGPT review of:

```text
kaonlt_review_e8_4_fix3_profile_repin.diff
```

returned:

```text
PASS
```

The review established, for source/provenance scope:

- starting/current committed HEAD was `29d7b7f9635db899939efeb3508e941e994e8928`;
- the profile re-pin changed `required_analysis_commit` only from
  `1aa1fd4184a6f8b20043e00ebb1ed3e9505a4935`
  to
  `29d7b7f9635db899939efeb3508e941e994e8928`;
- the focused test changed its exact reviewed-source expectation to the same Fix.3 commit;
- canonical-five setting inventory remained unchanged;
- artifact inventory remained unchanged;
- profile/test-only committed-range allowlist remained unchanged;
- `docs/memory/` remained the only non-analysis prefix;
- collector and `tcsh` wrapper remained unchanged;
- Fix.3 analysis source/test remained unchanged;
- F.6.3/E.8.4 scientific/runtime source remained unchanged;
- no farm, ROOT/PyROOT, procedure-PDF, production, Method-A promotion, or runtime claim was established.

Codex-reported deterministic tests remain **NOT RUN by ChatGPT**.

---

## 4. Required status reconciliation

### E.8.4.Fix.3

Keep:

```text
SOURCE REVIEWED
```

Do not upgrade it to runtime validated.

### E.8.4.Fix.3 validation-bundle profile re-pin

Advance from:

```text
ACTIVE
```

to:

```text
SOURCE REVIEWED
```

Record that independent ChatGPT actual-diff/source-provenance review returned PASS.

This status is source/provenance only.

It does not establish:

- ROOT/PyROOT validation;
- full `main.py` integration;
- procedure-PDF runtime correctness;
- farm validation;
- production-physics validation;
- Method-A promotion;
- final E.8 closure;
- F.6.4 closure.

### Preserved downstream state

Keep:

```text
E.8 — ACTIVE
E.8.4 — SOURCE REVIEWED
F.6.3 — SOURCE REVIEWED
Final E.8 — BLOCKED
F.6.4 — BLOCKED
```

The prior Left/lowe determinism observation remains blocker evidence, not runtime closure.

---

## 5. Exact NEXT after reconciliation

Update `docs/memory/CURRENT.md` so the sole exact NEXT is:

```text
NEXT — user-controlled commit/push of the independently reviewed E.8.4.Fix.3 validation-bundle profile re-pin.
```

Also preserve the project rule:

- do not authorize a farm command yet;
- after the user-controlled push, ChatGPT must review the pushed repository state;
- only after pushed-state review passes may the later farm gate be considered.

Do not write later-gate commands into CURRENT.

---

## 6. Allowed changes

Only the following files may change in this reconciliation task:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/e8-4-fix3-pion-alignment-determinism.md
docs/memory/phases/e8-4-fix3-bundle-profile-repin.md
docs/memory/phases/e8-4-fix3-bundle-profile-repin-final-pre-push-source-review-reconciliation-task-contract.md
docs/memory/roadmap/STATUS.md
```

No other path may change.

If another path appears necessary, stop and report the blocker.

---

## 7. Files that must remain byte-identical to the independently reviewed candidate

Do not alter:

```text
testing/pion_hgcer_validation_bundle_profile_e8_1.json
testing/test_pion_hgcer_validation_bundle_profile_e8_1.py
docs/memory/phases/e8-4-fix3-bundle-profile-repin-task-contract.md
src/cuts/pion_component_fits.py
testing/test_pion_component_dynamic_alignment.py
testing/collect_pion_hgcer_validation_bundle.py
testing/package_pion_hgcer_validation_bundle.tcsh
src/cuts/full_background_subtraction_plots.py
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
src/binning/calculate_yield.py
src/cuts/rand_sub.py
src/cuts/pion_component_subtraction.py
src/utility/background_config.py
src/main.py
run_Prod_Analysis.sh
```

No scientific or production behavior may change.

---

## 8. Required durable-memory edits

### `docs/memory/CURRENT.md`

Record:

- E.8.4.Fix.3 remains `SOURCE REVIEWED`;
- E.8.4.Fix.3 validation-bundle profile re-pin is now `SOURCE REVIEWED`;
- independent ChatGPT review PASS is recorded;
- no runtime/farm claim is added;
- exact NEXT is user-controlled commit/push of this reviewed candidate.

### `docs/memory/phases/e8-4-fix3-bundle-profile-repin.md`

Update status to:

```text
SOURCE REVIEWED
```

Add a concise final independent-review section recording:

- reviewed bundle name;
- independent ChatGPT result `PASS`;
- substantive profile/test change remained exact;
- collector/wrapper and scientific source remained frozen;
- Codex-reported tests were not run by ChatGPT;
- no farm/runtime claim;
- NEXT is user-controlled commit/push.

### `docs/memory/phases/e8-4-fix3-pion-alignment-determinism.md`

Only reconcile its stale downstream NEXT/reference so it no longer says the profile re-pin is awaiting independent review.

Keep Fix.3 itself `SOURCE REVIEWED`.

### `docs/memory/roadmap/STATUS.md`

Mark the Fix.3 profile re-pin `SOURCE REVIEWED` and preserve all downstream blocks/statuses.

Do not rewrite unrelated roadmap history.

### `docs/memory/manifest.json`

Regenerate after the final memory edits.

---

## 9. Deterministic reconciliation checks

Run:

```bash
python -B tools/update_memory_manifest.py --root . --write
python -B tools/update_memory_manifest.py --root . --check
python -B tools/check_memory_health.py --root .
python -B tools/memory_bootstrap.py --root . --json
python -B -m unittest testing.test_memory_health -v
git -c core.safecrlf=false diff --check
```

The known CURRENT soft-size warning is non-fatal only if the command exits 0. Record the actual result.

Do not rerun or alter the substantive profile/test implementation merely to perform this memory-only reconciliation.

---

## 10. Diff audit

Before stopping, verify:

```bash
git status --short
git diff --stat
git diff -- docs/memory/CURRENT.md
git diff -- docs/memory/phases/e8-4-fix3-pion-alignment-determinism.md
git diff -- docs/memory/phases/e8-4-fix3-bundle-profile-repin.md
git diff -- docs/memory/roadmap/STATUS.md
git diff -- testing/pion_hgcer_validation_bundle_profile_e8_1.json
git diff -- testing/test_pion_hgcer_validation_bundle_profile_e8_1.py
git diff -- src/cuts/pion_component_fits.py
git diff -- testing/test_pion_component_dynamic_alignment.py
git -c core.safecrlf=false diff --check
```

The substantive profile/test diffs must be byte-for-byte the same reviewed changes as before reconciliation.

The Fix.3 source/test diffs must remain empty.

---

## 11. Fresh cumulative review bundle

Create a **new timestamped cumulative review bundle** in the repository root.

Use the established `kaonlt_review(YYYYMMDD-HHMMSS).diff` naming pattern.

It must include:

1. starting branch and committed HEAD;
2. current branch and committed HEAD;
3. `git status --short`;
4. complete cumulative `git diff --stat`;
5. complete cumulative tracked diff;
6. complete `git diff --no-index /dev/null ...` sections for every intended new/untracked tracked candidate file, including:
   - `docs/memory/phases/e8-4-fix3-bundle-profile-repin-task-contract.md`
   - `docs/memory/phases/e8-4-fix3-bundle-profile-repin.md`
   - `docs/memory/phases/e8-4-fix3-bundle-profile-repin-final-pre-push-source-review-reconciliation-task-contract.md`
7. exact reconciliation-check commands and results;
8. final cumulative changed-path inventory.

Do not stage merely to create the review bundle.

The review bundle itself is temporary and must not be committed.

---

## 12. Acceptance criteria

This final pre-push reconciliation is ready for independent ChatGPT review only if:

1. committed HEAD remains exactly `29d7b7f9635db899939efeb3508e941e994e8928`;
2. the previously reviewed profile/test changes are unchanged;
3. only the allowed reconciliation memory paths plus this contract changed after the PASS;
4. E.8.4.Fix.3 remains `SOURCE REVIEWED`;
5. the Fix.3 bundle-profile re-pin is now `SOURCE REVIEWED`;
6. CURRENT exact NEXT is user-controlled commit/push;
7. E.8 remains `ACTIVE`;
8. E.8.4 and F.6.3 remain `SOURCE REVIEWED`;
9. final E.8 and F.6.4 remain `BLOCKED`;
10. no farm/runtime/production claim is added;
11. manifest is regenerated and checks pass;
12. memory health reports no hard failure;
13. `git diff --check` passes;
14. a fresh timestamped cumulative review bundle is produced;
15. no commit, push, packaging, rerender, or farm run occurs.

---

## 13. Hard stop

Stop after:

- final source-review memory/status reconciliation;
- manifest regeneration;
- memory/integrity checks;
- actual-diff audit;
- creation of one fresh timestamped cumulative review bundle.

Do not commit.
Do not push.
Do not package.
Do not rerender.
Do not run the farm.

The exact NEXT after this task is only:

```text
independent ChatGPT review of the final pre-push reconciliation bundle
```

After that independent PASS, the user alone performs the scoped commit/push.
