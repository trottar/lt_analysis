# KaonLT E.8.4.Fix.3 — Validation-Bundle Profile Re-pin

## Purpose

Re-pin the existing E.8.4 validation-bundle provenance gate to the pushed and independently source-reviewed E.8.4.Fix.3 source.

Required pushed analysis source:

`29d7b7f9635db899939efeb3508e941e994e8928`
`E8.4 Fix.3: repair pion alignment determinism`

Its parent is:

`9daca79051a4eb451a158dd7c34b7028c6b96456`
`Re-pin E8.4 validation bundle profile`

This task changes validation provenance only. It must not alter Fix.3 science, production physics, the collector, wrapper, renderer, Method A, Method B, yields, or cross sections.

## Exact starting state

Required branch: `test`

Required committed HEAD:

`29d7b7f9635db899939efeb3508e941e994e8928`

Before editing, verify:

```bash
git branch --show-current
git rev-parse HEAD
git log -1 --oneline
git status --short
```

Hard stop if committed HEAD differs. Do not reset, stash, clean, commit, push, package, or run the farm.

## Mandatory startup reading

Read in this order:

1. root `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read only the task-relevant records/source:

- `docs/memory/phases/e8-4-fix3-pion-alignment-determinism.md`
- `docs/memory/evidence/e8-4-left-lowe-pion-alignment-determinism-blocker.md`
- `docs/memory/phases/e8-4-bundle-profile-repin.md`
- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/phases/phase-f6-method-a-production-promotion.md`
- `docs/memory/roadmap/STATUS.md`
- this task contract
- `testing/pion_hgcer_validation_bundle_profile_e8_1.json`
- `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`
- `testing/collect_pion_hgcer_validation_bundle.py`
- `testing/package_pion_hgcer_validation_bundle.tcsh`

Do not reopen settled scientific architecture.

## Required profile change

In:

`testing/pion_hgcer_validation_bundle_profile_e8_1.json`

change only:

`source_identity.required_analysis_commit`

from:

`1aa1fd4184a6f8b20043e00ebb1ed3e9505a4935`

to:

`29d7b7f9635db899939efeb3508e941e994e8928`

Preserve exactly:

- schema version `pion_hgcer_validation_bundle_profile/v4`
- validation profile `phase_e8_1_full_background_procedure_pdf_farm_review/v1`
- collection mode `generic_artifacts`
- ordered five-setting inventory:
  - Left / lowe
  - Left / highe
  - Center / lowe
  - Center / highe
  - Right / highe
- global frozen F.6.2 acceptance-refinement JSON
- per-setting full-background procedure PDF
- per-setting full-background page-manifest JSON
- `allowed_committed_files` exactly:
  - `testing/pion_hgcer_validation_bundle_profile_e8_1.json`
  - `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`
- `allowed_non_analysis_path_prefixes` exactly:
  - `docs/memory/`

Do not add any Fix.3 source file to the allowlist. Fix.3 belongs in the required analysis commit itself.

## Required focused-test change

In:

`testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`

change only the reviewed-source constant/expectation from:

`1aa1fd4184a6f8b20043e00ebb1ed3e9505a4935`

to:

`29d7b7f9635db899939efeb3508e941e994e8928`

Preserve all existing coverage proving:

- exact five-setting inventory
- exact artifact inventory
- complete synthetic packaging
- missing/invalid required artifacts fail closed
- required source ancestry is mandatory
- only the profile/test pair and `docs/memory/` may follow the required analysis commit
- later analysis-source changes fail closed
- output collisions fail closed

Do not weaken provenance enforcement.

## Frozen implementation

Remain byte-unchanged from pushed HEAD `29d7b7f...`:

- `src/cuts/pion_component_fits.py`
- `testing/test_pion_component_dynamic_alignment.py`
- `testing/collect_pion_hgcer_validation_bundle.py`
- `testing/package_pion_hgcer_validation_bundle.tcsh`
- `src/cuts/full_background_subtraction_plots.py`
- `src/cuts/pion_hgcer_method_a_parallel_full_procedure.py`
- `src/binning/calculate_yield.py`
- `src/cuts/rand_sub.py`
- `src/cuts/pion_component_subtraction.py`
- `src/utility/background_config.py`
- `src/main.py`
- `run_Prod_Analysis.sh`

Also preserve all accepted production/scientific behavior, including random/dummy subtraction, slow-proton subtraction, baseline pion weights, canonical binning, SIMC normalization/templates, yields/errors, efficiencies, acceptance, L/T separation, frozen F.6.2 science, Method-B diagnostic-only status, and F.6.3 ownership of the private Method-A branch.

## Allowed changes

Substantive validation files only:

- `testing/pion_hgcer_validation_bundle_profile_e8_1.json`
- `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`

Warranted durable-memory/history only:

- `docs/memory/CURRENT.md`
- `docs/memory/manifest.json`
- `docs/memory/phases/e8-4-fix3-pion-alignment-determinism.md`
- `docs/memory/phases/e8-4-fix3-bundle-profile-repin-task-contract.md`
- `docs/memory/phases/e8-4-fix3-bundle-profile-repin.md`
- `docs/memory/roadmap/STATUS.md`

No other tracked path may change. If another source file appears necessary, stop and report the blocker rather than expanding scope.

## Status reconciliation

Record E.8.4.Fix.3 as:

`SOURCE REVIEWED`

with:

- independent ChatGPT actual-diff/source review PASS
- user-controlled pushed commit `29d7b7f9635db899939efeb3508e941e994e8928`
- independent ChatGPT pushed-state review PASS
- no ROOT/PyROOT, full-analysis, procedure-PDF, farm, or runtime validation
- existing Left/lowe observation remains blocker evidence, not closure

Create:

`docs/memory/phases/e8-4-fix3-bundle-profile-repin.md`

and keep this profile re-pin:

`ACTIVE`

pending independent ChatGPT review.

Preserve:

- E.8 `ACTIVE`
- E.8.4 `SOURCE REVIEWED`
- F.6.3 `SOURCE REVIEWED`
- final E.8 `BLOCKED`
- F.6.4 `BLOCKED`

The exact NEXT after Codex completes is only:

`independent ChatGPT review of the complete E.8.4.Fix.3 bundle-profile re-pin diff/review bundle`

Do not place or authorize a farm command yet.

## Deterministic local checks

Run:

```bash
python -B -m py_compile \
  testing/test_pion_hgcer_validation_bundle_profile_e8_1.py \
  testing/test_pion_component_dynamic_alignment.py

python -m json.tool \
  testing/pion_hgcer_validation_bundle_profile_e8_1.json >/dev/null

python -B -m unittest \
  testing.test_pion_hgcer_validation_bundle_profile_e8_1 -v

python -B -m unittest \
  testing.test_collect_pion_hgcer_validation_bundle -v

python -B -m unittest \
  testing.test_pion_component_dynamic_alignment -v
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

Report PyROOT skips honestly. The known `CURRENT.md` soft-size warning is non-fatal only if the command exits 0; record the actual result.

These checks do not establish farm/runtime validation.

## Diff audit

Before stopping:

```bash
git status --short
git diff --stat
git diff -- testing/pion_hgcer_validation_bundle_profile_e8_1.json
git diff -- testing/test_pion_hgcer_validation_bundle_profile_e8_1.py
git diff -- docs/memory/CURRENT.md
git diff -- docs/memory/phases/e8-4-fix3-pion-alignment-determinism.md
git diff -- docs/memory/phases/e8-4-fix3-bundle-profile-repin.md
git diff -- docs/memory/roadmap/STATUS.md
git diff -- src/cuts/pion_component_fits.py
git diff -- testing/test_pion_component_dynamic_alignment.py
git -c core.safecrlf=false diff --check
```

The final two diffs must be empty.

## Review bundle

Create exactly one temporary repository-root review artifact:

`kaonlt_review_e8_4_fix3_profile_repin.diff`

It must contain:

1. starting branch and HEAD
2. current HEAD
3. `git status --short`
4. `git diff --stat`
5. complete tracked diff
6. complete `git diff --no-index /dev/null <file>` sections for every intended new/untracked repository file
7. exact local test commands and results
8. final changed-path inventory

Do not stage merely to generate the review bundle. Do not add the review bundle to Git.

## Farm boundary

No Jefferson Lab farm command is authorized in this task.

Do not package a farm bundle.
Do not run `main.py`.
Do not rerender a procedure PDF.
Do not claim ROOT/PyROOT/full-runtime validation.

## Acceptance criteria

Ready for independent review only if:

1. starting committed HEAD is exactly `29d7b7f9635db899939efeb3508e941e994e8928`
2. profile required analysis commit is exactly `29d7b7f...`
3. focused test reviewed-source expectation is exactly `29d7b7f...`
4. profile/test-only committed-range allowlist is unchanged
5. `docs/memory/` remains the only non-analysis prefix
6. five-setting inventory is unchanged
7. artifact inventory is unchanged
8. collector and wrapper are unchanged
9. Fix.3 source and focused alignment test are unchanged
10. F.6.3/E.8.4 source is unchanged
11. required local checks pass, with skips reported
12. memory integrity passes apart from any explicitly reported non-fatal soft warning
13. E.8.4.Fix.3 is `SOURCE REVIEWED`, not runtime validated
14. profile re-pin remains `ACTIVE` pending independent review
15. no farm/runtime/production claim is introduced
16. no commit or push occurs
17. complete review bundle is produced

## Hard stop

Stop after:

- exact profile/test re-pin
- warranted status/history reconciliation
- manifest regeneration
- deterministic local checks
- diff audit
- creation of `kaonlt_review_e8_4_fix3_profile_repin.diff`

Do not commit.
Do not push.
Do not run the farm.
Do not provide later-gate commands.

NEXT is only independent ChatGPT review of the review bundle.
