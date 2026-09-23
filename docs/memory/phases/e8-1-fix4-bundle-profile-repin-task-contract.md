# E.8.1.Fix.4 — re-pin the E.8.1 validation profile to the Fix.3 layout source

## Objective

Re-pin the existing E.8.1 canonical-five-setting validation bundle profile to
the now-pushed, independently source-reviewed E.8.1.Fix.3 procedure-PDF layout
source.

This is a narrow validation/provenance packaging change only. It must not alter
analysis physics, the Fix.3 renderer, the generic collector, the bundle wrapper,
artifact declarations, canonical settings, or accepted F.6.2 evidence.

## Exact starting state

Start only from committed `test` HEAD:

`350c34c55b2de33ad01011559dc6d8ed84d9c8a7`

Commit message:

`Fix E8.1 procedure PDF layout`

Its parent is:

`9f7094b19f3cdb0887c95f1f972d98471f2ecfdc`

Before editing:

1. establish branch/HEAD/worktree;
2. read the normal five-file repository-memory startup sequence;
3. read the active E.8.1.Fix.3 phase record and its layout-blocker evidence;
4. confirm this task-contract file is the only expected user-created
   uncommitted file;
5. stop if the starting HEAD differs or unrelated worktree changes exist.

The pushed E.8.1.Fix.3 source commit is now the exact analysis/procedure source
that the E.8.1 validation profile must require:

`350c34c55b2de33ad01011559dc6d8ed84d9c8a7`

The current profile still pins the earlier Fix.1 reader source:

`8985d9a212799c021c4ad1a759a689ea3826e0ea`

The old pin must not be broadened by allowlisting the later renderer source.
Replace the required source identity with the exact pushed Fix.3 commit.

## Current status boundary

Preserve:

- F.6.2 — `CLOSED / RUNTIME VALIDATED`
- F.6.2.Fix.5 — `CLOSED / RUNTIME VALIDATED`
- E.8.1 — `ACTIVE`
- E.8.1.Fix.1 — `CLOSED / RUNTIME VALIDATED` only for the narrow reader defect
- E.8.1.Fix.2 — `CLOSED / RUNTIME VALIDATED` only for the narrow profile
  provenance defect that admitted the Fix.1 source
- E.8.1.Fix.3 — `SOURCE REVIEWED`, user-pushed at
  `350c34c55b2de33ad01011559dc6d8ed84d9c8a7`; no ROOT/PyROOT/farm/runtime
  validation yet
- E.8.1.Debug.1.Fix.1 — `SOURCE REVIEWED`
- F.6.3 — `BLOCKED` pending E.8.1

This task is the separate provenance re-pin required before a fresh Fix.3
Left/lowe farm/PDF gate.

## Scientific and architectural ownership

E.8.1 remains presentation-only.

Method A remains detached/non-production.
Method B remains diagnostic/cross-check only.

This task changes validation source provenance only. It must not change:

- accepted F.6.2 JSON, artifact fingerprint, validation fingerprint, or values;
- cuts, binning, normalizations, templates, priors, fits, subtraction formulas,
  production weights, SIMC, yields, efficiencies, cross sections, L/T
  separation, or uncertainties;
- Fix.3 text/header/grid rendering;
- ordinary D.6-D.9 + E.8 page semantics or page inventory;
- canonical-five-setting membership.

## Allowed implementation files

Only:

- `testing/pion_hgcer_validation_bundle_profile_e8_1.json`
- `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`

## Allowed durable-memory files

Update only as warranted:

- `docs/memory/CURRENT.md`
- `docs/memory/phases/e8-1-fix3-procedure-pdf-layout.md`
- `docs/memory/phases/e8-1-fix4-bundle-profile-repin.md`
- `docs/memory/manifest.json`

This task contract is also an intended tracked file:

- `docs/memory/phases/e8-1-fix4-bundle-profile-repin-task-contract.md`

Regenerate `docs/memory/manifest.json` after all versioned memory changes.

Do not change `MEMORY.md`, `USER.md`, `CODEX.md`, `TOOLS.md`,
`CURRENT_HANDOFF.md`, roadmap, decisions, or evidence records unless a concrete
blocker proves their owned durable knowledge changed. Stop rather than expand
scope automatically.

## Frozen files and interfaces

Do not modify:

- anything under `src/`;
- `run_Prod_Analysis.sh`;
- `src/main.py`;
- `src/cuts/rand_sub.py`;
- `src/cuts/full_background_subtraction_plots.py`;
- `testing/test_full_background_subtraction_plots.py`;
- `testing/collect_pion_hgcer_validation_bundle.py`;
- `testing/package_pion_hgcer_validation_bundle.tcsh`;
- any other validation profile;
- any F.6.2 analyzer/renderer/profile;
- any accepted farm artifact;
- `docs/memory/evidence/e8-1-fix2-left-lowe-layout-blocker.md`;
- E.8.1 artifact declarations;
- E.8.1 canonical setting declarations.

The generic collector and wrapper already implement the required fail-closed
provenance semantics. Do not edit them.

## Required profile change

In:

`testing/pion_hgcer_validation_bundle_profile_e8_1.json`

change only:

```json
"required_analysis_commit": "8985d9a212799c021c4ad1a759a689ea3826e0ea"
```

to:

```json
"required_analysis_commit": "350c34c55b2de33ad01011559dc6d8ed84d9c8a7"
```

Preserve exactly:

- `schema_version`:
  `pion_hgcer_validation_bundle_profile/v4`
- `validation_profile`:
  `phase_e8_1_full_background_procedure_pdf_farm_review/v1`
- `collection_mode`:
  `generic_artifacts`
- canonical five settings, in current order:
  - Left / lowe
  - Left / highe
  - Center / lowe
  - Center / highe
  - Right / highe
- one required global frozen F.6.2 JSON artifact;
- two required per-setting artifacts:
  - full-background-subtraction PDF
  - full-background-subtraction page-manifest JSON
- `allowed_committed_files`, exactly:
  - `testing/pion_hgcer_validation_bundle_profile_e8_1.json`
  - `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`
- `allowed_non_analysis_path_prefixes`, exactly:
  - `docs/memory/`

Do not bump the profile schema or `validation_profile` identifier. The evidence
package contract is unchanged; only its exact required analysis/procedure
source identity advances to Fix.3.

Do not add `src/cuts/full_background_subtraction_plots.py`, its test, or any
other source file to `allowed_committed_files`. The correct provenance model is
a new exact required source commit, not a broader exception list.

## Required focused-test change

In:

`testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`

change only the exact source identity constant:

```python
REVIEWED_SOURCE = "8985d9a212799c021c4ad1a759a689ea3826e0ea"
```

to:

```python
REVIEWED_SOURCE = "350c34c55b2de33ad01011559dc6d8ed84d9c8a7"
```

Preserve the existing source-identity tests and their fail-closed behavior.

The tests must continue to prove that:

- the profile declares exactly the canonical five settings;
- artifact declarations are unchanged;
- the profile `required_analysis_commit` equals `REVIEWED_SOURCE`;
- only the profile/test paths are allowed committed files after the required
  source;
- `docs/memory/` remains the only allowed non-analysis prefix;
- a complete synthetic package can be collected;
- missing/invalid required artifacts remain incomplete;
- missing required source ancestry remains incomplete;
- any analysis-source change after the required source, including
  `src/cuts/full_background_subtraction_plots.py`, remains rejected;
- output collisions still fail without overwrite.

Do not weaken or remove negative provenance coverage merely because the Fix.3
renderer now becomes the required source baseline.

## Before/after provenance behavior

Before:

`required analysis source = 8985d9a...`

which correctly described the prior Fix.2 gate but treats the later Fix.3
renderer commit as an unexpected analysis-source change.

After:

`required analysis source = 350c34c...`

and the future profile-repin commit may differ from that source only by the two
explicitly allowed profile/test files plus `docs/memory/`.

This preserves fail-closed provenance while admitting exactly the reviewed
Fix.3 layout source.

## Deterministic local validation

Run at minimum:

```text
python -B -m py_compile testing/test_pion_hgcer_validation_bundle_profile_e8_1.py
python -m json.tool testing/pion_hgcer_validation_bundle_profile_e8_1.json
python -B -m unittest testing.test_pion_hgcer_validation_bundle_profile_e8_1 -v
```

Also run the generic collector regression module if locally deterministic:

```text
python -B -m unittest testing.test_collect_pion_hgcer_validation_bundle -v
```

Report exact skips/failures.

Run the normal memory manifest/integrity/health/bootstrap checks required by
current repository instructions and:

```text
git diff --check
```

These are source/provenance checks only, not farm or ROOT validation.

## Required diff/provenance audit

Before stopping:

1. confirm committed HEAD remains
   `350c34c55b2de33ad01011559dc6d8ed84d9c8a7`;
2. show `git status --short`;
3. show `git diff --stat`;
4. show `git diff --name-only`;
5. confirm every changed path is allowlisted;
6. verify no `src/` file changed;
7. verify profile settings/artifacts/allowlists are unchanged except the exact
   required commit field;
8. verify the focused test changes only the expected `REVIEWED_SOURCE` identity
   unless a narrowly necessary assertion update is required;
9. trace the intended farm packaging path:
   future Fix.4 profile commit -> wrapper `--bundle-commit` -> detached worktree
   -> E.8.1 profile -> generic collector -> exact ancestry/range checks ->
   existing required artifacts only.

Create one complete temporary root-level:

`kaonlt_review.diff`

containing the tracked diff plus complete
`git diff --no-index /dev/null ...` additions for every intended new/untracked
file. Do not stage merely to create the review.

## Durable-memory update

Update `CURRENT.md` to record:

- live/pushed E.8.1.Fix.3 source commit:
  `350c34c55b2de33ad01011559dc6d8ed84d9c8a7`;
- E.8.1 remains `ACTIVE`;
- E.8.1.Fix.3 remains `SOURCE REVIEWED`, now user-pushed;
- E.8.1.Fix.4 profile re-pin is `ACTIVE` until ChatGPT independently reviews
  its actual diff;
- Fix.1 and Fix.2 retain their narrow `CLOSED / RUNTIME VALIDATED` statuses;
- no Fix.3 farm/runtime acceptance is claimed;
- F.6.3 remains `BLOCKED`;
- NEXT after source review is user commit/push of Fix.4, then pushed-state
  review, then a fresh targeted `Q4p4W2p74 / Left / lowe` analysis and bundle;
- that fresh PDF must be independently inspected for the exact Fix.3 clipping,
  mojibake, and overlay-header defects before returning to `Left / highe`.

Update the existing Fix.3 phase record only to append:

- actual pushed source commit `350c34c...`;
- pushed-state source review remains distinct from farm validation;
- separate Fix.4 profile re-pin is the next provenance step.

Create:

`docs/memory/phases/e8-1-fix4-bundle-profile-repin.md`

as the permanent chronology for this profile-only re-pin.

Do not create Fix.3 runtime-acceptance evidence yet.

## Farm boundary and future identities

Do not run the farm in this task.

Do not hard-code or guess the future Fix.4 profile commit into any file.

After user commit/push of this task:

- ChatGPT will inspect the pushed Fix.4 commit;
- that future pushed Fix.4 commit becomes the wrapper `--bundle-commit`;
- the profile `required_analysis_commit` remains the already known Fix.3
  analysis/procedure source:
  `350c34c55b2de33ad01011559dc6d8ed84d9c8a7`.

These two identities must remain distinct.

## Forbidden shortcuts

Do not:

- broaden `allowed_committed_files`;
- allow any `src/` path after the required analysis commit;
- change artifact declarations;
- change canonical settings;
- change the collector or wrapper;
- change Fix.3 renderer/test source;
- modify accepted F.6.2 artifacts;
- bump profile/version merely for the provenance change;
- insert the future profile commit as `required_analysis_commit`;
- commit, push, or run the farm.

## Acceptance criteria

Source-review candidate PASS requires:

- profile required source is exactly
  `350c34c55b2de33ad01011559dc6d8ed84d9c8a7`;
- `REVIEWED_SOURCE` matches exactly;
- profile/test are the only allowed committed files;
- canonical settings and artifact declarations are unchanged;
- source identity remains fail-closed;
- later analysis-source changes remain rejected;
- collector and wrapper are untouched;
- no `src/` or production file changed;
- deterministic local tests pass to local capability;
- memory distinguishes Fix.3 analysis source from future Fix.4 profile/bundle
  source;
- actual diff is narrow and allowlisted.

## Hard stop

Do not commit, push, or run Jefferson Lab farm validation.

Return to ChatGPT with:

- implementation summary;
- exact changed files;
- exact local test/check results;
- provenance/runtime-path trace;
- complete root-level `kaonlt_review.diff`;
- proposed user Git handoff commands, not executed.

ChatGPT must inspect the actual diff before the user commits or pushes.
