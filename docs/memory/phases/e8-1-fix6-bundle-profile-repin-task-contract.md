# E.8.1.Fix.6 — re-pin the E.8.1 validation profile to the pushed Fix.5 overlay renderer source

## Objective

Re-pin the existing E.8.1 canonical-five-setting validation bundle profile to
the now-pushed, independently source-reviewed E.8.1.Fix.5 persisted-overlay
ROOT-PDF geometry source.

This is a narrow validation/provenance packaging change only. It must not alter
analysis physics, the Fix.5 renderer, either full-background-subtraction test
module, the generic collector, bundle wrapper, artifact declarations, canonical
settings, accepted F.6.2 evidence, or any production behavior.

## Exact starting state

Start only from committed `test` HEAD:

`53fd262b730af8f1254e411a38231aebeb6a1da3`

Commit message:

`Repair E8.1 persisted overlay PDF geometry`

Its parent is:

`fdd368f2084c9b9508e7ce679b8f7391f5b556f1`

The pushed Fix.5 commit is the exact analysis/procedure renderer source that
the E.8.1 validation profile must now require:

`53fd262b730af8f1254e411a38231aebeb6a1da3`

Before editing:

1. establish branch/HEAD/worktree;
2. read the normal five-file repository-memory startup sequence in order;
3. read the current Fix.5 phase record and Fix.4 farm blocker evidence;
4. inspect the existing E.8.1 profile and focused profile test;
5. confirm this task-contract file is the only expected user-created
   uncommitted file;
6. stop if the committed HEAD differs or unrelated worktree changes exist.

## Pushed-state source review already established

Independent ChatGPT pushed-state review established:

- remote `test` HEAD is exactly
  `53fd262b730af8f1254e411a38231aebeb6a1da3`;
- its parent is exactly
  `fdd368f2084c9b9508e7ce679b8f7391f5b556f1`;
- the commit contains exactly the reviewed Fix.5 12-path change set;
- the renderer still uses the reviewed square `3600 x 3600` overlay canvas;
- the canonical `3 x 9` persisted-overlay grid remains intact;
- the grid is rendered before the separate top-level header pad;
- L/B/A legend semantics/colors remain unchanged;
- `Modified()` / `Update()` occur before PDF `Print()`;
- the focused renderer regression still preserves canonical zero-based
  `phi_index` values and derives the first title as
  `Left-lowe phi0 [-180, -140)`;
- no E.8.1 profile, collector, wrapper, `main.py`, `rand_sub.py`, accepted
  F.6.2 artifact, or unrelated production/scientific source changed.

Fix.5 remains source-reviewed only. The real-PyROOT PDF regression was not run
locally because PyROOT was unavailable. No Fix.5 farm/runtime acceptance is
claimed yet.

## Current profile state

The current E.8.1 profile still requires the older Fix.3 renderer source:

`350c34c55b2de33ad01011559dc6d8ed84d9c8a7`

The current focused profile test uses the matching constant:

`REVIEWED_SOURCE = "350c34c55b2de33ad01011559dc6d8ed84d9c8a7"`

That old pin correctly described the Fix.4 farm gate. It must now move to the
new pushed Fix.5 renderer source.

Do not admit Fix.5 by broadening `allowed_committed_files`. Replace the exact
required analysis source identity.

## Status boundary

Preserve/update durable status as follows:

- F.6.2 — `CLOSED / RUNTIME VALIDATED`
- F.6.2.Fix.5 — `CLOSED / RUNTIME VALIDATED`
- E.8.1 — `ACTIVE`
- E.8.1.Fix.1 — `CLOSED / RUNTIME VALIDATED`
- E.8.1.Fix.2 — `CLOSED / RUNTIME VALIDATED`
- E.8.1.Fix.3 — `SOURCE REVIEWED`; its context/handoff repair has farm evidence,
  while the old overlay subrepair was superseded by Fix.5
- E.8.1.Fix.4 — `CLOSED / RUNTIME VALIDATED` only for its narrow profile
  provenance re-pin to Fix.3 source
- E.8.1.Fix.5 — `SOURCE REVIEWED`, user-pushed at
  `53fd262b730af8f1254e411a38231aebeb6a1da3`; no ROOT/PyROOT/farm/runtime
  acceptance yet
- E.8.1.Fix.6 — `ACTIVE` until independent ChatGPT actual-diff review
- E.8.1.Debug.1.Fix.1 — `SOURCE REVIEWED`
- F.6.3 — `BLOCKED` pending E.8.1
- F.6.4 — `BLOCKED` pending F.6.3 evidence
- lifecycle-hook dispatch — `BLOCKED / DEFERRED`

Do not close E.8.1 or Fix.5 from this profile-only task.

## Scientific and architectural ownership

E.8.1 remains presentation-only.

Method A remains detached/non-production.
Method B remains diagnostic/cross-check only.

This task changes validation source provenance only.

It must not change:

- accepted F.6.2 JSON SHA-256, artifact fingerprint, validation fingerprint,
  or numerical content;
- cuts, normalization, binning, templates, priors, fits, subtraction formulas,
  production weights, SIMC, yields, efficiencies, cross sections, L/T
  separation, or uncertainties;
- E.8 reader/authority behavior;
- Fix.5 overlay renderer geometry/draw order;
- context, handoff, or map-page renderers;
- page IDs or page inventory;
- canonical-five-setting membership.

## Allowed implementation files

Only:

- `testing/pion_hgcer_validation_bundle_profile_e8_1.json`
- `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`

## Allowed durable-memory files

Update only as warranted:

- `docs/memory/CURRENT.md`
- `docs/memory/phases/e8-1-fix5-overlay-pdf-geometry.md`
- `docs/memory/phases/e8-1-fix6-bundle-profile-repin.md`
- `docs/memory/manifest.json`

This task contract is also an intended tracked file:

- `docs/memory/phases/e8-1-fix6-bundle-profile-repin-task-contract.md`

Regenerate `docs/memory/manifest.json` after all versioned-memory changes.

Do not modify `MEMORY.md`, `USER.md`, `CURRENT_HANDOFF.md`, roadmap, decisions,
or evidence records unless a concrete blocker proves their owned durable
knowledge changed. Stop rather than expand scope automatically.

## Frozen files and interfaces

Do not modify:

- `src/cuts/full_background_subtraction_plots.py`
- `testing/test_full_background_subtraction_plots.py`
- `run_Prod_Analysis.sh`
- `src/main.py`
- `src/cuts/rand_sub.py`
- `testing/collect_pion_hgcer_validation_bundle.py`
- `testing/package_pion_hgcer_validation_bundle.tcsh`
- any other validation profile;
- any F.6.2 analyzer/renderer/profile;
- accepted farm artifacts or evidence;
- `docs/memory/evidence/e8-1-fix4-left-lowe-overlay-blocker.md`
- Fix.5 implementation/task-review contracts except the allowed phase status
  record named above;
- E.8.1 artifact declarations;
- E.8.1 canonical setting declarations.

The generic collector and wrapper already provide the required fail-closed
source-provenance behavior. Do not edit them.

## Required profile change

In:

`testing/pion_hgcer_validation_bundle_profile_e8_1.json`

change only:

```json
"required_analysis_commit": "350c34c55b2de33ad01011559dc6d8ed84d9c8a7"
```

to:

```json
"required_analysis_commit": "53fd262b730af8f1254e411a38231aebeb6a1da3"
```

Preserve exactly:

- `schema_version`:
  `pion_hgcer_validation_bundle_profile/v4`
- `validation_profile`:
  `phase_e8_1_full_background_procedure_pdf_farm_review/v1`
- `collection_mode`:
  `generic_artifacts`
- settings, in existing order:
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

Do not add `src/cuts/full_background_subtraction_plots.py`,
`testing/test_full_background_subtraction_plots.py`, or any other source path
to `allowed_committed_files`.

Do not bump the profile schema or validation-profile identity.

## Required focused-test change

In:

`testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`

change only:

```python
REVIEWED_SOURCE = "350c34c55b2de33ad01011559dc6d8ed84d9c8a7"
```

to:

```python
REVIEWED_SOURCE = "53fd262b730af8f1254e411a38231aebeb6a1da3"
```

Preserve the existing source-identity and fail-closed tests.

The test suite must continue to establish that:

- the ordered canonical five settings are exact;
- artifact declarations are unchanged;
- `required_analysis_commit == REVIEWED_SOURCE`;
- only the profile/test pair is allowed as committed changes after the required
  analysis source;
- `docs/memory/` remains the sole allowed non-analysis prefix;
- a complete synthetic package can be collected;
- missing/invalid required artifacts fail completeness;
- missing source ancestry fails completeness;
- any later analysis-source change, including
  `src/cuts/full_background_subtraction_plots.py`, is rejected;
- output collisions still fail without overwrite.

Do not weaken negative provenance coverage.

## Before/after provenance behavior

Before:

`required analysis source = 350c34c55b2de33ad01011559dc6d8ed84d9c8a7`

which correctly supported the Fix.4 farm gate but treats the later Fix.5
renderer source as an unexpected analysis change.

After:

`required analysis source = 53fd262b730af8f1254e411a38231aebeb6a1da3`

and the future Fix.6 profile/bundle commit may differ from that source only by
the two explicitly allowed profile/test files plus `docs/memory/`.

This preserves fail-closed provenance.

## Required durable-memory update

### CURRENT.md

Update stale pre-push wording and record:

- Fix.5 is `SOURCE REVIEWED` and user-pushed at
  `53fd262b730af8f1254e411a38231aebeb6a1da3`;
- no Fix.5 ROOT/PyROOT/farm/runtime validation is claimed;
- Fix.6 profile re-pin is `ACTIVE` until independent ChatGPT actual-diff
  review;
- Fix.4 remains `CLOSED / RUNTIME VALIDATED` only for its narrow prior
  profile-provenance gate;
- E.8.1 remains `ACTIVE`;
- F.6.3 remains `BLOCKED`;
- lifecycle/debug continuity facts remain unchanged;
- no Fix.6 commit exists yet.

Set NEXT to:

1. independent ChatGPT actual-diff review of Fix.6;
2. user-controlled commit/push of the reviewed Fix.6 profile re-pin;
3. ChatGPT pushed-state review;
4. use that future pushed Fix.6 commit as wrapper `--bundle-commit`;
5. keep profile `required_analysis_commit` equal to distinct Fix.5 source
   `53fd262b730af8f1254e411a38231aebeb6a1da3`;
6. run fresh targeted `Q4p4W2p74 / Left / lowe`;
7. collect fresh bundle/PDF;
8. inspect pages 37--47, especially persisted-overlay pages 38/41/44;
9. only after Left/lowe visual PASS return to detailed `Left / highe`;
10. broaden only after that detailed gate passes.

### Fix.5 phase record

Update `docs/memory/phases/e8-1-fix5-overlay-pdf-geometry.md` only to append the
actual pushed source identity:

`53fd262b730af8f1254e411a38231aebeb6a1da3`

Record that pushed-state source review passed and remains distinct from
farm/runtime validation. State that the separate Fix.6 provenance re-pin is the
next required step.

Do not change Fix.5 to runtime validated.

### Fix.6 phase record

Create:

`docs/memory/phases/e8-1-fix6-bundle-profile-repin.md`

Record:

- starting HEAD `53fd262b730af8f1254e411a38231aebeb6a1da3`;
- status `ACTIVE` pending independent actual-diff review;
- exact required-source re-pin from Fix.3 `350c34c...` to Fix.5 `53fd262b...`;
- unchanged allowlists/settings/artifacts;
- no scientific or renderer change;
- local deterministic validation results when available;
- farm boundary;
- future profile/bundle commit remains unknown and distinct from the required
  Fix.5 source.

Do not create new farm evidence in this task.

## Deterministic local validation

Run at minimum:

```text
python -B -m py_compile \
  testing/test_pion_hgcer_validation_bundle_profile_e8_1.py

python -m json.tool \
  testing/pion_hgcer_validation_bundle_profile_e8_1.json

python -B -m unittest \
  testing.test_pion_hgcer_validation_bundle_profile_e8_1 \
  -v
```

Also run the generic collector regression:

```text
python -B -m unittest \
  testing.test_collect_pion_hgcer_validation_bundle \
  -v
```

Report exact tests, failures, and skips.

Run normal memory manifest/integrity/health/bootstrap checks and:

```text
git diff --check
```

These are source/provenance checks only, not ROOT or farm validation.

## Required actual-diff audit

Before stopping:

1. confirm committed HEAD remains
   `53fd262b730af8f1254e411a38231aebeb6a1da3`;
2. show `git status --short`;
3. show `git diff --stat`;
4. show `git diff --name-only`;
5. confirm every changed path is allowlisted by this task;
6. verify no `src/` path changed;
7. verify `testing/test_full_background_subtraction_plots.py` is untouched;
8. verify profile settings/artifacts/allowlists are unchanged except the exact
   required source commit;
9. verify the focused profile test changes only `REVIEWED_SOURCE` unless a
   narrowly necessary assertion update is required;
10. trace the intended packaging path:
    future Fix.6 profile commit -> wrapper `--bundle-commit` -> detached
    worktree -> E.8.1 profile -> generic collector -> exact ancestry/range
    checks -> existing required artifacts.

Create one complete temporary root-level:

`kaonlt_review.diff`

containing the full tracked diff plus complete
`git diff --no-index /dev/null ...` additions for every intended new/untracked
file.

Do not stage merely to create the review bundle.

## Farm boundary and future identities

Do not run the farm in this task.

Do not hard-code or guess the future Fix.6 profile commit.

After user commit/push of Fix.6:

- Fix.5 analysis/procedure source remains:
  `53fd262b730af8f1254e411a38231aebeb6a1da3`
- the future pushed Fix.6 commit becomes the wrapper `--bundle-commit`.

These identities must remain distinct.

## Forbidden shortcuts

Do not:

- broaden `allowed_committed_files`;
- allow any `src/` path after the required analysis commit;
- change artifact declarations;
- change canonical settings;
- change collector or wrapper source;
- change the Fix.5 renderer or its regression tests;
- change accepted F.6.2 artifacts;
- bump profile/version merely for this provenance update;
- insert the future Fix.6 commit as `required_analysis_commit`;
- commit, push, or run the farm.

## Acceptance criteria

Source-review candidate PASS requires:

- profile `required_analysis_commit` is exactly
  `53fd262b730af8f1254e411a38231aebeb6a1da3`;
- `REVIEWED_SOURCE` matches exactly;
- profile/test remain the only allowed committed files;
- canonical settings and artifact declarations are unchanged;
- source identity remains fail-closed;
- later analysis-source changes remain rejected;
- collector/wrapper and all `src/` files are untouched;
- deterministic local tests pass to local capability;
- memory distinguishes Fix.5 analysis source from future Fix.6 profile/bundle
  source;
- actual diff is narrow and allowlisted.

## Hard stop

Do not commit, push, or run Jefferson Lab farm validation.

Return to ChatGPT with:

- implementation summary;
- exact changed files;
- exact local tests/checks;
- provenance/runtime-path trace;
- complete root-level `kaonlt_review.diff`;
- proposed user Git handoff commands, not executed.

ChatGPT must inspect the actual diff before the user commits or pushes.
