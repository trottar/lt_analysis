# E.8.1.Fix.6 — post-review source-status reconciliation

## Purpose

Reconcile durable repository memory after independent ChatGPT actual-diff review
of the complete E.8.1.Fix.6 profile-provenance re-pin.

The Fix.6 profile/test candidate has now passed independent actual-diff review.
This task is memory-only. It must not alter the reviewed profile, focused
profile test, Fix.5 renderer/test source, collector, wrapper, artifacts,
analysis/runtime source, or accepted F.6.2 evidence.

## Exact committed base

Committed `test` HEAD remains:

`53fd262b730af8f1254e411a38231aebeb6a1da3`

The existing uncommitted E.8.1.Fix.6 worktree is intentional.

Before editing:

1. establish exact branch/HEAD/worktree;
2. read the normal five-file repository-memory startup sequence;
3. read the Fix.5 and Fix.6 phase records and current Fix.6 task contract;
4. confirm there are no unrelated worktree changes.

## Independent ChatGPT actual-diff result

ChatGPT inspected the complete refreshed Fix.6 review bundle and found the
profile/test provenance re-pin PASS.

The reviewed implementation:

- changes
  `testing/pion_hgcer_validation_bundle_profile_e8_1.json`
  only by moving
  `source_identity.required_analysis_commit`
  from
  `350c34c55b2de33ad01011559dc6d8ed84d9c8a7`
  to the pushed Fix.5 renderer source
  `53fd262b730af8f1254e411a38231aebeb6a1da3`;
- changes
  `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`
  only by moving `REVIEWED_SOURCE` to the same exact Fix.5 source;
- leaves the profile schema, validation-profile identity, collection mode,
  ordered canonical five settings, global artifact declaration, per-setting
  artifact declarations, `allowed_committed_files`, and
  `allowed_non_analysis_path_prefixes` unchanged;
- does not add any `src/` exception;
- preserves fail-closed rejection of later analysis-source changes;
- leaves the Fix.5 renderer and its regression tests unchanged;
- leaves the generic collector, bundle wrapper, accepted F.6.2 artifact, and
  all production/scientific source unchanged;
- correctly distinguishes the required Fix.5 analysis/procedure source
  `53fd262b730af8f1254e411a38231aebeb6a1da3`
  from the future Fix.6 profile/bundle commit, which is not known yet.

The actual changed-path set is exactly the seven intended Fix.6 paths:

- `docs/memory/CURRENT.md`
- `docs/memory/manifest.json`
- `docs/memory/phases/e8-1-fix5-overlay-pdf-geometry.md`
- `testing/pion_hgcer_validation_bundle_profile_e8_1.json`
- `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`
- `docs/memory/phases/e8-1-fix6-bundle-profile-repin-task-contract.md`
- `docs/memory/phases/e8-1-fix6-bundle-profile-repin.md`

No `src/` path changed.

## Test/check boundary

Codex-reported deterministic checks remain `NOT RUN by ChatGPT`:

- focused profile-test `py_compile`: PASS;
- profile JSON parse/format check: PASS;
- focused E.8.1 bundle-profile tests: PASS, 5 tests, no skips;
- generic collector tests: PASS, 28 tests, no skips;
- repository memory manifest/integrity/health/bootstrap checks: PASS;
- `git diff --check`: PASS.

These are source/provenance checks only. No ROOT/PyROOT/farm/runtime validation
is claimed for Fix.6.

## Status boundary

Preserve/update:

- F.6.2 — `CLOSED / RUNTIME VALIDATED`
- F.6.2.Fix.5 — `CLOSED / RUNTIME VALIDATED`
- E.8.1 — `ACTIVE`
- E.8.1.Fix.1 — `CLOSED / RUNTIME VALIDATED`
- E.8.1.Fix.2 — `CLOSED / RUNTIME VALIDATED`
- E.8.1.Fix.3 — `SOURCE REVIEWED`
- E.8.1.Fix.4 — `CLOSED / RUNTIME VALIDATED` only for its narrow prior
  profile-provenance gate
- E.8.1.Fix.5 — `SOURCE REVIEWED`, user-pushed at
  `53fd262b730af8f1254e411a38231aebeb6a1da3`, with no farm/runtime acceptance
- E.8.1.Fix.6 — now `SOURCE REVIEWED`
- E.8.1.Debug.1.Fix.1 — `SOURCE REVIEWED`
- F.6.3 — `BLOCKED` pending E.8.1
- F.6.4 — `BLOCKED` pending F.6.3 evidence
- lifecycle-hook dispatch — `BLOCKED / DEFERRED`

Do not upgrade E.8.1 or Fix.5/Fix.6 to runtime validated.

## Allowed files

Modify only:

- `docs/memory/CURRENT.md`
- `docs/memory/phases/e8-1-fix6-bundle-profile-repin.md`
- `docs/memory/manifest.json`

This reconciliation contract is also an intended tracked file:

- `docs/memory/phases/e8-1-fix6-source-review-reconciliation-task-contract.md`

Do not modify any other file.

In particular, the following must remain byte-for-byte unchanged during this
reconciliation:

- `testing/pion_hgcer_validation_bundle_profile_e8_1.json`
- `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`
- `docs/memory/phases/e8-1-fix5-overlay-pdf-geometry.md`
- `docs/memory/phases/e8-1-fix6-bundle-profile-repin-task-contract.md`
- `src/cuts/full_background_subtraction_plots.py`
- `testing/test_full_background_subtraction_plots.py`
- all collector/wrapper files;
- `run_Prod_Analysis.sh`;
- `src/main.py`;
- `src/cuts/rand_sub.py`;
- accepted F.6.2 artifacts and evidence.

## Required CURRENT.md reconciliation

Record one unambiguous current state:

- E.8.1 remains `ACTIVE`;
- Fix.5 remains `SOURCE REVIEWED`, pushed at
  `53fd262b730af8f1254e411a38231aebeb6a1da3`;
- Fix.6 is now `SOURCE REVIEWED` from independent ChatGPT actual-diff review;
- Codex-reported Fix.6 checks were `NOT RUN by ChatGPT`;
- no Fix.6 ROOT/PyROOT/farm/runtime validation is claimed;
- no Fix.6 commit exists yet because the user has not committed/pushed it;
- the profile now requires the distinct Fix.5 source
  `53fd262b730af8f1254e411a38231aebeb6a1da3`;
- the future Fix.6 profile/bundle commit remains unknown;
- F.6.3 remains `BLOCKED`;
- lifecycle/debug continuity facts remain unchanged.

Set exact NEXT to:

1. user-controlled commit/push of the reviewed Fix.6 profile/test/memory set;
2. ChatGPT pushed-state review;
3. use the resulting pushed Fix.6 commit as wrapper `--bundle-commit`;
4. keep profile `required_analysis_commit` equal to the distinct Fix.5 source
   `53fd262b730af8f1254e411a38231aebeb6a1da3`;
5. run fresh targeted `Q4p4W2p74 / Left / lowe`;
6. collect a fresh validation bundle and procedure PDF;
7. independently inspect pages 37--47;
8. require pages 38/41/44 to show the E.8 header, L/B/A legend, first canonical
   child title, and unclipped first row;
9. only after Left/lowe visual PASS return to detailed `Left / highe`;
10. broaden only after that detailed gate passes.

Do not create a new farm-evidence record yet.

## Required Fix.6 phase-record reconciliation

Update:

`docs/memory/phases/e8-1-fix6-bundle-profile-repin.md`

so that:

- status becomes `SOURCE REVIEWED`;
- it records that ChatGPT independently inspected the complete actual diff and
  found the profile/test provenance re-pin PASS;
- it states the exact required Fix.5 analysis/procedure source:
  `53fd262b730af8f1254e411a38231aebeb6a1da3`;
- it preserves the unchanged settings/artifacts/allowlists and fail-closed
  source boundary;
- Codex-reported local checks are explicitly `NOT RUN by ChatGPT`;
- no ROOT/PyROOT/farm/runtime validation is claimed;
- user commit/push remains pending;
- the future Fix.6 profile/bundle commit remains distinct from the required
  Fix.5 source;
- NEXT is user commit/push -> pushed-state review -> fresh Left/lowe farm/PDF
  validation.

## Validation

Regenerate:

`docs/memory/manifest.json`

Run applicable repository-memory manifest/integrity/health/bootstrap checks and:

```text
git diff --check
```

Verify explicitly that the already reviewed profile/test/Fix.5 phase files have
no new diff introduced during this reconciliation.

## Review bundle

Refresh one complete temporary root-level:

`kaonlt_review.diff`

It must contain:

1. the complete tracked diff for the entire current Fix.6 worktree;
2. complete `git diff --no-index /dev/null ...` representations for every
   intended new/untracked Fix.6 file, including this reconciliation contract.

Do not stage merely to create the review bundle.

## Hard stop

Do not commit, push, alter the profile/test implementation, or run Jefferson
Lab farm validation.

Return:

- exact changed paths;
- exact memory/check results;
- confirmation that profile/test/renderer/collector/wrapper files were
  untouched during this reconciliation;
- refreshed root-level `kaonlt_review.diff` ready for final ChatGPT review.
