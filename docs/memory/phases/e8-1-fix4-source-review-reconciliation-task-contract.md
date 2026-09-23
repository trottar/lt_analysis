# E.8.1.Fix.4 — post-review status reconciliation

## Purpose

Reconcile durable repository memory after independent ChatGPT actual-diff review
of E.8.1.Fix.4.

The profile/test provenance re-pin has now been independently reviewed and is a
PASS. This task is memory-only. It must not alter the E.8.1 profile, its focused
test, Fix.3 renderer/test source, any analysis/runtime source, the collector,
wrapper, or accepted F.6.2 artifacts.

## Starting committed state

Committed `test` HEAD remains:

`350c34c55b2de33ad01011559dc6d8ed84d9c8a7`

The existing uncommitted E.8.1.Fix.4 profile/test/memory changes are intentional.

Before editing:

1. establish branch/HEAD/worktree;
2. read the normal five-file repository-memory startup sequence;
3. read the current Fix.4 task contract and phase record;
4. confirm there are no unrelated worktree changes.

## Independent review result

ChatGPT inspected the complete actual `kaonlt_review.diff` and found the
E.8.1.Fix.4 source candidate PASS.

The reviewed change set:

- changes only
  `testing/pion_hgcer_validation_bundle_profile_e8_1.json`
  `source_identity.required_analysis_commit` from
  `8985d9a212799c021c4ad1a759a689ea3826e0ea`
  to
  `350c34c55b2de33ad01011559dc6d8ed84d9c8a7`;
- changes only the matching `REVIEWED_SOURCE` constant in
  `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`;
- leaves the profile schema, validation-profile identifier, collection mode,
  ordered canonical five settings, artifact declarations,
  `allowed_committed_files`, and `allowed_non_analysis_path_prefixes`
  unchanged;
- adds no `src/` exception and preserves fail-closed rejection of later
  analysis-source changes;
- leaves Fix.3 renderer/test source, the generic collector, bundle wrapper,
  accepted F.6.2 artifacts, and production/scientific source untouched;
- correctly distinguishes the pushed Fix.3 analysis/procedure source
  `350c34c...` from the future Fix.4 profile/bundle commit.

The actual changed paths are all within the Fix.4 task allowlist.

Codex-reported deterministic checks remain `NOT RUN by ChatGPT`:

- profile-test `py_compile`: PASS;
- JSON parse/format check: PASS;
- focused E.8.1 bundle-profile tests: PASS, 5 tests, no skips;
- generic collector tests: PASS, 28 tests, no skips;
- memory manifest/health/bootstrap and `git diff --check`: PASS.

No ROOT/PyROOT/farm/runtime validation is claimed for Fix.4.

## Allowed files

Modify only:

- `docs/memory/CURRENT.md`
- `docs/memory/phases/e8-1-fix4-bundle-profile-repin.md`
- `docs/memory/manifest.json`

This reconciliation contract is also an intended tracked file:

- `docs/memory/phases/e8-1-fix4-source-review-reconciliation-task-contract.md`

Do not modify any other file.

In particular, the following must remain byte-for-byte unchanged during this
reconciliation:

- `testing/pion_hgcer_validation_bundle_profile_e8_1.json`
- `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`
- `docs/memory/phases/e8-1-fix3-procedure-pdf-layout.md`
- `docs/memory/phases/e8-1-fix4-bundle-profile-repin-task-contract.md`
- all `src/` files;
- collector and wrapper files;
- accepted F.6.2 artifacts and evidence.

## Required CURRENT.md reconciliation

Record:

- E.8.1 remains `ACTIVE`;
- E.8.1.Fix.1 and E.8.1.Fix.2 retain their narrow
  `CLOSED / RUNTIME VALIDATED` statuses;
- E.8.1.Fix.3 remains `SOURCE REVIEWED`, user-pushed at
  `350c34c55b2de33ad01011559dc6d8ed84d9c8a7`, with no farm/runtime validation;
- E.8.1.Fix.4 is now `SOURCE REVIEWED` from independent ChatGPT actual-diff
  inspection;
- Codex-reported Fix.4 checks were `NOT RUN by ChatGPT`;
- no ROOT/PyROOT/farm/runtime validation is claimed for Fix.4;
- no Fix.4 commit exists yet because the user has not committed/pushed it;
- F.6.3 remains `BLOCKED`;
- NEXT is user-controlled commit/push of the reviewed Fix.4 change set;
- after push, ChatGPT must inspect the actual pushed Fix.4 commit;
- that resulting pushed Fix.4 commit, distinct from the profile's required
  Fix.3 source `350c34c...`, becomes the wrapper `--bundle-commit`;
- only after pushed-state review should the user rerun
  `Q4p4W2p74 / Left / lowe`, collect a fresh bundle/PDF, and return it for
  independent visual inspection of the Fix.3 clipping, mojibake, handoff, and
  overlay-header repairs;
- only after that Left/lowe gate passes should validation return to detailed
  `Left / highe`, then broaden.

Do not upgrade E.8.1 itself and do not create a new farm-evidence record.

## Required Fix.4 phase-record reconciliation

Update:

`docs/memory/phases/e8-1-fix4-bundle-profile-repin.md`

so that:

- status becomes `SOURCE REVIEWED`;
- it records that ChatGPT inspected the complete actual diff and found the
  profile/test provenance re-pin PASS;
- Codex-reported deterministic checks are explicitly `NOT RUN by ChatGPT`;
- no ROOT/PyROOT/farm/runtime validation is claimed;
- user commit/push remains pending;
- the future Fix.4 profile/bundle commit remains distinct from required
  Fix.3 analysis source `350c34c...`;
- NEXT remains user commit/push, pushed-state review, then fresh targeted
  Left/lowe farm/PDF validation.

Preserve the existing provenance and scientific-boundary description.

## Validation

Regenerate `docs/memory/manifest.json`.

Run the applicable repository-memory manifest/integrity/health/bootstrap checks
and:

```text
git diff --check
```

Verify explicitly that the already reviewed profile/test files have no new diff
introduced during this reconciliation.

## Review bundle

Refresh one complete temporary root-level review bundle:

`kaonlt_review.diff`

It must contain:

1. the complete tracked diff for the entire current Fix.4 worktree;
2. complete `git diff --no-index /dev/null ...` representations for every
   intended new/untracked Fix.4 file, including this reconciliation contract.

Do not stage merely to create the review bundle.

## Hard stop

Do not commit, push, alter any analysis source, or run the Jefferson Lab farm.

Return:

- exact changed paths;
- exact memory/check results;
- confirmation that profile/test/analysis/collector/wrapper files were untouched
  during this pass;
- the refreshed root-level `kaonlt_review.diff` ready for ChatGPT review.
