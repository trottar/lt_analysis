# E.8.1.Fix.2 — post-review status reconciliation

## Purpose

Reconcile repository memory after independent ChatGPT actual-diff review of
E.8.1.Fix.2.

The implementation/provenance diff has now been reviewed and is a PASS. This
task is memory-only. It must not alter the E.8.1 profile, its test, any analysis
source, the collector, wrapper, or accepted artifacts.

## Starting committed state

Committed `test` HEAD remains:

`8985d9a212799c021c4ad1a759a689ea3826e0ea`

The existing uncommitted E.8.1.Fix.2 profile/test/memory changes are intentional.

Before editing:

1. establish branch/HEAD/worktree;
2. read the normal five-file repository-memory startup sequence;
3. read the current Fix.2 task contract and phase record;
4. confirm there are no unrelated worktree changes.

## Independent review result

ChatGPT inspected the complete actual review diff and found:

- the E.8.1 profile changes only
  `source_identity.required_analysis_commit`, from
  `0bf1281ebf44868b68ff808090436653eea6ed60`
  to
  `8985d9a212799c021c4ad1a759a689ea3826e0ea`;
- canonical five-setting membership is unchanged;
- artifact declarations are unchanged;
- `allowed_committed_files` is unchanged and remains limited to the profile and
  its focused test;
- `docs/memory/` remains the only allowed non-analysis prefix;
- no `src/` file, collector, wrapper, analysis source, accepted F.6.2 artifact,
  or production physics changed;
- the focused test changes only the exact reviewed-source constant;
- the profile remains fail-closed for later analysis changes.

Therefore E.8.1.Fix.2 may now be recorded as `SOURCE REVIEWED`.

Codex-reported deterministic checks remain `NOT RUN by ChatGPT`. No
ROOT/PyROOT/farm/runtime validation is claimed.

## Allowed files

Modify only:

- `docs/memory/CURRENT.md`
- `docs/memory/phases/e8-1-fix2-bundle-profile-repin.md`
- `docs/memory/manifest.json`

This reconciliation contract is also an intended tracked file:

- `docs/memory/phases/e8-1-fix2-source-review-reconciliation-task-contract.md`

Do not modify any other file.

In particular, the following must remain byte-for-byte unchanged during this
reconciliation:

- `testing/pion_hgcer_validation_bundle_profile_e8_1.json`
- `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`
- `docs/memory/phases/e8-1-fix2-bundle-profile-repin-task-contract.md`
- `docs/memory/phases/e8-1-fix1-frozen-f6-2-setting-reader.md`
- all `src/` files
- collector and wrapper files.

## Required CURRENT.md reconciliation

Record:

- E.8.1 remains `ACTIVE`;
- E.8.1.Fix.1 remains `SOURCE REVIEWED`, user-pushed at
  `8985d9a212799c021c4ad1a759a689ea3826e0ea`;
- E.8.1.Fix.2 is now `SOURCE REVIEWED` from independent ChatGPT actual-diff
  inspection;
- Codex-reported tests/checks were `NOT RUN by ChatGPT`;
- no ROOT/PyROOT/farm/runtime validation is claimed for Fix.2;
- F.6.2 and F.6.2.Fix.5 remain `CLOSED / RUNTIME VALIDATED`;
- F.6.3 remains `BLOCKED`;
- no profile/bundle commit exists yet because the user has not committed/pushed
  Fix.2;
- NEXT is user-controlled commit/push of the reviewed profile re-pin;
- after that push, ChatGPT must inspect the pushed commit;
- the resulting pushed Fix.2 commit, not the Fix.1 analysis source commit,
  becomes the wrapper `--bundle-commit` for the fresh
  `Q4p4W2p74 / Left / lowe` bundle-only farm gate;
- only after that narrow Left/lowe gate passes should review return to detailed
  `Left / highe` and then broaden.

Do not upgrade E.8.1 itself or create a farm-evidence record.

## Required Fix.2 phase-record reconciliation

Update
`docs/memory/phases/e8-1-fix2-bundle-profile-repin.md` so that:

- status becomes `SOURCE REVIEWED`;
- it records that ChatGPT inspected the complete actual diff and found the
  profile/test provenance repair PASS;
- Codex-reported deterministic checks are explicitly `NOT RUN by ChatGPT`;
- no farm/runtime claim is made;
- user commit/push remains pending;
- the future pushed Fix.2 commit remains distinct from the profile's required
  analysis commit `8985d9a...`;
- NEXT remains user commit/push, pushed-state review, then Left/lowe bundle-only
  farm validation.

Preserve the existing scientific/provenance description.

## Validation

Regenerate `docs/memory/manifest.json`.

Run the applicable repository-memory manifest/integrity/health/bootstrap checks
and:

```text
git diff --check
```

Verify that the already reviewed profile/test files have no new diff introduced
during this reconciliation.

## Review bundle

Refresh one complete temporary root-level review bundle:

`kaonlt_review.diff`

It must contain:

1. the complete tracked diff for the entire current Fix.2 worktree;
2. complete `git diff --no-index /dev/null ...` representations for every
   intended new/untracked Fix.2 file, including this reconciliation contract.

Do not stage merely to create the review bundle.

## Hard stop

Do not commit, push, or run the farm.

Return:

- exact changed paths;
- exact memory/check results;
- confirmation that profile/test/source files were untouched during this pass;
- the root-level `kaonlt_review.diff` ready for ChatGPT review.
