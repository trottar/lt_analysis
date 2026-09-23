# E.8.1.Fix.4 — final CURRENT status wording repair

## Purpose

Repair one post-review durable-memory wording contradiction found by ChatGPT in
the refreshed E.8.1.Fix.4 review bundle.

The Fix.4 profile/test provenance re-pin itself remains independently source
reviewed and PASS. This task is memory-only and must not alter any implementation,
test, profile, analysis, collector, wrapper, evidence, or prior phase record.

## Starting state

Committed `test` HEAD remains:

`350c34c55b2de33ad01011559dc6d8ed84d9c8a7`

The existing uncommitted E.8.1.Fix.4 change set is intentional.

Before editing:

1. establish branch/HEAD/worktree;
2. read the normal five-file startup sequence;
3. confirm there are no unrelated worktree changes.

## Exact defect

`docs/memory/CURRENT.md` currently says, in the same Current Work Item
paragraph, both:

- Fix.4 is `ACTIVE` for the profile-only provenance re-pin; and
- Fix.4 is now `SOURCE REVIEWED`.

That is internally inconsistent after ChatGPT completed the actual-diff review.

## Required repair

Modify only:

- `docs/memory/CURRENT.md`
- `docs/memory/manifest.json`

This contract is also an intended tracked file:

- `docs/memory/phases/e8-1-fix4-final-current-status-repair-task-contract.md`

In `CURRENT.md`, remove the stale `Fix.4 is ACTIVE` wording from the Current Work
Item paragraph and leave one unambiguous state:

- E.8.1 remains `ACTIVE`;
- E.8.1.Fix.3 remains `SOURCE REVIEWED`, user-pushed at
  `350c34c55b2de33ad01011559dc6d8ed84d9c8a7`;
- E.8.1.Fix.4 is `SOURCE REVIEWED`;
- Codex-reported Fix.4 checks were `NOT RUN by ChatGPT`;
- no Fix.4 ROOT/PyROOT/farm/runtime validation is claimed;
- no Fix.4 commit exists yet;
- NEXT remains user commit/push -> ChatGPT pushed-state review -> fresh targeted
  `Q4p4W2p74 / Left / lowe` farm/PDF gate.

Do not change any other CURRENT semantics.

Do not modify:

- `testing/pion_hgcer_validation_bundle_profile_e8_1.json`
- `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`
- `docs/memory/phases/e8-1-fix3-procedure-pdf-layout.md`
- `docs/memory/phases/e8-1-fix4-bundle-profile-repin.md`
- either prior Fix.4 task contract
- any `src/` file
- collector or wrapper files
- evidence files
- accepted F.6.2 artifacts.

Regenerate `docs/memory/manifest.json`.

Run the applicable memory manifest/integrity/health/bootstrap checks and:

```text
git diff --check
```

Refresh the complete root-level `kaonlt_review.diff`, including this new
untracked contract via `git diff --no-index /dev/null ...`.

Do not stage merely for review.

## Hard stop

Do not commit, push, or run the farm.

Return the refreshed `kaonlt_review.diff` for final ChatGPT review.
