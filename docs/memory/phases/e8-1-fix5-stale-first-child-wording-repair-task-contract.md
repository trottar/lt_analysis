# E.8.1.Fix.5 — stale first-child evidence wording repair

## Purpose

Repair one remaining internal inconsistency found during independent ChatGPT
review of the refreshed E.8.1.Fix.5 diff.

The Fix.5 renderer implementation and repaired real-ROOT regression are not
being reopened. The current source/test candidate remains otherwise correct.

The only defect is in the tracked original Fix.5 task contract: its earlier
"Independent PDF review" section still calls
`Left-lowe phi1 [-180, -140)` the first canonical child, even though the same
contract now correctly states later that the persisted `phi_index` contract is
zero based (`0..8`) and the real first child is
`Left-lowe phi0 [-180, -140)`.

## Exact committed base

Committed `test` HEAD must remain:

`fdd368f2084c9b9508e7ce679b8f7391f5b556f1`

The existing uncommitted E.8.1.Fix.5 worktree is intentional.

Before editing:

1. establish branch/HEAD/worktree;
2. read the normal five-file startup sequence;
3. read the Fix.5 task contract and actual-diff repair contract;
4. confirm no unrelated changes exist.

## Independent review result already established

The refreshed diff now correctly:

- keeps the Fix.5 renderer patch byte-for-byte unchanged from the prior review;
- uses the validated persisted E.8 authority in the conditional real-ROOT test;
- preserves zero-based child indices `0..8`;
- derives the first display label with `_e8_child_label`;
- explicitly requires `Left-lowe phi0 [-180, -140)`;
- records the farm evidence as:
  - `phi0 [-180, -140)` absent;
  - `phi1 [-140, -100)` absent;
  - later titles including `phi2 [-100, -60)` extractable;
- restores CURRENT lifecycle `BLOCKED / DEFERRED`;
- restores the Debug.1.Fix.1 provenance caveat and wrapper/debug references.

Do not alter any of those already-correct repairs.

## Allowed modifications

Only:

- `docs/memory/phases/e8-1-fix5-overlay-pdf-geometry-task-contract.md`
- `docs/memory/manifest.json`

This repair contract is also intended to be tracked:

- `docs/memory/phases/e8-1-fix5-stale-first-child-wording-repair-task-contract.md`

Do not modify any other file.

In particular, do not modify:

- `src/cuts/full_background_subtraction_plots.py`
- `testing/test_full_background_subtraction_plots.py`
- `docs/memory/CURRENT.md`
- `docs/memory/evidence/e8-1-fix4-left-lowe-overlay-blocker.md`
- `docs/memory/phases/e8-1-fix3-procedure-pdf-layout.md`
- `docs/memory/phases/e8-1-fix4-bundle-profile-repin.md`
- `docs/memory/phases/e8-1-fix5-overlay-pdf-geometry.md`
- `docs/memory/phases/e8-1-fix5-actual-diff-review-repair-task-contract.md`
- any profile, collector, wrapper, production source, or accepted artifact.

## Required wording correction

In the original Fix.5 task contract, under:

`### Remaining blocker: pages 38/41/44`

replace the stale sentence:

```text
On page 38, the first canonical child title
`Left-lowe phi1 [-180, -140)` is also absent, while later child titles are
extractable.
```

with the direct corrected evidence:

```text
On page 38, `Left-lowe phi0 [-180, -140)` and
`Left-lowe phi1 [-140, -100)` are absent, while later titles including
`Left-lowe phi2 [-100, -60)` are extractable.
```

Do not rewrite any other part of that contract.

This wording must agree exactly in substance with:

`docs/memory/evidence/e8-1-fix4-left-lowe-overlay-blocker.md`

and with the zero-based regression contract later in the same Fix.5 task file.

## Status

Keep:

- E.8.1 — `ACTIVE`
- E.8.1.Fix.5 — `ACTIVE`

Do not mark Fix.5 `SOURCE REVIEWED` during this repair. That status change, if
warranted, follows a fresh independent ChatGPT review of the complete actual
diff.

## Validation

Regenerate:

`docs/memory/manifest.json`

Run the applicable memory manifest/integrity/health/bootstrap checks and:

```text
git diff --check
```

No source/test rerun is required by this memory-only wording repair, but report
if any checks are rerun.

## Required diff audit

Before stopping:

1. confirm committed HEAD remains
   `fdd368f2084c9b9508e7ce679b8f7391f5b556f1`;
2. confirm the Fix.5 renderer patch is unchanged;
3. confirm the Fix.5 test patch is unchanged;
4. confirm CURRENT and the farm-evidence record are unchanged;
5. confirm the original task contract no longer contains the stale
   `phi1 [-180, -140)` first-child claim;
6. confirm no profile/collector/wrapper/production source changed;
7. refresh the complete root-level `kaonlt_review.diff`, including this new
   untracked contract through `git diff --no-index /dev/null ...`.

Do not stage merely to create the review bundle.

## Hard stop

Do not commit, push, alter the E.8.1 validation profile, or run the Jefferson
Lab farm.

Return the refreshed root-level `kaonlt_review.diff` for independent ChatGPT
review.
