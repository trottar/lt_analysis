# E.8.1.Fix.3 — post-review status reconciliation

## Purpose

Reconcile durable repository memory after independent ChatGPT actual-diff review
of E.8.1.Fix.3.

The implementation/test diff has now been independently reviewed and is a PASS.
This task is memory-only. It must not alter the E.8 presentation implementation,
its tests, the accepted F.6.2 artifact, any analysis/runtime source, the E.8.1
profile, collector, or bundle wrapper.

## Starting committed state

Committed `test` HEAD remains:

`9f7094b19f3cdb0887c95f1f972d98471f2ecfdc`

The existing uncommitted E.8.1.Fix.3 implementation/test/memory changes are
intentional.

Before editing:

1. establish branch/HEAD/worktree;
2. read the normal five-file repository-memory startup sequence;
3. read the current Fix.3 task contract and phase record;
4. confirm there are no unrelated worktree changes.

## Independent review result

ChatGPT inspected the complete actual `kaonlt_review.diff` and found the
E.8.1.Fix.3 source candidate PASS.

The reviewed implementation:

- keeps the E.8 reader/authority/scientific payload path unchanged;
- adds an optional presentation-only text size to `_e8_text_page`;
- replaces the E.8 context-page Unicode em dash with renderer-safe ASCII;
- splits the context and handoff prose into explicit bounded display lines;
- retains the complete frozen-input SHA-256, artifact fingerprint, validation
  fingerprint, L/B/A semantics, non-production statement, and F.6.3 handoff;
- introduces a distinct overlay header `TPad` above a grid `TPad`;
- divides only the grid pad into the unchanged `3 x 9` canonical overlay grid;
- renders all 27 child/variable tiles through that grid pad;
- draws the existing L/B/A legend in the reserved header pad;
- preserves L/B/A colors and sparse-state proxy behavior;
- does not change map-page semantics, page IDs/order/inventory, persisted
  numerical values, production weights, yields, corrections, or cross sections.

The focused regression changes exercise the bounded text/provenance content,
ASCII-safe context presentation, separate header/grid geometry, complete
27-tile ordering, legend placement, and sparse proxy colors.

The actual changed paths are all within the Fix.3 task allowlist.

Codex-reported deterministic checks remain `NOT RUN by ChatGPT`:

- `py_compile`: PASS;
- focused E.8 tests: PASS, 11 tests;
- complete `testing.test_full_background_subtraction_plots`: PASS, 101 tests,
  17 expected skips because PyROOT is unavailable locally;
- memory manifest/health/bootstrap checks: PASS.

No ROOT/PyROOT/farm/runtime validation is claimed for Fix.3. The repaired visual
layout still requires a fresh Jefferson Lab farm PDF and independent visual
inspection.

## Allowed files

Modify only:

- `docs/memory/CURRENT.md`
- `docs/memory/phases/e8-1-fix3-procedure-pdf-layout.md`
- `docs/memory/manifest.json`

This reconciliation contract is also an intended tracked file:

- `docs/memory/phases/e8-1-fix3-source-review-reconciliation-task-contract.md`

Do not modify any other file.

In particular, the following must remain byte-for-byte unchanged during this
reconciliation:

- `src/cuts/full_background_subtraction_plots.py`
- `testing/test_full_background_subtraction_plots.py`
- `docs/memory/evidence/e8-1-fix2-left-lowe-layout-blocker.md`
- `docs/memory/phases/e8-1-fix1-frozen-f6-2-setting-reader.md`
- `docs/memory/phases/e8-1-fix2-bundle-profile-repin.md`
- `docs/memory/phases/e8-1-fix3-procedure-pdf-layout-task-contract.md`
- all profiles, collectors, wrappers, and accepted artifacts.

## Required CURRENT.md reconciliation

Record:

- E.8.1 remains `ACTIVE`;
- E.8.1.Fix.1 remains `CLOSED / RUNTIME VALIDATED` only for its narrow reader
  defect;
- E.8.1.Fix.2 remains `CLOSED / RUNTIME VALIDATED` only for its narrow profile
  provenance defect;
- E.8.1.Fix.3 is now `SOURCE REVIEWED` from independent ChatGPT actual-diff
  inspection;
- Codex-reported checks were `NOT RUN by ChatGPT`;
- no ROOT/PyROOT/farm/runtime validation is claimed for Fix.3;
- no Fix.3 commit exists yet because the user has not committed/pushed it;
- F.6.3 remains `BLOCKED` pending E.8.1;
- E.8.1.Debug.1.Fix.1 remains `SOURCE REVIEWED`;
- NEXT is user-controlled commit/push of the reviewed Fix.3 change set;
- after push, ChatGPT must inspect the actual pushed commit;
- after pushed-state review, create a separate narrow E.8.1 profile re-pin to
  the pushed Fix.3 source commit;
- only after that separately reviewed profile push should the user rerun
  `Q4p4W2p74 / Left / lowe`, collect a fresh bundle, and return the PDF/bundle
  for independent visual review;
- only after the repaired Left/lowe PDF passes should validation return to
  detailed `Left / highe`, then broaden.

Do not upgrade E.8.1 itself and do not create a Fix.3 farm-evidence record.

## Required Fix.3 phase-record reconciliation

Update:

`docs/memory/phases/e8-1-fix3-procedure-pdf-layout.md`

so that:

- status becomes `SOURCE REVIEWED`;
- it records that ChatGPT inspected the complete actual diff and found the
  presentation-only source/test repair PASS;
- it states that Codex-reported deterministic checks were `NOT RUN by ChatGPT`;
- it preserves the scientific/presentation boundary;
- it makes no ROOT/PyROOT/farm/runtime claim;
- user commit/push remains pending;
- NEXT remains user commit/push, pushed-state review, separate profile re-pin,
  then fresh Left/lowe farm PDF validation.

Preserve the blocker description and local-test chronology.

## Validation

Regenerate `docs/memory/manifest.json`.

Run the applicable repository-memory manifest/integrity/health/bootstrap checks
and:

```text
git diff --check
```

Verify explicitly that the already reviewed implementation/test files have no
new diff introduced during this reconciliation.

## Review bundle

Refresh one complete temporary root-level review bundle:

`kaonlt_review.diff`

It must contain:

1. the complete tracked diff for the entire current Fix.3 worktree;
2. complete `git diff --no-index /dev/null ...` representations for every
   intended new/untracked Fix.3 file, including this reconciliation contract.

Do not stage merely to create the review bundle.

## Hard stop

Do not commit, push, alter any E.8.1 profile, or run the Jefferson Lab farm.

Return:

- exact changed paths;
- exact memory/check results;
- confirmation that implementation/test/profile/collector/wrapper files were
  untouched during this pass;
- the refreshed root-level `kaonlt_review.diff` ready for ChatGPT review.
