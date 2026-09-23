# E.8.1.Fix.5 — post-review source-status reconciliation

## Purpose

Reconcile durable repository memory after independent ChatGPT actual-diff review
of the complete E.8.1.Fix.5 change set.

The Fix.5 source/test candidate has now passed independent actual-diff review.
This task is memory-only. It must not alter the reviewed renderer, regression
tests, farm evidence, original Fix.5 task contract, any validation profile,
collector, wrapper, analysis/runtime source, or accepted F.6.2 artifact.

## Exact committed base

Committed `test` HEAD remains:

`fdd368f2084c9b9508e7ce679b8f7391f5b556f1`

The existing uncommitted E.8.1.Fix.5 worktree is intentional.

Before editing:

1. establish exact branch/HEAD/worktree;
2. read the normal five-file repository-memory startup sequence;
3. read the Fix.5 phase record and the current Fix.5 contracts;
4. confirm there are no unrelated worktree changes.

## Independent ChatGPT actual-diff result

ChatGPT inspected the complete refreshed Fix.5 review bundle and found the
source/test repair PASS.

The reviewed implementation:

- changes only E.8 persisted-overlay presentation geometry in
  `src/cuts/full_background_subtraction_plots.py`;
- changes the overlay canvas from `1800 x 3600` to `3600 x 3600`;
- preserves the distinct header and grid pads;
- preserves the exact `3 x 9` grid and all 27 tile positions in canonical
  child/variable order;
- draws/populates the grid first and draws the header pad afterward as the final
  top-level sibling;
- keeps the existing L/B/A legend in that header pad;
- preserves L/B/A colors, sparse proxy behavior, persisted values, annotations,
  scientific definitions, page IDs, and page order;
- returns to the parent canvas and calls `Modified()` / `Update()` before
  `Print()`;
- does not change context, handoff, map-page, reader-authority, profile,
  collector, wrapper, production, yield, cross-section, or accepted F.6.2
  behavior.

The reviewed tests:

- extend fake ROOT only enough to expose canvas dimensions and
  Modified/Update state;
- verify the square overlay canvas;
- verify grid/header draw order and non-overlapping pad geometry;
- retain the exact 27-tile ordering and L/B/A color/proxy checks;
- add a conditional real-PyROOT + `pdftotext` multipage-PDF regression;
- validate the actual persisted E.8 authority before rendering;
- preserve the canonical zero-based child index contract `0..8`;
- derive the first title through `_e8_child_label()` and require
  `Left-lowe phi0 [-180, -140)`;
- require the overlay header and all three L/B/A legend labels in extracted PDF
  text.

The durable evidence now correctly records that the Fix.4 farm PDF lacked
`phi0 [-180, -140)` and `phi1 [-140, -100)` on page 38 while later titles
including `phi2 [-100, -60)` were extractable.

The original Fix.5 task contract was corrected to the same zero-based evidence
wording. The separate repair contracts may retain the old string only as
historical description of the defect they repaired.

The refreshed CURRENT state also correctly preserves:

- E.8.1.Debug.1.Fix.1 as `SOURCE REVIEWED`;
- the caveat that earlier incomplete debug output does not establish launcher
  provenance or the intentional full-high-epsilon skip;
- lifecycle-hook dispatch as non-required and `BLOCKED / DEFERRED`;
- the wrapper/debug references needed for the next farm gate.

## Test/check boundary

Codex-reported deterministic checks remain `NOT RUN by ChatGPT`:

- `py_compile`: PASS;
- focused `FullBackgroundSubtractionE8Tests`: PASS, 12 tests with 1 expected
  skip because PyROOT is unavailable locally;
- complete `testing.test_full_background_subtraction_plots`: PASS, 102 tests
  with 18 expected skips;
- memory manifest/integrity/health/bootstrap checks: PASS;
- `git diff --check`: PASS.

The conditional real-PDF regression did not run locally because PyROOT is
unavailable.

No ROOT/PyROOT/farm/runtime acceptance is claimed for Fix.5. A fresh Jefferson
Lab farm PDF remains mandatory.

## Allowed files

Modify only:

- `docs/memory/CURRENT.md`
- `docs/memory/phases/e8-1-fix5-overlay-pdf-geometry.md`
- `docs/memory/manifest.json`

This reconciliation contract is also an intended tracked file:

- `docs/memory/phases/e8-1-fix5-source-review-reconciliation-task-contract.md`

Do not modify any other file.

In particular, the following must remain byte-for-byte unchanged during this
reconciliation:

- `src/cuts/full_background_subtraction_plots.py`
- `testing/test_full_background_subtraction_plots.py`
- `docs/memory/evidence/e8-1-fix4-left-lowe-overlay-blocker.md`
- `docs/memory/phases/e8-1-fix3-procedure-pdf-layout.md`
- `docs/memory/phases/e8-1-fix4-bundle-profile-repin.md`
- `docs/memory/phases/e8-1-fix5-overlay-pdf-geometry-task-contract.md`
- `docs/memory/phases/e8-1-fix5-actual-diff-review-repair-task-contract.md`
- `docs/memory/phases/e8-1-fix5-stale-first-child-wording-repair-task-contract.md`
- all profiles;
- collector and wrapper files;
- `run_Prod_Analysis.sh`;
- `src/main.py`;
- `src/cuts/rand_sub.py`;
- accepted F.6.2 artifacts.

## Required CURRENT.md reconciliation

Record one unambiguous current state:

- E.8.1 remains `ACTIVE`;
- E.8.1.Fix.1 remains `CLOSED / RUNTIME VALIDATED`;
- E.8.1.Fix.2 remains `CLOSED / RUNTIME VALIDATED`;
- E.8.1.Fix.3 remains `SOURCE REVIEWED`; its context/handoff repair is
  farm-confirmed but the old persisted-overlay subrepair was superseded by
  Fix.5;
- E.8.1.Fix.4 remains `CLOSED / RUNTIME VALIDATED` only for its narrow
  profile-provenance re-pin;
- E.8.1.Fix.5 is now `SOURCE REVIEWED` from independent ChatGPT actual-diff
  inspection;
- Codex-reported Fix.5 checks were `NOT RUN by ChatGPT`;
- no Fix.5 ROOT/PyROOT/farm/runtime validation is claimed;
- no Fix.5 commit exists yet because the user has not committed/pushed it;
- E.8.1.Debug.1.Fix.1 remains `SOURCE REVIEWED` with its launcher-provenance
  caveat;
- F.6.3 remains `BLOCKED` pending E.8.1;
- F.6.4 remains `BLOCKED` pending F.6.3 evidence;
- lifecycle-hook dispatch remains `BLOCKED / DEFERRED`.

Set exact NEXT to:

1. user-controlled commit/push of the reviewed Fix.5 change set;
2. ChatGPT pushed-state review;
3. a separate narrow E.8.1 profile re-pin to the pushed Fix.5 source;
4. independent actual-diff review of that profile re-pin;
5. user commit/push of the profile re-pin;
6. ChatGPT pushed-state review;
7. fresh targeted `Q4p4W2p74 / Left / lowe` farm run and validation bundle;
8. independent inspection of pages 37--47, with pages 38/41/44 required to show
   the header, L/B/A legend, first canonical child title, and unclipped first
   row;
9. only after Left/lowe visual PASS return to detailed `Left / highe`;
10. broaden only after that detailed gate passes.

Do not upgrade E.8.1 itself.

Do not invent the future Fix.5 source commit or future profile/bundle commit.

## Required Fix.5 phase-record reconciliation

Update:

`docs/memory/phases/e8-1-fix5-overlay-pdf-geometry.md`

so that:

- status becomes `SOURCE REVIEWED`;
- it records that ChatGPT inspected the complete refreshed actual diff and
  found the narrow renderer/test repair PASS;
- it states that Codex-reported deterministic checks were
  `NOT RUN by ChatGPT`;
- it preserves the presentation-only/scientific ownership boundary;
- it explicitly states that the conditional real-PyROOT PDF test did not run
  locally because PyROOT is unavailable;
- it makes no ROOT/PyROOT/farm/runtime acceptance claim;
- user commit/push remains pending;
- NEXT is user commit/push -> pushed-state review -> separate profile re-pin ->
  fresh Left/lowe farm/PDF validation.

Preserve the farm blocker and repair chronology.

## Validation

Regenerate:

`docs/memory/manifest.json`

Run applicable repository-memory manifest/integrity/health/bootstrap checks and:

```text
git diff --check
```

Verify explicitly that the reviewed source/test/evidence/contracts have no new
diff introduced during this reconciliation.

## Review bundle

Refresh one complete temporary root-level:

`kaonlt_review.diff`

It must contain:

1. the complete tracked diff for the entire current Fix.5 worktree;
2. complete `git diff --no-index /dev/null ...` representations for every
   intended new/untracked Fix.5 file, including this reconciliation contract.

Do not stage merely to create the review bundle.

## Hard stop

Do not commit, push, alter any E.8.1 validation profile, or run Jefferson Lab
farm validation.

Return:

- exact changed paths;
- exact memory/check results;
- confirmation that reviewed source/test/evidence/profile/collector/wrapper
  files were untouched during this reconciliation;
- refreshed root-level `kaonlt_review.diff` ready for final ChatGPT review.
