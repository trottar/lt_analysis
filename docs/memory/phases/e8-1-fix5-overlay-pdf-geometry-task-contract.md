# E.8.1.Fix.5 — repair remaining farm-visible persisted-overlay page clipping

## Objective

Repair the one remaining E.8.1 farm-visible presentation blocker after the
Fix.3 text/layout repair and Fix.4 provenance re-pin.

The fresh `Q4p4W2p74 / Left / lowe` farm bundle proves that the E.8 context
page, handoff page, profile provenance, frozen F.6.2 identity, page inventory,
and renderer/checker gates now work. The remaining defect is isolated to the
three E.8 persisted-overlay pages: their reserved header/legend is absent from
the rendered PDF and their first plot row is clipped at the top page boundary.

This task is presentation-only. It must not change any accepted F.6.2 value,
E.8 reader authority, population definition, page inventory, map page,
production object, weight, yield, correction, or cross section.

## Exact starting state

Start only from committed `test` HEAD:

`fdd368f2084c9b9508e7ce679b8f7391f5b556f1`

Commit message:

`Re-pin E8.1 profile to layout repair`

Its parent is the pushed Fix.3 renderer source:

`350c34c55b2de33ad01011559dc6d8ed84d9c8a7`

Before editing:

1. establish exact branch/HEAD/worktree;
2. read the normal five-file repository-memory startup sequence;
3. read the Fix.3 and Fix.4 phase records;
4. read the fresh Fix.4 Left/lowe evidence record created by this task once its
   facts are entered;
5. confirm this task-contract file is the only expected user-created
   uncommitted file;
6. stop if HEAD differs or unrelated worktree changes exist.

## Fresh farm evidence that owns this repair

Bundle:

`KaonLT_E8_1_Fix4_Q4p4W2p74_Left_lowe_20260923-062807.zip`

Direct archive inspection established:

- ZIP SHA-256:
  `794c719d994444b49a58b0b9b326d5de019305a2e876615bf989f8c11bc01d59`
- bundle schema:
  `pion_hgcer_validation_bundle/v4`
- validation profile:
  `phase_e8_1_full_background_procedure_pdf_farm_review/v1`
- bundle `git_head`:
  `fdd368f2084c9b9508e7ce679b8f7391f5b556f1`
- required analysis/procedure source:
  `350c34c55b2de33ad01011559dc6d8ed84d9c8a7`
- required analysis commit is an ancestor:
  `true`
- unexpected committed files after required analysis commit:
  none
- manifest:
  `complete=true`
- manifest:
  `errors=[]`
- source-state worktree:
  clean
- frozen F.6.2 JSON SHA-256:
  `5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1`
- accepted F.6.2 artifact fingerprint:
  `ee713b70de898bad8fa61164cbc8a54886af1ea1eb1e712ed95de17df8890cd0`
- accepted F.6.2 validation fingerprint:
  `7edc73fce20ad7dc7622b8c23a7ba7e8986e367ace605884e010945595370b3b`
- fresh Left/lowe procedure-PDF SHA-256:
  `9a00fc6179891e936523613293dc40ebe5790f7ac49acf37aa942aea0333b98c`
- fresh Left/lowe page-manifest SHA-256:
  `d4bc89f49ce374f63055afd974d734bf782cea6e55cefcc26caf3c2f83028dc1`
- page manifest:
  47 pages
- `renderer_failures=[]`
- E.8 page inventory remains exactly:
  - page 37 `full_background.e8.context`
  - pages 38/41/44 `full_background.e8.persisted_overlays`
  - pages 39/42/45 `full_background.e8.delta_acceptance_maps`
  - pages 40/43/46 `full_background.e8.mm_acceptance_maps`
  - page 47 `full_background.e8.handoff`
- no `full_background.e8.unavailable` page exists.

Farm bundle source checks passed:

- `py_compile`: PASS
- `testing.test_full_background_subtraction_plots`:
  PASS, 101 tests, 14 skips
- Method-A acceptance contract:
  PASS, 7 tests
- Phase-F runtime contract:
  PASS, 3 tests
- validation-bundle collector:
  PASS, 28 tests
- `git diff --check`:
  PASS
- required-analysis-commit range diff check:
  PASS

These facts validate the narrow Fix.4 provenance re-pin. They do not validate
the remaining overlay-page visual layout.

## Independent PDF review

The fresh farm PDF was rendered independently with both PDFium and Poppler.

### Confirmed PASS: page 37 context

`full_background.e8.context` now renders cleanly:

- no right-edge clipping;
- no Unicode em-dash mojibake;
- complete 64-character frozen-input SHA-256 visible;
- complete artifact fingerprint visible;
- complete validation fingerprint visible;
- L/B/A definitions visible;
- non-production statement visible;
- F.6.3 handoff statement visible.

Independent `pdftotext` extraction also contains those complete strings.

### Confirmed PASS: page 47 handoff

`full_background.e8.handoff` now renders cleanly:

- the formerly clipped D.10-omission sentence is fully visible;
- the no-Method-A-production-promotion statement is visible;
- the F.6.3 handoff is visible;
- the independent-farm-review requirement is visible.

### Confirmed PASS: map pages

Pages 39/40/42/43/45/46 do not show the overlay-header blocker. Do not redesign
them in this task.

### Remaining blocker: pages 38/41/44

All three `full_background.e8.persisted_overlays` pages still fail the visual
gate in both PDFium and Poppler:

- the reserved E.8 header/title is absent from the rendered page;
- the three-column L/B/A legend is absent;
- the first plot row is clipped at the physical top page boundary.

`pdftotext` confirms that the expected header/legend text is absent from pages
38, 41, and 44. On page 38, `Left-lowe phi0 [-180, -140)` and
`Left-lowe phi1 [-140, -100)` are absent, while later titles including
`Left-lowe phi2 [-100, -60)` are extractable.

This is a real farm-output defect, not a viewer artifact.

## Source-level diagnosis boundary

Current reviewed source in
`src/cuts/full_background_subtraction_plots.py` already creates:

- a header `TPad` at NDC `(0.0, 0.89, 1.0, 1.0)`;
- a grid `TPad` at NDC `(0.0, 0.0, 1.0, 0.89)`;
- a `3 x 9` grid inside the grid pad;
- the L/B/A legend in the header pad.

The existing fake-ROOT regression verifies those object relationships, but the
farm PDF proves that this object-level check is insufficient to guarantee the
actual multipage ROOT-PDF output.

A concrete source difference is that the overlay page uses an unusually tall
`1800 x 3600` canvas, while the E.8 map pages use a square `3600 x 3600`
canvas and render without this top-page failure.

Treat the canvas-aspect explanation as the leading source hypothesis, not as a
scientific fact. The required behavior is the farm-visible result below.

## Status consequences

Preserve/update durable status as follows:

- F.6.2 — `CLOSED / RUNTIME VALIDATED`
- F.6.2.Fix.5 — `CLOSED / RUNTIME VALIDATED`
- E.8.1 — `ACTIVE`
- E.8.1.Fix.1 — `CLOSED / RUNTIME VALIDATED`
- E.8.1.Fix.2 — `CLOSED / RUNTIME VALIDATED`
- E.8.1.Fix.3 — `SOURCE REVIEWED`; fresh farm evidence confirms its context and
  handoff subrepairs but disproves completion of its persisted-overlay visual
  subrepair
- E.8.1.Fix.4 — `CLOSED / RUNTIME VALIDATED` only for the narrow profile
  provenance re-pin to Fix.3 source
- E.8.1.Fix.5 — `ACTIVE` until independent ChatGPT actual-diff review
- E.8.1.Debug.1.Fix.1 — `SOURCE REVIEWED`
- F.6.3 — `BLOCKED` pending E.8.1

Do not close E.8.1.

## Scientific ownership

This task owns only ROOT presentation geometry for
`full_background.e8.persisted_overlays`.

Frozen scientific content includes:

- accepted F.6.2 JSON;
- all F.6.2 fingerprints and authority checks;
- canonical t/phi geometry;
- L/B/A population definitions;
- persisted one-dimensional distributions;
- DeltaH, kappa, DeltaP^K, f_refine^K and all other stored values;
- sparse/unavailable child states;
- kaon missing-mass window;
- normalization semantics.

Method A remains detached/non-production.
Method B remains diagnostic/cross-check only.

No production object or correction may change.

## Allowed implementation files

Only:

- `src/cuts/full_background_subtraction_plots.py`
- `testing/test_full_background_subtraction_plots.py`

## Allowed durable-memory files

Update only as warranted:

- `docs/memory/CURRENT.md`
- `docs/memory/phases/e8-1-fix3-procedure-pdf-layout.md`
- `docs/memory/phases/e8-1-fix4-bundle-profile-repin.md`
- `docs/memory/phases/e8-1-fix5-overlay-pdf-geometry.md`
- `docs/memory/evidence/e8-1-fix4-left-lowe-overlay-blocker.md`
- `docs/memory/manifest.json`

This task contract is also intended to be tracked:

- `docs/memory/phases/e8-1-fix5-overlay-pdf-geometry-task-contract.md`

Regenerate `docs/memory/manifest.json` whenever versioned memory changes.

Do not update unrelated durable-memory owners.

## Frozen files and interfaces

Do not modify:

- `run_Prod_Analysis.sh`;
- `src/main.py`;
- `src/cuts/rand_sub.py`;
- E.8 reader/authority constants;
- `_e8_validate_artifact`;
- `_e8_validate_parent_setting`;
- E.8 parent/child scientific validation;
- context-page text/layout that just passed the farm visual gate;
- handoff-page text/layout that just passed the farm visual gate;
- `_e8_render_map_page` or map-page semantics;
- page IDs or page order;
- E.8.1 validation profile;
- E.8.1 profile test;
- generic collector;
- bundle wrapper;
- any F.6.2 source/artifact/profile;
- any production analysis source.

## Required overlay-page repair

Keep exactly one persisted-overlay page per canonical t parent.

Keep exactly:

- three variables:
  `analysis_MM`, `SHMS_xptar`, `SHMS_yptar`;
- nine canonical phi children;
- 27 tile positions;
- existing child order;
- existing L/B/A colors;
- existing sparse/unavailable displays;
- existing metrics and annotations;
- existing page ID:
  `full_background.e8.persisted_overlays`.

### Canvas geometry

Replace the current extreme `1800 x 3600` overlay canvas with a PDF-safe aspect
ratio.

Preferred exact repair:

- use `3600 x 3600`, matching the farm-successful E.8 map-page canvas family.

Do not change map-page dimensions.

Retain a distinct header pad and grid pad. The grid must not overlap the header.

Retain the current 3x9 grid.

The implementation may also explicitly return to the parent canvas and call
`Modified()` / `Update()` before `Print()` if needed for stable ROOT batch
layout, but do not add unrelated drawing behavior.

If draw order matters for ROOT output, prefer:

1. draw the grid pad;
2. draw/populate the header pad as a separate top-level sibling;
3. ensure the header is the final top-level overlay before printing.

Do not solve this by drawing the legend inside a data tile or by deleting
content.

## Required farm-visible output

The repaired persisted-overlay page must contain, visibly and extractably:

- header:
  `E.8 persisted L/B/A overlays - tN`
- L key:
  `L: upstream 0 < NPE <= 2 diagnostic reference`
- B key:
  `B: physical pion control, NPE > 2, baseline w0`
- A key:
  the existing same-B-population `w0*C` label
- first canonical child title on each t page;
- all remaining 26 tiles.

The first row must sit below the header region and must not touch or cross the
physical PDF page boundary.

## Required regression tests

Extend `testing/test_full_background_subtraction_plots.py`.

### Fake-ROOT structural regression

Keep the existing geometry/ordering assertions and add explicit proof that:

- overlay canvas dimensions are the new PDF-safe dimensions;
- header and grid remain separate top-level pads;
- grid remains `3 x 9`;
- all 27 tile positions remain in exact canonical order;
- header/title and legend remain attached to the header pad;
- colors and sparse proxies are unchanged.

### Real-ROOT PDF regression

The existing fake-ROOT test did not catch this farm failure.

Add a narrowly scoped conditional real-ROOT regression that runs when the local
environment has PyROOT and a usable PDF text extractor such as `pdftotext`.

The regression should render the persisted-overlay page in a multipage-PDF
context representative of the ordinary procedure PDF and assert that the
rendered overlay page contains at least:

- `E.8 persisted L/B/A overlays - t1`
- the L legend label;
- the B legend label;
- the A legend label;
- the first canonical child title, e.g.
  `Left-lowe phi0 [-180, -140)`.

The canonical persisted phi_index contract is zero based (0..8); regression
fixtures must not renumber it.

The test must fail with the current farm-defective output behavior.

If PyROOT or the extractor is unavailable locally, skip explicitly and report
that skip. Do not make local source acceptance depend on unavailable farm
software.

Do not use OCR.

## Required deterministic local validation

Run at minimum:

```text
python -m py_compile \
  src/cuts/full_background_subtraction_plots.py \
  testing/test_full_background_subtraction_plots.py
```

Run focused E.8 tests and then:

```text
python -m unittest testing.test_full_background_subtraction_plots -v
```

Report exact test and skip counts.

Run normal repository-memory manifest/integrity/health/bootstrap checks and:

```text
git diff --check
```

These checks do not establish farm acceptance.

## Required evidence record

Create:

`docs/memory/evidence/e8-1-fix4-left-lowe-overlay-blocker.md`

Record only directly established facts from the fresh bundle and independent
PDF review.

Include:

- bundle filename and ZIP SHA-256;
- bundle/profile source identities;
- complete/errors/provenance result;
- frozen JSON SHA-256;
- artifact/validation fingerprints;
- PDF SHA-256;
- page-manifest SHA-256;
- 47-page inventory and empty renderer failures;
- farm test counts;
- page 37 PASS;
- page 47 PASS;
- map-page PASS/no new blocker;
- pages 38/41/44 overlay FAIL in both PDFium and Poppler;
- absent overlay header/legend in `pdftotext`;
- absent first child title on page 38;
- E.8.1 remains ACTIVE;
- Fix.4 narrow provenance gate may close, but Fix.3 visual work does not.

Do not invent launcher-log facts that are not in the supplied ZIP.

## Required actual-diff audit

Before stopping:

1. confirm committed HEAD remains
   `fdd368f2084c9b9508e7ce679b8f7391f5b556f1`;
2. show `git status --short`;
3. show `git diff --stat`;
4. show `git diff --name-only`;
5. confirm every changed path is allowlisted;
6. confirm no profile/collector/wrapper/production source changed;
7. trace the real renderer path for the persisted-overlay page;
8. confirm context/handoff/map renderers are unchanged;
9. confirm page IDs/order/inventory are unchanged.

Create one complete temporary root-level:

`kaonlt_review.diff`

containing the tracked diff plus complete `git diff --no-index /dev/null ...`
representations for every intended new/untracked file.

Do not stage merely to create it.

## NEXT after source review

If ChatGPT independently reviews the Fix.5 actual diff and passes it:

1. user commits/pushes Fix.5;
2. ChatGPT reviews the actual pushed Fix.5 source;
3. create a separate narrow E.8.1 profile re-pin to the pushed Fix.5 source;
4. user commits/pushes that profile re-pin after source review;
5. ChatGPT reviews the pushed profile commit;
6. rerun targeted `Q4p4W2p74 / Left / lowe`;
7. collect a fresh Left/lowe bundle;
8. inspect pages 37-47 independently;
9. only after persisted-overlay pages 38/41/44 pass return to detailed
   `Left / highe`;
10. broaden only after that detailed gate passes.

Do not encode a future Fix.5 commit into the profile in this task.

## Forbidden shortcuts

Do not:

- regenerate or modify accepted F.6.2 JSON;
- change scientific values to improve plots;
- remove sparse/empty phi children;
- remove the header or legend;
- reduce the 27-tile inventory;
- split the overlay into extra pages;
- change page IDs/order;
- change map pages merely for style;
- alter context/handoff pages that now pass;
- broaden the profile allowlist;
- modify collector/wrapper behavior;
- change production physics;
- commit, push, or run the farm.

## Acceptance criteria

Source-review candidate PASS requires:

- new overlay canvas uses a PDF-safe aspect ratio;
- header and L/B/A legend remain separate from the 3x9 grid;
- all 27 tiles remain;
- no scientific content changes;
- context/handoff/map code remains frozen;
- fake-ROOT structural tests pass;
- real-ROOT/pdf-text regression exists and runs when dependencies are available;
- local deterministic tests pass to capability;
- memory accurately records Fix.4 provenance closure and the remaining overlay
  blocker;
- actual diff is narrow and allowlisted.

## Hard stop

Do not commit, push, alter the E.8.1 validation profile, or run Jefferson Lab
farm validation.

Return to ChatGPT with:

- implementation summary;
- exact changed files;
- exact local checks/tests and skips;
- renderer-path trace;
- complete root-level `kaonlt_review.diff`;
- proposed user Git handoff commands, not executed.

ChatGPT must inspect the actual diff before user commit/push.
