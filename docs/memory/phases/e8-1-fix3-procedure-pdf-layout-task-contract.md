# E.8.1.Fix.3 — repair farm-visible E.8 procedure-PDF layout defects

## Objective

Repair only the farm-visible E.8 presentation defects exposed after the
E.8.1.Fix.1 reader repair and E.8.1.Fix.2 provenance re-pin succeeded.

This task is presentation-only. It must not change the accepted F.6.2 artifact,
the E.8 reader's scientific/provenance validation, persisted values, page
inventory, production analysis physics, or any correction logic.

## Exact starting state

Start only from committed `test` HEAD:

`9f7094b19f3cdb0887c95f1f972d98471f2ecfdc`

Commit message:

`Re-pin E8.1 validation profile to repaired source`

Its parent is:

`8985d9a212799c021c4ad1a759a689ea3826e0ea`

Before editing:

1. establish branch/HEAD/worktree;
2. read the normal five-file repository-memory startup sequence;
3. read the active E.8.1.Fix.1/Fix.2 phase records;
4. confirm this task-contract file is the only expected user-created uncommitted
   file;
5. stop if HEAD differs or unrelated worktree changes exist.

## Fresh farm evidence that owns this repair

The fresh targeted bundle is:

`KaonLT_E8_1_Fix2_Q4p4W2p74_Left_lowe_20260923-040427.zip`

Direct archive inspection established:

- ZIP SHA-256:
  `52e0d283acdb57b56ed4c34c59f98dfe077f2a89e47fd0ced60fdcc91d5a5fd4`
- manifest `complete=true`
- manifest `errors=[]`
- bundle `git_head`:
  `9f7094b19f3cdb0887c95f1f972d98471f2ecfdc`
- required analysis commit:
  `8985d9a212799c021c4ad1a759a689ea3826e0ea`
- required analysis commit is an ancestor
- unexpected committed files after required analysis commit: none
- frozen F.6.2 JSON SHA-256 remained exactly:
  `5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1`
- fresh Left/lowe procedure-PDF SHA-256:
  `0745172508e2131e9256425b05bc44b19ffa6bcaa649224328679a75b7e18007`
- fresh Left/lowe page-manifest SHA-256:
  `d4bc89f49ce374f63055afd974d734bf782cea6e55cefcc26caf3c2f83028dc1`
- page manifest has `renderer_failures=[]`
- page manifest contains 47 pages:
  - 36 retained D.6-D.9 pages;
  - `full_background.e8.context`;
  - three E.8 pages for each of t1, t2, t3:
    `persisted_overlays`, `delta_acceptance_maps`, `mm_acceptance_maps`;
  - `full_background.e8.handoff`.
- there is no `full_background.e8.unavailable` page.

The bundle source checks also passed the repaired
`testing.test_full_background_subtraction_plots` farm module:
100 tests, 14 skips, return code 0.

This means the Fix.1 reader now accepts the actual frozen F.6.2 artifact on the
farm, and the Fix.2 profile provenance gate now collects it successfully.
Neither result authorizes E.8.1 closure because independent PDF review exposed
the presentation defects below.

## Independently observed PDF defects

The fresh PDF was rendered and inspected outside ROOT. The defects are visible
in both PDFium and Poppler where applicable.

### A. Context page text clipping

Page 37 (`full_background.e8.context`) contains multiple lines that run past the
right page boundary. The clipped content includes:

- the presentation-only statement;
- the frozen-input SHA line;
- the combined artifact/validation fingerprint line;
- the L/B/A semantic definitions;
- the B-to-A interpretation sentence;
- the common-scale statement.

The validation fingerprint is therefore not fully visible on the page.

### B. ROOT/UTF-8 mojibake in the context page

The source uses Unicode em dashes in the E.8 context title and L/B/A definitions.
ROOT 6.24/08 renders these as literal mojibake such as:

`â€”`

This is visible in both independent PDF renderers. The E.8 presentation must
use renderer-safe ASCII punctuation.

### C. Handoff page text clipping

Page 47 (`full_background.e8.handoff`) clips the long D.10-omission sentence at
the right page boundary.

### D. Overlay-page top/header clipping

Pages 38, 41, and 44
(`full_background.e8.persisted_overlays` for t1/t2/t3) place the 3x9 plot grid
against the top page boundary. The first-row plot titles/header area are visibly
clipped. The explicit L/B/A legend is currently drawn on the mother canvas at
NDC y=0.945-0.995 after `canvas.Divide(3, 9)`, with no reserved header region;
the farm PDF does not present that key cleanly.

The sparse canonical child rows themselves are correctly retained rather than
suppressed. Do not remove empty/unavailable rows to improve appearance.

## Status consequences

For durable memory after this task is implemented and independently source
reviewed:

- F.6.2 remains `CLOSED / RUNTIME VALIDATED`.
- F.6.2.Fix.5 remains `CLOSED / RUNTIME VALIDATED`.
- E.8.1 remains `ACTIVE`.
- E.8.1.Fix.1 may be recorded as `CLOSED / RUNTIME VALIDATED` for the targeted
  frozen-parent-setting reader defect: the fresh farm PDF contains the available
  E.8 section rather than the prior `parent_setting_identity_invalid` unavailable
  page.
- E.8.1.Fix.2 may be recorded as `CLOSED / RUNTIME VALIDATED` for the profile
  provenance re-pin: the fresh bundle is complete, fail-closed provenance
  checks pass, and no unexpected committed files were admitted.
- E.8.1.Fix.3 is `ACTIVE` until ChatGPT independently reviews its actual diff.
- E.8.1.Debug.1.Fix.1 remains `SOURCE REVIEWED`; this bundle does not by itself
  establish complete launcher provenance or prove the intentional full-high
  epsilon skip.
- F.6.3 remains `BLOCKED` pending E.8.1.

Do not create an E.8.1 closure record.

## Scientific ownership

This task owns only final E.8 visual layout/rendering.

Frozen scientific content includes:

- accepted F.6.2 JSON;
- accepted artifact and validation fingerprints;
- parent/child inventory;
- persisted one-dimensional distributions;
- persisted joint maps;
- persisted DeltaH, kappa, DeltaP^K, f_refine^K, support and bootstrap content;
- L/B/A population definitions;
- kaon missing-mass window;
- no renormalization/recalculation policy.

Method A remains non-production.
Method B remains diagnostic only.
No production weight, yield, correction, or cross section may change.

## Allowed implementation files

Only:

- `src/cuts/full_background_subtraction_plots.py`
- `testing/test_full_background_subtraction_plots.py`

## Allowed durable-memory files

Update only as warranted:

- `docs/memory/CURRENT.md`
- `docs/memory/phases/e8-1-fix1-frozen-f6-2-setting-reader.md`
- `docs/memory/phases/e8-1-fix2-bundle-profile-repin.md`
- `docs/memory/phases/e8-1-fix3-procedure-pdf-layout.md`
- `docs/memory/evidence/e8-1-fix2-left-lowe-layout-blocker.md`
- `docs/memory/manifest.json`

This task contract is also an intended tracked file:

- `docs/memory/phases/e8-1-fix3-procedure-pdf-layout-task-contract.md`

Regenerate `docs/memory/manifest.json` whenever versioned memory changes.

Do not change `MEMORY.md`, `USER.md`, `CODEX.md`, `TOOLS.md`,
`CURRENT_HANDOFF.md`, roadmap, decisions, or unrelated evidence unless a
concrete blocker proves their owned durable knowledge changed. Stop instead of
expanding scope automatically.

## Frozen files and interfaces

Do not modify:

- accepted F.6.2 farm artifacts;
- `run_Prod_Analysis.sh`;
- `src/main.py`;
- `src/cuts/rand_sub.py`;
- any pion/proton subtraction or HGCer scientific calculation source;
- binning, SIMC, yield, efficiency, cross-section, L/T or uncertainty source;
- E.8 reader authority constants;
- `_e8_validate_artifact`;
- `_e8_validate_parent_setting`;
- canonical parent/child inventory checks;
- E.8.1 validation profile;
- E.8.1 profile test;
- generic collector;
- bundle wrapper;
- canonical setting list;
- artifact declarations.

Do not alter page IDs or the 36 + 11 ordinary page inventory.

## Required context/handoff repair

Repair only visual text layout.

### Renderer-safe punctuation

Replace E.8 presentation-only Unicode em dashes with ASCII-safe punctuation in
the context page. At minimum:

- title:
  `E.8 - frozen F.6.2 acceptance-refinement atlas`
- L/B/A definitions use ` - ` rather than Unicode em dash.

Do not globally rewrite unrelated source strings.

### Explicit line layout

`TPaveText.AddText()` does not provide reliable automatic wrapping for these
ROOT pages. Do not rely on implicit wrapping.

Refactor the E.8 text-page interface narrowly so context and handoff pages can
supply renderer-safe, intentionally wrapped display lines and an appropriate
text size.

Requirements:

- no context or handoff display line may depend on clipping at the page edge;
- show the full 64-character frozen-input SHA-256;
- show the full artifact fingerprint;
- show the full validation fingerprint;
- preserve all current scientific/provenance statements;
- preserve the L/B/A definitions exactly in meaning;
- preserve the explicit F.6.3 handoff;
- preserve the statement that no Method-A production correction/promotion
  occurs in E.8;
- use ASCII punctuation where ROOT text rendering is known to be fragile.

A small optional `size` parameter on `_e8_text_page`, or a similarly narrow
presentation-only mechanism, is allowed.

Prefer explicit semantic line breaks over arbitrary breaking inside hash
strings. A label may occupy one line and its full hash/fingerprint the next.

Do not remove provenance text merely to make it fit.

## Required overlay-page layout repair

Keep the complete 3 variables x 9 canonical phi-child overlay inventory.

Do not suppress sparse or empty child rows.

Provide a real header/legend region that is separate from the 3x9 plot grid.

Preferred implementation:

- create a small header `TPad` at the top of the overlay canvas;
- create a grid `TPad` below it;
- divide only the grid pad into `3 x 9`;
- render all 27 overlay tiles into that grid pad;
- render the explicit three-column L/B/A legend into the header pad.

An equivalent narrow layout is acceptable only if it demonstrably reserves
non-overlapping vertical space and keeps the top plot row away from the page
boundary.

The legend must remain:

- L: upstream `0 < NPE <= 2` diagnostic reference;
- B: physical pion control, `NPE > 2`, baseline `w0`;
- A: same B population, `w0*C`.

Preserve the existing colors and sparse-display proxy behavior.

Do not reorder children, populations, variables, t parents, or page IDs.

## Map pages

The farm-rendered delta-acceptance and missing-mass-acceptance map pages are
not the blocker in this fix.

Do not redesign them merely for style.

Only touch shared helper behavior if required to keep the context/handoff and
overlay fixes internally consistent, and prove no map semantics or inventory
changes.

## Required regression tests

Extend `testing/test_full_background_subtraction_plots.py`.

### Text-page tests

Use the existing fake-ROOT text capture to prove:

- context title and L/B/A key contain no Unicode em dash and no mojibake literal;
- the complete frozen SHA, artifact fingerprint, and validation fingerprint are
  present in the captured text;
- context and handoff semantics remain explicit;
- display strings are split into bounded presentation lines rather than one
  known-overlong line;
- F.6.3 handoff remains present;
- no production-promotion language is weakened.

Do not write a test that simply searches for the exact bugfix source text.

### Overlay-layout tests

Extend the fake ROOT/canvas/pad plumbing as needed to prove:

- an explicit header region exists;
- the 3x9 grid is below that header;
- all 27 tile positions are still rendered;
- the L/B/A legend is drawn into the reserved header rather than over a data
  tile or against the top page edge;
- sparse legend proxies remain supported;
- child order and population colors are unchanged.

### Existing reader/provenance tests

All existing E.8 reader, exact-authority, malformed-payload, page-order,
sparse-state, shared-scale, and runtime-ownership tests must remain intact.

## Optional ROOT-aware regression

If the local environment has PyROOT, add or extend a focused real-ROOT render
test only if it can deterministically assert the repaired layout without
introducing fragile pixel-perfect thresholds.

Do not make local source acceptance depend on a renderer that is unavailable
off-farm.

Farm PDF inspection remains the final visual gate.

## Deterministic local validation

Run at minimum:

```text
python -m py_compile \
  src/cuts/full_background_subtraction_plots.py \
  testing/test_full_background_subtraction_plots.py
```

Run focused E.8 tests, then:

```text
python -m unittest testing.test_full_background_subtraction_plots -v
```

Report exact skips.

Also run the normal memory manifest/integrity/health/bootstrap checks required
by current repository instructions and:

```text
git diff --check
```

These are local/source checks only.

## Required actual-diff audit

Before stopping:

1. confirm committed HEAD remains
   `9f7094b19f3cdb0887c95f1f972d98471f2ecfdc`;
2. show `git status --short`;
3. show `git diff --stat`;
4. show `git diff --name-only`;
5. confirm every changed path is allowlisted;
6. confirm no frozen source/profile/collector/wrapper changed;
7. inspect the exact E.8 runtime rendering path;
8. confirm page IDs/inventory and scientific payload consumption are unchanged.

Create one complete temporary root-level:

`kaonlt_review.diff`

containing the tracked diff plus complete no-index additions for every intended
new/untracked file. Do not stage merely to create it.

## Durable evidence record

Create:

`docs/memory/evidence/e8-1-fix2-left-lowe-layout-blocker.md`

Record only directly established facts from the supplied farm bundle:

- bundle identity and SHA-256;
- complete manifest/provenance PASS;
- frozen F.6.2 JSON identity;
- PDF and page-manifest SHA-256 values;
- 47-page inventory and empty renderer-failure list;
- Fix.1 reader success;
- Fix.2 provenance success;
- independent PDF review findings:
  context clipping, em-dash mojibake, handoff clipping, overlay top/header
  clipping;
- E.8.1 overall remains ACTIVE;
- the evidence is a blocker record, not E.8.1 acceptance.

Do not invent terminal-log facts that are absent from the bundle.

## NEXT after this source repair

After ChatGPT independently reviews the actual Fix.3 diff:

1. user commits/pushes Fix.3;
2. ChatGPT reviews the pushed source;
3. create a separate narrow E.8.1 profile re-pin to the new Fix.3 source commit;
4. user commits/pushes that profile re-pin;
5. ChatGPT reviews the pushed profile commit;
6. rerun the targeted `Q4p4W2p74 / Left / lowe` debug analysis;
7. collect a fresh Left/lowe bundle;
8. independently inspect the repaired PDF;
9. only after Left/lowe passes return to detailed `Left / highe`;
10. broaden only after the detailed gate passes.

Do not try to encode the future Fix.3 commit into the profile in this task.

## Forbidden shortcuts

Do not:

- modify or regenerate the accepted F.6.2 JSON;
- change persisted values to improve plots;
- remove empty/sparse children;
- suppress unavailable-state labels;
- change L/B/A population definitions;
- alter common map scaling;
- remove full provenance hashes from the context page;
- hide clipping by deleting explanatory text;
- change page IDs or page order;
- broaden the validation-profile allowlist;
- modify collector or wrapper behavior;
- alter production physics;
- commit, push, or run the farm.

## Acceptance criteria

Source-review candidate PASS requires:

- context page uses renderer-safe punctuation;
- all context provenance and semantics are retained in bounded display lines;
- handoff text is bounded and retains its complete semantics;
- overlay pages reserve a visible non-overlapping header/legend region;
- all 27 overlay tiles per t parent remain;
- all 9 E.8 parent pages plus context/handoff remain exactly 11 E.8 pages;
- sparse children remain visible;
- no scientific payload or production logic changes;
- deterministic local tests pass to local capability;
- memory/evidence accurately records the successful reader/provenance gates and
  the new presentation blocker;
- actual diff is narrow and allowlisted.

## Hard stop

Do not commit, push, alter the E.8.1 validation profile, or run Jefferson Lab
farm validation.

Return to ChatGPT with:

- implementation summary;
- exact changed files;
- exact local checks/tests;
- runtime/presentation path trace;
- complete root-level `kaonlt_review.diff`;
- proposed user Git handoff commands, not executed.

ChatGPT must inspect the actual diff before user commit/push.
