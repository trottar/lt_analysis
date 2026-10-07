# KaonLT E.8.2 Left/lowe PDF-legibility repair — task contract

## 1. Task class and authority

This is a **narrow presentation-only source repair** for the active E.8.2
Left/lowe baseline full-analysis audit.

It is triggered by a fresh farm run at reviewed source
`8ae17e18e480527cf470a93a7a8fb6f0747308fd` that passed the tracked E.8.2
owner, isolation, provenance, freshness, artifact, page-manifest and ZIP gates,
but failed independent rendered-PDF visual acceptance.

The owner/isolation repair itself is accepted for its narrow scope. This task
does **not** reopen that runtime-isolation repair and does not modify KaonLT
physics.

Follow repository-root `AGENTS.md` and the exact startup core before editing:

1. `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read only the task-relevant records, including:

- `docs/memory/MAINTENANCE.md`
- `docs/memory/CODEX.md`
- `docs/memory/COMMUNICATION.md`
- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/phases/e8-2-left-lowe-runtime-owner-output-isolation-repair.md`

Authority remains:

```text
current source/diff
-> newest applicable farm/runtime evidence
-> newest applicable validation/handoff artifact
-> tracked durable memory
-> older chat/history
```

---

## 2. Exact starting identity and start gate

Required repository:

```text
https://github.com/trottar/lt_analysis/tree/test
```

Required branch:

```text
test
```

Required starting HEAD and local `origin/test`:

```text
8ae17e18e480527cf470a93a7a8fb6f0747308fd
```

Commit:

```text
Repair E8.2 owner OUTPUT isolation
```

Before editing, establish:

```bash
git branch --show-current
git rev-parse HEAD
git rev-parse refs/remotes/origin/test
git status --short --untracked-files=all
git diff --check
```

The only expected new task file before implementation is:

```text
docs/memory/phases/e8-2-left-lowe-pdf-legibility-repair-task-contract.md
```

Temporary root-level `kaonlt_review*.diff` files are permitted only as review
artifacts.

If branch, HEAD, local `origin/test`, or unrelated workstation state differs,
**STOP as `BLOCKED`**. Do not reset, stash, clean, overwrite, or repair unrelated
user state.

No farm command is authorized in this implementation task.

---

## 3. Accepted returned farm evidence

The user supplied the fresh E.8.2 runtime package with stem:

```text
KaonLT_E8_2_Left_lowe_scientific_audit_20261006-204510
```

ChatGPT independently reviewed the ZIP, gate-status, run-summary, log, page
manifest, structured artifacts, and rendered E.8.2 PDF pages.

External artifact SHA-256 identities:

```text
ZIP:
60451f4cbf2848c972fc62152bf3c300933a3f6b3e002e9288a9a4d6e6228826

gate-status JSON:
17c772994ea73c58bd748c3db6fbd0f1f9eaa7e54990b88f73f3b70afabc011a

run-summary JSON:
43d89180773dcedcfa4c061fab9de9904c7069d3b0da3e100c807e86d3d02f9f

child log:
fbb786b6b82cd2868e4b143bb30241134ba124289e97bdb698db622bb9d80e48
```

Packaged procedure PDF:

```text
Left_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction.pdf
SHA-256:
a9c761ae2d7e046d4f730d7df6afc374faab96b812d548ddd1c1907906d56594
bytes:
1555866
pages:
74
```

Packaged E.8.2 page manifest:

```text
SHA-256:
e51fe03cff97ce00e05c850b031206e97a9cb9d15b0c71279618fda1a9577f02
bytes:
18203
```

`RUNTIME VERIFIED` by the supplied evidence and ChatGPT review:

- owner gate status is `success`;
- owner stage is `complete`;
- failure reason is null;
- farm-evaluated source is
  `8ae17e18e480527cf470a93a7a8fb6f0747308fd`;
- child analysis return code is 0;
- repaired disposable-worktree OUTPUT link passed;
- child OUTPATH resolved to the configured external KaonLT artifact root;
- ordinary farm checkout preservation passed before and after collection;
- installed `ltsep` preservation passed before and after collection;
- artifact verification passed;
- collection passed;
- ZIP verification passed;
- the bundle manifest has `complete=true` and no errors;
- requested runtime scope is only `Q4p4W2p74 / Left / lowe`;
- all 18 required E.8.2 page IDs are structurally present;
- no invalid/unavailable canonical phi children were declared;
- owner/isolation infrastructure is therefore
  `CLOSED / RUNTIME VALIDATED` for this narrow Left/lowe execution scope.

The returned run **does not** establish E.8.2 scientific/visual acceptance
because rendered-page inspection failed.

Codex does not need local copies of these farm artifacts for this source task
and must not claim to have independently rerun or revalidated them. Repository
memory must attribute the runtime facts to the supplied artifacts and ChatGPT
review.

---

## 4. First failed invariant: rendered-page legibility/completeness

ChatGPT visually inspected all required E.8.2 pages, which occupy PDF pages
47 through 64 in the 74-page procedure PDF.

The owner/page-manifest checker passed structural inventory but did not and
cannot establish visual legibility.

### 4.1 Subtraction pages — phi 1 and phi 2 clipped

Affected pages:

```text
47–50  t1
53–56  t2
59–62  t3
```

These are the random, dummy, proton, and baseline-pion E.8.2 pages.

The page manifest declares all nine canonical phi children, but the actual PDF
visibly clips the first two child rows above the physical page. PDF text
extraction begins at phi 3 on these pages.

Current source in:

```text
src/cuts/full_background_subtraction_plots.py
```

uses a tall canvas:

```python
canvas = ROOT.TCanvas(..., 1800, max(600, 360 * len(children)))
canvas.Divide(len(columns), len(children))
```

for nine child rows. This farm rendering demonstrates that the physical
multi-page PDF does not preserve the full tall-canvas region.

### 4.2 Final-MM pages — header overlaps top-row plots

Affected pages:

```text
51
57
63
```

The current final-MM renderer divides the full canvas into the 3x3 child grid,
then draws `_draw_page_header(...)` over the same canvas NDC region. The farm
PDF visibly places the E.8.2 heading / t-range on top of the first-row plots.

### 4.3 Stage-yield pages — long lines truncate on the right

Affected pages:

```text
52
58
64
```

The current stage-yield renderer creates one long text line per phi containing
all persisted stage names and values. The actual farm PDF truncates those lines
on the right. Later stage values are therefore not fully readable.

### 4.4 Acceptance consequence

This is the first failed invariant. Therefore:

```text
E.8.2 owner/isolation infrastructure:
  CLOSED / RUNTIME VALIDATED

E.8.2 baseline scientific audit:
  BLOCKED at rendered-page visual acceptance

scientific interpretation of stage magnitudes/prune impact:
  NOT YET AUTHORIZED
```

Do not numerically interpret the E.8.2 stage progression in this source task.

---

## 5. Scientific and production boundaries

This repair is presentation-only.

The following are frozen and **must not change**:

```text
src/cuts/calculate_yield.py
src/cuts/rand_sub.py
src/cuts/particle_subtraction.py
src/cuts/pion_component_subtraction.py
src/cuts/pion_component_fits.py
src/cuts/pion_component_shapes.py
src/cuts/proton_contamination_weights.py
src/main.py
run_Prod_Analysis.sh
set_SymLinks.sh

testing/run_e8_2_left_lowe_scientific_audit_gate.py
testing/test_run_e8_2_left_lowe_scientific_audit_gate.py
testing/collect_pion_hgcer_validation_bundle.py
testing/pion_hgcer_validation_bundle_profile_e8_2_left_lowe.json

testing/run_e8_4_fix5_canonical_five_plot_gate.py
```

Also freeze all unrelated `src/`, `testing/`, `tools/`, `farm_env/`, and
`background_samples/` paths.

Do not change:

- event traversal;
- prompt/random subtraction;
- dummy subtraction;
- slow-proton weights or application;
- `prune_hist`;
- baseline pion weights `w0`;
- pion templates/components;
- final baseline clean-kaon spectrum;
- final extracted yields;
- statistical/total uncertainties;
- cuts;
- binning;
- efficiencies;
- acceptance;
- SIMC;
- normalization;
- L/T separation;
- cross sections;
- active `no_empirical_residual` profile;
- Method-A/Method-B ownership;
- E.8.2 source/payload schemas or scientific values;
- E.8.2 required page IDs or semantic stage identities.

Method A remains detached/non-production and is not required by E.8.2.
Method B remains diagnostic/cross-check only and numerically excluded.

Canonical-five provenance/identity repair remains `DEFERRED`.

---

## 6. Exact allowed versioned paths

Only these seven paths may change:

```text
src/cuts/full_background_subtraction_plots.py
testing/test_e8_2_baseline_stage_audit.py
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/evidence/e8-2-left-lowe-runtime-owner-pass-visual-blocker-2026-10-06.md
docs/memory/phases/e8-2-left-lowe-pdf-legibility-repair-task-contract.md
docs/memory/phases/e8-2-left-lowe-pdf-legibility-repair.md
```

No other versioned path may change.

The final audit must prove all frozen scientific/production/owner/profile/
collector paths are byte-identical to starting HEAD.

---

## 7. Repair objective

Repair only the E.8.2 presentation geometry so that the same persisted,
already-authoritative E.8.2 payload is visibly complete and readable in the
farm-generated procedure PDF.

Preserve exactly:

```text
3 canonical t parents
x
(
  random page
  dummy page
  proton/prune page
  baseline-pion page
  final-MM page
  stage-yield page
)
=
18 required E.8.2 pages
```

Do not add continuation pages and do not alter the required E.8.2 page IDs.

The repaired presentation must make all nine canonical phi children visible on
every applicable page.

---

## 8. Required renderer behavior

Modify only the E.8.2 presentation helpers in:

```text
src/cuts/full_background_subtraction_plots.py
```

Do not redesign unrelated D/F/E presentation pages.

### 8.1 Use bounded page geometry

The E.8.2 farm failure shows that a very tall ROOT canvas cannot be treated as
a reliably printable physical page.

Use a bounded, page-shaped or square canvas and explicit content regions.

The existing E.8 persisted-overlay renderer in the same file already provides
an established pattern:

```text
canvas
-> dedicated header TPad
-> dedicated grid TPad
-> grid_pad.Divide(...)
```

Prefer reusing that local presentation pattern or a narrowly shared E.8.2
equivalent rather than inventing an unrelated rendering framework.

A header must occupy its own NDC region and must not overlap the content grid.

### 8.2 Subtraction pages

Preserve the current physics-oriented column semantics exactly.

Random:

```text
prompt pre-proton
random component removed
after random pre-proton
```

Dummy:

```text
after random pre-proton
dummy component removed
after dummy pre-proton
```

Proton:

```text
after dummy pre-proton
proton component removed
production after proton, pre-prune
after existing production prune_hist; pion input state
```

Pion:

```text
proton-cleaned pion input
exact baseline pion component B_pi^0
final baseline clean kaon K_0
```

For each page:

- retain one row per canonical phi child;
- retain one column per listed stage;
- retain all nine phi children in canonical order 1 through 9;
- retain existing Lambda-window markers;
- retain signed normalized-yield axes;
- retain phi range labeling;
- retain invalid/unavailable behavior;
- retain shared per-child y-range semantics;
- do not overlay/recombine/recalculate stages merely to solve layout;
- do not omit phi 1 or phi 2;
- do not create a second page for the same stage.

Use a bounded canvas/grid that physically fits all rows.

### 8.3 Final-MM pages

Preserve the current 3x3 final-MM child grid and all existing plotted
information:

- final baseline `MM_0`;
- all nine phi children;
- Lambda window;
- phi index/range;
- authoritative `Y0`;
- statistical uncertainty;
- existing total uncertainty.

Move the E.8.2 page heading and t-context into a dedicated non-overlapping
header region. The heading must never cover any plot, note, axis label, or
top-row panel.

### 8.4 Stage-yield pages

Preserve all six persisted diagnostic Lambda-window stage integrals and the
authoritative final yield information for every valid phi child.

The six stage quantities are:

```text
prompt_pre_proton
after_random_pre_proton
after_dummy_pre_proton
after_proton_pre_prune
after_proton_post_prune
after_pion_final
```

Also preserve:

```text
final_yield
statistical_error
total_error
```

Do not drop a value to make text fit.

Replace the current single unbounded line per phi with an explicitly bounded
presentation. A narrow acceptable implementation is:

- dedicated header region;
- compact, documented display labels for the six stage names;
- at most two or three bounded text rows per phi child;
- all values readable within the page;
- a visible legend/heading that defines any abbreviated labels;
- retain the statement that the stage-window integrals are diagnostic and are
  **not** final extracted yields.

If using a formatting helper, make it pure/deterministic so unit tests can
verify that every required value is represented exactly once.

Do not recompute integrals or yields in the renderer.

---

## 9. Preserve manifest and scientific semantics

The repair must not alter:

```text
E8_2_SOURCE_SCHEMA_VERSION
E8_2_PRESENTATION_SCHEMA_VERSION
```

unless a concrete source blocker proves impossible to avoid. If a schema change
appears necessary, **STOP as `BLOCKED`** instead of broadening scope.

Preserve the exact required E.8.2 page IDs:

```text
full_background.e8_2.random.t1
full_background.e8_2.random.t2
full_background.e8_2.random.t3

full_background.e8_2.dummy.t1
full_background.e8_2.dummy.t2
full_background.e8_2.dummy.t3

full_background.e8_2.proton.t1
full_background.e8_2.proton.t2
full_background.e8_2.proton.t3

full_background.e8_2.pion.t1
full_background.e8_2.pion.t2
full_background.e8_2.pion.t3

full_background.e8_2.final_mm.t1
full_background.e8_2.final_mm.t2
full_background.e8_2.final_mm.t3

full_background.e8_2.stage_yields.t1
full_background.e8_2.stage_yields.t2
full_background.e8_2.stage_yields.t3
```

Each record must continue to represent canonical phi inventory
`0..8` with no invalid child for this accepted Left/lowe evidence.

The existing owner remains the farm authority; do not weaken or bypass it.

---

## 10. Deterministic local regression coverage

Extend only:

```text
testing/test_e8_2_baseline_stage_audit.py
```

Use lightweight/fake ROOT objects as needed. Do not add ROOT/PyROOT as a local
dependency.

At minimum add deterministic coverage for the exact visual-regression
mechanisms.

### 10.1 Subtraction geometry

Prove for a nine-child fixture that each subtraction renderer:

- creates a bounded canvas rather than height scaling with nine rows;
- creates a distinct header region and content grid;
- divides the content grid into exactly:
  - `3 x 9` for random;
  - `3 x 9` for dummy;
  - `4 x 9` for proton;
  - `3 x 9` for pion;
- visits every canonical phi child 1 through 9;
- draws every required stage column for every valid child;
- appends exactly one existing manifest page record.

The test should fail if the implementation regresses to direct
`canvas.Divide(columns, 9)` on an unbounded/tall canvas with the header sharing
the plot region.

### 10.2 Header separation

For subtraction and final-MM renderers, prove the header and content regions
have non-overlapping NDC bounds.

For final MM, prove the content grid remains 3x3 for nine children.

Do not treat canvas pixel size alone as sufficient; assert the separate
header/content region behavior.

### 10.3 Stage-yield text completeness

Create a valid nine-child fixture with distinguishable values for all six stage
integrals plus Y0/stat/total.

Prove the formatting/rendering path represents every required field for every
phi child.

If display abbreviations are introduced, assert the legend maps every
abbreviation unambiguously to the persisted stage name.

Assert the implementation uses bounded lines/rows and does not concatenate all
six long persisted stage names into the previous single unbounded phi line.

### 10.4 Scientific immutability

Retain or extend tests proving:

- source histograms are cloned only for display;
- no source histogram content/error is modified;
- Lambda-window values are consumed, not recalculated;
- final Y0/stat/total are consumed from payload;
- stage-window integral values are consumed from payload;
- manifest page IDs/semantic-stage identities are unchanged;
- proton page still visibly identifies the `prune_hist` boundary between
  pre-prune and post-prune columns;
- renderer failure/finalization recovery behavior remains unchanged.

### 10.5 Existing suite

All existing E.8.2 tests must continue to pass.

---

## 11. Local checks

Do not run ROOT, PyROOT, the production analysis, or the farm.

Discover the local Python interpreter according to tracked
`docs/memory/TOOLS.md`.

Run at minimum:

```text
py_compile:
  src/cuts/full_background_subtraction_plots.py
  testing/test_e8_2_baseline_stage_audit.py

testing.test_e8_2_baseline_stage_audit

testing.test_run_e8_2_left_lowe_scientific_audit_gate

directly relevant existing full-background-subtraction presentation tests
identified from source

git diff --check
```

If an existing deterministic fake-ROOT presentation suite directly covers the
same renderer helpers, run it rather than inventing a redundant broad suite.

Report exact commands, interpreter/OS, counts, skips, and results.

Local fake-ROOT geometry tests are `SOURCE VERIFIED` evidence only; actual PDF
legibility remains farm-only.

---

## 12. Repository evidence record

Create:

```text
docs/memory/evidence/e8-2-left-lowe-runtime-owner-pass-visual-blocker-2026-10-06.md
```

Record only facts accepted by the supplied farm evidence and ChatGPT review.

It must state:

```text
farm source:
8ae17e18e480527cf470a93a7a8fb6f0747308fd

scope:
Q4p4W2p74 / Left / lowe

owner/isolation:
CLOSED / RUNTIME VALIDATED

E.8.2 audit:
BLOCKED at rendered-page visual acceptance
```

Include the four external companion SHA-256 values and the packaged PDF/page
manifest identities listed in this contract.

Record the three visual blocker classes and exact page groups:

```text
subtraction clipping:
47–50, 53–56, 59–62
phi 1 and phi 2 not physically visible

final-MM header overlap:
51, 57, 63

stage-yield right-edge truncation:
52, 58, 64
```

State that structured/page-manifest validation passed but does not substitute
for rendered-page inspection.

Do not record scientific stage-magnitude conclusions.

---

## 13. Phase repair record and CURRENT

Create:

```text
docs/memory/phases/e8-2-left-lowe-pdf-legibility-repair.md
```

Record:

- starting source identity;
- accepted farm evidence identity;
- source-verified presentation cause;
- exact presentation-only repair;
- changed paths;
- deterministic local checks;
- no farm execution;
- runtime visual validation still pending.

Update `docs/memory/CURRENT.md` to remove stale pre-farm wording and reflect:

```text
E.8 remains ACTIVE.

The E.8.2 Left/lowe owner/isolation repair is
CLOSED / RUNTIME VALIDATED for the narrow owner/isolation/provenance scope
at farm source 8ae17e18...

The returned E.8.2 run passed owner/provenance/artifact/page-manifest/ZIP gates
but is BLOCKED at rendered-page visual acceptance.

No E.8.2 scientific magnitude/prune interpretation has been accepted yet.

A narrow E.8.2 presentation-only PDF-legibility repair is the active source
task.

Canonical-five provenance repair remains DEFERRED.

Final E.8 and F.6.4 remain BLOCKED.
```

After successful local implementation but before fresh farm validation, mark the
presentation repair:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

The sole ordinary NEXT after successful implementation must be equivalent to:

```text
NEXT — after ChatGPT actual-diff review, user commit/push and pushed-state
synchronization/farm-readiness review, rerun exactly one isolated Q4p4W2p74
Left/lowe E.8.2 owner gate; inspect the fresh E.8.2 rendered pages before any
scientific stage-magnitude/prune interpretation.
```

Do not make commit/push itself the scientific NEXT.

`CURRENT.md` was already at approximately the 8 KiB soft threshold before this
checkpoint. Consolidate stale wording rather than adding a long new block.
Prefer keeping it below 8192 bytes. Do not perform a broad memory rewrite.

Regenerate/check `docs/memory/manifest.json` according to tracked procedure.

---

## 14. Farm-validation boundary

No farm run occurs during Codex implementation.

After:

```text
implementation
-> deterministic local checks
-> ChatGPT actual-diff review
-> user commit/push
-> ChatGPT pushed-state synchronization
-> farm-readiness review
```

the user will run exactly one fresh isolated E.8.2 Left/lowe owner gate using
the existing reviewed owner.

The fresh review must again validate owner/provenance/artifact/page-manifest
gates and then visually inspect all 18 E.8.2 pages.

Visual acceptance requires:

- all nine phi children physically visible on every subtraction page;
- no top-row clipping;
- no header/content overlap;
- every final-MM plot and Y0/stat/total note readable;
- every stage-yield stage value readable;
- no right-edge truncation;
- no new rendering failures;
- page identities and scientific semantics unchanged.

Only after that visual gate passes may the E.8.2 stage magnitudes and prune
impact be scientifically interpreted.

---

## 15. Final integrity/health gate

At the end:

1. regenerate and check the memory manifest;
2. run ordinary memory health;
3. run bootstrap JSON if required by tracked final procedure;
4. run `git diff --check`;
5. confirm only the seven allowed versioned paths changed;
6. confirm frozen scientific/production/owner/profile/collector files are
   byte-identical to starting HEAD;
7. produce one complete root-level `kaonlt_review.diff` containing:
   - tracked diffs;
   - complete no-index additions for every new/untracked contract/evidence/phase
     file.

Do not stage files merely to create the review bundle.

Report:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
warning classification: <none or blocking/nonblocking>
manifest check: PASS | FAIL
bootstrap: PASS | FAIL
git diff --check: PASS | FAIL
```

---

## 16. Hard stops

Stop as `BLOCKED` rather than broadening scope if any of the following appears
necessary:

- changing `calculate_yield.py`;
- changing subtraction/proton/pion production logic;
- changing the E.8.2 source or presentation schema;
- changing E.8.2 page IDs;
- adding continuation pages;
- changing the E.8.2 owner/profile/collector;
- changing canonical-five owner/provenance machinery;
- changing final yields or uncertainties;
- activating Method A or Method B numerically;
- changing `no_empirical_residual`;
- changing SIMC or normalization;
- running ROOT/PyROOT locally;
- running the farm;
- interpreting scientific stage magnitudes before the fresh visual gate passes.

---

## 17. Acceptance criteria for source review

The candidate is ready for ChatGPT actual-diff review only when all are true:

1. exactly the seven allowed versioned paths changed;
2. all scientific/production files are byte-identical;
3. owner/profile/collector files are byte-identical;
4. E.8.2 source/presentation schemas are unchanged;
5. all 18 required E.8.2 page IDs are unchanged;
6. subtraction pages use bounded geometry with separate header/content regions;
7. all nine phi children and every current stage column remain represented;
8. final-MM pages retain 3x3 child coverage with a non-overlapping header;
9. stage-yield pages preserve all six diagnostic integrals plus Y0/stat/total
   in bounded readable rows;
10. no renderer recomputes scientific quantities;
11. deterministic layout/text/immutability tests pass;
12. existing E.8.2 owner tests pass unchanged;
13. memory evidence/status accurately distinguish owner runtime PASS from
    E.8.2 visual BLOCKED;
14. manifest/health/bootstrap/diff checks pass;
15. no farm command was run.

Status after successful local implementation:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

This applies only to the PDF-legibility repair. The owner/isolation remains
`CLOSED / RUNTIME VALIDATED`; E.8.2 scientific audit remains `BLOCKED` until
fresh rendered-page acceptance.

---

## 18. Codex stop point

Codex must stop after:

- narrow renderer implementation;
- deterministic local tests;
- evidence/phase/CURRENT memory update;
- manifest/health/bootstrap checks;
- complete `kaonlt_review.diff` creation.

Codex must **not**:

- stage merely for review;
- commit;
- push;
- update refs;
- run ROOT/PyROOT;
- run production analysis;
- run the farm;
- interpret E.8.2 scientific stage magnitudes.

Return:

- starting/ending branch;
- starting/ending HEAD and local `origin/test`;
- exact changed paths;
- concise renderer repair summary;
- exact local checks and results;
- byte-preservation result for frozen paths;
- memory-health report;
- `kaonlt_review.diff` size/path;
- explicit confirmation that no farm/ROOT/PyROOT command was run.
