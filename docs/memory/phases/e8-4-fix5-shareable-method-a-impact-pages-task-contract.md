# KaonLT E.8.4 Fix.5 — shareable Method-A impact pages with per-(t,phi) SIMC comparison and tracked Left/lowe farm gate

## 1. Objective

Implement the next **presentation-only E.8 extension** needed to decide and communicate whether the detached Method-A reweighting behaves as expected in the real full-analysis chain.

The immediate runtime scope remains:

```text
Q4p4W2p74 / Left / lowe
```

The new procedure-PDF output must make the following visually explicit for **all 27 canonical `(t,phi)` cells**:

1. baseline versus Method-A final clean-kaon missing mass;
2. final Method-A clean-kaon missing mass versus the **existing authoritative SIMC missing-mass prediction**;
3. baseline and Method-A extracted yield versus phi for each t parent;
4. absolute yield change `DeltaY = YA - Y0` versus phi for each t parent;
5. fractional yield change `DeltaY / Y0` versus phi for each t parent, only where defined;
6. parent-t preservation/closure.

The purpose is to make the E.8 result less diagnostic-looking and more directly shareable without changing any physics.

This task also performs the **milestone memory alignment warranted by the already-supplied Left/lowe runtime evidence** and creates a tracked one-command farm owner that runs the existing Left/lowe debug analysis, verifies the required new E.8 pages, and packages the fresh evidence. The user must not again be left with a completed multi-hour run and no tracked packaging path.

No production promotion is allowed. Method A remains detached. Method B remains diagnostic-only.

---

## 2. Exact starting state

Required branch:

```text
test
```

Required starting HEAD:

```text
da38444e7aa60efd62d6638780776344daf40276
```

Commit subject at that HEAD:

```text
F6.3: adopt detached current-baseline candidate lineage
```

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short --untracked-files=all
git log --oneline -1
```

Requirements:

- branch must be `test`;
- committed HEAD must be exactly `da38444e7aa60efd62d6638780776344daf40276`;
- the only intended pre-existing versionable change may be this task-contract file after the user copies it into the repository;
- root `AGENTS.md`, `.codex/`, and temporary review bundles remain local-only/untracked as applicable;
- unrelated tracked changes are a blocker;
- do not reset, stash, clean, commit, push, or run the farm.

If these conditions are not met, STOP.

---

## 3. Mandatory startup read and memory alignment audit

Read in order:

1. root `AGENTS.md`;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read the task-relevant canonical records:

- `docs/memory/CODEX.md`
- `docs/memory/COMMUNICATION.md`
- `docs/memory/MAINTENANCE.md`
- `docs/memory/LEARNINGS.md`
- `docs/memory/TOOLS.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/decisions/farm-validation-bundle-procedure.md`
- `docs/memory/phases/f6-3-current-baseline-candidate-lineage-adoption.md`
- `docs/memory/phases/f6-3-current-baseline-candidate-lineage-adoption-task-contract.md`

Then read the relevant source/tests:

- `src/cuts/full_background_subtraction_plots.py`
- `src/cuts/rand_sub.py`
- `src/cuts/pion_hgcer_method_a_parallel_full_procedure.py`
- `testing/test_full_background_subtraction_plots.py`
- `testing/test_e8_4_production_impact_audit.py`
- `testing/test_f6_3_parallel_full_procedure_method_a.py`
- `testing/test_run_prod_analysis_debug_left_low.py`
- `testing/collect_pion_hgcer_validation_bundle.py`
- `testing/package_pion_hgcer_validation_bundle.tcsh`
- current applicable bundle profiles/tests.

Do not reopen unrelated historical phase records.

### Known memory issue that must be repaired

`docs/memory/TOOLS.md` is stale relative to the memory-bracketed workflow: it still states that `--fail-on-warning` is the universal task-final health gate, while `CODEX.md`, `MAINTENANCE.md`, `USER.md`, `AGENTS.md`, and the contract template correctly make the non-strict health command the ordinary gate and reserve strict-warning mode for explicit zero-warning/milestone cases.

Repair `TOOLS.md` in this task. Do not create a separate memory-only cycle.

---

## 4. Already-supplied runtime evidence to record

The user supplied the fresh post-F.6.3 Left/lowe runtime package:

```text
KaonLT_F6_3_Left_lowe_runtime_evidence_Q4p4W2p74_20261001-222142.zip
```

Independent ChatGPT evidence review established:

```text
bundle SHA-256:
200fda66fe410274df1c9a8252b9e87114d8fd10a520fb7d71691ec6b3772874

farm/source HEAD:
da38444e7aa60efd62d6638780776344daf40276

kinematic:
Q4p4W2p74

setting:
Left / lowe

candidate F.3 SHA-256:
eeb480d9b7ddffab1e97c7f9f27c3099d0780c2bbb4d42a0ee68c077a5752b3d

candidate F.4 SHA-256:
1d545924eba89c7f9ffa28028e307aca9b434a89beec06863cf2893887b6b902
```

The supplied procedure-PDF manifest showed:

```text
87 pages
renderer_failures = []
```

and contained the complete existing E.8.4 block:

```text
full_background.e8_4.authority
full_background.e8_4.pion_consequence t1/t2/t3
full_background.e8_4.final_mm t1/t2/t3
full_background.e8_4.signed_difference t1/t2/t3
full_background.e8_4.yield_impact t1/t2/t3
full_background.e8_4.setting_summary
```

Direct evidence showed:

- F.6.3 branch role is `parallel_nonproduction_method_a_full_analysis`;
- live-cache parity passed;
- Method B is numerically excluded;
- production promotion is false;
- baseline production is not mutated;
- Method-A reweighting changes real `(t,phi)` yields;
- the largest observed fractional changes are concentrated in the lowest t parent;
- parent-t preservation closes to floating-point precision.

The exact yield-impact evidence includes, among others:

```text
t1 phi [-180,-140): DeltaY/Y0 about -25.8%
t1 phi [140,180):   DeltaY/Y0 about -45.1%

t2 phi [-180,-140): about -4.55%
t2 phi [100,140):   about +8.65%

t3 phi [-180,-140): about -1.83%
t3 phi [100,140):   about +2.84%
```

Do not encode these percentages as physics constants in source. They belong only in the evidence record.

### Provenance caveat to preserve

The packaging-time worktree snapshot contained:

```text
M src/models/xmodel_kaon_pl.f
?? src/kaon/functions/Q4p4W2p74.model
```

The package was created **after** the run, so that snapshot cannot prove the pre-run checkout was globally clean. Record this as a provenance caveat. Do not silently erase it, and do not reinterpret it as a demonstrated scientific failure because the runtime candidate identities, live-cache parity, and baseline outputs were independently checked.

### Status established by that evidence

Record narrowly:

```text
F.4.Refresh.2 detached candidate materialization:
CLOSED / RUNTIME VALIDATED

F.6.3 current-baseline candidate-lineage adoption,
Q4p4W2p74 / Left / lowe runtime gate:
CLOSED / RUNTIME VALIDATED

E.8.4 Q4p4W2p74 / Left / lowe production-impact presentation/runtime gate:
CLOSED / RUNTIME VALIDATED
```

This does **not** close canonical-five E.8, does not close F.6.4, and does not promote Method A.

---

## 5. Scientific and production ownership — frozen

Preserve all of the following exactly:

- random subtraction;
- dummy subtraction;
- slow-proton subtraction;
- pion-component models/fits/windows/amplitudes;
- active `no_empirical_residual` profile;
- baseline pion weight `w0`;
- Method-A correction factors `C`;
- F.4 parent-preserving mathematics;
- F.5 propagation mathematics;
- F.6.3 `w0 -> w0*C` mathematics;
- live-cache parity;
- canonical t/phi binning;
- existing yield extraction;
- Lambda integration window;
- data and SIMC normalization;
- SIMC event weights and normalization factors;
- efficiencies;
- acceptance;
- L/T separation;
- cross-section formulas;
- Method-B diagnostic-only status.

No fit, correction, normalization, yield, or SIMC scaling may be changed to improve visual agreement.

E.8 remains a **consumer/presentation layer**.

---

## 6. Source/runtime audit required before implementation

Codex must trace the real current path for the existing per-`(t,phi)` SIMC missing-mass objects before editing.

Trace:

```text
SIMC producer / already-normalized in-memory histogram
-> existing data-vs-SIMC consumer
-> E.8 render-state capture or F.6.3/E.8.4 handoff
-> E.8.4 payload
-> procedure-PDF renderer
-> page manifest
```

Establish:

1. which in-memory SIMC missing-mass histogram(s) are already used for the current data-vs-SIMC path;
2. whether there is already one authoritative SIMC MM object per canonical `(t,phi)` cell;
3. exactly which current normalization is applied before those objects are consumed;
4. whether E.8 can clone those existing objects directly at render-state capture time.

### Hard boundary

The E.8 renderer must **not**:

- reopen SIMC files;
- refill SIMC from events;
- rederive normalization;
- shape-normalize data or SIMC;
- introduce a new comparison normalization;
- reuse a setting-wide SIMC shape as if it were a canonical child unless the existing production path already does so and the identity is explicitly correct.

If authoritative per-cell SIMC MM is not available on the current runtime path, STOP and report the exact blocker rather than inventing a new SIMC calculation.

---

## 7. Required E.8 presentation additions

Keep every currently accepted E.8/E.8.2/E.8.3/E.8.4 page intact unless a small layout-only adjustment is required to avoid overlap in the newly appended pages.

Append the following **shareable E.8.4 pages**.

### 7.1 Final Method-A clean-kaon MM versus SIMC — all 27 cells

For each of the three canonical t parents, add one 3x3 page covering all nine phi children.

Required page IDs:

```text
full_background.e8_4.method_a_vs_simc.t1
full_background.e8_4.method_a_vs_simc.t2
full_background.e8_4.method_a_vs_simc.t3
```

Each populated pad must show:

- final Method-A clean-kaon `MM_A(t,phi)`;
- the authoritative already-normalized SIMC missing-mass histogram for the same canonical `(t,phi)` cell;
- identical MM axis/binning or an explicit fail-closed geometry check;
- the same Lambda signal/integration window used by the authoritative yield extraction;
- t index/range;
- phi index/range;
- a concise legend identifying `Method A data` and `SIMC`;
- existing stored Method-A final yield `YA` where available;
- explicit `EMPTY`/unavailable status for cells with no valid data, rather than silently dropping them.

Do not rescale one histogram to the other for display.

### 7.2 Baseline / Method A / SIMC direct comparison — all 27 cells

To make “did Method A move the data in a sensible direction?” visually answerable, add one 3x3 page per t parent overlaying:

- baseline final clean-kaon `MM_0`;
- Method-A final clean-kaon `MM_A`;
- the same authoritative SIMC MM object.

Required page IDs:

```text
full_background.e8_4.baseline_method_a_simc.t1
full_background.e8_4.baseline_method_a_simc.t2
full_background.e8_4.baseline_method_a_simc.t3
```

This is a display-only comparison. No goodness-of-fit score, fit, or renormalization is introduced in this task.

### 7.3 Shareable yield-impact summaries — all phi bins

For each t parent add one clean page with all nine canonical phi bins containing three panels:

1. `Y0` and `YA` versus phi;
2. `DeltaY = YA - Y0` versus phi;
3. `DeltaY / Y0` versus phi where defined.

Required page IDs:

```text
full_background.e8_4.yield_summary.t1
full_background.e8_4.yield_summary.t2
full_background.e8_4.yield_summary.t3
```

Requirements:

- consume stored authoritative `Y0`, `YA`, `DeltaY`, and defined fractional values from the E.8.4/F.6.3 payload;
- do not reintegrate MM histograms in the renderer;
- retain every canonical phi bin;
- undefined fractions must be visibly marked undefined, not set to zero;
- preserve sign;
- use existing stored uncertainty information when already available; do not invent a new uncertainty for `DeltaY` in this task.

### 7.4 Parent-t closure page

Add one clean setting-level page showing the stored parent-preservation closure for all three t parents:

```text
full_background.e8_4.parent_closure
```

Show:

- baseline signed parent sum;
- Method-A adjusted signed parent sum;
- stored closure residual / relative residual.

This page is a **sanity check**. Label it clearly so equality is understood as the required parent-preservation constraint, not the size of the reweighting effect.

### 7.5 Existing direct reweighting pages remain

Retain the existing E.8.3 and E.8.4 direct reweighting evidence:

- `B_pi^0` versus `B_pi^A`;
- pion redistribution per canonical child;
- `MM_0` versus `MM_A`;
- signed differences;
- original yield-impact pages.

Do not replace those with the new SIMC pages.

---

## 8. Presentation quality

These are E.8 analysis/presentation pages, not engineering diagnostic pages.

For the new pages:

- titles and legends must not overlap plotted content;
- use readable axis labels and line widths;
- keep a consistent legend convention across all t parents;
- use the same x-range and Lambda-window markers for directly compared MM pages;
- include phi range in each pad;
- include t range in the page title;
- no provenance dumps in plot titles;
- no giant technical text blocks on the shareable plot pages;
- detailed provenance remains on the authority/summary pages and manifest.

Do not remove provenance from the procedure PDF entirely; simply keep technical provenance separate from the shareable figures.

---

## 9. Allowed source files

Scientific/presentation runtime:

```text
src/cuts/full_background_subtraction_plots.py
src/cuts/rand_sub.py
```

`rand_sub.py` may change **only** if needed to hand already-existing authoritative SIMC per-cell MM objects into the E.8 render-state/payload. No calculation or normalization change is permitted.

Focused tests:

```text
testing/test_full_background_subtraction_plots.py
testing/test_e8_4_production_impact_audit.py
testing/test_run_prod_analysis_debug_left_low.py
```

Farm validation/profile ownership:

```text
testing/pion_hgcer_validation_bundle_profile_e8_4_fix5.json
testing/test_pion_hgcer_validation_bundle_profile_e8_4_fix5.py
testing/run_e8_4_fix5_left_lowe_plot_gate.py
testing/test_run_e8_4_fix5_left_lowe_plot_gate.py
```

If the generic collector needs no change, keep it frozen. If Codex believes the generic collector must change, STOP and report why rather than expanding scope automatically.

Memory:

```text
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/USER.md
docs/memory/LEARNINGS.md
docs/memory/TOOLS.md
docs/memory/roadmap/STATUS.md
docs/memory/decisions/e8-full-analysis-procedure-roadmap.md
docs/memory/manifest.json
docs/memory/phases/f6-3-current-baseline-candidate-lineage-adoption.md
```

Create exactly:

```text
docs/memory/evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md
docs/memory/phases/e8-4-fix5-shareable-method-a-impact-pages.md
docs/memory/phases/e8-4-fix5-shareable-method-a-impact-pages-task-contract.md
```

No other memory records unless an owned hard inconsistency makes one necessary. If so, STOP before expanding.

---

## 10. Frozen files

Do not modify:

```text
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
src/cuts/pion_hgcer_method_a_parent_preserving_correction.py
src/cuts/pion_hgcer_method_a_tphi_propagation.py
src/cuts/pion_component_subtraction.py
src/cuts/pion_component_fits.py
src/cuts/particle_subtraction.py
src/cuts/calculate_yield.py
binning/calculate_yield.py
run_Prod_Analysis.sh
main.py
```

Also freeze all F.1-F.6.2 accepted authority artifacts/identities and Method-B source.

No production physics change is allowed.

---

## 11. Validation profile and one-command farm owner

The new Left/lowe farm validation must be fully tracked **before** the user starts the multi-hour run.

### 11.1 Validation profile

Create:

```text
testing/pion_hgcer_validation_bundle_profile_e8_4_fix5.json
```

Use the existing generic collector schema. Scope it to:

```text
Q4p4W2p74 / Left / lowe
```

The profile must require at minimum:

- fresh `Left_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction.pdf`;
- fresh corresponding page-manifest JSON;
- final structured `kaon_FullAnalysis_Q4p4W2p74_lowe.json` or the narrowest structured artifact needed to independently audit the E.8.4 source;
- no-empirical-residual correction ledger JSON/CSV;
- current candidate F.3 JSON;
- current candidate F.4 JSON;
- applicable run summary/provenance output.

Do not require the enormous yield-data PDF unless it is actually needed to validate the new pages.

The profile must pin the **new reviewed analysis source commit after the user pushes this change** through the established source-identity mechanism. Because that SHA is not known during local implementation, the task may use the repository's existing profile pattern for candidate source identity but must leave a clear push-stable re-pin mechanism if the existing schema requires an exact pushed SHA.

Prefer a profile design that can be source-reviewed before push without requiring another scientific implementation cycle. Do not weaken provenance.

### 11.2 Tracked farm owner

Create:

```text
testing/run_e8_4_fix5_left_lowe_plot_gate.py
```

It must own the complete farm operation:

```text
preflight
-> invoke existing ./run_Prod_Analysis.sh -d 4p4 2p74
-> verify successful Left/lowe completion
-> verify fresh procedure PDF and manifest
-> verify renderer_failures == []
-> verify every required new page ID exists exactly once
-> verify all existing critical E.8.4 page IDs remain present
-> package through the reviewed generic collector/profile
-> verify ZIP integrity/manifest
-> print exactly one returned ZIP path
```

The user should not need a second manual bundle command after the multi-hour run.

The owner must:

- require exact pushed bundle/source commit;
- require the intended repository branch/HEAD;
- verify pre-run source provenance;
- record pre-run worktree status;
- preserve the existing `-d` semantics: paired low/high preflight, full analysis only Left/lowe, stop before high-epsilon full processing;
- not clean/reset/stash/delete farm-local files;
- fail closed if the run fails;
- fail closed if the expected new pages are missing;
- fail closed if renderer failures are nonempty;
- use a unique timestamped ZIP under `/volatile/hallc/c-kaonlt/trottar/globus`;
- never overwrite an existing ZIP.

Do not put the eight-hour analysis behind an unreviewed shell sequence.

---

## 12. Positive tests

Add deterministic tests proving at least:

1. E.8.4 SIMC comparison consumes passed/captured authoritative SIMC histograms and never opens a SIMC file.
2. SIMC normalization/scaling is not modified by E.8.
3. all three t parents retain all nine canonical phi children.
4. populated child panels contain Method-A MM and SIMC.
5. baseline/Method-A/SIMC page contains all three passed histograms.
6. Lambda-window markers are identical to the existing authoritative E.8.4 window.
7. empty/unavailable child is explicit.
8. yield-summary page consumes stored values and does not integrate a histogram.
9. undefined `DeltaY/Y0` remains undefined.
10. signed `DeltaY` is preserved.
11. parent-closure page consumes stored closure quantities.
12. existing E.8.4 page IDs remain unchanged and present.
13. new page IDs appear exactly once.
14. renderer failure accounting includes failures from the new pages.
15. Method B never enters any new numerical payload.
16. production objects are not mutated.
17. existing `MM_0` and `MM_A` payloads remain unchanged.
18. existing debug launcher behavior remains unchanged.
19. farm owner refuses wrong HEAD.
20. farm owner refuses dirty tracked source before run when the dirt touches gate-relevant source; harmless farm-local output/model artifacts must be handled according to the established farm-provenance policy, not deleted.
21. farm owner refuses missing new page IDs.
22. farm owner refuses nonempty renderer failures.
23. farm owner refuses an existing output ZIP.
24. farm owner prints one ZIP only after run + verification + collection succeed.
25. profile has only Left/lowe setting scope and exact required artifacts.

---

## 13. Negative / regression tests

Prove fail-closed behavior for:

- mismatched t/phi geometry between Method-A data and SIMC;
- missing SIMC child histogram;
- mismatched MM binning/axis where direct overlay would be invalid;
- malformed optional E.8.4 source;
- missing stored yield;
- undefined fractional yield;
- wrong candidate F.3/F.4 identity where already validated by the existing path;
- Method-B contamination;
- attempted SIMC shape renormalization in E.8;
- attempted child-by-child normalization;
- failed analysis return code;
- stale procedure PDF/manifest;
- missing manifest page;
- ZIP collection failure.

Run existing E.8.2/E.8.3/E.8.4/F.6.3 regressions.

---

## 14. Local validation

Use the repository-authoritative interpreter from `TOOLS.md`.

Run at minimum:

```bash
<PYTHON> -B -m py_compile \
  src/cuts/full_background_subtraction_plots.py \
  src/cuts/rand_sub.py \
  testing/run_e8_4_fix5_left_lowe_plot_gate.py \
  testing/test_full_background_subtraction_plots.py \
  testing/test_e8_4_production_impact_audit.py \
  testing/test_run_e8_4_fix5_left_lowe_plot_gate.py \
  testing/test_pion_hgcer_validation_bundle_profile_e8_4_fix5.py

<PYTHON> -B -m unittest testing.test_e8_4_production_impact_audit -v
<PYTHON> -B -m unittest testing.test_full_background_subtraction_plots -v
<PYTHON> -B -m unittest testing.test_f6_3_parallel_full_procedure_method_a -v
<PYTHON> -B -m unittest testing.test_run_prod_analysis_debug_left_low -v
<PYTHON> -B -m unittest testing.test_run_e8_4_fix5_left_lowe_plot_gate -v
<PYTHON> -B -m unittest testing.test_pion_hgcer_validation_bundle_profile_e8_4_fix5 -v
```

Run collector/profile regressions if the new profile exercises shared collector behavior.

Then:

```bash
<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/memory_bootstrap.py --root .
<PYTHON> -m unittest testing.test_memory_health
<PYTHON> -B tools/check_memory_health.py --root .
git -c core.safecrlf=false diff --check
```

`--fail-on-warning` is not the ordinary task-final gate. It may be run additionally for this milestone if zero-warning cleanliness is desired, but its success must not be rewritten back into the universal policy.

Do not run `main.py` locally and do not claim PyROOT/farm rendering validation from local tests.

---

## 15. Memory updates

### 15.1 CURRENT.md

Make CURRENT concise and accurate.

It must record:

- current pushed/starting source identity `da38444e7aa60efd62d6638780776344daf40276`;
- Refresh.2 materialization narrow closure;
- F.6.3/E.8.4 Left/lowe narrow runtime closure from the supplied package;
- Method A still detached;
- Method B still diagnostic-only;
- canonical-five E.8 not closed;
- this new E.8.4.Fix.5 shareable plot task as the active implementation.

Exactly one ordinary NEXT:

```text
NEXT — after user commit/push and pushed-state synchronization, run the tracked Q4p4W2p74 / Left / lowe E.8.4.Fix.5 farm owner to regenerate the full-analysis procedure PDF with all shareable Method-A impact pages and return its single fresh validation ZIP.
```

Do not make commit/push itself the NEXT.

### 15.2 MEMORY.md

Add only durable technical facts:

- Left/lowe current-baseline Method-A full-analysis path has direct runtime evidence;
- real `(t,phi)` yields change while parent-t normalization remains preserved;
- current evidence is narrow Left/lowe only;
- E.8 shareable comparison now requires final Method-A MM versus authoritative SIMC for every canonical child without renormalization.

Do not copy detailed phase chronology.

### 15.3 USER.md

Add the stable collaboration preference:

> When the user requests new KaonLT analysis/procedure plots that belong in tracked E.8 or another repository output, implement them through the normal source-changing Codex workflow and farm-render them from the repository. Chat-generated/extracted plots are inspection aids only when explicitly requested; they are not substitutes for tracked analysis deliverables.

Keep this concise.

### 15.4 LEARNINGS.md

Add the generalized lesson:

> A request for a new procedure-PDF analysis figure is a source/presentation task, not a request for ad-hoc plotting in ChatGPT. Trace authoritative inputs, implement the page in tracked source, test locally, then farm-render it.

### 15.5 TOOLS.md

Repair only the stale memory-health wording so it agrees with CODEX/MAINTENANCE/template:

- ordinary task-final health: `tools/check_memory_health.py --root .`;
- strict `--fail-on-warning` only for explicit milestone/zero-warning or blocking-warning cases.

Do not change unrelated tool/path facts.

### 15.6 Roadmap / E.8 decision

Update E.8.4 closure/presentation requirements to include:

- Method-A final clean-kaon MM versus authoritative SIMC per canonical `(t,phi)`;
- baseline/Method-A/SIMC direct overlay;
- shareable Y0/YA, DeltaY, and DeltaY/Y0 summaries versus phi;
- parent-closure sanity page.

This is presentation-only and does not change F.6.4 promotion ownership.

### 15.7 Evidence record

Create:

```text
docs/memory/evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md
```

Record the exact bundle/source/F.3/F.4 identities above, page-manifest result, narrow runtime conclusions, representative yield effects, parent closure, and packaging-time worktree caveat.

### 15.8 Phase record

Create:

```text
docs/memory/phases/e8-4-fix5-shareable-method-a-impact-pages.md
```

Record this task's source scope, presentation ownership, new page groups, tracked farm owner/profile, status, and exact NEXT.

---

## 16. Diff audit

Before stopping:

```bash
git status --short --untracked-files=all
git -c core.safecrlf=false diff --stat
git -c core.safecrlf=false diff --check
git -c core.safecrlf=false diff --no-ext-diff
```

Verify no frozen production/scientific file changed.

Create one fresh repository-root review bundle containing:

- branch/HEAD/status;
- complete tracked diff;
- complete `git diff --no-index /dev/null ...` for every new/untracked file;
- changed-path inventory;
- all local check outputs;
- memory-health report and byte counts;
- manifest result.

Use:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

Do not stage merely for review.

---

## 17. Acceptance criteria

PASS only if:

1. starting HEAD is exactly `da38444e7aa60efd62d6638780776344daf40276`;
2. supplied Left/lowe runtime evidence is recorded without overclaiming canonical-five or production closure;
3. TOOLS health policy is realigned;
4. new E.8 pages cover all 27 canonical `(t,phi)` cells;
5. final Method-A MM is directly compared to authoritative SIMC for every valid cell;
6. baseline/Method-A/SIMC overlay exists for every valid cell;
7. no SIMC normalization or shape normalization changes;
8. no renderer reloads/refills SIMC;
9. yield summary consumes stored Y0/YA/DeltaY/fraction values only;
10. parent closure consumes stored values only;
11. all existing E.8 pages remain available;
12. Method B remains excluded;
13. Method A remains non-production;
14. no production physics source changes;
15. debug Left/lowe semantics remain unchanged;
16. a tracked one-command farm owner exists and performs run -> verify -> package;
17. the user will not need a separate manual post-run bundle command;
18. new profile is Left/lowe scoped and provenance-safe;
19. deterministic local tests pass with exact skips reported;
20. memory manifest/health/diff checks pass;
21. one fresh complete review bundle is returned;
22. Codex does not commit, push, or run the farm.

---

## 18. Hard stop

After implementation, local deterministic tests, warranted memory updates, manifest/health checks, and creation of one review bundle:

**STOP.**

Do not:

- run the farm;
- commit or push;
- broaden to Left/highe or canonical five;
- make F.6.4 promotion decisions;
- modify scientific normalization/cuts/fits/yields;
- create ad-hoc plot files as substitutes for E.8;
- start another memory reconciliation phase.

Return the review bundle for independent ChatGPT actual-diff review.

After PASS, the user commits/pushes, ChatGPT performs one lightweight pushed-state synchronization review, and the user runs **one tracked Left/lowe farm owner** that both regenerates the new E.8 procedure PDF and packages the fresh evidence.
