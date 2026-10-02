# KaonLT E.8.4 Fix.5.1 — unblock authoritative per-(t,phi) SIMC handoff, then implement shareable Method-A pages

## Objective

Resolve the concrete pre-implementation blocker found while executing
`e8-4-fix5-shareable-method-a-impact-pages-task-contract.md`, then complete that
E.8.4 shareable-plot task without changing production physics.

The blocker is source-architectural, not scientific:

```text
find_yield_data(...)
-> current E.8.2/E.8.4 finalization
-> find_yield_simc(...)
```

The authoritative normalized per-`(t,phi)` SIMC missing-mass histograms are
created and attached to each setting histogram only inside `find_yield_simc`.
Therefore the current E.8 finalizer runs too early to consume them.

The narrow repair is to make the already-existing E.8 finalizer run only after
the already-existing SIMC-yield calculation has populated
`hist["_xsect_support_simc"]`, and then pass/read that existing support into the
presentation-only E.8.4 consumer.

After that handoff is available, implement the shareable E.8.4 pages from the
original Fix.5 contract and the tracked Left/lowe run -> verify -> package owner.

No SIMC recalculation, reload, renormalization, production correction, or
Method-A promotion is allowed.

---

## Exact starting state

Required branch:

```text
test
```

Required committed HEAD:

```text
da38444e7aa60efd62d6638780776344daf40276
```

Required prior Codex blocker result:

```text
E.8.4.Fix.5 PRE-IMPLEMENTATION HARD STOP / BLOCKED
```

The prior task made no tracked source/testing/memory edits.

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short --untracked-files=all
git log --oneline -1
```

Expected pre-existing versionable files may include:

```text
docs/memory/phases/e8-4-fix5-shareable-method-a-impact-pages-task-contract.md
docs/memory/phases/e8-4-fix5-1-simc-runtime-handoff-unblock-task-contract.md
```

and no unrelated tracked changes.

Do not reset, stash, clean, commit, push, or run the farm.

---

## Mandatory startup read

Read in this order:

1. root `AGENTS.md`;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read:

- `docs/memory/CODEX.md`
- `docs/memory/COMMUNICATION.md`
- `docs/memory/MAINTENANCE.md`
- `docs/memory/LEARNINGS.md`
- `docs/memory/TOOLS.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/decisions/farm-validation-bundle-procedure.md`
- `docs/memory/phases/f6-3-current-baseline-candidate-lineage-adoption.md`
- `docs/memory/phases/e8-4-fix5-shareable-method-a-impact-pages-task-contract.md`

Then inspect the live source/tests around the blocker:

- `src/main.py`
- `src/binning/calculate_yield.py`
- `src/cuts/full_background_subtraction_plots.py`
- `src/cuts/rand_sub.py`
- `src/utility/bg_optimization.py`
- `testing/test_e8_4_production_impact_audit.py`
- `testing/test_full_background_subtraction_plots.py`
- `testing/test_f6_3_parallel_full_procedure_method_a.py`
- `testing/test_run_prod_analysis_debug_left_low.py`

Do not trust the previous Codex blocker summary by itself; confirm the source path.

---

## Audited blocker facts to preserve

The repair is based on these source facts at the required starting HEAD.

### Current caller order

`src/main.py` currently does:

```text
find_yield_data(histlist, inpDict)
finalize_full_background_subtraction_e8_2(hist) for each kaon setting
find_yield_simc(histlist, inpDict)
```

So E.8 finalization currently precedes the SIMC-yield producer.

### Authoritative per-cell SIMC producer

`src/binning/calculate_yield.py` already creates canonical per-`(t,phi)`:

```text
H_MM_SIMC_<t>_<phi>
```

It applies the existing:

```text
Scale(normfac_simc)
```

before cloning the histogram into:

```text
support_hist_dict["mm"][t][phi]
```

and `calculate_yield_simc(...)` stores the result on the production setting
histogram as:

```text
hist["_xsect_support_simc"] = binned_dict[kin_type]["support_hist_dict"]
```

These are the authoritative already-normalized child SIMC MM objects required
for E.8 presentation.

### Why the prepass support is not the right handoff

`src/utility/bg_optimization.py` can construct SIMC support for temporary
candidate dictionaries, but the selected final low-epsilon analysis reruns
`rand_sub` on final histograms. Do not wire temporary optimization-candidate
support into E.8 and do not treat it as final production-setting authority.

### Required architecture

The intended path is:

```text
existing find_yield_data
-> existing find_yield_simc
   -> existing normfac_simc scaling
   -> existing canonical support_hist_dict["mm"]
   -> existing hist["_xsect_support_simc"]
-> presentation-only E.8 finalizer
   -> validate/clone existing per-cell SIMC MM
   -> render E.8 pages
```

No second SIMC producer is permitted.

---

## Scope resolution

The original Fix.5 contract stopped because `src/main.py` was outside its
allowlist. This repair explicitly allows one narrow caller-order change in
`src/main.py`.

### Allowed scientific/presentation source files

```text
src/main.py
src/cuts/full_background_subtraction_plots.py
```

`src/cuts/rand_sub.py` is allowed only if a strictly presentation-state plumbing
change is proven necessary after the caller-order repair. Prefer no change.

### `src/main.py` ownership restriction

The only permitted behavioral change in `src/main.py` is:

- keep `find_yield_data(...)` where it is;
- run the existing `find_yield_simc(...)` before the E.8 finalization loop;
- run the existing `finalize_full_background_subtraction_e8_2(hist)` after
  SIMC support has been attached;
- preserve the same single call to each calculation;
- preserve all yield dictionaries/results;
- preserve correction-ledger ordering after E.8 finalization unless source
  dependencies prove otherwise;
- do not alter arguments, normalization, loops, cuts, or physics.

This is a presentation-consumer ordering repair only.

Do not invoke `find_yield_simc` twice.

---

## Frozen scientific source

Do not modify:

```text
src/binning/calculate_yield.py
src/utility/bg_optimization.py
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
src/cuts/pion_hgcer_method_a_parent_preserving_correction.py
src/cuts/pion_hgcer_method_a_tphi_propagation.py
src/cuts/pion_component_subtraction.py
src/cuts/pion_component_fits.py
src/cuts/particle_subtraction.py
binning/calculate_yield.py
run_Prod_Analysis.sh
```

Do not change any:

- SIMC event loop;
- SIMC normalization;
- data normalization;
- yield integration;
- pion fit;
- proton subtraction;
- F.4/F.5/F.6.3 mathematics;
- Lambda window;
- canonical t/phi binning;
- efficiency/acceptance/cross-section logic.

Method B remains diagnostic-only.
Method A remains detached/non-production.

---

## SIMC presentation handoff

Extend the E.8.4 presentation consumer so that it can consume the existing
`hist["_xsect_support_simc"]`.

Requirements:

1. use only `support_hist_dict["mm"]`;
2. require exactly the canonical t/phi geometry expected by E.8.4;
3. clone display histograms before styling;
4. require finite bin contents/errors;
5. require compatible MM binning/axis with `MM_0`/`MM_A`;
6. never call `Scale`, `Normalize`, `Integral` for a display renormalization, or
   any equivalent shape-normalization operation on SIMC;
7. never reopen a SIMC ROOT file;
8. never refill a SIMC histogram;
9. never use the setting-wide Step-5 `H_MM_SIMC` as a child substitute;
10. fail closed with an explicit E.8.4 unavailable reason if required child
    SIMC support is missing or geometrically incompatible.

The E.8.4 payload may include detached SIMC display clones alongside the
existing `MM_0` and `MM_A` clones.

No authoritative scientific payload is rewritten.

---

## Required E.8.4 pages

Complete the original Fix.5 presentation requirement.

### A. Method-A final MM versus SIMC

One 3x3 page per canonical t parent:

```text
full_background.e8_4.method_a_vs_simc.t1
full_background.e8_4.method_a_vs_simc.t2
full_background.e8_4.method_a_vs_simc.t3
```

Each valid pad overlays:

- `MM_A(t,phi)`;
- authoritative already-normalized SIMC MM for the same child;
- Lambda-window markers;
- t/phi identity;
- concise legend.

No display rescaling.

### B. Baseline / Method A / SIMC

One 3x3 page per t parent:

```text
full_background.e8_4.baseline_method_a_simc.t1
full_background.e8_4.baseline_method_a_simc.t2
full_background.e8_4.baseline_method_a_simc.t3
```

Overlay exactly:

- `MM_0`;
- `MM_A`;
- authoritative SIMC MM.

No fit or goodness-of-fit statistic is introduced here.

### C. Yield summaries

One shareable page per t parent:

```text
full_background.e8_4.yield_summary.t1
full_background.e8_4.yield_summary.t2
full_background.e8_4.yield_summary.t3
```

Each page shows all nine canonical phi bins:

1. `Y0` and `YA` vs phi;
2. `DeltaY = YA - Y0` vs phi;
3. `DeltaY / Y0` vs phi when defined.

Use the already validated E.8.4/F.6.3 values. Do not reintegrate MM histograms
for these summary curves.

Undefined fractional values remain visibly undefined.

### D. Parent closure

Add:

```text
full_background.e8_4.parent_closure
```

Show baseline signed parent, adjusted signed parent, and closure residual for all
three t parents. Label this explicitly as a parent-normalization sanity check.

### E. Preserve existing E.8 evidence

Retain all existing E.8.2/E.8.3/E.8.4 pages and IDs, including:

- baseline subtraction pages;
- `B_pi^0` versus `B_pi^A`;
- final `MM_0` versus `MM_A`;
- signed difference;
- existing yield-impact pages;
- setting summary.

---

## Presentation quality

These new pages are analysis/shareable E.8 output, not engineering diagnostics.

Requirements:

- all 27 canonical cells represented;
- 3x3 layout for per-child MM pages;
- readable titles and legends;
- t range in page header;
- phi range in each pad;
- consistent MM x-range for direct comparisons;
- Lambda-window markers;
- no giant provenance blocks on shareable figure pages;
- no overlapping title/legend/pad labels;
- detailed provenance remains on authority/manifest pages.

---

## Tests allowed

Modify:

```text
testing/test_e8_4_production_impact_audit.py
testing/test_full_background_subtraction_plots.py
testing/test_run_prod_analysis_debug_left_low.py
```

Create:

```text
testing/test_e8_4_fix5_main_order.py
testing/pion_hgcer_validation_bundle_profile_e8_4_fix5.json
testing/test_pion_hgcer_validation_bundle_profile_e8_4_fix5.py
testing/run_e8_4_fix5_left_lowe_plot_gate.py
testing/test_run_e8_4_fix5_left_lowe_plot_gate.py
```

Keep the generic collector unchanged unless a concrete incompatibility is
demonstrated; if so STOP rather than silently expanding scope.

---

## Required caller-order tests

Add deterministic source/runtime-path checks proving:

1. `find_yield_data` still occurs before SIMC-yield calculation;
2. `find_yield_simc` occurs exactly once;
3. E.8 finalization occurs after `find_yield_simc`;
4. E.8 finalization occurs exactly once per kaon hist;
5. correction-ledger export remains after E.8 finalization;
6. no new analyzer/SIMC producer is added.

Prefer a narrow AST/source-order test over importing/running full farm `main.py`.

---

## Required E.8/SIMC tests

Prove:

1. a complete 3x9 SIMC `mm` support matrix is accepted;
2. each SIMC display histogram is cloned;
3. input SIMC histogram integrals/bin contents are unchanged after payload build
   and rendering helpers;
4. there is no `Scale`/shape normalization in the E.8 handoff;
5. missing child SIMC fails closed;
6. wrong matrix dimensions fail closed;
7. wrong MM binning fails closed;
8. nonfinite SIMC content fails closed;
9. setting-wide SIMC fallback is not accepted;
10. Method-B data is never consulted;
11. existing `MM_0`, `MM_A`, Y0, YA, and existing E.8.4 page payloads remain
    unchanged apart from the addition of SIMC display support;
12. all 27 cells are retained;
13. all new page IDs appear exactly once;
14. all existing critical E.8 page IDs remain present;
15. renderer failures include failures of new pages.

---

## Tracked Left/lowe farm owner

Complete the original operational requirement before the next long run.

Create:

```text
testing/run_e8_4_fix5_left_lowe_plot_gate.py
```

It must own:

```text
preflight
-> invoke existing ./run_Prod_Analysis.sh -d 4p4 2p74
-> verify Left/lowe analysis success
-> verify fresh full-background-subtraction PDF + manifest
-> require renderer_failures == []
-> require every new E.8.4 page ID exactly once
-> require existing critical E.8.4 page IDs still present
-> package through the generic collector with the reviewed Fix.5 profile
-> verify ZIP/manifest integrity
-> print exactly one fresh ZIP path
```

Preserve existing debug semantics:

```text
paired low/high canonical preflight
full analysis only Left / lowe
stop before high-epsilon full processing
```

The owner must not clean/reset/stash/delete farm-local files.

It must verify the exact pushed source identity before the expensive run and
record pre-run worktree provenance.

The user must not need a separate post-run bundle command.

---

## Validation profile

Create a Left/lowe-only generic-artifact profile:

```text
testing/pion_hgcer_validation_bundle_profile_e8_4_fix5.json
```

Require at minimum:

- `Left_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction.pdf`;
- matching page-manifest JSON;
- `kaon_FullAnalysis_Q4p4W2p74_lowe.json`;
- no-empirical-residual correction-ledger JSON and CSV;
- current-baseline candidate F.3 JSON;
- current-baseline candidate F.4 JSON;
- narrow run summary/provenance artifact needed by the owner.

Do not require the giant yield-data PDF.

Use the repository's established exact-source provenance mechanism. Do not
weaken it to avoid the post-push SHA.

If the profile needs a literal pushed commit that cannot be known until after
user commit, structure the owner/profile so the user-controlled pushed SHA is
supplied/validated through the existing supported mechanism rather than
requiring a second source implementation phase.

---

## Milestone memory alignment

The previous Fix.5 contract identified stale memory and the supplied runtime
evidence. Perform that alignment in this implementation, not as a separate
cycle.

Allowed memory modifications:

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

Create:

```text
docs/memory/evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md
docs/memory/phases/e8-4-fix5-shareable-method-a-impact-pages.md
docs/memory/phases/e8-4-fix5-1-simc-runtime-handoff-unblock-task-contract.md
```

Retain the original Fix.5 contract as the record of the pre-implementation
boundary that exposed this blocker.

### Runtime evidence to record

Use the already-reviewed package:

```text
KaonLT_F6_3_Left_lowe_runtime_evidence_Q4p4W2p74_20261001-222142.zip

bundle SHA-256:
200fda66fe410274df1c9a8252b9e87114d8fd10a520fb7d71691ec6b3772874

source HEAD:
da38444e7aa60efd62d6638780776344daf40276

F.3 candidate SHA-256:
eeb480d9b7ddffab1e97c7f9f27c3099d0780c2bbb4d42a0ee68c077a5752b3d

F.4 candidate SHA-256:
1d545924eba89c7f9ffa28028e307aca9b434a89beec06863cf2893887b6b902
```

Record narrow statuses:

```text
F.4.Refresh.2 detached candidate materialization:
CLOSED / RUNTIME VALIDATED

F.6.3 current-baseline candidate lineage,
Q4p4W2p74 / Left / lowe:
CLOSED / RUNTIME VALIDATED

E.8.4 existing Left/lowe Method-A production-impact audit:
CLOSED / RUNTIME VALIDATED

E.8.4.Fix.5 shareable SIMC/yield presentation:
ACTIVE until the new farm PDF is reviewed
```

Do not close canonical-five E.8 and do not promote Method A.

Preserve the packaging-time worktree caveat:

```text
M src/models/xmodel_kaon_pl.f
?? src/kaon/functions/Q4p4W2p74.model
```

as post-run provenance only.

### USER/LEARNINGS alignment

Record the durable workflow rule:

- repository analysis/procedure plots are implemented through Codex/source and
  farm-rendered;
- ChatGPT-generated plots are inspection aids only when explicitly requested,
  not substitutes for tracked E.8 deliverables.

### TOOLS alignment

Repair the stale health wording:

```text
ordinary task-final:
<PYTHON> -B tools/check_memory_health.py --root .

strict --fail-on-warning:
only explicit milestone/zero-warning or materially blocking-warning audits
```

Do not make strict-warning mode universal again.

### CURRENT exact NEXT

Leave exactly one ordinary NEXT:

```text
NEXT — after user commit/push and pushed-state synchronization, run the tracked Q4p4W2p74 / Left / lowe E.8.4.Fix.5 farm owner to regenerate the full-analysis procedure PDF with the authoritative per-(t,phi) SIMC comparisons and shareable Method-A yield-impact pages, verify them, and return its single fresh validation ZIP.
```

---

## Local validation

Use repository-authoritative Python from `TOOLS.md`.

At minimum:

```bash
<PYTHON> -B -m py_compile \
  src/cuts/full_background_subtraction_plots.py \
  testing/test_e8_4_fix5_main_order.py \
  testing/run_e8_4_fix5_left_lowe_plot_gate.py \
  testing/test_run_e8_4_fix5_left_lowe_plot_gate.py \
  testing/test_pion_hgcer_validation_bundle_profile_e8_4_fix5.py

<PYTHON> -B -m unittest testing.test_e8_4_fix5_main_order -v
<PYTHON> -B -m unittest testing.test_e8_4_production_impact_audit -v
<PYTHON> -B -m unittest testing.test_full_background_subtraction_plots -v
<PYTHON> -B -m unittest testing.test_f6_3_parallel_full_procedure_method_a -v
<PYTHON> -B -m unittest testing.test_run_prod_analysis_debug_left_low -v
<PYTHON> -B -m unittest testing.test_run_e8_4_fix5_left_lowe_plot_gate -v
<PYTHON> -B -m unittest testing.test_pion_hgcer_validation_bundle_profile_e8_4_fix5 -v
```

Run applicable generic collector/profile tests.

Then:

```bash
<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/memory_bootstrap.py --root .
<PYTHON> -m unittest testing.test_memory_health
<PYTHON> -B tools/check_memory_health.py --root .
git -c core.safecrlf=false diff --check
```

Tests do not establish farm/PyROOT behavior.

---

## Diff audit

Before stopping:

```bash
git status --short --untracked-files=all
git -c core.safecrlf=false diff --stat
git -c core.safecrlf=false diff --check
git -c core.safecrlf=false diff --no-ext-diff
```

Confirm:

- `src/main.py` change is only the narrow Step-6 order repair;
- `src/binning/calculate_yield.py` unchanged;
- `src/utility/bg_optimization.py` unchanged;
- no production physics files changed;
- no SIMC scaling changed;
- no Method-A mathematics changed;
- no Method-B source changed.

Create one fresh complete repository-root review bundle:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

Include branch/HEAD/status, complete tracked diff, complete diffs for new files,
test outputs, memory-health/manifest output, and changed-path inventory.

Do not stage merely for review.

---

## Acceptance criteria

PASS only if all are true:

1. exact starting HEAD is correct;
2. the pre-implementation blocker is confirmed from source;
3. `find_yield_simc` remains a single existing calculation;
4. only its order relative to E.8 finalization changes in `src/main.py`;
5. E.8 finalization now sees the authoritative existing per-cell SIMC support;
6. no SIMC producer/normalization changes;
7. no setting-wide SIMC substitute is used;
8. all 27 canonical cells are validated;
9. final Method-A MM vs SIMC pages exist for t1/t2/t3;
10. baseline/Method-A/SIMC pages exist for t1/t2/t3;
11. Y0/YA, DeltaY, DeltaY/Y0 summary pages exist for t1/t2/t3;
12. parent-closure page exists;
13. existing E.8 pages remain;
14. Method A remains detached;
15. Method B remains diagnostic-only;
16. production objects/corrections are unchanged;
17. tracked farm owner performs run -> verify -> package;
18. no separate manual post-run bundle step is required;
19. Left/lowe-only profile is provenance-safe;
20. milestone memory is aligned without overclaiming canonical-five closure;
21. ordinary memory health/manifest/diff checks pass;
22. one complete review bundle is returned;
23. Codex does not commit, push, or run the farm.

---

## Hard stop

After implementation, deterministic local validation, warranted memory updates,
and creation of one complete review bundle:

**STOP.**

Do not commit, push, run the farm, broaden to other settings, promote Method A,
change scientific normalization/yields, or start another memory-only cycle.

The next workflow is:

```text
ChatGPT actual-diff review
-> user commit/push
-> ChatGPT pushed-state synchronization
-> one tracked Left/lowe farm run
-> owner verifies and packages automatically
-> ChatGPT evidence/PDF review
```
