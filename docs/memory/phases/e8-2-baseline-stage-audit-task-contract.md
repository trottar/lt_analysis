# KaonLT E.8.2 — Current Baseline Stage Audit Implementation Contract

## Objective

Implement **E.8.2 — baseline full-analysis stage audit** for the accepted current
KaonLT analysis chain.

E.8.2 is presentation/diagnostic infrastructure only. It must make the authoritative
baseline kaon missing-mass chain visually readable in canonical `(t,phi)` bins:

```text
selected prompt kaon
    -> random subtraction
    -> dummy subtraction
    -> accepted slow-proton cleaning
    -> accepted baseline pion-background subtraction using w0
    -> final baseline clean kaon MM_0(t,phi)
    -> final baseline Y_0(t,phi)
```

For every subtractive stage, the procedure presentation must show:

```text
authoritative input spectrum
    -> exact component/treatment removed
    -> authoritative output spectrum
```

This task must **not** implement Method A numerically. E.8.3/F.6.3 will follow
E.8.2 and will add the explicitly required `w0 -> w0*C` reweighting and the direct
baseline-versus-reweighted pion-background comparison.

This task must **never** activate, rerun, tune, consume, or present the historical
empirical residual Fit 1 / Fit 2 machinery. The accepted active profile is
`no_empirical_residual`.

No production physics, cuts, weights, normalization, binning, pion model, proton
model, yield extraction, SIMC, efficiencies, or cross-section behavior may change.

---

## Exact starting state

Required branch:

`test`

Required starting HEAD:

`5b2ba7f92dc51b303257dbf209df687415364602`

Commit message:

`Correct E8 no-empirical-residual chain`

Parent:

`942444fb9586dc53e847039c97c6c3b52fe9a351`

Codex must begin with:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
```

The worktree must be clean except for this newly placed untracked task-contract file.

Hard stop on any different branch, different HEAD, or unrelated pre-existing worktree
change.

Do not pull, reset, stash, clean, switch branches, commit, push, or run the Jefferson
Lab farm.

---

## Required startup reading

Read repository memory in the required order:

1. `docs/memory/AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read:

- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/decisions/e8-no-empirical-residual-chain-correction.md`
- `docs/memory/phases/phase-f6-method-a-production-promotion.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/evidence/e8-1-fix6-left-lowe-runtime-closure.md`
- this task contract.

Before editing, inspect the exact live source/runtime path around:

- `src/utility/background_config.py`
  - active background profile;
  - resolved empirical residual scales.
- `src/binning/calculate_yield.py`
  - `_process_yield_data_tree`;
  - proton-cleaning factor lookup/application;
  - prompt/random/dummy/dummy-random event filling;
  - random and dummy normalization/subtraction;
  - baseline pion subtraction;
  - authoritative component-subtraction payload;
  - final wide and cut-window missing-mass histograms;
  - `calculate_yield_data`;
  - final `Y_0(t,phi)` and its existing uncertainty;
  - existing legacy `stage_window_yields`.
- `src/cuts/rand_sub.py`
  - construction of D.6-E.8 procedure payloads;
  - current full-background procedure render;
  - current page-manifest write;
  - returned per-setting histogram state.
- `src/cuts/full_background_subtraction_plots.py`
  - current payload ownership/cloning;
  - current E.8 renderer and handoff ordering;
  - page-manifest helpers.
- `src/main.py`
  - Step 3 `rand_sub(...)`;
  - `histlist` lifetime;
  - Step 6 `find_yield_data(histlist, inpDict)`;
  - point immediately after authoritative data-yield production.

Do not redesign adjacent architecture.

---

# Source facts already independently established

These facts define the implementation boundary.

## 1. Active analysis explicitly disables the historical empirical residual fits

At the starting HEAD, `src/utility/background_config.py` has:

```python
BG_STAT_SCALE1 = 0.0
BG_STAT_SCALE2 = 0.0
BG_OPT_ACTIVE_PROFILE = "no_empirical_residual"
```

and the active profile forces both empirical residual scales to zero.

The old Fit-1/Fit-2 blocks remain in `calculate_yield.py` only as dormant historical
source machinery behind positive-scale conditionals.

E.8.2 must not call, expose, label, snapshot, test as an active stage, or otherwise
resurrect those fits.

## 2. Current canonical yield-level data fills already apply slow-proton cleaning

`_process_yield_data_tree(...)` obtains the accepted event-level
`final_cleaned_factor` from the setting proton-cleaning result and fills the ordinary
yield histograms with that factor.

Therefore the current yield-level:

- raw prompt snapshot;
- after-random snapshot;
- after-dummy / `H_MM_nosub_DATA`;

are already **post-proton-cleaning** when the cleaner is accepted.

They must not be relabeled as pre-proton.

## 3. The existing event traversal already has the information needed to capture the
true pre-proton chain

For each accepted data/dummy/random event, the authoritative traversal already knows:

- selected event identity;
- shifted missing mass;
- canonical t index;
- canonical phi index;
- source role;
- accepted proton-cleaning factor;
- random-window normalization ownership;
- data/dummy normalization ownership.

A true pre-proton E.8.2 stage can therefore be captured during the **same event
traversal** without a new replay/tree physics loop and without modifying the production
fill.

## 4. The current pion payload already contains the authoritative before/component/after
objects needed for the baseline pion stage

The accepted component-subtraction path already retains the exact baseline:

- wide pre-pion kaon spectrum;
- exact applied pion subtraction template;
- wide post-pion kaon spectrum;

and corresponding cut-window objects.

E.8.2 must consume those exact objects. It must not reconstruct the pion background in
the renderer.

## 5. Under `no_empirical_residual`, the post-pion kaon spectrum is the final baseline
kaon spectrum

There is no accepted empirical residual stage after pion subtraction.

The wide post-pion spectrum is the final baseline spectrum for presentation, while the
existing cut-window `H_MM_DATA` is the authoritative histogram used by the existing
yield extraction.

## 6. Final yield extraction occurs later than the current procedure-PDF render

The current full-background procedure PDF is rendered in Step 3 inside `rand_sub(...)`.

Canonical data yields are produced later in Step 6:

```python
yieldDict.update(find_yield_data(histlist, inpDict))
```

`calculate_yield_data(...)` stores the authoritative per-cell processed dictionary on
the setting histogram and then computes the existing final yield.

Therefore E.8.2 final MM/yield pages must be finalized **after Step-6 data-yield
production**, without rerunning any scientific calculation.

---

# Scientific ownership and frozen interfaces

## Frozen production baseline

The following are frozen:

- event selection and cuts;
- shifted missing-mass definition;
- shifted t definition;
- canonical t and phi edges;
- random-window definition and `nWindows`;
- random subtraction coefficients;
- data normalization;
- dummy normalization;
- dummy subtraction coefficients;
- slow-proton PID model;
- slow-proton event-factor calculation;
- setting-wide proton gate/application;
- baseline pion-control population;
- baseline pion event weight `w0`;
- pion component definitions;
- pion-component fits and accepted parent policy;
- pion fit/control windows;
- pion amplitudes/scales;
- pion application/fallback policy;
- final baseline yield calculation;
- existing yield error calculation;
- SIMC;
- efficiencies;
- acceptance;
- L/T separation;
- cross-section formulas;
- all accepted F.1-F.6.2 artifacts/fingerprints.

E.8.2 may observe, clone, summarize, and render authoritative objects. It may not alter
them.

## No empirical residual fits

Freeze:

`BG_OPT_ACTIVE_PROFILE = "no_empirical_residual"`

E.8.2 must not modify `src/utility/background_config.py`.

E.8.2 must not:

- call `bg_fit`;
- resolve a positive residual-fit scale for presentation;
- create an empirical residual component;
- use `H_MM_fit1sub_*` as a scientific stage;
- label any page or yield row Fit 1 / Fit 2;
- expose dormant residual-fit uncertainty terms as an E.8 stage.

If the runtime active profile is not `no_empirical_residual`, the E.8.2 payload must
fail closed/unavailable rather than silently presenting the wrong chain.

## No Method A / Method B numerical use

E.8.2 is the **baseline audit only**.

It must not:

- obtain or apply `C_j`;
- multiply by `C_j`;
- construct `B_pi^A`;
- construct `K_A`;
- modify `w0`;
- use Method B numerically.

Existing E.8.1/F.6.2 pages remain frozen and may still appear in the procedure PDF as
their already accepted detached presentation.

---

# Allowed files

Source edits are limited to:

- `src/binning/calculate_yield.py`
- `src/cuts/rand_sub.py`
- `src/cuts/full_background_subtraction_plots.py`
- `src/main.py`

Focused test edits/additions are limited to:

- `testing/test_full_background_subtraction_plots.py`
- new `testing/test_e8_2_baseline_stage_audit.py`

Warranted repository-memory edits are limited to:

- `docs/memory/CURRENT.md`
- new `docs/memory/phases/e8-2-baseline-stage-audit.md`
- `docs/memory/manifest.json`
- this contract:
  `docs/memory/phases/e8-2-baseline-stage-audit-task-contract.md`

Do not modify:

- `src/utility/background_config.py`;
- validation profiles;
- collectors;
- farm wrappers;
- launchers;
- accepted F.6.2 JSON/artifacts;
- any other production/scientific source.

If implementation requires another source file, hard stop and report why.

---

# E.8.2 authoritative producer contract

## A. Capture the true pre-proton source spectra during the existing traversal

Do **not** add a second data-tree event loop for E.8.2.

Extend the existing authoritative `_process_yield_data_tree(...)` path so that, for the
same accepted event and canonical `(t,phi)` assignment, E.8.2 can retain detached
wide-MM diagnostic contributions for:

1. **pre-proton contribution**
   - event diagnostic factor `1.0`;

2. **production post-proton contribution**
   - existing accepted `proton_cleaning_factor`;

3. **proton-removed contribution**
   - diagnostic factor `1.0 - proton_cleaning_factor`.

The production histogram fill remains exactly the current fill with
`proton_cleaning_factor`.

The E.8.2 auxiliary fills are presentation diagnostics only.

When no accepted proton cleaning applies:

- pre-proton and post-proton must be identical;
- proton-removed must be identically zero.

Use the same wide missing-mass axis as the current pion-stage wide-MM audit so all
baseline stages can be compared consistently.

Do not divide by the proton factor or infer pre-proton content from a post-proton
histogram.

## B. Preserve the four authoritative source roles

Build E.8.2 early-stage spectra from the same existing source roles and normalizations:

- prompt data;
- random data;
- dummy prompt;
- dummy random.

Do not invent a new sign convention.

The resulting stage semantics must be:

### Stage 0 — prompt

`prompt_pre_proton`

Normalized prompt data before random subtraction and before proton cleaning.

### Stage 1 — random subtraction

Retain:

- `random_component_pre_proton`
- `after_random_pre_proton`

Require bin-by-bin closure:

```text
prompt_pre_proton
    - random_component_pre_proton
    = after_random_pre_proton
```

### Stage 2 — dummy subtraction

Retain:

- `dummy_component_pre_proton`
- `after_dummy_pre_proton`

The dummy component must use the established normalized dummy-prompt minus
dummy-random treatment.

Require:

```text
after_random_pre_proton
    - dummy_component_pre_proton
    = after_dummy_pre_proton
```

### Stage 3 — slow-proton cleaning

Retain:

- `proton_component_removed`
- `after_proton`

`after_proton` must be the actual production post-random/post-dummy proton-cleaned wide
MM spectrum **before pion subtraction**.

Require:

```text
after_dummy_pre_proton
    - proton_component_removed
    = after_proton
```

The removed proton component must be derived from the same authoritative source-event
traversal and accepted event factors. Do not fit a new proton shape. Do not substitute a
setting/t-level proton diagnostic for this canonical `(t,phi)` closure.

---

# E.8.2 baseline pion-stage contract

For every canonical `(t,phi)` child, consume the exact existing accepted pion-application
objects.

Retain:

- `pion_input = K_proton`
- `pion_component_removed = B_pi^0`
- `after_pion_final = K_0`

where:

```text
K_0 = K_proton - B_pi^0
```

and the baseline event contribution remains:

```text
b_j^0 = s_j * w0_j
```

Use the exact applied production pion template associated with the accepted child
application/fallback state.

Require bin-by-bin closure where subtraction is valid:

```text
pion_input
    - pion_component_removed
    = after_pion_final
```

For `zero`, `skip_bin`, rejected, or otherwise unavailable child states:

- preserve the exact authoritative status/reason;
- do not manufacture a pion spectrum;
- do not silently drop the child.

Do not reconstruct `B_pi^0` in the renderer as `input - output` when the exact applied
template exists.

---

# Final baseline MM and yield contract

## Final wide MM

The E.8.2 final baseline presentation spectrum is:

`after_pion_final = MM_0(t,phi)`

on the existing wide-MM axis.

There is **no** empirical residual-fit stage after this.

## Existing cut-window yield histogram

The existing cut-window `H_MM_DATA` remains the authoritative yield-extraction input.

Do not replace or modify it.

## Final extracted yield

Retain the existing authoritative:

`Y_0(t,phi)`

exactly as currently returned by `calculate_yield_data(...)`.

Do not change the public yield dictionary schema merely for presentation.

Add a private E.8.2 presentation record on the setting histogram or per-cell source
structure containing:

- final existing yield value;
- statistical error from the existing `integral_with_stat_error(final_hist)` call;
- existing total yield error already computed by the production code;
- validity/status.

The existing production yield and existing production yield error must remain
bit-for-bit/numerically unchanged for deterministic fixture inputs.

Do not mislabel the existing total yield uncertainty as statistical uncertainty.

---

# E.8.2 stage-window diagnostic integrals

Do **not** repurpose or rename the existing legacy `stage_window_yields` field. Its
existing keys/semantics must remain unchanged for compatibility.

Add a separate E.8.2-specific ordered summary computed producer-side from the correct
current-chain snapshots:

1. `prompt_pre_proton`
2. `after_random_pre_proton`
3. `after_dummy_pre_proton`
4. `after_proton`
5. `after_pion_final`

Integrate with the established Lambda signal window.

These are **diagnostic stage-window integrals**.

They are not a replacement for the final extracted yield.

The renderer must not perform these integrals itself.

---

# Private E.8.2 source/payload schema

Use a new explicit versioned schema.

Recommended source schema:

`e8_2_baseline_stage_source/v1`

Recommended presentation schema:

`full_background_subtraction_e8_2/v1`

The source/payload must include at minimum:

- `available`;
- literal reason when unavailable;
- `non_authoritative = True` for presentation;
- `production_objects_mutated = False`;
- active profile;
- setting;
- epsilon;
- canonical t edges;
- canonical phi edges;
- missing-mass edges;
- Lambda window;
- ordered child inventory;
- per-child validity/reason;
- exact stage histograms;
- stage-window diagnostic integrals;
- final `Y_0`;
- statistical uncertainty;
- existing total uncertainty;
- pion application status/reason.

Fail closed if:

- active profile is not `no_empirical_residual`;
- canonical geometry is malformed;
- a populated child claims a subtraction stage without the required authoritative
  source object;
- histogram binning is inconsistent across one stage triplet;
- final yield authority is missing for a valid populated child.

---

# Procedure-PDF integration

## Existing timing problem

The current procedure PDF is rendered in Step 3, before Step 6 produces final canonical
yield information.

Do not move or rerun scientific calculations.

## Required approach

Preserve the current Step-3 procedure render as the preliminary/legacy artifact path.

During Step 3, retain a **private transient detached render-state bundle** on each
setting histogram containing everything required to reproduce the current procedure
pages later without scientific recomputation.

The render state must include:

- the exact already-built D.6-E.8.1 presentation payloads;
- current E.8/F.6.2 frozen payload exactly as already consumed;
- page-manifest setting identity;
- final PDF path;
- final page-manifest path or information sufficient to derive it.

### Detachment requirement

The preserved Step-3 render state must not hold mutable aliases to production objects
that later Step-6 processing can mutate.

Use the repository's existing ROOT histogram cloning/ownership helpers.

Add a deterministic test that mutates an original source object after render-state
capture and proves the retained render state is unchanged.

The render-state bundle is transient/private:

- do not export it into normal analysis JSON;
- do not write it into production ROOT output;
- release/remove it after successful final E.8.2 rerender.

## Step-6 finalization point

Immediately after:

```python
yieldDict.update(find_yield_data(histlist, inpDict))
```

and before SIMC yield calculation, invoke a presentation-only E.8.2 finalizer.

For each applicable kaon setting:

1. consume the preserved detached Step-3 render state;
2. consume the authoritative E.8.2 source created by yield processing;
3. build the E.8.2 presentation payload;
4. rerender the **complete** full-background-subtraction procedure PDF from scratch to
   a private temporary path;
5. render all existing pages with unchanged existing payloads;
6. insert E.8.2 pages at the required position;
7. write a temporary complete page-manifest sidecar;
8. only after successful render/close/manifest construction, atomically replace the
   preliminary PDF and manifest;
9. release the transient Step-3 render state.

Do not:

- rerun random subtraction;
- rerun proton cleaning;
- rerun pion fits/subtraction;
- rerun yield extraction;
- reopen or append to a closed multipage ROOT PDF;
- use external PDF merge/concatenation tools.

## Failure behavior

Presentation failure must **not mutate production physics objects or yields**.

If finalization fails:

- remove incomplete temporary replacement files;
- preserve the existing preliminary PDF/manifest rather than corrupting them;
- record an explicit E.8.2 finalization failure in the setting's private presentation
  status/log;
- do not claim E.8.2 complete;
- do not silently create substitute data.

Do not change the ordinary production-analysis failure policy solely because an E.8.2
renderer failed.

---

# E.8.2 page contract

Retain every existing procedure page and its existing relative order.

Existing E.8.1/F.6.2 values and semantics are frozen.

Insert the new E.8.2 block:

**after the existing E.8.1 context/overlay/map pages and before the final E.8 handoff
page.**

The handoff remains the final page.

For each authoritative canonical t parent, render the following six pages.

## 1. Random-subtraction page

Rows:

canonical phi children in authoritative order.

Columns:

1. normalized prompt pre-proton MM;
2. normalized random component removed;
3. after-random pre-proton MM.

## 2. Dummy-subtraction page

Columns:

1. after-random pre-proton MM;
2. normalized dummy contribution removed;
3. after-dummy pre-proton MM.

## 3. Slow-proton-cleaning page

Columns:

1. after-dummy pre-proton MM;
2. accepted proton component removed;
3. production after-proton MM.

## 4. Baseline pion-subtraction page

Columns:

1. proton-cleaned MM entering pion treatment;
2. exact applied baseline pion background `B_pi^0`;
3. final baseline clean kaon MM `K_0`.

## 5. Final baseline `MM_0(t,phi)` page

Use a compact canonical-phi grid.

Every valid child must show:

- canonical t index/edges;
- phi index/range;
- final wide clean-kaon MM spectrum;
- Lambda signal window;
- authoritative final `Y_0(t,phi)`;
- statistical uncertainty;
- existing total uncertainty, clearly distinguished;
- explicit invalid/unavailable reason where applicable.

## 6. Baseline stage-yield page

For every canonical phi child, show the producer-computed diagnostic Lambda-window
integrals:

- prompt pre-proton;
- after random;
- after dummy;
- after proton;
- after pion/final baseline clean sample.

Also show final authoritative `Y_0(t,phi)` separately.

Preferred layout:

- final `Y_0` vs phi with statistical error bars;
- a compact stage-integral table or clearly separated stage display.

Do not place stage-window diagnostic integrals and final extracted yield on a common axis
unless their normalization/units are explicitly identical.

---

# Plot semantics

These are **actual normalized analysis spectra**, not shape-only diagnostics.

Therefore:

- no unit-area normalization;
- no arbitrary display rescaling of one stage relative to another;
- no clipping of negative signed bins;
- no zeroing of negative bins;
- no smoothing;
- no interpolation;
- no independent child normalization.

Within each before/component/after phi row:

- use a common MM x axis;
- choose comparison-safe y limits;
- preserve signed content;
- visibly mark the Lambda signal window;
- show phi index/range.

The component column may have signed content after random/dummy algebra. Render it
faithfully.

---

# Page-manifest requirements

Extend the existing ordinary page manifest additively.

Do not break existing fields required by accepted E.8.1 validation.

For each E.8.2 page record at minimum:

- page ID;
- E.8.2 schema version;
- semantic stage;
- setting;
- epsilon;
- t index/edges;
- represented phi inventory;
- invalid/unavailable child inventory;
- `authoritative = False`;
- presentation-only flag.

Use stable page IDs, for example:

- `full_background.e8_2.random.t1`
- `full_background.e8_2.dummy.t1`
- `full_background.e8_2.proton.t1`
- `full_background.e8_2.pion.t1`
- `full_background.e8_2.final_mm.t1`
- `full_background.e8_2.stage_yields.t1`

Do not renumber or rename existing page IDs.

---

# Production-regression requirements

The implementation must prove that E.8.2 observation/presentation does not change the
baseline analysis.

At minimum verify deterministically:

1. ordinary production post-proton histograms are unchanged by enabling E.8.2 capture;
2. proton-cleaning factors are unchanged;
3. random/dummy normalized production spectra are unchanged;
4. baseline pion component payload is unchanged;
5. baseline pion-subtracted production histograms are unchanged;
6. existing public `stage_window_yields` remains unchanged;
7. existing returned `Y_0(t,phi)` remains unchanged;
8. existing returned yield error remains unchanged;
9. existing E.8.1/F.6.2 payload values are unchanged;
10. no historical empirical residual fit is called or enabled;
11. no Method-A or Method-B numerical dependency enters E.8.2.

---

# Positive deterministic tests

Add focused tests covering at least:

## Early-chain event capture

- accepted proton factor `1.0`:
  - pre-proton == post-proton;
  - proton removed == zero.
- nontrivial accepted proton factor:
  - correct pre-proton fill;
  - correct post-proton production fill;
  - correct removed component;
  - no second data-tree event traversal.
- exact random normalization closure.
- exact dummy normalization/subtraction closure.
- exact proton closure.
- canonical t/phi identity preserved.

## Pion stage

- exact accepted pion template is used.
- before - component = after.
- wide final `after_pion_final` equals the accepted post-pion object.
- child fallback/invalid reason is preserved.
- no empirical residual stage exists in the E.8.2 payload.

## Yield stage

- producer-side stage integral sequence is exactly:
  - prompt;
  - random;
  - dummy;
  - proton;
  - pion/final.
- existing legacy `stage_window_yields` is unchanged.
- final `Y_0` matches the existing returned production value.
- statistical error is the existing histogram statistical error.
- existing total yield error remains unchanged and separately labeled.

## Active-profile guard

- `no_empirical_residual` is accepted.
- any non-`no_empirical_residual` profile makes E.8.2 unavailable/fail closed.
- a test spy/mock proves `bg_fit` is never called by any E.8.2 producer, builder, or
  renderer helper.

## Render-state detachment

- mutate original source histogram/payload after Step-3 render-state capture;
- retained render state remains unchanged.

## Full rerender orchestration

- existing pages retain relative order;
- E.8.2 pages occur after E.8.1 parent pages and before handoff;
- handoff remains last;
- every canonical phi child is represented or explicitly invalid/unavailable;
- successful temporary render atomically replaces preliminary PDF/manifest;
- failed temporary render preserves preliminary output and records failure;
- transient render state is released after success.

---

# Negative tests

Require fail-closed/unavailable behavior for:

- missing true pre-proton source support;
- malformed canonical t/phi geometry;
- inconsistent MM binning in a before/component/after triplet;
- missing exact pion template for a child that claims accepted pion subtraction;
- missing final wide post-pion spectrum;
- missing final yield for a valid populated child;
- malformed E.8.2 schema;
- missing Step-3 render state at finalization;
- active profile not equal to `no_empirical_residual`;
- attempted inclusion of Fit 1 / Fit 2 as E.8.2 stages;
- attempted Method-A or Method-B numerical application.

Do not add fallback reconstruction that hides missing authority.

---

# Existing regression tests to run unchanged

Run:

- `testing.test_full_background_subtraction_plots`
- `testing.test_proton_cleaning_runtime_status`
- `testing.test_pion_t_bin_particle_stage`
- `testing.test_pion_hgcer_phase_f_runtime_contract`

Do not edit the latter three merely to make the new implementation pass.

Expected environment-dependent ROOT/PyROOT skips are acceptable only when they are
already legitimate conditional skips.

---

# New focused test module

Add:

`testing/test_e8_2_baseline_stage_audit.py`

Use fake histogram objects and mocks where practical.

Do not build a duplicate fake scientific analysis merely for testing.

Prefer extracting narrow pure helpers for:

- stage closure;
- source-role algebra;
- payload validation;
- page inventory/order;
- finalization transaction logic.

Production scientific functions remain the authority.

---

# Local deterministic validation

Use the repository-appropriate local interpreter discovered by Codex.

At minimum run:

```bash
python -B -m py_compile   src/binning/calculate_yield.py   src/cuts/rand_sub.py   src/cuts/full_background_subtraction_plots.py   src/main.py   testing/test_e8_2_baseline_stage_audit.py   testing/test_full_background_subtraction_plots.py

python -B -m unittest testing.test_e8_2_baseline_stage_audit -v
python -B -m unittest testing.test_full_background_subtraction_plots -v
python -B -m unittest testing.test_proton_cleaning_runtime_status -v
python -B -m unittest testing.test_pion_t_bin_particle_stage -v
python -B -m unittest testing.test_pion_hgcer_phase_f_runtime_contract -v
```

If a different Python executable is required, report it exactly.

Also run:

- established memory manifest regeneration;
- memory manifest check;
- memory health/integrity;
- memory bootstrap check;
- `git diff --check`.

Local tests do not establish ROOT/PyROOT/full-main/farm validation.

---

# Memory update

Add:

`docs/memory/phases/e8-2-baseline-stage-audit.md`

Record:

- starting HEAD;
- exact scope;
- true pre-proton capture path;
- baseline pion authority;
- post-yield finalization path;
- explicit `no_empirical_residual` boundary;
- explicit absence of Method A/Method B numerical application;
- deterministic checks actually run;
- farm boundary.

Update `docs/memory/CURRENT.md` narrowly:

- E.8 remains `ACTIVE`;
- E.8.2 becomes `ACTIVE`;
- local implementation awaits independent actual-diff/runtime-path review;
- E.8.3 remains `BLOCKED`;
- F.6.3 remains `BLOCKED`;
- no farm/runtime claim is added.

Regenerate:

`docs/memory/manifest.json`

Do not alter historical runtime evidence.

Do not mark E.8.2 `SOURCE REVIEWED`; ChatGPT owns that independent status after actual
diff inspection.

---

# Farm boundary

No Jefferson Lab farm run is requested.

The user explicitly chose milestone-based farm validation rather than farm execution
after every source update.

This E.8.2 change may proceed through local deterministic tests, independent source
review, and user commit/push without an immediate farm run.

Do not modify a validation bundle/profile in this task.

---

# Diff audit

Before finishing:

```bash
git status --short
git diff --stat
git diff --check
git diff --name-only
```

Audit explicitly for absence of changes to:

- `src/utility/background_config.py`;
- random/dummy coefficients;
- proton factor calculation;
- pion weights;
- pion components/fits/windows;
- canonical binning;
- yield formulas;
- Method A;
- Method B numerical use;
- dormant empirical residual fit configuration.

Because this implementation is expected to be substantial, create:

`kaonlt_review.diff`

in the repository root containing:

1. complete tracked `git diff`;
2. complete `git diff --no-index -- /dev/null ...` entries for every intended
   untracked file.

At minimum include the new:

- `testing/test_e8_2_baseline_stage_audit.py`
- `docs/memory/phases/e8-2-baseline-stage-audit.md`
- `docs/memory/phases/e8-2-baseline-stage-audit-task-contract.md`

Do not add `kaonlt_review.diff` to Git.

---

# Acceptance criteria

The implementation is ready for independent ChatGPT review only if:

1. starting HEAD is exactly `5b2ba7f92dc51b303257dbf209df687415364602`;
2. only allowlisted paths change;
3. the production scientific chain remains unchanged;
4. no dormant empirical residual fit is activated, rerun, or presented;
5. a true pre-proton canonical `(t,phi)` wide-MM chain is captured during the same
   authoritative event traversal;
6. random before/component/after closes;
7. dummy before/component/after closes;
8. proton before/component/after closes;
9. baseline pion before/component/after uses the exact applied `B_pi^0`;
10. post-pion wide MM is explicitly the final baseline `MM_0(t,phi)`;
11. final public production `Y_0(t,phi)` and existing error are unchanged;
12. E.8.2 stores the statistical error separately without changing production output;
13. legacy `stage_window_yields` remains unchanged;
14. new E.8.2 stage integrals contain only prompt/random/dummy/proton/pion-final;
15. active profile is fail-closed to `no_empirical_residual`;
16. Step-3 render state is detached;
17. post-Step-6 finalization performs no scientific recalculation;
18. existing E.8.1/F.6.2 pages remain numerically/semantically unchanged;
19. E.8.2 pages appear before final handoff and handoff remains last;
20. every canonical child is represented or explicitly unavailable/invalid;
21. no Method-A or Method-B numerical application is introduced;
22. deterministic local and memory checks pass;
23. no farm/runtime claim is made.

---

# Hard stop

Stop and report rather than broadening scope if:

- true pre-proton canonical support cannot be captured without changing the accepted
  proton-cleaning calculation;
- the exact applied baseline pion template cannot be identified unambiguously;
- the final baseline yield cannot be mapped to the canonical child without changing
  yield logic;
- final procedure rerender requires rerunning a scientific calculation;
- an external PDF merge appears necessary;
- a non-`no_empirical_residual` active profile is discovered;
- another source file must be modified;
- any frozen scientific interface would have to change.

Codex must not commit, push, or run the farm.
