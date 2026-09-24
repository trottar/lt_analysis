# E.8.2 — Baseline full-analysis stage audit

## Status

`SOURCE REVIEWED` — independent ChatGPT inspection of the complete cumulative
E.8.2/Fix.1/Fix.2/Fix.3 actual diff passed. This is source review only, not a
ROOT/PyROOT, full-analysis, procedure-PDF farm-rendering, farm, or runtime
acceptance claim.

## Starting identity and scope

- Required/observed branch: `test`.
- Required/observed starting HEAD:
  `5b2ba7f92dc51b303257dbf209df687415364602`
  (`Correct E8 no-empirical-residual chain`).
- The only pre-existing worktree item was the user-created E.8.2 task contract.
- Edited implementation scope: `src/binning/calculate_yield.py`,
  `src/cuts/rand_sub.py`, `src/cuts/full_background_subtraction_plots.py`, and
  `src/main.py`; focused test: `testing/test_e8_2_baseline_stage_audit.py`.
- Memory scope: this record, `CURRENT.md`, and regenerated `manifest.json`.

## Implemented presentation-only baseline audit

The accepted current baseline chain is represented exactly as:

```text
prompt/random -> dummy -> slow proton -> baseline pion background w0
-> final baseline clean kaon MM_0(t,phi) -> Y_0(t,phi)
```

- `calculate_yield._process_yield_data_tree` retains detached wide-MM auxiliary
  contributions during its existing authoritative traversal for each canonical
  `(t,phi)` and each established prompt/random/dummy/dummy-random source role:
  factor `1.0` before proton cleaning, the existing production proton factor,
  and `1.0 - factor` removed. Ordinary production fills remain the exact
  existing factor-weighted fills; no second scientific tree loop was added.
- The producer composes the pre-proton prompt/random/dummy closures with the
  existing random-window and data/dummy normalizations. Its E.8.2 production
  snapshots are cloned from the actual post-random/post-dummy wide-MM histogram
  on each side of the existing production pruning treatment, never inferred
  from a diagnostic spectrum.
- The pion stage consumes only the exact accepted application payload's wide
  pre-pion spectrum, applied baseline pion template, and post-pion spectrum.
  `zero`, `skip_bin`, rejected, and unavailable child states remain explicit;
  no substitute pion spectrum is reconstructed.
- The private versioned source is
  `e8_2_baseline_stage_source/v1`. It records only the current baseline stages,
  producer-side Lambda-window diagnostic integrals, `Y_0`, the existing
  histogram statistical error, and the existing total yield error. It does not
  change the public yield dictionary or legacy `stage_window_yields`.
- The presentation reader is
  `full_background_subtraction_e8_2/v1`. It fails closed unless the active
  profile is literally `no_empirical_residual`, geometry/closures are valid,
  exact required objects exist, and valid children have authoritative yields.
- Step 3 preserves its existing preliminary PDF/manifest route and retains a
  private detached render-state bundle. Immediately after Step-6 data yields,
  a presentation-only finalizer rebuilds the complete PDF to temporary paths,
  inserts the six E.8.2 pages per canonical parent after E.8.1 context/overlay/
  map pages and before the terminal handoff. After a successful close and
  sidecar write, it creates recoverable copies of the preliminary pair and
  installs the two final names sequentially. A handled partial-installation
  failure restores the preliminary pair; preservation is claimed only when
  that recovery succeeds, and transient render state is released only after
  both final artifacts are installed. This is a pair-safe transaction, not
  impossible literal multi-file filesystem atomicity.

## Frozen boundaries preserved

- Active profile remains `no_empirical_residual`; legacy empirical residual
  Fit 1/Fit 2 remain dormant and are not called, displayed, or represented as
  E.8.2 stages.
- No Method-A or Method-B numerical data, factor, template, or comparison enters
  E.8.2.
- No change was made to production cuts, proton factor calculation, pion `w0`,
  pion components/fits/windows, random/dummy coefficients, canonical binning,
  normalizations, yield formula/error, SIMC, efficiencies, cross sections, the
  F.6.2 artifact, the E.8.1 profile, collector, or wrapper.

## E.8.2.Fix.1 source-review repair

The independent actual-diff review identified only local presentation-sidecar
authority and finalization gaps. This narrow repair leaves the E.8.2 staging and
all scientific ownership unchanged:

- A non-empty complete-renderer failure list now fails finalization before either
  temporary artifact can replace the preliminary PDF/page-manifest pair. The
  preliminary pair and detached Step-3 render state remain available for
  diagnosis/retry.
- The post-yield finalizer now copies both preliminary artifacts to private
  recovery paths before installation. If either replacement fails, it restores
  both preliminary artifacts, reports recovery errors explicitly, and claims
  `preliminary_artifact_preserved` only after successful recovery. It releases
  the transient render state only after both final artifacts are installed.
- A current-chain child with missing early-stage capture, accepted-pion exact
  object/clone authority, required stage object, or final `Y_0`/statistical/total
  authority now fails the entire E.8.2 source closed with a literal authority
  reason. Defined `zero`, `skip_bin`, rejected, and other non-applicable child
  states remain explicit child-level unavailable records.
- The focused suite now directly covers the one-pass event traversal at proton
  factors `1.0` and `0.4`, existing random/dummy/proton algebra, exact pion
  authority, yield authority failures, dormant-fit exclusion, detached state,
  E.8 parent/E.8.2/handoff ordering, renderer-failure handling, second-install
  rollback, recovery-failure truthfulness, and successful pair installation.

## Codex-reported deterministic checks

The following checks were run and reported by Codex. They were `NOT RUN by
ChatGPT`; they are source-review context only and do not establish ROOT/PyROOT,
full-analysis, procedure-PDF farm-rendering, farm, or runtime acceptance.

The following were run locally after implementation:

- `python -B -m unittest testing.test_e8_2_baseline_stage_audit -v` — PASS
  (15 tests after Fix.1): direct one-pass producer capture, normalization and
  closure algebra, exact pion/yield authority, current-chain fail-closed cases,
  dormant-fit exclusion, render-state detachment, ordering, renderer failure,
  rollback/recovery, and successful pair transaction.

## E.8.2.Fix.2 production-pruning boundary repair

Independent review of the cumulative local staging established that the normal
authoritative yield path applies its existing `prune_hist(...)` treatment to
`H_MM_nosub_DATA` after random/dummy subtraction and before baseline pion
application. Fix.2 exposes that frozen production boundary without moving,
calling, tuning, or otherwise changing it:

- `after_proton_pre_prune` is a detached clone of the actual factor-cleaned
  production wide-MM spectrum immediately before the existing
  `prune_hist(H_MM_nosub_DATA, event_threshold)` call. The exact slow-proton
  closure is therefore `after_dummy_pre_proton - proton_component_removed =
  after_proton_pre_prune`.
- `after_proton_post_prune` is a separately detached clone of that same
  production spectrum immediately after the unchanged call and before pion
  application. It is explicitly the production handoff to the pion stage; it
  may differ from the pre-prune snapshot and is never described as proton
  subtraction.
- The exact accepted pion payload remains authoritative. Fix.2 requires its
  `H_MM_nosub_before_pion_subtraction` to match the post-prune snapshot in
  binning and bin content, or fails the current-chain E.8.2 source closed.
  The exact accepted `pion_input - B_pi^0 = after_pion_final` closure remains
  unchanged.
- Valid presentation children now require exactly the eleven current-chain
  stage keys: prompt, random component, after random, dummy component, after
  dummy, proton component removed, post-proton/pre-prune, post-proton/post-
  prune, pion input, baseline pion component, and final baseline clean kaon.
  The reader fails closed on any missing or extra stage, including historical
  Fit 1/Fit 2, empirical-residual, Method-A, or Method-B numerical injection.
- The private producer-side Lambda-window diagnostics now show, in order,
  prompt, after random, after dummy, post-proton/pre-prune,
  post-proton/post-prune, and final baseline clean kaon. They remain detached
  diagnostics; existing public `stage_window_yields`, `Y_0`, its existing
  histogram statistical error, and total error are not changed or recomputed
  by the renderer.
- The existing slow-proton/procedure page retains its placement and count, but
  now has four columns per phi child: after dummy, proton component removed,
  actual pre-prune production, and actual post-prune production/pion-input
  state. Its header states that production `prune_hist` lies between columns
  three and four. Random, dummy, pion, final-MM, page placement, finalizer, and
  terminal handoff behavior are otherwise unchanged.
- Focused deterministic coverage now executes the detached boundary capture
  around a mutating production-prune fixture, verifies source-path ordering and
  unchanged production-prune call inventory, validates pre-prune proton closure
  while allowing the observed pruning transition, validates post-prune/pion
  binning/content handoff failures, rejects forbidden stage injection, checks
  malformed geometry/binning/schema and missing Step-3 state, retains
  producer-owned integrals/statistical error, and preserves all Fix.1 pair-safe
  finalization tests. No E.8.2 helper invokes dormant `bg_fit`.

All Fix.2 changes remain presentation-sidecar observation and validation only.
They do not change event traversal count, production pruning, proton cleaning,
pion `w0`/components/fits/windows, random/dummy normalization, canonical
binning, public yields/errors, Method A, Method B, empirical-residual behavior,
accepted F.6.2 evidence, or farm/runtime state.

## E.8.2.Fix.3 contract-closure regression and wording repair

Fix.3 makes no production or scientific source change. It closes the remaining
independent-review test gap with a narrow deterministic fixture that executes
the actual `calculate_yield_data(...)` public-return path around a stubbed
`bin_data(...)` result, not a duplicate yield calculation:

- The fixture runs the same public path with private E.8.2 source construction
  suppressed and enabled. It proves the returned public `Y_0`, existing public
  yield error, and the sentinel legacy `stage_window_yields` mapping are
  value-identical in both cases.
- With the sidecar enabled, the private source copies that same public `Y_0`
  and total error, while its statistical error is the exact value returned by
  the existing `integral_with_stat_error(final_hist)` call. Building the
  private source and then the presentation payload leaves the returned public
  container and processed-entry `stage_window_yields` object unchanged.
- The dormant-`bg_fit` guard now executes the early capture, exact pion-stage
  capture, source builder, and presentation reader under a raising spy; none
  reaches the dormant historical fit.
- The record now describes the Fix.1 PDF/page-manifest transaction accurately:
  completed temporary artifacts, preliminary recovery copies, sequential final
  installation, restoration on handled partial failure, truthful preservation
  status, and render-state release only after both installations. No literal
  multi-file filesystem atomicity is claimed.

The focused suite passed 23 tests after Fix.3. These local checks do not
establish ROOT/PyROOT, full-analysis, procedure-PDF farm rendering, or runtime
validation.
- `python -B -m unittest testing.test_full_background_subtraction_plots -v` —
  PASS (102 tests; 18 expected PyROOT/superseded-contract skips).
- `python -B -m unittest testing.test_pion_t_bin_particle_stage -v` — PASS
  (7 tests).
- `python -B -m unittest testing.test_pion_hgcer_phase_f_runtime_contract -v`
  — PASS (3 tests).
- `python -B -m unittest testing.test_proton_cleaning_runtime_status -v` —
  one unrelated unchanged fixture error (34 pass): its
  `test_application_coordinate_fix_leaves_committed_lookup_fields_unchanged`
  calls `fingerprint_histogram_content_error(raw_mm)` with `raw_mm=None` in
  `proton_contamination_weights.apply_kaon_proton_cleaning_to_targets`. The
  implicated test, proton source, and ownership helper are unchanged from
  `5b2ba7f...`; an isolated `git archive` source snapshot of that exact HEAD
  reproduces the identical 34-pass/one-error result. E.8.2 did not edit those
  paths; no out-of-scope repair was made.

- `python -B -m py_compile` for the four edited source files and two focused
  test modules — PASS.
- Manifest regeneration — WROTE; manifest check — PASS; memory-health check —
  PASS; bootstrap check — PASS; `git diff --check` — PASS.

After Fix.2, the following additional local deterministic checks passed:

- `python -B -m py_compile src/binning/calculate_yield.py
  src/cuts/full_background_subtraction_plots.py
  testing/test_e8_2_baseline_stage_audit.py` — PASS.
- `python -B -m unittest testing.test_e8_2_baseline_stage_audit -v` — PASS
  (22 focused tests).
- `python -B -m unittest testing.test_full_background_subtraction_plots -v` —
  PASS (102 tests; 18 expected PyROOT/superseded-contract skips).
- `python -B -m unittest testing.test_pion_t_bin_particle_stage -v` — PASS
  (7 tests); `python -B -m unittest
  testing.test_pion_hgcer_phase_f_runtime_contract -v` — PASS (3 tests).
- `python -B -m unittest testing.test_proton_cleaning_runtime_status -v` —
  unchanged out-of-scope 34-pass/one-error result: the same
  `test_application_coordinate_fix_leaves_committed_lookup_fields_unchanged`
  fixture reaches `fingerprint_histogram_content_error(raw_mm)` with
  `raw_mm=None` in unchanged proton source. No repair was made.
- Manifest regeneration/check, memory health, memory bootstrap, and
  `git diff --check` are rerun after the Fix.2 memory update.

Local results do not validate PyROOT rendering, full main execution, the Jefferson
Lab farm, or runtime acceptance.

## Post-review source-status reconciliation

Independent ChatGPT review of the complete cumulative E.8.2/Fix.1/Fix.2/Fix.3
actual diff is `PASS`. The review confirms the accepted baseline path: one
authoritative traversal; true pre-proton capture; random/dummy subtraction;
actual post-proton pre-/post-prune production snapshots; unchanged
`prune_hist(...)`; exact accepted pion input, `B_pi^0`, and `K_0`; and existing
public `Y_0` with its existing total error. The public-output non-regression,
pair-safe finalization, and fail-closed authority behavior remain part of that
reviewed source scope.

The review changes E.8.2 to `SOURCE REVIEWED` only. It does not claim or request
ROOT/PyROOT, full-analysis, procedure-PDF farm-rendering, farm, or runtime
validation. It preserves the active `no_empirical_residual` profile, dormant
Fit 1/Fit 2 status, detached/non-production Method A, numerically absent Method
B, frozen F.6.2 evidence, and all production ownership boundaries. The
Codex-reported deterministic checks above were `NOT RUN by ChatGPT`; the known
34-pass/one-error proton fixture remains pre-existing and out of scope.

## Next

The reviewed E.8.2 source/test/memory set was pushed at
`91bb7809d27d84d7709a6600dbb0dc9ab514a458`; pushed-state review passed. This
remains source review only and makes no E.8.2 ROOT/PyROOT, farm, or runtime
validation claim.

E.8.3 has since completed independent source review. F.6.3 is now the
dependency `NEXT`; the immediate repository action is the user-controlled
commit/push of the reviewed cumulative candidate, then ChatGPT pushed-state
review, followed by the F.6.3 source/runtime-path audit and standalone
implementation contract. Neither E.8.2 nor E.8.3 has ROOT/PyROOT,
full-analysis, farm, or runtime acceptance.
