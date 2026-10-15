# Baseline pion-control Center/lowe neutron-anchor regression — October 14

Status: `ACTIVE`. This bounded read-only investigation concerns the existing
production baseline at Q4p4W2p74 / Center / low epsilon / canonical `t2`,
approximately `0.567 <= |t| <= 0.733 GeV^2`; compare adjacent Center `t1` where
relevant. It starts at the first staged pion-control `B_n=0`, before downstream
transfer amplification. [CURRENT](../CURRENT.md) alone owns the exact NEXT.

## Evidence classes and recording boundary

The [memory-only contract](../phases/kaonlt-baseline-pion-memory-reanchor-2026-10-14.md)
records the user-directed October 13–14 prior-chat priority that could not be
checkpointed when that chat reached capacity. MEMORY ONLY: archived PDF
quantities and the earlier assistant's source comparison below are historical,
attributed evidence supplied through that contract. They are not newly executed
farm gates or independently re-inspected ROOT data/PDF pages in this task.
Current source and the two cited Git changes were read locally at observed
`test` HEAD `d66b215f30d9e43875bdfd6504c4392b49261122`; those limited checks
are SOURCE VERIFIED below. No fit was rerun and no root cause is accepted.

## Reported archived August 26 diagnostic

MEMORY ONLY (prior chat/user-provided PDFs, transcribed from the contract):

| Quantity | Reported Center t2 value |
| --- | --- |
| Staged pion-control `B_n` | `0` |
| Joint-refined `B_n` | `0.0002021` |
| Kaon-side `A_n` | `0.0369527` |
| Neutron amplitude ratio; maximum transfer weight | approximately `182.8`; approximately `183` |
| Staged `(B_n, B_sidis, B_delta)` | approximately `(0, 2.608, 1.692)` |
| Joint-refined `(B_n, B_sidis, B_delta)` | approximately `(0.000202, 1.689, 2.446)` |
| Final kaon-side `(A_n, A_sidis, A_delta)` | approximately `(0.036953, 0.215900, 0)` |
| Fitted background; event-weighted applied background | approximately `0.2529`; approximately `6.556` |
| Event/model integral ratio | approximately `25.9` |

The protected kaon-side delta fit was reported rejected/zeroed. Historic
neighboring Center event/model ratios were approximately `4.32` for `t1` and
`1.84` for `t3`; retain those historical labels rather than silently remapping
them to another canonical-bin convention.

The concurrent setting-wide fit was much better behaved, with integral ratio
around `1.01`. A separate archived 67-page Q4p4 Center setting-wide reference
reported maximum weight `0.0819`, fitted background `0.5967` and event-weighted
background `0.6010`. These differing scopes/PDFs are not an identical-input
regression test, and their agreement establishes no per-t correctness.

## Authority transition and chronology

SOURCE VERIFIED from the local read-only Git inspection:

- `f924a95c55f601546316de7b78185df653f3b0c0` (August 19) added the
  `pion_subtraction_scope="t_bin"` runtime default and retained
  `emit_setting_wide_pion_diagnostic=True` in `src/utility/background_config.py`.
- `ee9e29251865acd1d78b4cbbef1223bfae67565c` (August 20) introduced
  `src/cuts/pion_t_bin_parents.py` for authoritative parent-fit production.

MEMORY ONLY (attributed prior read-only source audit): independently fitted
canonical-t solutions became authoritative while the setting-wide fit remained
diagnostic; per-t fitting already existed in averaging code. Authority changed;
neither correctness of earlier per-t final yields nor the first onset of
instability is proven by that transition. The original staged single-template
solver, staged pass, joint refinement, regularization construction and
scope-template resolver were reported unchanged across it.

MEMORY ONLY (reported chronology): the original anomalous result was August 26,
before the explicit August 27 RF/no-RF tree switch. INFERENCE from that ordering:
the later switch cannot initiate the August 26 pathology and is not the lead
hypothesis. Neither protected kaon-side delta failure nor a final kaon-side
refit is established as the initiating cause of staged `B_n=0`.

## Source findings and bounded hypotheses

SOURCE VERIFIED: `src/cuts/pion_component_fits.py`,
`_solve_nonnegative_template_amplitude`, returns amplitude zero on invalid
template validation, no valid fit bins or invalid/nonpositive template support.
For valid least-squares support it calculates
`weighted_numerator=sum(template*y/sigma^2)` and
`weighted_denominator=sum(template^2/sigma^2)`, then clips their ratio at zero.
The successful path records both moments and `fit_bin_indices`, `fit_y`,
`fit_sigma`, `fit_template`. Early failure can lack these diagnostic arrays;
absence alone must not be interpreted as signed cancellation.

MEMORY ONLY (prior source audit): the neutron is fitted first from the signed
pion-control histogram in the neutron anchor around `0.88–0.94 GeV` with this
nonnegative least-squares estimator. NOT VERIFIED: which zero-producing branch
was taken in the original Center t2 fit. A signed weighted numerator explanation
requires the matching original inputs and solve result, not the final amplitude
alone.

SOURCE VERIFIED: `_build_regularization_width_map` uses a default `1e-3`
reference-amplitude floor, scaled by the effective prior scale; current
`background_config.py` also sets `joint_refinement_amplitude_floor=1e-3`.
This is a regularization-width construction, not a fitted-amplitude lower bound.
MEMORY ONLY (prior audit): joint refinement is centered on the staged zero.
INFERENCE: such a prior can discourage recovery from that initialization;
the original effective prior and its actual effect remain NOT VERIFIED.

SOURCE VERIFIED: current baseline config leaves maximum full-range chi2/ndf
and oversubtraction count/fraction limits unset (`None`). MEMORY ONLY (prior
audit): fitting/application gates can admit weak decompositions, large ratios
and event/model discrepancies. This source-level quality gap neither establishes
the original solve's cause nor validates an alternate physical model.

SOURCE VERIFIED: `src/utility/pion_component_shapes.py` independently
unit-normalizes binned `(t,phi)` templates. In
`src/cuts/pion_component_fits.py`, `_resolve_component_scope_hist` sums the phi
shapes for a parent through `_sum_hist_list_to_unit_area`, then renormalizes.
INFERENCE / testable input-shape hypothesis: populated unit-area phi slices may
receive equal mixture weight instead of their original event-weighted mixture.
Actual populated-phi inputs, normalization integrals and fallback selection
must be identified before testing its relevance. This is not an accepted cause
or an approved template change.

Signed cancellation, wrong normalization, alignment, SIMC template construction
and a particular commit remain candidate mechanisms, not findings of root cause.

## First-failure evidence boundary

The future read-only audit first searches existing applicable artifacts and
source-backed diagnostics for the original Center/lowe t2 solve, comparing t1
where relevant. Required identities include cuts, fit windows, histogram/bin
and coordinate conventions, signed-source coefficients, run/source lineage and
SIMC component/template construction. Recover:

- Prompt, random, dummy and dummy-random neutron-anchor contributions, retaining
  signed combination and uncertainties.
- `fit_bin_indices`, `fit_y`, `fit_sigma`, `fit_template`,
  `weighted_numerator`, `weighted_denominator`, and `success/message`.
- Staged initialization, effective joint prior/scale/width and populated-phi
  SIMC inputs, pre-normalization integrals and any support fallback.

Determine whether zero came from a nonpositive weighted numerator or
absent/invalid template/fit-bin support. If existing evidence cannot decide,
report the missing inputs and define a separate narrow diagnostic-only source
contract. This memory task authorizes no diagnostic implementation, farm
command, production correction, fit redesign or blocked-owner retry.

## Paired E.8 crosswalk and preserved gates

MEMORY ONLY (attributed prior chat): later Left/lowe + Center/lowe paired progress
was reported after tracked CURRENT was last published. No independent paired
closure is present in the relevant tracked evidence inspected here. Retain the
paired extension's recorded `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`
status; do not manufacture receipts, hashes or 16-page visual acceptance. Its
former farm NEXT is superseded as the immediate user priority, not failed.
The [October 12 evidence and October 13 addendum](../evidence/e8-current-snapshot-contribution-runtime-closure-2026-10-12.md)
retain their exact detached Left-only numerical/presentation scopes.

Accepted E.8.2 Left/lowe and [baseline-five numerical/procedure gates](../evidence/e8-2-canonical-five-baseline-yield-runtime-closure-2026-10-09.md)
stay `CLOSED / RUNTIME VALIDATED` within their accepted scopes. The latter's
135/135 data and SIMC checks establish producer/support consistency, not all
Center t2 component-decomposition physics. E.8.3, F.6.3/E.8.4 and all other
accepted closures are preserved without reinterpreting their evidence.
Baseline fit/transfer validity is a prerequisite to further interpreted full
E.8 production-impact claims and any F.6.4 promotion decision, not a new phase.

Private canonical-five F.4 remains `BLOCKED` at
`f6_3_scientific_equivalence_mismatch:f4:$.parents[0].absolute_baseline_parent_sum`.
Absolute-SIMC unit proof remains separately `BLOCKED` at
`SIMC_normfac_luminosity_and_charge_units_not_source_proven`; final E.8/F.6.4
remain `BLOCKED`. Historical F.1/F.4 forensics remain `DEFERRED`, a different
question from this baseline investigation. Method A stays detached/non-production;
Method B stays independent diagnostic/cross-check only and numerically excluded.
All random/dummy/slow-proton/pion handling, HGCer roles, SIMC inputs, cuts,
models, priors, templates, numerical ratios, weights, errors, yields, L/T,
cross sections, provenance and `no_empirical_residual` remain unchanged.

## Read-only-first audit addendum — October 14

### Temporary interruption and return condition

`ACTIVE`: this is an explicit **TEMPORARY interruption** of the normal
[roadmap](../roadmap/STATUS.md), not a new phase or replacement program.
The [read-only-first contract](../phases/baseline-pion-regression-readonly-first-task-contract.md)
owns this local audit. Once the discrepancy is explained and the existing
subtraction validated, or a necessary minimum correction is independently
reviewed, user-pushed and user-farm validated, return mandatorily to the
interrupted E.8 roadmap position. Preserve every closed, active, deferred,
blocked and pending item and existing scientific role. CURRENT alone owns NEXT.
The farm-owner output-isolation blocker remains real, but applies only to a
future operation requiring the normal launcher. The separate isolation proposal
is preserved; neither isolation engineering nor owner implementation proceeds
under this task. No source, fit, configuration, physics or acceptance changed.

### Local identity, available evidence and comparison limits

SOURCE VERIFIED: this audit inspected local `test` HEAD
`d860eeafbb889bcc04cfddb1ff6019cfcff50922`. Following the user's correction,
continuation uses local files and local Git objects only; remote checks are
outside the continued task. No farm, ROOT, fit or procedure rendering ran.
The existing unrelated `.M` entries have equal HEAD/index/Git-filtered-worktree
blob identities; their raw bytes and existing review bundles were preserved.

The local workspace search (including ignored files) found no ROOT file, log,
or applicable Center/full-analysis/pion procedure PDF containing the original
solve. Local Git history of ROOT, FullAnalysis/rand_sub PDF and log paths yielded
only unrelated historical log filenames, not a matching numerical record.
Exact quoted historical integrals occur in the existing investigation and the
two baseline memory/read-only contracts, not an independently available raw
fit artifact. No external artifact directory was supplied or guessed.

| Comparison | Supported identity and scope | Numerical conclusion |
| --- | --- | --- |
| Current source audit | Exact local HEAD above; configured Q4p4W2p74 runtime defaults and explicit canonical-t parent path | SOURCE VERIFIED code behavior only; no current numerical fit/application |
| Reported August 26 anomaly | MEMORY ONLY: Center/lowe, historical canonical t2 approximately 0.567–0.733 GeV^2, phi-integrated parent | Staged zero precedes tiny refined denominator and large application; original SHA, run lineage, exact cuts, normalization, binning and resolved configuration NOT VERIFIED |
| Adjacent Center t1 | MEMORY ONLY: reported event/model ratio approximately 4.32 | Comparison control, not a demonstrated good same-scope predecessor; exact inputs NOT VERIFIED |
| Archived 67-page Center/lowe reference | MEMORY ONLY: setting-wide, 0.5967 model versus 0.6010 event-weighted | Different fit scope; no proof of earlier per-t correctness or historical working SHA |

Thus a last demonstrably good / first bad **same-scope numerical pair is NOT
VERIFIED**. The earliest reported internal anomaly is staged `B_n=0`; the
earliest independently verified historical numerical divergence cannot be
located from these local records. A code authority transition is not that datum.

### Actual source/history comparison

SOURCE VERIFIED from actual diffs, not commit subjects:

- Immediately before `f924a95c55f601546316de7b78185df653f3b0c0` is
  `fe967c2eb08c25830a96525deee012a5c4e57fa8`. Per-t fitting/application
  already existed in `src/binning/ave_per_bin.py::process_hist_data`.
  `f924a95` adds the `t_bin` runtime default, makes the setting-wide result
  diagnostic-only in `rand_sub`, and builds RF-restored authoritative parents
  in `src/main.py` before child yield/average consumers. The averages path
  reuses those parent fit results instead of refitting them at that point.
  This is the earliest authority/consumer difference in the requested pair
  of transitions, not proof that it initiated the numerical pathology.
- Immediately before `ee9e29251865acd1d78b4cbbef1223bfae67565c` is
  `5ad24e044f76aff933f83a3b59ea4d2b12033910`. `ee9e292` moves normal parent
  construction/freezing to `src/cuts/pion_t_bin_parents.py` during `rand_sub`,
  still calling the shared integrated input builder. It adds an explicit
  half-open t-edge filter to the detached parent event-application diagnostic
  and makes downstream consumers validate the frozen collection. These are
  real producer/selection changes; no original numerical inputs establish an
  erroneous boundary choice or a first-bad result there.
- AST comparisons across both commits and their immediate predecessors find
  identical `_solve_nonnegative_template_amplitude`,
  `_run_staged_component_pass`, `_fit_staged_anchor_templates`,
  `_run_joint_template_refinement`, `_apply_joint_component_refinement`,
  `_build_regularization_width_map`, `_resolve_component_scope_hist` and
  `_resolve_single_scope_hist` in `pion_component_fits.py`, and identical
  `build_simc_shape_pion_control_weights` in `pion_component_subtraction.py`.
  `pion_component_shapes.py` has no diff across this bounded transition range.
  This excludes a change to those routines in that range, not changes to
  their callers, inputs, configuration or later history.

Current lineage was independently traced at the stated HEAD:
`pion_t_bin_parents.py:561–626` consumes producer-owned kaon and signed pion
inputs, resolves t-specific shapes/alignment and calls the component fitter.
It is no longer the original August shared-averages construction path.
`rand_sub.py:1261–1326` builds the parent event diagnostic from those frozen
records. Historical input equivalence to that current path is NOT VERIFIED.

### Model versus event application

SOURCE VERIFIED: `pion_component_subtraction.py:745–894` constructs
`P_model=sum B_i*T_i` and `K_model=sum A_i*T_i` over the active pion background
components, excluding signal templates. Final/staged amplitude dictionaries
and protected kaon applied amplitudes determine the selected model. In normal
supported, finite, unclipped bins it uses
`w(M)=K_model(M)/P_model(M)`, hence `w*P_model=K_model`.

The event product is instead `w(M)*P_signed_data(M)` under the same MM bins
and selection. `rand_sub.py:1217–1258` multiplies each authoritative record's
signed coefficient by the bin weight; `pion_component_subtraction.py:977–1010`
preserves prompt/random/dummy/dummy-random coefficients
`(+norm_data, -norm_data/nWindows, -norm_dummy, +norm_dummy/nWindows)`.
Its cached-event path (`1078–1094`, `1114–1187`) retains cached coefficients
and scope masks; baseline Method-A multiplier defaults to one. Source/mask,
MM-coordinate and optional downstream treatment identities must match before
comparing any particular histogram or integral.

For those matched, supported/unclipped bins the model/event difference is
`w(M)*(P_signed_data(M)-P_model(M))`. A small positive fitted denominator can
amplify a control fit/data mismatch even when the model-ratio closure passes.
This is a SOURCE VERIFIED mechanism, not a demonstrated historical root cause.
Denominator at/below the configured floor produces zero weight; nonfinite
weights become zero; clip bounds apply afterward. Flow bins and unsupported
or clipped bins do not inherit `w*P_model=K_model`. An amplitude ratio alone
does not determine the full multi-component weight or its integral.

### Earliest staged solve and template checks

SOURCE VERIFIED at this HEAD: `background_config.py:357–434` fits control
`pi_n` first, one staged pass, neutron anchor 0.88–0.94 GeV with configured
data-MM offset. Only `k_sigma0_signal` has the explicit `window_integral`
mode/window; neutron resolves to `least_squares` through
`pion_component_fits.py:3925–3998`. Initial other amplitudes are zero, so the
first neutron solve sees the signed control target, before joint refinement
or the later protected kaon delta fit. Actual historical overrides and shifted
window/bin identities are NOT VERIFIED.

`pion_component_fits.py:620–840` selects bins by center, windows and finite
target content; invalid/nonpositive errors use `max(sqrt(abs(y)),1)`.
The LS solution is `max(sum(T*y/sigma^2)/sum(T^2/sigma^2),0)` only after valid
template/bin/support checks. Missing/mismatched/nonpositive-integral template,
no valid bins or invalid/nonpositive denominator can also return zero with a
failure message. The alternative window-integral mode uses a different target/
template sum estimator. Therefore the original zero cannot be assigned to a
nonpositive LS numerator without its estimator, arrays and success/message.
The `1e-3` reference floor (`3222–3238`) constructs regularization widths;
it does not impose a lower fitted amplitude or explain the observed zero.

SOURCE VERIFIED: `pion_component_shapes.py:657–722` fills SIMC with
`iter_weight`, `-evt.t` and phi converted to degrees, using SIMC cuts and the
kaon hole exclusion, then independently normalizes setting/full and each
`(t,phi)` shape. `pion_component_fits.py:5893–5985` sums the normalized phi
shapes for a parent and renormalizes, falling back to setting-wide shape if
t support is unusable. Consequently this aggregation does not retain raw phi
population weights. Its relevance to the historical t2 solve is INFERENCE:
original populated slices, raw integrals, selected fallback and aligned/scaled
template are absent. No template reweighting, normalization or fit change is
authorized or validated by this observation.

### Findings matrix and one decisive missing artifact

| Evidence label | Finding |
| --- | --- |
| SOURCE VERIFIED | Authority/ownership changes and unchanged bounded solver routines; separate fitted-model and signed-event products; possible denominator amplification; multiple zero-producing solver paths; normalized-phi aggregation |
| MEMORY ONLY | August 26 staged zero, refined approximately 0.0002021, A_n approximately 0.0369527, maximum weight approximately 183 and 0.2529 versus 6.556; setting-wide and t1 comparisons |
| INFERENCE | Signed cancellation, template mixture/support, prior-constrained recovery or an amplified control-model mismatch may matter; none is an accepted cause |
| NOT VERIFIED | Same-scope working predecessor, first-bad numerical identity, original estimator/solver branch and input decomposition, historical root cause, current farm reproduction |

**One earliest missing decisive artifact:** the original August 26
Q4p4W2p74 / Center / lowe / t2 **first neutron-anchor solve record**, including
source/run/configuration and coordinate/selection identity, prompt/random/
dummy/dummy-random contributions and errors, selected SIMC shape/support,
fit bins/y/sigma/template, estimator, numerator/denominator (or integral-mode
sums), and success/message. This one record distinguishes clipped LS from
invalid/absent support before any downstream explanation. A matched good
predecessor remains necessary to date a regression, but is not substituted
for this earlier unresolved solve. Stay `ACTIVE`; no diagnostic implementation
or farm run follows from the evidence gap.

Local integrity note: the initial manifest mismatch consisted solely of the
two pre-existing untracked baseline contracts, with no changed listed-file
hashes. The ordinary generator inventories both tracked and nonignored
untracked memory files; regeneration therefore also records the preserved
isolation proposal. Its presence in integrity metadata is not execution or
publication authorization. Actual-diff review must account for that untracked
dependency before using the contract's publication block; do not silently
publish a manifest referencing an absent file. The existing CURRENT soft-size
warning is nonblocking and scheduled for the next meaningful memory checkpoint;
no accepted-status compaction is part of this task.
