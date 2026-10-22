# Full baseline pion-SIMC model-to-application discrepancy — October 14–15

Status: `ACTIVE`. This temporary roadmap interruption concerns the entire
existing baseline SIMC-derived pion-background construction, final model and
signed event application. Initial focus: Q4p4W2p74 / Center / lowe / canonical
`t2`, approximately `0.567 <= |t| <= 0.733 GeV^2`; adjacent t1/t3 and setting-wide
results remain scoped controls. The candidate staged `B_n=0` does not define
the scientific objective or an established first historical defect.
[CURRENT](../CURRENT.md) alone owns the exact NEXT; the historical filename
is retained for provenance/link stability.

## Scientific-memory correction / re-anchor — October 15

The [October 15 documentation contract](../phases/baseline-pion-full-model-memory-reanchor-2026-10-15-task-contract.md)
restores the user's full-model versus full-subtraction question. The recovered
Library excerpts from `Pasted text(9).txt` and `Pasted text(10).txt` are supplied
historical chat evidence, not locally available raw PDFs or new farm receipts.
At local `test` HEAD `06b4573e6879fd2cd619f179a703afb85d1e6aae`, this task changes
memory only; it neither reruns nor independently validates those artifacts.

The user's recovered words were:

> "So the issue is with baseline. Something fundamentally wrong has occured here. This is from Center_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction.pdf"
>
> "The second plot is from Center_kaon_rand_sub_Q4p4W2p74_lowe.pdf  this looks correct, but something breaks between the two. We have a good background model. I fear this may be happening across the board. We need a detailed audit"

The recovered procedure description was:

> "The way, it should, work is we have 3 pion channels (neutron, sidis, delta). Each simc model is used to build the background, in particular, pi-n and pi-sidis are used. Then using this two channel model it should be improved with the pi-delta look. Then this initial pion model should be checked and adjusted to account for decomposition."

The earlier Q4p4 Center/lowe reference was identified in that discussion as
`Center_kaon_rand_sub_Q4p4W2p74_lowe(8).pdf`. The user said:

> "See it use to work fine, something broke the baseline at some point"

Preserve these observations even without a last-good SHA. They do not prove
that the exact same canonical-t producer, cuts, population or normalization
previously passed. The concern about broader failure is a user-reported lead,
not a verified across-kinematics defect. The original component-model versus
full-background page contrast requires its own scope/run/coordinate crosswalk.

### Historical PDF evidence matrix

All quantities/page identifiers below are MEMORY ONLY, recovered historical
supplied-PDF descriptions. No raw PDF was newly rendered or farm validated.

| Reference | Fitted/model | Applied/event-weighted | Maximum weight; reported scope/pages |
| --- | ---: | ---: | --- |
| Q3p0W3p14 Center/lowe well-behaved procedure | ~0.3164 | ~0.3190 | ~0.0381; setting-wide, pp. 54–57; staged pp. 42–44, joint p. 35 |
| Earlier Q4p4W2p74 Center/lowe 67-page reference | ~0.5967 | ~0.6010 | ~0.0819; setting-wide, pp. 54–59 in recovered discussion |
| Later Q4p4W2p74 Center/lowe problematic parent | ~0.2529 | ~6.556 | ~183; canonical t2, later 115-page procedure, component pp. 78–91 |

The setting-wide and canonical-t rows are not an identical-input regression
pair. Adjacent Center t1/t3 event/model ratios ~4.32/~1.84, staged `B_n=0`,
joint `B_n≈0.0002021`, `A_n≈0.0369527`, rejected protected-delta candidate and
final `A_delta=0` remain historical reports/candidates. Recovered chat also
describes abrupt t1/t2 final-MM features and events in bins with zero inherited
`w0`, distinguishing stored-weight discontinuities from display-scale effects.
Their physical interpretation and cross-setting extent remain NOT VERIFIED.

### Full-chain priority and superseded narrowing

The progressive narrowing to a neutron-anchor regression was the memory error.
`pi_n` is one component; its fit window is an implementation detail. The
historical first-zero exploration, source/history findings, numeric reports
and uncertainties below remain preserved as candidate-level audit history.
In particular, **the conclusion that the first neutron solve is the single
decisive artifact is superseded as the sole NEXT**, not proved false as a
component observation. Neither that solve nor protected-delta failure is an
accepted root cause or authority for a new global refit. The later historical
instruction to define a diagnostic contract is not a current authorization.

The full read-only, PDF-first comparison establishes scope, run/input/cut,
normalization, MM-coordinate/bin and provenance identities, then traces:
all three SIMC templates (`pi_n`, `pi_sidis`, `pi_delta`) -> staged and any
joint/protected fits -> final pion-control and kaon-side decompositions ->
composite model and `w0(MM)` -> signed pion-control event application -> actual
applied background -> resulting full-subtracted kaon MM/yields. Keep protected
kaon-signal treatment separate. Locate the earliest supported discrepancy;
if original PDFs or numerical records are unavailable after searching accessible
sources, identify the single earliest decisive missing evidence for that chain.
This memory correction does not perform that subsequent scientific audit.

| Evidence label | Retained meaning |
| --- | --- |
| SOURCE VERIFIED (prior bounded review, not newly rerun here) | Multi-component configuration and distinct control/kaon models; supported finite unclipped `w0=K_model/P_model` versus applied `w0*P_signed_data`; zero-support/clipping behavior, potentially unset quality limits; protected kaon delta fallback can leave `A_delta=0` without final readjustment of the first two kaon amplitudes |
| MEMORY ONLY | User's PDF contrast, previously working observations, table/page reports, abrupt MM/zero-inherited-weight descriptions |
| INFERENCE | In common supported unclipped bins, `w0*(P_signed_data-P_model)` is a possible difference mechanism; small denominators, template mixture, decomposition or fallback may matter; no culprit accepted |
| NOT VERIFIED | Matched first-good/first-bad identity, original initiating defect, root cause, across-kinematics failure, current reproduction or physical correctness |

Distinct zero-producing solver paths prohibit inferring cancellation from
`B_n=0`. August 19/20 authority changes remain context, not a demonstrated
initiator; by reported chronology, the August 27 RF-tree change cannot initiate
the August 26 anomaly. Supported-bin identities do not extend universally to
masked, flow, unsupported or clipped bins.

This remains a **TEMPORARY** interruption, not a phase. Prospective Method-A
development/promotion and normal E.8/F.6 forward implementation/interpretation
are paused. Accepted completed scopes, pending local candidates, independent
F.4/absolute-SIMC/final-E.8/F.6.4 blockers and deferred F.1/F.4 archaeology remain
unchanged. The component diagnostic stays
`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`, without verified capture;
farm owner stays `BLOCKED` on isolation,
and the tracked isolation proposal remains preserved, not launched. No code,
diagnostic widening, plotting repair, farm run or new physical hypothesis is
authorized. After explanation and validation of the existing subtraction, or
independently reviewed, user-pushed and user-farm-validated minimum correction,
return mandatorily to the approved roadmap at its interrupted position.

## Historical October 14 exploration (retained; narrower priority superseded)

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

## Low-epsilon farm evidence and scientific checkpoint — October 17–18

`ACTIVE`: temporary baseline-pion investigation; CURRENT alone owns NEXT.
The [October 18 memory-only contract](../phases/kaonlt-low-epsilon-baseline-science-memory-checkpoint-2026-10-18-task-contract.md)
carries the supplied farm artifacts and prior ChatGPT scientific review.
The RUNTIME VERIFIED findings below are attributed to that supplied evidence
and review, not a local rerun, raw-artifact revalidation or new farm closure
by this checkpoint. Historical sections above retain their original dates
and uncertainty labels; the new capture supersedes their no-current-capture
statements only for the explicit scope below.

SOURCE VERIFIED locally at `test` HEAD
`f1d1984b5e0fca9ba3c5666477996af783423db7`: the published launcher `-d` clears
the legacy Left-only selector on its single full low-epsilon child, selects
ordinary available settings, and exits before full high epsilon. Plain `-d`
does not invoke the `-c` merge block. This debug selection change did not
alter analysis physics. The supplied checkpoint records GitHub `test` at that
same HEAD and explicit verification of farm checkout HEAD before execution;
no live remote query was performed by this memory-only task.

RUNTIME VERIFIED (supplied farm invocation/review):
`./run_Prod_Analysis.sh -d 4p4 2p74` ended with `DEBUG_EXIT_STATUS=0`, reporting
all available low-epsilon settings complete and full high epsilon skipped.
Center and Left lowe produced analysis and full-subtraction PDFs; no populated
Right/lowe result is established. One Center/lowe canonical t1/t2/t3 seven-stage
JSON/PDF pair was captured; the PDF has 19 pages. Producer status remains
`captured_pending_independent_validation`, completeness `captured`, PDF
`available`, and `observed_source_HEAD=not_provided`. Verified checkout identity
does not repair that missing producer identity or imply full independent
provenance/checker acceptance or broader production scientific validation.
The reviewed `KaonLT_loweps_20261018_093137.zip` has eight intact members: six
Center/Left analysis, full-subtraction and pion-fit-debug PDFs plus the Center
JSON/PDF pair. The desired future handoff is one ZIP per farm run, preserving
original evidence and checksums; no new hashes or receipts are invented here.

RUNTIME VERIFIED (supplied scientific review), retaining setting/t scope:

- Center/t1: staged SIDIS normalization describes its local anchor well, but
  joint full-range refinement heavily suppresses SIDIS and underpredicts pion
  control around Lambda. The staged fit itself overpredicts portions of the
  high-MM tail. This is transition-specific fit/model tension, not proof of
  a broken optimizer or acceptance of the staged fit over the entire range.
- Center/t2: one high-leverage negative signed random contribution near
  0.910–0.916 GeV contributes to a negative neutron least-squares numerator
  (~-61.674). Initial staged `B_n=0`, refined `B_n~0.000202`, and the ratio
  reaches `w0~182.82`. The `K_model` pion integral is ~0.2529 versus event-applied
  background ~6.5562. Model algebra and after-subtraction arithmetic close,
  while signed-event/model agreement fails in each Center t bin. Do not drop
  the negative event or identify it as the sole root cause. These full-range
  integrals are neither Lambda cross sections nor extracted yields.
- Left: setting-wide neutron amplitude stays positive, but fit shape quality
  is poor and production canonical-t plots show concerns. Setting-wide model
  closure is not canonical-t signed-event/model closure. Matching source-resolved
  Left canonical-t diagnostics are NOT VERIFIED.

SOURCE VERIFIED: `src/main.py` calls `load_setting_pion_component_shapes`
without `t_bins` or `phi_bins`; their loader defaults are `None`, and the scope
resolver can fall back to setting-wide component shapes. RUNTIME VERIFIED
(supplied Center artifact): no constructed t-binned SIMC shapes are reported.
This is a possible limitation, not a proven cause or permission to activate
new shapes. It does not validate the historical normalized-phi mixture
hypothesis for this run.

INFERENCE: statistical-support-gated fits per t followed by conditional
parent-anchored `(t,phi)` refinement could be useful. No such physics change
is accepted or implemented. First inspect existing Center/Left lowe evidence
for signed/absolute counts, Sumw2, cancellations, meaningful effective support,
parent-child event covariance, per-child SIMC occupancy, parent fit quality,
component support and under-Lambda residuals; identify actual missing evidence.
The audit tests justification, without implementing refinement or rerunning
the farm under this checkpoint.

Preserve the initially useful staged `pi_n`/`pi_sidis`/`pi_delta` model,
normalization, signed events, weights, cuts, priors, templates and all frozen
analysis source. No first-bad historical commit or sole cause is established.
Method A remains detached/non-production; Method B remains independent,
diagnostic/cross-check only and numerically excluded. Earlier accepted phase
closures and their original scopes, separate private canonical-five/F.4 and
absolute-SIMC blockers, and historical normal-launcher/attempt-isolation
blockers remain unchanged. Successful plain `-d` is distinct from those
broader blocked operations. After this checkpoint's independent actual-diff
review, user publication and pushed-state synchronization, return directly
to CURRENT's read-only support audit; no new phase, source contract or
farm-readiness project follows from this documentation update.

## Global-to-fine-bin procedure audit and user correction — October 19

`ACTIVE`: memory-only checkpoint under the
[October 19 contract](../phases/kaonlt-pion-hierarchy-audit-memory-checkpoint-2026-10-19-task-contract.md).
CURRENT alone owns NEXT. Source was inspected at local `test` HEAD
`1975c3d8656253e262ff7e18306a93ca800b6027`; no scientific source, fit,
configuration, farm run or accepted status changed. Earlier October 17–18
findings remain bounded diagnostics, not canonical-t model acceptance.

### PDF location and observed onset

MEMORY ONLY (user correction/supplied display comparison): the full historical
setting-wide component procedure and diagnostic plots are still present in
`Center_kaon_rand_sub_Q4p4W2p74_lowe_pion-fit-debug.pdf`. The earlier reference
had 67 pages and October's main rand-sub PDF had 30; this presentation routing
change does not show that global fit machinery vanished. SOURCE VERIFIED:
`src/cuts/rand_sub.py:6733,6980–7015` routes setting-wide component pages to
`pdf_destinations["pion_fit_debug"]`.

The user finds the setting-wide model substantially more sensible and places
visible deterioration during canonical-t construction/refinement, before final
pion subtraction. This corrects the earlier framing of final subtraction as
the apparent first failure. It is a visual observation, not independently
accepted global model quality, an exact first-bad operation or a historical
culprit. No raw PDF was newly inspected or rendered by this checkpoint.

### Live source audit

SOURCE VERIFIED at the HEAD above:

- `src/main.py:789–794` calls `load_setting_pion_component_shapes()` without
  `t_bins/phi_bins`. `src/utility/pion_component_shapes.py:313–330,679–690,801–806`
  constructs/fills binned shapes only with supplied edges. The component scope
  resolver (`src/cuts/pion_component_fits.py:5951–5987`) falls back to
  `setting_shape_full` when t support is unavailable. Separately, RUNTIME
  VERIFIED (supplied October Center artifact/prior review, not rechecked here):
  no constructed t-binned SIMC shapes were reported. This neither proves
  t-dependent SIMC defective nor authorizes activating new shapes.
- In active `t_bin` scope, `src/cuts/rand_sub.py:5149–5167` marks the setting-wide
  fit `diagnostic_only`, clears its subtraction payload and explicitly excludes
  it from later weight/spectrum authority. `src/cuts/pion_t_bin_parents.py:566–636`
  calls the same full `build_particle_subtraction_component_result()` with
  t-restricted pion-control and post-proton kaon inputs. It inherits parent
  alignment, without global fitted amplitudes/covariance as hierarchical
  constraints. Frozen parent identity and producer-fingerprint checks
  (`pion_t_bin_parents.py:503–534`) remain real and required.
- `src/cuts/pion_component_fits.py:3611–3628,3728–3740` passes that fit's local
  staged amplitudes to joint refinement. When enabled, regularization centers
  on those coefficients (`3442–3474`), not a transferred global expectation.
  Local-stage validation and `_promote_sparse_zero_solution_if_applicable`
  (`427–501,7333–7357`) are separate from parent-conditioned model adequacy.
- `src/binning/calculate_yield.py:2437–2466` assigns
  `component_fit_results[j][k] = scope_result` from the authoritative t parent
  and applies its accepted event weight to local `(t,phi)` control events;
  it does not run a complete local component refit. Frozen-parent application
  acceptance/fallback policy is retained. The separate free local-fit branch
  (`2510–2571`) is outside authoritative `t_bin` scope; merely enabling it
  unchanged is not the requested hierarchy or an approved correction.

### Scientific intent and evidence limits

User-requested future scientific intent: preserve the complete setting-wide
staged `pi_n`/`pi_delta`/`pi_sidis` procedure at t and, later, `(t,phi)`.
The global solution should anchor t without suppressing genuine potentially
large t dependence. Each child refinement should be anchored to its own t
parent and conditional on meaningful local support and covariance. Shared-event
parent/child correlations preclude treating the parent as independent Gaussian
evidence. No fixed percentage clamp, prior width, new covariance model or
automatic fallback algorithm is authorized or validated by this checkpoint.

The October 19 SIMC support/shape report establishes only its reported local
support/shape information, not model acceptance or a first-bad commit (MEMORY
ONLY here; report not independently revalidated). High parent SIMC totals,
local tests, frozen-source consistency and accepted model algebra/subtraction
arithmetic do not establish canonical-t scientific fit goodness. INFERENCE:
local selection/normalization, global-shape fallback versus candidate genuine
t-shapes, and unconstrained staged/joint amplitudes need comparison against
existing October 17–19 evidence before choosing an intervention. Keep global
and t acceptance/normalization conventions separate. Historical culprit,
changed-model root cause, hierarchy efficacy and full farm acceptance remain
NOT VERIFIED.

After actual-diff review, user publication and pushed-state synchronization,
proceed directly to CURRENT's read-only global-versus-t comparison (Center lowe
t1/t2/t3, Left lowe control), establishing what changes before subtraction and
the smallest directly evidenced intervention. Any warranted t-level correction
requires a separate scientific implementation contract; only after t-level
scientific fit agreement and farm validation consider conditional full child
refinement. No source edit, child-fitter promotion, final-subtraction/yield
change, farm command, owner/pathing cycle or new phase follows here. Preserve
all prior accepted scopes, independent blockers, frozen physics and Method-A/
Method-B roles, and return to the interrupted roadmap after baseline resolution.

## Controlled canonical-t template comparison implementation — October 19

`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING` under the
[single comparison contract](../phases/kaonlt-pion-t-template-controlled-comparison-2026-10-19-task-contract.md).
SOURCE VERIFIED against starting `test`/`origin/test`
`4007d8a43a604628e5b74aca54229d3e6bd48848` and the local proposed diff:
the prior read-only source comparison confirms the ordinary loader lacks t
edges and the t resolver can inherit global shapes. Summing its independently
normalized phi children would give an equal-phi mixture. This completes that
source audit, without independent numerical revalidation of October farm inputs
or a claim that template choice caused the observed deterioration.

The additive cuts-stage observer uses only
`LT_ANALYSIS_PION_T_TEMPLATE_COMPARISON=1`, inherited by the unchanged launcher,
and limits output to Q4p4W2p74 Center/Left lowe, three canonical t intervals.
Baseline A remains the frozen parent by reference. Candidate B uses cloned
original pion-control/post-proton kaon histograms, the same protected templates,
configuration and existing full fitter, and the exact retained parent alignment
map/windows. No alignment rescan/persistence, transferred global-amplitude prior,
production weight, child refit or final subtraction is introduced.

The detached raw loader reads original `iter_weight` (the same unit default
when absent), raw `missmass`, `-t` and `phipq*180/pi`; it reuses the active
SIMC cuts and existing HGCer cut without resetting selection globals or creating
a replacement named cut. It accumulates raw canonical-t parents and raw phi
children independently, with counts, signed/absolute sums, Sumw2 and flows.
Weighted child/parent and selected/accepted/excluded closure precede one parent
normalization. The raw selected stream must also reproduce the already-loaded
global template contents and normalized Sumw2, preventing source/selection drift
from confounding the t-shape comparison. Empty/nonpositive/negative-bin templates and normal-bin effective
support <= 1 are unavailable, not global fallbacks. That support check is a
minimal estimability guard, not a model-adequacy threshold or physics acceptance.
Source SHA-256/tree, cut state, hole points, phi coverage, raw/normalized arrays
and differences from global/baseline-selected shapes are retained.

One separate JSON and six-page diagnostic PDF per setting retain all parent
fit/solver diagnostics, staged/refined coefficients and uncertainties as provided
by the existing fitter, histogram fingerprints and copied data/model/component
arrays. The PDF compares staged/refined A/B with common data and display-only
residuals. Original data/dummy file identities and relevant source-file hashes
are recorded. Parent identity/content/error checks, candidate input/configuration/
alignment checks, detached clones and restoration of the ownership-debug ledger,
ROOT directory and drawing pad protect existing consumers/output routing.
Candidate alignment metadata deliberately retains baseline provenance: it is
reused fixed alignment, not a newly calibrated local map. Missing alignment,
source/support or candidate fit produces explicit unavailable evidence.

Local deterministic validation: 122 tests ran across the new comparison and
existing parent-stage, chain/anchor audit, parent-integrity, launcher and dynamic
alignment suites; 111 passed, 11 existing ROOT-dependent checks skipped.
Scalar/fake-hist/AST and mocked PDF-command checks do not establish ROOT/PyROOT,
real fitter integration, rendered page correctness, farm behavior or scientific
acceptance. Windows Python with existing NumPy was used; WSL Python lacked
NumPy. Runtime comparison JSON/PDF and any improvement/cause remain NOT VERIFIED.

CURRENT retains the sole NEXT: after actual-diff acceptance, user publication
and pushed-state synchronization, one controlled farm A/B comparison, without
providing a farm command here. The existing CURRENT soft-size warning is
nonblocking and scheduled for consolidation at the next memory checkpoint or
milestone audit. Preserve all accepted scopes, independent blockers, ordinary
PDFs, fits, weights, cuts, priors, subtraction and Method-A/Method-B roles.


## Full-MM reconstruction scientific checkpoint — October 20

`ACTIVE`: memory-only checkpoint under the
[October 20 contract](../phases/kaonlt-full-mm-reconstruction-scientific-memory-checkpoint-2026-10-20-task-contract.md),
which carries the supplied artifacts, histogram crosswalk and prior ChatGPT
review. Local `test` HEAD and local `origin/test` were independently checked at
`6d13abb4681b938c56def2a1f455244b0880096c`; no live remote query or farm run.
The following RUNTIME VERIFIED observations are attributed to those supplied
artifacts/prior review, not independent local raw-artifact revalidation.
Earlier sections retain their historical statuses; the former A/B farm-pending
next step has been consumed only within the scope below. CURRENT alone owns NEXT.

### Controlled comparison and bounded histogram crosswalk

SOURCE VERIFIED / RUNTIME VERIFIED (scope only): the reviewed, user-pushed
`6d13abb4681b938c56def2a1f455244b0880096c` gate produced the accepted October
19–20 Center/Left lowe canonical t1/t2/t3 detached comparison JSON/PDFs and
receipts. Exact comparison basenames are
`Center_kaon_Q4p4W2p74_lowe_pion-t-template-comparison.json` / `.pdf` and
`Left_kaon_Q4p4W2p74_lowe_pion-t-template-comparison.json` / `.pdf`, as declared
in the [execution contract](../phases/kaonlt-pion-t-template-science-execution-gate-2026-10-19-task-contract.md).
No unprovided receipt or archive hash is invented. Local source confirms the
raw-weighted, normalized-once t-local candidate and fixed-parent-alignment,
diagnostic-only/non-production boundary in `src/cuts/pion_t_template_comparison.py`.
Template-comparison A uses inherited/global shapes; B uses those t-local shapes,
with the same data/configuration and protected-signal references. Center/t2
reported pion-control joint chi2/ndf worsens ~5.72 -> ~6.53; kaon-side **staged
full-range validation** worsens ~18.33 -> ~41.66. The latter is not final
protected-background chi2. Both protected-delta candidates are rejected with
applied `pi_delta=0`; staged-fit totals are not finally applied spectra. B can
improve diagnostics in other bins; no overall scientific preference or adoption.

The earlier `KaonLT_loweps_20261018_093137.zip` contains the 19-page Center
canonical t1/t2/t3 seven-stage JSON/PDF. Its producer still says
`captured_pending_independent_validation`, `observed_source_HEAD=not_provided`.
The supplied crosswalk found exact Center/t2 input/model histogram fingerprints
and physical coordinate shifts matching newer baseline A. This bounded match
does not upgrade the older owner/checker/HEAD provenance or establish a broader
provenance bridge. Its earlier capture limitations remain intact.

RUNTIME VERIFIED (supplied numerical crosswalk): Center/t2 fitted kaon pion-model
integral ~0.25285 versus signed event-applied ~6.5562 are normal-bin full-MM
integrals, not cross sections. About 94.6% of their ~6.303 difference lies at
MM 0.94–1.03 GeV. Signed pion control there is largely prompt; the tiny fitted
control denominator drives `w0` to ~182.82. Post-proton kaon data in that region
are small: ~0.001661 versus full-range ~1.37076. This is a transfer-stability/
accounting concern, not evidence of a major physical low-MM background or proven
signal bias. The negative neutron-anchor LS numerator and tiny refined control
neutron coefficient remain observations. INFERENCE: denominator amplification
may matter; unique particle-level cause and distortion of the signal yield/error
remain NOT VERIFIED. No weight clipping or removal of signed events is authorized.

### Scientific evaluation and nomenclature

The user's correction restores the overall L/T-separation objective: reconstruct
full experimental MM after random/dummy subtraction from proton contamination,
pion backgrounds and K+Lambda/K+Sigma0 SIMC under explicitly common selection,
coordinate, binning and charge/efficiency/acceptance normalization. Inspect actual
spectra, signal-region zooms, Lambda/Sigma0 positions, widths, tails, relative
structures and residuals before choosing a model from background chi2. Then test
the properly normalized background-subtracted K+Lambda data/SIMC yield ratio in
matching `(t,phi,epsilon)` bins through established extraction. Movement toward
unity is a consistent-model goal, never a constraint imposed by area-matching or
arbitrary final-signal SIMC normalization. Diagnostic background amplitude and
physically justified MM-offset nuisance fits may study shapes if documented;
they authorize no production shifts/scales, template retuning to noise, changed
weights or suppressed uncertainty. No new fit is made in this checkpoint.

Template-comparison A/B are not HGCer Method A/B. HGCer Method A is a candidate
reweighting independently validated only in limited scopes, preserving signed
parent normalization rather than independently correcting it; it stays detached/
non-production without explicit promotion. HGCer Method B stays diagnostic-only
and numerically excluded. Kaon-side coefficient `A_n` and pion-control `B_n` are
also distinct: `w0=K_model/P_model` correctly transfers pion control to kaon
background and must not be inverted. Baseline remains authoritative.

`SIMC_normfac_luminosity_and_charge_units_not_source_proven` still blocks the
relevant absolute-SIMC interpretations/pages; full absolute data/MC closure or
unity is not established. Canonical-five Method-A/F.6.4 lineage blockers and
historical ordinary-launcher isolation blockers remain separately `BLOCKED`.
All earlier E.8/F.6 closures retain their exact scopes. Preserve active
`no_empirical_residual`, random/dummy -> proton -> pion ordering, signed events,
SIMC, yields/errors, acceptance, kinematic bins, cuts, templates, priors and paths.
Historical first-bad SHA, all-setting acceptance, a provenance bridge beyond the
bounded fingerprint/coordinate comparison, final K-Lambda yield agreement and
physical origin of low-MM pion control remain NOT VERIFIED.

After this checkpoint's actual-diff review, user publication and pushed-state
synchronization, return directly to the read-only existing-product full-MM
component/source/normalization inventory and reconstruction feasibility audit
specified in CURRENT: Center/lowe t2 first; Center t1/t3 and separate Left/lowe
controls only after matching input identities. Show admissible existing MM plots
and signal-region zooms, determine whether matched authoritative yield/SIMC
objects support a credible comparison, and identify the first missing/blocked
observable. Stop at the first failed evidence gate. This is the next substantive
audit, not another documentation cycle, new fit, diagnostic implementation
contract, promise of new farm artifacts or farm command.

## Detached local-t transfer-variance candidate C — October 20

`ACTIVE`: the user explicitly directed the
[candidate-C contract](../phases/kaonlt-pion-local-t-transfer-variance-candidate-c-2026-10-20-task-contract.md),
returning to pion modeling before the previously recorded existing-product
inventory. CURRENT alone owns the revised NEXT. This is template-comparison
A/B/C, not HGCer Method A/B or fitted A_n/B_n; no new scientific phase.

SOURCE VERIFIED at local `test` / local `origin/test`
`52d8b44335efde614111f727d90fd7987f86d941`: the new detached reader uses the
already-aligned B fitting columns, explicit signed-control Sumw2, original fit
bins, nonnegative bounds, staged amplitudes and original joint prior widths.
It holds `H_kaon_pion_bg_fit_total_refined` fixed, verifies protected-applied
component closure and retains A/B input snapshots unchanged. No production
source, configuration, authoritative fitting or event application changes.

The locally found accepted A/B ZIP SHA-256 is
`8bd71f9ab17fa456a86911b8b43eed9c2a1a9da3101cbbff776c56f66fa901ee`.
Its four artifact hashes/sizes and manifest identity were checked before use;
both JSON source HEADs are `6d13abb4681b938c56def2a1f455244b0880096c`.
Center JSON SHA-256:
`7e34ec134ef28378ffebf1686c47915a266478531d9910bbdd211f5392cb1f24`;
Left JSON SHA-256:
`f5c8cf2483ce77d74dc2e2a497d8b7c5f1b0aa3c96b4d76c3964d813e675b812`.
Original acceptance remains attributed to prior supplied farm evidence/review;
this local archive/array verification is not a new farm run.

Local deterministic evidence: all six reconstructed B joint fits and their
data/prior objective terms match recorded B within declared numerical tolerances.
Maximum coefficient discrepancy is `1.2671e-8`, below `1e-7`. The predeclared
lambda grid is `0, .01, .1, 1, 10, 100, 1000`; its added objective is
`lambda * V_full(beta)/V_full(B)`, with the normalization fixed before the sweep.
Full normal-bin MM includes 0.94–1.03 GeV. All 42 fits converge; no unsupported
positive-kaon/Sumw2 normal bins are found in these results. Strengths are
unselected; convergence is not model acceptance.

SOURCE VERIFIED (deterministic local arrays only): Center/t2 full-MM risk changes from
about `188.968` at lambda zero to `0.521317` at lambda 1000, while signal-window
1.10–1.16 GeV risk stays approximately `3.63e-6`; data chi2/ndf changes from
`6.52955` to `6.52510`. These values concern local-template B/C, not inherited
baseline A or the historical maximum weight near 183. Several other parents
show worse data chi2 at large strengths; for example Center/t1 changes from
`4.01346` to `5.69835`. Retain the complete sweep and shape residuals; no
strength is validated or selected, and no scientific preference follows.

Outputs are detached local JSON and a 24-page local PDF, separate from accepted
farm/analysis-procedure products. Risk, histogram-applied background integrals,
weight tails, component contributions, bound hits, coefficient sensitivities and
conditional covariance are retained by setting/parent/strength. Signal windows
use explicit bin centers, not authoritative extracted-yield integration.
All-cuts Lambda application, RF-restored yields, independent event replay,
held-out validation and detailed shared-event/fit covariance remain NOT VERIFIED
or explicitly unavailable. The fixed-denominator Sumw2 proxy and conditional
curvature are not independent uncertainties or confidence intervals near bounds.
No clipping, suppressed signed event, altered K_B or production adoption.

Local validation: 27 deterministic tests pass, including the supplied-archive
six-parent controls, analytic regularized solution, signed cancellation/Sumw2,
unsupported support, solver nonconvergence, fixed K_B, covariance/bound limits,
provenance/input/path rejection, unchanged A/B and PDF branch smoke checks.
Both new modules compile without importing ROOT. These checks establish only
detached local implementation behavior; independent actual-diff review remains.

The initial memory integrity discrepancy was solely the transferred new contract;
the authorized generator repairs its manifest inventory. The existing CURRENT
soft-size warning is nonblocking and remains scheduled for the next checkpoint.
Independent actual-diff acceptance, user publication and pushed-state review
precede independent six-parent feasibility evidence review; any farm-integrated
C or weight-rule change requires a later separately reviewed decision. Preserve
all earlier closures, normalization/isolation blockers and Method-A/Method-B roles.

## Detached signed-event validation implementation — October 21

`ACTIVE` under the [single event-validation contract](../phases/kaonlt-pion-adaptive-c-event-validation-2026-10-21-task-contract.md).
Starting local `test`, HEAD and local `origin/test` were
`1eac918a58bb700e03deb89411e99a38befd3f95`; no live remote query or farm run.
The supplied contract records the completed numerical/visual C feasibility
review and establishes independent event-level evidence before any selector.
That review attribution is not a new farm acceptance. CURRENT's consumed NEXT
is corrected within this change, with broader science and blockers preserved.

SOURCE VERIFIED against current source/proposed diff: the one default-off
`LT_ANALYSIS_PION_LOCAL_T_EVENT_VALIDATION=1` observer follows both existing
yield-cache assignments. It copies scalar canonical-t pion-control records
(including outside-phi support), verifies their physical child-cache projection,
and copies kaon-yield rows with original post-proton/RF factors. No ROOT event
loop, cache construction change, new fit, physical selection or normalization.
Data/dummy hashes, signed coefficients, source/tree-local entry identities,
coordinate/PID/parent/source fingerprints and both cut flags are retained.
Pion canonical endpoints are inclusive at the final edge; the existing yield
cache excludes that final edge. Null phi children remain explicit outside
acceptance. Duplicate source entries and population mixing fail closed.

The standalone capture owner reuses the source-verified disposable-worktree,
overlaid-ltsep, output-link and preservation helpers of the established A/B
owner. One unchanged `-d` child writes two declared Center/Left packets. Freshness,
source/data identities, six parents, finite values, input-histogram/source
closure, ZIP CRC/inventory/hash checks and preservation precede successful
delivery. Existing collector/profile/helper files are unchanged. Capture is
`captured_pending_independent_validation`, never event-replay acceptance.
The separate ROOT-free replay CLI consumes these packets and the accepted
original C JSON; the C artifact need not be on the farm. It pins the supplied
original C hash, accepted A/B hashes and capture HEAD, joins source/coordinate/
canonical-edge/parent/model identities and never refits C. It retains A, B and
all seven strengths; B/C hold the protected-applied K_B fixed. Source and phi
spectra, signed/absolute/positive/negative support, quadratic Sumw2 and effective
support are separate for actual allcuts and nommcuts. Normal full-MM,
0.7–1.3, 1.10–1.16, 1.105–1.135 and 0.94–1.03 GeV use explicit normal-bin
center windows. Flows are retained with zero application weight as in existing
code. Unclipped diagnostic ratios and the captured production floor/clip
convention are both labeled; no equivalence is inferred from rounding or
unsupported support. These integrals are not authoritative extracted yields.

Deterministic five-fold source-file-hash/label/entry partitions use a predeclared
seed and retain source-weighted training/evaluation balance. They are disjoint
record partitions, not held-out models: saved K_B, staged amplitudes/priors and
protected signal fits consumed full-data events. The first missing independent
fitted-estimator step is training-only staged/protected fits and priors; this
change supplies no replacement fitter or full-data-prior fallback. Local SIMC
shapes alone provide no independent kaon-PID leakage truth with source
normalization. Report `bias_truth_unavailable` and `MSE_selection_unavailable`;
no adaptive lambda choice. Fixed-weight Sumw2, parameter covariance and true
repeated-training variance remain distinct. Synthetic injected-transfer controls
are analytic test evidence only, not physical bias or runtime acceptance.

Local deterministic tests include six-parent synthetic replay, negative
provenance/support/population gates, default-off and ROOT-import isolation,
original-module AST preservation outside the guard, owner failures and the real
unchanged generic collector with synthetic packets. The prior local C JSON and
accepted A/B archive are source/schema/numerical controls, not fresh signed
events. Synthetic figure layouts were inspected as PNGs; actual PDF encoding/
rendering, ROOT/PyROOT, full main.py integration and farm capture remain NOT
VERIFIED. No original C artifact was altered, no farm executed and no strength
selected. Independent actual-diff review remains pending. The initial manifest
inventory discrepancy was solely this transferred contract; existing entries
matched. CURRENT's soft-size warning is nonblocking and scheduled for the next
checkpoint/milestone. All earlier science closures, absolute-normalization and
ordinary-launcher blockers, Method-A/Method-B roles and production physics stay
frozen. CURRENT alone owns the next evidence gate.

Local validation receipt: 156 deterministic tests passed with no skips across
the new capture/replay/owner suites and existing C, A/B comparison/owner, parent
integrity, baseline-stage and pre-particle-subtraction binning suites. Six
changed/new Python modules compile. This is local/source evidence only.

### Enabled capture interface correction — October 21

`ACTIVE`: the user's independent actual-diff review identified two capture
adapter blockers missed by the initial 156-test suite: treating a source's
GetName as a backing filename, and requiring parent_fit_configuration instead
of observing the existing resolved production weight policy. The prior local
checks did not establish enabled-observer viability. This correction remains
inside the unchanged October 21 contract and eleven-path allowlist; CURRENT,
its scientific NEXT, source observation HEAD and all science boundaries remain.

SOURCE VERIFIED against current source/proposed diff: the adapter resolves a
tree's GetCurrentFile, or an existing TFile's GetFile, checks the backing file's
open/non-zombie state, and hashes its actual filename. A tree name is never a
path fallback. Missing backing-file provenance fails before packet creation;
no file opening, tree iteration, source mutation or ROOT directory change.
The adapter captures only floor/clip values returned by the existing
background_config production resolvers. It reads no parent_fit_configuration
field, and checks any retained final/proposed application policy diagnostics
for disagreement. Absent optional diagnostic payloads do not replace resolver
values with guesses. Production configuration and application remain untouched.

Accepted A/B data/dummy SHA-256 and byte-count projections were independently
recovered locally from both original JSONs in the archive already identified
above. Archive SHA-256 and both JSON SHA-256 identities were checked first;
input source HEAD remains 6d13abb4681b938c56def2a1f455244b0880096c.
The hash/size-only pins carry no machine paths. They are checked immediately
after backing-file hashing, before parent serialization or packet output, and
again by the owner before collection. No C file or guessed farm path is needed
during capture; offline replay retains its independent accepted-C/input checks.

Local validation: 163 deterministic tests passed, no skips, and six Python
modules compile. Positive enabled-observer tests use realistic TTree/TFile and
TH1 interfaces, parent records without parent_fit_configuration, actual retained
histogram serializers, and actual production configuration resolvers. They
verify identity/policy output and unchanged caches, parents, input dictionaries,
source handles, backing-file state and bytes. Negative tests cover absent,
closed and zombie backing files, path-looking tree names without provenance,
accepted-input mismatch before parent serialization, retained floor/clip
configuration disagreement, and self-consistent packets with nonaccepted source
hashes rejected before owner collection. Synthetic tests explicitly substitute
synthetic input pins; a separate archive test checks the real source-owned pins
against both original accepted A/B JSONs. No synthetic test is farm evidence.

Actual ROOT/PyROOT, full main.py integration, farm capture and scientific
acceptance remain NOT VERIFIED. Independent actual-diff review is still the
current workflow gate; user publication and pushed-state review precede the
same owned six-parent capture/replay scientific NEXT. No commit, push, farm
execution, new contract, selector, model choice or production change follows.

### Event-capture prepass boundary failure and repair — October 21

`BLOCKED` at the first farm child for source
`04f1ec403a9555634e61c229aa347cb76d0c78d3`. RUNTIME VERIFIED only as
user-supplied failure evidence transcribed in the
[repair contract](../phases/kaonlt-pion-event-capture-prepass-boundary-fix-2026-10-21-task-contract.md):
owner `analysis_returncode=1`, `stage=completion_markers`,
`failure_reason=analysis_failed`; packet and collection checks not reached.
Ordinary-worktree and installed-ltsep preservation passed per that receipt.
The child log identifies `bg_optimizer_pre_parent`, scalar `single_scale`,
no required parent/free component fit, then
`accepted_AB_input_identity_mismatch:Center:data` from the observer. Failed
outputs/receipts remain immutable and inadmissible as physical C event evidence.

SOURCE VERIFIED: optimizer candidates legitimately calculate preliminary yields;
the opt-in-only hook incorrectly admitted that pre-parent stage. The narrow
repair requires kaon `simc_shape_components`/`t_bin` and excludes either prepass
flag or context before diagnostic import, hashing or output. The standalone
observer enforces the same boundary. Existing frozen-parent validation and cache
population remain unchanged; eligible missing/incoherent parents or caches fail
rather than disappearing. All accepted A/B input SHA-256/byte pins and their
fail-closed comparison remain unchanged. The mismatch is not established to be
solely a prepass effect: a canonical mismatch must stop for independent
scientific/provenance adjudication, with no override or replay.

SOURCE VERIFIED file-identity derivation: `rand_sub.py` obtains `OUTPATH` from
its existing ltsep Root configuration, and receives `InDATAFilename` and
`InDUMMYFilename` through the main input dictionary (CLI arguments 6 and 7).
For each Center/Left setting, the exact input paths are constructed as
`OUTPATH / (phi_setting + "_" + ParticleType + "_" + InDATAFilename + ".root")`
and the corresponding `InDUMMYFilename` expression. `open_root_file` returns
`TFile.Open`; those handles are retained in `hist["InFile_DATA"]` and
`hist["InFile_DUMMY"]`. Existing `GetFile`/`GetCurrentFile` resolution is
preserved, with no tree-name fallback. An exact read-only check can use
`file_identity` and `verify_accepted_input_identity` on those two derived paths
per setting, given the actual active ltsep path dictionary, exact invocation
filenames and readable unchanged farm files. Those external inputs are not
available locally; current farm hashes/bytes are NOT VERIFIED. No farm path,
pin replacement or physical-equivalence assumption is introduced.

Local repair is `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; validation
receipts below are source/synthetic evidence only. Production physics, optimizer,
parent fits, yields and published owner/profile/collector remain frozen.
CURRENT alone retains the same six-parent capture and detached A/B/C replay
NEXT, conditional on independent actual-diff review, user publication,
pushed-state synchronization and readiness checks. No farm command or automatic
retry is authorized; no lambda selection or production promotion.

Local validation receipt: 170 tests ran across capture/replay, capture owner,
optimizer prepass, frozen-parent integrity, pre-particle-subtraction binning,
baseline-stage audit and existing C/A/B comparison/owner suites: 161 passed,
nine skipped (two original-artifact controls not supplied and seven ROOT-dependent
checks). After tightening the test's file-I/O tripwire, the affected 28-test suite
was rerun: 26 passed, the same two original-artifact controls skipped. Three
changed Python modules compile. Actual dispatch AST tests execute both normal
cache assignments and the guard: repeated prepass flag/context/scalar candidates
perform no diagnostic import or file I/O; canonical Center/Left dispatch once
per setting and preserve the source handles, caches and parent histograms.
Canonical wrong data/dummy identity still raises its precise invariant before
parent serialization. Defective eligible parent/cache inputs fail with no
packet. Disabled dispatch is equivalent. These checks establish no ROOT/farm
integration, fresh event acceptance or scientific preference. CURRENT's soft-size
warning remains nonblocking, scheduled for the next checkpoint/milestone audit.
