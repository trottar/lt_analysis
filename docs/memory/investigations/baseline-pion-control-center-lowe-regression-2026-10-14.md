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
