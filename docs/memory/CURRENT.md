---
memory_schema: 3
---
# Current KaonLT development state

## Active Objective

E.8 is `ACTIVE`: correct the [Method-A/reweighted-baseline roles](decisions/e8-full-analysis-procedure-roadmap.md#live-four-way-producer-audit).
Existing private `YA/MM_A` is the Method-A-reweighted baseline, not a
Method-A-only final kaon result. Detached weight/response review has priority.
Microscopic F.4 discrepancy investigation and F.1 archaeology are `DEFERRED`;
accepted scopes, exact failed gates and frozen science remain intact.

## Current Work Item

E.8.2 Q4p4W2p74 / Left / lowe remains `CLOSED / RUNTIME VALIDATED`.
[Accepted closure](evidence/e8-2-left-lowe-runtime-scientific-closure-2026-10-06.md)
owns its narrow visual/scientific conclusions. SOURCE VERIFIED: the final
canonical plot consumes after-baseline-pion `after_pion_final` / baseline
`MM_0(t,phi)` and stored authoritative `Y0`. Diagnostic Lambda-window stage
integrals (including after pion) are not final `Y0` or a new subtraction stage.
No legacy residual-fit stage is active.

E.8.3 Left/lowe and its [authority-page repair](phases/e8-3-authority-page-provenance-legibility-repair.md)
are `CLOSED / RUNTIME VALIDATED`. [Closure](evidence/e8-3-left-lowe-runtime-closure-2026-10-07.md)
owns supplied/prior ChatGPT acceptance. Historical F.6.1 lineage remains
separate from baseline E.8.2 and private current F.6.3/E.8.4.
SOURCE VERIFIED: E.8.4 consumes explicit `MM_0/MM_A` and stored `Y0/YA`
from that private branch. F.3 independently derives relative response `r`, not
absolute mis-ID probability. F.4 uses signed baseline `b_j=s_j*w0_j`,
`B=sum b_j`, `U=sum b_j*r_j`, and valid-parent `C_j=r_j/(U/B)` to preserve
`sum b_j*C_j=B`; normalized `C` depends on baseline support, without child
renormalization. F.6.3 fills `B_pi_A` with `s_j*w0_j*C_j`, subtracts it from
`pion_input` to form `MM_A`, and extracts the stored yield with the existing error helper.
Thus `YA/MM_A` already owns the final reweighted-baseline result. Standalone
Method-A-only final yield/MM is NOT VERIFIED; four-category presentation intent
must distinguish response/weight diagnostics from final yield/MM curves.
INFERENCE: a cross-branch or diagnostic-integral
comparison could cause an apparent discrepancy; none is proved resolved.

The [fresh-F1 canonical-five preflight](evidence/e8-4-canonical-five-lineage-preflight-runtime-closure-2026-10-05.md)
is `CLOSED / RUNTIME VALIDATED` pre-analysis only (analysis_started=false).
Supplied receipt/prior review owns its lineage/isolation/preservation scope.

The historical [second full run](evidence/e8-4-canonical-five-second-full-run-provenance-blocker-2026-10-05.md)
is `BLOCKED` by post-analysis page/artifact provenance. Its reported
F.2/F.3/F.4 equality does not accept the newer F.4 mismatch.

The [F.2 handoff/equivalence](phases/e8-4-canonical-five-f2-handoff-equivalence-repair.md),
[F.2 staging](phases/e8-4-canonical-five-f2-staging-predecessor-repair.md) and
[regenerated-F.1 preflight](phases/e8-4-canonical-five-f1-lineage-preflight-repair.md)
repairs remain `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`.
Strict validators and reviewed/current provenance stay frozen; unknown staging
targets fail closed. Historical limits remain in the
[October 8 attempt](evidence/e8-4-canonical-five-f1-lineage-preflight-blocker-2026-10-08.md),
[accepted audit](phases/e8-4-canonical-five-provenance-identity-repair.md) and
[first-run diagnosis](evidence/e8-4-canonical-five-fresh-f1-lineage-materialization-2026-10-05.md).
Local checks add no farm acceptance.

The newer October 8 owner at `ebf3a0be2f50c788c22defe3601fd5b7fee48644`
is `BLOCKED` before analysis at
`f6_3_scientific_equivalence_mismatch:f4:$.parents[0].absolute_baseline_parent_sum`.
Supplied comparison: F.2/F.3 exact PASS; F.4 exact FAIL in three Left-lowe
parents. The [diagnostic contract](phases/e8-4-f1-baseline-reproducibility-diagnostic-task-contract.md)
owns the failed comparison/search limits. The
[detached diagnostic](phases/e8-4-f1-baseline-reproducibility-diagnostic.md)
retains `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`.
The pushed CLI at `6cdd2fca` has a
[supplied diagnostic review](evidence/e8-4-f1-forensic-reproducibility-diagnostic-review-2026-10-08.md):
three counts agree, all six signed/absolute sums differ exactly, reviewed F.1
is absent. Prior ChatGPT review owns these observations; forensic parsing
grants no canonical acceptance. Full F.1 availability/event-level cause are
NOT VERIFIED. Further F.1 archaeology is `DEFERRED` unless concrete evidence
makes it material or unavoidable for a required gate.

The [accepted detached Left/lowe diagnostic](evidence/e8-4-left-lowe-method-a-current-lineage-diagnostic-runtime-closure-2026-10-04.md)
remains `CLOSED / RUNTIME VALIDATED` for t0/t1/t2, primary t1. With accepted
F.6.2/comparator t1 continuity, no redesign is warranted (INFERENCE).
Absolute mis-ID, proxy validity, hardware cause and promotion are NOT VERIFIED.

Independent E.8.2 baseline-five `Y0/MM_0` is `CLOSED / RUNTIME VALIDATED`:
[evidence](evidence/e8-2-canonical-five-baseline-yield-runtime-closure-2026-10-09.md)
owns farm identity, prior ChatGPT acceptance, 135/135 data and SIMC table/MM
gates and representative PDF inspection limits. Private `YA/MM_A`/F.4 stays blocked.

[Baseline-five QA](evidence/e8-baseline-five-extracted-yield-scientific-qa-2026-10-11.md)
records detached inspection/caveats only, with no new phase or physics acceptance.

## Verified State

Gate-3 ownership: NPE=0 kaon PID, NPE>0 pion tree, NPE>2 physical control.
F.3 `hgcer3` is relative response; baseline `w0` owns transfer. F.4 preserves
signed canonical-t parents, never independent children/MM regions. The HGCer
hole is excluded; zero-response transfer is diagnostic-only; slow-proton
architecture is analogue-only. The Gate-3 reference owns source evidence.

Prior comparator t1 continuity does not extend to substantive t0/t2 changes;
F.6.2 validates acceptance/MM correlation, not absolute HGC.
F.4.Refresh.2 remains `CLOSED / RUNTIME VALIDATED`;
[F.6.3/E.8.4 closure](evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md)
is only Left/lowe execution, live-cache parity, real child changes and signed
parent preservation. Fix.5.7/Fix.5.8 retain narrow Left/lowe provenance and
[97-page legibility closure](evidence/e8-4-fix5-8-left-lowe-runtime-closure-2026-10-03.md).
Earlier closures retain [roadmap](roadmap/STATUS.md) scopes. Workflow hardening
is `SOURCE REVIEWED`; literature adds no runtime evidence or authority.

## Source / Evidence Identity

Linked records own farm sources/artifact identities; none asserts live HEAD.
Current reviewed F.4 SHA-256:
`79e7ceda7221cbeeead4ed5bc306b0e0e670741a27beaa980e22349c555e96d7`.

## Blockers

Private Method-A canonical-five runtime remains `BLOCKED` by the exact F.4 mismatch;
owner retry is unauthorized pending separately reviewed resolving evidence.
Final E.8/F.6.4 remain `BLOCKED`; production correctness is unproven.
Optional low-epsilon RF was not performed.

Absolute-SIMC interpretation remains separately `BLOCKED`:
`SIMC_normfac_luminosity_and_charge_units_not_source_proven`.
`iter_weight * normfac / Ncontribute` has unproven luminosity/effective-charge
units. Only two absolute-SIMC page families are unavailable; valid data/identity/
yield/parent-closure payloads remain available. No scale/amplitude conclusion.

## Next Action

NEXT — after ChatGPT actual-diff acceptance, user commit/push and pushed-state
synchronization, scientifically review the existing detached Method-A weights:
distinguish relative response `r` from parent-normalized `C`; trace input
populations, range, shape, uncertainties/support and connection to `w0`;
identify standalone weight diagnostics comparable with baseline without
inventing a final-yield producer. F.4 microscopic investigation and F.1
archaeology stay `DEFERRED`; exact owner and absolute-SIMC blockers remain
binding. No automatic correction implementation, owner retry or farm work.

## Success Criteria

Manifest, health, bootstrap, preservation and diff must pass.
Source audit is not four-way physics completion; no staging, commit, push or farm.

## Do Not Reopen Without New Evidence

Keep `no_empirical_residual`, zero legacy scales, detached/non-production
Method A and numerically excluded diagnostic Method B. Freeze subtraction,
SIMC, weights, yields, uncertainties, cuts, templates, priors, binning,
efficiencies, acceptance, L/T and cross sections. No child normalization,
authority replacement or promotion. E.8 consumes; F.6.3 owns the private
branch. Failed artifacts cannot close a gate; owner success is insufficient.

## Relevant References

- [Gate-3 source/science investigation](investigations/e8-4-left-lowe-method-a-detector-response-source-science-audit-2026-10-03.md)
- [Current-baseline comparator](evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md)
- [Prior Fix.5.7 evidence](evidence/e8-4-fix5-7-left-lowe-runtime-and-fix5-visual-gate-2026-10-03.md)
