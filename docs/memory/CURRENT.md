---
memory_schema: 3
---
# Current KaonLT development state

## Active Objective

E.8 is `ACTIVE`: audit actual Method-A-reweighted baseline contributions.
The [role correction](decisions/e8-full-analysis-procedure-roadmap.md#live-four-way-producer-audit)
identifies private `YA/MM_A` as the reweighted baseline, not Method-A-only.
Microscopic F.4/F.1 investigations are `DEFERRED`;
Accepted scopes, failed gates and frozen science remain intact.

[Current-snapshot debug](phases/e8-current-snapshot-method-a-contribution-debug-task-contract.md):
`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`. Supplied preflight: FAIL at
`snapshot_failed:f3_fingerprint_input_content_mismatch`, outputs_written=false.
[Detached materialization](phases/e8-current-f1-coherent-detached-materialization-task-contract.md)
is `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`, conditional on independent
actual-diff review and user publication.
Current farm coherence NOT VERIFIED; October 4 acceptance and historical
implementation remain unchanged. No canonical-owner equivalence granted.

## Current Work Item

E.8.2 Left/lowe is `CLOSED / RUNTIME VALIDATED` under its
[closure](evidence/e8-2-left-lowe-runtime-scientific-closure-2026-10-06.md).
SOURCE VERIFIED: final MM consumes `after_pion_final` / `MM_0(t,phi)` and
stored `Y0`. Lambda-window stage integrals are diagnostics, not final yields
or another subtraction; no legacy residual fit is active.

E.8.3 Left/lowe and its [authority-page repair](phases/e8-3-authority-page-provenance-legibility-repair.md)
are `CLOSED / RUNTIME VALIDATED` per [prior acceptance](evidence/e8-3-left-lowe-runtime-closure-2026-10-07.md).
F.6.1 remains separate from baseline E.8.2/current F.6.3/E.8.4.
SOURCE VERIFIED: F.3 fits response `r`, not absolute mis-ID.
F.4 uses `b=s*w0`, `B=sum b`, `U=sum b*r`, `C=r/(U/B)` and preserves
`sum b*C=B`, without child normalization. F.6.3 fills `B_pi_A` with `s*w0*C`,
forms `MM_A=pion_input-B_pi_A`, and stores yield YA via the existing error helper.
E.8.4 consumes `MM_0/MM_A` and `Y0/YA`; standalone Method-A-only final
kaon yield/MM is NOT VERIFIED. Background contributions are not final yields.

[October 11 inspection](evidence/e8-method-a-weight-aggregate-json-inspection-2026-10-11.md)
records prior ChatGPT aggregate recalculation from accepted Left/lowe JSON:
finite tails, parent closure, phi/MM/coordinate redistribution. Marginals do not
prove `w0*C` tail dominance, effective statistics or joint locality. No new
farm run/runtime acceptance.

The [fresh-F1 canonical-five preflight](evidence/e8-4-canonical-five-lineage-preflight-runtime-closure-2026-10-05.md)
is `CLOSED / RUNTIME VALIDATED` for lineage/isolation/preservation only,
pre-analysis (analysis_started=false), per supplied receipt/prior review.

The historical [second full run](evidence/e8-4-canonical-five-second-full-run-provenance-blocker-2026-10-05.md)
is `BLOCKED` by post-analysis page/artifact provenance. Its reported
F.2/F.3/F.4 equality does not accept the newer F.4 mismatch.

The [F.2 handoff/equivalence](phases/e8-4-canonical-five-f2-handoff-equivalence-repair.md),
[F.2 staging](phases/e8-4-canonical-five-f2-staging-predecessor-repair.md) and
[regenerated-F.1 preflight](phases/e8-4-canonical-five-f1-lineage-preflight-repair.md)
repairs remain `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`.
Strict provenance/staging gates stay frozen; unknown staging fails closed. Historical limits remain in the
[October 8 attempt](evidence/e8-4-canonical-five-f1-lineage-preflight-blocker-2026-10-08.md),
[accepted audit](phases/e8-4-canonical-five-provenance-identity-repair.md) and
[first-run diagnosis](evidence/e8-4-canonical-five-fresh-f1-lineage-materialization-2026-10-05.md).
Local checks add no farm acceptance.

October 8 owner `ebf3a0be2f50c788c22defe3601fd5b7fee48644`
is `BLOCKED` before analysis at
`f6_3_scientific_equivalence_mismatch:f4:$.parents[0].absolute_baseline_parent_sum`.
F.2/F.3 exact PASS; F.4 exact FAIL in three Left-lowe parents.
[Contract](phases/e8-4-f1-baseline-reproducibility-diagnostic-task-contract.md)
owns comparison/search limits; [diagnostic](phases/e8-4-f1-baseline-reproducibility-diagnostic.md)
is `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`.
[Prior review](evidence/e8-4-f1-forensic-reproducibility-diagnostic-review-2026-10-08.md)
of pushed CLI `6cdd2fca`: three counts agree, six signed/absolute sums differ exactly,
reviewed F.1 absent. No canonical acceptance; full F.1 availability/event cause
NOT VERIFIED. F.1 archaeology `DEFERRED` unless material to a required gate.

[Detached Left/lowe diagnostic](evidence/e8-4-left-lowe-method-a-current-lineage-diagnostic-runtime-closure-2026-10-04.md):
`CLOSED / RUNTIME VALIDATED` for t0/t1/t2, primary t1. F.6.2/comparator t1
continuity warrants no redesign (INFERENCE). Absolute mis-ID, proxy validity,
hardware cause and promotion are NOT VERIFIED.

E.8.2 baseline-five `Y0/MM_0` is `CLOSED / RUNTIME VALIDATED`:
[evidence](evidence/e8-2-canonical-five-baseline-yield-runtime-closure-2026-10-09.md)
owns identity, prior acceptance, 135/135 data/SIMC table/MM gates and PDF
inspection limits. Private `YA/MM_A`/F.4 stays blocked.

[Baseline-five QA](evidence/e8-baseline-five-extracted-yield-scientific-qa-2026-10-11.md)
records inspection/caveats, no phase or physics acceptance.

## Verified State

Gate-3: NPE=0 kaon PID, NPE>0 pion tree, NPE>2 physical control.
F.3 `hgcer3` is relative response; `w0` owns baseline transfer. Signed parent
preservation does not preserve MM/phi cells. HGCer hole excluded; zero-response
transfer diagnostic-only; slow-proton architecture analogue-only.

Prior comparator t1 continuity does not extend to substantive t0/t2 changes;
F.6.2 validates acceptance/MM correlation, not absolute HGC.
F.4.Refresh.2 remains `CLOSED / RUNTIME VALIDATED`;
[F.6.3/E.8.4 closure](evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md)
is only Left/lowe execution, live-cache parity, real child changes and signed
parent preservation. Fix.5.7/Fix.5.8 retain narrow Left/lowe provenance and
[97-page legibility closure](evidence/e8-4-fix5-8-left-lowe-runtime-closure-2026-10-03.md).
Closures retain [roadmap](roadmap/STATUS.md) scopes; workflow hardening
is `SOURCE REVIEWED`; literature adds no runtime authority.

## Source / Evidence Identity

Linked records own farm/artifact identities, not live HEAD.
Current reviewed F.4 SHA-256:
`79e7ceda7221cbeeead4ed5bc306b0e0e670741a27beaa980e22349c555e96d7`.

## Blockers

Private canonical-five Method A and final E.8/F.6.4 remain `BLOCKED`.
Exact F.4 mismatch forbids owner retry pending reviewed resolving evidence.
Production correctness is unproven; optional low-epsilon RF was not performed.

Absolute-SIMC interpretation remains separately `BLOCKED`:
`SIMC_normfac_luminosity_and_charge_units_not_source_proven`.
`iter_weight * normfac / Ncontribute` luminosity/effective-charge units are unproven.
Only two absolute-SIMC page families are unavailable; data/identity/yield/parent
closure remains available. No scale/amplitude conclusion.

## Next Action

NEXT — after ChatGPT actual-diff acceptance, user commit/push and pushed-state
synchronization, perform one detached coherent-current-F.1 materialization farm
gate. Review its manifest/ZIP before the separate current-snapshot seven-input
preflight. No owner retry or automatic farm work; owner/absolute-SIMC blockers
and `DEFERRED` work stay binding.

## Success Criteria

Manifest, health, bootstrap, preservation and diff must pass.
No physics completion, staging, commit, push or farm.

## Do Not Reopen Without New Evidence

Keep `no_empirical_residual`, zero legacy scales, detached Method A and
numerically excluded diagnostic Method B. Freeze subtraction, SIMC, weights,
yields, uncertainties, cuts, templates, priors, binning, efficiencies, acceptance,
L/T and cross sections. No child normalization, authority replacement or promotion.
E.8 consumes the F.6.3-owned private branch; failed artifacts/owner success
alone cannot close a gate.

## Relevant References

- [Gate-3 source/science investigation](investigations/e8-4-left-lowe-method-a-detector-response-source-science-audit-2026-10-03.md)
- [Current-baseline comparator](evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md)
- [Prior Fix.5.7 evidence](evidence/e8-4-fix5-7-left-lowe-runtime-and-fix5-visual-gate-2026-10-03.md)
