---
memory_schema: 3
---
# Current KaonLT development state

## Active Objective

`ACTIVE`: baseline-pion interruption, **Q4p4W2p74 / Center / highe / t3**.
Preserve iteration-zero subtraction. Local-t is the user-selected scientific
baseline for continued development and final analysis, not production-adopted
physics. Global and Local-(t,phi) remain historical comparisons. INFERENCE from corrected spectra: the narrow
Lambda peak is largely preserved; its background is mixture-sensitive.
Full-spectrum chi2 is not the primary adequacy criterion beneath Lambda.
High epsilon has no RF slow-proton subtraction. Missing configured Sigma0
input is a separate improvement, not this blocker. Method A remains detached,
B diagnostic/numerically excluded; normal [roadmap](roadmap/STATUS.md) forward
work stays paused. No accepted scope is reopened.

## Current Work Item

The independent joint anchor-mean gate remains blocked on the missing
independent component-mean evidence specified under Next Action.

`BLOCKED`: independent joint anchor-mean gate, **Outcome B**. SOURCE VERIFIED:
scalar empty-error handling is not unit-covariant; staged/gauge fits reuse the
kaon data. November 3 source-aware diagnostics improve some conditional laws,
not universal coverage or physical truth; they are not adopted. Acceptance/
fallback changes actual subtraction. No independent normalized constraint on
non-SIDIS contamination of the SIDIS anchor and pion contamination of the
Lambda gauge was found. SIMC shapes and same-data amplitudes cannot supply it.
A newly inspected retained Left/low zero-photoelectron diagnostic infers
response from positive-NPE controls; it supplies neither required Center/high/t3
mean nor their joint covariance. No speculative correction contract or implementation.

See the [November 3 checkpoint](investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md#scientific-memory-checkpoint-and-independent-anchor-gate--november-3)
for hashes, source trace and results. Historical details below are linked context, not the active NEXT.

## Verified State

The [relocated historical details](investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md#historical-status-details-retained-from-current)
retain the prior CURRENT statements, exact source/artifact limits, failed
attempts and contract links. They are historical context; no scope is reopened,
closed or promoted by this consolidation. Earlier acceptance is MEMORY ONLY
in this memory pass, attributed to its governing evidence records.

| Scope | Preserved state and boundary |
| --- | --- |
| E.8.2 and E.8.3 Left/lowe | CLOSED / RUNTIME VALIDATED within [E.8.2](evidence/e8-2-left-lowe-runtime-scientific-closure-2026-10-06.md) and [E.8.3](evidence/e8-3-left-lowe-runtime-closure-2026-10-07.md); historical F.6.1 is separate from current F.6.3/E.8.4. |
| E.8.2 baseline-five Y0/MM_0 | CLOSED / RUNTIME VALIDATED: [135/135 data and SIMC gates](evidence/e8-2-canonical-five-baseline-yield-runtime-closure-2026-10-09.md), original PDF limits; no added component-physics acceptance. |
| F.4.Refresh.2 and Left/lowe F.6.3/E.8.4 | CLOSED / RUNTIME VALIDATED within [accepted scopes](evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md); [Fix.5.7/5.8](evidence/e8-4-fix5-8-left-lowe-runtime-closure-2026-10-03.md) retains 97-page legibility closure. |
| October 12-13 snapshot/presentation | CLOSED / RUNTIME VALIDATED within [recorded limits](evidence/e8-current-snapshot-contribution-runtime-closure-2026-10-12.md); internal coherence only, canonical-owner equivalence BLOCKED; no Center MM/PDF closure. |
| Canonical-five lineage-only preflight | CLOSED / RUNTIME VALIDATED with [analysis_started=false](evidence/e8-4-canonical-five-lineage-preflight-runtime-closure-2026-10-05.md), not scientific-equivalence acceptance. |
| Paired lowe extension; full-chain/component diagnostics; F.1/F.2 staging; resolver/receipt/history repairs | DEVELOPMENT COMPLETE, FARM VALIDATION PENDING within their original contracts; no broader acceptance. |

## Source / Evidence Identity

Local `test`, `origin/test` and public `test` matched
`5d5c372007b926bd5245e0f0d8283c98a06da851` on November 3 before this roadmap
update. Only memory differs from farm scientific source
`9172aba853fd5b5b744c19add992fa4f27466180`.
RUNTIME VERIFIED: supplied October 30 `20261030T022205` high-epsilon acquisition
at the farm scientific source completed analysis, collection and nine-artifact verification.
Stored `evidence_admissible=false` and `scientifically_accepted=false` remain.
That acquisition consumes the old ledger/fresh-high NEXT, not scientific acceptance.

## Blockers

Separately BLOCKED: canonical-five Method A/final E.8/F.6.4 at
`f6_3_scientific_equivalence_mismatch:f4:$.parents[0].absolute_baseline_parent_sum`;
microscopic F.1/F.4 archaeology remains DEFERRED. Absolute SIMC remains BLOCKED
at `SIMC_normfac_luminosity_and_charge_units_not_source_proven`; this affects
absolute interpretation, not valid data-only/identity/yield/parent closure.
The [historical SIMC Work audit](investigations/simc-absolute-normalization-ledger-2026-11-03.md)
is retained with provenance, not adopted as a new prerequisite or demonstrated
normalization correction. The [November 3 accountability correction](investigations/kaonlt-assistant-roadmap-displacement-2026-11-03.md)
requires preserving mature functionality and the approved priorities; no repeated
audit or generator archaeology follows without concrete new defect evidence and
user authorization. Existing fitting/subtraction machinery is available for
the recorded bounded blind predictive closure once qualifying new evidence is
independently reviewed. Outcome B and the scientific NEXT are unchanged.
The earlier normal-launcher isolation blockers and matched-cohort ninth failure
(`external_input_verification/source_changed_during_hash`) remain in the
historical details and [incident record](investigations/kaonlt-matched-cohort-farm-readiness-failure-2026-10-22.md).
No fresh matched-event physical C or held-out fitted estimator is accepted;
no strength is selected. None of these separate scopes defines the present NEXT.

## Next Action

NEXT — obtain one independently derived, source-identified joint constraint
record for **Center/highe/t3**: non-SIDIS mean in the existing kaon SIDIS
anchor `[1.0601751288600272, 1.1101751288600272]` GeV and pion mean in the
Lambda gauge `[1.105, 1.125]` GeV, with joint uncertainty, covariance to the
evaluated signed counts and disjoint/calibrated provenance. It must not be
derived by relabeling the evaluated staged/gauge fits as truth. Require
independent scientific review of the new constraint before one blind predictive closure of those windows and the actual
iteration-zero Lambda-window subtraction with frozen inputs and conventions.
This remains the scientific NEXT after memory review/publication and
pushed-state verification; no farm acquisition or correction contract yet.

## Success Criteria

No production statistical correction or farm readiness is accepted.
Memory publication and pushed-state synchronization do not change scientific NEXT.

## Do Not Reopen Without New Evidence

Preserve `no_empirical_residual`, zero legacy scales and all current cuts,
templates, priors, binning, alignment, normalizations, weights, subtraction,
efficiencies, acceptance, yields/errors, L/T and cross sections. No child
renormalization. `YA/MM_A` is the reweighted-baseline branch; standalone
Method-A-only final yield/MM is NOT VERIFIED. Frozen reviewed F.4 SHA-256:
`79e7ceda7221cbeeead4ed5bc306b0e0e670741a27beaa980e22349c555e96d7`.

## Relevant References

The [investigation](investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md)
owns numerical evidence/history; [MEMORY](MEMORY.md) owns reusable knowledge;
the [roadmap](roadmap/STATUS.md) retains the interrupted approved progression.
