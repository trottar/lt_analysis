---
memory_schema: 3
---
# Current KaonLT development state

## Active Objective

`ACTIVE`: physically credible clean K+Lambda: close cleaned-data/K+Lambda-SIMC
shape agreement in the actual Lambda MM window (peak/relevant tail), matched
selection/binning/resolution/acceptance and documented normalization. Integrated
yields alone are insufficient; close shape agreement is necessary, not proof
of purity. Then data/SIMC ratios flat versus phi for every t/both epsilons;
data itself has physical angular dependence. Initial unity is not required.
Later SIMC iterations approach unity with iteration-zero subtraction fixed.

## Current Work Item

Q4p4W2p74: the November 5 meeting assessment, as reported by the user, finds
generally consistent Lambda peak/relevant-tail shapes in latest Local-t data
and unchanged iteration-zero KLambda SIMC, at parent-t and populated child
level. Local-t remains selected, not production-adopted; this qualitative
assessment proves neither purity, covariance validity nor final yields.
Principal reported exceptions: Center/low and Center/high t3. Center's restricted
high-minus-t coverage makes an acceptance-edge interpretation plausible, not
established. Test it before attributing the behavior to pion subtraction.

Check 1 retained-input assessment is `SOURCE REVIEWED` for the supported
canonical occupancy/MM arithmetic only; evidence remains insufficient for both
Center epsilons. Lambda t3 record shares are 9.8%/13.8% (low/high), versus SIMC
yield shares 7.4%/12.4%; broad-MM records are not similarly depleted. Counts do
not establish occupied t length. Continuous t and joint t-MM/applied-background
information are absent from the named retained inputs, so full accepted-t and
common-interior tests were not produced. See the
[Check 1 record](investigations/full-t-acceptance-edge-check1-2026-11-05.md).

The approved immediate sequence is full-t/acceptance-edge diagnostic, review of
the colleague's exclusion procedure when supplied, remaining beneath-Lambda
pion discrepancies, physically justified conditional Method-A assessment, and
existing final-yield gates; see MEMORY's November 5 meeting section.
Preserve separate unit-normalized shape and existing yield-normalized magnitude
diagnostics. No arbitrary scaling of yield spectra or repeated normalization.
Iteration-zero SIMC uses the old parameterization; initial D/M need not be unity.
Method A stays detached, optional pending a justified target; Method B remains
diagnostic/excluded. Outcome B and source-aware statistics remain separate.

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

Check 1 startup: local/origin/live HTTPS test matched
`6663874fae232729093a7dd51daffd6c710621a8`; published meeting-memory changes
were checked against the reviewed files, allowing the user-run LF conversion.
This Check 1 result/memory update is local; index/refs and pending offline
testing work remain preserved. No new farm source or runtime gate is claimed.
Production physics remains unchanged from scientific source
`9172aba853fd5b5b744c19add992fa4f27466180`.
RUNTIME VERIFIED: supplied October 30 `20261030T022205` high-epsilon acquisition
at the farm scientific source completed analysis, collection and nine-artifact verification.
Stored `evidence_admissible=false` and `scientifically_accepted=false` remain.
That acquisition consumes the old ledger/fresh-high NEXT, not scientific acceptance.

## Blockers

The corrected November 5 detached comparison now uses October 28 lowe /
October 30 highe Local-t experimental allcuts spectra and the existing unchanged
iteration-zero KLambda SIMC yield arrays. The explicit corrective handoff
authorizes retained October 9 SIMC, never its experimental arrays. A change in
experimental subtraction does not require regenerating unchanged simulation.
The previous October 9 experimental comparison remains rejected for this task.
Both existing normalizations stay unchanged: no additional factor, area matching,
fitted amplitude or unit-area template substitution. All 135 slots and 15 parents
are compared, with original signed data and empty states. Original allcuts
boundary contributions permit full-window integrals; earlier interior arrays
are checked unchanged. Nonaligned MM grids give binwise content/ratio/residual
bounds without fractional fills. The 22-page atlas and report are in
KaonLT_Latest_LocalT_Data_Existing_SIMC_20261105.zip. This is detached observer
evidence, not final producer or subtraction-physics acceptance. Conditional
counting errors omit unretained fit/training/shared-systematic covariance.

The selected-Local-t producer-evidence requirement remains unresolved for
final combined angular ratios/yield acceptance: same-cohort producer Y0/MM_0,
errors, validity/skip reasons, source/fit/cut/normalization lineage and available
shared covariance, matched to iteration-zero SIMC. Detached spectra and meeting
agreement do not complete it. This later gate does not block the bounded offline
full-t diagnostic using available data; no producer rebuild or farm run follows.

Production Step 1 remains `BLOCKED`, Outcome B: independently derive the
Center/highe/t3 non-SIDIS mean in the kaon SIDIS anchor
`[1.0601751288600272, 1.1101751288600272]` GeV and pion mean in the Lambda
gauge `[1.105, 1.125]` GeV, with joint uncertainty, covariance to evaluated
signed counts and disjoint/calibrated provenance. Same-data staged/gauge fits
cannot be relabeled truth. Independent review precedes the already scoped
blind predictive closure. This parallel study clears none of that gate.

Separately BLOCKED: canonical-five Method A/final E.8/F.6.4 at
`f6_3_scientific_equivalence_mismatch:f4:$.parents[0].absolute_baseline_parent_sum`;
microscopic F.1/F.4 archaeology remains DEFERRED. The scoped absolute-unit
provenance limitation remains; no factor-of-1000 correction was demonstrated.
Use existing normalization for the authorized matched MM comparison: the
limitation is not a new standalone study or artificial prerequisite.
The [historical SIMC Work audit](investigations/simc-absolute-normalization-ledger-2026-11-03.md)
is retained with provenance, not adopted as a new prerequisite or demonstrated
normalization correction. The [November 3 accountability correction](investigations/kaonlt-assistant-roadmap-displacement-2026-11-03.md)
requires preserving mature functionality and the approved priorities; no repeated
audit or generator archaeology follows without concrete new defect evidence and
user authorization. Existing fitting/subtraction machinery is available for
the recorded bounded blind predictive closure once qualifying new evidence is
independently reviewed. Outcome B remains unresolved; the authorized parallel Part 1 action does not bypass it.
The earlier normal-launcher isolation blockers and matched-cohort ninth failure
(`external_input_verification/source_changed_during_hash`) remain in the
historical details and [incident record](investigations/kaonlt-matched-cohort-farm-readiness-failure-2026-10-22.md).
No fresh matched-event physical C or held-out fitted estimator is accepted;
no strength is selected. None of these separate scopes defines the present NEXT.

## Next Action

NEXT — independent actual-diff review of the November 8 diagnostic-only
scientific-identity failure receipt; after source acceptance, user publication
and exact pushed-state synchronization, consider only a separately reviewed,
evidence-preserving narrow diagnostic execution. No second full acquisition
is authorized. Local implementation base: test
`912ca31d0afdd7b9757395d205d7c54bc49018e8`.

The supplied November 7 failed-attempt fields in the
[contract](phases/check1-scientific-identity-failure-diagnostics-2026-11-08-task-contract.md)
identify `KaonLT_Check1_Acquisition_Q4p4W2p74_20261107T163552Z-3952780`
at that source with frozen reference SHA-256
`02ff71998fd51d5e56f07484bb1c71185f78ed51d6ce409b0d34154bda7a6d7a`.
Analysis return code was 0 and analysis completed, but `verify_artifacts` failed
at `check1:scientific_identity:Left:lowe.cells[0].fit.pion_component_alignment.applied_components.pi_delta.parent_alignment_hash`.
Artifact verification, collection and ZIP verification did not complete;
evidence is inadmissible and scientifically unaccepted. These are supplied
failed-evidence facts, not an independent local parsing of the farm outputs.
The destination collision did not recur; the new hash operands/cause remain
NOT VERIFIED. Scientific acquisition stays BLOCKED at the original invariant.

The diagnostic candidate retains the comparator, projection, coordinates,
tolerances and hash producer unchanged. It records bounded actual/reference
values, alignment selection/lineage and exact parent-hash operands in the
existing owner status before cleanup, then re-raises the original exception.
The November 8 review found that compact JSON did not bound the indented
receipt. The local repair now bounds the real nested UTF-8/text-newline
representation; its oversized fixture falls from 90,582 to 21,621 persisted
diagnostic bytes while retaining first-difference hashes and operand evidence.
The v2 review additionally found lost short window operands during compaction.
The local continuation compacts ancillary context first and retains both exact
first-difference records whenever they fit; the saved window regression is
21,690 bytes with both short window lists intact. Summarization of a first
operand is deferred until its own remaining budget is exceeded. Replacement
actual-diff review remains required; no source acceptance follows.
See the [implementation record](investigations/full-t-acceptance-edge-check1-2026-11-05.md#november-8-scientific-identity-failure-diagnostics).
No farm/ROOT/main/PDF or scientific validation is claimed. Prior P1/P2,
both Center conclusions, Outcome B, Method A/B and accepted closures stand.
Stop before Check 2.

## Success Criteria

Pion subtraction -> Lambda-window data/SIMC peak/tail shape agreement -> six
phi-ratio distributions -> fixed-data SIMC iteration toward unity -> existing
L/T, with physical background support and available errors/covariance. No new
threshold, algorithm, normalization or promotion; exact NEXT/gates unchanged.

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
