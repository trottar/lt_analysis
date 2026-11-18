---
memory_schema: 3
---
# Current KaonLT development state

## Active Objective

`ACTIVE`: physically credible clean K+Lambda: close cleaned-data/K+Lambda-SIMC
shape agreement in the actual Lambda MM window (peak/relevant tail), matched
selection/binning/resolution/acceptance and documented normalization. Integrated
yields alone are insufficient; close shape agreement is necessary, not proof
of purity. Neither unity nor phi-flatness of data/SIMC is an acceptance
requirement at iteration zero. Later SIMC cross-section iterations address
ratio shape/magnitude with the initial pion subtraction fixed; data itself
retains physical angular dependence. Preserve existing SIMC and L/T machinery.
The completed D2.a constraint study awaits independent scientific review.
Source-kernel NumPy evidence is not ROOT or production acceptance.

## Current Work Item

November 12 user decision: adopt the four-criterion Lambda (t,phi) exclusion
procedure: >=30 total and >=20 signal-region raw entries, signed cleaned
signal yield >= each equal-width neighboring yield. The reviewed cohort mask
retains 67/135 physical cells; user adoption is distinct from pending production
integration and farm validation. Preserve signed observations, unchanged
setting-level normalization and identical retained data/SIMC support. Exact
criteria, windows, mask identity and ownership are in the
[governing decision](decisions/lambda-mask-and-pion-decomposition-priority-2026-11-12.md).

D1 and D1.a are SOURCE REVIEWED within their bounded scopes, per the November
16 user closure handoff. The [D1/D1.a record](phases/kaon-sidis-empty-bin-support-2026-11-16.md)
owns exact identities, findings and review attribution. The empty-bin defect is
established: do not repeat its diagnosis or a general D1 covariance study.
Full covariance and physical SIDIS absence remain unverified.

D1.a changes only its SIDIS fitting objective, retaining physical bins and
signed cancellations without a positive floor. Its reviewed
[contract/candidate](phases/kaon-sidis-empty-bin-support-2026-11-16-task-contract.md)
remain unchanged and default-off; source publication is verified at test
`e1f02933774d1c066154be49a77aa66340593bc3`, without farm validation.

D2/D2.a remain SOURCE REVIEWED, not corrected ROOT accepted; the bounded
mechanism is published at `55a1258791483b658729c13508c9aa2c49864e74`.
The detached acceptance/constraint study is numerically complete at source
`0ae7c5f6b0673df3345977c7d1568eccc3cc4ddc`; see the
[measured assessment](phases/event-population-authority-2026-11-16.md#november-17-completed-acceptance-and-model-dependence-sensitivity).
Policy-only passes6/6/7/8;150 central and19200 paired protected refits.
Continuous constraints are stable (max parent span0.04649%) without broad
variance gain. Window1.27 gives10 passes but up to14.31% parent/15.53% retained
cell yield response and unresolved omitted-region residuals. No better physical
full-cohort tradeoff or baseline adoption is established. STRICT unchanged;
FORCED diagnostic;fixed67,Local-t,six steps,D3-D5 DEFERRED,OutcomeB BLOCKED.

Planned after decomposition assessment and any separately authorized correction:
rerun the established Local-t chain through (t,phi) binning and the adopted
mask, then evaluate Lambda spectra/yields using unchanged SIMC and L/T machinery.
This records the sequence; it does not execute or authorize a farm run now.
The prior physical and exclusion assessments, their support losses and unresolved
physical/covariance limitations remain preserved in the
[Check 1 record](investigations/full-t-acceptance-edge-check1-2026-11-05.md).

Q4p4W2p74 Local-t remains the selected scientific baseline. The November 5
meeting identified Center/low and Center/high t3 as the principal Lambda-MM
exceptions. The November 9 Check 1 packet supplies continuous -t, occupancy
and matched-interior comparisons. Both Center settings have observed support
through nearly all canonical t3; this does not establish geometric acceptance
or a resolved Lambda signal. The earlier low-epsilon lower-t3 tail interpretation
and its use against an upper-edge explanation are withdrawn. Low-epsilon t3
has weak, poorly resolved Lambda signal; high-epsilon t3 has a potentially
substantial applied pion background of unverified physical accuracy. See the
[Check 1 record](investigations/full-t-acceptance-edge-check1-2026-11-05.md#november-9-continuous-t-evidence-and-pid-follow-up).

Check 1.B reproduced Center/low t2 broad-MM subtraction: the maximum transfer
was 3435.63 from a small fitted pion-control denominator. This is the known
denominator-amplification mechanism, with no cross-cohort regression established.
SHMS delta and HGCer position vary with -t. Low-positive HGCer fractions do not
show a common Center upper-t3 increase. Center/high t3 applied pion background
is concentrated at negative delta; the experimental source split does not
identify the pion-production channels. See the
[Check 1.B result](investigations/full-t-acceptance-edge-check1-2026-11-05.md#november-9-check-1b-results-and-check-1c-decision).

Check 1.C measured different positive-response trends at the two Center
settings on common joint detector support. Kaon-selected rows support kinematic
comparison; matched kaon detector values were not available in the reviewed
evidence. The retained channel SIMC shapes differ in MM and phi, with substantial
overlap beneath Lambda. The current fitted Center/high t3 SIDIS term contributes
to the pion-control denominator while its kaon-side amplitude is zero. Fitted
allocation cannot identify experimental species or establish channel delta
occupancy. See the [Check 1.C record](investigations/full-t-acceptance-edge-check1-2026-11-05.md#november-10-ten-bin-channel-constraint).

The November 5 sequence is historical context. The November 12 user decision
now governs mask adoption, the acceptance-informed pion-decomposition priority
and subsequent separately authorized reanalysis; no roadmap status is promoted.
Preserve separate unit-normalized shape and existing yield-normalized magnitude
diagnostics. No arbitrary scaling of yield spectra or repeated normalization.
Iteration-zero SIMC uses the old parameterization; initial D/M need not be unity.
Method A stays detached, optional pending a justified target; Method B remains
diagnostic/excluded. Outcome B and source-aware statistics remain separate.

The November 11 regenerated pion h10 projection is CLOSED / RUNTIME VALIDATED
within its accepted [15-source/45-parent diagnostic-acquisition boundary](evidence/regenerated-pion-h10-projection-runtime-closure-2026-11-11.md).
The detached joint MM/delta fit, source-aware 1024-replica bootstrap and
five-setting/15-parent uncertainty propagation with 116 coherent model cases
are SOURCE REVIEWED per the supplied Outcome B handoff. Publication and
independent pushed-state synchronization are recorded at test
d28ad6e7f620f0c4cc18c389373560f56fd9e84a; this checkpoint verifies live public/
local identity and reviewed source snapshots. Amplitudes, covariance, responses
and central predictions are unchanged. Conditional statistical covariance and
model sensitivity ranges are available; complete calibrated physical covariance
and external PID/cross terms remain unavailable, not zero.

Outcome B remains BLOCKED: current supplied source-bound evidence contains no
qualifying independent Center/high t3 joint anchor-component constraint. The
same-data staged/gauge fits, fixed-protected-A bootstrap and propagated model
responses do not provide independent kaon-side means or their count covariance.
The current [Outcome B audit](investigations/full-t-acceptance-edge-check1-2026-11-05.md#november-11-outcome-b-sourceevidence-audit-and-publication-synchronization)
records the exact missing input. The present assessment does not resolve that
formal limitation or commission another independent-evidence acquisition.

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

At the November 9 checkpoint, local/origin/live test and farm acquisition source matched
`6fb1a94ed8f6b7af61e130c9b42b9762164444b9`. The published deterministic
identity correction matches the seven independently reviewed blobs. The owner
reports success, including artifact verification and publication. Direct
readback verified the archive, log and summary hashes, 35 artifacts, both
scientific comparison paths and all ten numerical packets against the pinned
reference. These are scoped diagnostic acquisition checks. The stored
`evidence_admissible=false` and `scientifically_accepted=false` remain.
Production physics remains unchanged from scientific source
`9172aba853fd5b5b744c19add992fa4f27466180`.

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
PID/acceptance diagnostic using available data; no producer rebuild or farm run follows.

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

NEXT — Independent ChatGPT review of the completed numerical evidence and
actual four-file memory diff, followed by a user iteration-zero baseline decision.
No parameter,threshold or hybrid adoption follows automatically. Method A remains
separate; no additional acquisition,farm or Codex task is commissioned.

## Success Criteria

Physically supported pion decomposition/subtraction -> adopted bin exclusion
with matched retained data/SIMC support -> Lambda MM/producer-yield assessment
on identical bin edges -> existing SIMC iteration with initial subtraction fixed
-> existing L/T. Neither unity nor phi-flatness is an iteration-zero acceptance
gate. Preserve formal Outcome B, uncertainty limits and source/runtime boundaries.

## Do Not Reopen Without New Evidence

Except for the user-adopted exclusion procedure awaiting separate integration,
preserve `no_empirical_residual`, zero legacy scales and all current cuts,
templates, priors, binning, alignment, normalizations, weights, subtraction,
efficiencies, acceptance, yields/errors, L/T and cross sections. No child
renormalization. `YA/MM_A` is the reweighted-baseline branch; standalone
Method-A-only final yield/MM is NOT VERIFIED. Frozen reviewed F.4 SHA-256:
`79e7ceda7221cbeeead4ed5bc306b0e0e670741a27beaa980e22349c555e96d7`.

## Relevant References

The [investigation](investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md)
owns numerical evidence/history; [MEMORY](MEMORY.md) owns reusable knowledge;
the [roadmap](roadmap/STATUS.md) retains the interrupted approved progression.
