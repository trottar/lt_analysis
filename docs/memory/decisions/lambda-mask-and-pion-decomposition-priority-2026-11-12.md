# Lambda exclusion adoption and pion-decomposition priority — November 12, 2026

## Authority and status boundaries

Authority: the user's explicit November 12 memory-synchronization/scientific-priority
instruction. This records scientific decisions only. Public/local test at checkpoint:
`77c44eb60192f38c6979e71a2ae0cd8a920ada76`; pending source/memory changes remain
local. User adoption does not certify source review, publication, production
integration, final yield acceptance or farm/runtime validation. Those remain
separate; no phase is promoted here. CURRENT owns the single active NEXT.

## Adopted exclusion procedure and reviewed cohort mask

The user adopts the four Lambda (t,phi) criteria:

1. At least 30 total experimental entries.
2. At least 20 signal-region entries.
3. Signed cleaned signal-region yield at least the equal-width left-region yield.
4. Signed cleaned signal-region yield at least the equal-width right-region yield.

The reviewed diagnostic cohort retains 67/135 physical cells. Counts and signed
normalized yields are distinct. Its signal [1.114,1.126), left [1.102,1.114) and
right [1.126,1.138) GeV windows use the retained native experimental edges and
half-open inclusion. Total counts are retained kaon nommcuts source rows; signal
counts are allcuts rows in the signal interval, including signed-source streams
and zero-weight records. These are the reviewed definitions, not net significance.
Every failed criterion and empty/weak/nonpositive/finite-MC support state remains.
No extra negative-yield or poor-data/SIMC-agreement exclusion is adopted.

Preserve authoritative signed data, unchanged setting-level normalization,
no bin-level charge correction, and identical retained experimental/SIMC cells.
No renormalizing the surviving cells or replacing excluded measurements by
0 +/-10000. Preserve source/cohort/mask identity; a later correction must not
silently retune the mask against agreement or replace this reviewed support.
The 67-cell count describes this cohort, not guaranteed acceptance for other data.

Mask evidence: `kaonlt_lambda_bin_exclusion_20261112_review_bundle.zip`, SHA-256
`cf18308e93c9bc1444b2b7d73f139a97a2aff9e66c8bfa500e2a304d9182d9e5`,
member `report/primary_135_cell_decisions.csv`. The SHA-256 of the sorted retained
ID list encoded as compact JSON is
`504c22d76d71473ea7635c2b4e2e42b65c49edb2df707d179693b8078bfce09d`.
The review receipt includes the exact member hash and all 67 IDs. Existing support
losses and remaining MM discrepancies are retained, not erased by user adoption.
Production integration and farm validation of this mask are pending.

## Iteration-zero interpretation

Neither unity nor phi-flatness of data/SIMC ratios is an acceptance requirement
at iteration zero. Inspect these ratios diagnostically; do not reject an otherwise
supported initial result solely for nonunity or angular variation. The SIMC
cross-section parameterization is iterated subsequently with the initial pion
subtraction fixed. Do not tune background, area-match or alter normalization to
force either property. Preserve physical MM/yield requirements and uncertainties;
this is not a declaration of purity or final cross-section acceptance.
This supersedes earlier wording treating immediate constant-ratio agreement as
a prerequisite, without rewriting historical numerical results.

## All-15-parent MM presentation requirement

Every experimental/SIMC parent MM comparison uses identical bin edges. Never
mix coarse and fine histograms in the same comparison. Produce separate
unit-normalized shape and original yield-normalized comparisons; display shape
normalization never changes authoritative yields. Only exact bin aggregation or
source-event-backed binning can provide a common grid: no assumed fractional-bin
contents, interpolation or invented within-bin distribution. Report unsupported
mapping honestly; no new acquisition or replot is authorized by this memory task.
Earlier mixed-grid plots remain historical evidence, not fulfillment of this
requirement. No existing plot, histogram or analysis source is changed here.

## Immediate scientific priority and planned sequence

The fitted zero kaon-side SIDIS component is physically concerning and requires
an acceptance-informed assessment, particularly Center/high t3. Use the established
pion-control sample and existing neutron/Delta/SIDIS channel SIMC decomposition.
Examine accepted phase-space support, channel shapes/overlap, component degeneracies,
physically supported constraints and available uncertainties. Distinguish pion-control
amplitudes from kaon-side amplitudes and fitted zero from proven physical absence.
Do not impose a positive SIDIS lower bound without evidence or repeat completed
fit/bootstrap/model-response investigations as a new prerequisite.

Preserve accepted source-aware work, the 1024-replica conditional covariance and
116 model responses within their existing scopes. Complete physical covariance,
external PID/cross terms and the formal independently constrained Outcome B gate
remain unresolved, not zero or resolved by this decision. Neither Outcome B nor
unrelated normalization/Method-A gates displace this user's immediate assessment;
no replacement evidence-acquisition exercise is commissioned.

Following decomposition assessment and any separately authorized correction,
plan the established Local-t pion-subtraction chain -> (t,phi) binning -> adopted
exclusion mask -> Lambda spectra and yields. Preserve existing SIMC iteration,
L/T machinery and random/dummy/slow-proton ownership. Method A remains detached,
Method B diagnostic/excluded. This records the intended sequence, not present
permission to implement a correction or execute farm analysis.

## Review boundary

This memory-only update stops at independent ChatGPT actual-diff review. User-only
commit/push and subsequent published-state synchronization remain required. The
substantive priority remains the decomposition assessment, not another memory
reconciliation phase. Existing six roadmap statuses/dependencies and scoped
source/runtime closures are preserved. No physics source, Codex implementation,
farm run, refit, new bootstrap, commit or push occurs in this task.
