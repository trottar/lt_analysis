# Phase F.6 staged Method-A validation and production promotion

F.6.1 is `CLOSED / RUNTIME VALIDATED`; no production promotion is authorized
by this roadmap.
F.1 through F.5.2 remain `CLOSED / RUNTIME VALIDATED`; in particular, accepted F.5 scientific fingerprint
`d11b728d1089301a12c29e7f8b1798c6e0b6021ac47bd5d2afc2b62b47a1effa` remains
frozen. Method B is diagnostic/cross-check only and is numerically excluded
from every F.6 correction.

## Frozen baseline and Method-A boundary

The current accepted pion subtraction remains authoritative. For each
authoritative pion-control event, its baseline signed contribution is
`b_j = s_j * w0_j`, where `s_j` is the established signed
prompt/random/dummy source coefficient and `w0_j` is the accepted
MM-dependent pion-background weight. F.6 does not replace or refit `w0_j`.

The only intended Method-A variation is `b_j^A = s_j * w0_j * C_j`, using
exactly the accepted F.4 parent-preserving correction. F.6 may not redefine
pion components, fits, windows, random/dummy or proton treatment, canonical
`(t,phi)` binning, F.3 map mathematics, F.4 parent normalization, SIMC,
acceptance, efficiencies, L/T separation, or cross-section formulas. It may
not independently normalize children, clip/cap/winsorize factors, or introduce
Method B numerically.

## F.6.1 — detached Method-A reweighting validation

`CLOSED / RUNTIME VALIDATED` — detached validation closed from accepted bundle
`KaonLT_PhaseF6_1_validation_Q4p4W2p74.zip` at farm/bundle HEAD
`ffb7d4dc251610f6034eac11038a9841e3ef9f58`, with reviewed science source
`bfc4fe421f9fc9139a1992a0ec92e31aa101b86c`. It exactly reproduces accepted
F.4/F.5 closure through five settings, 15 parents, and 135 cells. The
27-page review was manually accepted. Its parent-dependent shape outcome is
evidence for F.6.2, not a uniform-improvement requirement; see
evidence/f6-1-runtime-closure.md.

For every canonical setting and canonical-t parent, compare the low-response
shape reference `0 < P_hgcer_npeSum <= 2`, the baseline physical control
population `P_hgcer_npeSum > 2` weighted by `w0`, and the transiently
reweighted physical control population weighted by `w0 * C`. The low-response
population is not an absolute leakage probability. Shape normalization must be
explicit on every comparison, either unit-area or scaled to the preserved
baseline parent integral.

Required comparisons cover the accepted `hgcer3` variables `SHMS_delta`,
`P_hgcer_xAtCer`, and `P_hgcer_yAtCer`; independent validation variables
`SHMS_xptar`, `SHMS_yptar`, `phi`, and missing mass; and, where consistently
available in frozen F.1 records, `Q2`, `W`, `HMS_xptar`, `HMS_yptar`, and
`theta_cm`. No validation variable becomes a Method-A input. Compare actual
signed baseline and reweighted pion-background shapes in missing mass, phi,
canonical `(t,phi)`, and any other frozen yield observable needed to show the
event-level effect.

F.6.1 is detached and shadow-only: production histogram mutation, production
yield mutation, and production-correction promotion are all false. Application
identities must exactly match F.4 transient factor identities and, when joined,
the authoritative child cache. Missing, extra, or duplicate records fail the
validated population. Training/Part-1 producer `evt.ph_q` is frozen in radians;
physical application `phi_degrees` is frozen as `evt.ph_q * 180/pi`; F.6.1
explicitly closes this rad-to-degree source contract without inferring units.
Aggregate event reweighting must reproduce accepted F.5
within its existing floating-point closure tolerance. Persist aggregate
validation summaries and provenance only; obtain factors transiently from the
accepted F.4 calculation and do not persist a new authoritative event-factor
table. Its review PDF must be physics-readable and cover population/authority,
model variables, independent variables, signed background shapes, F.5
`(t,phi)` continuity, and closure/provenance. Its exact page layout belongs to
the later implementation contract.

### F.6.1.Validation.1 — detached farm bundle/profile gate

`CLOSED / RUNTIME VALIDATED` — the generic v4 profile
`testing/pion_hgcer_validation_bundle_profile_f6_1.json` reuses the unchanged
collector and packages only the F.6.1 JSON/PDF, accepted F.5/F.4/F.3 JSONs,
five F.1 acceptance-contract JSONs, and manifest/source provenance. It pins
reviewed F.6.1 source `bfc4fe421f9fc9139a1992a0ec92e31aa101b86c`; the only
allowed committed files after that source are the profile and
`testing/test_pion_hgcer_validation_bundle_profile_f6_1.py`, plus the
`docs/memory/` prefix. Later F.6.1 validator/analyzer/science changes fail this
source gate until separately reviewed. The profile does not duplicate F.6.1
physics or make `complete=true` scientific acceptance. The accepted complete
bundle and runtime closure are recorded in evidence/f6-1-runtime-closure.md.

## F.6.2 — acceptance-correlated Method-A refinement validation

`NEXT` — before any full-procedure comparison, assess whether substantive
Method-A departures from the baseline are physically explainable and
statistically supported in each canonical `(setting, t, phi)` child. The
baseline remains the guide for broad physical missing-mass structure; Method A
is an acceptance-dependent refinement, not a replacement. This phase treats
low response as a reference shape, permits pion leakage in kaon-bearing
regions, and remains detached from production.

Persist raw normalized-shape discrepancy/refinement/alignment measurements,
independent acceptance and MM×acceptance diagnostics, fixed-kaon-window pion
change, support and diagnostic variance statistics, and paired bootstrap
intervals. Do not construct a composite score, freeze numerical case
thresholds, tune after inspection, introduce Method B numerically, or claim a
final-yield uncertainty reduction. The complete pre-implementation contract is
decisions/f6-2-acceptance-refinement-measurement-contract.md.

## E.8 — streamlined full-background-subtraction presentation

`BLOCKED` pending accepted F.6.2. E.8 is presentation-only and follows F.6.2
so the streamlined full-background-subtraction PDF can consume frozen accepted
F.6.1/F.6.2 results. It does not calculate or own Method-A science. It should
tell the background-subtraction story—random/dummy subtraction, slow-proton
removal, baseline pion determination, Method-A motivation, validated
refinement, and later F.6.3 impact—while detailed validation and Method-B
machinery remain in their dedicated artifacts.

## F.6.3 — full baseline procedure versus full procedure plus Method A

`BLOCKED` pending accepted F.6.2 + E.8. Preserve the current full yield
calculation unchanged as the baseline branch. The parallel Method-A branch
uses the same random/dummy subtraction, slow-proton treatment, pion
components/fits/windows/amplitudes, canonical binning, SIMC, background-fit
algorithms/configuration, efficiencies, and yield extraction. Its only
intended scientific change is `w0_j -> w0_j * C_j` while filling
pion-subtraction templates. Downstream empirical missing-mass background fits
rerun normally on the altered post-pion spectrum with their algorithms and
configuration frozen.

For every canonical `(t,phi)` bin compare baseline and Method-A pion background
removed, kaon yield after pion subtraction, residual-background Fit 1 and Fit 2
results, final kaon yield, and absolute/fractional final-yield differences.
Summarize versus phi per t parent, t per setting, and all five canonical
Q4p4W2p74 settings. The baseline path must be reproducible and unchanged when
Method A is disabled. Farm progression is narrow: Left-lowe, Left-highe,
Center-lowe, Center-highe, then Right-highe; do not broaden beyond Left-lowe
until its complete chain through final yield passes.

## F.6.4 — explicit production-promotion decision

`BLOCKED` pending F.6.2/F.6.3 evidence. This is the only phase that may decide
whether Method A becomes production pion-background treatment; no automatic
promotion follows from existing artifacts. A later promotion contract must
define supported kinematics, enable/disable behavior, fail-closed authority,
setting-wide atomicity, provenance, uncertainty treatment, downstream
validation, and regression against the unchanged baseline.

## Uncertainty and forbidden shortcuts

This roadmap defines no Method-A systematic uncertainty. In particular,
`|Method-A result - baseline result|` is a correction effect, not automatically
its uncertainty; later production closure requires separate justification.

This roadmap does not authorize production weights/yields, Method-A promotion,
or alteration of accepted F.1-F.5.2 evidence. F.6 remains detached. Its next
authorized step is design and implementation of F.6.2 acceptance-correlated
Method-A refinement validation; do not begin E.8 or alter ROOT or production
behavior.
