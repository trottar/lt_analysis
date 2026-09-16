# Phase F.6 staged Method-A validation and production promotion

F.6.1 is `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; no production
promotion is authorized by this roadmap.
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

`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING` — the detached validator and
analyzer were initially implemented at
`883ddf70af9935cd90e85c732c2d42d8765c783d` and source reviewed, including
F.6.1.Fix.1, at `bfc4fe421f9fc9139a1992a0ec92e31aa101b86c`. Before any
production or final-yield mutation, farm evidence must show whether accepted
Method-A reweighting moves the physical pion-control population toward the
observed low-response pion population while preserving accepted F.4/F.5
closure. Source review passed on 2026-09-16: py_compile; the F.6.1
validator/analyzer; F.4 correction; F.5 propagation/analyzer; and available
F.1 acceptance-contract tests, plus the `4928e145..bfc4fe42` range whitespace
check. Farm/runtime validation is not established.

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

## E.8 — streamlined full-background-subtraction presentation

`BLOCKED` until F.6.1 is accepted. E.8 is presentation-only and follows F.6.1
so the streamlined full-background-subtraction PDF can consume validated
Method-A reweighting comparisons, rather than only F.5 aggregate
redistribution. It does not calculate or own Method-A science. It should tell
the background-subtraction story—random/dummy subtraction, slow-proton
removal, baseline pion determination, Method-A motivation, validated
reweighting, and later F.6.2 impact—while detailed validation and Method-B
machinery remain in their dedicated artifacts.

## F.6.2 — full baseline procedure versus full procedure plus Method A

`BLOCKED` until F.6.1 is accepted and E.8 is updated in this order. Preserve
the current full yield calculation unchanged as the baseline branch. The
parallel Method-A branch uses the same random/dummy subtraction, slow-proton
treatment, pion components/fits/windows/amplitudes, canonical binning, SIMC,
background-fit algorithms/configuration, efficiencies, and yield extraction.
Its only intended scientific change is `w0_j -> w0_j * C_j` while filling
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

## F.6.3 — explicit production-promotion decision

`BLOCKED` pending accepted F.6.1 and F.6.2 evidence. This is the only phase
that may decide whether Method A becomes production pion-background treatment;
no automatic promotion follows from existing artifacts. A later promotion
contract must define supported kinematics, enable/disable behavior, fail-closed
authority, setting-wide atomicity, provenance, uncertainty treatment,
downstream validation, and regression against the unchanged baseline.

## Uncertainty and forbidden shortcuts

This roadmap defines no Method-A systematic uncertainty. In particular,
`|Method-A result - baseline result|` is a correction effect, not automatically
its uncertainty; later production closure requires separate justification.

This roadmap does not authorize production weights/yields, Method-A promotion,
or alteration of accepted F.1-F.5.2 evidence. F.6.1 remains detached. Its
next authorized step is a detached validation bundle/profile gate followed by
one fresh farm validation run; do not begin E.8 or alter ROOT or production
behavior.
