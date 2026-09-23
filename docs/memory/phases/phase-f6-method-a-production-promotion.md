# Phase F.6 staged Method-A validation and production promotion

F.1 through F.6.2, including F.6.1 and F.6.2.Fix.5, remain `CLOSED / RUNTIME VALIDATED`. Their accepted evidence, F.5 scientific fingerprint `d11b728d1089301a12c29e7f8b1798c6e0b6021ac47bd5d2afc2b62b47a1effa`, and frozen F.6.2 JSON/fingerprints remain unchanged. No production promotion is authorized by those closures.

## Frozen baseline and Method-A boundary

The accepted full pion-subtraction baseline remains authoritative. For each authoritative pion-control event, `b_j^0 = s_j * w0_j`, where `s_j` is the established signed prompt/random/dummy coefficient and `w0_j` is the accepted missing-mass-dependent pion-background weight. The only Method-A variation is `b_j^A = s_j * w0_j * C_j`, using exactly the accepted F.4 parent-preserving correction.

This roadmap cannot redefine pion components, fits, windows, random/dummy subtraction, slow-proton treatment, canonical `(t,phi)` binning, F.3 map mathematics, F.4 parent normalization, SIMC, acceptance, efficiencies, L/T separation, cross sections, or production yields. It cannot independently normalize a child, clip/cap/winsorize factors, or introduce Method B numerically. Method B remains diagnostic/cross-check only; Method A remains detached until the explicit F.6.4 decision.

The accepted active background profile is `no_empirical_residual`; it forces
both empirical residual-background scales to zero. Legacy empirical residual
Fit 1/Fit 2 are dormant historical source machinery and are excluded from the
current E.8/F.6.3 chain.

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

F.6.2 scientific validation is `CLOSED / RUNTIME VALIDATED`; see
[the scientific closure record](../evidence/f6-2-scientific-runtime-closure.md).
The accepted detached validation assessed whether substantive Method-A
departures from the baseline are physically explainable and statistically
supported in each canonical `(setting, t, phi)` child. The baseline remains the
guide for broad physical missing-mass structure; Method A is an
acceptance-dependent refinement, not a replacement. This phase treats low
response as a reference shape, permits pion leakage in kaon-bearing regions,
and remains detached from production.

The accepted scientific gate persists raw normalized-shape
discrepancy/refinement/alignment measurements, independent acceptance and
MM×acceptance diagnostics, fixed-kaon-window pion change, support and
diagnostic variance statistics, and paired bootstrap intervals. It constructs
no composite score, freezes no numerical case thresholds, introduces no Method
B numerical dependency, and claims no final-yield uncertainty reduction.

### F.6.2.Fix.5 — presentation-only rerender

`CLOSED / RUNTIME VALIDATED` — the accepted Fix.5 presentation bundle closed
the F.6.2 presentation gate without changing the accepted scientific JSON
SHA-256 or artifact/validation fingerprints. Its renderer source
`c88ed65cb18ba6a37358897292b77696016312d1` and bundle/profile commit
`b789d203e11f0927deb59ebcae9dc59fe8add4ae` are presentation-artifact
provenance only; the scientific source remains distinct. See the
[Fix.5 closure record](../evidence/f6-2-fix5-presentation-runtime-closure.md).
F.6.2 overall is therefore `CLOSED / RUNTIME VALIDATED`. This gate has no
scientific, yield, ROOT, Method-B, or production ownership.

## Accepted E.8.1 narrow runtime evidence

E.8.1.Fix.5 and E.8.1.Fix.6 are `CLOSED / RUNTIME VALIDATED` only for the fresh `Q4p4W2p74 / Left / lowe` persisted-overlay geometry and profile/provenance gate. The bundle/profile commit is `0ec29d4e1bb345eb37e8cca35b8b7e5cbe1b4d5b`, the required distinct Fix.5 analysis/procedure source is `53fd262b730af8f1254e411a38231aebeb6a1da3`, and the frozen F.6.2 authority is unchanged. The remaining E.8.1 canonical-five expansion is `DEFERRED` by user decision, not failed and not full E.8 closure.

## Approved dependency order

    E.8.2 baseline full-analysis audit
        -> E.8.3 detached Method-A reweighting audit
        -> F.6.3 parallel full procedure plus Method A
        -> E.8.4 production-impact audit
        -> final E.8 closure
        -> F.6.4 explicit production-promotion decision

## E.8.2 — baseline full-analysis stage audit

`NEXT` — E.8.2 exposes the complete existing baseline analysis as an authoritative before/component/after missing-mass and yield chain. It is presentation-only and does not change physics.

### E.8.2a — prompt/random subtraction

Show the prompt kaon missing-mass input, random contribution removed, after-random output, and a difference or overlay where useful.

### E.8.2b — dummy subtraction

Show the after-random input, normalized dummy contribution, after-dummy output, and a difference or overlay where useful.

### E.8.2c — slow-proton subtraction

Show the after-random/dummy kaon input, proton contamination estimate, proton-cleaned output, before/after comparison, and existing PID/weight diagnostics. `_process_yield_data_tree()` already applies accepted proton cleaning during initial `(t,phi)` filling; existing raw/random/dummy snapshots must never be relabeled as pre-proton. E.8.2 must persist or consume a true authoritative pre-proton snapshot, not reconstruct or guess one.

### E.8.2d — baseline pion treatment

Show the proton-cleaned kaon input, actual baseline pion background `B_pi^0`
using accepted `w0`, and baseline clean-kaon output
`K_0 = K_proton - B_pi^0`, with before/after and signed difference where
useful. The contribution remains `b_j^0 = s_j * w0_j`.

### E.8.2e — final baseline canonical `(t,phi)` missing mass

For every canonical t parent, show all nine canonical phi children with
explicit empty/invalid status. Each populated final panel exposes canonical-t
identity, phi index/range, final baseline clean-kaon `MM_0(t,phi)` used for
yield extraction, Lambda signal/integration window, final baseline extracted
yield, and statistical uncertainty. There is no empirical-residual stage
between baseline pion subtraction and this final spectrum.

### E.8.2f — baseline stage-yield audit

Expose existing authoritative per-cell yield progression: prompt, after random,
after dummy, after proton, after baseline pion/final baseline clean sample, and
authoritative final baseline `Y0(t,phi)`. A renderer never recalculates a yield
when its authoritative producer can persist or provide it.

## E.8.3 — detached Method-A reweighting audit

`BLOCKED` pending E.8.2 source-reviewed completion. It is detached/non-production and consumes accepted F.4/F.5/F.6.1/F.6.2 results only.

Method A operates on the accepted pion background after the proton-cleaned
baseline picture; it does not invoke legacy empirical residual fits.

### E.8.3a — reweighting operation

Show `b_j^0 = s_j * w0_j` to `b_j^A = s_j * w0_j * C_j`. Method A changes only the accepted pion event contribution by the accepted F.4 parent-preserving factor; no child-by-child renormalization is permitted.

### E.8.3b — baseline versus reweighted pion background

For every canonical t parent, show `B_pi^0(MM)`, `B_pi^A(MM)`, their overlay, `Delta B_pi = B_pi^A - B_pi^0`, and `R_pi = B_pi^A / B_pi^0` only where defined and presentation-safe. This direct baseline-versus-reweighted comparison is required, not a legend implication.

### E.8.3c — canonical redistribution and parent closure

For every canonical t, show all nine `B_pi,tphi^0`, `B_pi,tphi^A`, and `Delta B_pi,tphi` children. Explicitly display stored/accepted closure residual and relative residual for `sum_phi B_pi,tphi^A ~= sum_phi B_pi,tphi^0`; never independently normalize a child.

### E.8.3d — F.6.2 acceptance/HGCer explanation

Retain accepted explanatory diagnostics: L/B/A normalized shapes, `delta x xptar`, `delta x yptar`, missing-mass x acceptance maps, support/OOD, effective statistics, and kaon-window refinement metrics. They explain the observed `B_pi^0 -> B_pi^A` redistribution; they do not replace the direct before/after plots.

## F.6.3 — parallel full procedure plus Method A

`BLOCKED` pending source-reviewed E.8.2 and E.8.3, not final E.8 closure. F.6.3 alone constructs the actual parallel Method-A full-analysis branch. It preserves the accepted full yield calculation unchanged as baseline and changes exactly `w0_j -> w0_j * C_j` while filling the pion-subtraction template.

Both branches use the same proton-cleaned input. Random/dummy subtraction,
slow-proton treatment, pion component models/fits/windows/amplitudes except
that event-level factor, canonical binning, SIMC, ordinary current
non-empirical yield extraction, efficiencies, acceptance, L/T separation, and
cross-section formulas remain identical between branches. The Method-A-clean
kaon spectrum is produced directly after pion subtraction and carried through
the existing yield extraction. Legacy empirical residual Fit 1/Fit 2 remain
disabled: `BG_OPT_ACTIVE_PROFILE = "no_empirical_residual"` stays frozen unless
a separate future scientific contract changes it. F.6.3 must not activate,
rerun, tune, or compare those historical fits. Disabling Method A must reproduce
the unchanged baseline branch.

## E.8.4 — baseline-versus-Method-A production-impact audit

`BLOCKED` pending F.6.3. E.8.4 consumes, and never constructs, F.6.3 two-branch outputs.

### E.8.4a — actual pion-subtraction consequence

Using the same proton-cleaned input compare `B_pi^0`, `B_pi^A`,
`B_pi^A - B_pi^0`, baseline clean kaon `K_0`, Method-A clean kaon `K_A`, and
`K_A - K_0`.

### E.8.4b — final canonical `(t,phi)` missing-mass comparison

For every populated canonical cell show baseline `MM_0(t,phi)`, Method-A
`MM_A(t,phi)`, overlay, signed difference, and the identical signal/integration
window.

### E.8.4c — final yield comparison

For every canonical `(t,phi)` show `Y0(t,phi)`, `YA(t,phi)`, `DeltaY = YA - Y0`, and `DeltaY/Y0` where defined. Summarize versus phi for each t, versus t for each setting, and later across all five canonical Q4p4W2p74 settings. The Method-A minus baseline shift is a correction effect, not automatically a systematic uncertainty.

## Final E.8 closure and F.6.4

Final E.8 is `BLOCKED` pending E.8.4 and later runtime/visual validation. It
requires prompt/random, dummy, slow proton, baseline pion, final baseline
MM/yields and stage yields, detached reweighting, direct baseline-versus-
reweighted pion background with signed difference/ratio, canonical
redistribution/parent closure, F.6.2 explanatory diagnostics, actual F.6.3
baseline-versus-Method-A clean-kaon MM comparison, final `Y0` versus `YA`, and
absolute/fractional yield shifts.

F.6.4 remains `BLOCKED` pending completed F.6.3/E.8.4 production-impact evidence. It is the only phase that may decide whether Method A becomes production pion-background treatment; no detached validation, visualization, or source review promotes it automatically. A later promotion contract must define supported kinematics, enable/disable behavior, fail-closed authority, setting-wide atomicity, provenance, uncertainty treatment, downstream validation, and regression against the unchanged baseline.

## Uncertainty and forbidden shortcuts

This roadmap defines no Method-A systematic uncertainty. In particular,
`|Method-A result - baseline result|` is a correction effect, not automatically
its uncertainty; later production closure requires separate justification.

This roadmap does not authorize production weights/yields, Method-A promotion,
or alteration of accepted F.1-F.6.2 scientific evidence. F.6 remains detached.
Its active successor is the approved E.8.2 -> E.8.3 -> F.6.3 -> E.8.4 forward
sequence. Do not alter ROOT or production behavior.

## Farm milestone cadence

Memory reconciliation and coherent E.8.2/E.8.3/F.6.3/E.8.4 implementation may proceed through deterministic local checks and independent source review without a farm run after every edit. Farm validation remains mandatory for ROOT/PyROOT/full-runtime claims. After a coherent F.6.3 plus E.8.4 branch is source-reviewed, run one narrow `Q4p4W2p74 / Left / lowe` end-to-end farm gate, inspect fresh artifacts, then PASS or make one coherent repair; broaden only after that gate passes.
