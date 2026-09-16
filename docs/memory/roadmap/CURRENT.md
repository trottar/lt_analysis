# Approved KaonLT roadmap

This roadmap records approved Phase-F structure. It does not make a runtime
claim, authorize source changes, or promote a diagnostic into production.

## Frozen upstream program

Random subtraction, slow-proton PID cleaning, pion component subtraction,
Method-A diagnostics, Method-B closure, Phase-D comparison, and Phase-E
presentation retain separate ownership. Phase E consumes frozen records only.
Method B remains diagnostic/cross-check/historical comparison only; it never
numerically adjusts pion event weights.

## Phase F

### F.1 — detached Method-A acceptance event contract

CLOSED / RUNTIME VALIDATED for all five canonical Q4p4W2p74 settings. It
preserves separate Method-A training and physical application populations; no
estimator, correction, event adjustment, or production mutation. Runtime is
anchored to `126fa22c19bd29b9952f55b33ab43d59f9727ef6`; the later collector
reconciliation does not replace that runtime identity.

### F.2 — freeze probability-map representation

CLOSED / RUNTIME VALIDATED — detached global representation audit. The
accepted F.2 farm bundle ran at `8e919fc618cea900227db5090d65da728c3aa555`.
The owner accepted `hgcer3` as the unique supported reduced basis after all 15
groups passed its frozen information and support gates. `full5_reference`
remains diagnostic-only. See evidence/f2-fix1-runtime-closure.md.

### F.3 — detached Method-A map

CLOSED / RUNTIME VALIDATED — detached support-aware relative `hgcer3` map,
accepted at `5382cfc1994b078c620b32c043938134c33ffa39` with all 15 support
continuity gates passing. See evidence/f3-runtime-closure.md.

### F.4 — detached parent-preserving A-only correction

CLOSED / RUNTIME VALIDATED — accepted at
`67e0298c51759c7a5ba693464d2c2655bf39250d`; all 15 parent corrections passed
F.3-support and closure review. It remains detached from templates and
production. See evidence/f4-runtime-closure.md.

### F.5 — detached event-level (t,phi) propagation

CLOSED / RUNTIME VALIDATED — accepted at
`3c6a66b7df9bf17e5a428458a2281a80831f001a`; the detached signed aggregate
templates have scientific fingerprint
`d11b728d1089301a12c29e7f8b1798c6e0b6021ac47bd5d2afc2b62b47a1effa`. F.5
does not renormalize child bins independently or promote to production. See
evidence/f5-runtime-closure.md.

### F.5.2 — detached presentation-only F.5 review rerender

CLOSED / RUNTIME VALIDATED — the presentation-only rerender was accepted at
`6634e9cb470cf35f21f5d475ec6ce33b524cd233`. Its 12-page review uses physical
phi-interval labels while the accepted F.5 payload and scientific fingerprint
remain exactly unchanged. See evidence/f5-2-runtime-closure.md.

### F.6 — staged Method-A validation and explicit production promotion

ACTIVE — F.6.1 follows accepted F.5/F.5.2 under the detailed contract in
phases/phase-f6-method-a-production-promotion.md. The current full procedure
remains the accepted baseline until an explicit later promotion decision.

#### F.6.1 — detached Method-A reweighting validation

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING — the detached validator and
analyzer are source reviewed through F.6.1.Fix.1 at
`bfc4fe421f9fc9139a1992a0ec92e31aa101b86c`; initial implementation was
`883ddf70af9935cd90e85c732c2d42d8765c783d`. It validates accepted F.4
factors against physical pion-control populations without mutating production
weights or yields. The frozen producer phi contract is training `evt.ph_q` in
radians and application `phi_degrees = evt.ph_q * 180/pi`; farm/runtime
validation is not yet established.

#### F.6.1.Validation.1 — detached farm bundle/profile gate

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING — generic v4 profile
`testing/pion_hgcer_validation_bundle_profile_f6_1.json` reuses the unchanged
collector to package only F.6.1 JSON/PDF, accepted F.5/F.4/F.3 JSONs, five F.1
acceptance-contract JSONs, and source provenance. It pins
`bfc4fe421f9fc9139a1992a0ec92e31aa101b86c` and permits only the profile/test
after that source, plus `docs/memory/`; it has no F.6.1 semantic or production
ownership. Farm/runtime validation is not yet established.

#### E.8 — streamlined full-background-subtraction presentation

BLOCKED — presentation-only update after accepted F.6.1; it consumes frozen
validated reweighting results and does not recompute Method-A science.

#### F.6.2 — full baseline procedure versus full procedure plus Method A

BLOCKED — after F.6.1 and E.8, compare the unchanged full baseline procedure
with a parallel procedure whose only intended pion-template change is the
accepted Method-A factor.

#### F.6.3 — explicit production-promotion decision

BLOCKED — only after accepted F.6.1/F.6.2 evidence; no automatic promotion.

## NEXT

NEXT — One fresh detached F.6.1 farm validation for Q4p4W2p74, then inspect
the complete F.6.1 JSON/PDF/bundle evidence. E.8, F.6.2, and F.6.3 remain
BLOCKED.
