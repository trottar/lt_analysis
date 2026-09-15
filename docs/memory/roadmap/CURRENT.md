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

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING — construct one A-only candidate
per canonical-t parent, preserve the full signed parent sum, and remain
detached from templates and production.

### F.5 — detached event-level (t,phi) propagation

BLOCKED — build parallel baseline and A-adjusted templates only after F.4 is
accepted. Do not renormalize child bins independently or promote to production.

### F.6 — explicit production promotion

Only after F.5 scientific validation and a separate explicit promotion
decision.

## NEXT

NEXT — Run F.4's single detached farm analyzer/collector gate. F.5 remains
BLOCKED pending explicit F.4 review.
