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

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING — detached global representation
audit at C2. It compares only `delta_only`, `track3`, `hgcer3`, and diagnostic
`full5_reference`; it neither chooses a five-dimensional estimator nor freezes
a basis. Farm review and explicit human basis acceptance are required before
F.3.

### F.3 — detached Method-A map

Requires support and out-of-domain handling. No Method-B numerical input.

### F.4 — detached parent-preserving A-only correction

Construct an A-only candidate while preserving canonical-t parent
normalization. Remain detached from production.

### F.5 — detached event-level (t,phi) propagation

Build parallel baseline and A-adjusted templates. Do not renormalize child
bins independently and do not promote to production.

### F.6 — explicit production promotion

Only after F.5 scientific validation and a separate explicit promotion
decision.

## NEXT

NEXT — Run F.2's one detached farm analyzer gate and explicitly accept or
reject its reduced-basis recommendation. Do not begin F.3 before acceptance.
