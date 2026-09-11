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

SOURCE REVIEWED at live test. It preserves separate Method-A training and
physical application populations; no estimator, correction, event adjustment,
or production mutation. Fresh F.1.Fix.5 farm evidence is ACTIVE pending
source-identity reconciliation.

### F.2 — freeze probability-map representation

After F.1 support review, choose a supported representation; candidate reduced
bases may include (delta,xptar,yptar) or (delta,HGCer x,HGCer y). Do not
blindly choose a 5D histogram.

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

NEXT — Resolve F.1.Fix.5 source identity for the farm profile, then perform its
single fresh targeted farm gate before considering F.2.
