# KaonLT phase and fix history

This is a phase-oriented map, not a substitute for per-gate farm evidence.

## Pre-HGCer production program

Slow-proton cleaning and pion component subtraction have distinct scientific
owners. The former is an event-level PID-contamination treatment protected by a
setting-wide Lambda gate; the latter is the staged pion-control prediction.
Random/dummy correction and frozen binning precede both. HGCer diagnostics
consume the frozen pion baseline and must not redefine either production path.

## HGCer sequence

- Phase A, A.1, A.2, A.3 — frozen baseline/host event contract and provenance
  hardening before Method B. [HANDOFF]
- Phase B/B.1 — independent Method-A detector-response diagnostics; B.1
  source anchor 18b06ec1aa8ba42859dd0980705ce8374d4720ab. [HANDOFF]
- Phase C — independent Method-B same-canonical-t local closure. C.Fix.1
  produced the checkpoint round-trip rule; C.Fix.2 adaptive partition remains
  diagnostic. CLOSED / RUNTIME VALIDATED for the five recorded settings at the
  closure scope only; see evidence/phase-c-five-setting-closure.md.
- Phase D — frozen A/B comparison using legacy Method B. DEFERRED for final
  runtime status because a complete final farm manifest was not recovered.
- Phase E — presentation-only over stored upstream values. E.3.Fix.2 is
  CLOSED / RUNTIME VALIDATED only for Q4p4W2p74 Left-low; other four setting
  closures are DEFERRED. See evidence/e3-fix2-left-low-runtime.md.

## Phase F

F.1 is a detached Method-A acceptance event contract, not an estimator or
production correction. Fix.3 and Fix.4 farm repairs are limited mechanical
evidence, while the F.1 v1 data-population result is a rejected scientific
interpretation. F.1.Fix.5 separates Method-A response training from downstream
physical application. At live test it is SOURCE REVIEWED; fresh F.1.Fix.5 farm
evidence is absent.

Approved later sequence, without implied authorization to start it:

    F.2 representation freeze
    F.3 detached Method-A map
    F.4 detached parent-preserving A-only correction
    F.5 detached event-level (t,phi) propagation
    F.6 explicit production promotion after F.5 validation

