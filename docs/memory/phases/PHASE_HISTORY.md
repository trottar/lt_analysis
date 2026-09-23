# KaonLT phase and fix history

This is a concise phase-oriented map, not a replacement for per-gate evidence or `CURRENT.md`.

## Pre-HGCer production program

Random/dummy correction and frozen binning precede slow-proton cleaning and pion component subtraction; those paths retain separate scientific owners. HGCer diagnostics consume the frozen pion baseline and never redefine production.

## HGCer sequence

- Phase A/A.1/A.2/A.3 — frozen baseline, host event contract, and provenance hardening before Method B. [HANDOFF]
- Phase B/B.1 — independent Method-A detector-response diagnostics; B.1 source anchor `18b06ec1aa8ba42859dd0980705ce8374d4720ab`. [HANDOFF]
- Phase C — independent Method-B same-canonical-t local closure. C.Fix.2 adaptive partition remains diagnostic; recorded closure is `CLOSED / RUNTIME VALIDATED` only for its five-setting scope.
- Phase D — frozen A/B comparison, `DEFERRED` for final runtime status because a complete final farm manifest was not recovered.
- Phase E — presentation-only stored-upstream values. E.3.Fix.2 is `CLOSED / RUNTIME VALIDATED` only for Q4p4W2p74 Left-low; other setting closures remain `DEFERRED`.

## Phase F and E.8

F.1 through F.6.2, including F.6.1 and the accepted F.6.2.Fix.5 rerender, are `CLOSED / RUNTIME VALIDATED` within their detached/frozen scopes. The accepted F.6.2 scientific JSON and fingerprints remain unchanged, and no closure promotes Method A.

E.8.1.Fix.5/Fix.6 are `CLOSED / RUNTIME VALIDATED` only for the fresh Q4p4W2p74 Left/lowe persisted-overlay geometry and profile/provenance gate. The remaining E.8.1 canonical-five expansion is `DEFERRED`, not failed or closed.

E.8 is `ACTIVE` as the complete presentation-only full-analysis audit. Its approved sequence is E.8.2 baseline stage audit, E.8.3 detached Method-A reweighting audit, F.6.3 parallel full procedure plus Method A, E.8.4 production-impact audit, final E.8 closure, and F.6.4 explicit promotion decision. F.6.3 is `BLOCKED` pending E.8.2/E.8.3 source-reviewed prerequisites, E.8.4 is `BLOCKED` pending F.6.3, final E.8 is `BLOCKED` pending E.8.4/runtime evidence, and F.6.4 is `BLOCKED` pending full production-impact evidence.

Legacy empirical residual Fit 1/Fit 2 remain disabled under the accepted
`no_empirical_residual` profile and are excluded from the E.8/F.6.3 forward
chain.
