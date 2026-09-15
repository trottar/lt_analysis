# Current KaonLT development state

Last source-identity and workflow-memory reconciliation: 2026-09-14,
F.2.Fix.1 source implementation. Repository HEAD values in this record are
timestamped observations, not permanent claims about the live checkout.

## Repository identity at last reconciliation

- Repository: trottar/lt_analysis; branch: test.
- F.2 exact starting HEAD:
  02f21886cff8a721df46ebf673594856583fc3c5.
- F.2 implementation / F.2.Fix.1 exact source starting HEAD:
  5549098d2552b9c092b65e31aeb77edc0807ddda. F.2.Fix.1 is a local
  source-only authority/provenance repair from that clean checkout; it is not
  a farm-evaluated commit.
- Reviewed F.1.Fix.5 analysis-source commit:
  dc4fc6283001739a487ec80068f951b0e388cae6.
- F.1 runtime-evaluated analysis HEAD:
  126fa22c19bd29b9952f55b33ab43d59f9727ef6.
- F.1.Validation.Fix.1 collector commit:
  31dd034d8404e317863bf0933c51253e1d3deeb8.
- This docs/memory-only closure commit is the F.2 starting baseline. It is
  neither the F.1 analysis-source commit nor a KaonLT runtime commit.

## Supported state

CLOSED / RUNTIME VALIDATED — Phase-C Method-B diagnostic closure for Left
lowe, Left highe, Center lowe, Center highe, and Right highe, at C.Fix.2.3
9a66bc62d20a99172e326e915866877b65ae1e5d and later accepted pre-E.3
e3853655db0809923cbf2326e2f779219128eda9. Adaptive Method B is DO NOT
PROMOTE. See evidence/phase-c-five-setting-closure.md.

CLOSED / RUNTIME VALIDATED — E.3.Fix.2 independent Method-A presentation for
Q4p4W2p74 Left-low only: implementation
eb1710f4739ba6ef14f51419806e9fc5bd53c175 and runtime/bundle HEAD
bf53dac84e1396cfd7e3f4e0234749426bfbdcf4. See
evidence/e3-fix2-left-low-runtime.md.

CLOSED / RUNTIME VALIDATED — F.1.Fix.5 for Left-lowe, Left-highe,
Center-lowe, Center-highe, and Right-highe. The reviewed analysis source is
`dc4fc6283001739a487ec80068f951b0e388cae6`; runtime remains anchored to
`126fa22c19bd29b9952f55b33ab43d59f9727ef6`. The accepted clean post-hoc
collector reconciliation ran from C1 `31dd034d8404e317863bf0933c51253e1d3deeb8`
against the existing artifacts. The owner accepted the required complete source
gate and artifact-hash continuity with the supplied v4 runtime bundle. F.1 is
detached, Method-B independent, and makes no correction, estimator, or
production mutation. See evidence/f1-fix5-v4-bundle-inspection-2026-09-14.md.

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING — F.2 is a detached global
five-setting Method-A acceptance-representation audit. Its F.2.Fix.1 source
repair starts at `5549098d2552b9c092b65e31aeb77edc0807ddda` and reconstructs
the frozen F.1 v2 population/metadata/child/full fingerprints before any
metric, requires exact F.1 authority metadata and geometry, binds deterministic
filenames to payload identity before source hashes, and fails safely for output
path collisions and sparse zero-row non-prompt groups. The fixed probe,
candidates, thresholds, cKDTree criterion, and recommendation policy remain
unchanged. It writes one aggregate JSON plus one four-page matplotlib PDF, with
no ROOT import, map, correction, event probability, weight, normalization,
Method-B numerical input, or production side effect. Local F.2 and frozen F.1
tests pass. See phases/phase-f2-method-a-acceptance-representation.md.

## Frozen architecture

- Random subtraction, slow-proton treatment, pion subtraction, HGCer Method A,
  HGCer Method B, SIMC comparison, yields, cross sections, diagnostics/checkers,
  and presentation retain separate ownership.
- Proton cleaning precedes pion subtraction; K Lambda preservation remains a
  setting-wide production gate. Proposed proton quantities never become applied
  quantities without gate acceptance.
- Method A and Method B are independent. Method B is same-canonical-t relative
  closure and remains diagnostic/cross-check only; it never adjusts pion
  weights. Method A is positive-response relative leakage only, not absolute
  zero-NPE probability.
- Phase E presentation uses frozen upstream records and has no physics
  ownership. F.1 through F.5 remain detached; only an explicit validated F.6
  may promote a production change. Parent-level normalization forbids
  independent (t,phi)-child renormalization.

## Explicit gaps

DEFERRED — final Phase-D farm provenance; final four-setting E.3 closure;
detailed older proton farm artifacts; exact raw paths/hashes for recovered
historical bundles. See investigations/KNOWN_GAPS.md.

## NEXT

NEXT — Run the detached F.2 analyzer once on the farm against the accepted five
Q4p4W2p74 F.1 v2 artifacts. Inspect JSON/PDF provenance, every per-group
response/support gate, OOD metrics, and recommendation. Accept one reduced
basis or make one coherent F.2 repair; do not begin F.3 without explicit basis
acceptance.
