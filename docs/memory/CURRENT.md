# Current KaonLT development state

Last source-identity and workflow-memory reconciliation: 2026-09-15,
accepted F.5 farm-bundle closure and F.5.2 presentation implementation. Repository HEAD values in this record are
timestamped observations, not permanent claims about the live checkout.

## Repository identity at last reconciliation

- Repository: trottar/lt_analysis; branch: test.
- F.2 exact starting HEAD:
  02f21886cff8a721df46ebf673594856583fc3c5.
- F.2 implementation / F.2.Fix.1 exact source starting HEAD:
  5549098d2552b9c092b65e31aeb77edc0807ddda. F.2.Fix.1 is a local
  source-only authority/provenance repair from that clean checkout; it is not
  a farm-evaluated commit.
- F.2 profile-driven collector exact source starting HEAD:
  170e6fae3d2fed1949fc6932b8eac9ad83e3e01c. This local source-only repair
  does not alter the F.2 analyzer or establish farm evidence.
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

CLOSED / RUNTIME VALIDATED — F.2.Fix.1 detached global five-setting Method-A
acceptance-representation audit. The accepted farm bundle ran at
`8e919fc618cea900227db5090d65da728c3aa555`, is complete with no unexpected
committed files, and preserves the five F.1 source hashes. The owner accepted
the unique supported reduced basis `hgcer3 = (SHMS_delta,
P_hgcer_xAtCer, P_hgcer_yAtCer)`. The F.2 historical artifact remains
`basis_frozen = false`; the human decision freezes this basis only for F.3.
See evidence/f2-fix1-runtime-closure.md.

CLOSED / RUNTIME VALIDATED — F.3 detached support-aware relative `hgcer3`
response map. The accepted farm bundle at
`5382cfc1994b078c620b32c043938134c33ffa39` has 15 valid parents, exact F.2
support continuity, no sparse/support failures, and no correction/production
behavior. See evidence/f3-runtime-closure.md.

CLOSED / RUNTIME VALIDATED — F.4 detached parent-preserving Method-A
correction. The accepted farm bundle ran at
`67e0298c51759c7a5ba693464d2c2655bf39250d`, has 15 valid parents and exact
F.3 support continuity, and retains parent-only signed normalization. It has
no template, yield, production, or Method-B ownership. See
evidence/f4-runtime-closure.md.

CLOSED / RUNTIME VALIDATED — F.5 detached signed `(t,phi)` propagation. The
accepted bundle ran at `3c6a66b7df9bf17e5a428458a2281a80831f001a` with 135
canonical cells (115 occupied, 20 explicit empty), exact F.1/F.4 propagation,
and all 15 parent/five-setting closures passing. Its scientific fingerprint is
`d11b728d1089301a12c29e7f8b1798c6e0b6021ac47bd5d2afc2b62b47a1effa`. See
evidence/f5-runtime-closure.md.

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING — F.5.2 is a presentation-only
PDF terminology rerender from accepted F.5 evidence. It changes no F.5 JSON,
scientific fingerprint, template arrays, authority, calculation, or production
behavior. SOURCE REVIEWED — F.5.Fix.1 identities remain distinct: F.5
implementation/Fix.1 start `f9ce5aadbebc31f5163a093ab32244644ad2a088`; Fix.1
implementation/source review `46598878c102a67275b1600e73cba3b2f166dc26`.

SOURCE REVIEWED — The existing validation collector now has a generic
profile-declared artifact path in addition to its unchanged F.1-specialized
path. `testing/pion_hgcer_validation_bundle_profile_f2.json` declares only the
global F.2 JSON/PDF and five setting-scoped frozen F.1 JSON inputs. It pins
`170e6fae3d2fed1949fc6932b8eac9ad83e3e01c`, permits only the four
collector/profile/test repair files plus `docs/memory/`, and performs no F.2
analysis or physics work.

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

NEXT — Preserve the accepted F.5 v1 evidence, then rerun only the detached
F.5 analyzer for F.5.2 presentation review against the accepted F.1/F.3/F.4
inputs. Run the unchanged generic collector only if a new review bundle is
needed. F.6 remains BLOCKED pending F.5.2 visual acceptance.
