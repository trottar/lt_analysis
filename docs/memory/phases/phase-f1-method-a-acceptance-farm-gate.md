# Phase F.1 Method-A acceptance farm gate

## Current source review

SOURCE REVIEWED — the relevant F.1.Fix.5 analysis-source commit is
dc4fc6283001739a487ec80068f951b0e388cae6. The repository HEAD observed during
the 2026-09-14 archive inspection was
126fa22c19bd29b9952f55b33ab43d59f9727ef6; it is a timestamped observation, not
the permanent source-review identity. The supplied archive records that same
candidate farm HEAD, but not accepted clean-source closure. Later
documentation/memory commits and detached validation infrastructure do not by
themselves invalidate the source review.

The reviewed path is:

    Part-1 pion response records + Phase-E acceptance records
    exact join by (source_label, entry_index)
    -> prompt/noRF/nommcuts and NPE>0 Method-A training records
    -> separate training summary and fingerprint
    Phase-A pion records + authoritative parent/child physical cache
    -> NPE>2 application records
    -> separate application summary and fingerprint
    -> F.1 artifact and presentation payload

Source inspection confirmed:

- training is selected from prompt noRF/nommcuts positive-response Part-1 rows,
  with low 0<NPE<=2 and control NPE>2;
- application is restricted to authoritative physical NPE>2 records;
- source verifies exact paired identities, geometry/provenance, Method-A
  cell/setting closure, and separate training/application fingerprints;
- F.1 pages 1-4 use training rows, while page 5 uses only application rows;
- the F.1 builder occurs after Method A and takes no Method-B numerical input;
- unavailable and available contracts declare non-authoritative state, no
  correction/estimator/weight adjustment, and no production/event mutation;
- the renderer and static runtime contract preserve this detached boundary.

The current review is source inspection only. It is not a test execution, ROOT
render, artifact generation, farm run, or runtime claim.

## Farm status and source-identity gate

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING — A supplied five-setting v4
farm archive records candidate HEAD
`126fa22c19bd29b9952f55b33ab43d59f9727ef6` and directly validates every F.1
v2 artifact's hash, detached state, training/application population ownership,
fingerprints, and stored canonical-cell closure. The collector repair excludes
only allowed `docs/memory/**` paths from the historical range whitespace check
while retaining strict analysis-source and path identity gates. Local collector
and frozen F.1 tests pass. Per user direction, re-run the collector only from a
clean farm checkout against the existing outputs; do not mark F.1 runtime
validated until its complete manifest is inspected. See
evidence/f1-fix5-v4-bundle-inspection-2026-09-14.md.

## Boundaries

F.1 is a detached diagnostic event contract. It must not alter random
subtraction, slow-proton treatment, pion subtraction, Method B, yields, or
final analysis. Method B stays diagnostic-only; F.1 does not license a future
production correction. Any F.2-F.6 work requires its own narrow contract.

## Immediate gate

NEXT — Run the repaired collector only from a clean farm checkout against the
existing five `Q4p4W2p74` outputs, then inspect its complete manifest. This
does not authorize a correction, production mutation, or F.2 work.
