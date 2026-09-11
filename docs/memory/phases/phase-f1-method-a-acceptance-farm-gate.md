# Phase F.1 Method-A acceptance farm gate

## Current source review

SOURCE REVIEWED — the relevant F.1.Fix.5 analysis-source commit is
dc4fc6283001739a487ec80068f951b0e388cae6. The repository HEAD observed at the
start of this workflow-memory reconciliation (2026-09-11) was
7f8fcc02cde86f93e4e99042ba10c621f89dd159; it is a timestamped observation, not
the permanent source-review identity. Later documentation/memory commits and
the detached validation infrastructure do not by themselves invalidate the
source review. No F.1.Fix.5 farm-evaluated commit is recorded.

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

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING — No fresh F.1.Fix.5 farm
evidence was supplied. The detached profile/bundle v4 requires reviewed source
commit dc4fc6283001739a487ec80068f951b0e388cae6. It permits exactly its
collector/profile/test files and docs/memory/ after that commit; tracked root
AGENTS.md, unexpected analysis, and unrelated test changes remain identity errors. The
collector now accepts only the actual v2 dual-population artifact and verifies
separate training/application ownership, summaries, and fingerprints. Local
checks passed; no F.1.Fix.5 farm artifact is present. See
investigations/2026-09-11-phase-f1-source-identity-reconciliation.md.

## Boundaries

F.1 is a detached diagnostic event contract. It must not alter random
subtraction, slow-proton treatment, pion subtraction, Method B, yields, or
final analysis. Method B stays diagnostic-only; F.1 does not license a future
production correction. Any F.2-F.6 work requires its own narrow contract.

## Immediate gate

NEXT — The user first runs one fresh targeted v2 gate for Q4p4W2p74 Left lowe,
using the completed v4 collector/profile infrastructure. Collect fresh
checker/validation JSON and bundle provenance; inspect the F.1 v2 acceptance
artifact, relevant rendered F.1 pages, and traceback/log if anything fails;
then establish PASS or one coherent repair. Only after a targeted Left-lowe
PASS may validation broaden to Left lowe, Left highe, Center lowe, Center highe,
and Right highe. This gate does not authorize a correction, production mutation,
or F.2 work.
