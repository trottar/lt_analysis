# Phase F.1 Method-A acceptance farm gate

## Current source review

SOURCE REVIEWED — the relevant F.1.Fix.5 analysis-source commit is
dc4fc6283001739a487ec80068f951b0e388cae6. The repository HEAD observed at the
last active-memory reconciliation (2026-09-11) was
23bcbd2f8db05f0a2e6cf399557c70df44326121; it is a timestamped observation, not
the permanent source-review identity. Later documentation/memory-only commits
do not by themselves invalidate the source review. No F.1.Fix.5
farm-evaluated commit is recorded.

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

ACTIVE — No fresh F.1.Fix.5 farm evidence was supplied. The current validation
profile requires analysis commit d656e15761970d7d612bb028d2746d077795e9ad and
permits only detached collector/profile changes after it. It therefore rejects
the later Fix.5 analysis source. See
investigations/2026-09-11-phase-f1-source-identity-reconciliation.md.

## Boundaries

F.1 is a detached diagnostic event contract. It must not alter random
subtraction, slow-proton treatment, pion subtraction, Method B, yields, or
final analysis. Method B stays diagnostic-only; F.1 does not license a future
production correction. Any F.2-F.6 work requires its own narrow contract.
