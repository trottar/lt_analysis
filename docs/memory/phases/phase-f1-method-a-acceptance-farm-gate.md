# Phase F.1 Method-A acceptance farm gate

## Current source review

CLOSED / RUNTIME VALIDATED — the relevant F.1.Fix.5 analysis-source commit is
`dc4fc6283001739a487ec80068f951b0e388cae6`; the KaonLT runtime-evaluated
analysis HEAD is `126fa22c19bd29b9952f55b33ab43d59f9727ef6`; and the narrow
collector reconciliation commit is `31dd034d8404e317863bf0933c51253e1d3deeb8`.
The owner accepted the clean C1 collector result and artifact-hash continuity
with the supplied v4 runtime bundle. The later closure-memory commit is the
F.2 baseline only; it is not an F.1 runtime commit.

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

The source inspection above is distinct from the accepted farm/runtime closure.

## Farm status and source-identity gate

CLOSED / RUNTIME VALIDATED — The supplied five-setting v4 archive directly
validated every F.1 v2 artifact's hash, detached state, training/application
population ownership, fingerprints, and stored canonical-cell closure at
runtime HEAD `126fa22c19bd29b9952f55b33ab43d59f9727ef6`. C1 changed only the
committed-range whitespace check to the three explicit profile-owned
validation files; the global worktree check and committed-file identity audit
remain strict. The accepted clean C1 collector gate reconciled provenance
without a KaonLT/ROOT rerun. See
evidence/f1-fix5-v4-bundle-inspection-2026-09-14.md.

## Boundaries

F.1 is a detached diagnostic event contract. It must not alter random
subtraction, slow-proton treatment, pion subtraction, Method B, yields, or
final analysis. Method B stays diagnostic-only; F.1 does not license a future
production correction. Any F.2-F.6 work requires its own narrow contract.

## Next phase boundary

NEXT — F.2 may begin from the closure-memory baseline as a detached global
representation audit only. This closure does not authorize a map, correction,
event application, production mutation, or F.3 work.
