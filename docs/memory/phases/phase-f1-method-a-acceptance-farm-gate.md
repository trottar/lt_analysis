# Phase F.1 Method-A acceptance farm gate

## Scope and ownership

`SOURCE REVIEWED` — Current source provides a Method-A acceptance diagnostic
after Method A in `src/cuts/rand_sub.py`. It uses the frozen HGCer diagnostic,
Method-A result, Phase-A event contract, and pion-control cache to build a
separate acceptance event contract. The artifact is written with a
deterministic Method-A acceptance filename and its display payload is passed to
`src/cuts/full_background_subtraction_plots.py` for five setting-scope pages.

The source declares this path non-authoritative. Its unavailable fallback keeps
`production_objects_mutated`, `refinement_applied`,
`production_application_performed`, and `event_application_performed` false,
and declares no Method-B numerical dependency. This diagnostic must not change
random subtraction, slow-proton subtraction, pion subtraction, Method B, or
final analysis products.

## Evidence reviewed

- Current imports and post-Method-A invocation in `src/cuts/rand_sub.py`.
- Artifact construction/serialization in `src/cuts/rand_sub.py`.
- F.1 payload and terminal five-page renderer entries in
  `src/cuts/full_background_subtraction_plots.py`.
- Current contract implementation and focused source tests in
  `src/cuts/pion_hgcer_method_a_acceptance_contract.py`,
  `testing/test_pion_hgcer_method_a_acceptance_contract.py`, and
  `testing/test_pion_hgcer_phase_f_runtime_contract.py`.
- The detached profile/collector in `testing/pion_hgcer_validation_bundle_profile.json`
  and `testing/collect_pion_hgcer_validation_bundle.py`.

This review establishes source-path understanding only. It does not record a
test execution, a generated artifact, ROOT/PyROOT behavior, a rendered PDF, or
a farm outcome.

## Farm-gate boundary

`BLOCKED` — The checked-in profile requires analysis commit
`d656e15761970d7d612bb028d2746d077795e9ad` and allows only three detached
collector/profile paths after it. Live `test` contains later changes to F.1
analysis files, so its own identity guard rejects the current committed range.
See `../investigations/2026-09-11-phase-f1-source-identity-reconciliation.md`.

No `CLOSED / RUNTIME VALIDATED` claim is supported. A later narrow gate must
use fresh Phase-C and Phase-D checkpoints, parent-preserving correction,
Method-A acceptance artifact, procedure PDF/page manifest, and review of the
checker/provenance plus selected rendered PDF pages.

## Follow-up boundary

`NEXT` — Reconcile the reviewed analysis commit and the profile identity rule
in a narrow contract before any farm run. Do not use a profile archive or a
`complete=true` field as a substitute for raw-artifact and page review.
