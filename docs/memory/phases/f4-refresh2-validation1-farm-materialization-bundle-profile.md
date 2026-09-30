# F.4.Refresh.2.Validation.1 — farm materialization bundle/profile

## Status

`SOURCE REVIEWED` — independent ChatGPT actual-diff/source-provenance review of
`kaonlt_review(20260930-152818).diff` passed for the generic v4 profile and
focused synthetic test. The required analysis source is the pushed, independently reviewed
F.4.Refresh.2 materializer at `141a3d04f9e5d07be21dba14e0e63212c3990bf1`.
No farm materialization or ZIP has been produced for this gate.
The profile's `SOURCE REVIEWED` status does not establish operational farm
readiness. The farm execution gate is `BLOCKED` until a separate tracked,
reviewed, pushed, and pushed-state-reviewed owner covers validated inputs ->
materialization -> verification -> collection -> returned ZIP. The generic
wrapper packages only; it does not run the materializer. See the
[operational-readiness investigation](../investigations/f4-refresh2-validation1-operational-readiness-failure.md).

## Scope and source boundary

The profile uses the unchanged generic v4 collector with the canonical-five
top-level setting inventory, no setting-scoped artifacts, and exactly five
required global JSONs: the reviewed F.4.Refresh.1 comparison-input copy,
candidate F.2, candidate F.3, candidate F.4, and the F.4.Refresh.2
materialization manifest. It allows committed changes after the required
materializer source only in the new profile and its test, plus `docs/memory/`.
The canonical-five top-level setting inventory is required by the existing
collector schema, while the setting-artifact declaration is empty. Focused
synthetic coverage shows missing/malformed required JSONs, unreviewed
materializer/scientific commits, and output-ZIP collisions fail closed; an
undeclared extra source file is not packaged.
The later collector bundle's `complete=true` will attest only to packaging
and source checks; the materialization manifest and candidate science require
separate evidence review.

The materializer, comparator, collector, accepted F.2/F.3/F.4 artifacts,
scientific/runtime source, F.5/F.6.3/E.8.4 authority, production physics, and
Method-A promotion remain unchanged. F.4.Refresh.2/Fix.1 stay `SOURCE REVIEWED`.
F.4.Refresh.1/Fix.1 remain `CLOSED / RUNTIME VALIDATED` only for their detached
comparison gate. F.6.3/E.8.4 remain `SOURCE REVIEWED` and runtime blocked;
final E.8 and F.6.4 remain `BLOCKED`.

## Next gate

This profile remains source reviewed. The separately missing execution owner
must pass the source-changing workflow before any farm command; the
[hardening task](memory-health-operational-completeness-hardening.md) is active.
Earlier Codex-reported local checks passed: profile 4, collector 28,
materializer 13, memory 35 tests, all with 0 skips; manifest write/check,
bootstrap, and `git diff --check` passed. The then-existing CURRENT soft-size
warning is now a workflow blocker to be repaired by the hardening task;
ChatGPT did not run those unit suites.
Source review cannot establish ROOT/PyROOT, farm materialization,
candidate-authority acceptance, or production promotion.
