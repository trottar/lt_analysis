# E.8.4 Fix.5.4 — current-lineage numerical identity and SIMC normalization audit

**Status:** `ACTIVE`.
Memory checkpoint recorded 2026-10-02 at exact `test` HEAD
`f9d70732290ea461096374ca1270b47452644991`.

## Ownership and authorization

This phase owns the planned current-lineage identity audit. The
[checkpoint contract](e8-4-fix5-4-post-farm-identity-audit-memory-checkpoint-task-contract.md)
authorizes memory changes only. No audit code, plotting change or farm run is
implemented or authorized here. A separate implementation contract will be
written after independent review and user push of this checkpoint, followed by
pushed-state synchronization.

## Trigger and evidence boundary

The supplied Fix.5 farm observation at the above pushed source reports fresh
97-page Left/lowe output, zero renderer failures, all ten new pages and prior
E.8.4 pages retained. The owner returned no final ZIP; accepted bundle closure
and complete owner integration are not established. The
[post-farm investigation](../investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md)
records exact artifact names/timestamps, input identities, representative
scalar values, findings and uncertainties.

Confirmed findings: E.8.3 historical accepted F.6.1 aggregate lineage differs
from current F.6.3/E.8.4 candidate lineage; blue baseline curves are obscured
by magenta Method-A curves. Neither observation establishes which numerical
object, if any, is wrong.

## Planned numerical closure

Scientific interpretation of fresh Fix.5 pages is `BLOCKED` pending:

1. Binwise `MM_A - MM_0 = -(B_pi_A - B_pi_0)` for every populated child,
   with exact geometry and defined floating-point tolerance.
2. Histogram/scalar closure: `Integral(MM_0) = Y0`, `Integral(MM_A) = YA`,
   `Integral(MM_A - MM_0) = YA - Y0`, using exact producer-owned integration.
3. Separate signed integral, positive-bin, negative-bin and absolute support
   diagnostics before attributing large fractions to cancellation.
4. Same-cell authoritative SIMC object, normalization and unit identity,
   including its relevant integration-window integral and provenance.

All four are **NOT YET VALIDATED**. Trace producer -> sidecar -> payload ->
renderer and exact source/artifact/fingerprint identity. Preserve all 27
canonical cells, explicit empty status, signed values and uncertainties;
implement only warranted invariants or narrow source repairs under the future
contract. Do not assume an object mismatch, cancellation or SIMC scaling error.

## Preserved science and dependencies

F.4.Refresh.2, historical F.6.1/F.6.2 and the accepted narrow F.6.3 mechanics
retain their closures. E.8.3 remains `SOURCE REVIEWED` for its historical
lineage; it is not the current branch explanation. Method A remains detached,
Method B diagnostic-only and numerically absent. Final canonical-five E.8 and
F.6.4 remain `BLOCKED`. No accepted authority or production physics changes.

Visualization-only repair follows numerical closure: independently visible
curves, clearer styles/legends, or a permitted support panel using validated
objects. No content or normalization change for appearance. The missing ZIP
is a separate post-render operational issue; diagnosing it alone does not
warrant another expensive analysis run.

Sequence: checkpoint -> ChatGPT PASS -> user commit/push -> pushed-state
review -> separate Fix.5.4 numerical contract -> Codex numerical implementation
-> ChatGPT actual-diff review -> user commit/push -> pushed-state review ->
visualization-only contract and Codex implementation from that reviewed pushed
numerical source -> ChatGPT actual-diff review -> user commit/push ->
pushed-state review -> one narrow Q4p4W2p74 / Left / lowe farm run -> fresh
scientific and visual evidence review. CURRENT owns the sole
ordinary NEXT. This checkpoint stops before all implementation and farm work.
