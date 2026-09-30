# F.4.Refresh.2.Validation.1 operational-readiness failure

## Observed boundary

Validation.1 correctly source reviewed a generic bundle profile and focused
test, not a complete farm execution driver. ChatGPT then incorrectly treated
the reviewed profile and pushed state as sufficient to begin interactive farm
orchestration. An attempted invocation in a stale ordinary farm checkout failed
before Python could open `testing/materialize_method_a_current_baseline_authority.py`.
The materializer exists in pushed source; checkout staleness was observed, but
it was not the repository-design root cause. No candidate materialization,
returned ZIP, or accepted-authority mutation resulted from that attempt.

## Workflow defect

The full path from validated inputs through materialization, verification,
collection, and returned review artifact was not audited before farm handoff.
The existing generic `testing/package_pion_hgcer_validation_bundle.tcsh`
wrapper is bundle-only and cannot own candidate materialization. No tracked,
reviewed F.4.Refresh.2 operation owns the complete multi-step sequence. The
profile stays `SOURCE REVIEWED`; its farm execution gate is `BLOCKED` pending
a separately implemented and reviewed execution owner.

Repeated CURRENT soft-size warnings also persisted although
[MAINTENANCE.md](../MAINTENANCE.md) requires focused consolidation. Warning-only
CLI success was incorrectly allowed to stand in for a healthy next-gate handoff.

## Prevention

Before any farm command, audit input authority -> producer -> checker ->
collector -> invocation owner -> returned artifact. All required executable
steps must be tracked, locally deterministic where possible, independently
source reviewed, pushed, and pushed-state reviewed. Missing multi-step
orchestration returns to a source-changing workflow. Strict memory health and
the required Codex/ChatGPT health report block new gates while warnings remain,
unless the user approves a recorded maintenance exception. See the
[hardening phase](../phases/memory-health-operational-completeness-hardening.md).
