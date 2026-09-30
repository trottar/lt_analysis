# Memory-health and farm-gate operational-completeness hardening

## Status

`SOURCE REVIEWED` — independent ChatGPT actual-diff review of
`kaonlt_review(20260930-173355).diff` passed. The reviewed local candidate is
based on committed `test` HEAD `86590fa655512926f2e4d0c50bf12b57d5198da5`.
Codex-reported checks were not run by ChatGPT. This status does not establish
commit, push, farm execution, scientific runtime, or production validation.

## Scope

The [task contract](memory-health-operational-completeness-hardening-task-contract.md)
owns only repository memory, `tools/check_memory_health.py`, and its focused
test. CURRENT is compacted below 7 KiB while retaining schema 3 and one NEXT.
Strict `--fail-on-warning` makes warning-free health a task-final gate without
changing the existing thresholds or default diagnostic behavior. Codex and
ChatGPT must each report exact memory byte counts, warnings, and manifest state
after substantial gates. Unresolved warnings block the next phase, fix, or
farm gate unless the user approves a recorded maintenance exception.

Before a farm command, the tracked, source-reviewed, pushed, and pushed-state-
reviewed owners of input authority, production, checking, packaging,
invocation, and returned artifact must cover the full operation. Missing
multi-step orchestration blocks farm handoff; one reviewed CLI may own a
complete direct operation. See the [failure investigation](../investigations/f4-refresh2-validation1-operational-readiness-failure.md).

F.4.Refresh.2.Validation.1 stays `SOURCE REVIEWED` for its profile. Its farm
execution gate is `BLOCKED` until a separate tracked execution owner is built
and reviewed. F.4.Refresh.2/Fix.1 remain `SOURCE REVIEWED`; F.4.Refresh.1/Fix.1
retain their narrow `CLOSED / RUNTIME VALIDATED` detached comparator closure.
F.6.3/E.8.4 remain source reviewed and runtime blocked. No scientific/runtime
source or accepted authority changes here.
