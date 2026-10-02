# KaonLT farm communication

This record owns how a requested farm operation is communicated and what is
returned for independent review. See [TOOLS.md](TOOLS.md) for canonical shell
and path facts, and the [farm validation-bundle procedure](decisions/farm-validation-bundle-procedure.md)
for detailed packaging.

## Request and return pattern

Before providing a farm command, audit input authority -> producer/materializer/
analyzer/renderer (if any) -> verification/checker (if any) -> collector/
packager -> invocation owner -> expected returned artifact. Every required
executable step must be tracked, locally deterministic where possible,
independently source reviewed, pushed, and pushed-state reviewed. State
`Farm readiness: PASS` and name the tracked source owners of the complete
operation before a command. If an owner is missing or unreviewed, state
`Farm readiness: BLOCKED` and provide no command. Multi-step orchestration
needs a tracked, reviewed driver; do not assemble it interactively in the
farm shell. One reviewed CLI can directly own a complete single-command
operation. Hard failures and material active-state or provenance ambiguity
block farm readiness: current source identity, accepted evidence/status, frozen
scientific interfaces, active scientific ownership/blocker, or exact NEXT.
A nonblocking memory warning alone does not block an otherwise ready farm
gate; record it and batch correction at the next checkpoint or milestone audit.

Farm commands must be valid `tcsh`. Communicate a simple operation as:

1. the exact requested operation and artifact;
2. one concise caveat when needed;
3. one exact `tcsh` command block; and
4. the requested failure output or returned artifact/ZIP.

Return fresh artifacts or a fresh bundle for independent review. On failure,
return the complete command output needed to diagnose the requested gate.

## Requested-operation boundary

Identify the request before issuing a farm command:

- **Full scientific rerun:** use the applicable scientific contract.
- **Presentation-only rerender plus bundle:** retain the frozen scientific
  payload and rerender only the presentation output required by the gate.
- **Bundle-only:** package the specified existing artifacts only.

Do not silently escalate a bundle-only request into an analyzer rerun. Do not
disturb unrelated modified or untracked farm-local files merely to package
detached evidence.

## Evidence boundary

Successful ZIP creation or `complete=true` is not evidence acceptance. Review
provenance, checker results, payloads or logs, and rendered pages as applicable.
Do not claim farm, ROOT/PyROOT, full-runtime, or production validation without
direct supplied evidence.

Source-changing task workflow belongs to [CODEX.md](CODEX.md); this record
does not replace it, the command catalog, scientific contracts, or detailed
farm packaging procedure.

## Failed gates and verification before consumption

Before readiness or a farm command, perform the full health check in
[MAINTENANCE.md](MAINTENANCE.md), state `Farm readiness: PASS` or
`Farm readiness: BLOCKED`, and identify the supporting evidence and tracked
source owners. A previous successful run or a PDF/ZIP path alone proves no
readiness.

If an owner/checker/provenance/manifest/page/artifact gate fails, stop at that
gate. Request/inspect the exact failure evidence and determine whether partial
artifacts are admissible for diagnosis only. Label them failed/inadmissible
validation evidence. Do not move to collection, packaging, copying, scientific
interpretation or another run until the failed invariant is understood.
Child-analysis completion does not establish owner success.

A failed farm attempt does not authorize a rerun. Investigate the first failed
invariant and earliest valid repair/debug gate, according to failure class:

- Source/state failure -> source/state repair.
- Artifact-schema/provenance failure -> artifact/provenance diagnosis.
- Renderer/page failure -> renderer/presentation diagnosis.
- Scientific closure failure -> producer/payload/consumer closure audit.

Verification precedes consumption in this acceptance order:

```text
run
-> owner/checker PASS
-> provenance/freshness PASS
-> structured payload/page-manifest PASS
-> rendered-page inspection
-> packaging/handoff
-> scientific interpretation
```

This is the acceptance ordering, not authorization to bypass a tracked owner's
implementation. File existence, `complete=true` or analysis completion alone
is insufficient. Procedure-PDF acceptance includes required page IDs,
represented children, renderer failures, provenance, required comparison
objects and actual visual legibility, as applicable.

When scalar yields and displayed spectra appear inconsistent, or a required
comparison such as SIMC is absent, establish closure before explaining physics.
Trace producer -> serializer/sidecar/checkpoint -> payload -> consumer ->
renderer and test exact producer-owned closure/integration semantics.

Once a reviewed owner/wrapper exists, retain it as authoritative. Do not
improvise an alternate manual analysis/package workflow unless a concrete
source-level blocker requires an explicit repair contract. Partial outputs may
be inspected only to diagnose the failed invariant, never promoted as an
accepted gate or delivered as accepted validation artifacts.

## Execution authority

Codex may make allowlisted local changes and deterministic checks, but must not commit, push, update remote refs, or initiate Jefferson Lab farm execution.
ChatGPT audits the actual diff; the user commits/pushes accepted changes and runs farm validation.
Workflow: Codex local changes -> ChatGPT audit -> user commit/push
-> ChatGPT pushed-state synchronization review -> user farm run when required
-> ChatGPT evidence review. Commit/push and pushed-state review are mandatory
synchronization stages, not independent scientific phases. A matching pushed
candidate with materially accurate CURRENT proceeds directly to the substantive
gate; a push alone creates no reconciliation phase. See [CODEX.md](CODEX.md)
for the canonical source-changing workflow.
