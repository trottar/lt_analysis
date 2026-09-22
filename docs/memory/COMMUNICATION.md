# KaonLT farm communication

This record owns how a requested farm operation is communicated and what is
returned for independent review. See [TOOLS.md](TOOLS.md) for canonical shell
and path facts, and the [farm validation-bundle procedure](decisions/farm-validation-bundle-procedure.md)
for detailed packaging.

## Request and return pattern

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

## Execution authority

Codex may make allowlisted local changes and deterministic checks, but must not commit, push, update remote refs, or initiate Jefferson Lab farm execution.
ChatGPT audits the actual diff; the user commits/pushes accepted changes and runs farm validation.
Workflow: Codex local changes -> ChatGPT audit -> user commit/push -> user farm
run -> ChatGPT evidence review. See [CODEX.md](CODEX.md) for the canonical
source-changing workflow.
