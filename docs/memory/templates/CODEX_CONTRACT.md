# <phase or fix> Codex contract

- **Objective and exact starting HEAD:**
- **Allowed files / frozen files:**
- **Scientific ownership:**
- **Before and after behavior:**
- **Preserved source and runtime paths:**
- **Positive, negative, and regression checks:**
- **Forbidden shortcuts:**
- **Local validation:**
- **Farm-validation boundary:**
- **Operational readiness:** State whether the task leads to a farm gate.
  Audit input authority -> producer/materializer/analyzer/renderer (if any) ->
  verification/checker (if any) -> collector/packager -> invocation owner ->
  expected returned artifact. Name tracked source for every executable step,
  with local deterministic checks where possible, independent source review,
  push, and pushed-state review. State whether multi-step orchestration needs
  a driver/wrapper. Missing tracked, reviewed orchestration blocks farm handoff;
  one reviewed CLI may own a complete direct operation.
- **Memory health:** After memory changes, regenerate/check the manifest.
  Report exact CURRENT/MEMORY/CURRENT_HANDOFF byte counts and warnings. The
  ordinary task-final command is
  `<PYTHON> -B tools/check_memory_health.py --root .`.
  Reserve `--fail-on-warning` for explicit milestone/zero-warning audits,
  memory-hardening tasks owning warning elimination, or warnings already
  classified as materially blocking. Classify each warning as blocking or
  nonblocking. Hard failures and material ambiguity about current source
  identity, accepted evidence/status, frozen scientific interfaces, active
  scientific ownership/blocker, or exact NEXT block the next gate. With all
  hard checks passing, record and schedule nonblocking warnings for the next
  checkpoint or milestone audit; they alone create no reconciliation phase.
  Codex and ChatGPT each report health after substantial implementation, review reconciliation,
  closure, or pushed-state handoff.
- **Push-stable NEXT:** The task must leave CURRENT's sole
  ordinary NEXT the substantive next gate conditional on user-controlled
  commit/push and pushed-state review, never the push itself. Pushed-state
  synchronization review checks source identity, exact gate-relevant changed
  paths/blobs, and CURRENT/NEXT continuity. A matching pushed candidate with
  materially accurate CURRENT proceeds directly to the substantive gate; a
  push alone requires no separate final-pre-push reconciliation task or phase.
  Reconcile only concrete source/provenance/active-state blockers. Record and
  batch cosmetic/historical wording drift unless it changes active meaning.
- **Execution authority:** Codex performs local changes/checks only; ChatGPT
  audits the actual diff; the user alone commits/pushes and runs the farm.
  Workflow: Codex local changes -> ChatGPT audit -> user commit/push
  -> ChatGPT pushed-state synchronization review -> user farm run when required
  -> ChatGPT evidence review. Commit/push and pushed-state review are mandatory
  synchronization stages, not independent scientific phases.
- **Diff audit and acceptance criteria:**
- **Hard stop:**
