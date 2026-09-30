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
  task-final command is
  `<PYTHON> -B tools/check_memory_health.py --root . --fail-on-warning`.
  Do not hand off to another gate with unresolved warnings
  unless the user approves a recorded maintenance exception. Codex and ChatGPT
  each report health after substantial implementation, review reconciliation,
  closure, or pushed-state handoff.
- **Execution authority:** Codex performs local changes/checks only; ChatGPT
  audits the actual diff; the user alone commits/pushes and runs the farm.
  Workflow: Codex local changes -> ChatGPT audit -> user commit/push -> user
  farm run -> ChatGPT evidence review.
- **Diff audit and acceptance criteria:**
- **Hard stop:**
