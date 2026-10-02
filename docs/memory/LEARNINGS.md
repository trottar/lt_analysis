# KaonLT durable learnings

## Evidence and validation

- Source review and local unit tests do not establish farm/runtime behavior.
- A bundle reporting completion is insufficient by itself: inspect provenance,
  checker gates, structured payloads or logs, and rendered pages when relevant.
- Renderer, bundle, and profile provenance remain distinct from scientific
  source provenance.
- Farm acceptance requires direct fresh evidence.

## Scientific ownership

- Diagnostic and presentation quantities must never silently become production
  corrections.
- Preserve independent owners for random subtraction, slow-proton subtraction,
  pion treatment, Method A, Method B, SIMC, yields, and cross sections.
- Method B is diagnostic/cross-check only.
- Method A promotion requires an explicit, validated F.6 decision.

## Diagnostics and checkers

- When a checker disagrees with raw evidence or implementation behavior,
  inspect the raw evidence and implementation rather than forcing data to fit
  the checker.
- When applicable, trace persisted diagnostics through producer -> serializer
  or checkpoint -> checkpoint-first payload -> consumer -> renderer.
- Retain raw quantities needed to challenge classifications.

## Scientific linkage and visual comparison

- Sequential procedure pages can be individually correct yet misleading when
  they silently cross correction-lineage boundaries.
- Every displayed scalar yield attached to a histogram must be auditable
  against that exact histogram's producer-owned integral.
- Signed-integral percentages need positive-bin, negative-bin and absolute
  support diagnostics before cancellation can explain a large relative shift.
- Overlapping comparison curves must remain independently visible. Repair
  visibility after numerical identity closes; visualization must never change
  physics, histogram contents or normalization.

- Use redundant line/marker styles and baseline-last draw order for nearly
  overlapping comparisons; never change the values or ranges for visibility.
- Present stored signed-support/aggregate diagnostics directly. A renderer
  ownership test can use deliberately distinct stored values and prohibit
  histogram reads, without accepting that fixture as a scientific identity.
- Fake ROOT draw-order/style checks do not establish real PDF legibility.

## Workflow and failure handling

- Prefer one narrow gate -> one targeted run -> fresh evidence -> inspect ->
  one coherent repair.
- Preserve negative, failed, sparse, unavailable, or rejected outcomes rather
  than deleting them from history.
- Do not perform unrelated cleanup merely to make a gate pass.
- Do not make the user the first validator of deterministic logic that can be
  checked locally.
- A reviewed profile and collector do not prove the requested farm operation
  has a tracked execution owner. Audit the complete executable path before
  handoff; if multi-step orchestration is missing, stop for a source-changing
  contract instead of composing a farm-shell sequence.
- Hard integrity failures and material active-state/provenance ambiguities
  block progress. Record nonblocking memory warnings and batch consolidation
  at checkpoints/milestone audits. Report exact byte counts, warnings, and
  manifest state after substantial gates.
- Test earlier: process and speculative hardening cannot substitute for fresh
  narrow runtime evidence.
- Memory is a checkpoint mechanism bracketing substantive work, not a parallel
  deliverable stream. Batch nonblocking drift rather than recursively creating
  repair/reconciliation phases.
- Actual-diff review, user commit/push, and pushed-state synchronization are
  necessary. Keep synchronization lightweight: verify the same source and
  materially accurate CURRENT/NEXT, then continue directly to the substantive
  gate when no concrete blocker is exposed.
- Every scientific loop must produce visible evidence: a plot, yield table,
  validated numerical comparison, accepted runtime artifact, or one directly
  evidenced scientific/runtime blocker with one coherent repair. Documentation
  alone is not completion.

- A new procedure-PDF analysis figure is a source/presentation task: trace
  authoritative inputs, implement it in tracked source, test locally, then
  farm-render it. Ad-hoc chat plots do not fulfill that deliverable.

- Persist owner orchestration failures separately from a child analysis log:
  later verification/collection errors cannot appear in that child stream.
- Run deterministic detached collector/source checks before expensive analysis;
  keep final source rechecks, fresh artifact checks and ZIP verification intact.
- Publish per-attempt status atomically and refuse existing attempt files. Keep
  owner diagnostics out of a frozen scientific bundle inventory.

## Gate admissibility and chat health

- A failed downstream gate stops forward progression even if its analysis child
  completed. Child-process success and owner success are distinct.
- Artifact existence is not artifact admissibility. A generated PDF does not
  justify packaging or handoff after a failed manifest/page gate; partial
  outputs are failed-gate diagnostic evidence only.
- Never rerun by reflex. Investigate the first failed invariant and earliest
  valid repair/debug gate before another operation.
- Suspicious scalar/histogram disagreement requires exact producer-owned
  closure before physics interpretation.
- Full in-chat health checks, periodic health pulses and explicit evidence
  labels expose stale state and unverified assumptions before silent drift.

## Repository memory and provenance

- Current source and diff outrank stale summaries for implementation.
- Fresh runtime evidence outranks source review for runtime acceptance.
- A stored repository HEAD is a timestamped observation, not permanent live
  identity.
- Open the canonical record rather than recursively summarizing old summaries.
