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

## Workflow and failure handling

- Prefer one narrow gate -> one targeted run -> fresh evidence -> inspect ->
  one coherent repair.
- Preserve negative, failed, sparse, unavailable, or rejected outcomes rather
  than deleting them from history.
- Do not perform unrelated cleanup merely to make a gate pass.
- Do not make the user the first validator of deterministic logic that can be
  checked locally.

## Repository memory and provenance

- Current source and diff outrank stale summaries for implementation.
- Fresh runtime evidence outranks source review for runtime acceptance.
- A stored repository HEAD is a timestamped observation, not permanent live
  identity.
- Open the canonical record rather than recursively summarizing old summaries.
