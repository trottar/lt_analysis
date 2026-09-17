# KaonLT memory operating rules

This tracked file governs repository-owned memory only. It does not replace
the root `AGENTS.md` or machine-local instructions.

## Authority and startup

- For implementation, resolve current exact source and relevant diff before
  older records. For runtime acceptance, direct farm evidence is required and
  cannot be inferred from source review.
- When evidence competes, use current source/diff, then newer farm/runtime
  evidence, then newer validation or handoff artifacts, then tracked durable
  memory, then older chat summaries. State any source-versus-runtime distinction.
- For substantial work, establish live branch, HEAD, and worktree first. Read
  `CURRENT.md`, its directly linked task records, relevant sections of
  `MEMORY.md`, and `handoffs/CURRENT_HANDOFF.md`; do not eagerly load the full
  memory tree.
- Use only these work-state labels: `CLOSED / RUNTIME VALIDATED`, `SOURCE
  REVIEWED`, `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`, `ACTIVE`,
  `DEFERRED`, `BLOCKED`, and `NEXT`.

## Scientific and validation boundaries

- Preserve separate ownership of random subtraction, slow-proton treatment,
  pion subtraction, HGCer Method A, HGCer Method B, SIMC, yields, cross
  sections, diagnostics, and presentation.
- Method B remains a diagnostic/cross-check only. It never becomes a production
  pion-weight adjustment through memory, diagnostics, or presentation work.
- Method A remains detached until an explicit, validated production-promotion
  decision. Local source review is not Jefferson Lab farm validation.
- Do not reopen a closed phase without concrete regression evidence. Do not
  turn a commit message, local test, bundle creation, or documentation update
  into a farm/runtime claim.
