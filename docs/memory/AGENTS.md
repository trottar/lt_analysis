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
- For substantial work, establish live branch, HEAD, and worktree first, then
  follow the five-file startup contract below before task-directed expansion.
- Use only these work-state labels: `CLOSED / RUNTIME VALIDATED`, `SOURCE
  REVIEWED`, `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`, `ACTIVE`,
  `DEFERRED`, `BLOCKED`, and `NEXT`.

## Five-file startup contract

At substantial-work start, first establish actual branch, HEAD, and worktree
state. Read these five files in full, in this exact order:

1. `AGENTS.md`
2. `CURRENT.md`
3. `MEMORY.md`
4. `handoffs/CURRENT_HANDOFF.md`
5. `USER.md`

Only after those five are read in full may a session expand selectively from:

- CURRENT direct references
- exact active task
- required canonical evidence/decision/phase records

Do not eagerly load the whole memory hierarchy. Identify the active objective,
evidence boundary, blockers, and NEXT from the core; inspect current source or
diff whenever implementation is involved; and distinguish source proof from
farm/runtime proof. Schema 3 is active: `CURRENT.md` is the sole ordinary
active-state authority, the handoff is exceptional transfer state only, and the
roadmap preserves dependency/status structure without owning the exact next
action.

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

## Execution authority

Codex may make allowlisted local changes and run deterministic local checks.
Codex must not commit, push, update remote refs, initiate Jefferson Lab farm
execution, or claim farm/runtime validation. ChatGPT independently audits the
actual diff. The user alone commits/pushes accepted changes and runs farm
validation. Workflow: Codex local changes -> ChatGPT audit -> user commit/push
-> user farm run -> ChatGPT evidence review.
