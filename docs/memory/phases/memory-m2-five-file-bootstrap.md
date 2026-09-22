# Memory M2 five-file bootstrap ownership

## Status

`SOURCE REVIEWED` — M2 establishes the repository-owned five-file startup
contract and stable role records. It is local memory organization only, not
farm/runtime validation.

## Starting identity and scope

- Starting branch: `test`.
- Starting HEAD: `bc02e09728f2decb624e0101b5c53132190981b3`.
- Starting worktree: clean.

M2 changes only the tracked memory control layer. It preserves CURRENT as the
ordinary active-state authority while leaving its live schema-2 representation,
the handoff, roadmap, MEMORY, tools, tests, science, and farm boundary intact.

## Changed files

- `docs/memory/AGENTS.md`
- `docs/memory/README.md`
- `docs/memory/MAINTENANCE.md`
- `docs/memory/USER.md`
- `docs/memory/TOOLS.md`
- `docs/memory/LEARNINGS.md`
- `docs/memory/phases/memory-m2-five-file-bootstrap.md`
- regenerated `docs/memory/manifest.json`

## Startup ownership

The full fresh-session core is ordered: `AGENTS.md`, `CURRENT.md`,
`MEMORY.md`, `handoffs/CURRENT_HANDOFF.md`, then `USER.md`. After that full
read, task-directed expansion is limited to CURRENT references, the exact task,
and required canonical evidence/decision/phase records.

`USER.md` owns stable collaboration and operating context. `TOOLS.md` owns
concise Linux/JLab-oriented repository and memory-control references, pointing
to detailed farm-procedure and communication owners. `LEARNINGS.md` owns
generalized reusable lessons rather than active state or chronology.

## Transitional boundary and validation

The live repository remains schema 2. CURRENT, MEMORY, CURRENT_HANDOFF, and
roadmap/CURRENT were not migrated; their schema-2 mirrored metadata remains
temporary compatibility until M3 rather than independent semantic authority.

Using the discovered local Python interpreter with `-B`, the M1 33-test suite,
all three embedded self-tests, live memory health, manifest check, and bootstrap
summary pass. The direct M1 `check_schema3_bootstrap()` check against the live
M2 tree returns an empty error list. No science, farm/runtime behavior, or
status changed.

## NEXT

NEXT — Memory M3: activate schema 3, make CURRENT.md the sole active-state
representation, and normalize CURRENT_HANDOFF.md to exceptional transfer state
only.
