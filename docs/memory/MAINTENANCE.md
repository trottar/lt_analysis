# KaonLT memory maintenance

## Record roles

- **Active state:** `CURRENT.md` is the concise sole ordinary active-state
  authority and owns the exact next action.
- **Durable knowledge:** `MEMORY.md` preserves reusable cross-phase rules,
  identities, and decisions after evidence is secure.
- **Exceptional transfer:** `handoffs/CURRENT_HANDOFF.md` records exceptional
  transfer state only and cannot override CURRENT.
- **Chronology:** historical/chat indexes help a new session navigate; they do
  not outrank source or evidence.
- **Source and historical navigation:** source indexes retain provenance;
  history, chats, and import records retain historical/navigation context.
  None acquires active objective, blocker, or NEXT ownership.
- **Roadmap:** `roadmap/STATUS.md` preserves approved dependency/status
  structure and does not own the exact active next action.
- **Evidence:** `evidence/` contains farm/runtime gate records; `phases/`
  contains implementation and fix chronology; `decisions/` contains contracts.
- **User context:** `USER.md` contains stable collaboration context.
- **Operations:** `TOOLS.md` contains concise Linux/JLab command and path
  references; `COMMUNICATION.md` owns farm request/return communication.
- **Codex workflow:** `CODEX.md` owns source-changing workflow and contract
  policy.
- **Generalized lessons:** `LEARNINGS.md` contains reusable lessons, not phase
  chronology or current status.

Maintain one active objective and one authoritative next action. Do not create
competing live-state summaries or copy another record's detailed procedure.
At final-pre-push reconciliation, CURRENT's sole ordinary NEXT must remain
accurate before and immediately after the user's commit/push: name the actual
substantive next gate conditional on commit/push and pushed-state review, not
the push alone. Pushed-state review must check for consumed pre-push wording;
stale local/unpushed claims or a completed push as NEXT make memory continuity
`BLOCKED` until repaired. Do not start another scientific or farm gate from
stale CURRENT.

## Startup contract health

At substantial-work start, establish actual branch, HEAD, and worktree state.
Then read these five files in full, in this exact order:

1. `AGENTS.md`
2. `CURRENT.md`
3. `MEMORY.md`
4. `handoffs/CURRENT_HANDOFF.md`
5. `USER.md`

Only after the five-file core is read may task-directed expansion use CURRENT
direct references, the exact active task, and required canonical
evidence/decision/phase records. Do not eagerly load the whole hierarchy. Keep
`CURRENT.md` as the concise sole ordinary active-state authority under schema
3. The handoff is exceptional-only, and the roadmap contains no exact active
next action.

## Size and semantic triggers

| Record | Soft warning | Hard failure |
| --- | ---: | ---: |
| `docs/memory/CURRENT.md` | 8 KiB | 16 KiB |
| `docs/memory/handoffs/CURRENT_HANDOFF.md` | 6 KiB | 10 KiB |
| `docs/memory/MEMORY.md` | 30 KiB | 50 KiB |

`CURRENT.md` should remain very compact because it is the sole ordinary active
state. The handoff should remain exceptional and small. `MEMORY.md` may grow
more because it carries durable cross-phase knowledge. A soft-limit warning
calls for focused consolidation; a hard-limit failure requires maintenance
before further expansion.
An unresolved warning blocks beginning another phase, fix, or farm gate unless
the user approves a recorded maintenance exception. A warning is not itself a
scientific or runtime failure.

Perform maintenance after a meaningful implementation or review, accepted farm
evidence, regression, phase decision, changed blocker, or changed next action.
Do not recursively summarize summaries: recover the supporting source or
evidence first, retain its identity and path, then write the shortest accurate
statement needed for future work.

## Required maintenance sequence

1. Preserve or link authoritative evidence, source identity, and relevant diff
   before compacting prose.
2. Update only the record whose role changed; do not restate scientific status
   in unrelated chronology or chat records.
3. Keep accepted runtime, source review, active work, deferred work, and next
   actions explicitly distinct.
4. Discover `<PYTHON>` as described in [TOOLS.md](TOOLS.md). After memory
   changes regenerate and check the manifest, then run the task-final gate:
   `<PYTHON> -B tools/check_memory_health.py --root . --fail-on-warning`.
5. Inspect the resulting diff and report the exact remaining next action.

After every substantial implementation, source review reconciliation, closure,
or pushed-state handoff, Codex reports its deterministic local result and
ChatGPT independently verifies the available evidence. Both surface:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
manifest check: PASS | FAIL
```

Resolve warnings before handing off to a new gate, unless the user approves a
recorded maintenance exception. The non-strict health command remains useful
for diagnosis but cannot complete a substantial task.

The generated manifest is integrity metadata only: it contains no active-state,
date, or stored-HEAD authority. Regenerate/check it after versionable memory
changes. Bootstrap obtains repository identity dynamically; neither manifest
nor bootstrap can override CURRENT.md.
