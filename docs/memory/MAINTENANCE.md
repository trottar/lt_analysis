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
At each substantive checkpoint, CURRENT's sole ordinary NEXT must remain
accurate before and immediately after the user's commit/push: name the actual
substantive next gate conditional on commit/push and pushed-state review, not
the push alone. Pushed-state review must check for consumed pre-push wording;
source identity, exact gate-relevant changed paths/blobs, and CURRENT/NEXT
continuity are the required synchronization audit. If source matches the
reviewed candidate and memory is materially accurate, proceed directly to the
substantive gate. A push alone requires no new source-change contract or broad
memory reconciliation. Cosmetic/historical wording drift is recorded and
batched unless it makes active meaning ambiguous.

The default loop is opening memory checkpoint -> substantive science/
implementation -> deterministic local checks -> ChatGPT actual-diff review ->
user commit/push -> ChatGPT pushed-state synchronization review -> narrow
farm/runtime gate when required -> fresh evidence review -> closing memory
checkpoint. The middle synchronization stages are required, not independent
scientific phases; they do not recursively create maintenance work unless a
concrete blocker is exposed. See the [workflow decision](decisions/memory-bracketed-scientific-throughput.md).

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
calls for focused consolidation at the next checkpoint/milestone audit;
a hard-limit failure requires maintenance before further expansion.
A warning blocks only when it creates material ambiguity about current source
identity, accepted evidence/status, frozen scientific interfaces, active
scientific ownership/blocker, or exact NEXT, or violates a hard integrity requirement.
Hard failures remain blocking. Record nonblocking warnings, continue the active
milestone, and batch their correction; a warning is not runtime evidence.

Perform maintenance primarily at the opening checkpoint for a substantive
milestone, after accepted implementation/runtime evidence, when a concrete
blocker changes active state, and at major milestone audits. Batch accumulated
nonblocking corrections at completed scientific implementations, accepted farm
evidence, real blocker/NEXT changes, promotion decisions, and major E.8/F.6
gates. Do not start standalone repair cycles for cosmetic or historical drift
without material active-state ambiguity. Maintenance must not become an
indefinitely recursive sequence between scientific gates.
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
   changes regenerate and check the manifest, then run the ordinary task-final
   gate: `<PYTHON> -B tools/check_memory_health.py --root .`.
   Hard health failures return nonzero. Require `--fail-on-warning` only for an
   explicit milestone/zero-warning audit, a memory-hardening task that owns
   warning elimination, or an active warning classified as materially blocking.
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
warning classification: <none or each warning: blocking | nonblocking>
manifest check: PASS | FAIL
```

Repair hard failures and blocking warnings before handing off to a new gate.
The ordinary non-strict command can complete a substantial task when all hard
checks pass and each nonblocking warning is recorded and explicitly scheduled
for the next checkpoint or milestone audit. A nonblocking warning alone does
not require a new reconciliation phase.

The generated manifest is integrity metadata only: it contains no active-state,
date, or stored-HEAD authority. Regenerate/check it after versionable memory
changes. Bootstrap obtains repository identity dynamically; neither manifest
nor bootstrap can override CURRENT.md.
