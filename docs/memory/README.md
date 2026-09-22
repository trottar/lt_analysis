# KaonLT development memory

This directory is KaonLT's repository-owned development memory. It is
maintained with the code, not reconstructed from chat or commit titles.

## Startup contract

At substantial-work start, establish actual branch, HEAD, and worktree state,
then read these five files in full and in this exact order:

1. [AGENTS.md](AGENTS.md)
2. [CURRENT.md](CURRENT.md)
3. [MEMORY.md](MEMORY.md)
4. [handoffs/CURRENT_HANDOFF.md](handoffs/CURRENT_HANDOFF.md)
5. [USER.md](USER.md)

Only then expand selectively from CURRENT direct references, the exact active
task, and required canonical evidence/decision/phase records. Do not load the
whole hierarchy without task-specific need. `CURRENT.md` is the sole ordinary
active-state authority under schema 3; the handoff is exceptional transfer
state only, and the roadmap preserves dependency/status structure without an
exact active next action.

## Navigation

- [AGENTS.md](AGENTS.md) — high-level behavior, evidence, science, and
  execution boundaries.
- [CURRENT.md](CURRENT.md) — concise authoritative active state.
- [MEMORY.md](MEMORY.md) — curated durable cross-phase knowledge.
- [handoffs/CURRENT_HANDOFF.md](handoffs/CURRENT_HANDOFF.md) — exceptional
  transfer state only.
- [USER.md](USER.md) — stable collaboration and delivery preferences.
- [TOOLS.md](TOOLS.md) — Linux/POSIX/JLab commands and path facts.
- [COMMUNICATION.md](COMMUNICATION.md) — farm request/return interaction.
- [CODEX.md](CODEX.md) — source-changing workflow and contract policy.
- [MAINTENANCE.md](MAINTENANCE.md) — memory maintenance procedure.
- [LEARNINGS.md](LEARNINGS.md) — generalized reusable lessons.
- [chats/CHAT_INDEX.md](chats/CHAT_INDEX.md) — historical chat navigation.
- [templates/CODEX_CONTRACT.md](templates/CODEX_CONTRACT.md),
  [templates/MEMORY_UPDATE.md](templates/MEMORY_UPDATE.md), and
  [templates/RUNTIME_EVIDENCE.md](templates/RUNTIME_EVIDENCE.md) — concise
  templates for new records.
- `evidence/` — farm/runtime validation records.
- `decisions/` — scientific and architectural contracts.
- `investigations/` — active, closed, and deferred investigations.
- `phases/` — phase/fix implementation and validation chronology.
- `roadmap/CURRENT.md` — approved dependency/status structure only.

## Memory controls

Use [MAINTENANCE.md](MAINTENANCE.md) for maintenance, [TOOLS.md](TOOLS.md) for
commands, [COMMUNICATION.md](COMMUNICATION.md) for farm delivery, and
[CODEX.md](CODEX.md) for source-changing work. The generated
[manifest.json](manifest.json) is an integrity index only.

Bootstrap distinguishes dynamic repository observations from canonical
repository-memory/evidence facts; the exceptional handoff cannot override
CURRENT.
