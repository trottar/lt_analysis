# KaonLT development memory

This directory is the active, repository-owned engineering notebook for KaonLT
development. It is maintained with the code, not reconstructed from chat or
commit titles.

For substantial work, establish the live branch, HEAD, and worktree first;
then read [current state](CURRENT.md), the records it directly links, relevant
durable guidance in [memory](MEMORY.md), and the
[current handoff](handoffs/CURRENT_HANDOFF.md). Do not load the whole tree
without a task-specific reason. The tracked operating rules are in
[AGENTS.md](AGENTS.md), with [maintenance](MAINTENANCE.md),
[communication](COMMUNICATION.md), and [Codex workflow](CODEX.md) guidance.

- [CURRENT.md](CURRENT.md) — concise, authoritative live state.
- [MEMORY.md](MEMORY.md) — curated durable knowledge and recurring rules.
- [handoffs/CURRENT_HANDOFF.md](handoffs/CURRENT_HANDOFF.md) — continuation
  state for a fresh session.
- `evidence/` — farm/runtime validation records only.
- `decisions/` — significant architectural or scientific decisions.
- `investigations/` — active, closed, and deferred investigations.
- `phases/` — phase/fix implementation and validation summaries.
- `roadmap/CURRENT.md` — approved phase structure only.

## Memory controls

- [AGENTS.md](AGENTS.md) — repository-tracked operating rules.
- [MAINTENANCE.md](MAINTENANCE.md) — maintenance and integrity procedure.
- [COMMUNICATION.md](COMMUNICATION.md) — concise farm-delivery rules.
- [CODEX.md](CODEX.md) — Codex task and evidence boundaries.
- [chats/CHAT_INDEX.md](chats/CHAT_INDEX.md) — distilled chat-era navigation.
- [templates/CODEX_CONTRACT.md](templates/CODEX_CONTRACT.md),
  [templates/MEMORY_UPDATE.md](templates/MEMORY_UPDATE.md), and
  [templates/RUNTIME_EVIDENCE.md](templates/RUNTIME_EVIDENCE.md) — concise
  records for new work.
- [manifest.json](manifest.json) — generated integrity index only; regenerate
  it with `python tools/update_memory_manifest.py --write`.

Use `python tools/check_memory_health.py` to check this control layer and
`python tools/memory_bootstrap.py` for a narrow startup summary.

At bootstrap, records deliberately distinguish observed repository facts from
facts that need an authoritative project handoff. This prevents historical
commit messages from becoming accidental validation claims.
