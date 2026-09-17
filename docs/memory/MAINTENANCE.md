# KaonLT memory maintenance

## Record roles

- **Active state:** `CURRENT.md` is the concise current state and authoritative
  immediate next action.
- **Durable knowledge:** `MEMORY.md` preserves reusable rules, identities, and
  decisions after evidence is secure.
- **Chronology and continuation:** handoffs and historical/chat indexes help a
  new session navigate; they do not outrank source or evidence.
- **Evidence:** `evidence/` contains farm/runtime gate records and is preserved
  before any related summary is shortened.

Maintain one active objective and one authoritative next action. Do not create
competing live-state summaries.

## Size and semantic triggers

The active files `CURRENT.md`, `MEMORY.md`, and
`handoffs/CURRENT_HANDOFF.md` have a 16 KiB soft target and 32 KiB hard limit.
A soft-limit warning calls for focused consolidation; a hard-limit failure
requires maintenance before further expansion.

Perform maintenance after a meaningful implementation or review, accepted farm
evidence, regression, phase decision, changed blocker, or changed next action.
Do not summarize a summary recursively: recover the supporting source or
evidence first, retain its identity and path, then write the shortest accurate
statement needed for future work.

## Required procedure

1. Preserve or link the authoritative evidence, source identity, and relevant
   diff before compacting prose.
2. Update only the record whose role changed; do not restate scientific status
   in unrelated chronology or chat records.
3. Keep accepted runtime, source review, active work, deferred work, and next
   actions explicitly distinct.
4. Run `python tools/update_memory_manifest.py --write`,
   `python tools/update_memory_manifest.py --check`, and
   `python tools/check_memory_health.py` after structural memory changes.
5. Inspect the resulting diff and report the exact remaining next action.

The generated manifest is an integrity index, never a scientific, source, or
runtime authority.
