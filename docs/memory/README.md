# KaonLT development memory

This directory is the active, repository-owned engineering notebook for KaonLT
development. It is maintained with the code, not reconstructed from chat or
commit titles. Read the current state, durable memory, and handoff before
substantial work; see `AGENTS.md` for the mandatory protocol.

- `CURRENT.md` — concise, authoritative live state.
- `MEMORY.md` — curated durable knowledge and recurring rules.
- `handoffs/CURRENT_HANDOFF.md` — continuation state for a fresh session.
- `evidence/` — farm/runtime validation records only.
- `decisions/` — significant architectural or scientific decisions.
- `investigations/` — active, closed, and deferred investigations.
- `phases/` — phase/fix implementation and validation summaries.
- `roadmap/CURRENT.md` — approved phase structure only.

At bootstrap, records deliberately distinguish observed repository facts from
facts that need an authoritative project handoff. This prevents historical
commit messages from becoming accidental validation claims.
