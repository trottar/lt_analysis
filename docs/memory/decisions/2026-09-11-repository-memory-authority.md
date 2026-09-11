# Repository-owned KaonLT memory authority

Decision: `docs/memory/` is the authoritative development-state record for
KaonLT. Codex native `~/.codex/memories` is supplemental associative context
only.

Rationale: project state must travel with the repository, be reviewable beside
the source, and be reconciled with actual Git and farm evidence rather than
depending on a particular Codex session.

Consequences:

- Substantial work begins by reading the current state, durable memory, and
  handoff, then comparing the recorded accepted HEAD with the live `test` HEAD.
- Meaningful implementation, review, investigation, decision, and validation
  outcomes update the relevant memory record before work is considered complete.
- Runtime acceptance is entered only from supplied farm/runtime artifacts; a
  commit subject or local source test is insufficient.

This is a development-process decision. It makes no claim about analysis
runtime validation or the status of any historical physics phase.
