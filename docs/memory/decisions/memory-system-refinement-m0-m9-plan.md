# KaonLT memory-system refinement: M0-M9 architecture

## Decision and scope

`docs/memory/` remains KaonLT's repository-owned durable continuity layer.
Native assistant memory is supplemental associative context only; it never
overrides current source, repository memory, or supplied farm evidence.

This plan adapts the successful repository-memory ideas used by the
PrivyHub/fantasy-football pattern to KaonLT's needs.  It does not copy that
project's Windows, PowerShell, or `.ffpkg` infrastructure.  KaonLT's target
operational environment is Linux/JLab, and its scientific ownership and
farm-validation boundaries remain first-class requirements.

This is a staged representation refinement, not a scientific or analysis
change.  Each phase requires a narrowly reviewed implementation contract and
local validation appropriate to that phase.  Codex makes local changes and
checks; ChatGPT audits the actual diff; the user commits/pushes accepted work
and performs Jefferson Lab farm runs.  Project lifecycle-hook dispatch remains
`BLOCKED / DEFERRED` and is not a dependency of this plan.

The M0 baseline observations are recorded separately in
[the M0 audit](../evidence/memory-m0-baseline-audit.md).

## Target operating model

The eventual fresh-session bootstrap is deliberately small and ordered.  Read
these five files completely before expanding selectively from `CURRENT.md` or
the task's direct references:

1. `AGENTS.md`
2. `CURRENT.md`
3. `MEMORY.md`
4. `handoffs/CURRENT_HANDOFF.md`
5. `USER.md`

`CURRENT.md` will be the sole ordinary active-state authority.  A handoff
will carry only exceptional transfer state, not a second live checkpoint.
`MEMORY.md` will retain stable reusable knowledge, while phase chronology,
farm evidence, decisions, source indexes, and history keep their distinct
roles.  Actual branch, HEAD, and worktree status are dynamic startup facts;
stored SHAs retain explicit historical, reviewed-source, evaluated-source, or
artifact-provenance roles instead of posing as permanent live identity.

## Phased migration

### M0 — baseline audit and migration map

- Establish the exact local starting branch, HEAD, and worktree state.
- Execute the existing control tools against that untouched state.
- Inventory existing memory roles, duplicated ownership, and actual tool
  semantics.
- Record the migration map without rewriting active representation.
- Preserve the established scientific/runtime frontier exactly.

M0 is an audit and documentation checkpoint only.  It neither implements M1
nor redesigns active memory files.

### M1 — schema-v3 contract and compatibility checker

- Define the target representation and its compatibility boundary with the
  legacy representation.
- Teach health tooling the target semantics while retaining temporary legacy
  compatibility.
- Add deterministic health tests.
- Do not broadly migrate memory records yet.

### M2 — five-file bootstrap ownership

- Establish the five-file fresh-session core in the order stated above.
- Add and separate `USER.md`, `TOOLS.md`, and `LEARNINGS.md` with explicit
  canonical roles.
- Require full reading of the five core files before task-directed expansion.

### M3 — single active state and handoff normalization

- Make `CURRENT.md` the sole ordinary active-state authority.
- Restrict the handoff to exceptional transfer state.
- Remove mirrored ordinary active-state ownership.
- Eliminate self-referential current-checkpoint SHA semantics.

### M4 — durable MEMORY normalization

- Keep stable cross-phase knowledge in `MEMORY.md`.
- Retain detailed phase chronology in `phases/` and accepted runtime details
  in `evidence/`.
- Promote reusable failures and lessons to their appropriate durable roles.
- Preserve every unique scientific-provenance identity while normalizing.

### M5 — Linux/JLab operational-role separation

- Separate canonical ownership among `AGENTS`, `USER`, `TOOLS`,
  `COMMUNICATION`, `CODEX`, and `MAINTENANCE`.
- Keep the operational model Linux/JLab-oriented.
- Do not introduce Windows `.ffpkg` or PowerShell infrastructure.

### M6 — roadmap, source, and history normalization

- Limit the roadmap to phase/dependency/status structure rather than today's
  exact next action.
- Stop source indexes from describing a stored SHA as permanently live.
- Treat chat, import, and historical material as navigation/history rather
  than startup control.

### M7 — strict semantic memory-health enforcement

Enforce the five-file bootstrap contract, `CURRENT.md` sole authority, one
active objective, one exact next action, exceptional-only handoffs, MEMORY
role separation, canonical headings/order, local-link validity,
duplicate-role prevention, representation/rendering checks, and role-specific
size thresholds.

### M8 — manifest and bootstrap identity semantics

- Keep the manifest as integrity metadata only.
- Do not copy semantic active state into the manifest.
- Do not retain a stored field that purports to be permanently current HEAD.
- Determine branch, HEAD, and worktree dynamically at startup.
- Have bootstrap report repository facts without replacing reading the five
  core files.

### M9 — fresh-session integration audit

Prove that a genuinely fresh session can recover the repository/branch and
dynamic-HEAD requirement, active objective, exact next action, relevant
scientific evidence, Method-A/Method-B boundaries, production/diagnostic/
presentation separation, source-review versus farm-runtime distinction,
Codex/user/farm actor sequence, Linux/JLab environment, and an exceptional
handoff when one exists.  Only after this passes is the large
chat-continuation workflow superseded.

## Invariants throughout M0-M9

### Repository authority and fresh identity

Repository memory remains authoritative continuity.  For every substantial
task, establish actual branch, HEAD, and worktree state.  A stored SHA is a
historical/scientific/provenance identity with a stated role, never an eternal
claim about the current checkout.

### Evidence precedence

For implementation, precedence is current source/diff, applicable newer
source evidence, canonical durable records, then historical summaries.  For
runtime acceptance, precedence is fresh farm/runtime evidence, accepted
canonical runtime evidence, source review, then historical summaries.  Source
review never implies farm validation.

### Scientific ownership and Method boundaries

Keep random subtraction, slow-proton subtraction, pion-background treatment,
HGCer Method A, HGCer Method B, SIMC comparisons, yields, final cross sections,
diagnostics, validation/checkers, and presentation separate.  Diagnostics and
presentation never silently become production corrections.  Method B remains
diagnostic/cross-check only.  Method A remains detached until an explicitly
validated F.6 production-promotion decision.

### Farm and human actor boundaries

Jefferson Lab farm, ROOT/PyROOT, and full-runtime validation are farm-only
unless fresh artifacts are supplied.  Codex does not commit, push, or run the
farm; ChatGPT audits actual changes; the user commits/pushes accepted changes
and performs farm runs.

## Preserved frontier during this migration

The memory refinement does not reinterpret or alter the current scientific
frontier: F.6.2 and F.6.2.Fix.5 remain `CLOSED / RUNTIME VALIDATED`; E.8.1
remains `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; F.6.3 remains
`BLOCKED` pending E.8.1; and F.6.4 remains `BLOCKED` pending F.6.3 evidence.
The canonical-five-setting `Q4p4W2p74` E.8.1 farm gate remains the next
scientific continuation, with `Left / highe` inspected first.
