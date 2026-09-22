# KaonLT memory schema-v3 semantic contract

## Status and transition boundary

Schema 2 is the live legacy representation throughout Memory M1.  Schema 3 is
the target representation, validated only with deterministic fixtures and
tests in M1.  This record defines target semantics; it does not migrate any
live memory record or upgrade scientific/runtime status.

`docs/memory/` remains repository-owned durable continuity.  Native assistant
memory is supplemental only.  Actual branch, HEAD, and worktree state are
established dynamically at substantial-work startup; a stored SHA retains an
explicit provenance role and is never a permanently live HEAD claim.

## Schema-2 legacy representation

During M1, schema 2 preserves its exact current behavior:

- the full active-state frontmatter field set remains required on CURRENT;
- CURRENT, handoff, and roadmap active-state dictionaries remain identical;
- commit syntax and resolution, CURRENT heading order, local links,
  runtime-closure evidence links, authority markers, top-level role uniqueness,
  size thresholds, and manifest freshness remain enforced;
- the generated manifest remains schema version 2 and retains `active_state`,
  `files`, `generated_date_utc`, and `observed_git_head`.

M1 does not rewrite history to imply that schema 3 is already active.

## Schema-3 CURRENT

Target `CURRENT.md` begins with exactly this minimal frontmatter:

```
---
memory_schema: 3
---
```

No other schema-3 frontmatter field is permitted.  Active objective, current
work item, active status, exact next action, and baseline/source/bundle-profile
identities belong to rendered CURRENT content with explicit provenance roles,
not duplicate frontmatter.

CURRENT is the sole ordinary active-state authority.  It owns the active
objective, work item, verified state, active source/evidence identities,
blockers, exact next action, success criteria, do-not-reopen boundaries, and
direct canonical references.  No ordinary memory record may become an
alternative current-state authority.

Its H2 headings occur exactly once and in this order:

1. `Active Objective`
2. `Current Work Item`
3. `Verified State`
4. `Source / Evidence Identity`
5. `Blockers`
6. `Next Action`
7. `Success Criteria`
8. `Do Not Reopen Without New Evidence`
9. `Relevant References`

Every section is non-empty.  CURRENT contains exactly one active `NEXT —`
statement, and it occurs in `## Next Action`; it describes one coherent
immediate action.  Every `CLOSED / RUNTIME VALIDATED` claim in CURRENT remains
locally tied to a canonical `evidence/` link.  Source review never satisfies a
runtime-closure claim.

Stored SHAs are allowed only with unambiguous roles such as reviewed analysis
source, evaluated farm source, bundle/profile source, renderer source,
predecessor/baseline, or artifact provenance.  They must not be described as a
permanently live current HEAD.

## Schema-3 exceptional handoff

Target `handoffs/CURRENT_HANDOFF.md` is not a second CURRENT.  It has exactly:

```
# Current KaonLT handoff

## Transfer State

## Resume
```

Transfer State is never empty.  A normal stable state uses exactly `No
exceptional transfer state is recorded.`  The handoff explicitly states that
`CURRENT.md` is the sole authoritative resumable state and that the handoff
cannot override `CURRENT.md`.  Resume directs the next session to CURRENT and then
only task-relevant canonical references.  Schema-3 handoff has no mirrored
active-state frontmatter or routine CURRENT headings.

## Five-file bootstrap target

After M2, read this core completely and in exact order before selective
expansion:

1. `AGENTS.md`
2. `CURRENT.md`
3. `MEMORY.md`
4. `handoffs/CURRENT_HANDOFF.md`
5. `USER.md`

Only then expand from CURRENT direct references, the exact active task, and
required canonical evidence/decision/phase records.  Do not eagerly load the
whole memory hierarchy.  M1 fixtures may validate this target; M1 does not
modify live startup records to enact it.

M2 will create these additional records with distinct roles:

- `USER.md`: stable collaboration preferences and project operating context,
  never scientific evidence or active state.
- `TOOLS.md`: canonical Linux/JLab/repository commands and conventions, never
  a Windows packaging layer.
- `LEARNINGS.md`: reusable generalized lessons, never current status or phase
  chronology.

## M1 compatibility implementation

Schema detection is deterministic from valid CURRENT frontmatter only.  Missing,
malformed, duplicate, ambiguous, or unsupported markers fail; filenames, USER
presence, branch, and Git HEAD never select a schema.  Transitional manifest
parsing returns the unchanged full schema-2 state for schema 2 and exactly
`{"memory_schema": 3}` for minimal schema 3.  The manifest itself remains
schema version 2; M1 does not redesign its fields or bootstrap behavior.

## Deferred strictness

M1 deliberately defers full MEMORY chronology normalization to M4; Linux/JLab
role-content normalization to M5; roadmap/source/history normalization to M6;
role-specific thresholds, final duplicate-role prevention, rendered-text
checks, and strict MEMORY active-heading prohibition to M7; manifest/bootstrap
semantic redesign to M8; and fresh-session proof to M9.  Reusable validation
primitives may support M1 fixtures, but no M2-M9 live migration occurs here.

## Preserved scientific and actor boundaries

This memory contract preserves separate ownership of random subtraction,
slow-proton subtraction, pion-background treatment, HGCer Method A, HGCer
Method B, SIMC, yields, cross sections, diagnostics, validation, and
presentation.  Method B remains diagnostic/cross-check only; Method A remains
detached pending an explicitly validated F.6 promotion.  Farm/ROOT/PyROOT/full
runtime validation remains farm-only.  Codex performs local changes and checks,
ChatGPT audits actual diffs, and the user commits/pushes accepted work and runs
the farm.
