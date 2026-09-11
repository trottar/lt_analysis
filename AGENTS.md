# KaonLT development instructions

`docs/memory/` is active project state, not passive documentation. It is the
authoritative development memory for this repository. Codex native
`~/.codex/memories` is supplemental associative memory only and must never
override current repository source, `docs/memory/`, or newer supplied runtime
evidence.

## Required memory protocol

At the beginning of every substantial KaonLT task (implementation, review,
investigation, validation, scientific decision, or roadmap work), Codex must:

1. Read `docs/memory/CURRENT.md`.
2. Read `docs/memory/MEMORY.md`.
3. Read `docs/memory/handoffs/CURRENT_HANDOFF.md`.
4. Inspect the decisions, evidence, investigations, and phase records relevant
   to the requested scope.
5. Establish the actual current `test` HEAD and compare it with the recorded
   accepted HEAD.
6. Resolve any discrepancy from current source and newer runtime evidence;
   never silently trust stale memory.

During and after every meaningful task, update the appropriate repository
memory record when the work changes project understanding. This includes an
implementation or review, meaningful source discovery, farm validation result,
regression, diagnostic conclusion, scientific or architectural decision,
phase/fix status, deferred issue, roadmap change, and the exact next step.
Memory maintenance is part of completing the task, not optional cleanup.

Use only these status labels for work-state claims:

- `CLOSED / RUNTIME VALIDATED`
- `SOURCE REVIEWED`
- `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`
- `ACTIVE`
- `DEFERRED`
- `BLOCKED`
- `NEXT`

Do not infer an accepted HEAD, phase status, or runtime result from commit
messages, a source-level test, or the existence of a validation tool. When the
record is incomplete, say so explicitly and obtain the authoritative project
handoff or fresh runtime evidence.

## Validation and ownership rules

- The real KaonLT runtime exists only on the Jefferson Lab farm. Never claim
  farm, ROOT/PyROOT, full-analysis, or runtime validation without supplied
  runtime evidence.
- Always distinguish scientific calculation, production correction/subtraction,
  runtime integration, diagnostics, checkers, and presentation-only plotting.
- Preserve validated random subtraction, proton subtraction, pion treatment,
  HGCer Method-A/Method-B independence, SIMC handling, yield calculation,
  acceptance/binning, efficiencies, L/T separation, and uncertainties unless
  the task explicitly concerns them.
- Do not reopen a closed/runtime-validated phase without concrete regression
  evidence.
- Make source changes through narrow phase/fix implementation contracts. After
  implementation, inspect the actual diff and relevant runtime path rather
  than treating an implementation summary as proof.
- The farm workflow is: one narrow gate, targeted farm run, fresh artifacts,
  inspect evidence, then either `PASS` or one coherent repair. The user
  performs farm validation.
- A validation bundle contains only artifacts required for its current gate.
  Source-level `PASS` is not runtime `PASS`.

## Memory record responsibilities

- `CURRENT.md` is the concise authoritative current state and immediate next
  action.
- `MEMORY.md` is curated durable knowledge, not an append-only timeline.
- `handoffs/CURRENT_HANDOFF.md` must let a new Codex session continue without
  rediscovery.
- `evidence/` holds farm/runtime gate records with commit, setting, artifacts,
  and conclusion.
- `decisions/`, `investigations/`, and `phases/` hold durable supporting
  records. `roadmap/CURRENT.md` contains only approved phase structure.
