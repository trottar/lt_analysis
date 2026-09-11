# Durable KaonLT project knowledge

This is curated project knowledge, not a chronological log. It records source,
history, and runtime evidence separately so that future work cannot turn an
implementation or a checker into an unsupported physics claim.

## Authority and status semantics

- `docs/memory/` is the authoritative development-state record. Native Codex
  memory is supplemental context only and never overrides current source,
  repository memory, or newer farm evidence.
- Begin substantial work by reading `CURRENT.md`, this file, and
  `handoffs/CURRENT_HANDOFF.md`, then compare their recorded source identity
  with live `test`.
- Use only these work-state labels: `CLOSED / RUNTIME VALIDATED`, `SOURCE
  REVIEWED`, `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`, `ACTIVE`,
  `DEFERRED`, `BLOCKED`, and `NEXT`.
- Current source decides implementation state. Fresh farm artifacts decide
  runtime state. Commit messages, a clean tree, tests, a profile, or a bundle
  with `complete=true` do not by themselves decide either acceptance or
  runtime validity.

## Scientific ownership boundaries

- Preserve the separation of random subtraction, slow-proton contamination
  treatment, pion-background treatment, HGCer diagnostics, SIMC comparison,
  yield extraction, and final cross-section analysis. The slow-proton and pion
  paths address distinct physical backgrounds.
- Preserve existing cuts, normalizations, templates, priors, fit windows,
  component definitions, binning, production corrections, efficiencies,
  acceptance, L/T separation, and uncertainty propagation unless an explicit
  narrow contract owns the change.
- Diagnostics, checkers, and presentation pages are non-authoritative unless
  a reviewed contract explicitly says otherwise. They must not silently become
  a production correction or subtraction.

## HGCer architecture

- Method A and Method B are independent local pion-background diagnostics.
  Keep their inputs, calculations, and conclusions independent; never tune one
  to agree with the other. Their comparison is diagnostic evidence, not a
  production-normalization loop.
- Existing staged pion components/refinement machinery remains separate from
  slow-proton subtraction. Retain provenance, unavailable states, and frozen
  upstream records rather than rebuilding results in a later consumer.
- Phase-D work is comparison/closure/validation infrastructure. Phase-E
  presentation consumes frozen diagnostic records and remains
  presentation-only. In particular, a Method-A presentation page must not
  acquire a Method-B or production-correction dependency.
- The current Phase-F.1 source builds a non-authoritative Method-A acceptance
  event/artifact contract after Method A, renders five setting-scope acceptance
  pages, and explicitly declares no production mutation or Method-B numerical
  dependency. Its source record is
  `phases/phase-f1-method-a-acceptance-farm-gate.md`.

## Implementation and source-review discipline

- Make source changes through one narrow phase/fix contract at a time. State
  the starting `test` HEAD, allowed and frozen files/interfaces, scientific
  owner, before/after behavior, preserved runtime path, tests, forbidden
  shortcuts, validation commands, diff audit, acceptance criteria, and stop
  boundary.
- An implementation summary is not proof. Inspect the actual diff and trace
  the applicable producer -> serializer/checkpoint -> checkpoint-first payload
  -> consumer -> renderer path before recording `SOURCE REVIEWED`.
- Keep only the artifacts needed for the current validation gate. Never mix
  artifacts from different commits or settings without identifying the
  mismatch.

## Farm-validation discipline

- The real KaonLT runtime exists only on the Jefferson Lab farm. ROOT/PyROOT,
  full-analysis, checkpoint/runtime integration, PDFs, and production behavior
  require supplied farm evidence; the user performs those runs.
- Validate one narrow gate at a time: targeted farm run, fresh artifacts,
  source/provenance and checker inspection, rendered PDF-page inspection, then
  `PASS` or one coherent repair. Broaden only after that gate passes.
- Each runtime record must identify evaluated commit, setting/gate, fresh
  artifact paths and hashes where available, inspections, and conclusion.
  `complete=true` is insufficient without this review.

## Current Phase-F.1 gate trap

The checked-in Phase-F.1 profile pins required analysis commit
`d656e15761970d7d612bb028d2746d077795e9ad` and permits only detached
collector/profile files after it. Do not use that profile to bless a later
analysis commit containing other source changes: its collector intentionally
reports an identity error. Reconcile the identity rule and the reviewed
analysis commit first; this is a source-gate issue, not farm evidence. See
`investigations/2026-09-11-phase-f1-source-identity-reconciliation.md`.
