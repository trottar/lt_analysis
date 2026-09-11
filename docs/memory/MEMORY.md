# Durable KaonLT project knowledge

This is curated durable knowledge, not a chronological activity log. Update it
only when source review, a decision, or runtime evidence changes what a future
session needs to know.

## Authority and reconciliation

- Repository `docs/memory/` is the authoritative development memory. Native
  Codex memory is supplemental associative context only.
- Treat `CURRENT.md`, this file, and `handoffs/CURRENT_HANDOFF.md` as a set;
  then verify the actual current `test` HEAD before substantial work.
- If memory, source, and runtime evidence disagree, resolve implementation from
  current source and runtime status from the newest supplied farm evidence.
- Commit messages, source tests, and validation-bundle tooling are not by
  themselves proof of an accepted phase or runtime result.

## Scientific-development ownership boundaries

- Keep scientific calculation, production correction/subtraction, runtime
  integration, diagnostics, checkers, and presentation-only plotting distinct.
- Preserve established random subtraction, proton subtraction, pion treatment,
  HGCer Method-A/Method-B independence, SIMC handling, yield calculation,
  acceptance/binning, efficiencies, L/T separation, and uncertainty behavior
  unless a narrow requested phase/fix explicitly owns the change.
- A diagnostic/checker or presentation artifact is non-authoritative unless its
  contract explicitly says otherwise. Do not let it silently mutate a
  production correction/subtraction path.

## Runtime and validation rules

- The real KaonLT runtime is on the Jefferson Lab farm. ROOT/PyROOT,
  full-analysis, farm, and runtime claims require supplied runtime evidence.
- Source review and local tests may establish only `SOURCE REVIEWED` (or a
  development-pending-farm state); neither is `CLOSED / RUNTIME VALIDATED`.
- Use one narrow validation gate at a time: targeted farm run, fresh artifacts,
  evidence inspection, then either `PASS` or one coherent repair.
- The user runs farm validation. A bundle contains only artifacts necessary for
  the gate being evaluated.

## Current initialization boundary

As of the 2026-09-11 repository-only initialization, accepted `test` source
baseline `dc4fc6283001739a487ec80068f951b0e388cae6` was recorded from the
inspected clean checkout. Approved phase history and newest runtime acceptance
have not been entered from an authoritative handoff. Do not backfill them from
commit subjects. See `CURRENT.md` for the observed checkout and
`handoffs/CURRENT_HANDOFF.md` for the exact continuation requirement.

## Closed decisions

- Repository-owned memory authority is recorded in
  `decisions/2026-09-11-repository-memory-authority.md`.

## Recurring traps

- A clean worktree or a commit on `test` does not make that commit the accepted
  test HEAD.
- A collector/profile that can package farm artifacts is not a farm result.
- An implementation summary is not proof: inspect the actual diff and the
  relevant runtime path after every implementation.
