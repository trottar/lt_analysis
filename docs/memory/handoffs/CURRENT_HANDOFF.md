# Current KaonLT handoff

Prepared: 2026-09-11 during repository-memory initialization.

## Checkout to continue from

- Repository: `trottar/lt_analysis`
- Branch: `test`
- Observed current HEAD: `dc4fc6283001739a487ec80068f951b0e388cae6`
- Accepted `test` source baseline:
  `dc4fc6283001739a487ec80068f951b0e388cae6`
- Observed tracking point: `origin/test` at `eb253046f`; local `test` is one
  commit ahead.
- The worktree was clean before this memory-only change.

The accepted source baseline was initialized from the inspected clean `test`
checkout. Re-check `git status --short --branch` and `git rev-parse HEAD`
before work; source-baseline acceptance is not a claim of phase closure or
runtime validation.

## Closed/runtime-validated work

No `CLOSED / RUNTIME VALIDATED` phase is initialized from repository evidence.
Do not infer one from the historical commit graph.

## Source-reviewed or farm-pending work

No work is safely classified in this bootstrap record. The observed HEAD changes
these source/test files, but the diff alone does not establish their review or
farm state:

- `src/cuts/full_background_subtraction_plots.py`
- `src/cuts/pion_hgcer_method_a_acceptance_contract.py`
- `src/cuts/rand_sub.py`
- `testing/test_full_background_subtraction_plots.py`
- `testing/test_pion_hgcer_method_a_acceptance_contract.py`
- `testing/test_pion_hgcer_phase_f_runtime_contract.py`

The repository also contains a Phase-F.1 validation-bundle collector and
profile. They define a possible gate and requested artifact names; they are not
a runtime result.

## Active work and deferred issues

`ACTIVE` — Establish the remaining authoritative continuity baseline: current
phase/fix and newest supplied farm/runtime evidence.

`DEFERRED` — Historical phase classification is intentionally deferred until
the authoritative project handoff is available. This is a data-quality
boundary, not a judgment about the implementation.

## Exact next action and required evidence

`NEXT` — Request or receive the authoritative project handoff, then update
`CURRENT.md`, this handoff, `roadmap/CURRENT.md`, and any corresponding phase
or evidence record in the same change. Reconcile every supplied commit with the
accepted `test` source baseline first.

Required evidence for every farm/runtime record: evaluated commit, setting or
gate, fresh artifact paths (and hashes where available), inspection result, and
clear conclusion. Reconcile its commit with the live `test` HEAD before making
a status claim.
