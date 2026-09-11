# Current KaonLT handoff

Prepared: 2026-09-11 after durable-memory initialization and current-source
review. This handoff makes no farm-validation claim.

## Resume from this checkout

- Repository: `trottar/lt_analysis`
- Branch: `test`
- Live `HEAD` / `test` / `origin/test` at review:
  `cf9c804b55a501af6839b8fb3b0358c56853c276`
- Worktree was clean before the documentation-only memory update.
- The earlier memory bootstrap recorded `dc4fc6283001739a487ec80068f951b0e388cae6`.
  It is an ancestor, not the current head. Re-run `git status --short --branch`
  and `git rev-parse HEAD` before any work.

## Source-established architecture

- Preserve the independent random, slow-proton, pion-background, HGCer,
  SIMC, yield, and cross-section paths. Diagnostics/presentation do not own
  production correction.
- Preserve Method-A / Method-B independence. Phase-D is comparison/closure
  infrastructure; Phase-E is frozen-record presentation-only work.
- `SOURCE REVIEWED` — Phase F.1 is a non-authoritative Method-A acceptance
  diagnostic. `rand_sub.py` creates and serializes its acceptance contract
  after Method A, then routes a display payload to the full-background
  renderer. Its fallback reports unavailable state without production mutation,
  refinement, or Method-B numerical dependency. The renderer appends five
  setting-scope pages. Details are in
  `phases/phase-f1-method-a-acceptance-farm-gate.md`.

## Current active gate

`ACTIVE` — Reconcile the Phase-F.1 validation-profile source identity with the
live source before asking for a farm run.

`BLOCKED` — The profile requires analysis commit
`d656e15761970d7d612bb028d2746d077795e9ad` and allows only three detached
collector/profile files after it. The current branch has subsequent changes to
`src/cuts/full_background_subtraction_plots.py`,
`src/cuts/pion_hgcer_method_a_acceptance_contract.py`, and `src/cuts/rand_sub.py`.
The collector is designed to reject this committed range. This is an
intentional source-identity guard, not an indication that a farm run failed.

The existing profile expects five declared settings and these frozen artifact
classes: Phase-C checkpoint, Phase-D checkpoint, parent-preserving correction,
Method-A acceptance artifact, procedure PDF, and page manifest. A valid later
farm review must inspect provenance/checker gates and the selected rendered
pages; a collector archive alone is insufficient.

## Runtime evidence and historical status

No completed farm artifact record is present locally or supplied for this
initialization. Therefore no phase is recorded as `CLOSED / RUNTIME VALIDATED`,
including Phase F.1. Do not backfill older phase status from Git subjects,
source, tests, or this handoff.

`DEFERRED` — Import the authoritative history for previously completed farm
gates, if any. Each imported record needs the evaluated commit, setting, fresh
artifact paths/hashes, inspection, and conclusion.

## Exact next action

`NEXT` — Write one narrow Phase-F.1 source-identity reconciliation contract:
select the reviewed analysis commit to evaluate and make the profile/collector
identity rule agree with it. Only then perform one targeted farm gate using
fresh declared artifacts and PDF-page inspection. Do not modify Method A,
Method B, or production subtraction as part of that reconciliation unless a
separate approved contract explicitly requires it.
