# Current KaonLT development state

Last reconciled: 2026-09-11, durable-memory initialization and current-source
review only.

## Observed repository identity

- Repository: `trottar/lt_analysis`
- Branch: `test`
- Observed `HEAD` and `test` HEAD:
  `cf9c804b55a501af6839b8fb3b0358c56853c276`
- Observed upstream tracking point: `origin/test` at the same commit.
- Observed worktree: clean before this documentation-only change.

The previous bootstrap record named
`dc4fc6283001739a487ec80068f951b0e388cae6`. It is an ancestor of the live
checkout, not the current `test` HEAD. The `cf9c804...` commit adds the memory
infrastructure; the current Phase-F.1 analysis source is the source reviewed
in the ancestor range through `dc4fc628...`.

## Reconciliation baseline

`cf9c804b55a501af6839b8fb3b0358c56853c276` is the current repository
reconciliation baseline. It records source identity only. It is not an
acceptance, phase-closure, or farm-validation claim.

## Current source-established state

`SOURCE REVIEWED` — The current source implements a Phase-F.1 Method-A
acceptance diagnostic and its presentation/checker path. After Method A,
`rand_sub.py` builds a frozen acceptance event contract, serializes a
non-authoritative acceptance artifact, and sends its display payload to the
full-background renderer. The runtime fallback remains explicitly
non-authoritative and declares no production-object mutation, refinement, or
Method-B numerical dependency. The renderer appends five setting-scope
acceptance pages; it does not alter the subtraction.

`ACTIVE` — The Phase-F.1 farm-review gate needs source-identity reconciliation
before it can validly evaluate the live checkout. The declarative profile
requires analysis commit `d656e15761970d7d612bb028d2746d077795e9ad` and
permits only collector/profile files after it. Live `test` also contains later
changes to `src/cuts/full_background_subtraction_plots.py`,
`src/cuts/pion_hgcer_method_a_acceptance_contract.py`, and
`src/cuts/rand_sub.py`; the collector is designed to reject those as
unexpected committed files. See
`investigations/2026-09-11-phase-f1-source-identity-reconciliation.md`.

## Runtime-established state

No farm artifact record is present in `docs/memory/evidence/` and none was
supplied for this initialization. There is therefore no locally supported
`CLOSED / RUNTIME VALIDATED` phase, no runtime acceptance of Phase-F.1, and no
basis to classify older phases by farm result. The collector/profile and source
tests are tooling/source evidence only.

## Frozen architecture to preserve

- Random subtraction, slow-proton subtraction, pion-background treatment,
  HGCer diagnostics, SIMC comparison, yield extraction, and final
  cross-section analysis have separate scientific ownership.
- HGCer Method A and Method B are independent diagnostic/cross-check paths.
  Neither is a permission to tune the other or to alter production pion
  subtraction.
- Phase-D comparison/closure records and Phase-E presentation consume stored
  upstream diagnostics. Presentation is non-authoritative and must not
  recompute or mutate production physics.
- A persisted diagnostic must retain its producer, serializer/checkpoint,
  checkpoint-first payload, consumer, and renderer provenance.

## Deferred history and exact next step

`DEFERRED` — Import authoritative historical farm evidence before assigning
runtime status to the broad background program or named earlier phases.

`NEXT` — Create one narrow Phase-F.1 source-identity reconciliation contract:
choose and review the analysis commit that the farm profile will evaluate, then
make the profile/collector identity rule agree with that commit before a
targeted farm run. Keep the required Phase-C/Phase-D/correction/Method-A
artifacts and rendered PDF-page inspection in that later gate. Independently
import any authoritative historical farm records with their evaluated commit,
setting, fresh artifact paths, inspection, and conclusion.
