---
memory_schema: 3
---
# Current KaonLT development state

## Active Objective

`ACTIVE`: **TEMPORARY roadmap interruption** to explain the full existing
baseline SIMC pion-background model versus final signed-event subtraction
discrepancy. Focus: **Q4p4W2p74 / Center / lowe / canonical t2**; adjacent t1/t3
and setting-wide controls retain their distinct scopes. The user compared
`Center_kaon_rand_sub_Q4p4W2p74_lowe.pdf` (apparently reasonable model) with
`Center_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction.pdf` (problematic
final subtraction), and reported earlier working behavior. MEMORY ONLY:
no matched canonical-t regression pair, historical culprit or wider failure verified.

The [October 15 contract](phases/baseline-pion-full-model-memory-reanchor-2026-10-15-task-contract.md)
and [investigation](investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md)
restore the full `pi_n`/`pi_sidis`/`pi_delta` question. `B_n=0` and protected-delta failure are candidates, not an established
first historical defect, root cause or sole target.

Prospective Method-A development/promotion and normal E.8/F.6 forward work and
interpretation are paused pending baseline resolution; no new phase. After
explanation and validation of the existing subtraction, or a necessary minimum
correction independently reviewed, user-pushed and user-farm validated,
**resume [the roadmap](roadmap/STATUS.md) at its interrupted position**.
All existing statuses/scopes remain intact.

## Current Work Item

Prior evidence/reviews own acceptance.

- E.8.2 Left/lowe: `CLOSED / RUNTIME VALIDATED` per
  [closure](evidence/e8-2-left-lowe-runtime-scientific-closure-2026-10-06.md).
  Final baseline `after_pion_final/MM_0`, stored `Y0`; stage integrals are diagnostics.
- E.8.3 Left/lowe and authority-page repair: `CLOSED / RUNTIME VALIDATED` per
  [closure](evidence/e8-3-left-lowe-runtime-closure-2026-10-07.md).
  Historical F.6.1 is separate from current F.6.3/E.8.4.
- E.8.2 baseline-five `Y0/MM_0`: `CLOSED / RUNTIME VALIDATED` per
  [October 9](evidence/e8-2-canonical-five-baseline-yield-runtime-closure-2026-10-09.md):
  135/135 data and SIMC gates, original PDF limits. The
  [QA](evidence/e8-baseline-five-extracted-yield-scientific-qa-2026-10-11.md)
  adds no physics acceptance or proof of all component-decomposition physics.
- F.4.Refresh.2 and Left/lowe F.6.3/E.8.4: `CLOSED / RUNTIME VALIDATED` within
  [accepted scopes](evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md):
  live-cache parity, real child changes, signed parent preservation.
  Fix.5.7/Fix.5.8 retain provenance and
  [97-page legibility closure](evidence/e8-4-fix5-8-left-lowe-runtime-closure-2026-10-03.md).
  Earlier F-stage closures and source-reviewed hardening remain in the roadmap.
- [October 12 and October 13 addendum](evidence/e8-current-snapshot-contribution-runtime-closure-2026-10-12.md):
  `CLOSED / RUNTIME VALIDATED` only for five-F.1 candidate materialization,
  read-only seven-input preflight, eight-page Left/lowe t0/t1/t2 diagnostic at
  `75756d7b32767f18949386c4c5864c756d5e85f6`; internal coherence PASS,
  canonical-owner equivalence BLOCKED. Fingerprint/Neff failures resolved only
  here; no authority replacement. October 13 at
  `c8e830f38cf59a22b7b6d9753b1a847fe1124e70`: Left/lowe presentation only,
  eight pages inspected, 5–8 physical 0.8–1.4 GeV/distinct UF/OF; numerical
  objects/October 12 acceptance unchanged, no Center MM/PDF closure.
- [Paired extension](phases/e8-current-snapshot-paired-lowe-phi-mm-extension-task-contract.md):
  `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; separate Left/Center lowe
  aggregates, 16 pages, default eight-page Left-only path preserved. Later paired
  progress MEMORY ONLY; no independent closure/receipts/hashes/16-page acceptance.
  Left-only conclusions stay setting-local; priority interrupted, not failed.
- [Component diagnostic](phases/baseline-pion-center-t2-neutron-anchor-diagnostic-task-contract.md):
  `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; capture/cause NOT VERIFIED.
  [Farm-owner](phases/baseline-pion-center-t2-anchor-farm-owner-task-contract.md)
  stays `BLOCKED` before implementation on output isolation; the tracked
  [isolation proposal](phases/baseline-pion-center-t2-attempt-output-isolation-task-contract.md)
  is preserved, not launched, prerequisite only for a future normal-launcher run.
  These contracts do not define NEXT.

## Verified State

[Durable rules](MEMORY.md) retain Gate-3/HGCer/transfer roles, F.6.2 acceptance/MM
scope and t1-only continuity. Absolute mis-ID/proxy/hardware/RF remain NOT VERIFIED.
`YA/MM_A` is reweighted baseline; standalone Method-A-only final yield/MM is
NOT VERIFIED. Signed parents, no child normalization, existing errors and
`MM_0/MM_A`, `Y0/YA` consumers are frozen; no child/yield closure inferred.

## Source / Evidence Identity

Linked records retain source/farm/artifact identities. Frozen reviewed F.4 SHA-256:
`79e7ceda7221cbeeead4ed5bc306b0e0e670741a27beaa980e22349c555e96d7`.

## Blockers

Private canonical-five Method A and final E.8/F.6.4 remain `BLOCKED`:
`f6_3_scientific_equivalence_mismatch:f4:$.parents[0].absolute_baseline_parent_sum`.
October 8 owner `ebf3a0be2f50c788c22defe3601fd5b7fee48644`: F.2/F.3 exact PASS,
F.4 exact FAIL in three Left-lowe parents before analysis; no retry.
Roadmap retains the historical second-full-run page/artifact failure
and separate [lineage-only preflight](evidence/e8-4-canonical-five-lineage-preflight-runtime-closure-2026-10-05.md),
`CLOSED / RUNTIME VALIDATED` (`analysis_started=false`). F.2 handoff/staging, F.1 preflight repairs and detached
reproducibility diagnostic remain `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`;
microscopic F.1/F.4 archaeology `DEFERRED`. Production correctness unproven.

Absolute SIMC separately `BLOCKED`:
`SIMC_normfac_luminosity_and_charge_units_not_source_proven`.
`iter_weight * normfac / Ncontribute` units unproven; two absolute-SIMC page
families unavailable, not data/identity/yield/parent closure. No scale claim.

## Next Action

NEXT — after independent ChatGPT actual-diff review, user commit/push and
ChatGPT pushed-state synchronization, conduct a read-only **PDF-first** comparison
of earlier well-behaved component/model procedure and problematic full/final
subtraction pages. Establish scope, run/input/cut, normalization, MM-coordinate/bin
and provenance identities; trace all three SIMC templates -> staged/joint/protected
fits -> final control/kaon models -> MM transfer ratio -> signed pion-control event
application -> full subtracted MM/yields. Locate the earliest supported discrepancy.
`B_n=0`/protected-delta failure are candidates only. If PDFs/numerical records are
unavailable after searching accessible sources, identify the single earliest decisive
missing evidence. No default farm/implementation, plotting or isolation work. Remain `ACTIVE`
until resolved; return under the validation condition above.

## Success Criteria

Local integrity; no physics acceptance.

## Do Not Reopen Without New Evidence

Keep `no_empirical_residual`, zero legacy scales, detached/non-production Method A
and independent diagnostic/cross-check, numerically excluded Method B. Freeze
all subtraction, SIMC, cuts, templates, priors, binning, weights, efficiencies,
acceptance, yields/errors, L/T and cross sections; no promotion.

## Relevant References

See [roadmap](roadmap/STATUS.md).
