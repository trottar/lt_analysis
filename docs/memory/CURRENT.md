---
memory_schema: 3
---
# Current KaonLT development state

## Active Objective

`ACTIVE`: **TEMPORARY roadmap interruption** to explain the full existing
baseline SIMC pion-background model versus final signed-event subtraction
discrepancy. Focus: **Q4p4W2p74 / Center / lowe / canonical t2**; adjacent t1/t3
and setting-wide controls retain their distinct scopes. The user's visual
comparison places deterioration during canonical-t model construction/refinement,
before final pion subtraction; the setting-wide model appears substantially more
sensible. This is a user observation, not independently accepted model quality
or an exact first-bad operation. The full setting-wide procedure remains in
`Center_kaon_rand_sub_Q4p4W2p74_lowe_pion-fit-debug.pdf`: the older 67-page
reference versus October's 30-page main PDF reflects presentation routing,
not disappearance of global fit machinery. See the
[October 19 correction](investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md#global-to-fine-bin-procedure-audit-and-user-correction--october-19).
The user also reported earlier working behavior. The October 17–18
farm evidence and prior review are now recorded in the
[investigation checkpoint](investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md#low-epsilon-farm-evidence-and-scientific-checkpoint--october-1718).
No matched canonical-t regression pair, historical culprit or wider failure verified.

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

Prior evidence/reviews own acceptance. The user-requested
[memory correction](phases/kaonlt-science-gate-displacement-memory-correction-2026-10-17-task-contract.md)
records the [assistant's October 15–17 workflow failure](investigations/kaonlt-science-gate-displacement-2026-10-17.md).
It is detached from farm science and clears no runtime gate.

- [Full-chain diagnostic](phases/baseline-pion-full-chain-observability-2026-10-15-task-contract.md):
  `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING` for independent validation.
  RUNTIME VERIFIED (supplied artifacts/prior review, recorded through the
  [October 18 contract](phases/kaonlt-low-epsilon-baseline-science-memory-checkpoint-2026-10-18-task-contract.md)):
  published all-low debug source `f1d1984b5e0fca9ba3c5666477996af783423db7`
  ran on the farm with `DEBUG_EXIT_STATUS=0`; Center/Left lowe PDFs and a
  Center canonical t1/t2/t3 seven-stage JSON/PDF (19 pages) were captured.
  No populated Right/lowe result established. Producer status remains
  `captured_pending_independent_validation`, completeness `captured`, PDF
  `available`, `observed_source_HEAD=not_provided`. Farm checkout HEAD was
  explicitly verified before the run; that does not fill the artifact's
  missing identity. The reviewed eight-member `KaonLT_loweps_20261018_093137.zip`
  is not full independent provenance/checker or production acceptance.
  Normal-launcher isolation remains separately blocked.

- [Attempt-output isolation](phases/kaonlt-pion-full-chain-attempt-output-isolation-2026-10-16-task-contract.md):
  `BLOCKED`: the historical helper reports
  `post_clean_background_output_mkdir_outside_attempt` unconditionally for an
  external background directory. The revised owner audit below distinguishes
  existing input setup from actual writers; this negative helper grants no PASS.
  Source checks retain the frozen ancestor and record current HEAD.
  No safe normal-launcher proof or farm readiness.

- [Safe path interface](phases/kaonlt-pion-full-chain-safe-output-path-interface-2026-10-16-task-contract.md):
  `BLOCKED`: `effective_charge_default_root_outpath_not_verified`.
  [Resolver attestation](phases/kaonlt-pion-full-chain-ltsep-default-root-attestation-2026-10-16-task-contract.md):
  `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; installed behavior NOT VERIFIED;
  both isolation blockers remain.
  [Receipt repair](phases/kaonlt-pion-ltsep-receipt-canonical-globus-repair-2026-10-16-task-contract.md):
  `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; canonical-role policy preserved.
  [History compatibility](phases/kaonlt-pion-receipt-regression-history-compatibility-2026-10-16-task-contract.md):
  `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; explicit mutual helper/test
  successor lists repair the committed-history rejection without changing frozen
  source checks. Pre-existing local report: 47 WSL synthetic/history/isolation tests plus eight
  Windows Python owner regressions passed (55 total); initial WSL import lacked
  NumPy. Not rerun by this memory correction. The October 17 receipt/history
  repair is published at `dceaadb87f303c0c406978b2d85a2728b5bd2056`;
  publication establishes no farm readiness or completed execution owner.
  Supplied JLab attempt stopped at `receipt_directory_symlink`; no resolver acceptance.

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

Historical normal-launcher/attempt isolation remains separately `BLOCKED`:
`merge_error_log_outside_attempt` (`-c` merge writes through external `log`),
`effective_charge_default_root_outpath_not_verified`, and
`normal_launcher_downstream_write_coverage_unproven`. External SIMC/background
link repair can also write outside an attempt. Plain `-d` does not invoke
`-c`; its successful observed run does not clear those broader gates. Setup,
path configuration and all scientific source remain frozen. The prior
[debug contract](phases/kaonlt-center-lowe-signed-pion-chain-owner-task-contract.md)
and [workflow investigation](investigations/kaonlt-science-gate-displacement-2026-10-17.md)
retain historical scope; no owner/pathing/attestation cycle is authorized here.

## Next Action

NEXT — read-only comparison of setting-wide versus canonical-t production
input populations, template-source/shape fallbacks, amplitude scales and
fit-stage/validation outputs for Q4p4W2p74 Center lowe t1/t2/t3, with Left lowe
as control, using existing October 17–19 evidence first. Establish what changes
before subtraction; discriminate local selection/normalization, global-shape
fallback versus candidate genuine t-shapes, and unconstrained staged/joint
amplitudes. Quantify available support/covariance and identify missing evidence
and the smallest directly evidenced intervention. No source edit or farm command
is authorized. After independent ChatGPT actual-diff review, user publication
and pushed-state synchronization, proceed directly to this same audit.

RUNTIME VERIFIED (supplied October 17–18 evidence/prior review): Center/t1 has
staged-anchor versus full-range SIDIS tension; Center/t2 has a high-leverage
negative random contribution, negative neutron LS numerator (~-61.674), staged
`B_n=0`, refined `B_n~0.000202` and `w0~182.82`. Its model pion integral
~0.2529 versus event-applied background ~6.5562 illustrates failed signed-event/
model agreement in each Center t bin despite model algebra/subtraction arithmetic
closure. These are full-range integrals, not extracted yields or cross sections.
Left's setting-wide neutron amplitude is positive but fit quality is poor;
matching source-resolved Left canonical-t diagnostics remain NOT VERIFIED.
See the investigation for the transition-specific findings and limits.

SOURCE VERIFIED at `1975c3d8656253e262ff7e18306a93ca800b6027`: plain `-d` uses ordinary available lowe
settings, skips `-c` merge and exits before full highe. Analysis physics is
unchanged. `src/main.py` calls the pion-shape loader without `t_bins/phi_bins`;
canonical-t scopes can fall back to setting-wide shapes. The Center artifact
reports no constructed t-binned SIMC shapes (supplied evidence). This possible
limitation is not a proven cause or permission to activate new shapes.
SOURCE VERIFIED: the setting-wide fit is diagnostic-only in active `t_bin`
scope; canonical-t runs the full fitter on local inputs with parent alignment,
without transferred global amplitude/covariance constraints. Enabled joint
regularization centers on that fit's local staged coefficients. Authoritative
`(t,phi)` consumers apply the t-parent event weight without a full local refit.
User-requested future intent: preserve the complete staged three-component
procedure, anchor t to global while allowing genuine large t dependence, then
consider support/covariance-conditional children anchored to their own t parent.
INFERENCE: this may help; shared-event dependence and model adequacy require
audit. No hierarchy, prior width or fallback algorithm is approved. Numerical
cause, efficacy and full farm acceptance remain NOT VERIFIED. All existing
normalization/weights/cuts/priors/templates and signed events remain frozen.

## Success Criteria

Local integrity; no physics acceptance.

## Do Not Reopen Without New Evidence

Keep `no_empirical_residual`, zero legacy scales, detached/non-production Method A
and independent diagnostic/cross-check, numerically excluded Method B. Freeze
all subtraction, SIMC, cuts, templates, priors, binning, weights, efficiencies,
acceptance, yields/errors, L/T and cross sections; no promotion.

## Relevant References

See [roadmap](roadmap/STATUS.md).
