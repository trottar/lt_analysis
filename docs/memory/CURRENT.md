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

Prior evidence/reviews own acceptance. The user-requested
[memory correction](phases/kaonlt-science-gate-displacement-memory-correction-2026-10-17-task-contract.md)
records the [assistant's October 15–17 workflow failure](investigations/kaonlt-science-gate-displacement-2026-10-17.md).
It is detached from farm science and clears no runtime gate.

- [Full-chain diagnostic](phases/baseline-pion-full-chain-observability-2026-10-15-task-contract.md):
  `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; read-only companion JSON/PDF
  for Center/lowe t1/t2/t3, t2 primary. Local scalar/fake-ROOT checks only;
  no capture, procedure-PDF, farm or scientific acceptance. Observer source is
  published; the next gate is the execution-path source audit below.
  Normal-launcher isolation stays blocked.

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

## Next Action

NEXT — independent ChatGPT actual-diff review of the
[all-low debug correction](phases/kaonlt-center-lowe-signed-pion-chain-owner-task-contract.md),
then user publication and pushed-state synchronization, followed by the narrow
farm-readiness decision for fresh Center/lowe seven-stage companion JSON/PDF
(canonical t2 primary, separate t1/t3). No farm command before that decision.
The user-supplied replacement contract supersedes the abandoned, unimplemented
Center/lowe owner contract at the same path. `DEVELOPMENT COMPLETE, FARM
VALIDATION PENDING`: launcher `-d` clears the legacy Left-only selector on the
single full low-epsilon child, uses ordinary available-setting selection and
exits before full high epsilon. Paired low/high capture and shared preflight
remain unchanged; scientific source and normal mode are frozen.

SOURCE VERIFIED at `1c7c4714344a0d9245071b87d4abe74de1c50682`:
`run_Prod_Analysis.sh` normal low/high, Right/Left/Center route reaches
`src/main.py` -> `rand_sub.py` -> `build_setting_t_bin_pion_parents` ->
`observe_chain_audit` (including parent reuse), with `alignment_outpath=OUTPATH`.
At that HEAD, `-d` was Left/lowe only; the local correction above changes it.
Background manifest -> component file map ->
`pion_component_shapes.py` opens template ROOT inputs in read mode; the existing
background directory's `mkdir -p` can be idempotent. This traced consumer is
input-only, not proof of all consumers or safe launcher execution. The old
helper remains a negative historical diagnostic, not an input classifier.

Separate historical normal-launcher isolation blocker, not an all-low debug
acceptance gate cleared here: when `-c` is enabled, launcher merge calls
(first populated low-epsilon merge branch) pass
`${LTANAPATH}/log/<setting>_kaon_Analysed_<sample>_Q4p4W2p74_lowe.err` to
`src/setup/mergeRootFiles.py:44`, whose `open(err_fout, 'w')` truncates/creates
the log. After clean/setup, `log` resolves to `${VOLATILEPATH}/log`, outside
the attempt. Failure invariant: `merge_error_log_outside_attempt`.
Earlier setup is independently unsafe if external `simc_gfortran` children or
configured background `worksim` need link repair; those link paths are external.
That isolation path lacks an attempt-owned log destination that survives
launcher cleanup and is respected by setup; LTANAPATH-only overlay cannot
provide it. Plain `-d` does not enable `-c`: the existing merge/delete block is
skipped by its flag guard, not by detection of existing ROOT inputs. Combining
`-c` with `-d` still deletes/remerges populated inputs and writes merge logs.
Setup and path configuration remain unchanged.

Local synthetic checks executed unchanged setup with pre-existing background
inputs for both epsilons/all phi settings, valid external SIMC links and
idempotent background mkdir; input bytes and external links were preserved.
Executing the actual merge log-write AST then overwrote an external synthetic
log on success and before injected child failure. This is local filesystem/source
evidence only; actual farm realpaths, installed default/no-mode Root behavior,
full downstream write coverage and ROOT rendering remain NOT VERIFIED.
`effective_charge_default_root_outpath_not_verified` and
`normal_launcher_downstream_write_coverage_unproven` remain `BLOCKED`;
no historical shared-output safety gate is cleared by the corrected classification.

Plain `-d` still resolves ltsep paths, creates skim/output directories, validates
external Sigma0 paths, runs `set_SymLinks.sh` (including possible external link
repair/background mkdir), reads both epsilon run lists/efficiencies/effective
charges, runs `set_sig_fortran.py` model rewrites, and performs both raw captures
and shared canonical preflight before the full low analysis. Existing ROOT
inputs do not remove these dependencies. Local source/stub tests establish
selection, inherited-selector neutralization, failure propagation and normal
low/high preservation only; ROOT, farm realpaths, downstream write coverage and
companion rendering remain NOT VERIFIED. Historical isolation/effective-charge
blockers remain separate and `BLOCKED`; this change grants no farm readiness.

No owner implementation, pathing redesign, or attestation/history cycle is
authorized here. Requested evidence remains
fresh farm-generated **Q4p4W2p74 / Center / lowe seven-stage companion JSON/PDF**,
canonical t2 primary and separate t1/t3 controls, with reviewed source, logs,
provenance, checker and page evidence. Resolve the baseline discrepancy before
returning to the interrupted roadmap.

## Success Criteria

Local integrity; no physics acceptance.

## Do Not Reopen Without New Evidence

Keep `no_empirical_residual`, zero legacy scales, detached/non-production Method A
and independent diagnostic/cross-check, numerically excluded Method B. Freeze
all subtraction, SIMC, cuts, templates, priors, binning, weights, efficiencies,
acceptance, yields/errors, L/T and cross sections; no promotion.

## Relevant References

See [roadmap](roadmap/STATUS.md).
