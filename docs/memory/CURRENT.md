---
memory_schema: 3
---
# Current KaonLT development state

## Active Objective

`ACTIVE`: **TEMPORARY roadmap interruption** for the baseline-pion investigation,
now evaluated through full missing-mass reconstruction and signal-region relevance
for **Q4p4W2p74 / Center / lowe / canonical t2**. Establish common selection,
coordinates, binning and charge/efficiency/acceptance normalization for data,
proton contamination, baseline pion background and K+Lambda/K+Sigma0 SIMC.
Inspect full-MM spectra, signal-region peaks/tails/residuals, then authoritative
background-subtracted K+Lambda data/SIMC yields in corresponding `(t,phi,epsilon)`
bins. Unity is a model goal, never an area-matching constraint. Absolute
normalization remains separately blocked below.

The [October 20 checkpoint](investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md#full-mm-reconstruction-scientific-checkpoint--october-20)
records supplied farm findings, the histogram crosswalk and the user's scientific
correction. Large control/full-range disagreement alone does not establish large
contamination under Lambda or signal bias. Historical first-bad SHA and physical
root cause remain NOT VERIFIED. Earlier observations and procedure routing remain
in the investigation; no earlier accepted scope is reopened.

Prospective Method-A development/promotion and normal E.8/F.6 forward work remain
paused pending baseline resolution; no new phase. After explanation/validation
of the existing subtraction, or an independently reviewed, user-pushed and
user-farm-validated minimum correction, resume [the roadmap](roadmap/STATUS.md)
at its interrupted position.

## Current Work Item

`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`: local
[fast-fail correction](phases/kaonlt-matched-cohort-fast-fail-science-unblock-2026-10-24-task-contract.md)
from `c6c6963a00cfc82965eb2414968093a3da4a0fbb`; actual-diff review pending.
SOURCE VERIFIED: caller/cwd/module/config/private-output and input-availability
gates precede scientific hashing/copying. Authorized root aliases retain both
path identities. Full shared byte checks remain before collection and after
cleanup without a duplicate final scan. Status records stage timings. See the
[source/preservation review](investigations/kaonlt-matched-cohort-farm-readiness-failure-2026-10-22.md#fast-fail-source-and-preservation-review--october-23).
Scientific source, frozen resolver, helpers and launchers are unchanged.

Runtime acceptance remains `BLOCKED`: the seventh supplied pre-child failure at
`b05ad2a1d9af241527a7d6ce90196eb591188d19` was
`caller_path_invalid:src/setup/findEffectiveCharge.py`. Its exact farm alias,
status/timestamps, preservation and reported hours remain NOT VERIFIED.
No successful child, matched ZIP or fresh science is supplied. October 20
remains the accepted scientific checkpoint; local tests grant no farm readiness.

The [matched-cohort owner](phases/kaonlt-pion-matched-cohort-attempt-isolation-2026-10-21-task-contract.md)
retains private output/input copies, stable external roles and VOLATILEPATH,
default/Plot caller probes, both epsilon prepasses, diamond-selection parity,
optional unconfigured Sigma0 fallback and strict configured-file checks.
Only its child enables both observers. Original inputs, shared outputs,
external links, installed/copied packages and fresh six-parent A/B/event
joins/closure must pass before delivery. Launchers, isolation helpers, installed
`ltsep` and production source remain unchanged. Local synthetic checks grant
no farm readiness or scientific acceptance.

Frozen `shift_prep` writes `MM_shift`/`t_shift` only in private copies. Original
hashes/inodes and runtime-container identities remain separate; fresh A/B and
packets must match the latter. Real ROOT writes, current farm inputs and
preservation remain NOT VERIFIED here.
See the [October 22 source checkpoint](investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md#matched-cohort-private-acquisition-source-checkpoint--october-22).

The same CLI's explicit offline mode verifies the new bundle, reconstructs six
B controls, fits seven fresh C strengths per parent and replays signed events
with distinct A and fixed protected-applied `K_B` for B/C. Original/runtime
identities, event/cut/template/fit provenance, source/phi-child closure, MM/tail
statistics, production conventions and paired fixed-weight errors remain in the
detached six-parent JSON/PDF. Historical pins stay strict; no historical C reuse.

No matched farm bundle or fresh physical C result is available here. Historical
A/B/C remains immutable and cohort-specific. The earlier failed capture at
`04f1ec403a9555634e61c229aa347cb76d0c78d3` and published
[prepass boundary repair](phases/kaonlt-pion-event-capture-prepass-boundary-fix-2026-10-21-task-contract.md)
retain their recorded limits; failed artifacts supply no event acceptance.
The [independent event-validation objective](phases/kaonlt-pion-adaptive-c-event-validation-2026-10-21-task-contract.md)
remains: full-sample fits/priors/protected models do not provide a held-out fitted
estimator. Training-only refits, physical bias truth and MSE selection remain
unavailable. No strength is selected or promoted.

The [controlled template comparison](phases/kaonlt-pion-t-template-controlled-comparison-2026-10-19-task-contract.md)
and [execution gate](phases/kaonlt-pion-t-template-science-execution-gate-2026-10-19-task-contract.md)
were source reviewed, user-pushed and farm-run at
`6d13abb4681b938c56def2a1f455244b0880096c` (SOURCE VERIFIED / RUNTIME VERIFIED,
scope attributed to supplied artifacts and prior review in the October 20
contract). Center/Left lowe canonical t1/t2/t3 detached JSON/PDFs and receipts
were accepted for that diagnostic scope. Candidate t-local templates do not
establish a scientific preference or production adoption; baseline remains
authoritative. Center/t2 staged validation worsens, and both protected-delta
candidates are rejected with applied pi_delta zero. Detailed metrics, evidence
limits and nomenclature are in the checkpoint. No local farm/artifact rerun.

Prior evidence/reviews own acceptance. The user-requested
[memory correction](phases/kaonlt-science-gate-displacement-memory-correction-2026-10-17-task-contract.md)
records the [assistant's October 15–17 workflow failure](investigations/kaonlt-science-gate-displacement-2026-10-17.md).
It is detached from farm science and clears no runtime gate.

- [Full-chain diagnostic](phases/baseline-pion-full-chain-observability-2026-10-15-task-contract.md):
  `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING` for independent validation.
  RUNTIME VERIFIED only as supplied/prior-reviewed in the
  [October 18 contract](phases/kaonlt-low-epsilon-baseline-science-memory-checkpoint-2026-10-18-task-contract.md):
  all-low debug `f1d1984b5e0fca9ba3c5666477996af783423db7`,
  `DEBUG_EXIT_STATUS=0`, Center/Left lowe PDFs and Center t1/t2/t3 seven-stage
  JSON/PDF (19 pages). No populated Right/lowe result. Producer status remains
  `captured_pending_independent_validation`; `observed_source_HEAD=not_provided`.
  Checkout verification does not fill missing artifact identity. The reviewed
  eight-member `KaonLT_loweps_20261018_093137.zip` establishes neither full
  independent provenance/checker acceptance nor production/normal-launcher isolation.

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

NEXT — preserve the scientific target: fresh same-cohort Center/Left lowe
A/B and signed-event packets for canonical t1/t2/t3, then the existing fresh
seven-strength C fitting/replay and six-parent Lambda-window MM/variance
JSON/PDF comparison. First independently review the complete fast-fail diff and
CLI execution/preservation path, then user publication and pushed-state
synchronization before acquisition-readiness review. No farm command or new
mandatory smoke cycle follows from local tests. Hard-stop at the first failed invariant;
no automatic retry or failed-artifact acceptance. No farm/runtime acceptance,
held-out estimator, bias truth, optimal strength or promotion follows. Preserve
full-MM, normalization and authoritative signal-yield objectives without a
Lambda/SIMC unity constraint.

Hierarchical amplitudes and child refitting remain deferred; no prior width,
covariance model or fallback is approved. The investigation retains prior
construction, transfer and signed-application findings and their limitations.
The template diagnostic's single opt-in remains detached; it does not change
production templates, alignment, fits, subtraction or yields.

## Success Criteria

Local integrity; no physics acceptance.

## Do Not Reopen Without New Evidence

Keep `no_empirical_residual`, zero legacy scales, detached/non-production Method A
and independent diagnostic/cross-check, numerically excluded Method B. Freeze
all subtraction, SIMC, cuts, templates, priors, binning, weights, efficiencies,
acceptance, yields/errors, L/T and cross sections; no promotion.

## Relevant References

See [roadmap](roadmap/STATUS.md).
