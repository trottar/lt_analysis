# E.8 full-analysis procedure roadmap

## Governing decision

E.8 is the complete visual audit of how the kaon missing-mass spectrum and signal-region yield evolve through the analysis. Every substantive stage exposes the authoritative input spectrum, component or treatment applied, and authoritative output spectrum. The presentation culminates in final canonical `(t,phi)` missing-mass spectra and extracted yields.

E.8 is presentation-only. It consumes authoritative upstream runtime objects, persisted snapshots, accepted detached Method-A artifacts, or later F.6.3 branch outputs. It never recomputes a fit, factor, correction, normalization, or yield merely for plotting.

The baseline production branch remains authoritative and unchanged. Method B remains diagnostic/cross-check only and is numerically excluded from every Method-A correction/application. Method A uses exactly the accepted F.4 parent-preserving correction; no child `(t,phi)` is independently renormalized. The accepted F.6.2 JSON SHA-256 and fingerprints remain frozen.

The accepted active background profile is `no_empirical_residual`; it forces
both empirical residual-background scales to zero. Legacy empirical residual
Fit 1/Fit 2 remain dormant historical source machinery and are not an E.8,
Method-A, or current production/presentation stage.

## Accepted narrow E.8.1 evidence

E.8.1.Fix.5 and E.8.1.Fix.6 are `CLOSED / RUNTIME VALIDATED` only for the fresh `Q4p4W2p74 / Left / lowe` persisted-overlay geometry and profile/provenance gate. See [the Fix.6 Left/lowe runtime closure](../evidence/e8-1-fix6-left-lowe-runtime-closure.md). The remaining E.8.1 canonical-five expansion is `DEFERRED` by user decision. That deferral is neither a failure nor a canonical-five closure.

## Approved dependency order

    E.8.2 baseline full-analysis audit
        -> E.8.3 detached Method-A reweighting audit
        -> F.6.3 parallel full procedure plus Method A
        -> E.8.4 baseline-versus-Method-A production-impact audit
        -> final E.8 closure
        -> F.6.4 explicit production-promotion decision

## E.8.2 — baseline full-analysis stage audit

`CLOSED / RUNTIME VALIDATED` only for Q4p4W2p74 / Left / lowe at
`706ae708ae69be0585472982a4df3d87c63f186f`, from the
[accepted narrow closure](../evidence/e8-2-left-lowe-runtime-scientific-closure-2026-10-06.md).
SOURCE VERIFIED: final canonical MM consumes after-baseline-pion
`after_pion_final` / baseline `MM_0` with stored authoritative `Y0`.
Diagnostic Lambda-window stage integrals, including `PI`, are not final
extracted yields and do not identify a further subtraction stage.

### E.8.2a — prompt/random subtraction

Show prompt kaon missing-mass input, random contribution removed, after-random missing-mass output, and a difference or overlay where useful.

### E.8.2b — dummy subtraction

Show after-random input, normalized dummy contribution, after-dummy output, and a difference or overlay where useful.

### E.8.2c — slow-proton subtraction

Show after-random/dummy kaon input, proton contamination estimate, proton-cleaned kaon output, before/after comparison, and existing PID/weight diagnostics. `_process_yield_data_tree()` already applies accepted proton cleaning during initial `(t,phi)` filling. Existing per-cell raw/random/dummy snapshots must not be relabeled as pre-proton. E.8.2 must capture or consume a true authoritative pre-proton snapshot rather than reconstructing or guessing one.

### E.8.2d — baseline pion treatment

Show proton-cleaned kaon input, the actual baseline pion background `B_pi^0`
using accepted `w0`, and the baseline clean-kaon output
`K_0 = K_proton - B_pi^0`, with before/after and signed difference where
useful. The baseline contribution remains `b_j^0 = s_j * w0_j`.

### E.8.2e — final baseline canonical `(t,phi)` missing mass

For every canonical t parent show all nine canonical phi children, including
explicit empty/invalid status. Every populated final panel exposes canonical-t
identity, phi index/range, final baseline clean-kaon `MM_0(t,phi)` actually used
for yield extraction, Lambda signal/integration window, final baseline extracted
yield, and statistical uncertainty. There is no empirical-residual stage between
baseline pion subtraction and this final spectrum.

### E.8.2f — baseline stage-yield audit

Expose the existing authoritative per-`(t,phi)` progression: prompt, after
random, after dummy, after proton, after baseline pion/final baseline clean
sample, and authoritative final baseline `Y0(t,phi)`. A renderer does not
recalculate a yield when its producer can persist or provide it.

## E.8.3 — detached Method-A reweighting audit

`CLOSED / RUNTIME VALIDATED` only for Q4p4W2p74 / Left / lowe at farm source
`63eb2bf66d98baaffa39e88cf1db7ac37397d364`, from supplied artifacts/prior
independent ChatGPT review. See the
[closure evidence](../evidence/e8-3-left-lowe-runtime-closure-2026-10-07.md).
Owner/provenance/preservation passed; repaired authority page 65 shows all
seven full provenance values, and pages 66–72 are pixel-identical to the prior
accepted render. The detached/non-production audit still consumes accepted
F.4/F.5/F.6.1/F.6.2 results only, separately from current F.6.3/E.8.4.
No scientific redesign, canonical-five closure or Method-A promotion follows.

Method A operates on the accepted pion background after the proton-cleaned
baseline picture; it does not invoke legacy empirical residual fits.

### E.8.3a — reweighting operation

Explicitly show `b_j^0 = s_j * w0_j` to `b_j^A = s_j * w0_j * C_j`. Method A changes only the accepted pion event contribution by the accepted F.4 parent-preserving factor. No child-by-child renormalization is allowed.

### E.8.3b — baseline versus reweighted pion background

For every canonical t parent show baseline `B_pi^0(MM)`, Method-A-reweighted `B_pi^A(MM)`, overlay, signed difference `Delta B_pi = B_pi^A - B_pi^0`, and ratio `R_pi = B_pi^A / B_pi^0` only where numerically well-defined and presentation-safe. This comparison is a first-class E.8 result, not a legend implication.

### E.8.3c — canonical redistribution and parent closure

For every canonical t show all nine `B_pi,tphi^0`, `B_pi,tphi^A`, and `Delta B_pi,tphi` children. Display stored/accepted closure and relative residual for `sum_phi B_pi,tphi^A ~= sum_phi B_pi,tphi^0`. Never independently normalize a child.

### E.8.3d — F.6.2 acceptance/HGCer explanation

Retain accepted explanatory diagnostics: L/B/A normalized shapes; `delta x xptar`; `delta x yptar`; missing-mass x acceptance maps; support/OOD; effective statistics; and kaon-window refinement metrics. They explain the observed `B_pi^0 -> B_pi^A` redistribution and do not replace direct before/after reweighting plots.

## F.6.3 — parallel full procedure plus Method A

`CLOSED / RUNTIME VALIDATED` only for current-baseline candidate-lineage
adoption at Q4p4W2p74 / Left / lowe. The supplied independently reviewed
[runtime evidence](../evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md)
establishes live-cache parity, real child-yield changes and parent preservation
at farm HEAD `da38444e7aa60efd62d6638780776344daf40276`.
The preceding independent source review covered the private branch changing
exactly `w0_j -> w0_j * C_j` while filling the pion-subtraction template.
No canonical-five closure or production promotion follows.

The baseline full yield calculation remains unchanged. Both branches use the
same proton-cleaned input; random/dummy subtraction, slow-proton treatment,
pion component models/fits/windows/amplitudes except the event-level factor,
canonical binning, SIMC, ordinary current non-empirical yield extraction,
efficiencies, acceptance, L/T separation, and cross-section formulas remain
frozen and identical. The Method-A-clean kaon spectrum is produced directly
after pion subtraction and carried through the existing yield extraction.
Legacy empirical residual Fit 1/Fit 2 remain disabled:
`BG_OPT_ACTIVE_PROFILE = "no_empirical_residual"` stays frozen unless a
separate future scientific contract changes it. F.6.3 must not activate, rerun,
tune, or compare those historical fits. Disabling Method A must reproduce
baseline exactly.

## E.8.4 — baseline-versus-Method-A production-impact audit

`CLOSED / RUNTIME VALIDATED` only for the existing Q4p4W2p74 / Left / lowe
impact audit, from the supplied [runtime evidence](../evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md).
The preceding independent source review passed
`kaonlt_review(20260928-233040).diff`.
E.8.4 consumes authoritative F.6.3 two-branch outputs and does not construct
them. The reviewed candidate retains Fix.1's E.8.2-to-F.6.3 setting-token
identity, internal analysis-MM geometry, and common pion input, and Fix.2's
explicit malformed optional-sidecar Lambda-window/child-inventory fail-closed
boundary. Codex-reported deterministic checks were `NOT RUN by ChatGPT`; this
historical source review alone did not establish ROOT/PyROOT, full-analysis,
farm or runtime validation.

### E.8.4a — actual pion-subtraction consequence

Using the same proton-cleaned input compare `B_pi^0`, `B_pi^A`,
`B_pi^A - B_pi^0`, baseline clean kaon `K_0`, Method-A clean kaon `K_A`, and
`K_A - K_0`.

### E.8.4b — final canonical `(t,phi)` missing-mass comparison

For every populated canonical cell show baseline `MM_0(t,phi)`, Method-A
`MM_A(t,phi)`, overlay, signed difference, and the identical signal/integration
window.

### E.8.4c — final yield comparison

For every canonical `(t,phi)` show `Y0(t,phi)`, `YA(t,phi)`, `DeltaY = YA - Y0`, and `DeltaY/Y0` where defined. Summarize versus phi for each t, versus t for each setting, and later across all five canonical Q4p4W2p74 settings. The Method-A minus baseline shift is a correction effect, not automatically a systematic uncertainty.

### E.8.4d — shareable SIMC and yield-impact presentation

[Fix.5 / Fix.5.1 / Fix.5.2](../phases/e8-4-fix5-shareable-method-a-impact-pages.md)
records the historical `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`
checkpoint. The later [Fix.5.8 closure](../evidence/e8-4-fix5-8-left-lowe-runtime-closure-2026-10-03.md)
accepts only the narrow Left/lowe owner/presentation scope. For each canonical t parent,
show all nine children in final Method-A MM versus authoritative already-normalized
per-child SIMC and baseline/Method-A/SIMC overlays. Clone the existing SIMC
runtime support only after the single existing SIMC-yield producer; never reload,
refill, rescale or use a setting-wide substitute.

Shareable summaries show stored Y0/YA with existing statistical uncertainties,
signed DeltaY and defined DeltaY/Y0 versus phi. Undefined fractions remain
visible; no new DeltaY uncertainty is introduced. The parent-normalization
sanity page displays hash/fingerprint-pinned stored candidate F.4 signed parent
sums and residuals. Equality is the preservation constraint, not the magnitude
of child-yield changes. Existing E.8 pages remain.

The tracked narrow owner runs the existing Left/lowe debug analysis, verifies
fresh pages/artifacts, then collects only Left/lowe using the canonical-five
profile and a temporary clean detached collector/source-check worktree at the
exact reviewed/pushed SHA. Ordinary canonical farm artifacts are consumed
without ordinary-checkout cleanup or a separate manual packaging step.
The accepted Fix.5.8 evidence supplies narrow Left/lowe owner/presentation
validation; other-setting rendering and integration remain unaccepted.

### User-clarified final presentation requirement

User intent, clarified October 11, 2026, retains four named presentation
categories: baseline, Method A, reweighted baseline and K+Lambda SIMC.
The subsequent user correction establishes their distinct scientific roles;
these are not four established independent extracted-yield/MM pairs.

| Requested category | Established object / final-yield meaning | Availability and limits |
| --- | --- | --- |
| Baseline | `Y0(t,phi)` / `MM_0(t,phi)` after baseline pion subtraction with `w0` | Baseline-five numerical/procedure runtime closure remains accepted. |
| Method A | Independently derived relative response `r` and detached response/weight diagnostics; F.4 parent-normalized `C` | Not absolute mis-ID probability or standalone kaon yield. Standalone Method-A-only final yield/MM is NOT VERIFIED. |
| Reweighted baseline | Existing private `YA(t,phi)` / `MM_A(t,phi)` after applying `w0*C` | Already produced by F.6.3; narrow Left/lowe accepted scope only. `Y_rw/MM_rw` were request notation, not extra output keys or a missing producer. |
| K+Lambda SIMC | Existing stored SIMC yield/error and per-cell MM support | Baseline-five numerical-support matching accepted; absolute units remain blocked. |

SOURCE VERIFIED: F.3 derives relative HGCer response `r` from its own
positive-response training on accepted `hgcer3`. F.4 forms signed baseline
`b_j=s_j*w0_j`, `B=sum b_j`, `U=sum b_j*r_j`, and valid-parent
`C_j=r_j/(U/B)`, preserving `sum b_j*C_j=B`. Independent relative-response
derivation does not make normalized `C` independent of baseline `w0` and
signed source coefficients. No independent `(t,phi)` child normalization.
F.6.3 fills `B_pi_A` using `s_j*w0_j*C_j`, forms
`MM_A=pion_input-B_pi_A`, and extracts `YA` with the existing uncertainty
helper. `YA/MM_A` is the final Method-A-reweighted baseline in source/API
naming, not a separate Method-A-only final kaon result. `B_pi_A` is the
pion-background component, never an extracted final signal yield.

The [renderer](../../../src/cuts/full_background_subtraction_plots.py) still
implements baseline/private-reweighted-baseline/SIMC three-way MM overlays
and stored `Y0/YA` plus `DeltaY` summaries. Its existing Method-A labels and
API names are unchanged by this memory correction. E.8.4c/d retain their
implemented two-/three-series scope. Prior validated pages are neither
relabeled four-way nor failed by this intent; no new page or runtime claim.

Only genuine comparable final products may appear as final-yield/MM series,
with common canonical geometry and Lambda window and existing uncertainties.
Method-A relative response/weight diagnostics may occupy their own diagnostic
category without false yield units. Four independent final curves cannot be
promised until genuinely distinct final objects are defined; no unrelated
new reweighted-baseline formula is required. Retain unavailable/undefined
states, exact impact references and no fabricated `DeltaY` errors or SIMC
scaling. Absolute-SIMC units remain blocked; no amplitude agreement is
authorized. This correction creates no phase, correction, acceptance change
or farm gate. CURRENT owns the sole substantive NEXT; microscopic F.4
investigation and F.1 archaeology remain `DEFERRED` and the exact owner gate
remains binding.

### Live four-way producer audit

SOURCE VERIFIED against live `bfb0ae21ee0645f1c30cbec877a74979350f4548` under
[the implementation contract](../phases/e8-four-way-final-yield-mm-implementation-task-contract.md).
This historical bounded source audit is not new runtime acceptance. Its literal
key inventory was valid, but the inference from naming was wrong: `YA/MM_A`
already owns the final reweighted baseline. The role correction is traced
against live `78ecadbacf58ab25b6a3670137633b3124074fe7`, without changing
source, historical runtime identities or accepted scopes. Producer paths
below are in [calculate_yield.py](../../../src/binning/calculate_yield.py);
the [E.8.4 consumer](../../../src/cuts/full_background_subtraction_plots.py)
remains unchanged. No four-populated-series result or new page is implemented.

| Requested pair | Exact producer and storage/object authority | Final stage and existing errors | Units and accepted scope / missing dependency |
| --- | --- | --- | --- |
| `Y0/MM_0` | `calculate_yield_data()` stores `e8_2_yield_measurements[(j,k)]`: `yield`, `statistical_error`, `total_error`, `final_histogram`. `find_yield_data()` returns `binned_DATA[phi_setting]["yield"][(j,k)]["yield"/"yield_err"]`. `_build_f6_3_parallel_method_a_source()` copies the measurement to `children[].Y0` and clones `H_MM_after_pion_subtraction` to `MM_0`; its identity audit checks the actual measured final histogram. | Final baseline after pion subtraction; no active empirical-residual stage. `_measure_yield_with_current_uncertainty()` supplies statistical and total errors, copied as `Y0_statistical_error/Y0_total_error`. Wider E.8.2 `after_pion_final` stage snapshots/window diagnostics are not replacement extracted-yield authority. | Effective-data-charge normalized, with existing dummy target/charge normalization. Baseline-five numerical/procedure closure accepted at `8e68d2acae91fd97622793df563aba0296403a48`; all five settings. E.8.4's two-branch consumer still requires an available private source, not baseline acceptance alone. |
| `YA/MM_A` | `_build_f6_3_parallel_method_a_source()` returns the source attached by `calculate_yield_data()` at `hist["_f6_3_parallel_method_a_source"]`, schema `f6_3_parallel_method_a_source/v1`, branch `parallel_nonproduction_method_a_full_analysis`; `children[].YA/MM_A`. | Same proton-cleaned `pion_input` minus producer-owned `B_pi_A` filled with accepted `w0*C`; no child renormalization. Same measurement helper supplies `YA_statistical_error/YA_total_error`. `B_pi_A` itself is background, not another final clean spectrum. | Same data normalization as baseline. Narrow Left/lowe adoption/impact closure at `da38444e7aa60efd62d6638780776344daf40276` remains accepted. Newer private canonical-five owner is `BLOCKED` before analysis by exact F.4 equivalence; no five-setting availability follows. |
| Reweighted-baseline request `Y_rw/MM_rw` | Existing `children[].YA/MM_A` in `_build_f6_3_parallel_method_a_source()` already owns this role; request symbols are not extra keys. The preceding `YA/MM_A` row describes the same final branch, not an independent fourth pair. | Same `pion_input-B_pi_A`, with `B_pi_A` filled using signed `w0*C`, and existing stored `YA` errors. A standalone Method-A-only final pair is NOT VERIFIED; detached response/weights are diagnostics. | Same normalization and narrow accepted Left/lowe scope as `YA/MM_A`. No missing reweighted-baseline producer. Private canonical-five owner remains `BLOCKED` by exact F.4 equivalence; no four-independent-series claim. |
| `Y_SIMC/MM_SIMC` | `process_hist_simc()` fills cut-selected per-cell `H_MM_SIMC` from `h10.missmass` with `iter_weight`, scales once by `normfac_simc`, and detaches `support_hist_dict["mm"][j][k]`. `calculate_yield_simc()` attaches `hist["_xsect_support_simc"]["mm"]` and returns stored `yield/yield_err` via `find_yield_simc()` at `binned_SIMC[phi_setting]["yield"][(j,k)]`. E.8.4 currently clones the MM support under `histograms.SIMC`, not the stored yield. | Existing producer integrates normal bins with statistical Sumw2 error; negative-yield or nonfinite yield/error cases retain its existing zero convention. No separate total-error field. `normalization_audit.children[].SIMC_integral_in_the_same_lambda_window` is diagnostic, not the stored clamped yield/error authority. | `iter_weight * normfac / Ncontribute`; `get_eff_charge.find_events()` reads SIMC normfac/Ncontribute. `src/simc_ana/iter_weight.py:iter_weight()` changes SIMC model event weights, not data baseline yields. Five-setting numerical-support matching accepted by the baseline closure; absolute data/SIMC units remain `BLOCKED` by `SIMC_normfac_luminosity_and_charge_units_not_source_proven`. Stored scalar/error handoff to E.8.4 is not presently implemented. |

Geometry/error ownership: internal `(j,k) = (t_index,phi_index)` is zero-based;
producer processed keys `t_bin{j+1}phi_bin{k+1}` and printed bin numbers are
one-based. Canonical-five inventory is Left/lowe, Left/highe, Center/lowe,
Center/highe, Right/highe, each 3 t parents x 9 phi children. Final numerical
MM objects have 100 normal bins in `inpDict["mm_min"]` to `inpDict["mm_max"]`
(accepted baseline-five window 1.10–1.16 GeV). Wider stage displays remain
separate. `build_full_background_subtraction_e8_4_payload()` validates setting,
branch/authority/flags, complete coordinates, common Lambda window, identical
analysis-MM edges, finite contents/errors and stored baseline scalar identity;
it clones display objects and never supplies a fourth producer.

The data statistical error is the square root of summed normal-bin error
squares; total error additionally retains existing data-charge, dummy-charge
and residual-fit terms from the shared producer helper. Active residual scales
remain zero. E.8.4 yield summaries copy stored statistical errors, retain signed
`DeltaY = YA-Y0`, and leave `DeltaY/Y0` undefined when `Y0=0`; no new
DeltaY uncertainty. SIMC absolute-comparison-unavailable remains effective.
Existing page IDs, schemas, inventories, science and narrow closures are intact.

**Historical source stop and corrected interpretation:** The bounded audit
correctly found no separate `Y_rw/MM_rw` keys and no fourth independent final
pair. Its conclusion that an independent reweighted-baseline producer was
missing was incorrect: F.6.3 already applies Method A to baseline and produces
`YA/MM_A`. No standalone Method-A-only final kaon yield/MM has been established.
No additive page or source change was made by that audit or this correction.
The user's intent is retained by distinguishing weight/response diagnostics
from final-yield products, not by manufacturing another curve. CURRENT owns
the next substantive detached Method-A weight/response science review; this
memory correction is not physics completion, farm readiness or gate acceptance.

## Final E.8 closure and F.6.4

Final E.8 is `BLOCKED` by the newer October 8 canonical-five owner failure
before analysis at
`f6_3_scientific_equivalence_mismatch:f4:$.parents[0].absolute_baseline_parent_sum`
at `ebf3a0be2f50c788c22defe3601fd5b7fee48644`. F.2/F.3 exact PASS and
F.4 exact FAIL are diagnostic evidence, not a roundoff exemption or a proven
physics defect. The historical second full run failed page/artifact provenance
after analysis; its reported equality does not accept this newer mismatch.
Separate lineage preflight remains accepted pre-analysis only
(`analysis_started=false`). The owner repairs retain
`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; full-owner retry requires
a separately reviewed admissible resolution of the exact invariant.
See the [supplied forensic receipt](../evidence/e8-4-f1-forensic-reproducibility-diagnostic-review-2026-10-08.md).
Full-procedure extracted-yield readiness has priority; further F.1 archaeology
and further stage-yield/broader physics interpretation are `DEFERRED` under
the user decision until material gate evidence/proper yields respectively.
CURRENT owns the sole substantive next action. Accepted narrow Left/lowe
E.8.2/E.8.3 and F.6.3/E.8.4 scopes remain intact. Closure requires prompt/random, dummy, slow proton, baseline pion,
final baseline canonical MM/yields and stage yields, explicit `w0 -> w0*C`,
direct baseline-versus-reweighted pion background with signed difference/ratio,
canonical redistribution/parent closure, F.6.2 explanation, actual F.6.3
baseline-versus-Method-A clean-kaon MM comparison, final `Y0` versus `YA`, and
absolute/fractional yield shifts, per-child Method-A/SIMC and baseline/Method-A/SIMC
comparisons without renormalization, shareable yield-impact summaries and parent
closure sanity presentation. Accepted existing Left/lowe evidence does not
close the canonical-five program. Absolute-SIMC interpretation remains
separately `BLOCKED` by
`SIMC_normfac_luminosity_and_charge_units_not_source_proven`.

F.6.4 is `BLOCKED` pending completed F.6.3/E.8.4 production-impact evidence. It is the only phase that may decide whether Method A becomes production pion-background treatment. No automatic promotion follows from detached validation, presentation, or source review.

## Farm milestone cadence

Memory reconciliation is local/documentation only. E.8.2, E.8.3, and coherent F.6.3/E.8.4 development may proceed through deterministic local checks and source review without a Jefferson Lab farm run after every update. Farm validation remains mandatory for ROOT/PyROOT/full-runtime claims. After a coherent F.6.3 plus E.8.4 branch is source-reviewed, run one narrow `Q4p4W2p74 / Left / lowe` end-to-end farm gate, inspect fresh artifacts, then PASS or make one coherent repair; broaden only after that narrow gate passes.
