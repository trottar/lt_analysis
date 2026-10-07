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

`SOURCE REVIEWED` — exposes the complete existing baseline as a coherent
before/component/after spectrum and yield chain without changing physics.

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
is `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`. For each canonical t parent,
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
New-page rendering and owner integration still require fresh farm evidence.

## Final E.8 closure and F.6.4

Final E.8 is `BLOCKED` by the second canonical-five full runtime gate's
artifact/page provenance failure; provenance/identity repair remains `DEFERRED`
unless explicitly reopened by user decision. Accepted narrow Left/lowe
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
