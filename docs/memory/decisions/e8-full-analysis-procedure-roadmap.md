# E.8 full-analysis procedure roadmap

## Governing decision

E.8 is the complete visual audit of how the kaon missing-mass spectrum and signal-region yield evolve through the analysis. Every substantive stage exposes the authoritative input spectrum, component or treatment applied, and authoritative output spectrum. The presentation culminates in final canonical `(t,phi)` missing-mass spectra and extracted yields.

E.8 is presentation-only. It consumes authoritative upstream runtime objects, persisted snapshots, accepted detached Method-A artifacts, or later F.6.3 branch outputs. It never recomputes a fit, factor, correction, normalization, or yield merely for plotting.

The baseline production branch remains authoritative and unchanged. Method B remains diagnostic/cross-check only and is numerically excluded from every Method-A correction/application. Method A uses exactly the accepted F.4 parent-preserving correction; no child `(t,phi)` is independently renormalized. The accepted F.6.2 JSON SHA-256 and fingerprints remain frozen.

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

`NEXT` — expose the complete existing baseline as a coherent before/component/after spectrum and yield chain without changing physics.

### E.8.2a — prompt/random subtraction

Show prompt kaon missing-mass input, random contribution removed, after-random missing-mass output, and a difference or overlay where useful.

### E.8.2b — dummy subtraction

Show after-random input, normalized dummy contribution, after-dummy output, and a difference or overlay where useful.

### E.8.2c — slow-proton subtraction

Show after-random/dummy kaon input, proton contamination estimate, proton-cleaned kaon output, before/after comparison, and existing PID/weight diagnostics. `_process_yield_data_tree()` already applies accepted proton cleaning during initial `(t,phi)` filling. Existing per-cell raw/random/dummy snapshots must not be relabeled as pre-proton. E.8.2 must capture or consume a true authoritative pre-proton snapshot rather than reconstructing or guessing one.

### E.8.2d — baseline pion subtraction

Show proton-cleaned kaon input, baseline pion background actually removed, baseline pion-subtracted kaon output, before/after comparison, and difference. The baseline contribution remains `b_j = s_j * w0_j`.

### E.8.2e — residual background Fit 1

Show post-pion input missing mass, the actual production Fit-1 background component removed, post-Fit-1 output, before/after comparison, and existing fit/window information. Do not refit for presentation.

### E.8.2f — residual background Fit 2

Show Fit-1 output entering Fit 2, the actual production Fit-2 background component removed, final background-subtracted missing mass, and before/after comparison. Do not refit for presentation.

### E.8.2g — final baseline canonical `(t,phi)` missing mass

For every canonical t parent show all nine canonical phi children, including explicit empty/invalid status. Every populated final panel exposes canonical-t identity, phi index/range, the final missing-mass spectrum actually used for yield extraction, Lambda signal/integration window, relevant final fitted/background information, extracted yield, and statistical uncertainty.

### E.8.2h — baseline stage-yield audit

Expose the existing authoritative per-`(t,phi)` progression: prompt, after random, after dummy, after proton, after pion, after Fit 1, and after Fit 2/final. Summarize final baseline `Y0(t,phi)` versus phi for each canonical t. A renderer does not recalculate a yield when its producer can persist or provide it.

## E.8.3 — detached Method-A reweighting audit

`BLOCKED` pending E.8.2 source-reviewed completion. This remains detached/non-production and consumes accepted F.4/F.5/F.6.1/F.6.2 results only.

### E.8.3a — reweighting operation

Explicitly show `b_j^0 = s_j * w0_j` to `b_j^A = s_j * w0_j * C_j`. Method A changes only the accepted pion event contribution by the accepted F.4 parent-preserving factor. No child-by-child renormalization is allowed.

### E.8.3b — baseline versus reweighted pion background

For every canonical t parent show baseline `B_pi^0(MM)`, Method-A-reweighted `B_pi^A(MM)`, overlay, signed difference `Delta B_pi = B_pi^A - B_pi^0`, and ratio `R_pi = B_pi^A / B_pi^0` only where numerically well-defined and presentation-safe. This comparison is a first-class E.8 result, not a legend implication.

### E.8.3c — canonical redistribution and parent closure

For every canonical t show all nine `B_pi,tphi^0`, `B_pi,tphi^A`, and `Delta B_pi,tphi` children. Display stored/accepted closure and relative residual for `sum_phi B_pi,tphi^A ~= sum_phi B_pi,tphi^0`. Never independently normalize a child.

### E.8.3d — F.6.2 acceptance/HGCer explanation

Retain accepted explanatory diagnostics: L/B/A normalized shapes; `delta x xptar`; `delta x yptar`; missing-mass x acceptance maps; support/OOD; effective statistics; and kaon-window refinement metrics. They explain the observed `B_pi^0 -> B_pi^A` redistribution and do not replace direct before/after reweighting plots.

## F.6.3 — parallel full procedure plus Method A

`BLOCKED` pending source-reviewed E.8.2 and E.8.3, not final E.8 closure. F.6.3 is the only step here that constructs the actual parallel Method-A full-analysis branch. It changes exactly `w0_j -> w0_j * C_j` while filling the pion-subtraction template.

The baseline full yield calculation remains unchanged. Random/dummy subtraction, slow-proton treatment, pion components/fits/windows/amplitudes except the event-level factor, canonical binning, SIMC, Fit-1/Fit-2 algorithms/configuration, yield extraction, efficiencies, acceptance, L/T separation, and cross-section formulas remain frozen and identical between branches. Downstream fits rerun normally on the altered Method-A post-pion spectrum with algorithms/configuration frozen. Disabling Method A must reproduce baseline exactly.

## E.8.4 — baseline-versus-Method-A production-impact audit

`BLOCKED` pending F.6.3. It consumes authoritative F.6.3 two-branch outputs and must not construct them.

### E.8.4a — actual pion-subtraction consequence

Using the same proton-cleaned input compare `B_pi^0`, `B_pi^A`, `B_pi^A - B_pi^0`, baseline post-pion kaon `K_pi^0`, Method-A post-pion kaon `K_pi^A`, and `K_pi^A - K_pi^0`.

### E.8.4b — Fit-1 consequence

Compare baseline and Method-A Fit-1 input, fitted background component, output, and branch difference.

### E.8.4c — Fit-2 consequence

Compare baseline and Method-A Fit-2 input, fitted background component, final Fit-2 output, and branch difference.

### E.8.4d — final canonical `(t,phi)` missing-mass comparison

For every populated canonical cell show baseline final missing mass, Method-A final missing mass, overlay, signed difference, and identical signal/integration window.

### E.8.4e — final yield comparison

For every canonical `(t,phi)` show `Y0(t,phi)`, `YA(t,phi)`, `DeltaY = YA - Y0`, and `DeltaY/Y0` where defined. Summarize versus phi for each t, versus t for each setting, and later across all five canonical Q4p4W2p74 settings. The Method-A minus baseline shift is a correction effect, not automatically a systematic uncertainty.

## Final E.8 closure and F.6.4

Final E.8 is `BLOCKED` pending E.8.4 and later required runtime/visual validation. Closure requires the full visual chain: prompt, random, dummy, slow proton, baseline pion, Fit 1, Fit 2, final canonical missing mass/yields, stage yields, detached `w0 -> w0*C`, direct pion-background comparison/difference/ratio, canonical redistribution/parent closure, F.6.2 explanation, actual F.6.3 post-pion/Fit-1/Fit-2/final-MM/final-yield comparisons, and absolute/fractional yield shifts.

F.6.4 is `BLOCKED` pending completed F.6.3/E.8.4 production-impact evidence. It is the only phase that may decide whether Method A becomes production pion-background treatment. No automatic promotion follows from detached validation, presentation, or source review.

## Farm milestone cadence

Memory reconciliation is local/documentation only. E.8.2, E.8.3, and coherent F.6.3/E.8.4 development may proceed through deterministic local checks and source review without a Jefferson Lab farm run after every update. Farm validation remains mandatory for ROOT/PyROOT/full-runtime claims. After a coherent F.6.3 plus E.8.4 branch is source-reviewed, run one narrow `Q4p4W2p74 / Left / lowe` end-to-end farm gate, inspect fresh artifacts, then PASS or make one coherent repair; broaden only after that narrow gate passes.
