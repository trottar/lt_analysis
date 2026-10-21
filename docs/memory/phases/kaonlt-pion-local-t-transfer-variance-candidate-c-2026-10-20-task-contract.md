# KaonLT — Detached candidate C: local-t transfer-variance regularization

**Task class:** Narrow scientific diagnostic implementation, not production correction or a new phase.  
**Source-audited remote:** `test` at `52d8b44335efde614111f727d90fd7987f86d941` (audit anchor; always recheck live state).  
**Farm input evidence:** Accepted October 19–20 A/B four-artifact bundle, source `6d13abb4681b938c56def2a1f455244b0880096c`.  
**Scope:** Q4p4W2p74, Center/lowe t1–t3 and Left/lowe t1–t3, separately.  
**Work state:** `ACTIVE` — first detached A/B/C numerical feasibility comparison only; no promotion, farm authorization, or fitting-contract alteration in production.

## 1. Scientific objective and nomenclature

Introduce a third, **diagnostic-only** comparator:

- **A:** The authoritative canonical-t parent fit using inherited/setting-wide pion templates.
- **B:** The accepted genuinely local-t template fit with the existing fitter.
- **C:** **Local-t transfer-variance-regularized**, with the *identical aligned component templates and frozen kaon-side pion model as B*, but a soft variance-aware refinement of the pion-control amplitudes.

These are **template-comparison A/B/C**, not HGCer Method A/B and not the fitted `A_n`/`B_n` amplitudes. Neither physical component shapes nor protected Lambda/Sigma0 logic is altered. C is a test of estimator conditioning, not an approved cure.

The scientific question is whether variance-aware fitting can reduce excessive signed-event transfer variance **without** destroying measured MM structure or biasing the estimated pion background beneath Lambda. No constraint may target a chosen Lambda yield, data/SIMC ratio, or preconceived background fraction. C may fail; retain such results as explicit evidence.

## 2. Source/worktree and hard-stop gate

Before any edit: establish actual live `test` and `origin/test` identity, complete local tracked/staged/unstaged/untracked worktree and pending task contracts; preserve unrelated work exactly. Read repository-root `AGENTS.md`, then `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, `docs/memory/USER.md` in exact order. Then read `docs/memory/MAINTENANCE.md`, `docs/memory/CODEX.md`, the current pion investigation, the October 19 A/B comparison contracts, and the specific source cited below. Emit the full tracked health check.

Treat this contract as a newly proposed diagnostic task under the existing temporary baseline-pion interruption. The present `CURRENT.md` NEXT predates the user-directed return to pion modeling; reconcile it narrowly, preserving closed phases, other blockers, and the sole NEXT. Stop `BLOCKED` if the worktree overlaps, an intended new path exists, the farm JSON provenance differs, or any proposed operation requires a production-source change.

## 3. Frozen boundaries and existing evidence

Do **not** modify `src/`, launchers, fitting/weight/application code, `src/cuts/pion_t_template_comparison.py`, owner/collector/profile, input data/SIMC, fits, priors, alignment, templates, normalization, cuts, histograms, E.8/F.6/Method A/B interfaces, yields or cross sections. No farm run, ROOT/PyROOT claim, commit, or push by Codex.

Use only the two accepted comparison JSON objects from the October 19–20 A/B ZIP for initial C: `global/Center_kaon_Q4p4W2p74_lowe_pion-t-template-comparison.json` and `global/Left_kaon_Q4p4W2p74_lowe_pion-t-template-comparison.json`. Verify ZIP/manifest identity if supplied as an archive, source HEAD `6d13abb4681b938c56def2a1f455244b0880096c`, six available parents, frozen input content/error fingerprints, fixed alignment, valid genuinely-local B templates and all recorded fit-status limitations. Do not infer new farm or scientific acceptance from source availability.

Current source provides **all ingredients for a detached histogram-level feasibility test**: `pion_background_chain_audit.histogram()` retains normal/flow bins and explicit `stored_sumw2`; the accepted B snapshot contains `H_simc_shape_pi_n`, `H_simc_shape_pi_delta`, `H_simc_shape_pi_sidis` **after the fit's alignment**, `H_kaon_pion_bg_fit_total_refined`, fitted B coefficients, staged prior centers/regularization widths, and `diagnostics.pion.joint_refinement.fit_bin_indices`. Do not substitute the earlier *raw* `templates[*].normalized_values` for aligned fitting columns. `stored_sumw2` is a statistical variance proxy for the signed histogram, not event-level full covariance.

No assumption that six-bin A/B JSON contains independent training/validation partitions or exact all-cuts event-level source identities. Any in-sample C sweep is **feasibility only**. Do not select a production regularization strength without independent validation.

## 4. Controlled mathematical test

Reproduce B's pion joint fit **first**, for every Center/Left t1/t2/t3 parent. Use the identical valid fit-bin list, data vector and `stored_sumw2`, already aligned B SIMC-template columns, nonnegative coefficient bounds, original stage amplitudes and original `joint_refinement.regularization_widths`. Require exact modeled B data chi-squared and regularization contribution within deterministic numerical tolerances, and all coefficients within `1e-7` absolute or tighter of the supplied B solution. Stop on any failure; C must reduce to B when the new penalty strength is zero.

Define the prospective *fixed* kaon-side model `K_B(i)` from the **protected-applied** B kaon-pion background, not from an unprotected staged total or signal-containing data. Do not refit `K_B` as part of C: otherwise the variance penalty could artificially shrink the kaon numerator to zero. Only the three pion-control amplitudes `beta_j >= 0` vary. Use B's aligned templates `T_ij`, `P_beta(i)=sum_j beta_j*T_ij`, and the existing joint data plus staged-prior least-squares objective.

Define the statistical transfer-risk diagnostic on bin support `D` by

`V_D(beta) = sum_{i in D} S2_i * [K_B(i)/P_beta(i)]^2`,

where `S2_i` is the explicit signed-control histogram Sumw2. It incorporates squared prompt/random/dummy normalization coefficients as retained by the histogram producer; it is not the variance of plotted weights across MM bins, and not the complete fitted-coefficient covariance. Treat `K_B(i)>0` and `S2_i>0` with `P_beta(i)<=0` or nonfinite values as **unsupported**, not as free zero-cost weights. Never exploit the production floor's `weight=0` behavior to minimize the penalty; do not introduce a post-hoc hard clipping threshold or modify fitted/model data to evade unsupported support.

Use a *predeclared* normalized soft penalty of the above transfer risk in a **detached** objective, in addition to (not replacing) B's data-fit and existing staged-prior objective. Its mathematical normalization, numerical support policy and lambda scan must be printed into the diagnostic provenance before comparing candidates. Scan `lambda=0` and multiple log-spaced positive values, with no outcome-driven tuning of the scanned grid. The normal-bin full MM domain and the accepted Lambda yield window must both be reported separately; do not secretly exclude the 0.94–1.03 GeV region because it is numerically troublesome, or use it to mislabel a signal-window yield failure.

A candidate with unsupported positive kaon model/control variance must be flagged `unavailable` or `unstable`, with exact offending MM bins. A candidate that resolves transfer variance by destroying pion-control shape must be identified as such, not accepted. Quantify fit-parameter sensitivity and bound hits; covariance-derived estimates are conditional and are not automatically valid confidence intervals near bounds.

The fitting algorithm may be a deterministic standalone SciPy numerical optimizer acting on the recorded arrays; it must **not** change the authoritative ROOT fitter or silently call it with new settings. Use exact B baseline reconstruction at lambda zero as the numerical control. If solver nonconvergence, inconsistent templates or absent Sumw2 prevent a comparison, fail closed and report the first invariant.

## 5. A/B/C comparative evidence for all six bins

Produce a single **detached local** JSON plus a concise, readable **detached local** PDF (or explicitly report rendering unavailable). Neither is an accepted farm product or a replacement for source-tracked analysis procedure pages. Keep the two setting populations separate, with all three canonical-t bins in each. Record source ZIP/JSON hashes, candidate formulas, solver/versions, all parameter values and complete failed attempts.

Per bin, retain and display:

1. A and B **unchanged**, plus C for **each scanned lambda**; clearly distinguish any illustrated, unselected C from a validated choice.
2. Pion-control signed MM vs model totals, component contributions and uncertainty-aware residuals, **full spectrum and Lambda-region zoom**, with common comparison scales.
3. The kaon-side pion model `K_B` held fixed for B/C; measured post-proton kaon spectrum; A/B and prospective C transfer-applied backgrounds, showing all-cuts Lambda window separately from extended no-MM-cut diagnostics.
4. Model-applied/background-integral closure, maximum and tail transfer weights, `V_D`, the in-window component/statistical uncertainty proxy, maximum binwise amplification, unsupported support bins, staged/joint coefficients and correlations, and fit chi-squared/p-values with correct labels.
5. Full-range `0.7–1.3 GeV` and signal `1.10–1.16 GeV` risk summaries plus the narrower `1.105–1.135 GeV` visual comparison. Provide a Pareto-style sweep (MM shape fit quality versus transfer variance), not a claimed optimal lambda without held-out evidence.
6. Explicit evidence limits: no event-by-event independent application replay or detailed correlated covariance from the A/B JSON alone; historical Center A replay is independent prior evidence, not new farm verification. Do not claim full RF-restored yields, actual C event application, physical pion purity, or corrected cross sections.

Reject silent relabeling of diagnostic candidate C as production or as HGCer Method C.

## 6. Narrow edit scope and tests

New files (verify absent):

- `testing/pion_local_t_transfer_variance_comparison.py` — detached source-agnostic A/B/C array reader, optimizer, verified comparison output and optional PDF generation. No import-time ROOT or farm file operations.
- `testing/test_pion_local_t_transfer_variance_comparison.py` — synthetic deterministic tests and exact A/B/C invariant tests.
- `docs/memory/phases/kaonlt-pion-local-t-transfer-variance-candidate-c-2026-10-20-task-contract.md` — this contract after user transfer.

Existing files, only if required:

- `docs/memory/CURRENT.md` — short scoped work-state and sole NEXT correction, no reopening E.8/F.6 closure.
- `docs/memory/investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md` — attributable variance-aware candidate rationale and audit status; no invented runtime acceptance.
- `docs/memory/manifest.json` — generated integrity update.

No other path is authorized. In particular do not edit `src/cuts/pion_component_fits.py`, `src/cuts/pion_component_subtraction.py`, `src/cuts/rand_sub.py`, A/B diagnostic producer, or execution/collector owner. A subsequent contract may propose farm-integrated C **only after** the detached six-bin C feasibility evidence is independently reviewed and found scientifically promising.

Minimum tests: all six existing B fits match at lambda zero when accepted JSON is available; analytic synthetic positive/support-absent cases; zero numerator/denominator; independently signed sources with nonzero Sumw2 despite signed cancellation; variance-scaled event weights; perturbation/monotonic risk sensitivity; solver nonconvergence; correlation/uncertainty limits; no post-fit clipping; isolated provenance rejection; A/B content not mutated; output schema, numeric finiteness, reproducible JSON hashing; exact input/path scope. No local farm or ROOT behavior claims. Run relevant deterministic unit tests, compilation, diff/whitespace and memory manifest/health checks. If a suitable local SciPy/NumPy interpreter is not available, report the precise test as unexecuted; do not install packages or claim fake success.

## 7. Review and next gate

Codex performs only the allowlisted local implementation and tests, preserving all existing user files/index/HEAD/refs. Provide the **complete actual diff**, including complete new files, for independent ChatGPT review. Do **not** commit/push, run the farm, install ROOT, or use a farm path. A missing required artifact is an explicit evidence limitation, not a reason to manufacture a numerical result.

After independent actual-diff approval the **user** publishes; ChatGPT verifies pushed blobs/HEAD and NEXT. The next scientific gate is to run/review the **six-bin detached A/B/C variance-versus-shape comparison** on the accepted A/B JSON, with all λ values and meaningful MM plots. A farm-integrated C candidate or weight-rule change is a *later* separately reviewed decision, not authorized by this task.
