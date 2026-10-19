# KaonLT — Controlled canonical-t pion-template comparison

**Status:** `ACTIVE` — diagnostic-only implementation contract  
**Scientific setting:** Q4p4W2p74, Center/lowe, canonical t1/t2/t3; Left/lowe control  
**Source-audited starting HEAD:** `4007d8a43a604628e5b74aca54229d3e6bd48848`

## 1. Objective

Isolate the numerical effect of using genuine canonical-t SIMC pion-component shapes instead of the setting-wide shapes currently inherited by the authoritative t fits.

This is a controlled scientific comparison, **not** authorization to correct, replace, or promote the production background model.

The source audit established that:

- The setting-wide and canonical-t fits call the same full component-fitting builder.
- The normal SIMC loading path does not construct t-resolved pion shapes.
- The t-shape resolver can fall back to the setting-wide shapes.
- Current joint refinement regularizes toward the local staged solution, not the global fitted solution.
- The authoritative (t,phi) path inherits the t-parent fit without a new local component fit.

The user reports that the setting-wide model appears substantially more sensible than the later t-resolved models. The precise cause remains unverified.

**One variable is under examination in this contract: pion SIMC template shape.** Hierarchical amplitude constraints and (t,phi) refitting are later tasks.

## 2. Authority and startup

Before editing:

1. Establish live `origin/test`, local branch/HEAD, complete staged/unstaged/untracked state, relevant pending patches, and repository identity. Preserve all unrelated bytes and index state.
2. Read root `AGENTS.md`, then `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, and `docs/memory/USER.md`, in that order.
3. Read the applicable `MAINTENANCE.md`, `CODEX.md`, `COMMUNICATION.md`, `TOOLS.md`, October 19 hierarchy memory checkpoint, and baseline-pion investigation.
4. Recheck the live production path before implementing. If a material assumption is wrong or a task path conflicts with existing work, stop `BLOCKED` and report the first invariant.

The starting SHA identifies the audited source, not permission to overwrite a newer checkout.

## 3. Frozen scientific boundaries

Do not alter:

- Existing setting-wide or authoritative canonical-t fit results.
- Staged or joint fit mathematics, configured priors, postfit scales, windows, or validation thresholds.
- Random/dummy subtraction, slow-proton cleaning, physical pion-control selection, or signed-event normalization.
- SIMC source files, event weights, cuts, or production normalization.
- Canonical t/phi bin boundaries and missing-mass coordinates.
- Protected Lambda/Sigma0 procedures.
- The existing pion-weight transfer, child application, yield, or cross-section chain.
- Method A/B ownership or accepted roadmap status.

The comparison must not change the existing main PDF, `pion-fit-debug` PDF, authoritative parent outputs, or production subtraction.

## 4. Required comparison

For each Center t1/t2/t3, with Left t1/t2/t3 as a separate control:

**Baseline A:** Existing authoritative t-parent fit with its currently selected pion templates.

**Candidate B:** Repeat `build_particle_subtraction_component_result()` on the same frozen t-parent pion-control and post-proton kaon-side histograms, replacing only the three pion component templates with correctly constructed t-local SIMC shapes.

Both fits must use:

- Identical pion-control histogram contents and Sumw2.
- Identical kaon-side histogram contents and Sumw2.
- Identical component-fitting configuration and numerical solver.
- Identical protected-signal templates.
- Identical coordinate convention, alignment parameters and fitting windows.

For this first comparison, hold the existing t-parent alignment fixed. Do not rescan or recalibrate alignment as a second uncontrolled change. Record that candidate B is a *fixed-alignment template-shape test*, not a completed local optimization.

Retain baseline A by reference, never overwrite it.

## 5. Correct t-local SIMC construction

The current SIMC loader independently normalizes each populated `(t,phi)` shape, while `_sum_hist_list_to_unit_area()` combines those already normalized children.

That procedure generally constructs an equal-phi mixture rather than an event-weighted t-integrated distribution.

**Do not use that mixture as candidate B.**

Instead, accumulate original selected SIMC contributions directly into each canonical t interval, preserving the original event weights and histogram Sumw2. Normalize each completed t-integrated component histogram **once**.

Use the existing SIMC selection and coordinate conventions. Do not substitute the detached support-audit's raw event counts for the production SIMC weights without establishing their equivalence.

For each component, record:

- ROOT source identity and selected tree.
- Canonical t boundaries, coordinate identity and phi coverage.
- Selected event count, signed/absolute weighted sums, Sumw2 and effective support.
- Raw and normalized histogram integrals.
- Underflow/overflow and events excluded by canonical acceptance.
- Difference from the global template and baseline-selected template.
- Whether a candidate shape is genuinely local, globally inherited, empty, or unavailable.

Prove weighted child-to-parent and canonical-selection closure before calling a candidate shape valid.

A sparse or invalid candidate must be recorded as unavailable. Do not silently substitute the global template while labeling it t-local.

## 6. Diagnostic integration

Prefer an opt-in, additive diagnostic beside the existing cuts-stage authoritative t-parent construction, where the original parent inputs and fits already exist.

The implementation must be disabled by default. Use one explicit diagnostic enable mechanism, verified against the existing launcher; do not introduce a change to ordinary analysis selection.

Reuse the existing component fitter rather than writing a second fitter.

Any new template-building helper must preserve the current loader's default return values, cache behavior and production consumers.

If the candidate cannot be constructed or fitted without modifying authoritative inputs, mutating frozen parents, or writing unintended alignment artifacts, stop and report the blocker rather than weakening these boundaries.

## 7. Required evidence

Produce one detached, source-identifiable comparison JSON and one concise detached diagnostic PDF for the selected settings.

The JSON must include:

- Diagnostic source HEAD and relevant source-file identities.
- Original analysis input identity where available.
- Setting, epsilon, canonical interval pair, t index and boundaries.
- Parent fit identity and frozen input histogram fingerprints.
- Full component template provenance and normalization.
- Baseline A and candidate B staged/refined component coefficients and uncertainties.
- Fitting modes, status, chi2/ndf, fit p-value, validation metrics, component correlation and bound-hit diagnostics.
- Component/template availability and explicit rejection reasons.
- A before/after comparison of model construction **without applying pion subtraction**.
- A check that authoritative parent identity and fitted results are unchanged.

The PDF must show, for each t bin, the same pion-control and kaon-side data with baseline A and candidate B staged/refined predictions. It must clearly distinguish data, component shapes, total fitted model, and residuals.

Do not label either model scientifically accepted merely because its optimizer or existing software gate returned success.

## 8. Minimum deterministic checks

Add focused tests for:

1. Correct raw weighted t integration versus the incorrect equal-phi average, including deliberately unequal child populations.
2. Exactly one normalization after t integration, with preserved weighted bookkeeping.
3. Canonical boundary handling, acceptance coverage, empty bins and missing inputs.
4. Baseline/candidate input identity and fixed-alignment enforcement.
5. Diagnostic disabled by default and no changes to authoritative objects or existing output routing.
6. Explicit failure when a candidate's provenance or support is invalid.
7. JSON serialization and deterministic comparison summaries.

Run applicable focused tests, compilation, whitespace/diff checks, tracked memory-manifest checks, and ordinary memory health. Missing local ROOT/PyROOT must be reported as unexecuted farm-only validation, not replaced by a claimed full integration test.

## 9. Allowed edit scope

Existing files, only where required:

- `src/cuts/rand_sub.py` — narrow opt-in diagnostic attachment.
- `src/utility/pion_component_shapes.py` — additive raw t-template support only, without changing default loader behavior.
- `docs/memory/CURRENT.md` — record completion of the prior read-only comparison and the present diagnostic gate.
- `docs/memory/investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md` — source-audit observations, explicitly distinguished from farm evidence.
- `docs/memory/manifest.json` — generated integrity update.

Proposed new files:

- `src/cuts/pion_t_template_comparison.py`
- `testing/test_pion_t_template_comparison.py`
- `docs/memory/phases/kaonlt-pion-t-template-controlled-comparison-2026-10-19-task-contract.md` — this task contract, transferred into the checkout before Codex execution.

Confirm proposed new paths are absent before creation (except this already-delivered contract). No change to `main.py`, `calculate_yield.py`, `pion_t_bin_parents.py`, `pion_component_fits.py`, production fit configuration, or other physics files is authorized.

A necessary expansion is a contract blocker requiring review, not implicit permission.

## 10. Acceptance and stopping rules

**Local implementation acceptance** requires all focused deterministic checks and unchanged-default-path tests to pass, with an actual-diff review showing no production change.

**Farm scientific acceptance** requires subsequent user-run validation of the detached comparison and independent evidence review. It is not established by local tests, source review, successful PDF generation, or an optimizer's `success` status.

If t-local templates improve a fit, record which component and missing-mass region changes, not merely an aggregate goodness-of-fit statistic. If they do not improve it, retain that negative result. Do not adjust priors, tune thresholds or activate parent anchoring to force a desired outcome within this comparison.

## 11. Deliverable and next gate

Return the complete actual source/memory diff, new-file contents, deterministic test results, memory-health receipt and exact unresolved runtime limitations.

Do not commit, push, or execute Jefferson Lab farm code.

After ChatGPT actual-diff acceptance, user publication, and pushed-state synchronization, the next gate is **one controlled farm comparison of baseline A versus candidate B**.

Only its scientific evidence can justify the following hierarchical global-to-t fitting contract. The `(t,phi)` extension remains deferred until t-level model construction has been validated.
