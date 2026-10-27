# KaonLT — Local (t,phi) accepted-fit reporting and one private Center/Left evidence owner

**Task class:** narrowly bounded farm-facing implementation, **not** a new analysis phase or production promotion.
**Scientific NEXT:** obtain actual child-level Global / Local t / Local (t,phi) comparisons for Q4p4W2p74, Center/Left, low epsilon.
**Observed remote `test` HEAD at audit:** `60f1f4f0e5bc11667b2fecde328c477cbd83926e`. Recheck actual HEAD, `origin/test`, worktree, pending edits, and GitHub before changing anything.
**Status:** `ACTIVE` contract; no farm readiness or ROOT acceptance is asserted.

## 1. Startup, existing work, and stop rules

Read tracked root `AGENTS.md`, then `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, `docs/memory/USER.md`, in that exact order; then `docs/memory/MAINTENANCE.md`, this contract, and only directly relevant owner/comparison contracts. Emit the prescribed full health check before the source-change gate.

**Preserve the substantial *unpublished* prior Work changes.** The supplied October 27 `KaonLT_local_t_phi_20261027_review.zip` is review evidence, not a replacement checkout: it records both a pre-existing dirty worktree and newly implemented default-off `src/cuts/pion_t_phi_template_comparison.py`, `src/binning/calculate_yield.py` hook, tests, renderer and a previous task contract. Independently inspect the **actual local** state and compare changed paths/bytes with that receipt before edits. The live remote has not yet incorporated these files. If identity or pending-work preservation is ambiguous, stop and report the **first** conflict, without resetting, cleaning, staging, deleting, or overwriting anything.

Existing detached Local (t,phi) observer uses the genuine `build_particle_subtraction_component_result` fitter on separate selected (t,phi) child histograms; the existing signed and proton-cleaned caches provide observations, and ROOT/farm fits have **not** run. **Do not revisit feasibility, refit Global/Local t from scratch, revisit DiamondPlot, run C regularization, or introduce a new optimizer.**

## 2. Correction: fitted result is not necessarily an accepted fit

**Target:** `src/cuts/pion_t_phi_template_comparison.py` and `testing/test_pion_t_phi_template_comparison.py`.

The new observer currently assigns `fitted_local_templates` / `fitted_with_template_fallback` immediately after the fitter returns, even when the existing fit acceptance gate would reject the result. Fix only this reporting defect:

1. Call `evaluate_particle_subtraction_component_fit_result(child_fit, inp)` from the established `pion_component_subtraction.py` (or the exact existing gate that the weight builder invokes) after each returned fit.
2. Record **both** independently: `fit_returned`, `fit_accepted`, complete gate reasons/diagnostics, and per-component shape scope (`local_t_phi`, `local_t_aggregate_fallback`, `setting_fallback`, `unavailable`). Distinguish `fit_rejected` from fitter failure, empty child, and accepted with/without fallback; preserve full fit evidence and source support on rejection.
3. Keep `cell_status_counts` truthful and deterministic. No rejected result is labeled a successful Local (t,phi) fit. The existing transfer-weight gate remains authoritative; do not apply rejected child weights or change its rules.
4. Extend tests for returned-but-rejected, accepted-with-genuine-local, accepted-with-fallback, exception/empty, and correct renderer status. Rejected/sparse cells are **expected observations** in this experiment, not a reason to manufacture a coarse fit.

**No changes to fitting physics:** preserve exact component shapes, fit ranges, alignment procedure (including parent-context child alignment), stages, priors, protection of Lambda/Sigma0, weight clipping, source coefficients, cuts, binning, ROOT ownership and statistical conventions. Global refers to historical inherited-template A; Local t to existing B; Local (t,phi) is the new independent child fit. Method A is separate.

## 3. Minimal runtime owner for the actual evidence (not a new launcher)

**Start from the known successful private A/B owner** `testing/run_pion_t_template_comparison_gate.py` and its test `testing/test_run_pion_t_template_comparison_gate.py`. Its currently verified output is four A/B files, and it rejects unknown `LT_ANALYSIS_*` flags. It also strictly constrains committed descendants of historical `REVIEWED_HEAD` and its archive/profile inventory; simply setting a flag will **not** make this farm-ready.

Provide **one** new isolated Center/Left low-epsilon acquisition, extending that owner's already-established single private analysis child. No extra full analysis execution, no alternate farm owner, no bypass of source, isolation, provenance, output, and preservation checks.

1. Enable **both** `LT_ANALYSIS_PION_T_TEMPLATE_COMPARISON=1` and `LT_ANALYSIS_PION_T_PHI_TEMPLATE_COMPARISON=1` **only** inside the authenticated private analysis child, and reject other unauthorized opt-ins. Preserve the existing prepass exclusion and original env guard.
2. Keep the exact four required A/B files. Also require **exactly two fresh child JSONs**, named `Center_kaon_Q4p4W2p74_lowe_pion-t-phi-template-comparison.json` and `Left_kaon_Q4p4W2p74_lowe_pion-t-phi-template-comparison.json`. Do not require ROOT to produce a child PDF: the already implemented `testing/render_pion_t_phi_template_comparison.py` renders captured JSON offline after artifact acceptance.
3. Update the owner's exact artifact declaration, freshness/no-collision checks, complete-profile contract, strict collector manifest and ZIP inventory, byte-size/SHA identity and returned verification. Prefer a **versioned new profile** instead of silently redefining the historical A/B v1 profile and its archives. The same six objects must be collected from the same private run.
4. For each child JSON, independently verify schema, source commit, source-file identities against actual committed blobs (read the exact required source inventory from the observer), Center/Left identity, canonical `3 × 9` cell inventory, unique `(t_index,phi_index)` keys, edges, child/status accounting, `parent_input_closure=PASS`, signed source/phi/`Sumw2` closure, protected-model status and explicit fit acceptance. Reconcile the corresponding A/B JSON's original **runtime-container** data/dummy identity with this child JSON's `runtime_inputs`; never substitute original pre-shift files for owner-private modified containers. Fail on any mismatch.
5. Verify rejected or unsupported child fits are honestly reported as unavailable for transfer, **without** treating their scientific non-acceptance as a file/provenance failure. Separate collection success, fit success, and physics adequacy. No fallback may be represented as a genuinely local component.
6. Expand the owner's **exact allowlist** to cover only the genuinely reviewed, necessary unpublished scientific and test paths that this new committed acquisition consumes. Maintain strict Git ancestry, committed-blob identity, independent source-worktree inspection, frozen cleanup-repair pins, checked profile and exclusion of unrelated source changes. **Do not** replace guards with wildcards (`src/`, `testing/`) or weaken checks just to pass on a dirty checkout.
7. Preserve owner-private staging of pion low-epsilon Data/Dummy inputs, all five standard setting directories needed by the existing prepasses, original input hash/inode preservation, trusted volatile placement, no changes to ordinary checkout or installed `ltsep`, output isolation, collector identity, and bounded cleanup. No new scans of existence-only skim ROOT files or additional unnecessary large I/O.
8. Accept no old/stale output by filename existence alone. Failure must stop before packaging or interpretation and report the **first failed invariant**. Do not automatically rerun the farm.

**Likely owned paths:**
- `src/cuts/pion_t_phi_template_comparison.py`
- `testing/test_pion_t_phi_template_comparison.py`
- `testing/run_pion_t_template_comparison_gate.py`
- `testing/test_run_pion_t_template_comparison_gate.py`
- New versioned `testing/*t_phi*gate_profile*.json` **only if strictly needed**
- This standalone task contract, `docs/memory/CURRENT.md`, `docs/memory/manifest.json`.

The previously unpublished `src/binning/calculate_yield.py` guarded observer hook and renderer are already implementation inputs, not a license to rewrite them. Touch one of them only to fix a proven contract mismatch, document why, and include it in the exact actual-diff review. Treat all other pending files as preserved, not automatically owned. Do not edit `USER.md`; preserve the pending Work-first proposal byte-for-byte.

## 4. Deterministic validation before actual-diff review

- Demonstrate that the default-off production AST / runtime path remains unchanged; no observer/fit import in ordinary runs or prepasses.
- Cover the observer's **genuine** 54-child fixture paths (27 cells each for Center/Left), positive/negative fit gates, rejected-weight exclusion, fallback attribution and parent/child reconstruction; distinguish synthetic ROOT stubs from farm ROOT validation.
- Cover exact owner dual-flag behavior, unknown opt-in rejection, no second analysis child, duplicate/stale/missing/wrong-cohort child JSONs, fake/changed source hash, wrong private input identity, wrong child count/coordinates/status, rejected-but-recorded child outcomes, ZIP-profile mismatch, preserved original A/B collection and historical source-guard negatives.
- Rerun affected old comparison/owner/event regression tests and new tests. Report environment, executed/passed/skipped/failing counts accurately. Run compile/grammar checks, `git diff --check`, memory-manifest update/check and memory-health check.
- Verify normal configuration and all accepted source/evidence/phase statuses remain unchanged except accurate ACTIVE/NEXT synchronization for this task. Preserve historical task contracts and artifacts.
- Include **complete tracked diffs plus complete new/untracked files** in one independent review bundle, with a before/after full `git status --short`, hashes of initially pending files, and list of changed paths. Do not stage merely to prepare the review.

**Stop conditions:** failed source identity, unexpected dirty worktree changes, failed tests, changed frozen physics, incomplete provenance, unsafe output/path/copy changes, or any rejected checker. Report precisely which invariant failed; do not claim farm-ready or offer a farm command.

## 5. Handoff, evidence and physics output

Codex performs source changes and local deterministic checks only; **no commit, push, farm execution or environment-dependent ROOT claims**. Then ChatGPT independently reviews the **actual diff**. User alone commits/pushes after acceptance; ChatGPT verifies pushed HEAD/blob identity and complete readiness before giving any farm command.

**Once one accepted private Center/Left execution exists**, return its full verified six-artifact bundle for plot-first analysis. Use the tracked JSON-only renderer to generate real full-MM, signal-window, high-MM shoulder, cell residuals/fit components, signed A/B/Bphi transfer overlays, per-cell statuses and t×phi heatmaps. Label conditional fixed-weight errors versus fit-parameter covariance; do not invent inherited child control-fit curves or claim rejected kaon subfits are accepted. No acceptance, promotion, physical bias/MSE or L/T physics conclusion is authorized by fitted curves alone.

**Sole `NEXT`:** return Codex's actual owner-plus-status-fix diff and deterministic test receipts for independent review. No separate speculative investigation or preparatory farm attempt.


## Bounded independent-review repair - accepted application evidence

User authorized this repair on October 27 at unchanged test/origin/test
60f1f4f0e5bc11667b2fecde328c477cbd83926e. Reproduced review finding:
verify_child_payload accepted an available child transfer with no arrays, and
one with negative variance/missing allcuts. This appendix extends the same
contract; it creates no new phase or owner.

Own only testing/run_pion_t_template_comparison_gate.py and its existing test,
this appendix, CURRENT and generated manifest. Preserve all initial pending
bytes outside these five paths, especially observer/fitter/renderer/profile,
USER and previous review artifacts. Do not adjust physics or fit acceptance.

Require the exact Global/Local_t/Local_t_phi application inventory. An honest
unavailable result has a nonempty reason and no weights/cuts. Every available
result needs finite axis-sized weights with zero flows, both allcuts/nommcuts
spectra, finite signed contents and finite nonnegative variances. Check each
against captured child pion observations: H_applied=H_observed*w and
V_applied=V_observed*w*w, using the existing numerical tolerance. Do not refit
or introduce clipping, thresholds, scaling or a new estimator. Preserve the
rejected-child prohibition. No claim of fitted covariance or physical adequacy.

Regress both reproduced defects and missing/extra methods or cuts, arrays of
wrong sizes/types, nonfinite numbers, invalid flow weights, signed content and
quadratic variance mismatches. Include valid signed/cancellation/zero-transfer
examples and an application generated by the actual scalar producer with only
ROOT-dependent weight construction stubbed. Run owner tests in WSL plus the
existing affected observer/numerical/event-owner regressions. Verify syntax,
complete actual incremental/cumulative diffs, manifest/health and preservation.
Stop at any required frozen-source change. No stage, commit, push, farm command
or retry. NEXT remains independent actual-diff review of this bounded repair;
the previously unsettled publication dependency scope remains explicit.


### Bounded repair local outcome

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING. The old gate failed 70 new
negative subcases across two regression methods before the repair. After it,
71 WSL owner tests pass. Existing Windows observer/comparison/alignment/
statistics/renderer/event-owner suites: 131 cases, 118 passed and 13 skipped
(11 require PyROOT; two require historical artifacts). Combined distinct
inventory: 202 cases, 189 passed, 13 skipped, zero failures.

Six new owner test methods cover valid applications from the actual scalar
producer (only ROOT-dependent weight construction stubbed), both cut masks,
all three method names, exact method/cut inventories, finite/axis-sized
arrays, zero flows, negative variances, signed/quadratic reconstruction,
honest unavailable states, rejected fits, cancellations and zero weights.
The old owner fixture now includes all three explicit method statuses, just
as the real observer always did. The repair does not infer fitted-parameter
covariance, calibrate uncertainties, refit physics or modify scientific source.

The exact five-path incremental repair and complete cumulative unpublished
diff accompany the preservation/health receipt. Previous review artifacts,
USER and all initially pending paths outside this scope are preserved.
No index/ref mutation, publication or farm execution is part of the repair. NEXT is external
actual-diff review; this author's audit does not grant independent acceptance.
The prior publication dependency/hunk scope is still unsettled and must not
be bypassed by broadening the committed-source allowlist. No farm readiness.
