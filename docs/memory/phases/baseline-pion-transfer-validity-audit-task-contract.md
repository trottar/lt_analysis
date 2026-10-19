# KaonLT — setting-independent baseline pion-transfer validity audit

Status: `ACTIVE` — narrowly scoped **detached, read-only scientific diagnostic**. No new phase, source-physics adjustment, or farm run is authorized by this contract. Its purpose is to make the existing model-to-signed-event transfer defect measurable, reproducible, and setting independent before any corrective production proposal.

## Authority and startup gate

Repository: `trottar/lt_analysis`, `test`. Authoring-time remote HEAD observed: `c360bf591881d017bd0b5bf0c89d73dd780b1d9a` (October 18, 2026); a starting observation, **not** permission to assume the current local HEAD. Prior to edits, independently establish the actual remote `test`, local branch/HEAD, complete index/worktree and relevant untracked paths. Read root `AGENTS.md` and its five-file startup core in exact specified order. Read `docs/memory/MAINTENANCE.md`, `docs/memory/CODEX.md`, CURRENT's baseline investigation and linked October 18 evidence, and the exact source owners below. If identities or dirty state conflict materially, STOP `BLOCKED`; never reset, clean, overwrite, commit, push, or execute the farm.

## Scientific observation requiring this task

October 17 Center/low-epsilon canonical-t evidence exhibits closed component-model algebra but poor transfer from signed pion-control events to the applied subtraction. The existing source:

- `src/cuts/pion_component_subtraction.py`: builds final ratio `w_i=K_model_i/P_model_i` if denominator exceeds its numerical floor; otherwise writes zero, and applies this weight to frozen signed source events. Its present checks establish finite/template/algebra support, **not** scientific model-to-observed-control agreement.
- `src/cuts/pion_component_fits.py`: staged and regularized joint fits and fit validation; common configuration in `src/utility/background_config.py` leaves some rejection limits `None`, which may report validation accepted despite large reduced chi-square. Preserve this configuration for this contract.
- `src/cuts/pion_background_chain_audit.py`: already collects full signed source spectra, Sumw2, transfer histogram and coefficients, event-applied products, closures, and exact parent identities. Its `is_target` is intentionally hardcoded for the historical Center investigation; **do not change this producer's scope in this task**.

The existing Center seven-stage output has `schema_version=pion_background_chain_audit/v1` and `status=captured_pending_independent_validation`. This permits diagnostic study only, not claims of accepted farm production or an original historical root cause. Existing Center/t1, t2 and t3, and independently observed Center/Left pion/SIMC occupancy, are context, **not special-case logic**.

## Allowed changes

1. New `testing/analyze_pion_transfer_validity.py`: self-contained pure-Python validator/CLI accepting an explicit existing audit-JSON path and optional JSON report path. No ROOT, PyROOT, farm-specific imports, `main.py`, or input ROOT files. It must operate on any compatible `audit.parents[]` collection, irrespective of `Q2`, `W`, `EPSSET`, `phi_setting`, or the number and labels of `t` bins.
2. New `testing/test_analyze_pion_transfer_validity.py`: deterministic pure-Python tests with compact generated fixture dictionaries and a temporary JSON fixture. Tests must not depend on private farm artifacts.
3. This one Markdown contract under `docs/memory/phases/` and minimal synchronized `docs/memory/CURRENT.md`/existing baseline investigation entries **only if** necessary to document an actual reviewed work-state/evidence change. No other production, utility, `src/`, farm launcher, existing tests, or analysis-configuration file may change.

4. `docs/memory/manifest.json` solely for required regeneration with the canonical manifest generator (user-approved allowlist addition); confirm this contract is the only added inventory entry and existing file hashes remain unchanged.

If a narrow source helper outside this allowlist becomes necessary, STOP and report, rather than broadening ownership silently.

## Diagnostic definitions and validation contract

For each stored canonical parent, use only coherent same-axis *normal* missing-mass bins and the existing flow conventions. Read:

- `P_obs[i]` and axis from `parent.control.values/edges`, independently cross-checked against `parent.source_spectra.combined.signed`; use signed contributions, **never absolute-replace negative random/dummy events**.
- `V_obs[i]=parent.source_spectra.combined.sumsq[i]`, plus `absolute` and `count` support. Check finite and nonnegative Sumw2 and source-control closure; preserve underflow/overflow distinction.
- `P_model[i]`, `K_model[i]`, `w[i]` from `parent.transfer.histograms`, `parent.transfer.weights` and `parent.transfer` numerical policy, with per-bin denominator support and clipping state; do **not** fabricate a ratio where the denominator is unsupported.
- `C_applied[i]` from `parent.final.histograms.H_pion_subtraction_template_MM_nosub` and its stored independent signed-event reproduction checks. Never substitute `proposal` values as production acceptance.

On **qualified, un-clipped, denominator-supported** bins test the identity
`C_applied[i]-K_model[i] = w[i]*(P_obs[i]-P_model[i])` at scale-aware *floating-point-only* tolerance; a failed algebra/identity/source/axis/provenance invariant is a **hard stop for downstream interpretation**. Do not call arithmetic PASS scientific agreement.

For any requested region assembled exclusively from native histogram axes or metadata-defined windows (no assumed `0.7–1.3`, `1.10–1.16`, fixed Q2/W or other magic coordinates), report:

- Signed sums `sum(P_obs)`, `sum(P_model)`, `sum(K_model)`, `sum(C_applied)`; source counts, `sum |source weights|`, `sumw2`, and meaningful signed/absolute effective support and cancellation fraction (handling zeros without division-by-zero).
- `D_R=sum_i (C_applied[i]-K_model[i])`, and `L_R=sum_i |C_applied[i]-K_model[i]|` **without signed cancellation**. Report these both directly and as the reconstructed `w*(P_obs-P_model)` identity on qualified bins only; report unsupported/clipped bins separately.
- Conditional **fixed-transfer-weight** uncertainty `sigma_R=sqrt(sum_i w_i^2*V_obs[i])` on qualified bins, and `D_R/sigma_R` when defined, marked explicitly **diagnostic conditional pull, not a formal fit p-value**. It omits pion-fit covariance, SIMC-template uncertainty, fitted-event correlation, and transfer-weight uncertainty. Do not use this statistic alone to grant scientific acceptance.
- Normalized event/model discrepancy measured against `sum |K_model|` where nonzero, and an absolute discrepancy measure; keep magnitudes visible when the kaon-side pion model is very small, so ratios alone cannot cause false scientific decisions.
- Counts of bins with `K_model>floor` but unsupported `P_model`, denominator-floor touches, clipping, boundary-hit components, fit amplitude uncertainties, joint fit quality and producer validation flags (report `accepted` as an *existing flag*, not a newly earned result).

Do not introduce a numerically chosen weight clip, neutron floor, rescue amplitude, empirical residual, or any automatic fallback. All true numerical tolerances must be scale-aware; any optional warning threshold must be universal, named, and strictly presentation-only. When necessary covariance, windows, signed source support, or provenance are missing, report `NOT VERIFIED` for that dimension; do not mark scientific PASS. The report must state producer provenance limitations and that Method A/B have no authority here.

## Scope and dynamic invariance tests

Tests must exercise: exact identity; a model-algebra-PASS/event-model-FAIL fixture; positive numerator with unsupported/zero denominator; tiny positive denominator/high leverage; negative random/dummy signed cancellation; insufficient/zero Sumw2; clipping; underflow/overflow; missing or mismatched axes/identity and provenance; nonzero but very small kaon-side model; parent count other than 3; arbitrary kinematic tokens and epsilon/phi labels; variable missing-mass edges and fit windows; independent rescaling of all *physical event yields and corresponding variances* without changing dimensionless diagnostic conclusions. No setting-dependent branch or frozen-input mutation is allowed.

The acceptance criterion for this implementation is **correct read-only risk attribution and numerical closure**, not a new production validity threshold or correction. The future production decision would need an independently calibrated uncertainty budget including fit/template covariance, appropriately defined statistical criteria, actual-diff review, user publication and applicable farm validation. Do not use historical accepted E.8/F.6 status as authorization to revise their physics.

## Deterministic local checks and handoff

Execute `python -m py_compile testing/analyze_pion_transfer_validity.py testing/test_analyze_pion_transfer_validity.py`, `python -m unittest testing.test_analyze_pion_transfer_validity`, `git diff --check`, and inspect full actual diff and new-file contents. A detached invocation on the supplied October 17 JSON is optional only if its file is accessible; label resulting observations diagnostic-only and never smuggle the farm paths or values into code/tests. Show source SHA, worktree status, changed-path allowlist, test receipts, sample synthetic output, and failure/unknown semantics. Codex must **not** commit, push, run the farm, or promote production.

After ChatGPT independently reviews the actual diff, the user owns commit/push. ChatGPT then synchronizes pushed state; the substantive scientific NEXT remains interpretation of the baseline pion-control transfer and whether uncertainty-calibrated criteria justify a minimum dynamic correction.
