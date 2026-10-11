# KaonLT E.8 — Method-A versus reweighted-baseline scientific-role memory correction

## Classification and objective

**Task class:** Narrow repository-memory scientific-ownership correction after a user-identified interpretation error. **No new scientific phase, computation, physics implementation, four-way plot, farm run, F.1 archaeology, or failed F.4 owner retry.**

**Observed remote `test` HEAD at contract preparation:** `78ecadbacf58ab25b6a3670137633b3124074fe7`. Before editing, independently establish actual local `test` HEAD, `origin/test`, index and complete `git status --short --branch`. The last user-reported worktree contained three unrelated modified files (`docs/memory/LEARNINGS.md`, `docs/memory/evidence/e8-4-canonical-five-fresh-f1-lineage-materialization-2026-10-05.md`, `docs/memory/phases/e8-4-canonical-five-fresh-f1-candidate-lineage-refresh-task-contract.md`) and four untracked `kaonlt_review*.diff` files. **Preserve all of these exactly.** If branch/HEAD/ref unexpectedly differs, or an allowed edit would overwrite unexpected local changes, stop `BLOCKED` with the exact identity/path difference; do not reset, pull, stash, stage, clean, commit or push.

Follow root `AGENTS.md`: first read the full startup core, in this exact order:
1. `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read `docs/memory/MAINTENANCE.md`, `docs/memory/CODEX.md`, `docs/memory/TOOLS.md`, the relevant E.8 decision/roadmap records, and only source/evidence needed for this correction. Emit the mandatory full health check before work.

## User clarification and exact scientific ownership

The user explicitly corrected the previous handoff's misunderstanding: **Method A independently derives a refinement response/weight, and the purpose of the reweighted baseline is to apply that Method-A refinement to the existing baseline analysis to improve its pion-background estimate.** Do not describe the reweighted baseline as a missing independent *third data-analysis producer* or demand an unrelated new formula.

Preserve the distinction, confirmed by live tracked source:

- **Baseline:** Existing pion-background transfer `w0` and its baseline after-pion-subtraction `Y0/MM_0`. Baseline production remains unchanged.
- **Detached Method A:** F.3 derives a relative HGCer detector-response function `r` on the accepted `hgcer3` basis from its own positive-response training inputs. `r` is **not** an absolute pion-mis-ID probability, a stand-alone extracted kaon yield, or a complete replacement background normalization.
- **F.4 parent-preserving Method-A correction:** From signed baseline contributions `b_j = s_j*w0_j`, define `B = sum_j b_j`, `U = sum_j b_j*r_j`, `C_j = r_j / (U/B)` for valid parent sums; the accepted invariant is `sum_j b_j*C_j = B`. This means the **relative response is independently determined**, while the **parent normalization of applied `C` depends explicitly on baseline `w0` and signed source coefficients**. Do not claim the final normalized correction is numerically independent of baseline. No independent `(t,phi)` child normalization.
- **Method-A-reweighted baseline:** The F.6.3 parallel non-production full-analysis branch fills `B_pi_A` using signed source coefficient times `w0*C`, forms `MM_A = pion_input - B_pi_A`, and extracts `YA` using the existing uncertainty helper. **Existing `YA/MM_A` are already the final reweighted-baseline result** in source/API naming, not a separate "Method-A-only" final kaon result.
- **SIMC:** Existing independent K+Lambda SIMC per-cell MM and stored yield support. The absolute-normalization/unit blocker remains in force.

Do not conflate **Method-A relative weight/response**, **baseline plus Method-A correction**, and **extracted final kaon signal yield**. A standalone "Method-A-only" final kaon `Y/MM` pair is **NOT VERIFIED** as defined/produced; it must not be conjured by relabeling `B_pi_A` or `YA/MM_A`. The previously requested four-category presentation may distinguish Method-A weight/diagnostic evidence from final yield/MM curves until genuinely distinct comparable final objects are defined. Preserve the user's original intent, but do not promise four independent physical final-yield/MM series on the current evidence.

**Source anchors to trace, without altering:**
- `src/cuts/pion_hgcer_method_a_acceptance_map.py` (independent F.3 relative response)
- `src/cuts/pion_hgcer_method_a_parent_preserving_correction.py` (F.4 `b`, `r`, `C`, signed parent preservation)
- `src/cuts/pion_component_subtraction.py` (`fill_simc_shape_pion_subtraction_templates()` with `method_a_event_multipliers`)
- `src/binning/calculate_yield.py` (`_build_f6_3_parallel_method_a_source()`, `B_pi_0/B_pi_A`, `MM_0/MM_A`, stored `Y0/YA`)
- `src/cuts/full_background_subtraction_plots.py` (E.8.4 downstream consumer)

**Evidence limits:** Accepted F.6.3/E.8.4 Left/lowe execution at farm source `da38444e7aa60efd62d6638780776344daf40276` remains `CLOSED / RUNTIME VALIDATED` for that narrow scope, not full canonical five or Method-A production promotion. Baseline-five E.8.2 remains `CLOSED / RUNTIME VALIDATED`. The canonical-five owner failed before analysis at `f6_3_scientific_equivalence_mismatch:f4:$.parents[0].absolute_baseline_parent_sum`; preserve `BLOCKED`. Microscopic F.4 mismatch investigation and further F.1 archaeology remain `DEFERRED` by user priority. Absolute SIMC units remain `BLOCKED` by `SIMC_normfac_luminosity_and_charge_units_not_source_proven`. Method B is numerically excluded.

## Allowed memory files only

1. `docs/memory/CURRENT.md` — correct active misunderstanding, true blocker distinction, and one push-stable NEXT.
2. `docs/memory/MEMORY.md` — record concise, reusable definitions of `r`, `C`, and `w0*C`, and warn against treating `YA/MM_A` as "Method-A-only" or inventing a reweighted-baseline producer.
3. `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md` — correct four-way table/producer interpretation and old "missing reweighted baseline" conclusions; retain the user request as an **intent** and distinguish established output from not-yet-owned standalone Method-A final output. Record that prior bounded source-audit observations were valid as literal key inventories but the inference from naming was wrong; do not erase history, rewrite past runtime identities or say current source has fourth result.
4. `docs/memory/roadmap/STATUS.md` — correct only the false reweighted-baseline missing-ownership wording while preserving phase order and all accepted statuses.
5. `docs/memory/manifest.json` — generated integrity metadata only.
6. `docs/memory/phases/e8-method-a-reweighted-baseline-role-memory-correction-task-contract.md` — this one standalone contract.

**Frozen:** All `src/`, all `testing/`, all evidence records, old phase/implementation contracts (especially `e8-four-way-final-yield-mm-implementation-task-contract.md`), `AGENTS.md`, user/handoff records, `LEARNINGS.md`, F.1 archaeology, profiles/checkers/collectors/owner scripts and artifacts. Existing uncommitted/untracked files are not included. If an additional source or evidence record must change to avoid material conflict, STOP and report it, rather than broadening this task.

## Exact before and after

**Before:** `CURRENT` and E.8 roadmap assert an independent `Y_rw/MM_rw` producer is missing; the source audit labels `YA/MM_A` solely "Method A", leading to the false requirement for a second independently produced reweighted baseline. Actual F.6.3 applies `w0*C` and already owns the reweighted-baseline final outputs.

**After:** The active and durable memory unambiguously states:
- Method A's *relative* detector-response refinement is separate from baseline weight construction.
- Parent-normalized `C` uses signed baseline support and is not independent of `w0` at that step.
- F.6.3 private `YA/MM_A` equals the **Method-A-reweighted baseline final branch** in the existing source nomenclature.
- A standalone Method-A-only **final yield/MM** is not established; Method-A detached response/weights/diagnostics do exist.
- Only genuine final-yield products may be plotted as final-yield series; Method-A-only weight diagnostics may be shown as diagnostics without false yield units.
- All closed/blocked scopes and scientific freeze stay unchanged.

**Push-stable `CURRENT` NEXT (conditional on actual-diff PASS, user commit/push and pushed-state synchronization):** Review the existing detached **Method-A weight/response** scientifically: distinguish `r` and parent-normalized `C`; trace their input populations, range, shape, uncertainties/support and connection to `w0`, and identify exactly which standalone weight diagnostics can be compared with baseline without inventing a new final-yield producer. This is the requested next substantive physics task. Do not automatically start another correction implementation or farm run. The deferred F.4 numerical discrepancy archaeology remains deferred and the failed private owner remains blocked.

## Positive, negative and regression checks

**Positive:**
- Every active mention of an absent independent *reweighted-baseline* producer is corrected in the allowlisted records; `YA/MM_A` is explicitly mapped to `w0*C`.
- The four-category *presentation intent* is retained with current object availability/semantics distinguished.
- Method-A relative response independence and F.4 baseline-dependent signed parent normalization are BOTH explicit.
- `CURRENT` has precisely one substantive NEXT focused on Method-A weights.
- Preserve evidence links, statuses and prior accepted narrow closures.
- No stale statement claims an independent final Method-A-only `Y/MM` exists.

**Negative:**
- No relabeling `B_pi_A` as a final kaon yield; no `r` as absolute probability or `C` as baseline-independent.
- No new output key, formula, correction, ROOT object, standalone normalized Method-A yield, production branch or SIMC amplitude claim.
- No automatic F.4 owner retry, F.1 archaeology, equality-tolerance relaxation or F.6.4 promotion.
- No accidental staging, overwrite or deletion of unrelated local changes/review diffs.
- No unsanctioned rewrite of the old task contract or evidence as though prior accepted farm data changed.

**Regression:** No production, simulation, plotting, tests, runtime/profile or physical numerical behavior is changed. Retain `no_empirical_residual`, Method B exclusion, parent preservation, frozen binning, yields and errors. Existing page names/API `YA/MM_A` remain unchanged; memory may explain semantics without code alias.

## Required local validation and review return

Discover `<PYTHON>` from the actual local environment. After memory edits:
1. `<PYTHON> -B tools/update_memory_manifest.py --root . --write`
2. `<PYTHON> -B tools/update_memory_manifest.py --root . --check`
3. `<PYTHON> -B tools/check_memory_health.py --root .`
4. `<PYTHON> -B tools/memory_bootstrap.py --root . --json`
5. `git diff --check` and scoped changed-path inspection.

Report exact interpreter, command outcomes, manifest PASS/FAIL, actual `CURRENT.md`, `MEMORY.md`, `CURRENT_HANDOFF.md` byte counts, all health warnings with blocking/nonblocking classifications, and any material contradiction. Nonblocking size/history warnings can be batched at the next milestone; hard integrity and material ownership warnings STOP.

Generate **one complete temporary root-level review bundle** under a non-conflicting name, e.g. `kaonlt_review_e8_method_a_role_memory.diff`, including every tracked allowlisted diff and **full additions for this new contract even if untracked** (`git diff --no-index /dev/null ...`). Do not overwrite existing review bundles. Supply paths, bundle size/SHA-256, `git status --short --branch`, and the exact user-controlled scoped Git handoff commands **without executing them**. Do not stage, commit, push, or run JLab farm.

## Acceptance boundary / hard stop

Codex's summary and local tests do not establish ChatGPT actual-diff acceptance. The user commits/pushes only after that review; ChatGPT then independently synchronizes the pushed state. Only afterward start the detached Method-A weight science task. This is a memory **semantic correction**, not new science/runtime acceptance, a completed four-way yield plot or a farm authorization.

**Hard stop:** HEAD/ref mismatch, conflict with unexpected local changes, loss of accepted evidence/state, new independent final-yield claims unsupported by source, changed scientific code, unresolved blocking manifest/health results, or need to edit outside the allowlist -> `BLOCKED` with exact first failed invariant, no workaround.
