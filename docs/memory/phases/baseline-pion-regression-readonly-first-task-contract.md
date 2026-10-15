# KaonLT — temporary baseline pion-subtraction regression priority and read-only first-divergence audit

Status: `ACTIVE` — one narrow memory correction plus read-only source/evidence investigation; **no production source changes, farm execution, or physics acceptance**.

## Authority and startup

Repository: `trottar/lt_analysis`, live `test` HEAD at authoring `d860eeafbb889bcc04cfddb1ff6019cfcff50922`. Independently verify live remote `test`, exact local HEAD and complete index/worktree/untracked state before editing; stop and re-anchor on unexpected drift. Read root `AGENTS.md` followed *in exact order* by `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, and `docs/memory/USER.md`. Follow `docs/memory/MAINTENANCE.md`, `CODEX.md`, `COMMUNICATION.md`, and the existing October 14 regression investigation. Keep `CURRENT.md` the sole active-state/NEXT authority. Preserve all unrelated local files and review bundles, including any prior untracked isolation-task contract; do not delete, reset, stash, clean, or stage for review.

## User's corrected scientific objective (must be unmistakable)

The normal, already approved KaonLT scientific roadmap is **temporarily interrupted**, not replaced, to explain a discrepancy in the **existing baseline pion-background subtraction**. The procedure worked previously in some scope; identify and verify which *same-scope* behavior worked rather than asserting that a setting-wide fit proves a per-t fit was correct. The immediate target is Q4p4W2p74 / Center / low epsilon / canonical t2; adjacent Center t1 and the setting-wide fit are comparison controls with scope distinctions stated.

Historical August 26 (attributed **MEMORY ONLY**, not newly reproduced): staged pion-control `B_n=0`; refined `B_n≈0.0002021`; kaon-side `A_n≈0.0369527`; transfer weight max ≈183; fitted pion-background integral ≈0.2529 versus event-weighted applied ≈6.556, ratio ≈25.9. The 67-page available Center/lowe PDF is *setting-wide*, reporting 0.5967 fitted versus 0.6010 event-weighted; **not an identical-input per-t regression**. Do not invent the nonexistent “(4)” PDF, an original numerical fit branch, or a historical working SHA. The actual cause of staged `B_n=0` remains **NOT VERIFIED**.

**Primary objective:** locate the *first numerical divergence* between earlier genuinely comparable working behavior and the problematic subtraction: signed pion-control inputs -> per-t SIMC templates and fitted component amplitudes -> fitted control/kaon model ratio -> event-weighted application. No change to fits, PID, cuts, input population, signed subtraction, template mixing, normalizations, shifts, priors, weights, acceptance, statistics, yields, L/T, SIMC, or `no_empirical_residual`. Do not force positive `B_n`, clip weights, or redesign physics without a demonstrated source/evidence defect and a later explicitly authorized narrow fix.

**Return-to-roadmap requirement:** once the discrepancy is explained and the existing subtraction validated, or a necessary minimum correction is independently reviewed, user-pushed and **user-farm validated**, resume the normal `docs/memory/roadmap/STATUS.md` at the interrupted position. Preserve all closed, active, deferred, blocked and pending roadmap items, scientific roles and existing acceptance. This investigation is **not a new phase**, does not promote Method A or Method B, and must not permanently displace the normal roadmap. If evidence does not yet establish the root cause, remain `ACTIVE` and report the single missing decisive datum rather than claiming resolution.

## Current operational boundary

The normal-launcher output-isolation blocker is **real** but **only a prerequisite for a future run that actually needs the normal launcher**. Existing `baseline-pion-center-t2-anchor-farm-owner-task-contract.md` remains `BLOCKED`; the separate attempt-output isolation proposal is preserved, not launched by this task. Do not prioritize further isolation engineering before the read-only comparison of available evidence and source history. The user alone commits/pushes and runs all Jefferson Lab farm validation. Codex never runs farm, commits, or pushes.

## Allowlist

1. `docs/memory/phases/baseline-pion-regression-readonly-first-task-contract.md` — this exact delivered standalone contract.
2. `docs/memory/CURRENT.md` — change only the **Active Objective**, relevant blocker wording if needed, and the **single NEXT**, preserving the rest of accepted status.
3. `docs/memory/investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md` — narrow addendum documenting the corrected task framing, source-evidenced observations, uncertainties, and first-divergence plan/results.
4. `docs/memory/manifest.json` — regenerate/check after the memory edits.

No source, tests, plotting, scripts, existing operation owners, physics configuration, `roadmap/STATUS.md`, `MEMORY.md`, `USER.md`, `CURRENT_HANDOFF.md`, or accepted evidence modifications. A blocked implementation or lack of historical comparators does not authorize broadening this allowlist.

## Read-only analysis to perform *before* completing memory text

1. Search existing available records, logs, PDFs and commits for an independently identifiable **last demonstrably good and first bad same-scope per-t** fit/application. State exact `HEAD`, setting, epsilon, t/phi selection, cuts, normalization, fit scope, and source identities. If none exists, say `NOT VERIFIED`, don't infer from the setting-wide PDF. Do not ask user for artifacts already present.
2. Trace current `src/cuts/pion_component_fits.py`, `src/cuts/pion_component_subtraction.py`, `src/cuts/rand_sub.py`, `src/cuts/pion_t_bin_parents.py`, `src/utility/pion_component_shapes.py` and `src/utility/background_config.py`; inspect the August 19 `f924a95` and August 20 `ee9e292` transitions and the immediately preceding commit. Compare **real changed source**, not commit messages or summaries, and name the earliest differing producer/consumer/selection if supportable.
3. Record explicitly that the code's fitted-model weight is `w(M)=K_model(M)/P_model(M)` while the event product is `w(M)*P_signed_data(M)`; the model-event difference is `w(M)*(P_signed_data(M)-P_model(M))` in supported/unclipped bins. A near-zero fitted denominator can amplify a fit/data mismatch; this is a **SOURCE VERIFIED mechanism**, **not a demonstrated historical root cause**. Independently recheck the formula and lineage before recording. Keep non-supported/clipped cases separate.
4. Determine from actual available source and evidence whether first staged `B_n=0` is attributable to nonpositive LS weighted numerator, invalid/no template support, or unavailable original inputs. Inspect fit mode (`least_squares` versus `window_integral`), fit bins, signed source decomposition, SIMC shape and amplitude scaling and any populated-phi aggregation if evidence exists. Do not interpret zero amplitude alone as a negative numerator, and do not confuse per-t diagnostic with setting-wide results.
5. Finish with a compact findings matrix: `SOURCE VERIFIED`, historical `MEMORY ONLY`, hypotheses `INFERENCE`, and critical `NOT VERIFIED`; identify exactly **one earliest missing artifact or comparison** if unable to resolve. No new diagnostic implementation or farm run authorized here.

## Checks and handoff

Run manifest regeneration/check, ordinary memory health, `git diff --check`, `git status --short --untracked-files=all`, verify no science source changes, preserve existing unrelated status/index/refs/bytes and note nonblocking CURRENT size warning. A new untracked contract missing from the initial manifest is expected before regeneration; never bypass final integrity checks. Return a complete actual review diff including the full untracked contract addition plus the tracked memory changes. Report the source audit separately, grounded to exact live source identity. ChatGPT independently reviews the **actual diff**. Codex does not stage, commit/push, or run farm.

Only after ChatGPT accepts the complete diff, the **user** executes this *unexecuted* scoped publication block:

```bash
git add -- docs/memory/phases/baseline-pion-regression-readonly-first-task-contract.md docs/memory/CURRENT.md docs/memory/investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md docs/memory/manifest.json
git commit -m "Prioritize read-only baseline pion subtraction regression and roadmap return"
git push origin test
```

ChatGPT then checks pushed `test` HEAD, exact changed blobs, CURRENT/NEXT and manifest; continue the read-only discrepancy investigation immediately, not an administrative or isolation milestone. Any new farm gate requires separate source-reviewed, pushed readiness and is executed by the user only.
