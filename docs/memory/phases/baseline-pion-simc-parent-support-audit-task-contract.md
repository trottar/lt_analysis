# KaonLT — detached pion-SIMC canonical-parent support audit

Status: `ACTIVE` — one narrowly scoped, read-only scientific measurement. This contract does **not** authorize production correction, full analysis execution, new SIMC generation, or a farm run. `docs/memory/CURRENT.md` remains the only authority for active NEXT.

## 1. Scientific question and evidence boundary

Determine whether the *existing setting-wide* `pi_n`, `pi_sidis`, and `pi_delta` SIMC missing-mass templates have adequate statistical and shape support in the canonical-t parents of **Q4p4W2p74 / Center and Left / lowe**, concentrating on Center t2 and using Center t1/t3 as controls. This is the first missing measurement in the October 17–18 baseline pion-model versus signed-event subtraction audit.

The supplied `KaonLT_loweps_20261018_093137.zip` contains a Center seven-stage JSON/PDF and Center/Left analysis/fit PDFs, but **not** the component SIMC event trees. The Center JSON retains a large t2 discrepancy (model ~0.2529 versus applied signed-event background ~6.5562), near-zero refined pion-side neutron amplitude (~0.000202), and high transfer weight (~182.82). The original signed data, staged/joint fits, and transfer identities are investigated already. These values are context, **not new acceptance criteria or fixed test fixtures**. Existing SIMC shape diagnostics show setting-wide populations and zero constructed t-binned shapes. The absence of parent-specific SIMC event-level evidence is `NOT VERIFIED`, not proof that the setting-wide shapes are wrong.

Do not repeat the prior transfer-validity validator, redesign fit regularization, or adopt a hierarchical fit here. The deliverable is parent-resolved SIMC occupancy, weighted support, and missing-mass shape comparison, with transparent limits.

## 2. Authority and startup hard gate

Authoritative repository: `trottar/lt_analysis`, branch `test`. Authoring-time remote `test` was observed at `ef2bb15517e6b1c271c5189d487e8a71c042e79c`; this is **not** a pinned implementation HEAD. Before any modification, independently establish live remote `test`, local branch/HEAD, cached remote, complete index/worktree and relevant untracked paths. Read tracked root `AGENTS.md`, then `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, and `docs/memory/USER.md`, in that exact order. Next read `docs/memory/MAINTENANCE.md`, `docs/memory/CODEX.md`, CURRENT's baseline investigation and the exact current source owners. Stop `BLOCKED` on material drift, competing pending contracts/edits, missing authority, or unexpected dirty state; do not reset, stash, clean, overwrite, commit, or push.

Source owners to audit before implementation:

- `src/utility/pion_component_shapes.py`: `load_pion_simc_component_shape`, `load_setting_pion_component_shapes`, SIMC source/tree resolution, full-range MM histogram axes, fill/weight/normalization conventions.
- `src/cuts/apply_cuts.py`: `_compute_simc_cut_state`, `apply_simc_sub_cuts`, `apply_simc_cuts` and prerequisite cut configuration.
- `src/main.py`: existing setting-wide component-shape loader call.
- `src/cuts/pion_component_fits.py`: `_resolve_component_scope_hist` and actual parent/template fallback.
- `src/cuts/pion_t_bin_parents.py` and authoritative input/axis configuration: real canonical-t boundaries and parent identity.
- `testing/analyze_pion_transfer_validity.py`: detached prior evidence only; do **not** modify it.

Relevant baseline scientific record: `docs/memory/investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md` and the current October 18 checkpoint. Retain full scientific and provenance limitations.

## 3. Allowed changes — no others

1. **NEW** `testing/audit_pion_simc_parent_support.py`: one explicit, read-only, farm-runnable detached diagnostic/CLI with JSON output. Keep numerical accumulation/report validation importable without ROOT so deterministic checks work off-farm; import ROOT and existing farm source only for the requested event scan. Input paths and settings must be explicit/resolved from existing authoritative configuration, not machine-specific hard-coded paths.
2. **NEW** `testing/test_audit_pion_simc_parent_support.py`: deterministic ROOT-free unit tests using small synthetic events/histograms/axis metadata; no real farm fixtures.
3. **NEW** this single Markdown task contract, `docs/memory/phases/baseline-pion-simc-parent-support-audit-task-contract.md`.
4. `docs/memory/manifest.json` **only** for regeneration by the tracked canonical generator to inventory this new contract. Confirm no pre-existing manifest entry changes.

No `src/` change, `main.py` wiring, existing test change, cut change, SIMC/template/normalization change, fit/refinement change, subtraction, yield, error, or physics configuration change. No broad CURRENT/MEMORY/HANDOFF reconciliation: record any nonblocking memory warning for a later meaningful checkpoint. If source behavior cannot be reproduced without an out-of-allowlist change, stop and report the first specific blocker.

## 4. Required scientific measurement

Operate on the **same existing** `pi_n`, `pi_sidis`, and `pi_delta` SIMC trees and setting configuration used for each selected Center/Left lowe production setting. Check file existence, tree identity, required branch availability and configuration identity before scanning. No fallback to other SIMC generations, nearby kinematics or invented paths. The diagnostic must not mutate trees, analysis files, ROOT objects in production ownership, or production caches.

Reproduce the authoritative full-template event selection exactly: `set_val`/cut-state configuration as needed; `apply_simc_sub_cuts` (acceptance, diamond, t range and t-minimum resolution, without an extra MM cut); the *existing* kaon HGCer hole exclusion when appropriate; finite `iter_weight` with the same default as the loader; MM=`evt.missmass`, |t|=`-evt.t`, and phi from `evt.phipq` in degrees when requested. Separately account for the additional `apply_simc_cuts` population. Do not silently change which population constructs `setting_shape_full`.

Use **source-authoritative** canonical t edges and full-range normal MM axes. Cross-check Center parent edges/identity against the existing seven-stage JSON when it is available as explicit optional evidence; do not require it for other settings. Where Left canonical identity or phi edges cannot be established from authoritative source/configuration, fail that requested comparison as `NOT VERIFIED` rather than inventing bins or borrowing Center's edges. Use the existing `find_2d_bin_indices` convention where phi binning is available. The current component loader constructs binned histograms only when **both** t and phi edges are supplied: the new detached measurement must not mislabel absent cells as populated or introduce those arrays into production.

For each component and setting, retain the exact *unnormalized* setting-wide MM histogram and, separately, per-t-parent histograms plus optional `(t,phi)` children. For **each** population, region and bin report:

- event count; signed `sum_w`; absolute `sum_abs_w`; `sum_w2`; signed effective support `sum_w**2/sum_w2` and absolute effective support `sum_abs_w**2/sum_w2` where defined; cancellation fraction with safe zero handling;
- full-range MM edges, normal-bin contents and Sumw2, **plus** underflow and overflow separately; normalization integral and whether normalization is actually meaningful;
- each component's source-defined neutron, SIDIS and delta anchor-window occupancy/support and the under-Lambda region only when applicable window bounds can be sourced or supplied as explicitly labeled diagnostic metadata;
- selected-event t/phi coverage, out-of-parent and unassigned events, zero-support children and nonfinite/invalid weights, without silently dropping their accounting;
- independent parent-versus-setting and **parent-versus-disjoint-complement** shape comparisons where the latter is supported, with explicit weighted variances and shared-event covariance caveat. No independent-sample p-value for parent versus inclusive setting.

The detached report must demonstrate per-MM-bin and integral bookkeeping closure: setting-wide selected events = assigned parent events + explicitly unassigned/out-of-range events, with signed, absolute, and Sumw2 consistency; child-to-parent closure if child analysis is enabled. If the authoritative t-parent coverage is genuinely partial, report that limit rather than forcing parent totals to equal the setting total. Normalize only after retaining raw moments; never unit-normalize empty/invalid histograms. Compare reconstructed setting-wide template shape against the **actual existing production setting-wide shape** when same-input provenance permits; a mismatch is a hard stop, not a license to replace it.

Report per-setting/component/parent `SOURCE VERIFIED`, `RUNTIME VERIFIED`, or `NOT VERIFIED` only to the extent the actual evidence permits; scientific adequacy is **not** automatically PASS based on occupancy or good closure. Retain fit-data/SIMC distinction, unresolved template uncertainty, and source provenance including source HEAD and input file identity. Use generic settings/parent iteration where feasible; no Center-t2 numerical constants in code.

## 5. Local gates, farm boundary, and handoff

Implement deterministic ROOT-free tests covering histogram axis/flow boundaries, the exact chosen t/phi edge rules, signed weights and random-like negative contributions, Sumw2 and effective counts, zero/low support, disjoint parent/complement bookkeeping, missing branches/provenance, finite-value rejection, and supported/unsupported shape-comparison cases. Require compilation, targeted unit tests, `git diff --check`, a clean allowlist audit, and the canonical manifest and memory-health integrity commands, reporting any pre-existing CURRENT soft-size warning as nonblocking. These are **local source and arithmetic checks**, not farm validation.

Output from Codex: actual live/local identities and full worktree inventory; source trace confirming exact selection and axes; complete diff including full new files; deterministic test receipts; representative synthetic JSON fields; manifest/health results; one **proposed** narrow farm validation invocation with required explicit input roles and expected outputs, **not executed**. Do not create a second contract, phase, or unrelated remediation.

The normal gates remain: ChatGPT reviews **actual diff**, user commits/pushes, ChatGPT verifies published state, **user alone** runs the farm diagnostic if authorized, and ChatGPT independently reviews fresh JSON/evidence. Do not commit, push, or run any farm program as Codex.

**Scientific decision gate:** The output should let us say for each component and parent whether the setting-wide shape is reasonably supported, demonstrably inconsistent, or `NOT VERIFIED`. Do not automatically alter templates or fits, clip high weights, add neutron floors, change signed-event subtraction, or claim accepted production physics. If sufficient evidence is absent, identify the exact missing measurement and stop.
