# KaonLT — October 14 Baseline Pion-Subtraction Scientific Memory Re-anchor

## Task class and starting identity

Memory-only reconciliation/closure-of-handoff task, **not** a scientific-source or farm task.
The remotely observed `test` HEAD when this contract was authored was
`d66b215f30d9e43875bdfd6504c4392b49261122` (2026-10-14).

Before editing, establish the actual local `test` HEAD, complete worktree status (tracked, staged, untracked), and remote identity. Read tracked root `AGENTS.md`, then in order `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, and `docs/memory/USER.md`; then read `docs/memory/MAINTENANCE.md`, `docs/memory/CODEX.md`, relevant current roadmap/evidence records, and this contract. If branch/HEAD differs or the local worktree has unrelated unexpected changes, **stop BLOCKED** and report the mismatch. Do not reset, clean, stash, checkout, or overwrite. An externally delivered copy of this contract is permitted pre-existing state, together with only the user-authorized startup exception below.

### User-authorized startup exception - 2026-10-14

The following seven untracked review bundles are expected pre-existing state.
This exception does not expand the implementation allowlist and authorizes no
modification, staging, deletion or commit of these bundles. Record their startup
SHA-256 values below and verify byte-for-byte preservation at completion.

| Exact filename | Startup SHA-256 |
| --- | --- |
| `kaonlt_review.diff` | `d2e045b40b8f209db6486bdbb3391c0cfd99636ed97d4282500a9ccd2044c7ac` |
| `kaonlt_review_e8_current_f1_coherent_materialization_2026_10_11.diff` | `71915e12e494cc10f262b88cced581e30d7e99290906d4b8b7ade6a0f5949474` |
| `kaonlt_review_e8_current_snapshot_checkpoint_2026_10_12.diff` | `c2d1cb53914aaa600471eea3db686c7bdd976956bc501387591b7b4c68f535f5` |
| `kaonlt_review_e8_four_way_BLOCKED.diff` | `9b672ad0448c90681349813c6d47b56670568bc5203ac0c02477923631c92b33` |
| `kaonlt_review_e8_four_way_final_source_stop.diff` | `d494144419becba8ad6508264af577e054a0b6d4047167dc8b65aa440a7e4e99` |
| `kaonlt_review_e8_four_way_repaired.diff` | `9f86da18c49cca113be0f9031a5a04a8046caea813f85166dae84841b9a7985c` |
| `kaonlt_review_e8_reweighted_contribution_debug_BLOCKED_2026_10_11.diff` | `494284ab8c7e40bb5c6c0b634392a7c79e2577d3c1d9af47efb9751867ef3789` |

The existing tracked `.M` entries below are also expected and must remain
untouched; HEAD/index/worktree Git content hashes were verified identical:

- `docs/memory/LEARNINGS.md`
- `docs/memory/evidence/e8-4-canonical-five-fresh-f1-lineage-materialization-2026-10-05.md`
- `docs/memory/phases/e8-4-canonical-five-fresh-f1-candidate-lineage-refresh-task-contract.md`

All other startup invariants and hard stops remain in force.

## Objective

Record the user-directed scientific priority reached in the October 13–14 prior chat, which could not be checkpointed because that chat ran out of capacity. The current tracked `CURRENT.md` still points to paired E.8 diagnostics, whereas the discussion advanced to a **read-only audit of the existing production baseline pion-component subtraction**. Recover the correct single NEXT and dependency ordering for a future handoff without changing any scientific calculation or inventing acceptance.

This is a **narrow memory update**, not a new phase or a rewrite of the approved F/E program. Keep accepted closures, their scopes, and unrelated blockers intact. The October 14 investigation takes priority **before any attempt to interpret or promote Method-A impact through final E.8/F.6.4**; the independently blocked canonical-five F.4 and SIMC-unit issues stay blocked.

## Allowlist

1. `docs/memory/phases/kaonlt-baseline-pion-memory-reanchor-2026-10-14.md` — this single contract (already supplied; preserve its content).
2. `docs/memory/CURRENT.md` — one active objective and one exact scientific NEXT, with concise transition note.
3. `docs/memory/roadmap/STATUS.md` — a narrowly scoped dependency/prerequisite annotation for baseline scientific correctness; do not rewrite phase history.
4. `docs/memory/investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md` — one new bounded investigation record containing the evidence/chronology and limits below.
5. `docs/memory/manifest.json` — regenerate/check using the tracked memory tooling.

Do **not** change `MEMORY.md`, `handoffs/CURRENT_HANDOFF.md`, `USER.md`, `AGENTS.md`, `MAINTENANCE.md`, `CODEX.md`, `LEARNINGS.md`, implementation source, tests, plots, data, validation profiles, owner scripts, or any other path. The handoff is ordinary `CURRENT.md` state, not an exceptional handoff override. If another file is truly required for a hard consistency check, stop and report it rather than widening scope.

## Historical evidence to record accurately

Evidence classes must be explicit; the prior chat, user-provided PDFs, and previous assistant's source audit are **historical evidence**, not a newly executed farm gate in this memory task.

**Scientific setting:** Q4p4W2p74 / Center / low epsilon / canonical `t2`, approximately `0.567 <= |t| <= 0.733 GeV^2`; adjacent Center `t1` is also concerning. This is the **existing production baseline**. Do not switch to Method A or Method B.

**Reported archived August 26 parent diagnostic:**

- Staged pion-control `B_n = 0`; joint-refined `B_n = 0.0002021`.
- Kaon-side `A_n = 0.0369527`; neutron amplitude ratio approximately `182.8`, maximum transfer weight approximately `183`.
- Pion-control staged (`B_n`, `B_sidis`, `B_delta`) approximately `(0, 2.608, 1.692)` and joint-refined approximately `(0.000202, 1.689, 2.446)`.
- Final kaon-side (`A_n`, `A_sidis`, `A_delta`) approximately `(0.036953, 0.215900, 0)`; protected kaon-side delta fit was rejected/zeroed.
- Fitted background approximately `0.2529`; event-weighted applied background approximately `6.556`; event/model integral ratio approximately `25.9`.
- Adjacent Center `t1` event/model integral ratio approximately `4.32`, and `t3` approximately `1.84` according to the historic diagnostic.
- The concurrent setting-wide fit was much better behaved, with integral ratio around `1.01`; a separate archived 67-page Q4p4 Center setting-wide reference had maximum weight `0.0819`, fitted background `0.5967`, event-weighted `0.6010`. Do not treat differing scopes/PDFs as an identical-input regression test.

**Verified source-level chronology from the read-only Git comparison (recheck if needed):**

- `f924a95` (August 19) defaulted `pion_subtraction_scope` to `t_bin`, making independently fitted canonical-t solutions authoritative; setting-wide fit remained diagnostic.
- `ee9e292` (August 20) moved authoritative parent-fit production to `src/cuts/pion_t_bin_parents.py`.
- Per-t fitting existed previously in averaging code. Thus authority changed; it is **not proven** that earlier per-t final yields were correct, or that the algorithm first became unstable then.
- The original August 26 anomalous result preceded the explicit August 27 RF/no-RF tree switch. That later switch **cannot be the initiating cause** of the August 26 pathology; do not keep it as the lead hypothesis.
- The original staged single-template solver, staged pass, joint refinement, regularization construction, and scope-template resolver were reported unchanged across the transition. Neither the protected kaon-side delta failure nor a final kaon-side refit is established as the initiating cause of `B_n=0`.

**Source findings and bounded hypotheses:**

- The neutron component is fitted first from a **signed** pion-control histogram in the neutron anchor around `0.88–0.94 GeV`, using a nonnegative least-squares estimator; zero can result from a nonpositive weighted numerator or an unsuccessful/no-support solve. Which occurred in the original Center `t2` fit remains `NOT VERIFIED`.
- Joint refinement is regularized about the staged zero with a width floor of `0.001` (not a fitted-amplitude lower bound), discouraging recovery of `B_n` once initialized at zero.
- Baseline fitting/application quality gates can accept fits with weak decomposition, large ratios and model-versus-event discrepancies; absence of configured max quality limits is a source-level gap, **not** proof that any proposed alternate model is physically correct.
- `src/utility/pion_component_shapes.py` individually unit-normalizes `(t,phi)` SIMC templates before `src/cuts/pion_component_fits.py` sums and renormalizes them into a parent-t shape. Thus parent shapes may equally weight populated phi slices instead of preserving original event-weighted mixture; this is a **testable input-shape hypothesis**, not an accepted root-cause finding.
- Do not declare signed cancellation, wrong normalization, alignment, SIMC template construction, or a specific commit to be the root cause without matching actual neutron-anchor data/template/solver evidence.

**Paired E.8 continuation:** The previous conversation reported later paired Left/lowe + Center/lowe progress after tracked CURRENT was last published. Preserve any independently accepted paired evidence/closure already recorded. Where only previous-chat descriptions are available, mark them `MEMORY ONLY` / attributed prior review; do not invent a new farm receipt, hashes, 16-page visual acceptance, or silently upgrade a status. The old paired-farm NEXT is **superseded as the user's immediate scientific priority** but is not declared failed. Prior Left-only and canonical-five closures remain scoped exactly as recorded.

## Exact desired memory behavior

**Before:** `CURRENT.md` says paired diagnostic farm validation is the sole NEXT; October 14 baseline investigation is absent.

**After:** `CURRENT.md` describes one active priority: a narrow, read-only **baseline pion-control neutron-anchor regression investigation** for Q4p4W2p74 Center/lowe `t2` (compare Center `t1` where relevant), anchored to the first staged `B_n=0`. State the crosswalk to older paired E.8 separately, with no contradictory active NEXT. Keep an unambiguous distinction between (i) accepted E.8/F.6 status, (ii) baseline scientific anomaly, (iii) candidate mechanisms, and (iv) blocked/deferrred unrelated work.

Write **one push-stable NEXT** in `CURRENT.md`: determine from existing applicable records/source-backed diagnostics whether the original Center/lowe `t2` neutron-anchor nonnegative solve reached zero because of the signed weighted numerator or because of absent/invalid template/fit-bin support, with consistent cut, bin, coordinate, signed-source and SIMC-template identities. Record prompt/random/dummy/dummy-random anchor contributions, `fit_bin_indices`, `fit_y`, `fit_sigma`, `fit_template`, `weighted_numerator`, `weighted_denominator`, `success/message`, effective joint prior and populated-phi SIMC inputs; first audit existing artifacts, then, only if missing, define a **separate narrow diagnostic-only source contract**. No new farm command or production correction from this memory task.

The roadmap gets one concise note making **scientific baseline fit/transfer validity an explicit prerequisite** to further interpreted full E.8 production-impact claims and any F.6.4 promotion decision. This **does not** recategorize F.1/F.4, reopen closed statuses, create new Phase G, or imply that prior baseline-five owner/identity checks were invalid; they proved their stated limited gates, not all component-decomposition physics.

The new investigation record owns the detailed evidence and open questions; `CURRENT.md` links to it briefly. Avoid copying the investigation into multiple memory files. Do not attribute newly inspected ROOT data, rerun fit numbers, locally validated PDF pages, or fresh farm acceptance to this task.

## Invariants / frozen behavior

No scientific or production code modifications; preserve all random/dummy/slow-proton/pion handling, HGCer roles, SIMC inputs, cuts, fit models, priors, templates, numerical ratios, event weights, errors, yields, L/T, cross sections, `no_empirical_residual`, and provenance. Method A stays detached/non-production; Method B independent and numerically excluded. Canonical-five F.4 provenance, absolute-SIMC unit proof, final E.8 and F.6.4 remain `BLOCKED` where recorded. The separate historical F.1/F.4 forensics remain `DEFERRED`; baseline pion-control investigation is a different, user-authorized question. Keep `E.8.2 Left/lowe` and five-setting baseline numerical gates `CLOSED / RUNTIME VALIDATED` within their exact accepted scopes, without interpreting them as proof of Center `t2` fit quality.

## Local positive, negative, and regression checks

- Positive: `CURRENT.md` has exactly one active scientific objective and one NEXT (not a commit/push or farm command); the new investigation record matches the source/history and has explicit evidence-class limits.
- Positive: roadmap annotation states prerequisite ordering without inventing phases or changing accepted statuses.
- Negative: search the diff for invented `RUNTIME VERIFIED` claims, unsupported completed-source/farm status, production algorithm suggestions presented as decisions, changed numeric source/config or changes outside allowlist.
- Regression: compare preserved Left/lowe, E.8.2 baseline-five, E.8.3, F.6.3, F.6.4, Method A/B, canonical-five and SIMC blocker text/meaning with their previous source.
- Confirm only one standalone task contract and no second prompt `.md` were added.
- Discover local Python using `docs/memory/TOOLS.md`. Regenerate/check manifest with the tracked tooling; run `<PYTHON> -B tools/check_memory_health.py --root .` (ordinary non-strict gate). Report `CURRENT`, `MEMORY`, `CURRENT_HANDOFF` byte counts, manifest result, exact warnings, and classify each warning as blocking/nonblocking. Stop on any hard failure or blocking warning. `git diff --check` and scoped/full-status review must pass.

## Diff audit and acceptance

No farm, ROOT, plot, or integration validation is part of this task. Codex must report full `git status --short --untracked-files=all`, modified file list, before/after source-state meaning, validation output, and a complete reviewable diff including full additions for this new contract and investigation file; do not stage merely for diff review. ChatGPT independently reviews actual source and actual diff. **Codex must not commit or push.** Only after ChatGPT accepts the diff will the user run scoped Git publication; ChatGPT then performs pushed-state synchronization and returns to the baseline first-failure gate.

**Hard stop:** branch/HEAD mismatch, unexpected dirty worktree, missing authoritative file, unsupported status upgrade, substantive data/science source change, added path outside allowlist, failed manifest or memory health gate, or inability to keep a sole coherent NEXT. Report `BLOCKED` and the first invariant, with no speculative fallback.
