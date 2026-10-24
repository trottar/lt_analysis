# KaonLT — Seventh pre-child failure: durable incident and review-discipline correction

**Task class:** Repository-memory correction only, under the existing root `AGENTS.md` workflow. No new analysis phase, runtime gate, source implementation, or farm command.

**Source anchor (reverify):** GitHub `test` HEAD `b05ad2a1d9af241527a7d6ce90196eb591188d19` on 2026-10-23. The farm checkout and local WSL worktree/index are **NOT VERIFIED**. Before editing, inspect live/local branch, remote HEAD, staged/unstaged/untracked files, and any pending overlaps. Read repository-root `AGENTS.md`, then `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, and `docs/memory/USER.md` in that exact order; next read `docs/memory/MAINTENANCE.md`, existing October 22–23 matched-cohort failure investigation, and source at the failed call sites.

## Facts to record without upgrading evidence

The user operated the published full matched-cohort owner on ifarm2402 with `--source-commit b05ad2a1d9af241527a7d6ce90196eb591188d19`. It stopped with:

`Matched-cohort gate failed: caller_path_invalid:src/setup/findEffectiveCharge.py`

The printed status location was:

`/group/c-kaonlt/USERS/trottar/lt_analysis/OUTPUT/Analysis/KaonLT/KaonLT_MatchedCohort_Q4p4W2p74_20261023-gate-status.json`

The user reported that the attempt consumed hours. The gate-status file itself, its stage timestamps, actual filesystem resolution, and whether the output/checkout preservation checks passed were **NOT VERIFIED** here. Do not invent those values. **RUNTIME VERIFIED:** the supplied CLI failure and status path; there is no supplied successful child completion, matched cohort ZIP, or new scientific result.

**SOURCE VERIFIED:** In `testing/run_pion_local_t_matched_cohort_gate.py`, `probe_resolver()` rejects `not path.exists()` or `resolver.real(path) != path` under `caller_path_invalid`; `resolver.real()` resolves absolute paths and symlinks, including farm `/group` versus `/u/group` aliases. Whether this was the precise trigger is **NOT VERIFIED**, not a proven alias diagnosis. The full owner calls recursive `shared_state()` inventory/hashing of shared analysis output before `resolver_preflight()` reaches this check; `shared_state()` can be reevaluated during preservation/cleanup. This ordering was overlooked in source review; exact runtime attribution of the hours is **NOT VERIFIED**.

This is the **seventh** failed pre-child owner/smoke attempt in the existing incident series. The user correctly identified repeated assistant review/readiness failures and displaced science, time, and cost. Codex executed tasks and supplied local tests, but ChatGPT approved their scope, missed path/topology assumptions, failed to trace the entire actual runtime call chain and expensive ordering, and improperly promoted local/synthetic passes into farm readiness. Record assistant accountability, not fault in the user's command. The accepted scientific checkpoint remains October 20. **No fresh Center/Left low-epsilon canonical t1/t2/t3 same-run A/B + signed-event evidence, no seven-strength C replay, no six-parent Lambda-window MM/variance result.**

The latest further workflow violation was an assistant-produced direct `git apply` patch and manual repository mutation instructions. That **was not** the KaonLT process. It was explicitly withdrawn. Do not incorporate or apply that patch. The normal tracked flow is audit -> **one standalone Markdown task contract** -> short Codex launch -> Codex implementation/deterministic checks -> independent actual-diff review -> user commit/push -> pushed-state synchronization -> farm only if proven ready. No assistant direct patching and no unreviewed user changes.

## Memory-only changes

1. `docs/memory/CURRENT.md`: replace stale sixth-failure/current-gate prose with a concise seventh-failure `BLOCKED` summary. Preserve the scientific objective and latest accepted science, frozen roles and source identity. State that a source/path and pre-analysis cost audit, **not a new farm attempt**, is the first gate. Keep CURRENT under the tracked 16-KiB hard limit; compress obsolete material instead of appending.
2. `docs/memory/investigations/kaonlt-matched-cohort-farm-readiness-failure-2026-10-22.md`: append one tightly factual seventh-failure entry, source evidence, unknowns, repeated review accountability, and withdrawal of the direct-patch instruction. Preserve six historical entries unchanged.
3. `docs/memory/MEMORY.md`: add one compact durable cross-project review rule: prefer mature working source and the shortest science-producing path; do not approve farm readiness merely because a diff and synthetic suites pass; trace CLI -> preflight -> child -> validation/collection with actual farm path aliases, symlink/cwd/import geometry, cheap-to-expensive check ordering, and expected work, before authorizing farm; expose and stop at unverified hard invariants, not repeated speculative infrastructure fixes. Explicitly preserve the existing user-only farm/publication boundary and forbid direct assistant patch/application workflows.
4. `docs/memory/manifest.json`: regenerate with the existing tracked memory tooling after changes.
5. Add this single contract under `docs/memory/phases/kaonlt-seventh-farm-failure-memory-task-contract-2026-10-23.md`.

**No other file changes.** Do not touch `testing/`, `src/`, frozen helpers, analysis launchers, farm files, or previous scientific artifacts. No revised farm commands. Do not start a second support phase.

## Checks and handoff

Run the tracked memory manifest write and check, memory health and bootstrap, Markdown/file/size checks and `git diff --check`; report warnings separately from blockers. Review the full staged/unstaged diff for exactly the allowlisted paths. Produce one actual-diff review bundle containing the full new contract, changed records, and test receipt. **Do not commit, push, or run the farm.** The user publishes only after ChatGPT independent actual-diff review.

**Authoritative scientific NEXT remains:** fresh same-cohort Center/Left low-epsilon A/B and signed events, then existing C and six-parent MM/variance comparison. The seventh failed owner gate provides no accepted science; this memory correction does not authorize a new farm run.
