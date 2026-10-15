# KaonLT — Center/low-ε neutron-anchor attempt-output isolation proof

Status: `ACTIVE` — source-changing isolation-only task contract; no farm execution, no physics acceptance, no owner resumption.

## Authority and starting identity

The source-audited starting remote `test` HEAD is `d860eeafbb889bcc04cfddb1ff6019cfcff50922` (parent `474cbefec2a719c6b19ba06ae55a3c7488b1d28b`). Before any changes, independently verify local `test`, full HEAD, remote `test`, origin URL, index and complete staged/unstaged/untracked state. For read-only live identity, HTTPS `git ls-remote --heads https://github.com/trottar/lt_analysis.git refs/heads/test` is an acceptable fallback if `git ls-remote origin` cannot authenticate by SSH; never alter Git credentials/configuration to resolve a read-only check. Stop if the remote SHA differs or if any task-owned path was unexpectedly changed.

Read tracked root `AGENTS.md` and, in that order, `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, `docs/memory/USER.md`. Then read `docs/memory/MAINTENANCE.md`, `CODEX.md`, `COMMUNICATION.md`, `TOOLS.md`, `docs/memory/decisions/farm-validation-bundle-procedure.md`, `docs/memory/investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md`, and the existing `docs/memory/phases/baseline-pion-center-t2-anchor-farm-owner-task-contract.md`. Keep `CURRENT.md` as sole authoritative NEXT.

The previous owner task stopped correctly before implementation. It found that the normal launcher and `set_SymLinks.sh` can route output through shared volatile storage even from a detached worktree; therefore, the owner/checker/collector are not yet implemented and farm readiness is `BLOCKED`. Existing neutron-anchor capture source (`src/cuts/pion_control_anchor_audit.py`, `src/cuts/pion_t_bin_parents.py`) is `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING` and must remain unchanged. The original Center/low-ε canonical-t2 `B_n=0` cause is NOT VERIFIED. The 67-page setting-wide PDF does not prove the per-t cause.

## Single objective

Implement **only** a reusable, deterministic, source-backed **attempt-output isolation preflight/helper** for the future narrow neutron-anchor farm owner. Demonstrate, first in disposable local fixtures, that the **unchanged normal** `run_Prod_Analysis.sh` would be unable to delete, overwrite, or relink accepted/shared outputs, the ordinary checkout, protected external inputs or installed `ltsep` during a full `Q4p4W2p74` normal (low + high) analysis. If this cannot be proved without touching frozen source or production interfaces, stop `BLOCKED` with the first exact invariant. Do not run the analysis or farm to discover the answer.

This contract is an isolation proof and helper only. **Do not implement or resume the blocked farm owner, checker, collection profile, analysis execution, source promotion or any physics change in this task.** It must be consumable by that later separately reviewed owner, not be a second independent runner.

## Exact allowlist

1. `docs/memory/phases/baseline-pion-center-t2-attempt-output-isolation-task-contract.md` — this one new standalone task contract.
2. `testing/pion_anchor_attempt_output_isolation.py` — one new user-level pure-Python/standard-library helper providing preflight, temporary-attempt-path setup and post-setup verification; **no** analysis invocation, collector invocation, or ROOT import.
3. `testing/test_pion_anchor_attempt_output_isolation.py` — deterministic tests in disposable filesystem and Git fixtures; must not contact farm, network, production paths or real user data.
4. `docs/memory/CURRENT.md` — only necessary task-state/NEXT update, retaining accepted closures and all other blockers.
5. `docs/memory/manifest.json` — regenerate/check only after allowlisted memory changes.

Do **not** edit `run_Prod_Analysis.sh`, `set_SymLinks.sh`, `src/main.py`, `farm_env/print_ltsep_path_fields.py`, installed `ltsep`, `src/cuts/`, `src/utility/`, any model or data, existing owner/checker/collector/profile, `AGENTS.md`, `MEMORY.md`, `CURRENT_HANDOFF.md`, `roadmap/STATUS.md`, existing evidence, or unrelated review bundles. If this allowlist cannot accommodate a safe solution, STOP; request an explicit later contract rather than quietly expanding scope. No broad memory-maintenance cycle.

## Preservation and first source audit

Snapshot every pre-existing unrelated tracked/index/untracked path, including its exact bytes where possible, before touching task files; do not assume counts from prior Codex messages. The delivered contract is an expected untracked addition. The initial manifest will naturally lack that new contract: record this expected baseline state; regenerate/check the manifest as an allowed task step after implementation rather than treating it as an excuse to bypass health checks. Unexpected startup change to any task path or different live HEAD is a hard stop.

Read the full real entry path and resolve each write target:

- `run_Prod_Analysis.sh` calls `git clean -fdx` unless an existing mode changes that branch, then invokes `./set_SymLinks.sh`; the **normal** mode required to reach Center/low-ε also runs low + high and is not `-d` (which is Left/low-ε only). It performs explicit `rm -f` on merged output ROOT targets and writes additional artifacts. Do not use a different runtime mode as a workaround.
- `set_SymLinks.sh` derives `LTANAPATH`, `VOLATILEPATH`, `SIMCPATH`, `OUTPATH` and other paths from `farm_env/print_ltsep_path_fields.py`/`ltsep`. It can create or replace `LTANAPATH/OUTPUT`, `LTANAPATH/log`, `LTANAPATH/src/data`, `LTANAPATH/simc_gfortran` and its `OUTPUTS`, `input`, `worksim` links, and `background_samples/OUTPUTS`/SIMC links; it also creates the shared-volatility-derived background output directory. **Every reachable link, parent and destructive target must be included.** Resolve symlink chains to actual destinations rather than reasoning only from lexical `LTANAPATH`/`OUTPATH` strings.
- Existing reviewed isolation primitives in `testing/run_e8_4_fix5_canonical_five_plot_gate.py` include detached temporary user-owned worktree, copied `ltsep` overlay, source/path probes and ordinary-checkout preservation. They are *not* proof that the normal launcher output remains isolated. Reuse their audited semantics where applicable, but do not change their code or accept their present `OUTPATH` lexical check in place of realpath confinement.
- The active Project `KaonLT_Environment.md` owns concrete farm/workstation roles. The ordinary farm repo, accepted artifact directory and Globus destination must remain untouched. Do not insert private machine-specific paths into tracked code or memory. Derive any **new** attempt path only from an independently verified user-writable, owner-created temporary workspace, with strict parent/ownership/collision checks; do not assume a generic `/tmp` or arbitrary volatile path is allowed.

## Isolation invariant and implementation boundary

A proposed setup is acceptable only if **both** the path probes **and** a disposable Linux simulation prove the full normal-launcher cleanup/link-setup ordering cannot escape attempt ownership. It must not work merely before `git clean -fdx` and then lose its aliases when cleanup runs.

1. Create and identify a unique owner-created attempt directory, separate from the ordinary checkout and every accepted/shared artifact directory; require absent/new identity, exact parent, no unexpected symlinks, no pre-existing files, and ordinary user access. Use the established temporary-worktree parent only when independently proved safe. Never delete or replace a path not created and owned by this attempt.
2. Model the actual sequence: clean disposable *detached Git worktree* -> prepare `ltsep` overlay/path identities -> establish any candidate disposable aliases -> **execute the same `git clean -fdx` semantics in a wholly synthetic local Git fixture** -> execute or faithfully model only `set_SymLinks.sh` setup under synthetic fake external directories -> resolve all resulting output, log, data, background and SIMC links. No invocation of the real normal launcher in this task.
3. Require actual resolved `OUTPATH`, all produced `OUTPUT/Analysis/...` and log write targets, `src/data`, background output and any other write-capable path to be within the owner-created attempt path. A mere symlink named inside the worktree is insufficient. Shared `VOLATILEPATH` must remain available as an external *input* if needed, without granting it output writes.
4. Preserve read-only external inputs byte-for-byte: SIMC, background samples, replay/skim sources, installed `ltsep`, environmental config and ordinary accepted outputs. If `set_SymLinks.sh` would create/replace any link or directory inside an external input tree (`SIMCPATH`, background SIMC, shared volatile), reject the run setup. Do not synthesize or modify missing external inputs, or silently repoint scientific input trees to the attempt path.
5. Ensure the standard launcher cannot remove the candidate output aliases before it runs `set_SymLinks.sh`. Any proposed index/alias strategy is allowed **only in a disposable analysis worktree**, never in the ordinary Git index/checkout, and must be proved under `git clean -fdx`; if index state is temporarily changed, explicitly attest its worktree-only scope and inspect reset/cleanup without touching the ordinary checkout. Do not presume a specific workaround succeeds.
6. Demonstrate that repeated setup, stale attempt paths, unexpected symlink types, broken/hostile links, collisions, permission failures and cleanup interruption fail closed. If any preflight is uncertain, return a structured `BLOCKED` result and make no mutations beyond own new temp files. Cleanup must be bounded to exact owner-created paths and never recurse through external symlinks. Do not use `sudo`, `su`, chmod of external directories, or user-facing `exit`/`logout`.
7. The helper must return a JSON-safe *plan/receipt* containing source HEAD input, path realpaths, link targets, identities, verified writable/read-only roles, cleanup safety, untouched external role fingerprints/links, a list of every anticipated write surface, and an explicit `passed` flag and first failed invariant. It may not assert actual farm execution or runtime acceptance. The later owner will independently re-probe right before execution and after cleanup.

Do not alter `VOLATILEPATH`, `SIMCPATH` or other global input-role values in a copied `ltsep` configuration merely to make output paths appear isolated. If a solution necessarily changes them, explicitly prove identical read-only input source identity at each actual producer/consumer and that no shared write remains, otherwise STOP. Production source and physics remain frozen.

## Deterministic validation requirements

Tests must include: a disposable Git repository with `git clean -fdx` deleting an untracked output symlink (negative control); an exact source-backed safe attempt candidate only if found; output/log/data alias realpath checks; normal clean-then-link ordering; existing/nonmatching and broken links; `background_samples/OUTPUTS`, external SIMC `OUTPUTS/input/worksim` and shared-directory write hazards; archive/analysis path collision rejection; external sentinel hashes unchanged; symlink-escape and symlink-cycle rejection; interrupted setup/cleanup, permission failures, and reentry; unchanged ordinary repo index/head/status; all other settings and science untouched. A test of a workaround must prove it survives real `git clean -fdx` on the fixture and remains owner-confined after actual `set_SymLinks.sh` behavior is modeled or executed **only in a synthetic fixture**. If `set_SymLinks.sh` is invoked in tests, use a fully controlled sandbox with fake `ltsep` path resolver so no real external locations are reachable. Tests must never rely on JLab ROOT or actual farm paths.

Run `python -m py_compile` on any new Python; run the complete new test module and relevant source-isolation regression tests; `git diff --check`; manifest regenerate and check (`tools/update_memory_manifest.py`); ordinary memory health (`tools/check_memory_health.py`) and bootstrap (`tools/memory_bootstrap.py`). Report exact interpreter, tests/results, changed paths, source/manifest state, CURRENT/MEMORY/HANDOFF byte counts, nonblocking CURRENT soft-size warning, and byte-for-byte preservation of unrelated paths.

Hard stop at any failed scientific, identity, external-write, source/manifest, symlink, cleanup or test gate. Do not turn partial-failure results into farm-readiness. This task authorizes **no** farm command, analysis run, owner/checker/collector resumption, production correction, changes to accepted closures, or Method-A/B promotion.

## Codex handoff and user publication

Codex returns **one uniquely named complete actual diff** (tracked diffs plus full new Markdown/Python files, not just summaries), all checks, preservation receipts, and the first failed invariant if BLOCKED. Codex does not stage, commit, push or run the farm. ChatGPT independently reviews actual source/diff before publication.

After ChatGPT explicit diff acceptance, **the user** performs the following exact scoped publication steps (unexecuted here):

```bash
git add -- docs/memory/phases/baseline-pion-center-t2-attempt-output-isolation-task-contract.md testing/pion_anchor_attempt_output_isolation.py testing/test_pion_anchor_attempt_output_isolation.py docs/memory/CURRENT.md docs/memory/manifest.json
git commit -m "Prove isolated Center neutron-anchor attempt outputs"
git push origin test
```

If Codex must stop before creating every listed new source file, it must return the actual changed-file subset and a correspondingly narrower explicit user-controlled Git handoff **only after** independent review. No `git add -A`. After publication, ChatGPT verifies pushed remote HEAD/paths/manifest/CURRENT, then resumes only the previously blocked neutron-anchor owner in a **separate** authorized task. Before any farm execution, review the resulting tracked owner/checker/collection chain and issue an exact single `tcsh`-compatible user-level command **without `exit`, `logout`, `sudo`, shell replacement or implicit session termination**. A verified source-level isolation helper is necessary but not sufficient for farm readiness.
