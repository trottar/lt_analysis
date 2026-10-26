# KaonLT A/B owner — volatile-only disposable runtime placement repair

## Authority and starting point

- Authoritative repository: `https://github.com/trottar/lt_analysis/tree/test`.
- Reviewed pushed starting SHA on October 25, 2026: `942fab7ed06ba53bff18ce053093a89a2cd55896`.
- First verify actual local `test`/`origin/test` HEAD, status/index and tracked root `AGENTS.md`. Read, in order, `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, `docs/memory/USER.md`, followed by `docs/memory/MAINTENANCE.md` and the October 25 A/B owner preservation contract. Preserve unrelated local changes. Stop at unexpected HEAD, an overlapping edit, or an unproven isolation invariant.
- This is a **single narrow runtime-placement correction** to the already source-reviewed A/B owner, not a new method, science phase, acquisition mechanism or general worktree redesign.

## Observed defect and supplied farm evidence

**SOURCE VERIFIED:** `testing/run_pion_t_template_comparison_gate.py` calls `isolation.owned_worktree(repo, sha, "analysis")`. The shared helper in `testing/run_e8_4_fix5_canonical_five_plot_gate.py` creates its directory with `dir=repo.resolve().parent`. `private_output()` creates `parent / "ab-output"`; `stage_copies()` writes to `ab-output/Analysis/KaonLT`. `isolation.prepare_runtime_overlay()` independently calls the frozen shared `cleanup_path`, which also requires a repository-sibling worktree. Merely changing `TMPDIR`, `--outdir`, or one symlink cannot fix this.

**RUNTIME VERIFIED on ifarm2402:** `pwd -P` is `/u/group/c-kaonlt/USERS/trottar/lt_analysis`; the intended installed Python is `/u/group/c-kaonlt/USERS/trottar/replay_lt_env/bin/python3`, installed `ltsep` under `lib64/python3.9/site-packages`. Owner preflight passed, 32 owner tests passed, the memory manifest passed, and memory health had only its known nonblocking CURRENT size warning. `copy_plan` found **19 input files, 8.85 GiB**, including Center/Left/Right high-epsilon archived SIMC. `df` showed `/group` **11 GiB** free and `/volatile` **4.4 PiB** filesystem-wide free (not proof of user quota). The ordinary analysis `OUTPUT` path is already symlink-backed onto volatile; preserve it.

**Hard invariant:** Never place copied or newly written analysis data in `/group`. The actual owner-private analysis worktree, its copied ROOT/hist inputs, derived ROOT data, runtime caches, copied `ltsep` tree, intermediate files and private science outputs must all be beneath a fresh user-owned 0700 directory under the **existing resolver-provided `VOLATILEPATH`**. The code-only independent clean source worktree may retain its existing location if it provably writes no analysis data. Do not invent a replacement path or hard-code an unverified scratch mount.

## Exact repair

1. **Owner-local placement and confinement only.** In `testing/run_pion_t_template_comparison_gate.py`, derive the actual user volatile base from the unchanged installed baseline `ltsep` resolver (`isolation.path_fields(repo, repo)` / the established baseline field `VOLATILEPATH`). Before creating anything, require the configured base is absolute, exists, writable, has correct containment under the real `/volatile` mount after resolution, and is not an alias for `/group` or the ordinary repository. Use canonical paths and filesystem identity (`stat`/mount device as appropriate), not string-prefix matching alone. Reject surprising symlinks/escape and ambiguous base configurations. Do not modify installed `ltsep`, `PATH_TO_DIR`, `VOLATILEPATH`, `OUTPUT`, or shared scientific artifacts.

2. Create a **unique, private, mode-0700 disposable analysis parent** directly beneath that verified base. Place the detached `analysis` Git worktree in `parent/analysis`; retain the already reviewed `parent/ab-output`, `parent/runtime-python`, `parent/ab-tmp`, `parent/ab-cache`, `parent/ab-mpl`, `parent/ab-pycache` design. Keep every writable runtime path confined to this volatile parent. Do not reuse the repo-sibling `isolation.owned_worktree` for the analysis role: its explicit placement contract is incompatible. Implement a small owner-local context manager for the analysis worktree, with pre- and post-cleanup containment, ownership, identity and path-bounded checks. Cleanup may touch **only paths created by that owner in the verified disposable parent**. No `git reset`, `git clean`, stash, ordinary-checkout deletion, global prune, or symlink target removal.

3. **Account for the second incompatible guard.** `isolation.prepare_runtime_overlay()` calls the repo-sibling `isolation.cleanup_path`, so simply swapping the analysis-worktree context manager would still fail before execution. Make the narrowest **owner-local** overlay adapter that preserves the existing copied-installed-`ltsep` behavior *byte-for-byte except for worktree-location validation*: baseline import identity, unique matching PATH_TO_DIR config selection through the installed SetPath logic, only copied `LTANAPATH` substitution, untouched stable roles including VOLATILEPATH, original installed-package hashes, copied-Python import checks and `verify_ltsep_preservation`. Reuse existing pure utilities; do not monkeypatch the shared helper or broaden its safety gate. If a safe owner-local implementation cannot be made while preserving this contract, stop rather than modifying unrelated frozen owners.

4. Update the A/B `execute_gate` analysis-worktree creation and overlay call only. Keep the independent `"source"` collector worktree, profile, SHA/origin check, `source_preflight`, `private_output`, `private_callers` (all actual default/Plot/cwd observations), external symlink no-repair checks, mandatory Center/Left lowe inputs, optional SIMC availability, exact Diamond namespace and tie selection, 19-file **source-derived** `copy_plan`, byte-verifying independent copies, ordinary snapshot, accepted-artifact identities, one `-d` child, four A/B artifact checks, independent collection, ZIP/provenance/freshness/CRC, cleanup-before-delivery and failure receipts. Keep all scientific code and mathematics unchanged.

5. **Fail closed before large copying or analysis** if the actual resolved private output, private source data destination, temporary/cache destination, or copied package is on `/group`, outside verified `VOLATILEPATH`, overlaps shared original OUTPUT, points through a mutable symlink, has incorrect UID/mode, or cannot be safely cleaned up. Report the verified actual lexical/canonical private parent, filesystem device, source-resolver VOLATILEPATH, intended staging byte total and output locations in the existing gate-status/run-summary diagnostics. Do not introduce recursive output hashing, bulk copying from unrelated directories, long acquisition wrapper, or any test that runs farm/ROOT.

6. Add narrowly focused regression tests in `testing/test_run_pion_t_template_comparison_gate.py`. In temporary synthetic directories, model a repo under `group`, an ordinary OUTPUT symlink to `volatile/OUTPUT`, and a separately validated volatile runtime root; assert the actual *owner* analysis worktree, `ab-output`, copied ROOT fixtures, cache/temp and copied `ltsep` all stay within the volatile root. Reject repo-sibling private output even when enough space exists, swapped/escaped symlinks, a VOLATILEPATH alias to group, invalid volatile root, and unsafe cleanup. Include success and injected-failure cleanup cases, no accepted-output mutation, untouched ordinary inputs, detached science source, exactly one child in success, no child on preflight failure. Preserve all existing 32 owner tests and other relevant frozen regression tests.

## Strict changed-path allowlist

- `testing/run_pion_t_template_comparison_gate.py`
- `testing/test_run_pion_t_template_comparison_gate.py`
- `docs/memory/phases/kaonlt-ab-owner-volatile-runtime-placement-2026-10-25-task-contract.md` (this contract)
- `docs/memory/CURRENT.md` (minimum accurate blocker/NEXT wording only)
- `docs/memory/manifest.json` (generated)

Do **not** modify `testing/run_e8_4_fix5_canonical_five_plot_gate.py`, `testing/pion_t_template_comparison_gate_profile.json`, production `src/`, launcher, `farm_env/`, unrelated testing, or other phases without a new explicitly reviewed narrow justification. No manual changes on the JLab farm.

## Deterministic completion / review gates

- Run the 32+ owner tests including new Linux filesystem and symlink topology fixtures, existing other directly relevant owner/collector tests and `testing.test_pion_t_template_comparison` pure-Python suite, without ROOT/farm.
- Verify py_compile of touched Python files, `git diff --check`, memory manifest regeneration/check and ordinary memory-health check; record nonblocking warnings without an extra maintenance phase.
- Review actual diff versus the exact starting SHA. List every edited path, explain that all production/scientific source bytes and the frozen shared helper are unchanged, and show a concrete pre-child rejection of `/group` private data placement.
- Present the complete actual diff and new contract file for **independent ChatGPT review**. Codex must not commit, push or run the farm; user alone publishes reviewed files and runs later farm validation. The previous successful farm preflight and 19-input inventory remain useful historical evidence but do not upgrade the corrected source to RUNTIME VERIFIED.

## One substantive NEXT

After independent actual-diff review and user-pushed-state synchronization, execute **one** owner-controlled private A/B farm validation on the correctly configured volatile runtime, inspect the actual gate-status/log/ZIP and four A/B JSON/PDF objects, then return to the seven-strength C scientific comparison. Do not add another unrelated acquisition or infrastructure phase.


## Local implementation outcome - October 25, 2026

`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`. Independent ChatGPT
actual-diff review is NEXT. Starting branch, HEAD and local origin/test matched
942fab7ed06ba53bff18ce053093a89a2cd55896; staged content was empty.
Pre-existing documentation edits were preserved. Startup manifest drift was
exactly this new contract, with no changed or missing existing entries.

SOURCE VERIFIED: the owner obtains baseline VOLATILEPATH through the unchanged
installed resolver. It rejects relative/missing/unwritable/aliased bases,
non-mount or wrong-device roots, group/repository overlap and invalid UID.
A unique 0700 parent beneath that verified base owns analysis, ab-output,
runtime-python and the four cache/temp roots. Directory canonical paths,
UID/mode, device/inode and base/parent identities are checked before staging,
before/after the child and before bounded cleanup. Swapped/escaped directories
stop cleanup; symlink-safe removal is required before worktree creation.
No global prune, ordinary-checkout cleanup or symlink-target deletion is used.
Failure receipts retain placement observations; success summary/status include
lexical/canonical locations, device, resolver base and intended copy bytes.

The owner-local overlay adapter preserves the frozen helper's configuration
selection, SetPath lookup, installed/copy import identities, copied config
bytes except LTANAPATH, stable external roles and preservation checks. An AST
comparison confirms identical overlay algorithm after excluding only the
location guard, additional confinement checks and private directory mode.
Only execute_gate changed among pre-existing owner definitions; copy_plan,
stage_copies, private callers, accepted artifacts and scientific/payload/ZIP
validators are unchanged. The independent code-source worktrees still use the
frozen helper and collect only after private artifacts are verified/alive.

### Deterministic evidence

Local WSL Ubuntu / Linux Python 3.12.3:

```text
python3 -B -m unittest testing.test_run_pion_t_template_comparison_gate testing.test_run_e8_2_left_lowe_scientific_audit_gate testing.test_run_prod_analysis_debug_left_low testing.test_collect_pion_hgcer_validation_bundle -q
```

110 passed: 40 owner (all original 32 retained), 27 debug-owner, 15 launcher,
28 collector. No skips. Fixtures use actual Linux directories, independent
copies, symlinks and real harmless installed/copied-package subprocess probes;
Git worktree operations and analysis children are mocked. The actual owner
context/layout checks run in both success and injected-failure flows. Success
has exactly one child, independent collection while private artifacts exist,
and cleanup before delivery. Group placement is concretely rejected with
volatile_base_invalid before copying or child execution. Tests also reject a
shared group device even with ample space, repo-sibling output, changed
base/worktree identity, swapped parent/cache/output links, wrong UID/mode and
unsafe cleanup. Ordinary inputs/accepted outputs are preserved. The synthetic
copy scope remains 21 files / 1017 bytes; the supplied 19 files / 8.85 GiB farm
inventory is historical, not remeasured or replaced by fixtures.

Windows Python 3.12.10:

```text
python -B -m unittest testing.test_pion_t_template_comparison -q
python -B -m unittest testing.test_run_e8_4_fix5_canonical_five_plot_gate -q
```

23 template tests passed. The additional frozen-helper suite ran 32 tests:
31 passed; one could not create its Windows symlink fixture (WinError 1314),
including on an unsandboxed retry. Its exact unchanged test method and two
fixture-helper ASTs were then executed in WSL/Linux: all three symlink/nonfile/
alias subcases passed (one unittest, zero skips). No frozen test was edited or
weakened. This provides all 32 frozen-helper cases across the two platforms;
it is not a clean full-suite Windows result. The initial new fixture import
also exposed unavailable WSL NumPy; extracting only the existing harmless
fixture builder removed that unnecessary scientific import. Final required
owner/debug/launcher/collector and template runs have no failures or skips.

Both touched Python files pass py_compile with bytecode directed to temporary
storage. Whitespace checks pass. Memory manifest is regenerated/checked and
ordinary memory health reports only the pre-existing nonblocking CURRENT
soft-size warning, deferred to the next checkpoint/milestone audit. All 985
frozen/unrelated versionable file byte identities and the empty staged contents
match startup. Branch/HEAD and user/local remote refs are preserved; Codex's
automatic internal turn-diff refs changed independently of agent Git actions.
No staging, commit, push, remote update or farm operation was performed.

This establishes local/source behavior only. Installed farm resolver geometry,
quota/capacity, ROOT writes, actual PDFs and runtime preservation remain NOT
VERIFIED for the corrected owner. No science or production acceptance follows.
After independent review and user publication/pushed-state synchronization,
retain the single private A/B farm validation and then seven-strength C NEXT
specified above. No farm command is supplied by this outcome.
