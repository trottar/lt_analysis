# KaonLT A/B owner — minimal Lustre setgid-directory compatibility repair

## Authority and failure evidence

**Start from the actual current local `test` and `origin/test`; the expected reviewed pushed base is `e628f727d780a4d87c2ce928fedfcdc62958eb00`.** First follow repository-root `AGENTS.md` and read `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, and `docs/memory/USER.md` in exactly that order. Read `docs/memory/MAINTENANCE.md`, the October 25 A/B volatile placement/alias contracts, and the relevant owner and tests. Inspect local branch, HEAD, index and all tracked/untracked changes before editing. Stop on overlapping edits or unexpected HEAD; do not reset, clean, stash, commit, push, or run the farm.

**RUNTIME VERIFIED on ifarm2402, October 25 local / October 26 UTC:** The single owner invocation at source SHA `e628f727d780a4d87c2ce928fedfcdc62958eb00` failed at `runtime_isolation` with `runtime_directory_owner_or_mode:/lustre24/expphy/volatile/hallc/c-kaonlt/trottar/kaonlt-ab-analysis-c16r9oxs`. A separate read-only `lstat` of that **already-created** directory showed **EUID = 8874, owner UID = 8874, mode = `0o2700`**. The failed receipt reports `status=failed`, `stage=runtime_isolation`, and no preservation-failure reason. The ordinary source preflight, configured resolver, aliases, package, 19 required inputs (8.85 GiB), and filesystem capacity passed the earlier bounded read-only farm gate. No real analysis child, accepted fresh ROOT/PDF, cleanup success, or production preservation is established by this failed attempt. The failed private parent may remain; **do not delete or alter it in this task**, and do not claim cleanup occurred.

**SOURCE VERIFIED:** `testing/run_pion_t_template_comparison_gate.py::directory_identity(path, private=True)` rejects every private directory unless `stat.S_IMODE(path.stat().st_mode) == 0o700`. A `0o2700` Lustre directory is still inaccessible to group/others: `0o2000` is the inherited **setgid directory bit**, not group read/write/execute permission. The owner calls this exact-mode guard immediately after `tempfile.mkdtemp()` and *before* entering the current `analysis_worktree()` cleanup-protected `try`. This is the first failed invariant and an early-cleanup gap, not a reason to move data, change aliases, chmod existing directories, or redesign the farm runtime.

## Exact smallest correction

1. **Fix only the strict private-directory mode check** in `testing/run_pion_t_template_comparison_gate.py::directory_identity`. Continue to require canonical absolute non-symlink paths, actual current EUID ownership for owner-managed paths, and the complete existing identity / device / inode / canonical containment checks. For `private=True`, accept **only `0o700` and `0o2700`** from `stat.S_IMODE`, equivalently owner `rwx`, group and other **no permissions**, with no special bits other than optional `S_ISGID`. Reject `0o770`, `0o750`, `0o755`, `0o777`, `0o1700`, `0o4700`, `0o3770`, wrong UID, replaced directories, and links. Preserve the actual observed mode in the returned identity record. **Do not call chmod to remove the inherited setgid bit**, change the group, loosen access, or substitute a new directory topology. Existing owner-created `ab-output` and temporary/cache paths must use this same safe predicate where `private=True` is already passed.

2. **Close the pre-`try` cleanup gap narrowly.** In `analysis_worktree`, establish the cleanup-protected scope immediately after creating the unique private parent, so failures during its first identity/permission check do not silently leave an otherwise safely identifiable empty temporary directory. Preserve strict fail-closed semantics: a cleanup attempt may affect **only the unique parent just created by this invocation**, and only after checking its expected canonical location, ownership, device/inode and non-symlink identity. If any check is unsafe or ambiguous, retain the path and fail with an explicit diagnostic; never traverse a swapped link, clean an existing directory, or recursively delete a different path. Use the existing bounded cleanup for validated normal runs. This is not authorization to clean up the October 26 failed parent, or to change source-worktree handling.

3. Preserve **all** established path and analysis behavior: raw resolver `VOLATILEPATH`, trusted infrastructure aliases and canonical volatile base; the already implemented private analysis worktree, `analysis/OUTPUT` symlink, copied `ltsep` algorithm/SetPath selection and `LTANAPATH`-only modification; imported-package and caller probes, external link no-repair checks, root/hist byte-verified staging, Diamond selection, optional SIMC handling, debug child, four A/B files, independent source collector, ZIP/CRC/manifest, receipts and every preservation gate. **No alternative storage path, directory framework, source code physics change, production output modification, or new farm readiness stage.**

4. Tests in `testing/test_run_pion_t_template_comparison_gate.py` must **reproduce the actual rejected mode before the correction**: a real Linux owner-created directory with `chmod(0o2700)` / inherited `S_ISGID` and UID equal to EUID must fail the old exact-mode predicate, then pass the corrected predicate and the actual `analysis_worktree` / `verify_runtime_parent` / `verify_runtime_layout` paths. Use an appropriate disposable Linux synthetic filesystem where setgid inheritance can be exercised; if the fixture cannot reproduce it, report that limitation rather than silently passing. Check both `0o700` and `0o2700`, all denied group/other/special combinations, other-owner mock, symlink swaps, and changed-mode identity after snapshot. Test that all owner-private TMP/XDG/MPL/pycache and `ab-output` directories still reject group/other access. Inject a **first identity-check failure after successful mkdtemp**, assert that only the newly owned empty parent is removed when provably safe, with no unrelated write or analysis child. Retain tests of unsafe-cleanup refusal and prior alias/group-device behavior; do not weaken older tests.

5. Run the same applicable deterministic 139-test Linux group (46 owner, 27 debug-owner, 15 launcher, 28 collector, 23 detached comparison), the relevant frozen helper regressions, `py_compile` with Python 3.9 syntax check, `git diff --check`, generated memory manifest write/check and ordinary health check. Existing nonblocking CURRENT length warning is not a new task. Record direct source expectations separately from tests. Do **not** claim real PyROOT, PDFs, environment, or farm preservation verified by local tests.

## Edit scope and work-state continuity

Strict five-path allowlist:

- `testing/run_pion_t_template_comparison_gate.py` — only the `private` mode predicate and **first-created-parent cleanup lifetime**.
- `testing/test_run_pion_t_template_comparison_gate.py` — focused regressions above.
- `docs/memory/phases/kaonlt-ab-owner-lustre-setgid-mode-fix-2026-10-25-task-contract.md` — this contract plus implementation/test outcome.
- `docs/memory/CURRENT.md` — minimally record this concrete failed farm invariant and preserve substantive scientific NEXT.
- `docs/memory/manifest.json` — regenerate/check.

Everything else is frozen, especially `src/`, `run_Prod_Analysis.sh`, `farm_env/`, profile, frozen shared helper, installed `ltsep`, and already accepted science. Scientific work remains the *fresh Q4p4W2p74 Center/Left lowe canonical t1/t2/t3 A/B diagnostic*, followed by the unchanged seven-strength C comparison. Method A stays detached/non-production; Method B stays independent/diagnostic and numerically excluded. The current gate is **BLOCKED**, not farm-ready.

Provide a full actual diff (including the new contract), unchanged paths/refs/index evidence, focused test results, an exact explanatory comparison to mode `0o2700`, and the narrow remaining runtime evidence. **Do not commit or push.** After independent ChatGPT actual-diff review, only the user publishes; ChatGPT verifies the published state; the farm is user-only. Never automatically chain a source review or readiness observation into another A/B run. This repair itself does **not** authorize a farm run, cleanup of the old failed directory, or any scientific interpretation.


## Local implementation outcome - October 25, 2026

`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING` for the local repair; the farm
execution gate remains `BLOCKED`. Independent actual-diff review is NEXT.
Starting branch test, HEAD and local origin/test all matched
`e628f727d780a4d87c2ce928fedfcdc62958eb00`; the index was empty. Pre-existing
tracked documentation edits and this supplied untracked contract were captured
before editing. Startup manifest failure was exactly the absent new-contract
entry: no existing entry changed or disappeared. No farm operation occurred.

SOURCE VERIFIED: the private predicate accepts exactly 0700 or 02700 and
records the actual mode. 02700 is 0700 plus S_ISGID, with zero group/other
access; it does not require chmod, group changes or a different topology.
Cleanup protection now begins after mkdtemp, before the first lstat/identity
check. If initial validation fails, only an empty newly created parent may be
removed with rmdir after canonical location, unchanged base/aliases, current
EUID ownership, directory type, device/inode and captured mode checks. A missing
snapshot, changed identity/mode, link, wrong owner/device or nonempty directory
fails closed with a diagnostic. Validated-run bounded cleanup is unchanged.
The old farm parent was neither accessed nor cleaned by this task.

Direct source comparison against the starting HEAD finds changes in only
`directory_identity` and `analysis_worktree`; no owner definitions were added.
All other owner definitions, including the copied-ltsep algorithm, execute_gate,
resolver aliases, staging, child, collector, receipts and preservation gates
are byte-identical. Scientific source, scripts, symlinks, frozen shared helper,
profile, installed packages and farm files are outside this repair.

### Deterministic local evidence

WSL Ubuntu Linux / Python 3.12.3. Real /dev/shm and /tmp fixtures retain distinct
devices, actual symlinks and harmless copied-package subprocess probes; Git
worktree operations and the analysis child are mocked. No real analysis runs.
The exact starting-HEAD functions reject a genuinely inherited 02700 directory
owned by EUID and leave the empty parent, reproducing both observed defects.
The corrected predicate accepts the same directory. Both 0700 and inherited
02700 pass actual analysis_worktree, verify_runtime_parent, overlay and layout
checks and bounded cleanup. Failure to reproduce the exact bits would fail
these tests, rather than skip or silently substitute another mode.

Five new test methods retain all prior 46 owner tests. The mode matrix rejects
all 510 disallowed combinations of owner rwx, group/other bits and special bits,
plus missing owner execute, with explicit observed-mode assertions. Wrong UID,
changed accepted-mode identity, all private TMP/XDG/MPL/pycache and ab-output
access checks, and runtime-python's existing private check are covered. Injected
first-check failures on both 0700/02700 remove only the new empty parent, without
recursive removal or a child; unrelated and accepted bytes remain unchanged.
Initial symlink/directory replacement, mode change, wrong UID and nonempty
cleanup refuse removal. Prior alias, group/device and unsafe-cleanup tests remain.

```text
python3 -B -m unittest testing.test_run_pion_t_template_comparison_gate testing.test_run_e8_2_left_lowe_scientific_audit_gate testing.test_run_prod_analysis_debug_left_low testing.test_collect_pion_hgcer_validation_bundle testing.test_pion_t_template_comparison -q
```

144 passed, zero skips/failures (final run 13.124 s): 51 owner, 27 debug-owner,
15 launcher, 28 collector, 23 detached comparison. This retains the specified
139-test group and adds five focused owner regressions.

Frozen helper regressions: WSL full-module import failed because NumPy is not
installed. Existing Windows Python 3.12.10 ran 32 tests in 109.453 s: 31 passed,
one errored at Windows symlink creation (WinError 1314), including elevated
execution. The exact unchanged failing test method and its two fixture methods
were extracted by AST and executed on Linux: one test, all three symlink/nonfile/
alias subcases passed, zero skips. Thus all 32 cases have split-platform coverage;
this is not a full-suite single-platform PASS. No frozen test was weakened.

Both edited Python files pass py_compile with temporary bytecode and Python 3.9
syntax parsing; actual Python 3.9 execution is NOT VERIFIED. git diff --check,
generated manifest write/check and ordinary memory health accompany the review.
The known CURRENT soft-size warning is nonblocking and deferred to the next
checkpoint/milestone audit. CURRENT retains the fresh Center/Left lowe t1/t2/t3
A/B target and unchanged seven-strength C follow-up, without authorizing a rerun.

This is local/source evidence only. Corrected farm behavior, real PyROOT, PDF
rendering, actual runtime preservation and old-parent cleanup remain NOT VERIFIED.
Independent actual-diff review, user publication and pushed-state review remain
required; no automatic farm run follows. No commit, push, staging or user-ref
update was performed. The temporary review bundle contains the complete tracked
diff (including separately identified pre-existing edits) and full new contract.

Final preservation receipt: only the five allowlisted file hashes changed from
the captured baseline; all 987 frozen/unrelated versionable files are byte-identical.
HEAD and user branch/remote/tag refs match startup. Raw .git/index bytes changed
during the session; the cause was not independently established. Every indexed
mode/blob still equals HEAD and the staged diff remains empty, as at startup.
No index restoration or staging was attempted. Existing unrelated worktree edits
are preserved byte-for-byte and separately identified in the review bundle.
Final health: CURRENT 16242 bytes, MEMORY 29898 bytes, CURRENT_HANDOFF 323 bytes;
manifest PASS, ordinary health WARN only for CURRENT > 8192 bytes (nonblocking,
next checkpoint/milestone audit). Whitespace and syntax checks pass.
