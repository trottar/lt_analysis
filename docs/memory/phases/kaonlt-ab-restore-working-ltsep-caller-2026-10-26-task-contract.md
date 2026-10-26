# KaonLT A/B — restore the proven group-side `ltsep` caller with volatile-only private data

## Authority, entry gate, and exact failure

- Repository: `https://github.com/trottar/lt_analysis/tree/test`.
- Expected published starting `test` HEAD: **`f77b73786960f8a4136267372efa2b607f064706`**. Codex must establish *actual local* `test`, `origin/test`, staged/unstaged/untracked worktree state and stop on an incompatible identity or overlapping edit. Never reset/clean/stash, commit/push, or execute farm work.
- Follow root `AGENTS.md`, then read in this exact order: `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, `docs/memory/USER.md`. Then read `docs/memory/MAINTENANCE.md`, CURRENT's linked active task contracts, the October 19 A/B owner at `6d13abb4681b938c56def2a1f455244b0880096c`, and current source/tests. CURRENT's active objective and scientific NEXT remain authoritative.
- **RUNTIME VERIFIED from supplied `ifarm2402` traceback, October 26 UTC:** the owner at `f77b737...` reaches an isolated copied `ltsep` import but its `isolation.probe_paths(worktree, ...)` fails *before ROOT analysis* with `ltsep_probe_failed` / `ltsep.pathing.InvalidPath` / `PATHFILE not defined`; the rejected installed lookup shows `/lustre24/expphy/volatile/hallc/c-kaonlt/${USER}`. The worktree was `/lustre24/expphy/volatile/hallc/c-kaonlt/trottar/kaonlt-ab-analysis-0o8v_67p/analysis`. The frozen `farm_env/ltsep_paths.py::create_ltsep_root` deliberately passes `os.path.realpath(caller_path)` into `ltsep.Root`. The installed `ltsep` environment-selection contract is **path-dependent**; moving the caller's *analysis checkout* to the volatile filesystem was incompatible with the established farm lookup. Do not edit `ltsep`, the installed `.path` file or the `farm_env` resolver to conceal this error.
- **Previously RUNTIME VERIFIED diagnostic precedent:** the October 19 A/B owner successfully ran a source-controlled detached *analysis code* worktree created using the frozen `testing/run_e8_4_fix5_canonical_five_plot_gate.py::owned_worktree(repo, sha, "analysis")` under the ordinary repository's `/group` sibling. Its established pathing and copied-`ltsep` baseline worked. The ordinary analysis `OUTPUT` is already a symlink to volatile. This historical success proves that *code placement/caller identity*, not private data storage, should follow the group-side convention. It does not prove the new split-private-data owner.
- Earlier gates established `VOLATILEPATH=/volatile/hallc/c-kaonlt/trottar`, canonical `/lustre24/expphy/volatile/hallc/c-kaonlt/trottar`, trusted `/volatile` aliases, the installed Python 3.9 package, 19 source-derived ROOT/hist inputs totaling 8.85 GiB, and adequate filesystem-level capacity (user/project quota limits remain separately qualified). A prior mode `02700` failure was repaired by `f77b737`. These accepted checks and the canonical external symlinks must **not** be reimplemented or discarded.
- The failed `f77b737` run is inadmissible as new A/B science. Its private parent and the preceding failed `e628f727` parent must **remain untouched**; no cleanup instruction is authorized. No ROOT/PDF/preservation claim follows from the traceback.

## Scientific scope and non-negotiable safety boundary

The requested artifact is **one fresh private Q4p4W2p74 Center/Left lowe canonical t1/t2/t3 pion-template A/B comparison** (four JSON/PDF files and existing owner/collection/provenance evidence), then *separately* the unchanged seven-strength C physics comparison. Method A remains detached/non-production; Method B independent, diagnostic and numerically excluded; the production subtraction and `no_empirical_residual` profile remain frozen.

**Restore exactly the previously working division of responsibilities:**

1. **Code/installed-path caller role on `/group`:** create the detached analysis *code* checkout using the **existing frozen** `isolation.owned_worktree(repo, source_commit, "analysis")` or prove an exactly equivalent established group-sibling geometry. The existing `source` collection worktree stays independent and unchanged. Use the proven copied-`ltsep` algorithm under the code worktree's sibling, including its unique installed config matching, byte-for-byte *LTANAPATH-only* substitution, `SetPath`, import identities, three `probe_paths` and actual Root/Plot/cwd caller checks. Prefer calling the unchanged `isolation.prepare_runtime_overlay(repo, worktree)` where its preconditions are now satisfied; do not revise `farm_env/ltsep_paths.py`, the frozen helper, package contents, or PATH_TO_DIR matching rules. The actual `ltsep.Root` caller and patched LTANAPATH must resolve to this `/group`-side code worktree as in the successful October 19 case.
2. **Private data role under existing resolver-provided `/volatile`:** continue using one unique 0700/02700-safe owner-created private **data** parent beneath the already verified **canonical `VOLATILEPATH`**. Keep `ab-output/Analysis/KaonLT`, 19 exact source-derived independently verified ROOT/hist copies, caches, TMP/XDG/MPL/pycache and all sizeable intermediate/scientific outputs **on `/volatile`, never `/group`**. Preserve the already-established original shared `OUTPUT` symlink untouched. Inside the `/group` analysis worktree, create its `OUTPUT` symlink **only to that owner-private volatile `ab-output`**—not to the accepted shared `OUTPUT` directory. The debug child must resolve OUTPATH and findEffectiveCharge's default writer to private volatile storage. Existing literal baseline `VOLATILEPATH`, SIMCPATH/UTILPATH/SKIMPATH and all other stable `ltsep` roles stay unchanged; only private copied LTANAPATH differs.
3. **Separate code from data checks.** Update only the A/B owner-local placement and helper call sites presently assuming `analysis`, `runtime-python`, `ab-output` and caches all share `worktree.parent` on volatile. Now check *code worktree and copied package* against the previously working trusted `/group` sibling, and *private output/copy/cache/temp* against the existing attested volatile parent. A code-only copied installed `ltsep` under the proven group-side worktree parent is permitted, but **no ROOT/hist data or large generated science data on `/group`**. Preserve strict provenance, canonical owner/UID/mode/device/inode and symlink containment appropriate to each role. Never relax to unsafe shared writes or invent a third filesystem root. Snapshot and recheck the two role identities, and check every writable data directory before bulk staging and before the child.
4. **Preserve the owner behavior:** identical source ancestry/commit and profile gate, debug environment and **one** `-d` analysis child, Center/Left lowe required inputs, high-epsilon raw prepass, Diamond filename namespace and tie-breaking, optional SIMC handling, external-link *no-repair* checks, byte-verified independent inputs, initial/final ordinary checkout and accepted-artifact identities, installed/copied-`ltsep` identity checks, skim **existence-only** inventory (no giant hashes), four A/B artifact validators, fresh six-page PDFs, independent clean-source collection, ZIP/provenance/CRC/manifest, companion status/log/summary, and fail-closed cleanup-before-delivery. `run_Prod_Analysis.sh`, science and subtraction are immutable. Do not copy 8.85 GiB to `/group` as an alternative or use the ordinary shared analysis outputs for staging.
5. **Failure cleanup:** the code-only worktree uses only the frozen bounded group-worktree cleanup; the volatile data parent uses the already reviewed strict `0700`/`02700` identity and empty-first-check cleanup, plus safe owned-tree cleanup after successful validation. Track the two lifetimes explicitly so both are checked and removed before successful delivery, but neither cleanup ever follows unrelated aliases, touches accepted outputs, or cleans an earlier failed attempt. On cleanup uncertainty fail closed, preserve failed status, do not delete outside the fresh two owner-created areas.

## Required pre-child review and regression evidence

- **Directly test the regression that was missed:** the actual frozen `farm_env/ltsep_paths.py` uses `os.path.realpath()` for Root callers. Synthetic installed-package fixtures must enforce `ltsep` environment selection using a **group-side path** and must reject the same canonical volatile code-worktree shape supplied in the traceback (not merely accept any absolute directory). Demonstrate old `f77b737` volatile-code routing fails this representative selector, whereas the restored group-side code worktree matches; record that these are synthetic tests and the real installed Python 3.9 `ltsep` remains farm-only.
- Exercise the **actual** A/B owner entry through source preflight -> group analysis worktree -> copied package baseline / caller probes -> volatile-private link -> exact stage copy -> one mocked child -> four artifacts -> independent source collection -> preservation/cleanup/delivery with real temporary symlinks and distinct Linux filesystem devices. Explicitly confirm *no copied ROOT/hist file, ab-output, TMP/XDG/MPL/pycache, analysis result or bulk cache* is placed on the group-side parent, while it holds only the detached code, harmless copied package and small configuration files.
- Negative/failure cases must cover: accidental `/group` copy or cache, private link to ordinary/shared accepted OUTPUT, resolver VOLATILEPATH alias escape, writable/untrusted alias, wrong device, swapped parent/child/output/cache, wrong UID/mode (including `02700` allowed and actual identity preserved), installed package or copied PATH_TO_DIR mismatch, Root/Plot/default caller disagreement, accidental output collision, copied-input mismatch, cleanup failure. No child or accepted ZIP on a pre-child failure. Reuse existing negative science/provenance cases rather than weaken or drop them. Retain an explicit test that both historically failed parent names are *never enumerated/removed* by current cleanup.
- Audit the exact unchanged debug launcher path for any writer **outside the volatile-private output/cache**; code-local lightweight files may follow the historical analysis worktree, but a *large ROOT/hist/simulation artifact to `/group`* is a hard stop. Do not claim global no-group writes merely from symlink layout; document any unverified writer rather than silently weaken this invariant.
- Run the affected Linux A/B owner suite plus the relevant debug-owner, launcher, collector and detached template tests; establish no regressions relative to the recorded 144-test baseline (update tests only for intentionally changed code/data location semantics). If NumPy prevents the frozen helper suite on WSL, use the existing split-platform evidence carefully and report what was actually run; do not invent an installed-package PASS. Python 3.9 syntax, `py_compile`, whitespace, memory-manifest write/check and ordinary health check are mandatory. The long-standing CURRENT soft-size warning is nonblocking unless a new material contradiction is introduced.

## Edit scope and handoff

**Strict five-path allowlist** (one new standalone contract):

1. `testing/run_pion_t_template_comparison_gate.py` — narrow group-code/volatile-private-data split and corresponding existing owner-local callers/guards/cleanup; reuse frozen helper.
2. `testing/test_run_pion_t_template_comparison_gate.py` — representative path-selection and split-placement regressions; retain all prior relevant cases.
3. `docs/memory/phases/kaonlt-ab-restore-working-ltsep-caller-2026-10-26-task-contract.md` — this contract, plus concise implementation results.
4. `docs/memory/CURRENT.md` — minimally record the actual new `f77b737` farm failure and exact next step while preserving the scientific objective and closed statuses.
5. `docs/memory/manifest.json` — regenerate/check.

**Forbidden:** any `src/`, launcher, `farm_env/`, frozen shared `testing/run_e8_4_fix5_canonical_five_plot_gate.py`, profiler JSON, collector, installed `ltsep`/`.path` or external symlink, external farm directory, existing failed parent, accepted artifacts, scientific weights/formulas/binning/fits, Git publication, farm execution. Stop with one precise blocker if the narrow split cannot be proven without editing outside the allowlist.

**Required return:** entire actual diff including the new contract, evidence that all unrelated user changes/index/refs were preserved, exact tests and errors, reproduced old path-selection failure and new code-path success, frozen source assertions, memory health, remaining actual farm limitations, and a short user-controlled Git handoff. Do not stage/commit/push. ChatGPT reviews the *actual* diff before user publication; then independent pushed-state verification. **No automatic farm command or another production/pathing redesign** follows. Farm execution is user-only after a separate evidence-based approval; successful local tests prove no ROOT/PDF or farm preservation behavior.

**NEXT:** One bounded Codex implementation of the proven group-code/volatile-data split, followed by actual-diff review. Return to the fresh A/B scientific artifact, then seven-strength C; no new phase.


## Local outcome - October 26, 2026: BLOCKED

The split implementation is **not delivered**. The mandatory frozen code-worktree
helper cannot satisfy the requested swapped-parent cleanup refusal in the
representative Linux fixture. This is the contract's hard stop, not farm readiness
or a replacement architecture. Independent actual-diff review of this diagnosis
is NEXT; no implementation publication or farm command follows automatically.

Starting branch test, HEAD and local origin/test matched
f77b73786960f8a4136267372efa2b607f064706. The startup core was read in order.
Pre-existing documentation bytes and the empty staged contents were captured.
Startup manifest drift was exactly the newly supplied contract, with no changed
or missing existing entry. CURRENT's prior mode-repair wording was reconciled
minimally to the newer supplied caller failure and this local cleanup blocker.

### First non-bypassable invariant

SOURCE VERIFIED: frozen cleanup_path resolves parent before checking
not parent.is_symlink(). A symlink from a fresh group-side parent to another
canonical sibling can therefore pass its lexical/containment check.
owned_worktree unconditionally invokes this cleanup in finally, then calls
Git worktree remove using the original pathname. An owner-local guard rejecting
the swap does not prevent that inner finally from executing. TemporaryDirectory
later raises Cannot call rmtree on a symbolic link, after the removal call.

The retained Linux regression uses the unchanged frozen helper and the existing
fixture's mocked Git. It creates only disposable ordinary/volatile/unrelated
fixture directories, replaces the fresh group parent with a real symlink to the
unrelated sibling, detects the swap, and requires no removal request and an
unchanged unrelated sentinel. Actual result:

```text
expected: {'removals': 0, 'sentinel_preserved': True}
observed: {'removals': 1, 'sentinel_preserved': False}
OSError during helper cleanup: Cannot call rmtree on a symbolic link
```

This proves the helper issues a removal request through the swapped path; the
synthetic Git callback performs real deletion of the unrelated fixture sentinel.
It does NOT establish how actual farm Git would handle that request. No real Git
worktree add/remove, analysis child, farm path or previous attempt was mutated.
The original fresh code directory remains under its renamed fixture path until
the outer disposable fixture's teardown.

The failed invariant is refusal to clean outside the newly owned areas. Fixing
only the caller layout while keeping this mandatory helper unchanged would not
pass that invariant. The helper is explicitly frozen/outside the allowlist; no
monkeypatch, guard weakening, alternate cleanup mechanism or helper edit was
introduced. The attempted split edits were restored to exact captured startup
bytes. The delivered owner source is unchanged. The two new tests and minimal
records are diagnostic evidence, not a deployable implementation candidate.
Independent review must resolve this exact contract/helper conflict before
implementation resumes; no second support contract was generated here.

### Caller-selection reproduction

The other new regression uses real Linux symlinks/devices, unchanged farm_env
create_ltsep_root and its realpath conversion, real harmless copied-package
subprocess imports/SetPath, and a synthetic installed-package selector requiring
the group-side ANALYSISPATH. It reproduces the canonical volatile code-worktree
shape under expphy/volatile/.../kaonlt-ab-analysis-.../analysis and raises
ltsep_probe_failed / ltsep.pathing.InvalidPath / PATHFILE not defined before
analysis. The unchanged group-side owned_worktree and prepare_runtime_overlay
then pass all three actual path probes in the same fixture. Copied LTANAPATH is
the group-side code worktree; literal VOLATILEPATH and installed package bytes
are unchanged. This validates the representative caller diagnosis only, not a
corrected complete split owner or installed farm Python 3.9.

Disposable stand-ins for both historical failed parent names retain their KEEP
bytes. The new test forbids glob enumeration of the volatile base. These are
synthetic preservation tests; the actual old farm parents were never accessed.

### Deterministic checks and limits

```text
python3 -B -m unittest testing.test_run_pion_t_template_comparison_gate testing.test_run_e8_2_left_lowe_scientific_audit_gate testing.test_run_prod_analysis_debug_left_low testing.test_collect_pion_hgcer_validation_bundle testing.test_pion_t_template_comparison -q
```

WSL Ubuntu / Python 3.12.3: 146 tests in 50.029 s, 145 passed, one safety
regression failed, zero skips. All prior 144 tests passed unchanged: 51 original
owner tests, 27 debug-owner, 15 launcher, 28 collector, 23 detached comparison.
Two new owner tests were added: caller reproduction passes; cleanup refusal
fails with the exact assertion above. All existing test/fixture function source
is unchanged. The earlier exploratory split run failed (51 tests: 4 failures,
21 errors); its incomplete fixture/signature changes and implementation were
removed, not treated as accepted evidence. The initial diagnostic regex missed
newlines in the traceback and was corrected; this changed no guard behavior.

The frozen helper is exercised directly by the two new regressions. Its full
scientific-dependency suite was not rerun at this hard stop; the previous
split-platform helper results are historical, not revalidated here. No WSL
NumPy installation or fabricated installed-package success was attempted.
Both owner/test files pass py_compile with temporary bytecode and Python 3.9
syntax parsing. Actual Python 3.9 execution remains NOT VERIFIED. Elevated
local git diff --check passes; a sandbox-only attempt returned 128 due local
access limitations. No Git configuration was changed.

SOURCE VERIFIED writer audit: debug mode excludes merge and iterative launcher
execution. main.py derives foutroot, generated SIMC ROOT/hist paths and its
cache_transfer archive from OUTPATH. findEffectiveCharge derives its CSV sink
from LTANAPATH/OUTPUT. set_sig_fortran modifies lightweight code-local model
files; binning/parameter/iteration text may follow the historical code worktree.
No large group-side ROOT/hist writer was established by this bounded source
inspection. Full transitive/runtime writer coverage, native dependency caches,
ROOT temporary behavior, corrected split placement and all runtime preservation
remain NOT VERIFIED. The audit is not global no-group-write acceptance.

Scientific source, launcher, farm_env, frozen helper, installed package/config,
external symlinks, accepted artifacts and previous failed farm attempts are
unchanged. No staging, commit, push, remote-ref update or farm execution occurred.
The complete actual diff includes the new contract and failing diagnostic test;
it must be reviewed as BLOCKED/inadmissible for an implementation launch.


## Independent review amendment - October 26, 2026

This is the approved narrow amendment to this existing contract, not a second
contract or a new phase. It supersedes only the shared-cleanup freeze and
five-path limit above. The recorded BLOCKED outcome remains accurate history.
No completed implementation, publication or farm approval is implied.

### Independent findings and decision

Local test, origin/test and GitHub test were independently verified at
f77b73786960f8a4136267372efa2b607f064706. The review reproduced 146 Linux
cases: 145 pass, one cleanup-safety failure, zero skips. The production A/B
owner and shared helper have no diff from HEAD. The blocked patch matches the
worktree. Compilation, Python 3.9 syntax, manifest and ordinary health checks
pass apart from the existing nonblocking CURRENT size warning.

The group-restricted caller fixture supports restoring the established group
code-worktree/copied-ltsep geometry while retaining private volatile data.
The shared cleanup conflict is concrete: cleanup_path resolves away a swapped
parent link before checking it, and owned_worktree requests Git removal through
the original path. Sentinel deletion occurred in a mocked-Git fixture, not on
the farm.

An additional independent disposable-fixture test rejected a replaced parent
before Git removal (zero Git removal calls), yet TemporaryDirectory still
deleted the replacement sentinel. A guard-only patch is therefore insufficient:
both explicit removal and implicit context cleanup need identity protection.

APPROVED FOR IMPLEMENTATION: complete the original group-code/volatile-data
repair using the established helper with only its cleanup safety corrected.
Preserve the resolver, copied-ltsep algorithm, launcher and scientific path.
The current diagnostic package is not a completed repair or a farm candidate.

### Exact scope exception

The seven allowed paths are the existing five plus:

- testing/run_e8_4_fix5_canonical_five_plot_gate.py: only cleanup_path,
  owned_worktree and strictly necessary standard-library imports.
- testing/test_run_e8_4_fix5_canonical_five_plot_gate.py: focused cleanup
  regressions and minimally corrected harmless fixture creation.

All other shared-helper functions, including prepare_runtime_overlay and
scientific/lineage/collector behavior, remain frozen. This operational exception
reopens no closed scientific phase or accepted historical runtime scope.

### Required repair and adjacent checks

1. Preserve the helper API, analysis/source roles, detached source identity and
   group-sibling code geometry. Capture the freshly created parent and worktree
   identities. Reject lexical symlinks before canonicalization; verify canonical
   containment, owner, device/inode and captured identity before any Git or
   recursive removal. Preserve safe inherited 0700/02700 and observed modes.
2. Remove the implicit destructive cleanup route. Use an explicitly managed
   unique parent or an equally proven mechanism so TemporaryDirectory cannot
   recursively remove a replacement after a guard rejects it. Initial validation
   failure permits only identity-verified empty-parent removal. Retain swapped,
   rebound, foreign, uncertain or nonempty unvalidated paths with diagnostics.
   Never restore/delete replacements or scan/clean earlier attempts.
3. Cover both analysis and source roles. Retain the current failing assertion
   and make it pass without skips, expectedFailure, weaker assertions or an
   optimistic Git mock. Add directory replacement as well as symlink replacement,
   guard rejection followed by context exit, first-check/add/body/cleanup
   failures and successful cleanup. Require zero unsafe removal requests,
   unchanged unrelated sentinels and no success delivery after failure.
4. Reconcile the A/B committed-source guard with these exact two newly authorized
   files: its current three-file allowance would otherwise reject publication
   as unexpected_committed_source. Accept only the reviewed successor and pinned
   helper/test bytes; reject later edits to those frozen blobs and unrelated
   paths. Preserve historical scientific pins, exact requested SHA, dirty-source
   rejection, committed producer hashes, detached clean-source collection and
   the unchanged collector checks. Add positive and negative source-gate tests.
5. Complete the original split with the group-restricted selector fixture and
   unchanged copied-ltsep algorithm. Retain canonical volatile, alias, setgid,
   input-copy, Diamond, optional-SIMC, caller and preservation checks. Group-side
   code/package/config is permitted; bulk ROOT/hist, science outputs and writable
   cache/temp remain private on volatile. Trace the reached writers and both
   lifetimes before calling implementation complete.
6. Run the affected owner/debug/launcher/collector/template cases, both cleanup
   roles and relevant shared-helper regressions. Report exact platform/dependency
   limits. Verify Python 3.9 syntax, compilation, whitespace, manifest/health,
   complete seven-path diff and unrelated-state preservation. Do not stage,
   commit, push or run the farm.

NEXT: Codex updates this existing contract's implementation record, implements
only the amended repair and returns the complete actual diff for independent
review. No additional supporting contract or automatic smoke/rerun is authorized.
After source acceptance and user publication, verify the pushed state and reuse
applicable established farm evidence when deciding the single scientific run.
Corrected split runtime, ROOT/PDF behavior and preservation remain NOT VERIFIED.


## Amended implementation outcome - October 26, 2026

`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`. This supersedes the local
BLOCKED implementation status above, preserving its diagnosis and approval
history. Independent actual-diff review is NEXT; no publication, farm readiness,
ROOT/PDF acceptance or scientific promotion follows from local checks.

Starting test HEAD and local origin/test remain
f77b73786960f8a4136267372efa2b607f064706. The amended startup captured the prior
diagnostic test/records, approval amendment and blocked review artifact. No
incompatible source edits existed. Startup manifest drift was exactly the
amended contract bytes. Only the seven authorized paths are implementation
changes; prior unrelated documentation bytes and review artifact are preserved.

SOURCE VERIFIED: shared cleanup_path rejects lexical symlinks before resolution,
checks canonical sibling geometry, actual EUID ownership, private parent
0700/02700 modes and device. owned_worktree captures repository, parent and
worktree identities, including observed modes and device/inode. Explicit mkdtemp
ownership replaces implicit TemporaryDirectory deletion. Parent/worktree swaps,
changed identities, foreign ownership or unsafe cleanup capability stop before
Git removal and recursive deletion. Initial check/add failure may remove only
the verified empty parent; partial-add/nonempty/unvalidated paths are retained.
No earlier-attempt enumeration or cleanup occurs. Both analysis and source APIs,
detached HEAD checks, original sibling geometry and normal Git worktree removal
are retained. All other shared-helper functions, particularly the copied-ltsep
algorithm, are byte-identical to HEAD; only standard-library stat was added.

The original failing symlink sentinel assertions now pass unchanged. The added
directory-replacement regression also retains replacement/unrelated sentinels
with zero removal requests. Five shared test methods exercise both roles:
parent/child symlink and directory replacement, guard rejection/context exit,
first-check/add/partial-add/nonempty/body/removal failures, success, wrong UID,
mode changes and real inherited 02700. Harmless existing fixtures now create a
private parent and model actual Git directory creation/removal; assertions and
scientific fixture semantics were not weakened.

SOURCE VERIFIED: the A/B owner reuses the corrected owned_worktree and unchanged
prepare_runtime_overlay for group-side code/package placement. Its separate
volatile data context retains the already reviewed alias/device/UID/mode/inode
checks and first-check empty cleanup. OUTPUT points only to the fresh volatile
ab-output/Analysis/KaonLT; staged ROOT/hist copies and TMP/XDG/MPL/pycache remain
under the canonical volatile parent. Code/package snapshots and data snapshots
are checked separately before copying/child execution and at cleanup. Both
areas must disappear before delivery. Failed code cleanup retains volatile data
and prevents success. No additional root or installed config change was added.

The source guard retains all historical scientific pins and original three-file
scope, extending it only to the two exact cleanup files. It requires ancestry
from f77b737 and these candidate Git blob identities for independent review:

```text
a05f041ceaf2a1ef5798002272a1e29a22083b83 testing/run_e8_4_fix5_canonical_five_plot_gate.py
48cc342829933da42cf925f453ad19fe352a97c7 testing/test_run_e8_4_fix5_canonical_five_plot_gate.py
```

The test checks these against actual normalized file bytes. Positive successor,
wrong ancestry, dirty helper, unrelated path and changed helper/test blob cases
are covered. Requested SHA, branch/origin observation, committed producer hashes
and independent clean-source collector gates remain unchanged.

### Final deterministic results

WSL Ubuntu / Linux Python 3.12.3:

```text
python3 -B -m unittest testing.test_run_pion_t_template_comparison_gate testing.test_run_e8_2_left_lowe_scientific_audit_gate testing.test_run_prod_analysis_debug_left_low testing.test_collect_pion_hgcer_validation_bundle testing.test_pion_t_template_comparison -q
```

148 passed, zero skips/failures, final run 17.612 s: 55 A/B owner, 27 debug-owner,
15 launcher, 28 collector and 23 detached comparison. All prior cases remain;
location assertions now distinguish disposable group code from volatile data.
The first amended owner run had one diagnostic-regex mismatch after nested
cleanup retained the unsafe code parent; the diagnostic now explicitly reports
code identity change or incomplete cleanup. The final safety assertions pass.

The group-restricted installed-package fixture executes the exact old f77b737
volatile routing and reproduces InvalidPath / PATHFILE not defined. Live code
uses group-side callers and passes actual default/Plot/cwd probes through the
unchanged realpath wrapper. A full owner-entry fixture uses real harmless copied
package imports, SetPath/config substitution and private caller probes, genuine
Linux symlinks and distinct /tmp versus /dev/shm devices. It independently stages
21 synthetic ROOT/hist files (1017 bytes), invokes one mocked child, validates
four fixture artifacts, runs the real collector with mocked source observations,
checks CRC/provenance and cleans both parents before delivery. Its group code
parent contains no staged ROOT/hist, private output or runtime-cache directories.
Both previous-attempt-name stand-ins retain their bytes. Git worktree operations,
source preflight observations, external link setup, ROOT analysis and PDF page
backend are mocked; no farm or actual scientific analysis is executed. The real
19-file / 8.85 GiB inventory remains supplied historical evidence.

Shared-helper Linux coverage: 36 tests passed in 5.600 s, comprising 31 existing
OwnerTests and all five new CleanupSafety tests. The exact test module AST was
executed without its unrelated scientific page-fixture import, using the existing
accepted_pages() fixture output generated by Windows Python. Test methods were
not rewritten or mocked away. One unchanged scientific-lineage reconstruction
test was explicitly not run in this Linux route because NumPy is unavailable.
Direct full-module WSL import independently fails with ModuleNotFoundError:
numpy. A native Windows OwnerTests run (32 tests, 10.695 s) failed with 23
subtest failures and 10 errors: required POSIX symlink-safe cleanup/geteuid are
unavailable, plus the known Windows symlink privilege limitation. Those Windows
failures are not a full-suite PASS, and no guard was weakened for Windows.
The one scientific-lineage test remains NOT VERIFIED in this implementation;
its scientific/lineage source is frozen. This limitation is distinct from the
passing Linux operational and cleanup tests.

Four edited Python files pass py_compile with temporary bytecode and Python 3.9
syntax parsing; execution under Python 3.9 remains NOT VERIFIED. git diff --check
passes. Manifest write/check and ordinary memory health accompany the review.
CURRENT's existing soft-size warning remains nonblocking, for the next normal
checkpoint/milestone audit; no maintenance phase is created.

### Direct writer audit and preservation boundary

Current unchanged source derives main ROOT output from OUTPATH, SIMC archived/
model ROOT and hist paths from OUTPATH, and cache_transfer from OUTPATH.
findEffectiveCharge's default CSV sink is LTANAPATH/OUTPUT/Analysis/KaonLT.
The debug launcher skips git clean, excludes merge/iterative execution, runs
both canonical prepasses and the all-low child, then exits before high analysis.
set_sig_fortran and binning/parameter/iteration text may write lightweight files
inside the disposable code worktree. No reached large group-side ROOT/hist
writer was established. Full transitive native-library initialization/caches,
ROOT temporaries, real package selection and output preservation are NOT VERIFIED
until applicable farm evidence. The source trace is not a global no-group-write
or corrected farm acceptance claim.

All 987 frozen/unrelated captured files, including the earlier blocked review
artifact, retain their byte hashes. HEAD, user branch/remote/tag refs and indexed
modes/blobs match startup; staged contents remain empty. Raw .git/index bytes
changed during local Git inspection; cause is NOT VERIFIED, and no index
restoration or staging was attempted. No commit, push, fetch, farm execution,
external-link mutation, scientific-source change or previous-attempt cleanup
occurred. The temporary complete actual-diff bundle contains all six tracked
changes and the full untracked contract, for independent review before any
user-controlled publication.

Final memory receipt: CURRENT 16244 bytes, MEMORY 29898 bytes, CURRENT_HANDOFF
323 bytes. Manifest PASS; ordinary health WARN only for CURRENT > 8192 bytes
(nonblocking, next normal checkpoint/milestone audit). Final required Linux
148-test and applicable 36-test shared runs both pass with the lexical-first
helper and final pinned blobs. No checks establish actual Python 3.9 execution,
corrected farm runtime, ROOT/PDF production or scientific acceptance.
