# E.8.4 Fix.5.6 — owner farm readiness and failure provenance

**Status:** `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`.

## Source and scope

Started on exact `test` HEAD `df957a6414fc9c515d1f82228517cb801dc90350`,
parent `761fbb6c03d2d7a10bb911cf84e9ba898496fab6`; committed HEAD remains
unchanged. The [contract](e8-4-fix5-6-owner-farm-readiness-and-failure-provenance-task-contract.md)
supplies independent ChatGPT actual-diff and pushed-state review PASS for the
12-path Fix.5.5 source. Fix.5.4/Fix.5.5 remain `DEVELOPMENT COMPLETE, FARM
VALIDATION PENDING`; source acceptance is not farm numerical or PDF acceptance.
Fix.5.6 was `ACTIVE` during implementation and awaits actual-diff review.

Only the tracked Left/lowe owner, its tests and warranted memory change. Frozen
scientific/presentation source, launcher, collector/tests, profile scope,
candidate hashes, artifact definitions and eight-artifact inventory are preserved.
Method A is detached/non-production; Method B numerically absent; absolute-SIMC
interpretation remains blocked with the existing page-level provenance semantics.

## Confirmed gap and repair

Source trace: preflight -> effective in-memory profile binding -> candidates ->
child debug analysis/log -> completion markers -> fresh artifact/page checks ->
summary -> exact detached collector/source checks -> ZIP verification. The old
log records only the child launcher stdout/stderr. Later owner failures appear
only on owner stderr; absence of an owner-failure string from that log cannot
exclude a post-analysis failure. The previous exact failed stage is unknown.

Before analysis, the owner now reuses `clean_collection_worktree` and
`collection_module`, reads the profile with the exact detached collector, binds
only the effective required commit to the supplied SHA, checks profile equality,
then invokes existing `collect_source_checks` and `_committed_identity`. Nonzero
checks include their literal stderr/stdout and return code in the failure reason.
Unexpected committed paths, wrong detached identity, import/profile failures or
source-check failures stop before analysis and create no final ZIP. Bounded
worktree cleanup and ordinary-checkout dirt handling are unchanged. Final
collection still runs independently from an exact clean detached tree with the
ordinary canonical artifacts, and repeats source checks before packaging.

Per-attempt `e8_4_fix5_owner_gate_status/v1` uses the intended ZIP's unique stem
plus `-gate-status.json` in the canonical artifact directory. It records exact
source/setting, ZIP/log paths, UTC timestamps, status/stage, literal reason,
analysis started/completed, source preflight, artifact verification, collection
and ZIP-verification completion flags. Initial JSON publication uses an exclusive
atomic hard link from a flushed/fsynced same-directory temporary file; later
updates use atomic replacement. Existing status attempts are refused untouched.
Filesystem publication failure prevents analysis; no status can be guaranteed
when the artifact directory itself is unwritable.

Stages: preflight, profile, candidate_identity, collector_source_preflight,
analysis, completion_markers, verify_artifacts, write_run_summary, collection,
verify_zip, collection_cleanup, complete. Caught failures persist failed state,
current stage and literal reason; success is complete only after verification
and owned temporary cleanup. Status failure diagnostics go to stderr. Success
stdout remains exactly one POSIX ZIP path; failure stdout remains empty. The
status is excluded from the unchanged declarative profile and ZIP. The child log
is never appended after its run-summary hash is recorded.

## Deterministic verification

Python 3.12.10: **66 tests passed, no skips** across:

- `testing.test_run_e8_4_fix5_left_lowe_plot_gate` (14 tests);
- `testing.test_collect_pion_hgcer_validation_bundle`;
- `testing.test_e8_4_fix5_5_visualization`;
- `testing.test_e8_4_fix5_4_identity_audit`;
- `testing.test_pion_hgcer_validation_bundle_profile_e8_4_fix5`.

Coverage includes pre-run source failure (no analysis/log/ZIP), analysis failure,
completion-marker failure, stale artifact, collector failure, ZIP verification
failure, atomic/exclusive status publication and complete success flags. Real
unchanged generic collection on synthetic artifacts retains inventory/hashes/
source/profile/Left-lowe provenance. Tests check analysis starts only with a
persisted preflight-complete status, final source rechecks still occur, detached
profile/commit/unexpected-path rejection, import failure persistence, frozen log
hash, stderr status path, one-path success/empty-failure stdout, wrong ordinary
branch/head/origin/dirt, preserved known model dirt, existing ZIP refusal,
path-bounded cleanup and exact `-d 4p4 2p74` launcher argv.

In-memory syntax checks for both changed Python files and `git diff --check`
pass. AST audit confirms all existing owner functions/validators remain exact
except the execution wrapper and main failure reporting. Removing only added
status/preflight statements leaves the existing execution steps exactly intact,
including post-analysis checks/summary/collection/ZIP verification. SHA-256
checks pass for all 164 frozen files (src, launcher, collector/tests and profile). Local tests do not establish farm permissions/filesystem
atomicity, ROOT/PyROOT, analysis, delivery, actual PDF visual quality, observed
cancellation or new numerical closure.

## Memory and review scope

Regenerate the intended versioned manifest with unchanged tools in a temporary
candidate view containing tracked memory plus this phase/contract, then run
write -> check -> ordinary health serially and verify every resulting hash/byte
count against the worktree. Preserve/exclude unrelated untracked
`workflow-continuity-hardening-task-contract.md`; its preservation SHA-256 is
`0ee396c059be03e0a528bef18a1fba1128be12e75b63ef725e9cf9a02dffe421`.
Scoped manifest: 199 files; write/check and ordinary health PASS with zero
warnings/hard failures. CURRENT/MEMORY/CURRENT_HANDOFF bytes: 6992/12489/323.
No index, ignore rule or memory tool change. Complete raw `kaonlt_review.diff`
is relative to starting HEAD and includes both intended new memory additions;
it is temporary, excluded from the intended commit.

## NEXT and hard stop

NEXT — Fix.5.6 actual-diff review -> user commit/push -> pushed-state review ->
one narrow Q4p4W2p74 / Left / lowe tracked-owner farm gate -> fresh
ZIP/status/log/PDF/manifest evidence review. CURRENT owns the exact next action.
Hard stop after local checks, memory/manifest and complete diff. No commit,
push, farm command or farm execution occurs in this task.
