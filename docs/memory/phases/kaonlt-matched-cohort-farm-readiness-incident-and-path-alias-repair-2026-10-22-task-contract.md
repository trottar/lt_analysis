# KaonLT — October 22 farm-readiness incident and bounded path-alias repair

**State:** `BLOCKED`. This is ONE combined repair-and-accountability task, not a new physics phase, a separate memory project, or authorization to rerun the farm. **The desired scientific artifact remains** fresh `Q4p4W2p74` Center/Left low-epsilon canonical t1/t2/t3 same-run A/B diagnostics and signed-event packets, followed by the already-tracked seven-strength C fits/replay and six-parent missing-mass/variance comparison under Lambda. There is NO new matched-cohort science today.

**Audit anchor:** Live public `test` HEAD `f266045e9a535970f54a189a5eecb009401b8adb`, October 22, 2026. Verify local `test`, full HEAD, origin/test, staged/unstaged/untracked tree, pending contracts, review bundles and index before editing. Preserve unrelated state; stop on conflicts. Read repo-root `AGENTS.md` and then `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, `docs/memory/USER.md`, in that exact order. Read task-relevant `MAINTENANCE.md`, `CODEX.md`, `COMMUNICATION.md`, `TOOLS.md`, `LEARNINGS.md`, Oct 15–17 science-displacement record, Oct 21 matched-cohort owner and Oct 22 three earlier repairs. Emit the full tracked health receipt and apply the mandatory user-correction and science-displacement re-anchor before consequential edits.

## 1. Record the actual failure and responsibility accurately

The user explicitly identifies October 22 as an unacceptable, repeated assistant-approved failure of readiness and a day with no new science. **Attribute responsibility to inadequate assistant source-review and readiness judgment, not to user farm operation.** Synthetic tests were repeatedly treated as sufficient despite nonrepresentative filesystem/link and Python dependency behavior in a long-established JLab environment. Existing maintenance rules were not enforced effectively. Avoid minimizing this as a routine Python glitch.

Farm evidence supplied by user/earlier owner excerpts (all **failed and inadmissible for science**):

1. Source `30e029c744ca2638733a1f4c070a71968641f157`, `133634` attempt: owner status `source_preflight`, unresolved irrelevant HeeP symlink target, `analysis_started=false`, `evidence_admissible=false`.
2. Source `842dab7a957a9823a51dcb484a765b4d11c6d9c1`, `144029` attempt: unrelated historical `Q3p0W3p14/Trial_44/Trial_49` symlink loop in broad shared inventory. Actual caller NOT VERIFIED without full status/traceback; no ZIP.
3. Source `a14e740ab68531aca1c60c48d7a73e13e8922aff`, `155418` attempt: owner resolver subprocess failed in `findEffectiveCharge.py` during copied `ltsep -> uproot` import; frozen resolver audit hook rejected import source as `resolver_data_read_forbidden`.
4. Source `f266045e9a535970f54a189a5eecb009401b8adb`, `180522` attempt: new owner-only resolver subprocess failed in `code_path()` with `resolver_code_path_escape` on real Python 3.9 virtualenv module import expressed as `.../lib64/python3.9/site-packages/uproot/__init__.py`. The terminal output DOES NOT prove the canonical target of that path or enumerate the complete installed dependency closure. Owner status mentioned, not yet supplied in full. No child analysis or matched ZIP acceptance is evidenced.

All four failures occurred before requested new scientific evidence; previous October 20 scientific observations do not substitute for fresh matched signed events. Do not assert that a failed receipt or local source tests verified the real farm runtime. In incident memory state explicitly: repeated **SOURCE REVIEWED**/farm-readiness approvals were not supported by representative tests for the actual class of farm environmental dependency; the resulting day of lost scientific throughput is an assistant workflow failure.

## 2. Exact first invariant and one owner-scoped repair

**SOURCE VERIFIED:** In `testing/run_pion_local_t_matched_cohort_gate.py::RESOLVER_PROBE`, `code_path(path)` unconditionally rejects an otherwise valid absolute source path if `p.resolve(strict=True) != p`. This treats *every* legitimate symlinked/aliased dependency location as a path escape. `readonly(open)` contains a related unconditional lexical/canonical comparison; other nearby gates include `Finder`, `Source`, `Extension` loaders, namespace package paths, source identity and config checks. Audit the **entire connected preflight path** so the correction cannot merely move the same failure to the next line.

Change **only** `testing/run_pion_local_t_matched_cohort_gate.py` and `testing/test_run_pion_local_t_matched_cohort_gate.py`:

- A trusted Python code import under a genuine interpreter-controlled `lib64 -> lib` or other *explicitly verified, bounded* dependency-root alias may be admitted when both lexical provenance and canonical containment in the authorized stdlib/site-packages code roots are established. Obtain interpreter roots through active `sysconfig` and validated import search, not hardcoded workstation/farm paths or an unrestricted `sys.path`/`PYTHONPATH` allowlist. Do not treat a symlink from arbitrary external input to site-packages as authorized code merely because its referent is inside.
- Distinguish lexical imported path from resolved source identity; record both where needed for tamper detection. Accept only `SourceFileLoader`-owned source and explicitly checked extension code, never arbitrary package data or `.pyc` authority. Preserve identity/hash and before/after original/copied-package verification.
- Keep actual `ltsep.Root` **default/no-mode and Plot** calls, caller context/cwd, matching config, stable input roles and private `OUTPATH` authority intact. Do not replace resolver calls by static config or fake Root; the separately frozen `testing/attest_pion_chain_ltsep_default_root.py`, `farm_env` and production source must stay unchanged.
- Strictly prohibit scientific `.root`, `.csv`, arbitrary `.json` and other data reads, unauthorized code and package-data reads, ROOT bindings, `ctypes`/subprocess/writes/mutations, path escape, and any changed required input or protection. Don't relax or bypass any real scientific/checkout/installed-package gate. If a safe fix requires an unallowlisted change or depends on genuinely unavailable environment evidence, stop at `BLOCKED` and report the first exact invariant.

**Required local regression coverage before readiness review** (no fabricated farm integration):

1. Reproduce the **prior exact function's** `resolver_code_path_escape` using a real fixture with `lib64` symlink to `lib` and an imported `lib64/python3.9/site-packages/uproot/__init__.py` source. Test `ltsep -> uproot -> numpy` in a subprocess, including baseline and copied-package imports, no-mode and Plot roles, correct caller/cwd and private output checks.
2. Prove the corrected entire owner resolver subprocess and associated `readonly(open)`, `Finder`, `Source`, `Extension` and namespace path checks can pass legitimate code aliases rather than just mocking the failing function. No child or ROOT success claim.
3. Prove rejection of symlinks to outside-code roots, arbitrary external aliases INTO allowed roots, untrusted `PYTHONPATH`, forbidden module/data loads inside site-packages and elsewhere, `.pyc`, ROOT/root_numpy, write/mkdir/link/rename/subprocess operations, wrong/changed package/config identities, modified input roles and nonprivate OUTPATH.
4. Re-run owner and applicable neighbor suites including preservation/collection failure behavior, document skip/fail details honestly, and leave unrelated historical launcher-pinned tests classified separately. Include the earlier missing-target, symlink-cycle and resolver-data-read negative controls so regression coverage explicitly maps to **all four supplied farm failure classes**.
5. `py_compile`, `git diff --check`, review all modified import/provenance guards, repository manifest and memory bootstrap/health. Distinguish covered synthetic invariants from **NOT VERIFIED real installed farm import closure**.

## 3. A durable accountability record AND governing controls (same task)

Create one incident `docs/memory/investigations/kaonlt-matched-cohort-farm-readiness-failure-2026-10-22.md`: sourced chronology of all four failed attempts, source SHA, failure class, evidence limits, missing scientific deliverable, assistant responsibility, test-review shortcomings, preventive conditions and honest unverified details. Do not publish personal `/u/...`, farm worktree, machine or Downloads paths; use relative path roles and `lib64/python3.9/site-packages/...` as an observed path *shape*.

Also update the governing records, not just the new incident:

- `docs/memory/CURRENT.md`: replace superseded preflight status with the `180522` **BLOCKED** real-runtime alias finding and link incident; preserve sole scientific NEXT conditional on actual-diff review, user publication, pushed-state sync and one guarded acquisition; keep **below 16,384 bytes** by replacing/compacting, not appending.
- `docs/memory/MAINTENANCE.md` and `docs/memory/CODEX.md`: an enforceable readiness-review requirement for path-sensitive code: compare the *observed* farm failure class (interpreter/version, lexical/canonical location, symlink topology, installed package import roots and exact caller) against realistic positive and negative tests. Before another farm command after repeated pre-child failures, require regression reproductions of all relevant supplied tracebacks and an audit of adjacent equivalent checks. Synthetic tests never establish **RUNTIME VERIFIED** or unqualified farm compatibility; state exactly what remains unknown and stop if a known failure is not addressed. Don't invent a separate process/phase.
- `docs/memory/USER.md`: make the user's explicitly identified failure priority durable: science throughput, accountability, no repetition of synthetic-only farm-ready claims, no reflexive run or recursive support contracts.
- `docs/memory/LEARNINGS.md` and `docs/memory/MEMORY.md`: concise reusable lesson and reference to incident; avoid duplicating the whole chronology.
- `docs/memory/investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md`: short link/attribution to the new incident and status of absent new event evidence; preserve all prior accepted scientific findings.
- Regenerate `docs/memory/manifest.json`, check memory health/bootstrap (nonstrict), record exact byte counts and nonblocking warnings; no standalone memory-cleanup phase.

Do **not** claim these tracked memories have been updated or the issue has been eliminated until the full actual diff is reviewed and user-pushed. The standalone contract is an instruction, not publication.

## 4. File allowlist, review, and earliest valid science

Allowed files, **and no others**:

- `testing/run_pion_local_t_matched_cohort_gate.py`
- `testing/test_run_pion_local_t_matched_cohort_gate.py`
- `docs/memory/CURRENT.md`
- `docs/memory/MEMORY.md`
- `docs/memory/USER.md`
- `docs/memory/MAINTENANCE.md`
- `docs/memory/CODEX.md`
- `docs/memory/LEARNINGS.md`
- `docs/memory/investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md`
- `docs/memory/investigations/kaonlt-matched-cohort-farm-readiness-failure-2026-10-22.md` (new)
- `docs/memory/manifest.json`
- `docs/memory/phases/kaonlt-matched-cohort-farm-readiness-incident-and-path-alias-repair-2026-10-22-task-contract.md` (this contract)

Codex checks precise worktree and protected paths; preserves unrelated bytes and refs; produces one complete review diff including **full new/untracked files** without staging for review. Codex must NOT commit, push, or touch/run JLab farm. The user alone publishes accepted changes and runs the farm after ChatGPT actual-diff and pushed-state review.

**Single next evidence gate:** one existing owner-operated user farm acquisition of fresh Center/Left lowe t1/t2/t3 A/B and event packets once source review, publication, state sync and readiness have actually passed. On any failure stop, retain failed status/traceback; no automatic retry or scientific interpretation. Only an accepted bundle licenses the later seven-strength C replay and six-parent MM/variance comparison.
