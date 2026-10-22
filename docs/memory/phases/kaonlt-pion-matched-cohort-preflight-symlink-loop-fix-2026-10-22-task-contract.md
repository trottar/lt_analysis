# KaonLT — Finish the matched-cohort symlink preflight repair and acquire physics evidence

**Task class:** One minimal fix to the existing matched-cohort acquisition owner after a fresh, demonstrated farm `source_preflight` failure. Not a new scientific phase or permission to investigate/repair historical caches.

**Audit anchor:** Live public `test` at `842dab7a957a9823a51dcb484a765b4d11c6d9c1`, October 22, 2026. Verify actual local `test`, `HEAD`, `origin/test`, staged/unstaged/untracked status and existing review files before edits. Preserve unrelated work, bytes, index and refs; stop at an overlapping dirty source or unexplained identity.

**Required scientific result:** One owned, single-run `Q4p4W2p74` acquisition for **Center/lowe and Left/lowe**, canonical `t1/t2/t3`: fresh matched A/B diagnostics plus signed-event packets, then already-tracked offline seven-strength C fits and six-parent missing-mass/variance comparisons under Lambda. These results are currently **unavailable**. No fit-selection claim, Lambda yield acceptance, production promotion or alternate manual farm procedure.

## 1. Read authority; identify the precise failed invariant

Read repository-root `AGENTS.md`, then `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, and `docs/memory/USER.md`, in that exact order. Then read `docs/memory/MAINTENANCE.md`, `docs/memory/CODEX.md`, `docs/memory/TOOLS.md`, `docs/memory/COMMUNICATION.md`, and the **existing** October 21 matched-cohort attempt-isolation and October 22 dangling-link repair contracts. Emit the repository full health receipt. Follow only their relevant source references.

**RUNTIME VERIFIED from the user's terminal:** On `ifarm2402`, at pushed `test` HEAD `842dab7a957a9823a51dcb484a765b4d11c6d9c1`, the owner failed with `Symlink loop from '/lustre24/expphy/volatile/hallc/c-kaonlt/trottar/OUTPUT/Analysis/KaonLT/cache_transfer/kaon/Q3p0W3p14/Trial_44/Trial_49'`. The prior failure was a missing HeeP referent. Both returned no accepted ZIP or fresh A/B/C science. The **exact current farm caller/stage is NOT VERIFIED** because the October 22 `144029` owner status JSON/traceback was not provided; do not claim otherwise.

**SOURCE VERIFIED:** `testing/run_pion_local_t_matched_cohort_gate.py::link_inventory()` calls `Path.resolve(strict=True)` on every symlink in broad metadata snapshots but catches only `FileNotFoundError`. Python releases differ: loops may raise `RuntimeError("Symlink loop from ...")` or `OSError` with `errno.ELOOP`. `execute_gate` invokes `inventory(shared)` via `shared_state()` before spawning the child. `external_input_inventory()` invokes `inventory()` on broad `CUTPATH` and `PARAMPATH`, and separately records only `missing_target` entries for later comparison. Neither source branch currently handles cycle metadata consistently. **Source-level cause is strong, but the precise farm caller remains unverified.**

## 2. One robust and conservative owner repair

Change the **existing** `testing/run_pion_local_t_matched_cohort_gate.py` only as necessary:

- On *broad metadata-only* symlink inventory, preserve each link's lexical path/key and exact `os.readlink()` literal target. Record a deterministic distinct `cyclic_target` state for both Python exception representations (`RuntimeError` specifically indicating symlink loop, and `OSError` with `errno.ELOOP`). Keep existing `missing_target` representation for genuine `FileNotFoundError`. Do not read/traverse a cyclic referent; do not attempt to resolve it again while recording that state.
- Preserve the existing resolution and byte-identity hashes of all normally resolvable regular-file referents, and all ordinary regular-file hashes and directory metadata. Continue `os.walk(..., followlinks=False)`: no traversal into symlinked directories. A cyclic link is not an input and must not be treated as a valid real directory or regular file. A symlink's *literal target* is metadata, not proof of an actual physical nested trial.
- Update `external_input_inventory()` so **both** missing-target and cyclic metadata records under `CUTPATH`/`PARAMPATH` are stored as unresolved links and rechecked by `verify_external_inputs()` after child and on failure. Shared-output inventory must also retain these records in the existing before/after equality checks.
- No generic catch-all `except Exception`, no swallowing permission/read/stat failures or unrelated `OSError` / `RuntimeError`. `os.readlink` errors, unrelated IO/permission errors, and every required-input verification failure still stop. Do not modify strict `identity()`, `stage_input_files()`, protected path/link checks, required background neutron/delta/sidis inputs, configured Sigma0, external SIMC links, installed/copied ltsep checks, checkout/output preservation, owner failure-to-collect semantics, or collection/provenance logic.
- **Prevent another opaque failure**: if compatible with the existing owner status schema, record a bounded failure traceback identifying the exact source caller in the owner status JSON for preflight exceptions. Keep `failure_reason`, `evidence_admissible=false`, and preservation exceptions authoritative; never mark failure as success. Do not create a new diagnostic runner or extra farm command.
- Do not repair/relink/delete any `Trial_*`, cache, HeeP or shared path. Do not alter installed `ltsep`, launchers, source-controlled physics, observer dispatch, cuts, fits, SIMC, event weights, ROOT containers, normalizations, uncertainties, yields or cross sections.

**Hard stop:** If an actual required KaonLT input or frozen external authority link contains a cycle, leave it a hard failure and identify the exact role/path; do not convert required scientific data to metadata, invent a replacement path or relax file identities. If a safe fix requires modifying an unallowlisted frozen file, return `BLOCKED` with the first exact invariant.

## 3. Thorough, finite deterministic validation (not farm validation)

Expand existing `testing/test_run_pion_local_t_matched_cohort_gate.py` tests only. Required cases:

1. Reproduce the current failure on the **old implementation** with a synthetic `cache_transfer/kaon/Q3p0W3p14/Trial_44/Trial_49` cyclic link. Include self-cycle, two-link cycle and relative ancestor cycle. Verify `os.walk` does not descend via the links or recursively invent nested trials.
2. Verify metadata state, equality and change detection for **all three call sites**: shared `KaonLT` output, `CUTPATH`, `PARAMPATH`. Unchanged cycles permit simulated owner preflight to reach **exactly one** child; they do not imply successful real analysis.
3. Simulate both `RuntimeError("Symlink loop from ...")` (farm Python) and `OSError(errno.ELOOP, ...)` (newer Python); only true loop conditions are admitted as metadata. Other `OSError`, `PermissionError`, or unrelated `RuntimeError` must propagate.
4. Verify failure when a cycle is removed, retargeted, replaced by a file, repaired/becomes resolvable, or appears where an ordinary intact link existed during the child. Verify ordinary regular-file and live regular-link referent byte mutation detection remains unchanged. Before/after preservation must run after a simulated child failure, without allowing collection.
5. Required configured Sigma0 and a required neutron/delta/sidis ROOT link that point through cycles must still fail **before child launch**. All existing source, path, write-boundary, package, observer and hash protections stay strict.
6. Test that failure status still records `failure_reason` and (if implemented) the first caller traceback, never claims admissibility, and does not create/deliver ZIPs. Check that multiple independent preservation checks still run.
7. Re-run the existing targeted owner suites (including event/cohort tests as supported), deterministic synthetic end-to-end preflight, `py_compile`, `git diff --check`, manifest/check, and standard nonstrict memory-health/bootstrap checks. No artificial ROOT/farm success claim. If a suite cannot import NumPy, mark that suite unexecuted, not passed, and do not initiate unrelated dependency repair.

Before handing off, inspect the **complete actual diff and test output** and explicitly list each remaining symlink dereference along source preflight, classifying it as required-authority, required-input, or broad-metadata. No other *unhandled broad-metadata* symlink-cycle route may remain in this owner. Do not change other owners to address speculative issues.

## 4. Source/change allowlist and memory

**Source:** Only
- `testing/run_pion_local_t_matched_cohort_gate.py`
- `testing/test_run_pion_local_t_matched_cohort_gate.py`

**Records:** Only
- `docs/memory/CURRENT.md`
- `docs/memory/investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md`
- `docs/memory/manifest.json` (generated)
- `docs/memory/phases/kaonlt-pion-matched-cohort-preflight-symlink-loop-fix-2026-10-22-task-contract.md` (this one contract)

Keep CURRENT well below the **16,384-byte hard limit** (it is already approximately 16,189 bytes at the audit anchor): replace/summarize superseded preflight wording rather than appending. Record the October 22 `144029` failure as inadmissible, followed by this one pending repair. Do not copy personal farm paths into tracked memory, create another phase, alter a closed result, or start broad memory consolidation. Preserve the exact conditional NEXT: once reviewed/pushed/synchronized, acquire **fresh** Center/Left A/B/C science.

Regenerate/check manifest using the repository's documented interpreter selection, run standard `tools/check_memory_health.py --root .` and bootstrap, report actual CURRENT/MEMORY/CURRENT_HANDOFF bytes, warnings and classification. Nonblocking CURRENT soft-size warnings do not displace science.

## 5. Acceptance / earliest next operation

Provide one complete review diff in the repository root that includes tracked changes and the full new contract addition, without staging solely for review. Report actual branch/HEAD/worktree, exact root cause and Python-version coverage, frozen paths, tests, health, and first hard stop if encountered. **Codex must not commit/push, modify farm filesystem, or run the farm.**

Required chain: Codex localized implementation + tests → independent ChatGPT **actual diff** review → user commit/push → ChatGPT pushed-state verification → **one user-controlled guarded farm acquisition** → successful owner/provenance/ZIP review → fresh offline seven-strength C and event-applied six-parent MM/variance science. No farm command from this task until readiness is independently passed. Any new failure stops at its first invariant; do not automatically retry or package partial evidence.

**Fix success criterion:** Under synthetic preflight, irrelevant dangling/cyclic historical links no longer block the science owner and their presence/change is preserved, while actual KaonLT inputs and all preservation gates remain strict. This establishes source readiness only; there is still no new scientific acceptance.
