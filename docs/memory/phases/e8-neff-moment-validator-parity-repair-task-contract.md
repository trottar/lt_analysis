# KaonLT E.8 — Detached contribution N_eff producer/checker arithmetic parity repair

## 1. Scope and authority

**Task class:** narrowly scoped source-changing repair following the first failed detached current-F.1 materialization farm gate. Authoritative repository: `https://github.com/trottar/lt_analysis/tree/test`.

**Live remote `test` HEAD at contract preparation:** `9145cb287c2ad7eccdd2e5b2f38840abea920f15`. Codex must independently verify live `test`, local `HEAD`, `origin/test`, the index, and complete `git status --short --branch --untracked-files=all` before changing anything. Do not reset, stash, clean, rebase, stage unrelated paths, or replace dirty files to force that gate. A different branch/HEAD, conflicting target-file edits, or unpreservable unrelated work is `BLOCKED` pending re-anchor.

Read the repository-root `AGENTS.md`, then **in exactly this order and in full**: `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, `docs/memory/USER.md`. Next read `docs/memory/MAINTENANCE.md`, `docs/memory/CODEX.md`, `docs/memory/TOOLS.md`, `docs/memory/COMMUNICATION.md`, this task contract and the task-relevant CURRENT references. Emit the full repository-defined KaonLT health check before implementation.

## 2. Directly supplied failed farm gate

The user ran the **pushed and pushed-state-reviewed** detached materializer at `test` HEAD `9145cb287c2ad7eccdd2e5b2f38840abea920f15`:

`testing/materialize_method_a_current_f1_contribution_snapshot.py`

The farm process returned:

```text
materialization_failed:snapshot_failed:neff_moment_mismatch
```

**RUNTIME VERIFIED:** that failure string was supplied from `ifarm2401`. **SOURCE VERIFIED:** in the current owner, after building candidate F.2/F.3/F.4, `manifest_for()` calls `verify_stages()` -> `src/cuts/pion_hgcer_method_a_current_snapshot_contribution_debug.py:reproduce_snapshot()` -> its `validate_parent()` -> `src/cuts/pion_hgcer_method_a_reweighted_contribution_debug.py:validate_summary()`; the final validator enforces exact equality of a reconstructed `neff` metric. **NOT VERIFIED:** the precise failing setting/t/source/weight tuple and real-farm successful F.2/F.3/F.4 acceptance; the error occurred during detached measurement/validation, and no successful materialization ZIP or downstream preflight was accepted. Do not interpret any failed-attempt partial files as accepted evidence, automatically rerun the farm, or reopen the canonical-five owner.

The frozen descriptive `effective(values)` producer computes:

```python
total = math.fsum(values)
numerator = total * total
```

But the `validate_summary(s)` persisted-moment checker recomputes:

```python
expected = {"count": s["count"], **metric(s[total]**2, s[square])}
require(s["neff"][name] == expected, "neff_moment_mismatch")
```

Floating-point multiplication and exponentiation can differ by 1 ULP even though both represent the same mathematical square. A local CPython example, `x = 55.83414573673653`, yields `x*x == 3117.4518301511343`, `x**2 == 3117.451830151134`, so the computed `N_eff` values differ (`1.0` vs `0.9999999999999999`) and exact dict equality rejects a valid result. This is a directly reproducible *source-level arithmetic inconsistency*, not a verified identification of the farm's actual offending record. Scientific definition remains:

\[
N_{\mathrm{eff}} = (\sum_i w_i)^2 / \sum_i w_i^2,
\]

for each explicitly nonnegative measure (`w0`, `w_rw`, `abs_b`, `abs_b_rw`).

## 3. Sole executable repair

Change only the checker arithmetic expression in `src/cuts/pion_hgcer_method_a_reweighted_contribution_debug.py::validate_summary()`:

```python
# BEFORE
expected = {"count": s["count"], **metric(s[total]**2, s[square])}

# AFTER
expected = {"count": s["count"], **metric(s[total] * s[total], s[square])}
```

Preserve the **exact equality** check of the metric dictionary, `neff_exceeds_count`, all signed/absolute moment fields, `math.fsum` semantics, `metric` availability/zero-support policy, producer `effective()` implementation, historical pins, all serializations, and all scientific and geometry logic. This is **operation-order parity only**, not a tolerance relaxation or a new scientific formula. Do not modify `neff` persisted values, alter manifests/ZIP schemas, change weight/correction factors, or backfill/accept previously failed output.

If the exact source differs, or this one-expression parity change is insufficient to establish the described failure's local regression, **stop** and request a fresh narrowly scoped review rather than broadening this contract.

## 4. Exact path allowlist

The only paths eligible to be added or edited for this task are:

1. `src/cuts/pion_hgcer_method_a_reweighted_contribution_debug.py` — the single checker expression above only.
2. `testing/test_pion_hgcer_method_a_reweighted_contribution_debug.py` — add regression tests for the corrected numerical parity and fail-closed tamper behavior.
3. `docs/memory/phases/e8-neff-moment-validator-parity-repair-task-contract.md` — **this single standalone contract** (new file).
4. `docs/memory/CURRENT.md` — minimum active-state/failed-gate wording and **one** push-stable NEXT conditional on source review/publication; preserve all accepted statuses/other content. A soft-size warning alone is nonblocking and may be consolidated at a meaningful later checkpoint.
5. `docs/memory/manifest.json` — regenerate only using the tracked manifest owner.

**All other paths are frozen**, specifically the F.1/F.2/F.3/F.4 scientific builders and calculators; `src/cuts/pion_hgcer_method_a_current_snapshot_contribution_debug.py`; `testing/materialize_method_a_current_f1_contribution_snapshot.py`; the current-snapshot and materializer tests; the historical/current CLI owners; canonical-five owner and profiles; `main.py`, production inputs/plots, Method-B, roadmap, `MEMORY.md`, `USER.md`, handoff, accepted closure evidence, and unrelated existing local worktree state. If another executable edit becomes necessary, stop `BLOCKED` for a newly reviewed contract; do not silently change the allowlist.

## 5. Deterministic regression gates

**Tests using the real unchanged producer and checker** must demonstrate:

- For a positive one-event sample with `x=55.83414573673653`, `s=summary([x], [1.0], [x])` is accepted by `validate_summary(s)` with `s["neff"]["w0"]["value"] == 1.0`. Explicitly prove on the actual interpreter that `x*x != x**2`; if this fixed sample is not discriminating on the discovered interpreter, choose and document a verified finite, positive alternate *in the test only*, without weakening any check or modifying physics.
- At least one nontrivial positive-support/reweighted example using real `summary()` and `validate_summary()`, including `w_rw` and `abs_b_rw`, passes. Include zero support and signed cancellation/negative contribution cases, preserving the distinction between nonnegative effective-statistics measures and signed baseline contributions.
- Deliberate changes to a `neff` value or its stored numerator/denominator continue to fail at `neff_moment_mismatch` (or the strict existing relevant invariant), with no new tolerance; invalid/negative/nonfinite moments and `N_eff > count` remain rejected.
- A persisted diagnostic validation through the current-snapshot and historical diagnostic test paths still executes all existing moment, joint `(t,phi,MM)`, marginal, tail, provenance and F.4 parent gates. Preserve previously accepted historical artifacts' definition and authority pins; do not rewrite accepted reference JSON.

Run available deterministic suites for:

```text
testing/test_pion_hgcer_method_a_reweighted_contribution_debug.py
testing/test_pion_hgcer_method_a_current_snapshot_contribution_debug.py
testing/test_materialize_method_a_current_f1_contribution_snapshot.py
```

and any focused tests necessary to ensure the existing current-snapshot preflight remains strict. Use the locally discovered interpreter per `docs/memory/TOOLS.md`. Tests may use synthetic fixtures and must not mock scientific checks to manufacture acceptance. Report exact pass/fail/skip counts and commands. Local tests **do not** constitute farm or production validation.

## 6. Preservation, integrity and review

- Confirm only the five allowlisted paths changed, no scientific algorithm or production path changed, and pre-existing unrelated local files, index and refs remain preserved. Report relevant status before/after and exact changed-file inventory. `git diff --check` and syntax compilation must pass.
- Run `<PYTHON> -B tools/update_memory_manifest.py --root . --write` then `--check`, `<PYTHON> -B tools/memory_bootstrap.py --root . --json`, and `<PYTHON> -B tools/check_memory_health.py --root .`. Record hard-gate results, CURRENT/MEMORY/CURRENT_HANDOFF byte counts, warnings and their material/nonblocking classification.
- CURRENT must record the failed detached materialization at `snapshot_failed:neff_moment_mismatch` (not as science-invalidating evidence), and preserve exactly one **push-stable substantive NEXT**: after ChatGPT actual-diff acceptance, user commit/push and pushed-state synchronization, perform **one user-controlled detached current-F.1 materialization farm gate** using the existing tracked owner; review its ZIP/manifest before any separate current-snapshot contribution preflight. The historical failed F.3 mismatch remains recorded, and canonical-five F.4/absolute-SIMC blockers stay binding. No new phase closure or Method-A promotion.
- Produce one complete repository-root temporary **actual-diff review bundle**, comprising tracked diffs plus complete `git diff --no-index /dev/null ...` for the new contract file. Report its exact name, byte count and SHA-256. Do not stage merely for review; no commit, push, or farm execution by Codex.
- After ChatGPT reviews and accepts the real diff, the **user alone** will run scoped Bash/WSL `git add` of the five paths, `git commit`, `git push origin test`. ChatGPT then independently verifies published HEAD, changed Git blobs, CURRENT/NEXT, and farm readiness **before** issuing a farm command.

## 7. Hard stops and expected state

Immediately stop at the first concrete blocker if source/worktree identity has drifted, unrelated files cannot be preserved, the exact parity repair cannot be implemented without touching a frozen file, regression/health/manifest checks fail, or the repair unexpectedly alters persisted scientific results, authority policy, floating tolerances, parent closure or accepted historical artifacts. Do not reopen F.1 archaeology or rerun a canonical-five owner.

Until new accepted evidence: E.8 is `ACTIVE`; detached materialization is `BLOCKED` at the reported farm gate; the parity implementation is `ACTIVE` until local checks and actual-diff review, then at most `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; canonical-five F.4 equivalence, final E.8/F.6.4, and absolute-SIMC interpretation remain `BLOCKED`; F.1 archaeology remains `DEFERRED`. **No farm command** is authorized from this contract before publication and pushed-state verification. The scientific objective—measure whether Method-A reweighting concentrates or redistributes the baseline pion contributions—does not change.
