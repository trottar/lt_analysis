# KaonLT — Unblock matched-cohort resolver dependency imports without weakening read/write provenance

**Task class:** One narrowly evidenced source-preflight repair in the *existing* matched-cohort acquisition, not a new phase, model change, or general resolver redesign.

**Live source audit anchor:** `test` at `a14e740ab68531aca1c60c48d7a73e13e8922aff` on 2026-10-22. Codex must verify its **own** local `test`, `HEAD`, `origin/test`, staged/unstaged/untracked state, and any review bundles, preserving all unrelated local content. Stop if those conflict with this one repair.

**Immediate runtime evidence (user's farm terminal):** The reviewed owner invoked at the audit SHA failed before running its analysis child:

```
resolver_subprocess_failed:src/setup/findEffectiveCharge.py
  copied runtime-python/ltsep/__init__.py
    -> copied runtime-python/ltsep/ltsep.py
    -> import uproot
    -> importlib source load
    -> PROBE.readonly(...)
RuntimeError: resolver_data_read_forbidden
```

Owner status: `KaonLT_matched_cohort_Q4p4W2p74_20261022-155418-gate-status.json`. **No ZIP, fresh A/B/C event packets, six-parent MM plot or scientific acceptance.** `RUNTIME VERIFIED` is confined to the user-supplied exception and command; unprovided status details and actual installed-package import closure remain **NOT VERIFIED**.

**Retained scientific NEXT:** One owned Q4p4W2p74 Center/Left low-epsilon canonical t1/t2/t3 matched acquisition of same-run A/B and signed-event packets, followed by existing offline seven-strength C fit/replay and Lambda-region full-MM/variance comparison. Do not replace it with historical C, select a lambda optimum, promote production, or start a general infrastructure phase.

## 1. Required startup and first failure audit

Read repository-root `AGENTS.md`, then `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, `docs/memory/USER.md` **in precisely that order**. Then read `docs/memory/MAINTENANCE.md`, `CODEX.md`, `TOOLS.md`, `COMMUNICATION.md`, the active October 21 matched-cohort attempt-isolation contract and October 22 link repairs. Emit the full KaonLT health receipt and science-displacement re-anchor. Audit *actual* source before coding.

**SOURCE VERIFIED at the audit anchor:**
- `testing/run_pion_local_t_matched_cohort_gate.py::probe_private_callers()` invokes `resolver.probe(worktree, caller, modes, env, cwd)` and `resolver.probe(repo, ...)` for actual no-mode and Plot resolver callers **before** `run_owned_child`.
- `testing/attest_pion_chain_ltsep_default_root.py::PROBE` inserts an audit hook, then imports `ltsep` and `farm_env.ltsep_paths`. Its allowed read locations are the copied/installed `ltsep` package, repository `farm_env`, and stdlib. The actual `ltsep/__init__.py -> ltsep.py -> uproot` dependency import attempts source-file reads from outside those roots and fails `resolver_data_read_forbidden`.
- The same frozen `PROBE` also prohibits direct NumPy imports alongside ROOT/root_numpy, which could trigger the next failure after the first issue is addressed. **Do not claim that simply adding site-packages to `allowed_reads` is sufficient.**
- The delegate is a separately owned, source-frozen attestation interface with independent synthetic tests; it must not be silently rewritten or globally monkeypatched. `farm_env/ltsep_paths.py::create_ltsep_root()` uses the real installed/copied `ltsep.Root`, so a fabricated fake Root or mode-independent config reading cannot attest actual default versus Plot behavior.

Read `testing/test_attest_pion_chain_ltsep_default_root.py` negative tests, actual matched-owner tests, `isolation.prepare_runtime_overlay()`, `probe_paths()`, and `reached_resolver_callers()`. Identify precisely which resolver safety invariants are essential to the matched owner.

## 2. Small owner-specific repair; safe code imports are not data reads

Modify only the matched owner to give its actual default/Plot `ltsep.Root` probes an **owner-scoped, testable, read-only** dependency import policy. Do not edit the separately frozen `testing/attest_pion_chain_ltsep_default_root.py` or `farm_env/` source, and do not patch its global `PROBE` string or monkeypatch its module at runtime.

Implement an explicit locally owned subprocess probe or a clearly scoped wrapper using the same real copied/installed `ltsep.Root` and caller context. Preserve or strengthen:
- Correct caller path, mode (`[]` versus `Plot_LTSep`), inherited cwd, source-module identity, exact copied-ltsep import identity and matching unique `PATH_TO_DIR` config.
- Every stable input field, original-vs-private `OUTPATH` comparison, baseline-versus-child equivalence and physical private output root checks.
- The read-only hook's denial of **scientific-data or configuration reads outside the explicitly trusted owned/copied source/config roots**, including ROOT, `.root`, `.csv`, `.json` data, arbitrary external input paths and unexpected package-data files. A scientific data read must remain an explicit hard failure, even if the path lies within a Python site-packages directory.
- Denial of attempted writes, processes, symlink mutations, arbitrary command execution and path escapes; original and copied package/ordinary checkout preservation; no scientific module execution.
- Strict source-only import behavior (do not accept preexisting `.pyc` as an unreviewed code authority). The *minimum code import closure* needed by the actual ltsep (e.g., `uproot` and possibly NumPy as an import-time library dependency) may read Python source/approved binary extension code under actual, checked interpreter dependency roots. Do **not** broadly whitelist site-packages reads or remove all NumPy/source safeguards. If NumPy is unavoidably a dependency, distinguish that controlled import from running analysis code, and retain ROOT/ROOT binding rejection.
- No data reads at all from external scientific `VOLATILEPATH`, `UTILPATH`, `PARAMPATH`, `CUTPATH`, `SIMCPATH`, etc., in the resolver-only probe. Required actual scientific inputs are still inventoried independently and strictly before the child.

Inspect testable interpreter dependency locations through the executing Python and verified import machinery; do not hardcode a farm interpreter path, guess package locations or assume `uproot` alone is the entire closure. Fail closed when a needed import is outside a defensible code-only boundary, or requires an unsafe read/write, and report its exact denied event/path in bounded status traceback.

A minimal adaptation that cannot prove the *real* no-mode `ltsep.Root` behavior, silently skips the probe, substitutes static configuration for `Root`, changes frozen helpers, or broadly disables audit controls is **not accepted**. If preserving these safeguards cannot be achieved within this scope, stop `BLOCKED` and identify the first required extra authority rather than claiming readiness.

## 3. Specific negative and positive tests before acceptance

Extend `testing/test_run_pion_local_t_matched_cohort_gate.py` (and only it, not frozen attestation suites) with realistic isolated dependency fixtures. Do not use mock return values to claim the guarded import path works.

1. **Reproduce the old failure** with a synthetic copied `ltsep/__init__.py -> ltsep.py -> uproot` dependency chain whose `uproot` source is outside the copied-ltsep and stdlib roots. Check exact `resolver_data_read_forbidden` before the fix.
2. **Positive import closure:** a controlled source-only third-party import (`uproot` and, if needed, `numpy`) succeeds without reading `.pyc`, physics data or external configuration, and returns the actual fixture `Root` path fields for both `[]` and `Plot_LTSep`, retaining distinct caller/cwd and copied-package identity. Confirm the entire `probe_private_callers()` path, not just a helper. No analysis child need run to verify this narrow positive gate.
3. **Negative data and mutation:** attempts to read scientific `.root`, `.csv`, arbitrary `.json`, outside-allowed code, or an unexpected package-data file fail; attempts to write a file, mkdir, unlink/symlink/rename, subprocess execution, import ROOT/root_numpy, package-byte mutation or change package/config realpath fail. Include negative tests for a third-party *code* file trying to read an external physics file. Source `.pyc` loading must not bypass source identity.
4. **Baseline/private equivalence:** altered no-mode Root OUTPATH, mismatched plot mode, changed stable input roles, missing/multiple config matches and unexpected imported package/module must stop before child.
5. **Owner isolation:** on a failed probe `analysis_started=false`, no child, no ZIP or collection, with a bounded traceback and all applicable preservation checks; on a successful synthetic probe the owner advances to the previously established protected child gate. Do not fake full ROOT, actual analysis, or collector acceptance.
6. Re-run unchanged existing matched-owner, event capture/replay and relevant resolver tests where supported; report suites unexecuted for missing NumPy rather than asserting PASS. `py_compile`, `git diff --check`, manifest write/check, ordinary memory health/bootstrap and complete diff review required.

**Finite scope:** Confirm the actual first dependency-load failure is resolved in synthetic execution; enumerate remaining known probe imports/caller guards without launching another broad speculative audit. No new platform or global packaging infrastructure.

## 4. Allowed changes and checkpoint

**Only source:**
- `testing/run_pion_local_t_matched_cohort_gate.py`
- `testing/test_run_pion_local_t_matched_cohort_gate.py`

**Only memory:**
- `docs/memory/CURRENT.md`
- `docs/memory/investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md`
- `docs/memory/manifest.json` (regenerated)
- `docs/memory/phases/kaonlt-matched-cohort-resolver-import-preflight-fix-2026-10-22-task-contract.md` (this contract)

Do **not** change `testing/attest_pion_chain_ltsep_default_root.py`, old attestation tests, production `src/`, launchers, external package install, copied-ltsep source beyond the existing authorized config-value-span overlay, physics models, cuts/templates/priors, event selection, signal definitions, yields, SIMC normalization or archival scientific evidence.

CURRENT is close to its **16,384-byte hard limit**. Replace/compact superseded local preflight paragraphs and NEXT in place; never append past hard size. Record this `155418` failure as *failed/inadmissible* and note exact remaining dependency uncertainty. Preserve the scientific NEXT. No personal farm/workstation paths in tracked memory.

After changes run the existing manifest write/check, memory health in ordinary nonstrict form and bootstrap; report CURRENT/MEMORY/CURRENT_HANDOFF bytes and warnings. Nonblocking warnings must not create another memory-maintenance cycle.

## 5. Handoff and evidence boundary

Codex: report exact changed paths, first invariant, realistic positive/negative test results, remaining unverified farm/third-party import behavior, protection of source and memory; create **one complete review diff in repo root**, including full new contract addition (`git diff --no-index /dev/null contract` when untracked), without staging for review. Do not commit/push, alter farm filesystem, run JLab, or deliver artifacts as accepted.

Sequence: **Codex narrow fix → ChatGPT actual-diff review → user commit/push → ChatGPT pushed-state synchronization → ONE user-controlled guarded farm acquisition → owner/provenance/ZIP review → fresh offline seven-strength C and six-parent signed-event MM/variance science.** Any failed farm gate stops without automatic retry or packaging. A positive source-only fake importer does **not** establish real farm dependency loading. No farm command until the preceding gates pass.
