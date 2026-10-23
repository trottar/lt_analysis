# KaonLT — Python 3.9 ctypes import bootstrap failure: one owner-only repair

**Work state: BLOCKED.** The science target remains fresh Q4p4W2p74 Center/Left lowe canonical t1/t2/t3 same-run A/B + signed-event packets, then seven-strength C replay and six-parent MM/variance diagnostics. No new accepted science and no full acquisition authorized by this task.

**Source anchor:** live public `test` = `abff8a49ead8c185b601af7e00fa5a89903bde52` on October 22, 2026. Before edits, verify actual local `test`/HEAD/`origin/test`, index and complete dirty/untracked state. Stop on overlapping work or unexpected identity. Read root `AGENTS.md`, then CURRENT, MEMORY, CURRENT_HANDOFF, USER in exact mandated order; read MAINTENANCE, CODEX, the current October 22 failure incident, and current resolver-smoke contract. Emit full health receipt and invoke the science-displacement alarm. Do not create a second support phase.

## Exact observed failure and source cause

User's **first** real `--resolver-smoke` on farm at `abff8a49e` failed before child analysis:

```
resolver_subprocess_failed:src/setup/findEffectiveCharge.py
... ltsep -> uproot -> numpy.core._internal -> /usr/lib64/python3.9/ctypes/__init__.py
... ctypes/__init__.py:464 -> ctypes/__init__.py:374
RuntimeError: resolver_mutation_attempt:ctypes.dlopen:(None,)
```

**RUNTIME VERIFIED only for supplied traceback:** this is Python 3.9's standard-library module-level `pythonapi = PyDLL(None)` initialization, not evidence of an analysis mutation. The owner `testing/run_pion_local_t_matched_cohort_gate.py::RESOLVER_PROBE` rejects **all** `ctypes.dlopen` in its `readonly` audit hook, installed before importing actual `ltsep`; the synthetic tests imported fake uproot/numpy and did not exercise Python 3.9 ctypes startup. The failure is ours, not the user's command. Status is FAILED / science-inadmissible, and the actual installed full dependency closure remains NOT VERIFIED.

## Narrow fix: safely bootstrap the genuine stdlib dependency

Change only the existing owner's `RESOLVER_PROBE` and its test suite. Prefer a **two-stage import gate**, not an arbitrary ctypes allowlist:

1. **Before** installing the full code/data-restrictive import hook, verify the interpreter-selected `ctypes` module's import spec resolves to the interpreter's trusted stdlib/platstdlib root (including a verified interpreter-specific lexical alias). Reject module shadowing via `PYTHONPATH`/untrusted filesystem, nonstandard loader, symlink escape or unexpected source.
2. Install an early **write/process/network/mutation-denial audit hook** *before* actually importing verified stdlib `ctypes` (and whichever explicitly justified stdlib import dependencies it requires). This narrow bootstrap may read verified standard-library code; it must never read arbitrary scientific data, import `ltsep`, launch child programs or mutate any files. The bootstrap hook stays installed when the full existing readonly hook is added.
3. Import **only** that genuine `ctypes` stdlib module under the bootstrap guard. Verify its actual `__file__`, source/code location and interpreter provenance; do not create arbitrary `ctypes.CDLL` instances except those occurring intrinsically during Python 3.9 module import.
4. Install the unchanged/full (or stronger) resolver `readonly` and `Finder` controls *after* this bootstrap. Keep the unconditional rejection of **subsequent** `ctypes.dlopen`, `ctypes.dlsym`, arbitrary native library loads and ROOT, as well as external scientific reads, writes, network and subprocess activity. Strictly preserve current source-module provenance/hash verification, actual baseline/private no-mode/Plot `Root` calls, copied/installed ltsep preservation, symlink authorities, and failure barriers. If a broad preimport violates required read-only/provenance guarantees, STOP and explain the first constraint; do not silently weaken it.
5. Audit adjacent startup operations and the code path for known Python 3.9 differences; fix ONLY source-demonstrated contradictions in this same owner. Do not add another independent subsystem or suppress failures.

**Explicit anti-regression:** A naked `if event == "ctypes.dlopen" and args == (None,): return` is **not accepted**; it permits arbitrary post-import access to the process symbol table. Nor is globally disabling audits or bypassing `probe_private_callers`.

## Deterministic verification

- Reproduce the **exact old hook behavior** on controlled `ctypes.dlopen(None)` event and assert denial prior to the fix. Demonstrate the corrected real imported stdlib `ctypes` succeeds with the bootstrap guard active. Python 3.9 may NOT be available locally: report truthfully; do not call Python 3.12 a Python 3.9 test.
- Probe actual code and test fixture modules through the existing full `probe_resolver` and `resolver_preflight` functions. Confirm real import of `ctypes` before the full audit; all callers/default/Plot path identities remain checked.
- Negative tests: shadowed `ctypes`, untrusted site-packages masquerading as stdlib, second `ctypes.CDLL(None)` after initialization, other shared-library names, `ctypes.dlsym`, filesystem reads of physics `.root`/`.csv`/`.json`, external writes/mkdir/symlinks, subprocess/ROOT, input/config aliases and package mutations must fail closed. Test no scientific child, no input staging, no collection/ZIP and preservation on both smoke success and failure.
- Run owner, frozen-resolver, relevant no-child/negative suites, py_compile, diff whitespace, manifest write/check, ordinary memory health/bootstrap; report unrelated historical launcher-pin failures separately. Do not claim farm compatibility, full ROOT or science validation.
- Verify that smoke and acquisition share **the identical** fixed resolver path and complete source identity gate.

## Allowlisted files / publication workflow

- `testing/run_pion_local_t_matched_cohort_gate.py`
- `testing/test_run_pion_local_t_matched_cohort_gate.py`
- `docs/memory/CURRENT.md` (replace superseded words; stay below 16,384 hard bytes)
- `docs/memory/investigations/kaonlt-matched-cohort-farm-readiness-failure-2026-10-22.md` (short factual fifth-failure addendum and exact first invariant)
- `docs/memory/manifest.json` (regenerate)
- **This new contract:** `docs/memory/phases/kaonlt-python39-ctypes-bootstrap-fix-2026-10-22-task-contract.md`

No production `src`, frozen resolver, farm_env, launchers, settings, Method A/B, templates, kaon/pion cuts, yields or analysis scientific source. Preserve source and external state; provide ONE complete actual diff in the repository root including full addition of this contract, tests and receipt. Codex does NOT commit/push or use the JLab farm. After independent actual-diff review, user-only publication, ChatGPT pushed-state synchronization, the one permitted next farm gate is a **single no-child resolver smoke** using the actual interpreter. No science run before receiving and reviewing its actual failure/success evidence. This task must not be turned into another speculative infrastructure hardening cycle.
