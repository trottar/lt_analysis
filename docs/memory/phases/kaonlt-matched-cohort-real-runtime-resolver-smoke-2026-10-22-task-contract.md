# KaonLT — Prove the existing resolver path on actual farm Python before one scientific acquisition

**Task class:** One narrowly scoped source audit/implementation for the already-published matched-cohort acquisition owner. No new physics phase, no scientific change, no historical repair, and no authorization for a full farm rerun.

**Source audit anchor:** public `test` HEAD `bd308954c59ab10c528055ee111af20f18b0a04a` (October 22, 2026). Establish live local HEAD/branch/origin-test, index, dirty/untracked worktree and existing contract/bundles before edits; preserve unrelated state. Read `AGENTS.md`, then in order `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, `docs/memory/USER.md`; follow task-relevant MAINTENANCE/CODEX/COMMUNICATION/TOOLS, October 22 incident and current acquisition contract. Emit full health receipt.

## 1. Verified source call chain and actual blocker

In `testing/run_pion_local_t_matched_cohort_gate.py::execute_gate`, the farm already passes `isolation.import_identity(repo)` and `isolation.prepare_runtime_overlay(repo, worktree)`, then `isolation.probe_paths(worktree,env,paths)` before hitting `probe_private_callers(...)`. The existing `testing/run_e8_4_fix5_canonical_five_plot_gate.py::prepare_runtime_overlay()` imports the real installed `ltsep`, calls real `SetPath`, byte-copies the package, patches only copied `LTANAPATH`, imports the copied `ltsep`, checks import identity, and preserves installed bytes. This chain passed far enough for the reported `resolver_code_path_escape` in the later owner-only restricted resolver probe. Thus ordinary imports had already succeeded for that attempt; **the failing step is redundant newly added source-loader/path policy**, not proof that uproot was missing.

Four user-supplied farm attempts on October 22 failed before scientific evidence. The latest at `f266045e9...` was `resolver_code_path_escape` when importing real `uproot` from an interpreter `lib64/python3.9/site-packages` alias. Current `bd308954c...` repairs that exact alias handling, but has **no farm execution evidence**; 3.12 host fixture with a 3.9 path shape is not 3.9 execution.

**Main audit conclusion:** The owner has no published **no-child diagnostic CLI** to evaluate its actual real-Python-3.9 `probe_private_callers` import/caller semantics independently of the expensive acquisition and unrelated stage input scans. Avoid another premature readiness approval.

## 2. Minimal source change — one no-child resolver check

Keep the current published alias fix and the production capture path. Add an explicit `--resolver-smoke` entry to the **same tracked owner** (no second production driver), with source identity and exact relevant path checks, that:

- Accepts the published `--source-commit` and environment-configured established repo (and any existing owner path arguments that are necessary); rejects unrelated capture/offline argument combinations. Its output is a bounded JSON receipt to stdout and failure stderr, **not a ZIP, diagnostic packet or scientific result**.
- Uses the same real `probe_private_callers` implementation, the existing `ltsep` import and copied LTANAPATH overlay, baseline/private `Root` no-mode and Plot calls, inherited cwd, actual unchanged caller list, and protected output-path comparison that `execute_gate` uses. No mocks, static config substitutions, alternate interpreter, or separate code path which might accidentally pass while the real owner fails.
- Stops **before** `stage_input_files`, external input reading, manifest/observer installation, `run_owned_child`, collection, rendering or bundle delivery. It must not touch shared analysis outputs or create owner artifact destinations, and must not import/execute scientific production modules.
- Permits only the existing narrowly owned *disposable* worktree, copied ltsep overlay, and private scratch output under existing isolation context; checks original installed/copied package, worktree and external path/links preservation on both success and failure, with no external writes, source pollution, or stale disposable directories. Fail closed if the current primitives cannot enforce this scope.
- Provides the resolved Python executable/version and **actual import/source identities** for `ltsep`, `uproot` and dependent modules as captured by the existing resolver probe, with no private workstation/farm paths written to tracked memory. The receipt must state `runtime_not_science=true`, `analysis_started=false`, `archive_created=false`, `scientifically_accepted=false`.
- Preserve strict no-mode/Plot mismatch, external-link/input authority, exact source commit, and read/write provenance rules. Do **not** merely suppress the current resolver read hook or remove the path protections. Do not alter the independent frozen resolver, frozen helpers, launchers, installed ltsep or science.

**Do not redesign import auditing as speculative preventative work.** If a source-confirmed issue in the current alias code prevents the smoke from being safe, correct *that exact invariant* within the same tightly constrained source scope and add a regression. No catch-all exceptions and no silent fallback.

## 3. Deterministic evidence required before handoff

- Tests prove the new smoke route invokes the **actual** resolver subprocess and the same caller list as `execute_gate`; normal `ltsep` import checks and wrong-mode/default private output mismatch still fail.
- Tests prove hard stop on unknown commit, package/read guard exception, mutated original/copied package, wrong cwd/Root, unexpected external path, symlink alias outside allowed root, modified index or source, and temporary output collision.
- On both successful/failed synthetic smoke: assert **no** child, ROOT source import, stage copy, input artifact hash/traversal, ZIP, log, collector or external writes; verify cleanup/preservation.
- Run existing owner resolver fixture tests, 3.9-shape alias tests, relevant historical no-child/negative tests, `py_compile`, `git diff --check`, manifest check, standard memory health/bootstrap. Classify older launcher-pin historical failures separately. Report that **real Python 3.9 and real installed uproot are NOT VERIFIED** until the smoke is actually run by the user on the farm.
- Audit source references for preflight path equivalence: identify any difference between the smoke and `execute_gate` up to the last resolver return. A smoke that uses less stringent input/caller checks cannot be advertised as readiness for the actual resolver gate.

## 4. Allowlist, memory and handoff

Allowed source:
- `testing/run_pion_local_t_matched_cohort_gate.py`
- `testing/test_run_pion_local_t_matched_cohort_gate.py`

Allowed documentation:
- `docs/memory/CURRENT.md` (replace obsolete wording, keep <16384 bytes, preserve science NEXT)
- `docs/memory/investigations/kaonlt-matched-cohort-farm-readiness-failure-2026-10-22.md` (one short addendum on source audit and smoke boundary, no false validation)
- `docs/memory/manifest.json` (regenerate)
- `docs/memory/phases/kaonlt-matched-cohort-real-runtime-resolver-smoke-2026-10-22-task-contract.md` (this one standalone contract)

No other source, memory, production, ROOT, SIMC, model, binning, cuts, weights, normalizations, artifacts or external environment changes.

Codex does local implementation and deterministic checks only, produces a single complete actual-diff review bundle (including full new contract) for independent ChatGPT review, and must **not** commit, push, or run JLab. Then user publishes accepted source and ChatGPT verifies pushed state. Only then one user-operated safe resolver smoke; if it passes on **real Python 3.9**, review its receipt and decide readiness for the existing **single** matched-cohort science owner. If it fails, stop at the first real invariant without launching the science child. Never claim that a smoke pass itself is full analysis or physics validation.

**Scientific NEXT remains fresh Q4p4W2p74 Center/Left lowe t1/t2/t3 same-run A/B and signed-event packet evidence followed by seven-strength C and six-parent missing-mass/variance comparison under Lambda.**
