# KaonLT Check 1 — Scientific-identity failure diagnostics (November 8, 2026)

## Decision and authority

**Task class:** narrow, diagnostic-only source correction following a failed Check 1 owner gate. **Work status:** `ACTIVE` for implementation; Check 1 scientific acquisition remains `BLOCKED` pending the original invariant.

Authoritative repository: `trottar/lt_analysis`, branch `test`. Audited live HEAD on contract preparation: `912ca31d0afdd7b9757395d205d7c54bc49018e8`. The current Git worktree and index are **NOT VERIFIED** by ChatGPT; Codex must establish them before editing. The farm checkout previously reported two explicitly permitted local exceptions: modified `src/models/xmodel_kaon_pl.f` and untracked `src/kaon/functions/Q4p4W2p74.model`. Do not touch either.

Read tracked root `AGENTS.md`, then `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, and `docs/memory/USER.md` **in exactly that order**; then `docs/memory/MAINTENANCE.md`, `docs/memory/CODEX.md`, `docs/memory/COMMUNICATION.md`, and the existing `docs/memory/investigations/full-t-acceptance-edge-check1-2026-11-05.md`. Reconcile relevant intervening source changes if HEAD has moved; stop for a scientific, reference, ownership, or gate change. Do not reset, stash, clean, stage, commit, push, or execute on the Jefferson Lab farm.

### Fixed failed evidence

The 2026-11-07 owner attempt `KaonLT_Check1_Acquisition_Q4p4W2p74_20261107T163552Z-3952780` used the reviewed source SHA above and frozen reference ZIP `KaonLT_Check1_Reference_Inputs_20261105.zip`, SHA-256 `02ff71998fd51d5e56f07484bb1c71185f78ed51d6ce409b0d34154bda7a6d7a`. Its `analysis_returncode=0`, `analysis_completed=true`, but `failed_stage=verify_artifacts`, `artifact_verification_completed=false`, `collection_completed=false`, `zip_verification_completed=false`, and both scientific acceptance fields are false. The **first failed invariant** is exactly:

`check1:scientific_identity:Left:lowe.cells[0].fit.pion_component_alignment.applied_components.pi_delta.parent_alignment_hash`

The previous `check1:destination_exists` producer failure did not recur. The November 7 failed artifacts and log are diagnosis-only, not accepted science. Preserve the original frozen reference and all accepted previous phases. The current log does not print the generated 64-character hash or the actual/expected operands of the failed comparison. Do not guess the underlying scientific cause.

## Exact outcome

When the existing Check 1 owner encounters a **scientific-identity mismatch**, it must retain the **identical failing invariant** while also capturing a compact, structured, deterministic **actual-versus-frozen-reference diagnostic** in the existing owner gate-status JSON and emitting a concise matching stderr line. This must occur **before the private analysis output is removed** and before any downstream collection, ZIP, publication, or scientific interpretation. A passing verification must behave exactly as before.

The diagnostic should identify what differed, not merely restate a field path. The first target is the observed Left/lowe cell-0 `pi_delta.parent_alignment_hash`, but the implementation must be reusable for other existing Check 1 identity mismatch paths without selecting or dropping scientific fields opportunistically.

## Scope and edit allowlist

Modify only the following existing files (bounded additions; first verify they exist):

- `src/utility/check1_acceptance_evidence.py` — pure-Python, ROOT-independent helpers for first-difference reporting and compact alignment identity records, if necessary. **Do not modify the behavior or tolerance** of `scientific_projection`, `validate_coordinate_contract`, or `compare_science`.
- `testing/run_pion_t_template_comparison_gate.py` — capture failure context around the existing `module.compare_science(...)` invocation in `_verify_check1_packets`, feed it to the existing owner status machinery, print bounded stderr diagnostics, and re-raise the **original exception unchanged**. Extend the existing Check 1 call-chain with only the minimum status/callback plumbing needed. Leave passing execution, preflight, artifact inventory, collector and release paths unchanged.
- `testing/test_check1_acceptance_evidence.py` and, only if needed for the owner-level failure receipt, `testing/test_run_pion_t_template_comparison_gate.py` — deterministic regression tests using existing fixtures and the real serializer/owner logic without ROOT or farm execution.
- `docs/memory/CURRENT.md` and the existing Check 1 investigation — concise factual correction to the dated failed-attempt evidence and the next scientific gate. Preserve other status, closed phases, roadmap, Method A/B and Outcome B.
- `docs/memory/manifest.json` — regenerate/check with the tracked memory tool after the changes.
- **This one new task contract** at `docs/memory/phases/check1-scientific-identity-failure-diagnostics-2026-11-08-task-contract.md`.

Do not change `src/cuts/pion_component_fits.py` (the actual alignment-hash producer), `src/cuts/pion_t_phi_template_comparison.py`, `src/cuts/pion_t_bin_parents.py`, fitting/configuration/cuts/templates, `main.py`, any yield/normalization/SIMC logic, reference inputs, profile thresholds, collector, installation/launcher, or any other file. If a true fix requires another path, **stop** and report the missing scope rather than expand this contract.

## Required failure diagnostic

For the first identity mismatch, include bounded fields sufficient to answer:

1. **Context:** source SHA; frozen reference SHA; setting, epsilon, cell index and `t_index`/`phi_index` when available; original failure path and category; current/expected type and value (hash strings in full).
2. **Provenance chain:** for the affected cell, inspect both the generated child JSON and the frozen child JSON, plus `children[setting:epsilon].parents[t_index].baseline_A.pion_component_alignment` in each. Record the parent alignment's hash, child alignment's hash, and the affected applied component's inherited `parent_alignment_hash`, on both sides. Mark missing or inconsistent links explicitly.
3. **Alignment selection:** for parent and child on both sides, record `accepted`, `source`, `analysis_scope`, `persistence_status`, `persistence_rejection_reasons`, relevant applied-component source and fallback reason, without implying any state is scientifically valid merely because it was serialized.
4. **Exact parent-hash operands:** compare `parent_setting_key`, `common_setting_shift_gev`, `resolved_configuration_hash`, and each `components.{pi_n,pi_delta,pi_sidis}` entry restricted to `residual_shift_gev`, `window_expansion_gev`, `effective_windows`, `immutable_source_hist_identifier`, and `immutable_source_hist_checksum`. These are the operands of `_alignment_parent_hash` in the **unchanged** producer. Report the **first differing nested operand** with both values; report a limited additional list or count if useful. For large identifiers/windows, include short typed summaries and stable SHA-256 digests of the *complete* canonicalized values, with enough exact differing scalar/short-window values to explain the difference. Do not silently round or tolerance-match hash preimage inputs.
5. **Consistency:** when a parent record is available, verify in diagnostic-only mode whether its stored parent hash equals the hash of the producer-defined parent preimage (use identical JSON normalization in a pure helper, or report `not_verified` if faithful reproduction is impossible). This is explanatory information only; **never substitute a recomputed hash into scientific output or acceptance**.
6. **Disposition:** same `failed_stage=verify_artifacts`, same `first_failed_invariant` and `failure_reason`, unchanged `evidence_admissible=false`, `scientifically_accepted=false`, no release/ZIP/collection. Ensure the diagnostic is visible in the existing owner gate-status JSON even if the private worktree is cleaned. Bound serialized diagnostics (e.g., up to 64 KiB, fixed first-difference ordering, no uncontrolled full-array/histogram dumps).

Do not opportunistically add untrusted arbitrary paths, private file content, full event records, or unbounded traceback dumps to the receipt. Do not replace the first owner failure with an exception from the diagnostic itself. If diagnostic construction fails, retain the original invariant and a bounded `diagnostic_error` secondary note; never turn a failed scientific check into success.

## Source seams verified on preparation

- `src/utility/check1_acceptance_evidence.py::compare_science` recursively calls `require` and presently emits only `scientific_identity:<path>` (or a corresponding keys/length/value error), with no operand details. Its float closeness remains `rel_tol=1e-8, abs_tol=1e-12`, **unchanged**.
- `testing/run_pion_t_template_comparison_gate.py::_verify_check1_packets` reads the newly generated child JSON and `reference["children"][key]`, validates the coordinate contract, and calls `module.compare_science(module.scientific_projection(child), module.scientific_projection(reference["children"][key]), key)` **before** reading the separate parent comparison and validating the Check 1 packets. Add failure diagnostics there without moving or skipping any gates.
- The owner calls `_verify_check1_packets` inside `verify_check1_artifacts` after completed ROOT analysis, and catches the exception at `execute_gate` to populate the status. The owner retains the first failure even through preservation/cleanup. Use its existing status-update pattern rather than implementing a second gate owner or artifact publication path.
- The existing reference ZIP's Left/lowe `global/Left_kaon_Q4p4W2p74_lowe_pion-t-phi-template-comparison.json` has `cells[0].t_index=0`, `phi_index=0`, with child `parent_alignment_hash=c03775871b5033a2d9682e174293ecefe599ebefb1a2e9a10cbc102dd223c82a`. Its `parents[0].baseline_A.pion_component_alignment.parent_alignment_hash` has the **same hash**. The parent was `accepted=true`, `source=bin_local_scan`, with `common_setting_shift_gev=0.0012744878163235551`; the child `accepted=false`, `source=current_common_shift_fallback`, with `common_setting_shift_gev=0.0` and reason `parent alignment unavailable, stale, incompatible, or invalid`. Preserve these distinctions: **the child's fallback state is not evidence that the parent was unaccepted**.
- `_alignment_parent_hash` in `src/cuts/pion_component_fits.py` hashes setting, common shift, resolved configuration hash and each component's specified shift/window/template identifiers and checksums via canonical JSON/SHA-256. Its implementation is frozen.

## Required deterministic tests

1. A controlled owner/child fixture whose `pi_delta.parent_alignment_hash` differs from its frozen counterpart must fail at the **same exact invariant path**, produce both full hash values, identify the corresponding parent record and first differing parent-hash input, and serialize those details into the **owner's** failed status. Compare against the provided immutable October reference structure where practical without importing ROOT.
2. A changed parent hash with identical child scalar values otherwise must remain a hard failure. A mismatch due only to a different component template checksum or a float/shift operand must be identified at the correct nested field. Include a deliberate parent/child lineage inconsistency case, and a missing-parent-record case that reports insufficiency without masking the failure.
3. On an unrelated scientific difference (e.g., a float beyond tolerance, keys/length mismatch, different setting), retain the same strict comparator semantics. No data-dependent exclusions or fallback to a weaker comparison.
4. A success case must not create diagnostics or change status/artifacts. In an injected diagnostic-construction failure, retain the original first_failed_invariant and failed owner status.
5. Check the existing Check 1 packet and owner tests, at least relevant standard-library and available NumPy deterministic suites; Python 3.9 grammar/compile for changed files; tracked memory manifest/health. Report exact interpreter, tests, skips and reasons. No synthetic test may be represented as farm validation.

## Review / output contract

Produce a **single review bundle** containing a changed-path inventory, full tracked diffs (including this contract), complete added-file contents, before/after diagnostics for the controlled failing fixture, run logs and test/skip summary, memory manifest checks, and a clear hard-stop statement. No auto-commit/staging/push, no farm run, no new scientific plots, no generated acquisition ZIP.

Return `BLOCKED` if actual repository state, scientific ownership, external reference, or a required edit falls outside this contract. Otherwise return `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING` for Codex's implementation claim, **subject to ChatGPT's independent actual-diff review**; do not claim `SOURCE REVIEWED` until that review has passed.

**After ChatGPT source acceptance:** user alone commits/pushes; ChatGPT synchronizes exact published blobs/state; only then may a *narrow, evidence-preserving diagnostic execution* be considered under separately reviewed farm authority. A second full 2-hour acquisition is **not** authorized by this contract.
