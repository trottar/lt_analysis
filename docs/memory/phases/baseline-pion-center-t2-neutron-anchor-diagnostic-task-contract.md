# KaonLT — Center/lowe canonical-t2 baseline neutron-anchor first-failure diagnostic

Status: `ACTIVE` (source-changing diagnostic-only contract; no farm acceptance).

## Starting authority and evidence

- Starting remote `test` HEAD: `c99f1b8981903c36ffca26d5926e3afd54c39f40`. Before editing, independently confirm the actual local branch, HEAD, origin, complete staged/unstaged/untracked state, and read tracked root `AGENTS.md`, then `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, `docs/memory/USER.md`, in that order. Follow `docs/memory/MAINTENANCE.md` and current investigation references. Stop if HEAD differs or an unexpected local change touches a task path.
- Pushed `docs/memory/CURRENT.md` owns NEXT: diagnose the original Q4p4W2p74 / Center / low-epsilon / canonical t2 pion-control staged `B_n=0` before revising any baseline physics.
- The **available user-supplied 67-page Q4p4 Center/lowe PDF** is setting-wide, not a canonical per-t parent; its setting-wide pion-control neutron LS numerator is positive (`+87.77`, amplitude `0.6646`). The separately reported August-26 canonical t2 `B_n=0`, weight near 183, and bad model/application closure are **MEMORY ONLY** historical evidence. Do **not** invent a second PDF, identical-input regression, original solver numerator, or farm validation.
- Current source, not a reconstructed history, defines the current capture. `src/cuts/rand_sub.py::_build_authoritative_pion_control_source_cache` retains four classes of signed per-event records in `hist["_authoritative_pion_control_source_cache"]["by_t"]`; `src/cuts/pion_t_bin_parents.py::build_setting_t_bin_pion_parents` obtains current canonical parent `fit_result`; `fit_result["H_pion_fit_step_overlays"]` retains per-component/pass `fit_diagnostics` (numerator, denominator and arrays), and `fit_result["diagnostics"]["pion"]["staged_pass_history"]` retains `success/message`. The scalar parent summary alone is insufficient.

## Single objective

Add **read-only, diagnostic-only, same-run extraction** of the exact existing canonical t2 neutron-anchor input and staged solver decision, for Q4p4W2p74 / Center / low epsilon (plus immediately adjacent canonical t1 only for internal comparison if it already exists in the same parent collection). Emit one small, directly attributable JSON artifact after parent fitting, without refitting, event re-selection, template reshaping, production acceptance changes, or extra runtime input populations. If current source does not reproduce the historical `B_n=0`, report that fact and provenance mismatch; do not manipulate it to reproduce the old result.

The output is a **first-failure diagnostic**, not a scientific correction or an E.8/F.6 gate. It should decide, when the live fitted arrays permit, whether a zero staged neutron amplitude came from a nonpositive finite weighted numerator with positive support denominator, or from an early invalid-template/no-bin/no-support failure. Explicitly label cases where the evidence is insufficient.

## Allowed paths

1. This single standalone task contract: `docs/memory/phases/baseline-pion-center-t2-neutron-anchor-diagnostic-task-contract.md` (new).
2. `src/cuts/pion_t_bin_parents.py` — one minimal post-fit observation hook, no fitting changes.
3. `src/cuts/pion_control_anchor_audit.py` — new pure-Python read-only audit builder/JSON output helper; avoid importing ROOT at module import time.
4. `testing/test_pion_control_anchor_audit.py` — new deterministic unit tests with fake histograms and scalar records.
5. `docs/memory/CURRENT.md` — only narrowly necessary active-state/NEXT wording, preserve accepted closures.
6. `docs/memory/manifest.json` — regenerate using tracked repository procedure.

If a necessary input cannot be reached from these allowlisted paths **without changing scientific/production interfaces**, stop `BLOCKED` and report the exact missing producer/consumer boundary. Do not silently expand the allowlist or edit `src/main.py`, `src/cuts/rand_sub.py`, `src/cuts/pion_component_fits.py`, `src/utility/background_config.py`, `src/utility/pion_component_shapes.py`, existing owner/checker/profile/collector code, or data.

## Strict preservation and startup worktree policy

Keep all production and other stages frozen: random/dummy normalization, RF/noRF and pion PID, proton cleaning, coordinate conventions, shifts, canonical bins, SIMC cuts/weights/templates/normalization, alignment, pion/kaon stage fits, priors, `w0`, weighted events, fallback/acceptance, yields and errors, E.8/F.6, Method A/Method B, manifest/page contracts, PDFs, and the `no_empirical_residual` profile. The diagnostic may read existing objects and write one **new** named JSON; it cannot add an input gate or alter an exception/fallback path.

Previous Codex report identified pre-existing unrelated tracked `.M` entries in `docs/memory/LEARNINGS.md`, `docs/memory/evidence/e8-4-canonical-five-fresh-f1-lineage-materialization-2026-10-05.md`, and `docs/memory/phases/e8-4-canonical-five-fresh-f1-candidate-lineage-refresh-task-contract.md` with unchanged Git content hashes; and eight untracked unrelated root review files:

- `kaonlt_review.diff`
- `kaonlt_review_baseline_pion_memory_reanchor_2026_10_14.diff`
- `kaonlt_review_e8_current_f1_coherent_materialization_2026_10_11.diff`
- `kaonlt_review_e8_current_snapshot_checkpoint_2026_10_12.diff`
- `kaonlt_review_e8_four_way_BLOCKED.diff`
- `kaonlt_review_e8_four_way_final_source_stop.diff`
- `kaonlt_review_e8_four_way_repaired.diff`
- `kaonlt_review_e8_reweighted_contribution_debug_BLOCKED_2026_10_11.diff`

These are **expected only if actually present and unchanged from startup**; their presence is not permission to overwrite/remove/stage them. Snapshot SHA-256 and Git content state of every pre-existing unrelated path before editing and verify byte-for-byte/identity preservation at end. The newly delivered contract is also expected untracked. An additional unexpected change, changed expected review bundle, or changed task path without authorization stops the task. The *new review bundle* for this task must use a unique filename and must not overwrite any pre-existing review bundle.

## Capture semantics: exact fitted source and template

Build a read-only JSON object from the existing `hist`, frozen parent result, and pion-control cache **after** fitting, only for the exact Q4p4/W2p74/Center/low-epsilon kaon selection and canonical t2 (`t_bin_index==1`); optional t1 (`t_bin_index==0`) is comparative. Inspect source to validate actual `Q2`, `W`, `EPSSET` tokens (`low` versus `lowe`) and the exact moment all inputs exist; do not guess. No event loop beyond iterating already-cached scalar records; no input ROOT tree reopen.

Required JSON schema/fields, with explicit `unavailable_reason` where impossible:

- `schema_version`, `diagnostic_only=true`, `production_objects_mutated=false`, scope/setting/epsilon/t index/edges, observed source HEAD or `not_provided` if runtime cannot attest it, existing parent ID, canonical coordinate fingerprint and coordinate contract, input/source fingerprints, histogram bin axis/edges and fitted analysis window.
- The original **one-pass, first pi-n staged** solve: pass number, fit component, amplitude, uncertainty, `success`, `message`, estimator/mode, anchor and exclusion windows; `fit_bin_indices`, `fit_y`, `fit_sigma`, `fit_template`, `weighted_numerator`, `weighted_denominator` exactly as frozen in `H_pion_fit_step_overlays` / staged history. Include every pass if the configured count exceeds one; separately identify the *first* pass. Recompute the numerator/denominator from retained arrays independently in pure Python and compare within justified floating-point tolerance; record `moment_reconstruction` PASS/FAIL and reason. **Do not recalculate or replace any fit result.**
- Explicit decision classification: `zero_due_to_nonpositive_numerator`, `zero_due_to_invalid_or_unsupported_template`, `positive_staged_amplitude`, or `not_determined`, based only on the actual solve diagnostics and branch outcomes; never infer negative numerator merely from zero amplitude. Include active joint-prior width and stage/refined neutron amplitudes if already present, otherwise mark unavailable.
- For the identical fitted neutron anchor **MM bins**, compute a four-source decomposition from existing `by_t[t_index]["records"]` with source labels `prompt`, `rand`, `dummy`, `dummy_rand`, `nommcuts`, `adj_MM`, and the stored **signed** `coefficient`. Record per-source per-bin count, signed sum, sum of squared signed coefficients, plus combined sum/error; do not convert negative support to nonnegative or use proton factors on pion records. Include source tree/selection/PID/normalization identities already present in source audit. Compare the combined reconstructed histogram against `parent["H_pion_control"]` in the same bins, using the original histogram's ROOT contents/errors **read-only**; record binwise closure, tolerance, and clear failure reasons, not only an integrated total.
- Record aligned neutron template support in the *same exact fit bins* using frozen solver arrays. Record per-t SIMC `pi_n` populated-phi template provenance and fallback if already accessible from existing payload/state. If original pre-normalization per-phi integrals are absent, say `NOT AVAILABLE`; do not reread SIMC or manufacture event-weighted reference shapes in this contract.
- Make all numbers JSON-safe (`allow_nan=False`); represent unavailable/nonfinite values explicitly rather than silently dropping records.

## Output and invocation boundary

Use `alignment_outpath` already supplied to `build_setting_t_bin_pion_parents` (the verified existing analysis-output role), without guessing another filesystem path. Emit a uniquely distinguishable `Center` / Q4p4W2p74 / lowe / canonical-t2 / parent-fingerprint JSON filename, written atomically and never overwriting an existing mismatched file. The output must state that it describes the **current executed parent**, not automatically the August 26 parent.

Prefer the narrowest source-owned hook after the frozen parent/result exists. Do not add CLI flags, environment variables, extra farm driver commands, bundle/checker/collector bypasses, or a new analysis selection in this task. The targeted diagnostic output may be automatic for the exact target setting only; for all other settings and all non-kaon analyses, preserve byte-for-byte normal behavior except this entirely absent artifact. Reuse existing parent/input and coordinate fingerprints; if the targeted state is unavailable, write an explicit diagnostic-unavailable JSON when possible and do not silently claim successful extraction.

**No farm command is authorized by this contract.** Before a later farm handoff, separately audit actual analyzed-source identity, startup owner, output path, packaging, owner/checker/provenance and collection path using the tracked repository procedure. A current-run result cannot retroactively validate the historical August 26 fit; compare against prior reported values only with source/population/cut/bin/coordinate provenance clearly distinguished.

## Deterministic local validation and hard stop

- Pure-Python tests (fake histogram/records, no ROOT) must cover: positive numerator; finite negative numerator clipped to zero; zero support/failed solve; missing diagnostics/unavailable; negative signed bins with prompt/random/dummy/dummy-random; exact bin-edge inclusion; `Sumw2` variance closure; mutated-record/hist detection; multi-pass first-pass preservation; target-setting gating; JSON nonfinite handling; collision-safe atomic output and no production mutation.
- Confirm `python -m py_compile` for each changed Python file, run the new tests and applicable existing `testing/test_pion_t_bin_particle_stage.py` or narrowly relevant deterministic tests; do **not** claim full `main.py`, ROOT/PyROOT or PDF runtime checks.
- Review producer/consumer and preserve all frozen fit, alignment, template, PID, random/dummy and proton inputs. Audit full actual changed paths including the new files. `git diff --check` must pass.
- Update CURRENT, if needed, with one push-stable scientific NEXT conditional on review/push. Avoid phase-status promotion from local tests; no new runtime closure. Preserve the original first-failure question.
- Regenerate/check memory manifest and run ordinary `<PYTHON> -B tools/check_memory_health.py --root .`; report byte counts for CURRENT/MEMORY/CURRENT_HANDOFF and classify warnings. Existing `CURRENT.md` 8KiB soft-limit warning is nonblocking if no material ambiguity; do not undertake unrelated compaction.
- At hard stop, return `BLOCKED` with first violated invariant and unchanged state. No speculative fit correction, forced positive amplitude, new acceptance threshold, changed template mixing, promoted Method A/B, commit/push, or farm execution.

## Handoff

Codex stops uncommitted and unpushed with deterministic test receipts, all intended changed paths, preserved-state hashes, and **one uniquely named complete diff** containing tracked changes and full `git diff --no-index /dev/null` additions. ChatGPT independently reviews that actual diff before user-controlled scoped `git add`, `git commit`, `git push origin test`. After publication, ChatGPT checks pushed source/manifest/CURRENT continuity before any user-controlled farm validation is considered.
