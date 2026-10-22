# KaonLT — Candidate C event-capture prepass boundary repair

**Task class:** Narrow farm-evidenced failure repair, not a new phase or physics change.
**Observed public `test` HEAD:** `04f1ec403a9555634e61c229aa347cb76d0c78d3` (audit anchor; recheck current identity).
**Status:** `BLOCKED` at the first Center/Left low-ε event-capture farm child run.
**Scientific NEXT:** The original six-parent signed-event capture and detached A/B/C replay, with no chosen regularization strength.

## 1. Governing authority and initial gate

Establish the actual live `test`/`origin/test`/local HEAD and full local staged/unstaged/untracked state. Read the tracked root `AGENTS.md` and, in exact order, `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, `docs/memory/USER.md`; then `docs/memory/MAINTENANCE.md`, `docs/memory/CODEX.md`, the October 21 event-validation contract, the active baseline-pion investigation, and only task-relevant source. Emit the full repository-defined health check. Preserve unrelated local state; stop on conflicting edits, unexpected paths, or changed authority. Codex must not commit, push or run the Jefferson Lab farm.

## 2. Direct evidence and first invariant

The user-supplied farm owner receipt for source `04f1ec403…` reports `analysis_returncode=1`, `stage=completion_markers`, `failure_reason=analysis_failed`, collection/packet checks not reached, and passed ordinary-worktree and installed-ltsep preservation. The user's corresponding failed-child log ends in:

```
src/utility/bg_optimization.py::_evaluate_phi_candidate
  -> find_yield_data([candidate_hist], candidate_inp)
  -> src/binning/calculate_yield.py::calculate_yield_data
  -> testing/pion_local_t_event_validation_capture.py::observe
  -> verify_accepted_input_identity
ValueError: accepted_AB_input_identity_mismatch:Center:data
```

The preceding log explicitly declares `context=bg_optimizer_pre_parent`, `subtraction_mode=single_scale`, `parent_required=false`, `free_component_fit=false`. Source confirms `build_bg_optimization_prepass_config` sets `bg_opt_prepass_mode=True`, `analysis_subtraction_context=BG_OPT_PREPASS_CONTEXT`, and scalar `particle_subtraction_mode`; `_evaluate_phi_candidate` legitimately runs preliminary yield calculations before the authoritative pion parents. The new observer hook at `calculate_yield.py` currently checks only `LT_ANALYSIS_PION_LOCAL_T_EVENT_VALIDATION=1`, and `observe` filters merely setting/Q²/W/particle/ε. That is an incorrect *capture-stage* integration boundary.

**Do not infer that the observed data SHA-256 mismatch is solely a prepass artifact.** The mismatch might recur at the canonical parent stage. Preserve all accepted A/B input SHA-256+byte pins and fail-closed identity semantics. No physical C event evidence exists from this failed run; do not package, replay or interpret partial outputs or rerun automatically.

## 3. Narrow required implementation

1. In `src/binning/calculate_yield.py`, keep the opt-in default-off, but invoke the observer only for the *authoritative* kaon `simc_shape_components`/`t_bin` application, not `bg_optimizer_pre_parent`, `bg_opt_prepass_mode=True`, scalar/prepass candidate evaluations, or preparent caches. Reuse established authoritative mode/scope/context resolvers and existing validated frozen-parent collection. Place the stage eligibility test **before diagnostic import, file hashing or output writes**, after the normal cache assignments. Do not change the production calculation, results, cache population, optimizer, parent fitting, or yield handling. The supported authoritative stage with missing/incoherent parent/cache must fail closed; do not silently skip a defective eligible stage.
2. Where useful, harden `testing/pion_local_t_event_validation_capture.py::observe` with the same authoritative-stage condition as a defensive standalone boundary. Keep its existing file-identity checks **unchanged in authority**. Do not downgrade or catch `accepted_AB_input_identity_mismatch:*` and do not accept a different source merely because its normalized yields or source labels look similar.
3. Add deterministic positive/negative tests to `testing/test_pion_local_t_event_validation.py`: prepass (`bg_opt_prepass_mode=True`, `analysis_subtraction_context=bg_optimizer_pre_parent`, `single_scale`) must neither invoke the observer nor import its diagnostic module, open/hash input ROOT files, or create a packet; repeated optimizer candidates remain untouched. The correct authoritative Center and Left modes invoke once per setting with real-like parents/cache and leave all original arrays/histograms/ROOT objects unchanged. A **canonical** wrong data/dummy hash must still raise its precise invariant. Test inconsistent/missing authoritative parents fail rather than silently disappear, and test disabled opt-in equivalence. Use fake/mocked ROOT interfaces as needed, not farm/ROOT acceptance claims.
4. Independently examine whether `hist["InFile_DATA"]`/`hist["InFile_DUMMY"]` are `TFile` handles (they are assigned from `open_root_file()` in `rand_sub.py`) or trees before altering any source-file identification code. Preserve verified backing-file resolution; do not replace actual file identity with tree names or generated names.
5. If an exact read-only check of *current* Center/Left data/dummy file SHA-256 and byte count against the already accepted A/B pins is possible without running the full analysis, document its prerequisites and exact source-owned path derivation. Do not invent a farm path, alter an accepted pin, treat raw file differences as physically equivalent, or block the stage-boundary repair on unavailable farm files. On a persistent mismatch, record the source identities and **stop** for independent scientific/provenance adjudication; no automatic override/replay.

## 4. Scoped files and frozen boundaries

Allowlist only:

- `src/binning/calculate_yield.py` (default-off observer dispatch only)
- `testing/pion_local_t_event_validation_capture.py` (optional defensive stage guard/error-only diagnostics)
- `testing/test_pion_local_t_event_validation.py` (regression tests)
- `testing/test_run_pion_local_t_event_validation_gate.py` (only if required for precise fail-closed owner regression; no owner source or profile change)
- `docs/memory/phases/kaonlt-pion-event-capture-prepass-boundary-fix-2026-10-21-task-contract.md` (this single task contract)
- `docs/memory/CURRENT.md` (short accurate failed-runtime BLOCKED state and *same* scientific NEXT; do not rewrite unrelated closures)
- `docs/memory/investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md` (bounded October 21 observed failure and cause)
- `docs/memory/manifest.json` (regenerated integrity file)

Keep `src/utility/bg_optimization.py`, `src/cuts/rand_sub.py`, all physics and simulation source, all fitter/transfer code, `run_Prod_Analysis.sh`, `src/main.py`, the published owner, collector, profile, A/B and C comparators and existing original data/MC unchanged. Preserve random/dummy/proton/pion ordering, selection, timing, templates, priors, signed event weights, clipping/floor semantics, acceptance, parent/phi-child ownership, yields, L/T and cross sections. Method A remains detached, HGCer Method B numerically excluded; no C lambda promotion.

## 5. Local checks and handoff

Run available deterministic tests for the affected observer and owner, appropriate background optimizer/yield prepass regressions, Python compile and `git diff --check` (including all new-file additions). No regression may require ROOT/farm to be called PASS locally. Regenerate and verify the memory manifest, run `tools/check_memory_health.py --root .` using a discovered working interpreter, classify CURRENT soft-size warning as blocking/nonblocking under `MAINTENANCE.md`. Preserve HEAD, refs, staged index and unrelated files.

Provide one complete actual-diff review bundle, including the full new-contract addition; report exact changed paths, source/manifest/health/test results, any actual blocking failures, and the narrow next gate. ChatGPT independently reviews the diff; **only the user** commits/pushes after acceptance. ChatGPT then verifies published HEAD and changed blobs. **Only after pushed-state and farm-readiness checks** may one user-controlled, new-attempt-owned Center/Left low-ε capture gate run; the failed attempt and its receipts remain immutable. No farm command is authorized by this contract alone.
