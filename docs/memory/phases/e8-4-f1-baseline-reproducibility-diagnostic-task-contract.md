# E.8.4 — Detached F.1 baseline-input reproducibility diagnostic: task contract

## 1. Authority, starting gate, and purpose

Authoritative repository: `https://github.com/trottar/lt_analysis/tree/test`.

Observed published `test` HEAD during ChatGPT's source audit on 2026-10-08:

`ebf3a0be2f50c788c22defe3601fd5b7fee48644`

**Codex must establish its actual local branch, HEAD, `origin/test`, full porcelain status, index, and relevant worktree state independently.** This identity is a *start gate*, not a license to overwrite unrelated modifications or a claim that the user's local worktree is clean. If the source identity has moved, stop for a fresh source review; do not reset, clean, stash, force-checkout, or silently adapt this contract.

Read in the repository-root `AGENTS.md` mandated full startup order:

1. `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read `docs/memory/MAINTENANCE.md`, `CODEX.md`, `TOOLS.md`, `COMMUNICATION.md`, CURRENT's E.8.4 direct references, the 2026-10-05 fresh-F.1 lineage/materialization and canonical-five preflight evidence, the October 8 F.1 owner repair contract/record, and only further records needed for this audit. Emit the **full KaonLT health check** at startup and before the implementation gate.

**Task class:** narrow detached diagnostic and test implementation, not a new scientific correction, provenance identity repair, runtime-owner repair, production algorithm change, candidate repin, or a broad E.8 reopening. The live October 8 failed farm owner gate and detached F.4 comparison supersede CURRENT's stale instruction to retry the owner. Correct CURRENT's sole NEXT at the required checkpoint without reopening closed phases.

## 2. Supplied evidence and the exact unresolved question

The isolated canonical-five owner at farm source `ebf3a0be2f50c788c22defe3601fd5b7fee48644` failed **before starting analysis** at:

`f6_3_scientific_equivalence_mismatch:f4:$.parents[0].absolute_baseline_parent_sum`

The diagnostic comparison input was `KaonLT_F4_FailedPreflight_Comparison_20261008-044305.json` (user-supplied farm output; not a successful owner artifact). It establishes:

- F.2 scientific projection: exact PASS.
- F.3 scientific projection: exact PASS.
- F.4 scientific projection: exact FAIL; deviations confined to the three Left/low-ε canonical-t parents in the reported core parent comparisons.
- F.4 `Left-lowe, t=0` reviewed `absolute_baseline_parent_sum = 0.12072578658633845`; observed/current reconstruction `0.12072578661296547` (difference about `2.6627e-11`). Small magnitude is **not** scientific-equivalence acceptance.
- Current `Left-lowe` F.1 raw SHA-256: `b5d0b10cb888b87915b3b4090e14d2c03436f93616fd5c42e79c729ecc0748e9`.
- Frozen reviewed `Left-lowe` F.1 raw SHA-256: `eb6f659da0511f0f6ea867420fda0ee6feb77e6698509ccd31f80a92cc541c07`.
- Current-input five-setting raw hashes are in the supplied comparison; the other settings' F.4 core parent comparisons were reported equal.

The user searched loose Left/low-ε F.1 files in the configured farm output and transferable-bundle directories and checked matching ZIP members. **No found file/member had the reviewed `eb6f659d...` identity.** Observed archival hashes included `f16c5e89...`, `dc3af0ba...`, and `bece5900...`; none may be renamed or promoted to the reviewed input. The raw reviewed F.1 payload and exact event-level cause remain **NOT VERIFIED**. Other filesystem locations were not comprehensively searched.

The existing detached comparator `testing/compare_method_a_current_baseline_authority.py` proves stage-level comparison but does not give an event-level coefficient-versus-weight forensic decomposition. This task must **not** duplicate or replace that canonical comparator's scientific gate.

## 3. SOURCE VERIFIED: protected chain and existing evidence

Read the actual current source before coding, especially:

- `src/cuts/pion_hgcer_event_contract.py`: `_weight_contract`, `_build_pion_side`, F.1 precursor event construction, source-algebra/parent weight-provenance fingerprints;
- `src/cuts/pion_hgcer_method_a_acceptance_contract.py`: `_application_population`, F.1 payload/writer/schema and event assignment;
- `src/cuts/pion_component_subtraction.py`: `build_simc_shape_pion_control_weights`, `simc_shape_pion_weight_from_value`, source coefficient specs;
- `src/cuts/pion_component_fits.py`: accepted parent fit/amplitude source (including the joint fit path);
- `src/normalize/get_eff_charge.py`: `normfac_data` and `normfac_dummy` authority;
- `src/cuts/pion_hgcer_method_a_parent_preserving_correction.py`: F.4 inclusion, identity validation, `math.fsum`, parent/source/child baseline sums;
- `src/cuts/pion_hgcer_method_a_parallel_full_procedure.py`: exact scientific projection bridge;
- existing F.1 strict validators, comparator, and relevant tests.

For each eligible F.1 application event, the frozen Phase-A owner formed the signed baseline contribution as `b_i = signed_source_coefficient_i * baseline_pion_weight_w0_i`. The original `w0` was selected by analysis-missing-mass bin from the **frozen final diagnostic application** weight array. Coefficients originate from the established prompt/random/dummy source algebra and data/dummy effective-charge normalization. F.4 consumes *F.1's persisted* signed contributions and applies its existing inside-phi/parent-t ownership. Do not change any of those calculations.

An F.1 JSON carries event `source_label`, `entry_index`, `t_index`, `phi_status`, `phi_index`, `analysis_MM`, `signed_source_coefficient`, `baseline_pion_weight_w0`, and `signed_baseline_event_contribution` as applicable, plus F.1/Phase-A identity fingerprints. **It does not itself expose full upstream ROOT parent component-template objects, fitted amplitude arrays, the full final weight-array/bin-edge object, or the reviewed missing raw F.1 file.** A source-only comparison cannot conclusively isolate fit/template numerical variation from a missing historical input. The diagnostic must expose these limits explicitly; never fabricate fit-amplitude/template fingerprints or claim original parent-MM bin assignments where histogram edges are unavailable.

## 4. Allowed files and freeze boundary

The implementation may add exactly these **new** executable/test files:

- `testing/diagnose_e8_4_f1_baseline_reproducibility.py`
- `testing/test_diagnose_e8_4_f1_baseline_reproducibility.py`

The single standalone task contract is this file:

- `docs/memory/phases/e8-4-f1-baseline-reproducibility-diagnostic-task-contract.md`

For durable recording, only the following additionally may change or be added:

- `docs/memory/CURRENT.md` (one accurate objective/BLOCKED/NEXT checkpoint; consolidate rather than expand if needed)
- `docs/memory/manifest.json` (generated by the existing owner)
- `docs/memory/phases/e8-4-f1-baseline-reproducibility-diagnostic.md` (brief implementation/validation record, if needed)

**All preexisting `src/`, other `testing/` files, scientific code, F.1 producer/serializer, F.2/F.3/F.4/F.5 calculators, full-analysis owner, profiles, collection/checker, launcher, candidate materialization/authority hashes, frozen scientific policies, and production outputs are forbidden to modify.** No event-content restamping, source artifact rewriting, ROOT imports, numpy/scipy fit reruns, repairs to the scientific-equivalence bridge, rounding, tolerance relaxation, or new farm invocation belongs here. If real upstream fit/template capture is necessary for the diagnostic, stop `BLOCKED` and propose a *separate* reviewed observational instrumentation contract; do not widen this one.

## 5. Implement one deterministic, read-only CLI

Provide a pure-standard-library Python 3 CLI with explicit arguments (no ambient path inference):

- `--current-f1 PATH` required: existing Left/low-ε F.1 JSON, never modified.
- `--expected-current-sha256 SHA256` required: exact expected raw bytes.
- `--reference-f1 PATH` optional, **only** with `--expected-reference-sha256 SHA256`; requires distinct bytes/path and independent strict validation. If absent, do **not** label event-level drift resolved.
- `--reviewed-f4 PATH` optional with required `--expected-reviewed-f4-sha256 SHA256`: independently pinned reviewed F.4 wrapper, used **only** for parent baseline scalar comparison, not factor reconstruction or acceptance.
- `--output PATH` required: fresh output JSON path; fail if already exists, is a symlink, or aliases an input (including resolved paths/inodes); never overwrite inputs or accepted outputs. Require the output's existing parent directory, write atomically only after all preflight validations; no partial artifacts on failure.

Do not encode a farm path or a particular local checkout into defaults. Restrict the scientific claim to `Q4p4W2p74` / `Left-lowe`, validate the complete setting identity and known F.1 schema/contract values, and refuse other settings. For raw input reads use duplicate-key and nonfinite-JSON rejection, strict SHA-256 matching, object/schema validation, finite numeric checks, unique `(source_label, entry_index)` application IDs, complete canonical t inventory, and exact phi geometry/identity semantics required by the source. Validate the frozen `b_i = c_i * w0_i` identity under **existing** source semantics; do not invent a competing tolerance or weaken a strict gate. Prefer reusing a pure-Python existing validation helper where possible; if unavailable, be explicit that the CLI is a forensic parser, not a substitute for the canonical F.1/F.4 validators.

Output one deterministic machine-readable JSON object with `schema_version`, `non_authoritative=true`, `failed_owner_diagnostic_only=true`, `production_objects_mutated=false`, input SHA-256s, setting, source-verified lineage role, parser/identity checks, data availability and scientific limits. Include at least:

1. **Self-inventory** for the current file: all application counts and `inside_phi`-eligible counts separately; per canonical t and per source label, counts, `math.fsum(b)`, `math.fsum(abs(b))`, zero-weight count, coefficient distribution/distinct values, baseline-weight distribution (min/max/zero/distinct counts), finite MM range, and canonical-phi grouping. Clearly label group totals that F.4 actually consumes; do not assume outside-phi rows are F.4 applications.
2. **Recomputed baseline parent sums** using F.4's eligible rows and `math.fsum`. If reviewed F.4 is provided, verify its pinned bytes/schema/setting and compare `application_event_count`, `baseline_parent_sum`, and `absolute_baseline_parent_sum` for each of the three parents with exact equality plus reported absolute/relative diagnostics. These comparisons never grant runtime acceptance.
3. **Reference F.1 mode, when available:** validate exact two-input provenance and pair rows by `(source_label, entry_index)` with no positional joins. Summarize missing/added identities, changed `t_index`, phi assignment, analysis MM, `signed_source_coefficient`, `baseline_pion_weight_w0`, and `signed_baseline_event_contribution`; report first changed IDs by stable sort and counts/absolute-max/signed/absolute deltas for each variable and parent/source. Separate exact equality from magnitude. No rounding/rewrite or tolerance-based PASS. Preserve exact reference identity independently of current lineage. An alternative archive hash may be used only as explicitly labeled *historical comparison*, not as the reviewed `eb6f...` authority.
4. **Evidence-completeness matrix:** observable from F.1 now; observable only if a properly identified reference F.1 exists; unavailable without full upstream parent fit/templates/weight-array/bin edges. Tag classifications `SOURCE VERIFIED`, `RUNTIME VERIFIED` (only where input supplied), `INFERENCE`, or `NOT VERIFIED` accurately. Never declare fit numerical drift, source coefficient drift, or event-MM reassignment as fact if evidence is insufficient.
5. **Bounded, non-sensitive output:** no raw full event arrays or factor arrays, no HTML/PDF generation. For detailed sample differences retain at most a fixed small number (e.g. 10) of sorted `(source_label, entry_index)` identifiers per category. No user-specific filesystem path embedded except explicit input/output provenance where necessary; no machine-wide scanning or implicit archive retrieval.

This tool must not import or call full KaonLT runtime, ROOT/PyROOT, SciPy optimization, the analysis launcher, or the canonical-five owner. It must never mutate an F.1 artifact, reviewed F.4, analysis directory, repository refs, or environment variables. Use `main(argv=None)` plus a callable pure comparison builder to permit deterministic tests.

## 6. Required deterministic tests and audit evidence

Use synthetic but **contract-shaped** Left/low-ε F.1 fixtures and reviewed F.4 parent aggregates. Do not patch a scientific acceptance function to force PASS. Test at minimum:

- Valid self-inventory and three parent signed/absolute sums match explicit `math.fsum` calculations (including cancellation and zero values); inside/outside phi handling matches existing source ownership.
- Same events reordered yield identical scientifically keyed summaries, while each input raw hash and role remain correctly distinct.
- `--reference-f1` shows a change in coefficient alone, weight alone, missing-mass/phi/t assignment, one missing and one added identity; distinguish categories and first offending stable IDs. In particular, a common coefficient change must not be falsely classified as a changed weight.
- The observed ~`2.66e-11` reviewed/current t0 F.4 mismatch can be represented and reported **without** passing an exact scientific-equality gate.
- Without reviewed F.1, output explicitly says event-level source attribution **NOT VERIFIED**, even if reviewed F.4 is supplied.
- Wrong hash, wrong setting, malformed/duplicate/nonfinite JSON, missing/duplicated IDs, invalid canonical-t/phi geometry, missing fields, invalid `b_i`, wrong reviewed F.4 schema, same resolved input/output, preexisting output and symlinks fail before any input mutation.
- Input file byte hashes, indexed event identity, no output on negative path, and deterministic byte-equivalent JSON on repeated clean runs with different fresh output paths; no hidden timestamps/cwd-dependent fields in science summary.

Run scoped tests and an applicable isolated regression subset with the actual local Python interpreter; record pass/fail/skips exactly. Use `python -m py_compile` or equivalent only for source syntax; do not report ROOT or farm validation. Do not synthesize an accepted F.1 or pass the exact owner gate locally.

## 7. Integrity, checkpoint, and ChatGPT review handoff

- Record initial and final local status, HEAD, `origin/test`, index/ref state, and all unrelated preexisting changes; preserve them byte-for-byte. Audit `git diff --check` and ensure executable edits are limited to the two new `testing/` files. No `src/` changes.
- Update `CURRENT.md` narrowly to include the *newer* October 8 failed F.4 pre-analysis gate and its diagnostic-only comparison, retain closed E.8.2/E.8.3 and narrow accepted preflight scopes, and replace the stale canonical-five retry NEXT with **the detached F.1 reproducibility diagnostic evidence gate**, conditional on actual-diff review, user commit/push and pushed-state synchronization. Keep canonical-five runtime, final E.8 and F.6.4 `BLOCKED`; preserve no-empirical-residual, Method-A detached and Method-B excluded boundaries. Reconcile CURRENT's existing >8 KiB soft warning by focused consolidation if safe; hard integrity failures stop.
- Write the brief diagnostic implementation record only as evidence of source/local checks, never as farm acceptance.
- Run `tools/update_memory_manifest.py --root . --write`, `--check`, `tools/check_memory_health.py --root .`, `tools/memory_bootstrap.py --root . --json` under the interpreter specified by tracked TOOLS and detected environment. Report health, soft/hard violations, manifest/checker outcomes, and full KaonLT health receipt.
- Create one temporary **complete** `kaonlt_review.diff` in the repository root including tracked diffs plus complete new-file additions using `git diff --no-index /dev/null ...` for each intended untracked file (including this contract), without staging merely for review. Preserve any existing file of that name rather than overwriting it.
- Stop at **ChatGPT actual-diff review**. Do not stage, commit, push, alter remote refs, run a farm operation, package failed evidence as accepted, or announce runtime success. The user alone performs subsequent accepted Git publication and Jefferson Lab farm commands.

## 8. Success, hard stops, and next gate

**Local implementation success** is only: one read-only deterministic CLI and its tests, consistent checkpoint/manifest, untouched scientific/runtime source, complete reviewed diff available, and explicit completeness limits. Work-state is `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING` only if a genuine future farm/source-dependent gate is required; do not call a local parser test runtime-validated science. No scientific-promotion decision follows.

**Stop `BLOCKED`** if source/HEAD/worktree conflicts with the start gate; existing frozen validators cannot be respected; a new source/producer instrumentation path is necessary; accepted authorities would be modified; a test, manifest or preservation gate fails; original F.1 availability is invented; the missing template/fit provenance is falsely treated as reconstructed; or the exact comparison semantics would change.

**NEXT after accepted code review and the user's publication:** pushed-state synchronization, then one narrow review of a diagnostic JSON built from the SHA-pinned current Left/low-ε F.1 and SHA-pinned reviewed F.4 (optional distinct historical F.1 only for separately labeled context), with source/owner semantics preserved. A new diagnostic invocation/farm-readiness gate must be independently reviewed before prescribing an ifarm command. The full canonical-five owner retry remains unauthorized until its scientific-equivalence blocker is resolved through separately reviewed evidence and any necessary follow-on contract.
