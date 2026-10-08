# E.8.4 canonical-five regenerated-F.1 lineage-preflight repair — task contract

## 1. Authority and scope

This is a narrow **owner/checker source-changing task**, not a modification of KaonLT scientific analysis. Authoritative repository: `https://github.com/trottar/lt_analysis/tree/test`.

Observed published source at contract creation:

```text
branch: test
HEAD / origin/test: 75662969db261f2f6cf6387da09590710d871ae7
```

Codex must independently establish the live branch, exact HEAD, `origin/test`, complete worktree state, and root `AGENTS.md` instructions before proceeding. This SHA is the required start gate unless the user explicitly supplies a reviewed successor contract. Do not reset, clean, stash, stage, overwrite, or disturb unrelated work to force the gate.

Read in the root `AGENTS.md` prescribed order, in full:

1. `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read `docs/memory/MAINTENANCE.md`, `CODEX.md`, `COMMUNICATION.md`, `TOOLS.md`; CURRENT direct references relevant to E.8.4; the preceding F.2/F.3/F.4 equivalence repair and F.2 predecessor repair records; and the exact code listed below. Emit the full in-chat/repository health receipt. Do not reopen closed scientific phases.

## 2. Failed farm gate: directly supplied evidence

At farm source `75662969db261f2f6cf6387da09590710d871ae7`, the attempt
`KaonLT_E8_4_Fix5_CanonicalFive_Q4p4W2p74_20261008-040259_3914085` failed at `f6_3_lineage_preflight` with `lineage_f1_identity_mismatch` inside `testing/run_e8_4_fix5_canonical_five_plot_gate.py:validate_candidate_lineage`.

The supplied owner status records `status=failed`, `analysis_started=false`, `analysis_completed=false`, `artifact_verification_completed=false`, `collection_completed=false`, and `zip_verification_completed=false`. The status verifies the following narrower gates:

- Reviewed materialization verification PASS.
- F.2 candidate replacement PASS: exact predecessor `182433ccd2d13d3f7af80b963e69190ad8f2eb4e295c8216b36585e8ee7409b2` replaced by reviewed `2fa715b2b5c2d5077e38416fb9eceb814e013e72f308c43526fddff36ddd962e`.
- Reviewed F.3 `c5b864...` and F.4 `79e7ceda...` both `no-op`, exact identities preserved.
- Collector source checks PASS; ltsep import isolation PASS.
- Farm checkout `HEAD` and remote-tracking `origin/test` were both `75662969db261f2f6cf6387da09590710d871ae7`; its ordinary worktree contained unrelated modified `src/models/xmodel_kaon_pl.f` and untracked `src/kaon/functions/Q4p4W2p74.model`. Neither may be altered.

The five F.1 files were present but differed from the reviewed raw SHA-256 pins:

| Setting | Reviewed F.1 raw SHA-256 | Observed farm F.1 raw SHA-256 |
| --- | --- | --- |
| Left-lowe | `eb6f659da0511f0f6ea867420fda0ee6feb77e6698509ccd31f80a92cc541c07` | `b5d0b10cb888b87915b3b4090e14d2c03436f93616fd5c42e79c729ecc0748e9` |
| Left-highe | `544ea08f71e74b6b59bc33d05458d4acc051e96fb0327f54ec092a331245c01e` | `3648346330ab01cd52e7d58fb3002d12f758e4f9f14634017ec1e5a02e5b4a9c` |
| Center-lowe | `2b193dec46b10aebc74dfb51634d19a7fef29362fd0a252d940f2894d18a0a16` | `4fd5cc5c691d2f181ffa5134fcf803328e824b7bb94ae04fb8e91a7c8144c724` |
| Center-highe | `c857911396bed03f9e418ca45a609509ebf1eadf8084d4c0ff2f9cc3a9980941` | `fc5f340db724a0bd4dd9b2e527dd67003b088a8751bc53cce71956f28652c251` |
| Right-highe | `e03245494d54a6f812c426d189b220417fdc0888ba1f58a961fbe250ff097652` | `81a09c6dca80cbe5c5dc477ebcece19f297aadced33c7cb30cbf8492e3d7775a` |

These are **raw file hashes**. Scientific equivalence of the observed farm F.1 payloads to the reviewed materialization is **NOT VERIFIED**. Do not claim equivalence, infer a harmless provenance-only change, or refresh pins from this result. The returned failed-attempt status is diagnostic, not accepted runtime closure.

## 3. Source-level diagnosis and narrow objective

In `testing/run_e8_4_fix5_canonical_five_plot_gate.py`, `validate_candidate_lineage` currently checks:

```python
expected_f1 = module.F6_3_CANDIDATE_F1_SOURCE_FILE_SHA256
require(expected_f1 == F1_SHA256 and
        {k: hashes.get(k) for k in expected_f1} == expected_f1,
        "lineage_f1_identity_mismatch")
```

It later also requires `provenance["accepted_f1_source_file_sha256"] == expected_f1` and `scientific_equivalence["mode"] == "exact-lineage"`. This preflight is exact-raw-lineage-only even though the already-reviewed frozen source `src/cuts/pion_hgcer_method_a_parallel_full_procedure.py` implements **both** `exact-lineage` and `equivalent-new-lineage`. For changed raw F.1 lineage, the latter rebuilds F.2, then F.3, then F.4 using existing public scientific builders, requiring exact scientific-projection equality at every stage before factor consumption. That bridge is the scientific decision authority; the owner must not duplicate, bypass or weaken it.

Objective: make the **canonical-five owner preflight** validate regenerated F.1 via the existing exact-scientific-equivalence bridge, while retaining the reviewed candidate pins, strict F.1 artifact validators, 15-parent inventory, factor identity, full provenance, fail-closed scientific comparison, and all later live-cache/parent-preservation gates. The owner must never accept a raw mismatch by itself.

## 4. Allowed and frozen paths

Executable edits are limited to:

```text
testing/run_e8_4_fix5_canonical_five_plot_gate.py
testing/test_run_e8_4_fix5_canonical_five_plot_gate.py
```

The task may update exactly these memory/integrity paths, according to their roles:

```text
docs/memory/CURRENT.md
docs/memory/roadmap/STATUS.md
docs/memory/manifest.json
docs/memory/evidence/e8-4-canonical-five-f1-lineage-preflight-blocker-2026-10-08.md
docs/memory/phases/e8-4-canonical-five-f1-lineage-preflight-repair.md
```

This contract itself is the sole new task contract:

```text
docs/memory/phases/e8-4-canonical-five-f1-lineage-preflight-repair-task-contract.md
```

All `src/` scientific/runtime-analysis files, their tests/strict validators, profile, collector, launcher, shell scripts, ltsep, farm environment, materialized candidates, F.1 outputs and production artifacts are frozen. Never modify scientific projections/exclusion lists, floating tolerances, candidate materialization hashes, reviewed F.1 SHA-256 pins, Method-A factor algorithm, accepted Left/lowe owner or F.2 staging-predecessor policy. If this repair requires anything else, stop `BLOCKED` for fresh review.

## 5. Required owner behavior

1. Keep source definition consistency `module.F6_3_CANDIDATE_F1_SOURCE_FILE_SHA256 == F1_SHA256`, but do **not** require observed F.1 raw hashes to equal frozen historical pins. Require observed raw hashes to be a complete exact five-setting inventory of nonempty valid SHA-256 strings (as provided by the strict loader), with no additional, missing or duplicated settings.
2. Continue requiring F.2/F.3/F.4 files to match their exact reviewed SHA-256 pins, and use the existing strict source-owned wrapper validation and reconstruction. Do not broaden the allowed F.2 predecessor state.
3. For each canonical setting, invoke the unchanged `module.reconstruct_transient_factor_map(...)` using **actual loaded F.1 artifacts and their observed hashes**. Do not change the input scientific objects, select reviewed F.1 rows, or use stale candidate factors on mismatch.
4. Require every returned `scientific_equivalence` to report `all_stages_passed is True`, exact expected schema and exclusion policy, and the correct mode: `exact-lineage` iff *all five* observed raw hashes match the reviewed pin map, otherwise `equivalent-new-lineage`. Require per-stage F.2/F.3/F.4 PASS and no first mismatch, and retain the existing source-authenticated F.3/F.4 identity and parent/inventory checks. Review the existing decision/provenance fields rather than inventing aliases.
5. Require `accepted_f1_source_file_sha256` and `current_runtime_lineage.f1_source_file_sha256` to reflect **observed current F.1 hashes**, while `reviewed_candidate.f1_source_file_sha256` remains the frozen original pin map. Assert these roles remain distinct when the hashes differ. Retain reviewed candidate F.2/F.3/F.4 identities and the current rebuilt lineage in provenance.
6. Validate factor positivity/finite values, exact setting and 15-parent review/inventory, transient factor identity fingerprints, nonproduction/Method-B-exclusion flags, and all existing runtime owner preservation/cleanup and artifact checks.
7. The owner preflight receipt must record both reviewed and observed F.1 hashes and the selected equivalence mode per setting; do not persist transient event-factor arrays or replace historical candidate authorities.
8. A failure at any source/semantic/provenance/inventory gate must occur before analysis and must record the first exact failed invariant. Never silently skip the preflight, retry, restage unrestricted candidates, or infer that changed raw hashes imply equivalence.
9. Preserve acceptance of the exact old raw F.1 case and maintain current owner operation ordering: materialization verification -> collector source preflight -> candidate staging -> path/ltsep isolation -> lineage preflight -> isolated analysis only if PASS -> artifact/page checks -> collector/ZIP/status/preservation.

If available fields make any assertion above impossible without modifying frozen code, stop and document the concrete blocker rather than weakening the check.

## 6. Deterministic local validation

Add tests in `testing/test_run_e8_4_fix5_canonical_five_plot_gate.py` with real public scientific calculators and strict validators (no mocked science acceptance):

- Exact original F.1 hashes retain `exact-lineage` and previous result.
- A regenerated F.1 fixture with **different raw hashes but exact F.2/F.3/F.4 scientific projections** succeeds in `equivalent-new-lineage` for all five settings; provenance separately retains original reviewed and actual observed F.1 hashes. The synthetic fixture must explicitly distinguish benign serialization/provenance variation from scientific change.
- A scientifically changed but internally valid F.1 fixture fails the correct first F.2/F.3/F.4 scientific-equivalence or F.1 invariant *before analysis*, not merely at missing bytes.
- Missing/extra setting, malformed F.1, raw identity mismatch in F.2/F.3/F.4, duplicate-key/nonfinite JSON, invalid F.4 authority, altered factor identities, failing parent inventory and intentionally manipulated/stale provenance remain fail-closed.
- Both preflight-only and full owner gate paths must use the same reconciled lineage preflight and not bypass it.
- The accepted F.2 predecessor replacement test and unchanged F.3/F.4 policy must continue passing.
- Run the existing owner, parallel-full-procedure, F.2/F.3/F.4 and collector regression suites with available local dependencies. Do not claim ROOT or farm success from deterministic tests.

Do not mock or patch the scientific comparison/builders to force acceptance. It is acceptable to patch immutable fixture pin constants for synthetic science fixtures as existing tests do, while retaining real builder and comparison execution. Show assertions that distinguish provenance-only and physics changes. A failure should name the first failed scientific stage and mismatch path where supported.

## 7. Integrity, memory and review bundle

Preserve unrelated local state and the index/refs. Report exact changed paths and verify that `src/` and unrelated bytes are untouched. Require source check, deterministic tests, actual diff and `git diff --check`.

Discover `<PYTHON>` under `docs/memory/TOOLS.md` and run:

```text
<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/check_memory_health.py --root .
<PYTHON> -B tools/memory_bootstrap.py --root . --json
```

Report full memory health, CURRENT/MEMORY/CURRENT_HANDOFF byte counts, warnings and classification, manifest/bootstrap, deterministic test counts/skips and preservation checks. CURRENT’s existing 8-KiB soft warning alone is nonblocking if it creates no material ambiguity; batch consolidation at the next meaningful checkpoint.

Produce one repository-root temporary `kaonlt_review.diff` comprising the complete tracked diffs and complete `git diff --no-index /dev/null ...` additions for every new file, including this task contract. Do not stage merely for review. Stop for **ChatGPT actual-diff acceptance** before user commit/push. Codex must not commit, push, run the Jefferson Lab farm, or touch its files. The user alone publishes accepted changes.

## 8. State and next gate

The newly failed 2026-10-08 attempt is diagnostic only. Preserve prior Left/lowe and isolated lineage-preflight closures with their exact scopes. The new F.1 lineage-preflight failure leaves canonical-five full runtime, final E.8 and F.6.4 `BLOCKED`; the F.2 staging and earlier scientific-equivalence repairs remain `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`. Method A remains detached/non-production, Method B diagnostic-only and numerically excluded. Absolute-SIMC normalization remains separately `BLOCKED`.

When local implementation is complete, record `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`, conditional on actual-diff acceptance. CURRENT must have one push-stable substantive NEXT: **after ChatGPT actual-diff review, user commit/push and pushed-state synchronization, perform a fresh farm-readiness audit of the same single isolated Q4p4W2p74 canonical-five owner gate**. Only `Farm readiness: PASS` authorizes another user-controlled run. Do not issue a run command from this task.

## 9. Hard stops

Return `BLOCKED` without workaround if: live source identity is wrong; unrelated local work cannot be preserved; a scientific/runtime analysis source needs modification; existing bridge cannot accept an independently validated exact scientific projection; invalid current F.1 can pass; stale candidate factors can be used; either reviewed or current lineage is hidden/lost; finite/positive, provenance, parent/inventory, live-cache or strict hash gates are relaxed; a test, manifest, health hard gate or relevant source diff check fails; or farm execution would be required for local proof.
