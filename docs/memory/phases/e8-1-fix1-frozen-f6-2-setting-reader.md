# E.8.1.Fix.1 — frozen F.6.2 persisted-setting reader

## Status and farm finding

`CLOSED / RUNTIME VALIDATED` — only for the narrow frozen persisted
parent-setting reader defect. This repair started from `test` / `origin/test`
committed HEAD `1c1a0be5da997043378fb49263b1de6ffc0d473c`; ChatGPT independently
inspected the actual complete diff and found the reader/test repair PASS.
Codex-reported deterministic tests were `NOT RUN by ChatGPT`.

The user pushed the reviewed repair as
`8985d9a212799c021c4ad1a759a689ea3826e0ea` (`Fix E8 frozen F6.2 setting
reader`). The separate Fix.2 profile re-pin subsequently targeted this exact
analysis/procedure source.

The fresh, complete Fix.2 `Left / lowe` bundle now establishes this narrow
runtime result: the E.8 section was available with its ordinary 47-page
inventory, rather than ending at `full_background.e8.unavailable` with
`parent_setting_identity_invalid`. Its direct provenance and PDF-layout blocker
record is [E.8.1.Fix.2 Left/lowe evidence](../evidence/e8-1-fix2-left-lowe-layout-blocker.md).
That success neither accepts E.8.1 nor claims closure for the later layout repair.

Fresh targeted Jefferson Lab evidence exposed the blocker in
`KaonLT_E8_1_Debug1_Fix1_Q4p4W2p74_Left_lowe_20260923-012755.zip`:

- bundle SHA-256:
  `490036e8d6c53d6e025da99a5c6c4ddcd529b6573007534efb7aaceac1acd002`;
- bundle manifest: `complete=false`;
- frozen F.6.2 JSON SHA-256 remained the accepted
  `5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1`;
- Left/lowe procedure-PDF SHA-256:
  `d1bb7f945aee7b04b4c9baa0e25ba9b9e0d38956ad42086737b938a6b97ade01`;
- Left/lowe page-manifest SHA-256:
  `6bebd1b8b0382b946c7e9e4ef5e0eca5613dc2871d5e31f5c9020ddc09574828`;
- the manifest ended at `full_background.e8.unavailable`, with literal reason
  `frozen_f6_2_authority_rejected:parent_setting_identity_invalid`.

The output is runtime evidence that the targeted execution reached downstream
ordinary procedure rendering. Its incomplete supplied provenance does not
establish the full launcher path or prove the intentional full-high-epsilon
skip. It does not reopen F.6.2, establish E.8.1 acceptance, or create an
accepted runtime-evidence record.

## Frozen scientific boundary

F.6.2 and F.6.2.Fix.5 remain `CLOSED / RUNTIME VALIDATED`. The accepted JSON,
SHA-256, artifact fingerprint, validation fingerprint, numerical content,
populations, geometry, bootstrap results, and Method-A interpretation remain
immutable. Method A remains detached/non-production and Method B remains
diagnostic only.

This repair changes no production weights, yields, cuts, canonical binning,
normalization, templates, priors, subtraction, SIMC, efficiency, cross section,
L/T separation, uncertainty, E.8.1 profile, collector, bundle wrapper, or
accepted artifact. The ordinary procedure PDF still retains D.6-D.9 and appends
E.8; it does not restore D.10/D.11 or E.2-E.7.2 pages.

## Reader repair and regression coverage

The E.8 reader remains fail-closed and still requires the exact accepted input
SHA-256, artifact schema/fingerprint, validation schema/fingerprint, and full
15-parent / 135-child canonical inventory. The defect was only its synthetic
parent-setting check: it treated semantic `epsilon_setting` (`low`/`high`) as
the canonical filename token (`lowe`/`highe`) and required a non-persisted
`ordinal`.

The repaired reader validates each persisted parent setting directly against
its `setting_id`: exact `Q2=4p4`, `W=2p74`, `kinematic_token=Q4p4W2p74`,
`particle_type=kaon`, matching `phi_setting`, matching
`epsilon_filename_token`, and the exact `lowe -> low` or `highe -> high`
semantic relation. It neither adds nor infers `ordinal`, and does not alter the
loaded artifact.

Focused regression coverage now uses the actual persisted setting schema,
proves `Left-lowe` selects exactly its three canonical parents, and fails closed
for mismatched phi, filename token, semantic epsilon, Q2, W, kinematic token,
particle type, missing setting, and non-mapping setting. The stale ROOT-aware
test now expects canonical D.6-D.9 per-t pages followed by the explicit E.8
unavailable setting page when no E.8 payload is supplied, while preserving
source ROOT objects. Production rendering was not changed for that assertion.

## Deterministic local validation

- `python -B -m py_compile src/cuts/full_background_subtraction_plots.py testing/test_full_background_subtraction_plots.py` — PASS.
- `python -B -m unittest testing.test_full_background_subtraction_plots.FullBackgroundSubtractionE8Tests -v` — PASS, 10 tests.
- `python -B -m unittest testing.test_full_background_subtraction_plots -v` — PASS, 100 tests with 17 expected skips. The stale real-ROOT assertion is skipped locally because PyROOT is unavailable; no ROOT/farm result is claimed.

These deterministic results were reported by Codex and were `NOT RUN by
ChatGPT` during the actual complete-diff review. The required memory
manifest/health/bootstrap checks and `git diff --check` are run after this
record and manifest are updated.

## Scope, chronology, and next action

The intended change set is limited to:

- `src/cuts/full_background_subtraction_plots.py`;
- `testing/test_full_background_subtraction_plots.py`;
- `docs/memory/CURRENT.md`;
- this phase record;
- `docs/memory/manifest.json`; and
- the user-created task contract
  `docs/memory/phases/e8-1-fix1-frozen-f6-2-setting-reader-task-contract.md`.

Next: independently review the active E.8.1.Fix.3 presentation-layout diff;
after its user-controlled push and a separate profile re-pin to its new source,
rerun `Q4p4W2p74 / Left / lowe` and inspect the repaired PDF before returning to
detailed `Left / highe` review or broader coverage.
