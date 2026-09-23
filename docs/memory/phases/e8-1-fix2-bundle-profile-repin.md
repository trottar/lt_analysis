# E.8.1.Fix.2 — bundle-profile provenance re-pin

## Status and source identity

`CLOSED / RUNTIME VALIDATED` — only for the narrow E.8.1 bundle-profile
provenance re-pin. This profile-only repair started from `test` / `origin/test`
committed HEAD `8985d9a212799c021c4ad1a759a689ea3826e0ea`. ChatGPT inspected the
complete actual diff and found the profile/test provenance repair PASS;
Codex-reported deterministic checks were `NOT RUN by ChatGPT`.

The required analysis/procedure source is the user-pushed, `SOURCE REVIEWED`
E.8.1.Fix.1 reader/test commit:

`8985d9a212799c021c4ad1a759a689ea3826e0ea`.

It replaces the prior E.8.1 procedure-source pin
`0bf1281ebf44868b68ff808090436653eea6ed60` exactly. The user pushed the
reviewed profile re-pin as `9f7094b19f3cdb0887c95f1f972d98471f2ecfdc`
(`Re-pin E8.1 validation profile to repaired source`).

Its fresh `Q4p4W2p74 / Left / lowe` bundle is complete, passes its fail-closed
provenance checks, and admits no unexpected committed files after the required
analysis source. This closes the profile-provenance defect only; it does not
accept E.8.1 because its rendered procedure PDF exposed a separate layout
blocker. See [E.8.1.Fix.2 Left/lowe evidence](../evidence/e8-1-fix2-left-lowe-layout-blocker.md).

## Narrow profile change and preserved boundary

`testing/pion_hgcer_validation_bundle_profile_e8_1.json` changes only
`source_identity.required_analysis_commit`. It retains profile schema v4,
`phase_e8_1_full_background_procedure_pdf_farm_review/v1`, generic-artifacts
mode, the canonical ordered five settings, the one global frozen F.6.2 JSON,
and the PDF/page-manifest pair per setting.

The exact source-range allowlists remain fail-closed: only the profile and its
focused test are allowed committed files after the required analysis source,
and only `docs/memory/` is allowed outside analysis. No `src/` exception is
added. The generic collector, tcsh wrapper, artifact declarations, canonical
settings, accepted F.6.2 artifact/evidence, and all analysis/production physics
remain untouched.

F.6.2 and F.6.2.Fix.5 remain `CLOSED / RUNTIME VALIDATED`; E.8.1 remains
`ACTIVE`; E.8.1.Fix.1 is closed only for its reader defect; and F.6.3 remains
`BLOCKED`. This record does not claim E.8.1 runtime acceptance.

## Deterministic local validation

- `python -B -m py_compile testing/test_pion_hgcer_validation_bundle_profile_e8_1.py` — PASS.
- `python -m json.tool testing/pion_hgcer_validation_bundle_profile_e8_1.json` — PASS.
- `python -B -m unittest testing.test_pion_hgcer_validation_bundle_profile_e8_1 -v` — PASS, 5 tests.
- `python -B -m unittest testing.test_collect_pion_hgcer_validation_bundle -v` — PASS, 28 tests, no skips.

The required memory manifest/integrity/health/bootstrap checks and
`git diff --check` are run after this record and manifest are updated. These
deterministic results were reported by Codex and were `NOT RUN by ChatGPT`
during the complete actual-diff review.

## Provenance path and next action

Next: independently review the active E.8.1.Fix.3 presentation-layout diff;
after its user-controlled push, create and review a separate profile re-pin to
that new source. Then rerun `Q4p4W2p74 / Left / lowe`, inspect the repaired PDF,
and only after that pass return to detailed `Left / highe` and broader coverage.

The intended change set is limited to the profile, its focused test,
`CURRENT.md`, the Fix.1 and Fix.2 phase records, the memory manifest, and the
user-created Fix.2 task contract. Do not create a farm-evidence record until
fresh applicable farm artifacts are supplied and independently reviewed.
