# E.8.1.Fix.2 — bundle-profile provenance re-pin

## Status and source identity

`SOURCE REVIEWED` — this profile-only provenance repair starts from `test` /
`origin/test` committed HEAD `8985d9a212799c021c4ad1a759a689ea3826e0ea`.
ChatGPT inspected the complete actual diff and found the profile/test provenance
repair PASS. Codex-reported deterministic checks were `NOT RUN by ChatGPT`.
This is not a ROOT/PyROOT, farm, or runtime acceptance claim.

The required analysis/procedure source is the user-pushed, `SOURCE REVIEWED`
E.8.1.Fix.1 reader/test commit:

`8985d9a212799c021c4ad1a759a689ea3826e0ea`.

It replaces the prior E.8.1 procedure-source pin
`0bf1281ebf44868b68ff808090436653eea6ed60` exactly. The future profile-repin
commit is not known and is intentionally not written into the profile; after
user commit/push and ChatGPT pushed-state review, it will become the wrapper
`--bundle-commit` for the next bundle-only farm collection.

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
`ACTIVE`; E.8.1.Fix.1 remains `SOURCE REVIEWED`; and F.6.3 remains `BLOCKED`.
No ROOT/PyROOT, farm, or runtime validation is claimed here.

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

Next: user-controlled commit/push of this `SOURCE REVIEWED` re-pin, then
ChatGPT pushed-state review. The resulting profile-repin commit supplies wrapper
`--bundle-commit`; the wrapper creates a detached worktree at that commit,
invokes this E.8.1 profile with the generic collector, and the collector applies
its exact source-ancestry/range checks before archiving only existing declared
artifacts. The next farm gate is a fresh `Q4p4W2p74 / Left / lowe` bundle; only
after it passes does review return to detailed `Left / highe` and broader
coverage.

The intended change set is limited to the profile, its focused test,
`CURRENT.md`, the Fix.1 and Fix.2 phase records, the memory manifest, and the
user-created Fix.2 task contract. Do not create a farm-evidence record until
fresh applicable farm artifacts are supplied and independently reviewed.
