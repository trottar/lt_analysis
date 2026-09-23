# E.8.1.Fix.4 - bundle-profile provenance re-pin

## Status and source identity

`SOURCE REVIEWED` — this profile-only provenance repair starts from committed
`test` HEAD `350c34c55b2de33ad01011559dc6d8ed84d9c8a7` (`Fix E8.1 procedure PDF
layout`). It requires that exact, user-pushed `SOURCE REVIEWED` Fix.3
analysis/procedure source; it is not the future Fix.4 profile/bundle commit.
ChatGPT independently inspected the complete actual diff and found the
profile/test provenance re-pin PASS. Codex-reported deterministic checks were
`NOT RUN by ChatGPT`. No ROOT/PyROOT, farm, or runtime validation is claimed,
and no Fix.4 commit exists yet.

## Narrow profile change and preserved boundary

`testing/pion_hgcer_validation_bundle_profile_e8_1.json` changes only
`source_identity.required_analysis_commit`, from Fix.1 reader source
`8985d9a212799c021c4ad1a759a689ea3826e0ea` to pushed Fix.3 source
`350c34c55b2de33ad01011559dc6d8ed84d9c8a7`. Its focused test changes only the
matching `REVIEWED_SOURCE` identity.

The profile retains schema v4,
`phase_e8_1_full_background_procedure_pdf_farm_review/v1`, generic-artifacts
mode, the ordered five settings, one global frozen F.6.2 JSON, and the PDF plus
page-manifest pair per setting. Its committed-range allowlist remains only the
profile/test pair and its sole non-analysis prefix remains `docs/memory/`.
Later analysis-source changes remain fail-closed; no `src/` exception is added.

The Fix.3 renderer/test, collector, wrapper, artifact declarations, canonical
settings, accepted F.6.2 evidence, and all production/scientific boundaries are
unchanged. E.8.1 remains presentation-only; Method A remains detached and
non-production, Method B remains diagnostic only. F.6.2 and F.6.2.Fix.5 remain
`CLOSED / RUNTIME VALIDATED`; E.8.1 remains `ACTIVE`; Fix.1/Fix.2 retain their
narrow closures; F.6.3 remains `BLOCKED`.

## Deterministic local validation

- `python -B -m py_compile testing/test_pion_hgcer_validation_bundle_profile_e8_1.py` — PASS.
- `python -m json.tool testing/pion_hgcer_validation_bundle_profile_e8_1.json` — PASS.
- `python -B -m unittest testing.test_pion_hgcer_validation_bundle_profile_e8_1 -v` — PASS, 5 tests, no skips.
- `python -B -m unittest testing.test_collect_pion_hgcer_validation_bundle -v` — PASS, 28 tests, no skips.
- Manifest write/check, memory health, memory bootstrap, and `git diff --check`
  — PASS.

These are local source/provenance checks only, not ROOT/PyROOT or farm results.

## Next action

Next: user-controlled commit/push of this `SOURCE REVIEWED` profile/test/memory
change, followed by ChatGPT pushed-state review. That future pushed Fix.4
identity becomes the wrapper `--bundle-commit`, but the profile keeps requiring
the distinct pushed Fix.3 source above. Only then run the fresh `Q4p4W2p74 /
Left / lowe` gate, collect its bundle/PDF, and independently inspect whether
the exact Fix.3 context clipping, mojibake, handoff clipping, and overlay-header
defects are repaired before returning to detailed `Left / highe` or broader
coverage.
