# E.8.1.Fix.4 - bundle-profile provenance re-pin

## Status and source identity

CLOSED / RUNTIME VALIDATED — narrow profile-provenance closure only. This
profile-only repair started from committed test HEAD
350c34c55b2de33ad01011559dc6d8ed84d9c8a7, the distinct user-pushed
SOURCE REVIEWED Fix.3 analysis/procedure source that it requires. ChatGPT
independently inspected the complete profile/test diff and found it PASS;
Codex-reported deterministic checks were NOT RUN by ChatGPT.

The user pushed Fix.4 as fdd368f2084c9b9508e7ce679b8f7391f5b556f1. Its fresh
Q4p4W2p74 / Left / lowe farm bundle passed the profile provenance gate: the
required Fix.3 source is an ancestor, no unexpected committed files are
admitted, the captured source worktree is clean, and the bundle manifest is
complete with no errors. This is a runtime validation of Fix.4's narrow
provenance boundary, not E.8.1 acceptance or an overlay-PDF visual pass.

## Narrow profile change and preserved boundary

testing/pion_hgcer_validation_bundle_profile_e8_1.json changes only
source_identity.required_analysis_commit, from Fix.1 reader source
8985d9a212799c021c4ad1a759a689ea3826e0ea to pushed Fix.3 source
350c34c55b2de33ad01011559dc6d8ed84d9c8a7. Its focused test changes only the
matching REVIEWED_SOURCE identity.

The profile retains schema v4,
phase_e8_1_full_background_procedure_pdf_farm_review/v1, generic-artifacts
mode, the ordered five settings, one global frozen F.6.2 JSON, and the PDF plus
page-manifest pair per setting. Its committed-range allowlist remains only the
profile/test pair and its sole non-analysis prefix remains docs/memory/.
Later analysis-source changes remain fail-closed; no src/ exception is added.

The Fix.3 renderer/test, collector, wrapper, artifact declarations, canonical
settings, accepted F.6.2 evidence, and all production/scientific boundaries are
unchanged. E.8.1 remains presentation-only; Method A remains detached and
non-production, Method B remains diagnostic only. F.6.2 and F.6.2.Fix.5 remain
CLOSED / RUNTIME VALIDATED; E.8.1 remains ACTIVE; Fix.1/Fix.2 retain their
narrow closures; F.6.3 remains BLOCKED.

## Deterministic local validation

- python -B -m py_compile testing/test_pion_hgcer_validation_bundle_profile_e8_1.py — PASS.
- python -m json.tool testing/pion_hgcer_validation_bundle_profile_e8_1.json — PASS.
- python -B -m unittest testing.test_pion_hgcer_validation_bundle_profile_e8_1 -v — PASS, 5 tests, no skips.
- python -B -m unittest testing.test_collect_pion_hgcer_validation_bundle -v — PASS, 28 tests, no skips.
- Manifest write/check, memory health, memory bootstrap, and git diff --check
  — PASS.

These local source/provenance checks and the narrow farm provenance closure do
not claim ROOT/PyROOT validation of the Fix.3/Fix.5 renderer repair.

## Current successor

The same fresh Fix.4 bundle passes context, handoff, and map pages but blocks
persisted-overlay pages 38/41/44. The details are in [the Fix.4 overlay blocker
evidence](../evidence/e8-1-fix4-left-lowe-overlay-blocker.md). E.8.1.Fix.5 is
ACTIVE for that sole geometry repair. After source review and user push, it
requires a separate narrow profile re-pin to its future pushed source before a
fresh Left/lowe farm gate; detailed Left / highe and broader coverage remain
gated on passing that repaired overlay inspection.
