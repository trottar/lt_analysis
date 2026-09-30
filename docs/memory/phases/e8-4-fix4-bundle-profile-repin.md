# E.8.4.Fix.4 — validation-bundle profile re-pin

## Status

`SOURCE REVIEWED` — independent ChatGPT actual-diff/source-provenance review
of `kaonlt_review(20260929-222602).diff` passed. This is not ROOT/PyROOT,
full `main.py`, procedure-PDF, farm, runtime, production, or
Method-A-promotion evidence.

## Starting source and scope

- Required and observed starting `test` HEAD:
  `6e2adf7a37ac9e79cad99242686804cf51701644` (`E8.4 Fix.4: version
  persisted alignment semantics`), parent
  `8d62dcfdf2fd8d08298c075940b4ab1b28d7f079`.
- Independent pushed-state review passed for that exact Fix.4 source and its
  ten reviewed committed paths. Fix.4 remains `SOURCE REVIEWED`; the fresh
  Left/lowe E.8.4 farm gate remains `BLOCKED`.
- Only `testing/pion_hgcer_validation_bundle_profile_e8_1.json` and
  `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py` change
  substantively.

## Source-provenance change

The profile's `source_identity.required_analysis_commit` and the focused
test's `REVIEWED_SOURCE` now require pushed Fix.4 analysis source
`6e2adf7a37ac9e79cad99242686804cf51701644`, replacing Fix.3
`29d7b7f9635db899939efeb3508e941e994e8928`.

The profile schema, profile name, collection mode, canonical-five setting
order, F.6.2 frozen JSON basename, global/setting artifact inventory,
exact profile/test-only `allowed_committed_files`, and sole
`docs/memory/` non-analysis prefix remain unchanged. Collector, `tcsh`
wrapper, Fix.4 source/test, other scientific/runtime source, and accepted
artifacts are unchanged. No farm bundle was collected.

## Deterministic local checks

- `py_compile` passed for the focused profile test and generic collector.
- `json.tool` passed for the E.8.1 validation profile.
- `testing.test_pion_hgcer_validation_bundle_profile_e8_1`: 5 tests OK.
- `testing.test_collect_pion_hgcer_validation_bundle`: 28 tests OK. Its
  synthetic fixtures do not establish farm collection.
- `testing.test_pion_component_dynamic_alignment`: 19 tests OK, 11
  PyROOT histogram-path tests skipped because PyROOT is unavailable.
- Manifest, memory-health/bootstrap, memory-health unit tests, and
  `git diff --check` results are in the fresh cumulative review bundle.

These are Codex-reported local deterministic checks, `NOT RUN by ChatGPT`.
They establish no ROOT/PyROOT, full-analysis, procedure-PDF, farm, or runtime
validation.

## Independent source-provenance review

The passing ChatGPT review confirmed that the profile's
`source_identity.required_analysis_commit` and the focused test's
`REVIEWED_SOURCE` alone change to pushed Fix.4
`6e2adf7a37ac9e79cad99242686804cf51701644`. The two reviewed diffs remain
byte-identical. Canonical-five settings, global/setting artifacts, the exact
profile/test committed-file allowlist, and sole `docs/memory/` non-analysis
prefix are unchanged. Collector, wrapper, Fix.4 analysis/test source, and
scientific/runtime behavior remain unchanged. The fresh Left/lowe E.8.4 gate
remains `BLOCKED` by
`f6_3_f4_shared_reproduction_failed:f3_fingerprint_input_content_mismatch`.
User-controlled commit/push and pushed-state review precede any new narrow
farm evidence; source review does not close the runtime blocker.
