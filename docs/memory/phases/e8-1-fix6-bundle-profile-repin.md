# E.8.1.Fix.6 - bundle-profile provenance re-pin to pushed Fix.5 source

## Status and source identity

CLOSED / RUNTIME VALIDATED — this narrow validation-profile repair started
from committed `test` HEAD 53fd262b730af8f1254e411a38231aebeb6a1da3, the
distinct Fix.5 analysis/procedure source. The user pushed Fix.6 as
0ec29d4e1bb345eb37e8cca35b8b7e5cbe1b4d5b. Its fresh reviewed `Q4p4W2p74 /
Left / lowe` farm bundle passed source provenance, frozen-artifact, checker,
ordinary-PDF/page-manifest, real-PyROOT regression, and pages 37--47 visual
gates. This closure is limited to Fix.6 profile/provenance; it is not canonical-five
E.8.1 acceptance or a production change.

The profile and its focused test re-pin the exact required analysis source from
the pushed Fix.3 procedure renderer
350c34c55b2de33ad01011559dc6d8ed84d9c8a7 to the pushed Fix.5 overlay renderer
53fd262b730af8f1254e411a38231aebeb6a1da3. The future Fix.6 profile/bundle
commit is unknown and remains distinct from this required Fix.5 source.

## Preserved fail-closed boundary

The profile retains schema v4, its existing validation-profile identity and
generic-artifacts collection mode, the ordered canonical five settings, one
global frozen F.6.2 JSON artifact, and the ordinary procedure-PDF plus
page-manifest pair for each setting. Its committed-range allowlist remains
only the profile/test pair; `docs/memory/` remains the sole allowed
non-analysis prefix. Later analysis-source changes, including the Fix.5
renderer path, remain fail-closed rather than being admitted to the allowlist.

This task changes validation source provenance only. It does not change the
Fix.5 renderer or either of its tests, the E.8 reader/authority, artifacts,
canonical settings, generic collector, wrapper, accepted F.6.2 evidence,
production physics, Method-A/Method-B ownership, cuts, normalization,
templates, binning, or persisted numerical content.

## Deterministic local validation

- `python -B -m py_compile testing/test_pion_hgcer_validation_bundle_profile_e8_1.py` — PASS.
- `python -m json.tool testing/pion_hgcer_validation_bundle_profile_e8_1.json` — PASS.
- `python -B -m unittest testing.test_pion_hgcer_validation_bundle_profile_e8_1 -v` — PASS, 5 tests, no skips.
- `python -B -m unittest testing.test_collect_pion_hgcer_validation_bundle -v` — PASS, 28 tests, no skips.

These Codex-reported local source/provenance checks were NOT RUN by ChatGPT and
cannot validate ROOT/PyROOT, the Jefferson Lab farm, or the rendered
persisted-overlay pages.

## Accepted runtime evidence

The bundle/profile commit `0ec29d4e1bb345eb37e8cca35b8b7e5cbe1b4d5b` was the
wrapper `--bundle-commit`; its profile required the distinct Fix.5 source
above. The detached-worktree collector found required-source ancestry, no
unexpected committed files, a clean captured worktree, and complete artifacts
with no errors. The accepted frozen F.6.2 JSON SHA-256 and both fingerprints
were preserved. The procedure PDF has 47 pages with `renderer_failures=[]`.
Independent inspection accepted context, overlays 38/41/44, all map pages, and
handoff; details are in the [Fix.6 Left/lowe runtime closure](../evidence/e8-1-fix6-left-lowe-runtime-closure.md).

## Next action

NEXT — the remaining E.8.1 canonical-five expansion is DEFERRED by user
decision. E.8.2 baseline full-analysis stage audit is the active successor;
Fix.6 closure does not authorize broader farm coverage or production changes.
