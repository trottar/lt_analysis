# E.8.1.Fix.6 - bundle-profile provenance re-pin to pushed Fix.5 source

## Status and source identity

SOURCE REVIEWED — this narrow validation-profile repair starts from committed
`test` HEAD 53fd262b730af8f1254e411a38231aebeb6a1da3, the user-pushed,
independently SOURCE REVIEWED E.8.1.Fix.5 overlay-renderer source. ChatGPT
independently inspected the complete Fix.6 actual diff and found the
profile/test provenance re-pin PASS. Codex-reported deterministic checks were
NOT RUN by ChatGPT. No Fix.6 commit exists yet because the user has not
committed/pushed it, and no ROOT/PyROOT, Jefferson Lab farm, or runtime
validation is claimed.

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

## Farm boundary and successor

After user-controlled commit/push and ChatGPT pushed-state review, the future
pushed Fix.6 commit must be passed to the reviewed wrapper as `--bundle-commit`.
The wrapper uses its detached worktree with this E.8.1 profile; the generic
collector then applies exact required-source ancestry and committed-range checks
against the distinct Fix.5 source above before it collects the existing
required artifacts.

Only then may a fresh Q4p4W2p74 / Left / lowe gate be run and its bundle/PDF be
inspected on pages 37--47, especially persisted-overlay pages 38/41/44. A
Left/lowe visual PASS is required before returning to detailed Left / highe;
broader coverage remains gated on that later detailed result. No farm evidence
is created by this task.

## Next action

NEXT — user-controlled commit/push of the reviewed Fix.6 profile/test/memory
set, then ChatGPT pushed-state review, then the fresh Q4p4W2p74 / Left / lowe
farm/PDF validation. The future pushed Fix.6 commit must be the wrapper
`--bundle-commit`, while the profile continues to require the distinct Fix.5
source 53fd262b730af8f1254e411a38231aebeb6a1da3.
