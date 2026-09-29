# E.8.4.Fix.3 — validation-bundle profile re-pin

## Status

`SOURCE REVIEWED` — final independent ChatGPT actual-diff/source-provenance
review of `kaonlt_review_e8_4_fix3_profile_repin.diff` `PASSED`. The reviewed
candidate is not ROOT/PyROOT, procedure-PDF, full-analysis, farm, runtime,
production, or Method-A-promotion evidence.

## Starting identity and scope

- Required and observed branch: `test`.
- Required and observed starting HEAD:
  `29d7b7f9635db899939efeb3508e941e994e8928` (`E8.4 Fix.3: repair pion
  alignment determinism`).
- The pushed Fix.3 source is `SOURCE REVIEWED`: its independent actual-diff/
  source review and pushed-state review passed. Its Left/lowe observation
  remains blocker evidence, not runtime closure.
- Only `testing/pion_hgcer_validation_bundle_profile_e8_1.json` and
  `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py` change outside
  warranted durable memory and manifest records.

## Re-pin and preserved boundary

The v4 generic E.8.1 procedure-PDF profile and focused test now require pushed
Fix.3 analysis source `29d7b7f9635db899939efeb3508e941e994e8928` rather than
the prior E.8.4 source `1aa1fd4184a6f8b20043e00ebb1ed3e9505a4935`.

Schema, profile identity, collection mode, ordered five-setting inventory,
frozen F.6.2 JSON, per-setting procedure-PDF/page-manifest artifacts, exact
profile/test-only committed-file allowlist, and `docs/memory/` non-analysis
prefix remain unchanged. The generic collector and `tcsh` wrapper are
byte-unchanged. No artifact is packaged and no farm action is performed.

## Deterministic local checks

Codex ran the following local checks; they do not establish ROOT/PyROOT,
procedure-PDF, full-analysis, farm, or runtime validation:

```text
PASS  py_compile testing/test_pion_hgcer_validation_bundle_profile_e8_1.py
      testing/test_pion_component_dynamic_alignment.py
PASS  json.tool testing/pion_hgcer_validation_bundle_profile_e8_1.json
PASS  testing.test_pion_hgcer_validation_bundle_profile_e8_1 — 5 tests OK
PASS  testing.test_collect_pion_hgcer_validation_bundle — 28 tests OK
PASS  testing.test_pion_component_dynamic_alignment — 15 tests OK; 10 existing
      histogram-path tests skipped because PyROOT is unavailable
PASS  update_memory_manifest.py --write and --check
WARN  check_memory_health.py exit 0: docs/memory/CURRENT.md exceeds soft limit
      (10419 > 8192)
PASS  memory_bootstrap.py --json and testing.test_memory_health — 35 tests OK
PASS  git -c core.safecrlf=false diff --check
```

## Final pre-push source-review reconciliation

Independent ChatGPT actual-diff/source-provenance review of
`kaonlt_review_e8_4_fix3_profile_repin.diff` returned `PASS`. It establishes
`SOURCE REVIEWED` for this profile/test provenance re-pin only: the required
analysis source is exactly the pushed Fix.3 commit
`29d7b7f9635db899939efeb3508e941e994e8928`, while the generic collector,
`tcsh` wrapper, and scientific/runtime source remain frozen. Codex-reported
deterministic tests were `NOT RUN by ChatGPT`. This review does not establish
ROOT/PyROOT, procedure-PDF, full-analysis, farm, runtime, production, or
Method-A-promotion evidence.

For this final pre-push reconciliation, Codex regenerated and checked the
memory manifest, ran memory health and bootstrap, ran
`testing.test_memory_health` (35 tests OK), and ran
`git -c core.safecrlf=false diff --check`. Memory health exited 0 with the
non-fatal warning `docs/memory/CURRENT.md exceeds soft limit (10795 > 8192)`.

## Next

`NEXT` — user-controlled commit/push of the independently reviewed E.8.4.Fix.3
validation-bundle profile re-pin. Do not infer runtime closure, production
acceptance, or Method-A promotion from source review or local deterministic
checks.
