# E.8.4 — validation-bundle profile provenance re-pin

## Status

`SOURCE REVIEWED` — independent ChatGPT reviewed the complete
`kaonlt_review(20260929-100154).diff` bundle-profile re-pin candidate and
returned `PASS`. This is source/provenance work only and does not claim
ROOT/PyROOT, procedure-PDF runtime, Jefferson Lab farm, or production
validation.

## Starting identity and narrow scope

- Required/observed branch: `test`.
- Required/observed committed HEAD:
  `1aa1fd4184a6f8b20043e00ebb1ed3e9505a4935` (`Add E8.4 Method-A
  production-impact audit`).
- The existing E.8.4 source remains `SOURCE REVIEWED`; this record does not
  reopen its renderer, F.6.3 producer, baseline production, or scientific
  ownership.
- Only `testing/pion_hgcer_validation_bundle_profile_e8_1.json` and
  `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py` change outside
  warranted memory/manifest records.

## Re-pin and preserved validation boundary

The E.8.1 generic procedure-PDF profile previously required stale Fix.5
analysis source `53fd262b730af8f1254e411a38231aebeb6a1da3`. It now requires the
pushed E.8.4 source `1aa1fd4184a6f8b20043e00ebb1ed3e9505a4935`. The focused
test has the same exact required-source expectation.

The profile remains schema `pion_hgcer_validation_bundle_profile/v4`, identity
`phase_e8_1_full_background_procedure_pdf_farm_review/v1`, and
`generic_artifacts` mode. Its ordered canonical-five settings, global frozen
F.6.2 acceptance-refinement JSON, and per-setting ordinary full-background
procedure-PDF/page-manifest pair are unchanged. The committed-range rule still
allows only the profile/test pair after the required analysis source, with only
`docs/memory/` as a non-analysis prefix; later analysis changes remain fail
closed. The generic collector and `tcsh` wrapper are unchanged.

No accepted artifact, analysis/renderer source, production physics, Method-A/
Method-B ownership, or farm runtime state changed. No farm command was
prepared or run.

## Codex-reported deterministic local checks

Interpreter: `Python 3.12.10`. These checks were run locally by Codex; they
are not ROOT/PyROOT, procedure-PDF runtime, Jefferson Lab farm, or production
validation.

```text
PASS  python -B -m py_compile testing/test_pion_hgcer_validation_bundle_profile_e8_1.py
PASS  python -m json.tool testing/pion_hgcer_validation_bundle_profile_e8_1.json
PASS  python -B -m unittest testing.test_pion_hgcer_validation_bundle_profile_e8_1 -v
      Ran 5 tests; OK
PASS  python -B -m unittest testing.test_collect_pion_hgcer_validation_bundle -v
      Ran 28 tests; OK
```

## Fix.1 — provenance restoration and memory-integrity evidence

Independent review of `kaonlt_review(20260929-010558).diff` found that the
re-pin had compressed the pushed E.8.3/F.6.3/E.8.4 source-review provenance in
`CURRENT.md`. Fix.1 restores that exact detailed provenance, including the
three reviewed-diff identities, while retaining this separate `ACTIVE` re-pin
state and its no-farm/runtime boundary. It does not alter the reviewed profile
or focused test, the collector/wrapper, USER/roadmap records, or any
analysis/runtime path.

The following commands were run locally by Codex after the restoration; they
do not establish ROOT/PyROOT, procedure-PDF runtime, Jefferson Lab farm, or
production validation:

```text
PASS  python -B tools/update_memory_manifest.py --root . --write
      regenerated docs/memory/manifest.json
PASS  python -B tools/update_memory_manifest.py --root . --check
WARN  python -B tools/check_memory_health.py --root .
      exit 0; CURRENT.md soft-size warning: 8349 > 8192 bytes
PASS  python -B tools/memory_bootstrap.py --root . --json
      exit 0; reported branch test and HEAD 1aa1fd4184a6f8b20043e00ebb1ed3e9505a4935
PASS  python -B -m unittest testing.test_memory_health -v
      Ran 35 tests; OK
PASS  git -c core.safecrlf=false diff --check
```

## Final independent source-review reconciliation

Independent ChatGPT reviewed `kaonlt_review(20260929-100154).diff` and
returned `PASS`. The bundle-profile re-pin is therefore `SOURCE REVIEWED` for
source/provenance scope: the substantive profile/test re-pin passed, and Fix.1
preserved the detailed E.8.3/F.6.3/E.8.4 provenance while recording the
required local checks.

`check_memory_health.py` exited 0 with the known CURRENT soft-size warning
(`8349 > 8192` bytes), which is non-fatal and distinct from a source or runtime
failure. Codex-reported focused and memory checks were `NOT RUN by ChatGPT`.
No ROOT/PyROOT, procedure-PDF, Jefferson Lab farm, runtime, production-physics,
or production-promotion claim exists.

Final reconciliation commands run locally by Codex:

```text
PASS  python -B tools/update_memory_manifest.py --root . --write
      regenerated docs/memory/manifest.json
PASS  python -B tools/update_memory_manifest.py --root . --check
WARN  python -B tools/check_memory_health.py --root .
      exit 0; CURRENT.md soft-size warning: 8685 > 8192 bytes
PASS  python -B tools/memory_bootstrap.py --root . --json
      exit 0; reported branch test and HEAD 1aa1fd4184a6f8b20043e00ebb1ed3e9505a4935
PASS  python -B -m unittest testing.test_memory_health -v
      Ran 35 tests; OK
PASS  git -c core.safecrlf=false diff --check
```

## Next

`NEXT` — user-controlled commit/push of the reviewed bundle-profile re-pin. Do
not provide, prepare, or execute farm-run or packaging commands until this
exact bundle/profile is pushed and pushed-state reviewed; advance one workflow
gate at a time.
