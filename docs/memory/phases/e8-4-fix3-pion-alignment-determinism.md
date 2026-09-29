# E.8.4.Fix.3 — pion-alignment determinism repair

## Status

`SOURCE REVIEWED` — independent ChatGPT actual-diff/source review and
pushed-state review passed for the user-controlled pushed source
`29d7b7f9635db899939efeb3508e941e994e8928` (`E8.4 Fix.3: repair pion
alignment determinism`). This is not ROOT/PyROOT, full-analysis, procedure-PDF,
farm, or runtime validation.

## Starting identity and narrow scope

- Required and observed branch: `test`.
- Required and observed starting HEAD:
  `9daca79051a4eb451a158dd7c34b7028c6b96456` (`Re-pin E8.4 validation bundle
  profile`).
- Substantive source scope: `src/cuts/pion_component_fits.py` and
  `testing/test_pion_component_dynamic_alignment.py` only.
- The baseline production branch, F.6.3 producer, E.8.4 consumer, validation
  profile, `background_config.py`, and frozen F.6.2 artifacts remain unchanged.

## Implemented repair

- Pion-control cache compatibility retains generated `hist_name` as diagnostic
  provenance but compares only a valid SHA-256 content checksum and complete,
  valid axis specification. Missing or malformed semantic identity fails
  closed. Alignment schema version remains v2; all configuration, physical-bin,
  parent, immutable-template, setting, phi, and epsilon compatibility checks
  remain authoritative.
- The candidate scan uses an explicit `rel_tol=abs_tol=1e-12` comparison only
  when deciding whether a shifted-template integral is materially below the
  configured minimum. The configured threshold remains `1.0`; interpolation,
  renormalization, scan grids, windows, scoring, ordering, and production
  pion-subtraction behavior are unchanged.

## Deterministic local coverage

The focused suite includes a ROOT-independent persisted-cache regression for a
renamed but semantically identical pion-control histogram, checksum/axis/
malformed-identity rejection coverage, and direct coverage of the actual scan
boundary helper at exact, `1e-13` below, and materially below `1.0`. Existing
ROOT histogram-path coverage remains separately skipped when PyROOT is absent.

Codex ran these local deterministic checks with `python -B`; they establish no
ROOT/PyROOT, full-analysis, farm, or runtime result:

```text
PASS  py_compile src/cuts/pion_component_fits.py testing/test_pion_component_dynamic_alignment.py
PASS  testing.test_pion_component_dynamic_alignment — 15 tests OK; 10 existing
      histogram-path tests skipped because PyROOT is unavailable
PASS  testing.test_t_bin_pion_parent_integrity — 15 tests OK
PASS  testing.test_binning_pre_particle_subtraction — 16 tests OK
PASS  testing.test_f6_3_parallel_full_procedure_method_a — 21 tests OK
PASS  testing.test_e8_4_production_impact_audit — 13 tests OK
PASS  update_memory_manifest.py --write and --check
WARN  check_memory_health.py exit 0: docs/memory/CURRENT.md exceeds soft limit
      (9742 > 8192)
PASS  memory_bootstrap.py --json and testing.test_memory_health — 35 tests OK
PASS  git -c core.safecrlf=false diff --check
```

## Next

`NEXT` — user-controlled commit/push of the independently reviewed E.8.4.Fix.3
validation-bundle profile re-pin. Do not infer runtime closure from source
review or local implementation.
