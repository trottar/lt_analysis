# E.8.4.Fix.4 — persisted alignment semantic-version gate

## Status

`CLOSED / RUNTIME VALIDATED` — fresh Left/lowe farm evidence directly shows a
current-semantics setting-wide alignment record with
`persistence_status = rejected_stale_then_created` and rejection reason
`alignment_semantics_version mismatch`. This closure covers only stale-cache
rejection and recomputation. F.6.3 and E.8.4 remain unvalidated at runtime;
their fresh gate remains blocked by
`f6_3_f4_shared_reproduction_failed:f3_fingerprint_input_content_mismatch`.
See [fresh Fix.4 evidence](../evidence/e8-4-fix4-left-lowe-runtime-closure-and-f4-baseline-divergence.md).

Prior source review: independent ChatGPT actual-diff/source-runtime-path review
of `kaonlt_review(20260929-215203).diff` passed for the local candidate at
fixed `test` starting HEAD `8d62dcfdf2fd8d08298c075940b4ab1b28d7f079`.
The reviewed source/test diffs remain byte-identical. This review is not
ROOT/PyROOT, full `main.py`, procedure-PDF, farm, runtime, production, or
Method-A-promotion acceptance.

## Source change

`pion_component_fits.py` now owns deterministic resolver semantics version
`pion_component_dynamic_alignment_semantics/v2`, separate from alignment
schema v2 and scientific configuration. Expected metadata, normal and
fallback/disabled resolver results, persisted JSON, and CSV diagnostics carry
it. Cache compatibility requires the exact current version; missing, wrong,
or malformed values reject the old record and invoke the existing resolver.
A current-semantics record still reuses on full compatibility, including a
generated histogram-name change with unchanged pion-control checksum/axis.
Direct fine-bin parent validity also requires current semantics.

Only `src/cuts/pion_component_fits.py` and
`testing/test_pion_component_dynamic_alignment.py` have substantive changes.
The Fix.3 minimum-integral predicate, scientific configuration, scan grids,
thresholds, component order, templates, alignment schema v2, F.3/F.4/F.6.3,
E.8.4 consumer, validation profile, and baseline physics remain untouched.
Independent review confirmed the narrow source/test boundary, exact persisted
semantics rejection, current-semantics reuse, direct parent fail-closed check,
and preserved Fix.3 pion-control checksum+axis identity and minimum-integral
predicate. No Method A or Method B mathematics changed.

## Evidence boundary and next gate

The [fresh Left/lowe blocker](../evidence/e8-4-fix3-left-lowe-stale-alignment-cache-runtime-blocker.md)
combines supplied farm artifact/fingerprint observations with a separate source
audit; it does not directly prove farm cache-reuse status. Deterministic local
checks exercise cache rejection, one-time recomputation, subsequent reuse,
and the existing fail-closed contracts. Histogram-path tests require PyROOT
and are skipped locally when unavailable.

## Deterministic local checks

- `py_compile` passed for the two changed Python files.
- `testing.test_pion_component_dynamic_alignment`: 19 tests OK, 11 existing
  PyROOT histogram-path tests skipped because PyROOT is unavailable. The
  ROOT-independent parent-semantics fallback test passed.
- `testing.test_t_bin_pion_parent_integrity`: 15 tests OK.
- `testing.test_binning_pre_particle_subtraction`: 16 tests OK.
- `testing.test_f6_3_parallel_full_procedure_method_a`: 21 tests OK.
- `testing.test_e8_4_production_impact_audit`: 13 tests OK.
- Memory manifest/health/bootstrap checks, memory-health unit tests, and
  `git diff --check` are recorded in the fresh cumulative review bundle.

These are Codex-reported deterministic results from the reviewed bundle;
they were `NOT RUN by ChatGPT` and establish no runtime validation.

The fresh Left/lowe E.8.4 gate remains `BLOCKED` by
`f6_3_f4_shared_reproduction_failed:f3_fingerprint_input_content_mismatch`.
After final pre-push reconciliation review, the user controls commit/push;
pushed-state review, any required separately reviewed profile re-pin, and a
new narrow farm gate remain subsequent steps.

## Pushed-state provenance

Independent GitHub pushed-state review passed for user-controlled remote
`test` commit `6e2adf7a37ac9e79cad99242686804cf51701644` (`E8.4 Fix.4:
version persisted alignment semantics`), parent
`8d62dcfdf2fd8d08298c075940b4ab1b28d7f079`. The pushed commit contains
exactly the ten reviewed Fix.4 candidate paths, including only the two
substantive source/test files named above; temporary review bundles are absent.
This is source/provenance evidence only. The validation-bundle profile re-pin
is a distinct `ACTIVE` gate and no ROOT/PyROOT, farm, runtime, production, or
Method-A promotion conclusion changes.
