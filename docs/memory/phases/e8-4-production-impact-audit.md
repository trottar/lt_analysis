# E.8.4 — Baseline-versus-Method-A production-impact audit

## Status

`SOURCE REVIEWED` — independent ChatGPT actual-diff/source-runtime-path review
passed the complete cumulative `kaonlt_review(20260928-233040).diff` candidate.
This is source review only. It does not claim ROOT/PyROOT, full-analysis,
procedure-PDF runtime, Jefferson Lab farm, or production validation.

## Starting identity and scope

- Required/observed branch: `test`.
- Required/observed committed starting HEAD:
  `d84f9610427cf38ccaaa8d6d8b14a8a995ecbac9` (`Add F6.3 parallel Method-A
  full procedure`).
- The sole pre-existing worktree item was the user-created E.8.4 task contract.
- Source scope: `src/cuts/full_background_subtraction_plots.py` only.
- Focused regression: `testing/test_e8_4_production_impact_audit.py`.
- Warranted history scope: `CURRENT.md`, the E.8 decision roadmap, the Phase-F
  roadmap, roadmap status, this record, and the regenerated integrity manifest.

## Implemented presentation-only consumer

E.8.4 reads only the already-produced post-yield
`hist["_f6_3_parallel_method_a_source"]` in the existing historically named
`finalize_full_background_subtraction_e8_2(...)` transaction. The finalizer
first rebuilds the existing E.8.2 presentation payload, then passes the F.6.3
source plus that payload into the E.8.4 consumer, and sends the resulting
available or unavailable payload through the existing temporary-PDF/page-
manifest render, recovery, and install path. No second finalizer, PDF lifecycle,
tree traversal, factor lookup, template fill, fit, or yield calculation was
added. `main.py` and `rand_sub.py` remain unchanged.

The consumer accepts only `f6_3_parallel_method_a_source/v1` with an available
current setting, exact parallel non-production branch role, present authority
whose live-cache parity explicitly passed, and all required non-production flags.
It cross-checks the selected setting, Lambda window, exact canonical `(t,phi)`
inventory and geometry, and F.6.3 `Y0`/statistical/total-error values against
the already-built E.8.2 baseline presentation. Any missing, duplicate, stale,
malformed, nonfinite, or mismatched input is a literal unavailable E.8.4
payload; it is never replaced by E.8.3, F.5/F.6.1, unity factors, or baseline
placeholders.

Every retained `pion_input`, `B_pi_0`, `B_pi_A`, `MM_0`, and `MM_A` is a
detached display clone. Display helpers make fresh clones to show the stored
direct signed consequences `Delta B_pi = B_pi_A - B_pi_0` and
`Delta K = MM_A - MM_0`; neither is written into the F.6.3 source or used
downstream. The payload copies producer-owned `Y0`, `YA`, and their existing
statistical/total errors and derives only `DeltaY = YA - Y0` and, for finite
nonzero `Y0`, `DeltaY/Y0`. A zero baseline yield carries the explicit
`y0_zero_fraction_undefined` reason. No epsilon denominator, cap, clipping,
normalization, rebinning, child renormalization, interpolation, or DeltaY/
fractional-shift uncertainty is introduced.

## Procedure-PDF route

The optional keyword-only E.8.4 renderer argument preserves the historical
page inventory for callers that omit it. On the post-yield route its pages are
inserted after E.8.3 and before the terminal E.8 handoff:

- setting authority/context;
- one pion-consequence page per canonical t parent;
- one final-MM `MM_0`/`MM_A` page per t parent, including the common Lambda
  window and per-child stored yields;
- one direct signed-difference page per t parent;
- one stored yield-impact page per t parent, with points-only fractional plots
  so an undefined denominator is never bridged; and
- one setting-level canonical t-by-phi `DeltaY`/`DeltaY/Y0` summary without a
  constructed t-integrated observable.

The stable manifest IDs are `full_background.e8_4.authority`, `.pion_consequence`,
`.final_mm`, `.signed_difference`, `.yield_impact`, `.setting_summary`, and
`.unavailable`. The unavailable page preserves the literal sidecar reason and
does not invalidate the baseline result, E.8.2, or E.8.3.

The E.8/E.8.3/handoff wording now states the settled ownership accurately:
E.8.3 remains detached evidence; F.6.3 alone constructs the private branch;
E.8.4 consumes it for the actual impact audit; Method A is not promoted; and
F.6.4 remains the sole promotion decision.

## Frozen boundaries

The following are byte-unchanged from the starting source audit:
`src/main.py`, `src/binning/calculate_yield.py`, `src/cuts/rand_sub.py`,
`src/cuts/pion_component_subtraction.py`,
`src/cuts/pion_hgcer_method_a_parallel_full_procedure.py`, and
`src/utility/background_config.py`. The baseline production branch, F.6.3
producer, public `Y0`, random/dummy and slow-proton treatment, pruning,
pion components/fits/windows/templates and baseline `w0`, Method B, canonical
binning, normalizations, SIMC, efficiencies, cross sections, accepted F.6.2
artifacts, active `no_empirical_residual` profile, and dormant Fit 1/Fit 2 all
remain frozen.

## Codex-reported deterministic local checks

Interpreter: `python` (Python 3.12 on this host). These are local deterministic
checks run by Codex; they were **NOT RUN by ChatGPT** and are not ROOT/PyROOT,
full-analysis, procedure-PDF runtime, farm, or runtime validation.

```text
PASS  python -B -m py_compile src/cuts/full_background_subtraction_plots.py testing/test_e8_4_production_impact_audit.py
PASS  python -B -m unittest testing.test_e8_4_production_impact_audit -v
      Ran 9 tests; OK
PASS  python -B -m unittest testing.test_e8_2_baseline_stage_audit -v
      Ran 23 tests; OK
PASS  python -B -m unittest testing.test_e8_3_detached_method_a_reweighting_audit -v
      Ran 21 tests; OK
PASS  python -B -m unittest testing.test_f6_3_parallel_full_procedure_method_a -v
      Ran 21 tests; OK
PASS  python -B -m unittest testing.test_full_background_subtraction_plots -v
      Ran 102 tests; OK (18 expected PyROOT/superseded-tail skips)
PASS  python -B -m unittest testing.test_pion_hgcer_phase_e_runtime_contract -v
      Ran 3 tests; OK
PASS  python -B -m unittest testing.test_pion_hgcer_phase_f_runtime_contract -v
      Ran 3 tests; OK
```

```text
PASS  python -B tools/update_memory_manifest.py --root . --write
PASS  python -B tools/update_memory_manifest.py --root . --check
PASS  python -B tools/check_memory_health.py --root .
PASS  python -B -m unittest testing.test_memory_health -v
      Ran 35 tests; OK
PASS  python -B tools/memory_bootstrap.py --root . --json
      reports test at d84f9610427cf38ccaaa8d6d8b14a8a995ecbac9 and only
      the intended E.8.4 candidate paths as dirty
```

No farm action was run.

## Next

`NEXT` — independent ChatGPT review of the complete E.8.4 actual diff/source/runtime path.
Do not request a farm run until that source review passes. Final E.8 and F.6.4
remain `BLOCKED`; no Method-A production promotion is implied.

## E.8.4.Fix.1 — independent source-review repair

`ACTIVE` at implementation time — this historical narrow repair followed
independent review of the original local candidate. Its findings are now
incorporated into the cumulative E.8.4 `SOURCE REVIEWED` closure below; it does
not remain a separate active work item.

The reviewed candidate had three presentation-consumer defects, repaired only
in `src/cuts/full_background_subtraction_plots.py` with focused coverage in
`testing/test_e8_4_production_impact_audit.py`:

- E.8.2 retains semantic runtime epsilon `low`/`high`, while the frozen F.6.3
  selected-setting token uses `lowe`/`highe`. The consumer now explicitly maps
  only `low -> lowe` and `high -> highe` for the comparison, retains the
  semantic display epsilon, and fails closed for every other token.
- E.8.2's wide diagnostic-MM geometry remains distinct from F.6.3's narrow
  analysis/yield-MM geometry. The consumer now validates finite, strictly
  ordered, matching F.6.3 edges inside every child and across all canonical
  children, then exposes that analysis-MM edge vector without rebinning either
  source.
- The pion-consequence page now draws the common proton-cleaned `pion_input`
  first and overlays both `B_pi_0` and `B_pi_A` with explicit `same` semantics;
  E.8.4 renderer absence checks use explicit `is None` logic.

The strengthened focused suite covers low/high setting-token conversion,
unsupported epsilon failure, wide-E.8.2/narrow-F.6.3 availability, within- and
cross-child analysis-MM edge failures, baseline statistical/total-error
identity failures, detached source contents/errors through every E.8.4
renderer, common-input draw order, available and unavailable optional-sidecar
finalization, rollback/recovery, zero-Y0 gap safety, and existing static
ownership guards.

Codex ran these local deterministic checks with `Python 3.12.10`; they were
**NOT RUN by ChatGPT** and are not ROOT/PyROOT, full-analysis, procedure-PDF
runtime, farm, or runtime validation:

```text
PASS  python -B -m py_compile src/cuts/full_background_subtraction_plots.py testing/test_e8_4_production_impact_audit.py
PASS  python -B -m unittest testing.test_e8_4_production_impact_audit -v
      Ran 12 tests; OK
PASS  python -B -m unittest testing.test_e8_2_baseline_stage_audit -v
      Ran 23 tests; OK
PASS  python -B -m unittest testing.test_e8_3_detached_method_a_reweighting_audit -v
      Ran 21 tests; OK
PASS  python -B -m unittest testing.test_f6_3_parallel_full_procedure_method_a -v
      Ran 21 tests; OK
PASS  python -B -m unittest testing.test_full_background_subtraction_plots -v
      Ran 102 tests; OK (18 expected PyROOT/superseded-tail skips)
PASS  python -B -m unittest testing.test_pion_hgcer_phase_e_runtime_contract -v
      Ran 3 tests; OK
PASS  python -B -m unittest testing.test_pion_hgcer_phase_f_runtime_contract -v
      Ran 3 tests; OK
```

No frozen producer/runtime file or F.6.3 scientific behavior changed. Final E.8
and F.6.4 remain `BLOCKED`; Method A remains unpromoted.

`NEXT` — independent ChatGPT review of the refreshed complete cumulative
E.8.4 + Fix.1 diff/source/runtime path. Do not request a farm run until that
source review passes.

## E.8.4.Fix.2 — malformed-sidecar fail-closed repair

`ACTIVE` at implementation time — this historical narrow repair was pending
independent ChatGPT review of the refreshed cumulative E.8.4 + Fix.1 + Fix.2
diff/source/runtime path. It is incorporated into the cumulative E.8.4 `SOURCE
REVIEWED` closure below and does not remain a separate active work item.

Independent source review found that a nominally available F.6.3 sidecar with
an integer or mapping in `lambda_integration_window` or `children` could escape
the existing E.8.4 `_E8PayloadError` boundary while the builder coerced the raw
external value to a tuple. The consumer now first requires each external
container to be a non-string `Sequence`; malformed Lambda-window containers
return literal `e8_4_lambda_window_invalid`, and malformed child inventories
return literal `e8_4_child_inventory_invalid`, through the existing explicit
unavailable-payload flow. No blanket exception catch was added.

Focused coverage proves both integer and iterable mapping forms fail closed for
each container. It also passes a malformed but schema-valid available sidecar
through the existing finalizer: the optional E.8.4 unavailable page renders,
the baseline transaction completes with `status=available`,
`e8_4_available=False`, and the literal malformed-sidecar reason, while the
same E.8.2 payload reaches the renderer unchanged. The accepted Fix.1 setting
mapping, distinct wide/narrow MM geometry, common pion-input draw order,
rollback/recovery, zero-Y0, and static ownership coverage remain intact.

Codex ran these local deterministic checks with `Python 3.12.10`; they were
**NOT RUN by ChatGPT** and are not ROOT/PyROOT, full-analysis, procedure-PDF
runtime, farm, or runtime validation:

```text
PASS  python -B -m py_compile src/cuts/full_background_subtraction_plots.py testing/test_e8_4_production_impact_audit.py
PASS  python -B -m unittest testing.test_e8_4_production_impact_audit -v
      Ran 13 tests; OK
PASS  python -B -m unittest testing.test_e8_2_baseline_stage_audit -v
      Ran 23 tests; OK
PASS  python -B -m unittest testing.test_e8_3_detached_method_a_reweighting_audit -v
      Ran 21 tests; OK
PASS  python -B -m unittest testing.test_f6_3_parallel_full_procedure_method_a -v
      Ran 21 tests; OK
PASS  python -B -m unittest testing.test_full_background_subtraction_plots -v
      Ran 102 tests; OK (18 expected PyROOT/superseded-tail skips)
PASS  python -B -m unittest testing.test_pion_hgcer_phase_e_runtime_contract -v
      Ran 3 tests; OK
PASS  python -B -m unittest testing.test_pion_hgcer_phase_f_runtime_contract -v
      Ran 3 tests; OK
```

No frozen producer/runtime file or F.6.3 scientific behavior changed. Final E.8
and F.6.4 remain `BLOCKED`; Method A remains unpromoted.

`NEXT` — independent ChatGPT review of the refreshed complete cumulative
E.8.4 + Fix.1 + Fix.2 diff/source/runtime path.

## Final pre-push source-review reconciliation

`SOURCE REVIEWED` — independent ChatGPT actual-diff/source-runtime-path review
of `kaonlt_review(20260928-233040).diff` returned **PASS** for the complete
cumulative E.8.4 + Fix.1 + Fix.2 candidate. Fix.1's semantic `low/high` to
`lowe/highe` comparison identity, separate wide E.8.2/narrow F.6.3 analysis-MM
geometry, and visible common pion input, plus Fix.2's malformed Lambda-window
and child-inventory fail-closed boundary, are resolved in that reviewed
cumulative candidate.

The review confirms that E.8.4 remains a presentation-only consumer of the
already-produced F.6.3 `_f6_3_parallel_method_a_source`; the public baseline
branch is unchanged; F.6.3 alone owns the private parallel Method-A branch;
Method B is numerically absent; and legacy empirical residual Fit 1/Fit 2
remain dormant under `no_empirical_residual`. It also confirms detached retained
histograms, stored producer-owned `Y0`/`YA` display comparisons without an
invented Method-A-shift uncertainty, no child renormalization/reconstruction/
tree traversal/template refill/second yield calculation, and the existing
single post-yield pair-safe PDF/page-manifest transaction with E.8.4 pages after
E.8.3 and before the terminal handoff. No production promotion is performed.

Codex-reported local deterministic checks were **NOT RUN by ChatGPT**. This
review claims no ROOT/PyROOT, rendered-PDF, full-analysis, Jefferson Lab farm,
or runtime acceptance. Final E.8 and F.6.4 remain `BLOCKED`; Method A remains
unpromoted.

`NEXT` — user-controlled commit/push of the reviewed cumulative E.8.4 + Fix.1
+ Fix.2 + final source-review memory reconciliation
-> ChatGPT pushed-state review
-> one narrow Jefferson Lab farm gate for Q4p4W2p74 / Left / lowe
-> fresh artifacts
-> ChatGPT evidence review.
