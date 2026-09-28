# F.6.3 — parallel full procedure plus Method A

**Status:** `SOURCE REVIEWED` — independent ChatGPT actual-diff/source-runtime-
path review passed the complete cumulative F.6.3 + Fix.1 + Fix.2 candidate.
This record makes no ROOT/PyROOT, full-analysis, farm, or runtime acceptance
claim.

## Scope and starting identity

Implemented against pushed `test` HEAD
`c9ed0d6b4d0013f7475eedbfca57a096730ea840`, under
`f6-3-parallel-full-procedure-method-a-task-contract.md` and the governing
E.8 full-analysis procedure roadmap.

Only the private F.6.3 Method-A side branch was added. E.8.4 is `NEXT`; its
rendering and finalization remain out of scope for this record.

## Actual source path

- `src/cuts/pion_hgcer_method_a_parallel_full_procedure.py` loads the exact
  accepted F.1/F.3/F.4 files through existing filename helpers, validates them
  with existing F.4/F.5 validators, reproduces the F.4 aggregate through the
  shared F.4 calculator, and constructs a selected-setting transient
  `(source_label, entry_index) -> C_j` map. It does not import ROOT or traverse
  analysis trees.
- Before any private fill, the current authoritative yield cache must exactly
  match accepted F.1 identity, canonical `(t,phi)`, `analysis_MM`, `analysis_t`,
  signed source coefficient, `w0`, and signed baseline contribution. Missing,
  extra, duplicate, non-finite, or mismatched rows make the full selected
  setting unavailable. The transient event factors are never persisted.
- `fill_simc_shape_pion_subtraction_templates(...)` has one optional
  `method_a_event_multipliers` argument. With it omitted, the public baseline
  event algebra and output are unchanged. With it supplied, the private path
  changes exactly `w0_j -> w0_j * C_j` for both all-cuts and no-MM-cuts pion
  template filling; missing or invalid factors fail rather than substituting
  `C_j = 1`.
- After the existing public baseline yield extraction,
  `calculate_yield_data(...)` attaches the private
  `_f6_3_parallel_method_a_source` sidecar. It reuses the already-completed
  authoritative cache and existing pion-control inputs; it adds no second
  ROOT/tree traversal and does not alter public return values, `Y0`, E.8.2
  captures, component payloads, or baseline ROOT objects.
- The private source contains detached `B_pi^0`, `B_pi^A`, `MM_0`, `MM_A`,
  `Y0`, `YA`, and their current statistical/total errors with the identical
  Lambda integration window for later E.8.4 consumption. An authority,
  profile, mode, cache, child, or clone failure marks the whole selected setting
  unavailable while leaving baseline production intact.

## Frozen scientific and runtime boundaries

- The accepted F.4 shared calculation and F.1 identities are reused; F.3/F.4
  mathematics is not reimplemented.
- The baseline public procedure, production weights, pion control fits,
  normalizations, canonical binning, yields, SIMC, efficiencies, and accepted
  F.1–F.6.2 artifacts remain unchanged.
- No canonical child normalization, factor clipping/smoothing/interpolation,
  Method-B numerical dependency, empirical Fit 1/Fit 2 activation, or Method-A
  production promotion was added. The active profile remains
  `no_empirical_residual`.

## Codex deterministic checks

Passed locally with `python -B`:

- `py_compile` for the new helper, changed filler, changed yield producer, and
  focused F.6.3 test.
- `testing.test_f6_3_parallel_full_procedure_method_a` — 17 passed.
- `testing.test_e8_2_baseline_stage_audit` — 23 passed.
- `testing.test_e8_3_detached_method_a_reweighting_audit` — 21 passed.
- `testing.test_pion_hgcer_method_a_parent_preserving_correction` — 8 passed.
- `testing.test_pion_hgcer_method_a_tphi_propagation` — 8 passed.
- `testing.test_pion_hgcer_phase_f_runtime_contract` — 3 passed.
- `testing.test_pion_component_dynamic_alignment` — 2 passed and 10 skipped
  because PyROOT is unavailable on this host.

These are deterministic local checks only. They do not provide independent
source review or ROOT/PyROOT/farm/runtime validation.

## Next

`NEXT` — independent ChatGPT actual-diff/source-runtime-path review of this
F.6.3 candidate. Keep F.6.3 `ACTIVE` unless that review passes; do not begin
E.8.4 or run the Jefferson Lab farm as part of this review.

## Fix.1 — branch-safety and public-yield repair

Independent actual-diff review found five narrow candidate defects, not a
scientific-ownership change: an optional branch-local ROOT-like exception could
escape; live F.1 parity used a raw cache coefficient rather than the filler
coefficient; parity used an absolute rather than F.4 scaled tolerance; public
`Y0` and private `YA` had duplicated final-yield arithmetic; and retained
provenance derived a digest from event-level factors.

Fix.1 preserves the F.6.3 architecture and repairs only those defects:

- `_build_f6_3_parallel_method_a_source(...)` has one exception-only outer
  boundary around its private work. It catches ordinary `Exception` failures
  (not `BaseException`), returns a fresh unavailable source with exception
  class context, and retains no partial children. The completed public baseline
  remains outside that boundary.
- Live parity now invokes the existing private
  `_component_cache_event_coefficient(source_spec, index)` resolver used by the
  accepted template filler. It compares the effective applied `s_current`, the
  current `w0`, and `s_current * w0` to F.1; it does not modify a baseline
  coefficient to make an authority mismatch pass.
- The live cache gate uses the frozen F.4 scaled close rule exactly:
  `abs(left-right) <= 1e-12 * max(1, abs(left), abs(right))`.
- The public baseline loop and private Method-A branch both call the same
  producer-owned final yield/error helper. A real mocked
  `calculate_yield_data(...)` regression snapshots public groups, baseline
  yield/error, E.8.2 yield sidecar, stage-window yields, scale factors, and
  baseline histogram contents/errors under both unavailable and branch-runtime-
  failed F.6.3 conditions.
- Retained F.6.3 authority keeps only factor population count and identity
  fingerprint with accepted F.4 authority. The factor-value-derived population
  fingerprint was removed; event-level `C_j` remains transient only.

Codex ran the following deterministic local checks after Fix.1:

- `py_compile` for the F.6.3 helper, pion filler, yield producer, and focused
  test — PASS.
- `testing.test_f6_3_parallel_full_procedure_method_a` — 21 passed.
- `testing.test_e8_2_baseline_stage_audit` — 23 passed.
- `testing.test_e8_3_detached_method_a_reweighting_audit` — 21 passed.
- `testing.test_pion_hgcer_method_a_parent_preserving_correction` — 8 passed.
- `testing.test_pion_hgcer_method_a_tphi_propagation` — 8 passed.
- `testing.test_pion_hgcer_phase_f_runtime_contract` — 3 passed.
- `testing.test_pion_component_dynamic_alignment` — 2 passed and 10 skipped
  because PyROOT is unavailable on this host.

These are local deterministic checks only. Fix.1 is `ACTIVE`, is **NOT farm
validated**, and has no ROOT/PyROOT, full-analysis, or runtime acceptance
claim. `NEXT` remains independent ChatGPT review of the refreshed cumulative
F.6.3 diff; E.8.4 remains `BLOCKED` pending that review.

## Fix.2 — public-regression contract closure

Independent review found the Fix.1 source repairs correct. One deterministic
public-regression coverage gap remained: the existing real
`calculate_yield_data(...)` test reduced the public `groups` result to two
numbers and did not preserve/compare a processed component payload. Fix.2 is
test-only: it invokes the actual public `"yield"` mode; serializes and asserts
the complete one-child `groups` container/type, child-key inventory,
per-child field-key inventory, and values; and attaches a nested serializable
component-payload sentinel before each public call. Both the ordinary
unavailable branch and injected branch-local `RuntimeError` now prove that
sentinel is value/byte-equivalent to its pre-call snapshot and is the original
unreplaced object. The prior E.8.2 yield/stage-window/scale-factor/final-
histogram assertions and the private all-one Method-A regression remain in
place.

Codex ran the following deterministic local checks after Fix.2:

- `python -B -m py_compile testing/test_f6_3_parallel_full_procedure_method_a.py`
  — PASS.
- `testing.test_f6_3_parallel_full_procedure_method_a` — 21 passed.
- `testing.test_e8_2_baseline_stage_audit` — 23 passed.
- `testing.test_e8_3_detached_method_a_reweighting_audit` — 21 passed.
- `testing.test_pion_hgcer_method_a_parent_preserving_correction` — 8 passed.
- `testing.test_pion_hgcer_method_a_tphi_propagation` — 8 passed.
- `testing.test_pion_hgcer_phase_f_runtime_contract` — 3 passed.
- `testing.test_pion_component_dynamic_alignment` — 2 passed and 10 skipped
  because PyROOT is unavailable on this host.

These are deterministic local checks only. Fix.2 is **NOT farm validated** and
does not claim ROOT/PyROOT, full-analysis, or runtime acceptance. F.6.3 remains
`ACTIVE`; E.8.4 remains `BLOCKED`. `NEXT` is independent ChatGPT review of the
refreshed cumulative F.6.3 diff.

## Independent source-review closure

Independent ChatGPT actual-diff/source-runtime-path review of
`kaonlt_review(20260928-161516).diff` passed. The reviewed committed base
remained `c9ed0d6b4d0013f7475eedbfca57a096730ea840`. Fix.1 source
implementation passed independent review, and comparison with the prior
candidate found no substantive source change in Fix.2: it changed only the
focused test plus warranted memory, contract, and manifest state.

The Fix.2 public regression now verifies the actual `"yield"` mode; complete
one-child public `groups` structure and values; unavailable/runtime-failed
branch equality; E.8.2 final yield/statistical/total values; stage-window
yields; scale factor; baseline histogram contents/errors; and processed
component-payload object identity plus deep/serialized value preservation. The
private all-one Method-A baseline-reproduction regression remains covered.

The deterministic checks recorded above are Codex-reported and were **NOT RUN
by ChatGPT**. This establishes `SOURCE REVIEWED` only; no ROOT/PyROOT,
full-analysis, farm, or runtime acceptance is claimed. F.6.3 is `SOURCE
REVIEWED`; E.8.4 is `NEXT` and has not started.

`NEXT` — user-controlled commit/push of the reviewed cumulative F.6.3 set ->
ChatGPT pushed-state review -> E.8.4 source/runtime-path audit and standalone
implementation contract.
