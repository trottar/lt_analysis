# KaonLT F.6.3 Fix.1 — Branch Fail-Closed Safety, Live Baseline Parity, and Public Yield Regression

## Purpose

Repair the current local F.6.3 candidate only. Do not redesign F.6.3, begin
E.8.4, promote Method A, or change production physics.

The current F.6.3 architecture is retained:

```text
authoritative baseline:
    b_j^0 = s_j * w0_j

private parallel Method-A branch:
    b_j^A = s_j * w0_j * C_j
```

The independent actual-diff review found that the candidate is scientifically
well scoped but is not yet acceptable because the optional parallel branch does
not fully satisfy the required fail-closed/public-baseline contract and its
live baseline parity check is not guaranteed to use the exact coefficient that
the production filler actually applies.

This is a **narrow repair**.

---

## Exact committed base and current worktree

Authoritative committed base:

```text
branch: test
HEAD: c9ed0d6b4d0013f7475eedbfca57a096730ea840
subject: Add E8.3 detached Method-A reweighting audit
```

The worktree already contains the uncommitted cumulative F.6.3 candidate
reviewed in:

```text
kaonlt_review(20260924-185653).diff
```

Do not reset, stash, clean, checkout over, or discard that candidate.

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
git log -1 --oneline
```

Hard stop if committed HEAD differs from
`c9ed0d6b4d0013f7475eedbfca57a096730ea840` or if unrelated dirty work is
present.

---

## Required startup reading

Read, in order:

1. local root `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read:

- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/phase-f6-method-a-production-promotion.md`
- `docs/memory/phases/e8-2-baseline-stage-audit.md`
- `docs/memory/phases/e8-3-detached-method-a-reweighting-audit.md`
- `docs/memory/phases/f6-3-parallel-full-procedure-method-a-task-contract.md`
- `docs/memory/phases/f6-3-parallel-full-procedure-method-a.md`
- this Fix.1 contract

Status remains:

```text
E.8   — ACTIVE
E.8.2 — SOURCE REVIEWED
E.8.3 — SOURCE REVIEWED
F.6.3 — ACTIVE
E.8.4 — BLOCKED
```

Do not mark F.6.3 `SOURCE REVIEWED`.

---

# Review findings that Fix.1 must repair

## Finding 1 — branch-local ROOT/runtime exceptions can escape the optional branch

The current builder creates/reset/clones/fills/adds ROOT-like objects inside:

```python
_build_f6_3_parallel_method_a_source(...)
```

but its outer failure boundary catches only:

```python
(MethodAParallelFullProcedureError, KeyError, TypeError, ValueError)
```

The existing repository ownership helper can itself raise `RuntimeError` when
a required ROOT clone fails. Other branch-local ROOT-like operations may also
raise Python exceptions.

F.6.3 is a non-promoted private side branch invoked only after the
authoritative baseline yield has already been constructed. Therefore an F.6.3
branch-local failure must never abort or invalidate the public baseline run.

### Required repair

Install one narrow outer exception boundary around **only**
`_build_f6_3_parallel_method_a_source(...)`'s optional F.6.3 work.

Requirements:

- baseline calculation remains outside this boundary;
- `KeyboardInterrupt`, `SystemExit`, and other `BaseException` classes are not
  swallowed;
- ordinary branch-local `Exception` failures are converted into one fresh
  `unavailable_parallel_source(...)`;
- do not return a partially constructed Method-A source;
- preserve a precise F.6.3 failure reason, including enough exception class
  context to diagnose the failure;
- do not silently substitute `C_j = 1`;
- do not weaken accepted authority checks.

A broad `except Exception` is permitted **only at this detached F.6.3 outer
boundary**, because its sole purpose is to prevent optional side-branch
failures from changing the already-completed baseline result. Do not introduce
broad exception swallowing anywhere in baseline physics.

Add deterministic tests injecting at least:

- `clone_root_histogram -> RuntimeError`;
- Method-A template fill -> `RuntimeError`;
- Method-A histogram `Add` or equivalent branch-local operation ->
  `RuntimeError`.

Each must produce:

```text
F.6.3 source available = false
public baseline calculation/result unchanged
no partial Method-A children retained
```

---

## Finding 2 — live parity uses the raw cached coefficient, not necessarily the exact coefficient the filler applies

The current candidate builds live F.1 parity rows with:

```python
coefficient = float(section["coefficient"][index])
```

but the actual accepted pion-template filler computes the event coefficient via:

```python
_component_cache_event_coefficient(source_spec, idx)
```

and only then multiplies by `w0`.

F.6.3 must validate the **actual current baseline event contribution that the
filler uses**, not a nearby raw cache field whose equivalence is merely
assumed.

### Required repair

Use exactly the same coefficient resolver for both:

1. live parity construction; and
2. baseline/Method-A template filling.

Do not duplicate the coefficient formula.

The live row must be built from:

```text
s_current = exact _component_cache_event_coefficient(source_spec, idx)
w0_current = exact simc_shape_pion_weight_from_value(...)
b0_current = s_current * w0_current
```

and those values must be compared to accepted F.1:

```text
signed_source_coefficient
baseline_pion_weight_w0
signed_baseline_event_contribution
```

If current coefficient rescaling differs from accepted F.1, Method A is
unavailable for the full selected setting. Do not modify the current baseline
coefficient to force agreement.

A private import/reuse of the existing coefficient helper is acceptable if it
is the narrowest implementation. Do not rename or redesign the production
filler merely for aesthetics.

Add a regression in which cached coefficient and effective current
coefficient differ through the existing `coefficient/base_coefficient`
semantics. The validator must compare the **effective applied coefficient** and
fail if it no longer matches F.1 authority.

---

## Finding 3 — live parity tolerance is not the accepted F.4 tolerance semantics

The candidate currently uses:

```python
math.isclose(..., rel_tol=0.0, abs_tol=1.0e-12)
```

for live F.1 parity.

Accepted F.4 uses its frozen scaled closure rule:

```text
abs(left - right) <= tolerance * max(1, abs(left), abs(right))
```

with the accepted `1e-12` signed-baseline identity tolerance.

The new F.6.3 parity gate must not invent a different numerical acceptance
rule.

### Required repair

Use the exact accepted F.4 close semantics for:

- `analysis_MM`;
- `analysis_t`;
- effective signed source coefficient;
- `w0`;
- signed baseline contribution.

Do not loosen the tolerance beyond accepted F.4 behavior.

Prefer direct reuse of the existing F.4 close semantic/helper when practical;
otherwise a tiny validation-only helper with the exact formula is acceptable.
Do not create a new configurable tolerance.

Add tests at values where a pure absolute `1e-12` check and the accepted
scaled check differ, proving the accepted F.4 rule is used.

---

## Finding 4 — YA uses a new yield helper, but the public Y0 path still owns separate arithmetic

The candidate adds:

```python
_measure_yield_with_current_uncertainty(...)
```

and calls it for `YA`.

The authoritative public `Y0` path still executes the historical yield/error
arithmetic separately. That means baseline and Method A do not yet invoke one
shared implementation even though F.6.3 requires identical yield-extraction
mathematics.

### Required repair

Refactor only the common final-yield/error arithmetic so both:

```text
public baseline Y0
private Method-A YA
```

call the same helper.

The helper must preserve the historical public baseline behavior exactly,
including:

- `integral_with_stat_error(...)`;
- dummy integral treatment;
- data charge uncertainty;
- dummy charge uncertainty;
- Fit1/Fit2 fractional terms;
- `ZeroDivisionError` behavior;
- non-finite yield/error sentinel behavior;
- statistical error behavior.

Do not change the formulas, operation meanings, uncertainty ownership, Lambda
histogram/window geometry, or public return schema.

Under `no_empirical_residual`, the existing Fit1/Fit2 terms remain zero and
dormant. This refactor must not activate those fits.

### Mandatory public regression

Add an executable regression through the real public:

```python
calculate_yield_data(...)
```

using the existing detached test infrastructure/mocks.

Snapshot the historical baseline/public result before the helper refactor and
prove after the refactor that the following are unchanged:

- returned `groups` structure and values;
- public `Y0`;
- public total error;
- E.8.2 source and its final baseline yield/statistical/total values;
- `stage_window_yields`;
- scale factors;
- component payloads;
- baseline ROOT-like histogram contents/errors used by the fixture.

Also prove:

- an unavailable F.6.3 branch does not alter any of those outputs;
- an injected F.6.3 branch-local runtime failure does not alter any of those
  outputs;
- an all-one Method-A branch gives `YA == Y0`,
  `YA_statistical_error == Y0_statistical_error`, and
  `YA_total_error == Y0_total_error`.

Do not satisfy this only by testing the private builder.

---

## Finding 5 — retained provenance hashes event-level factor values

The candidate retains:

```python
transient_factor_population_fingerprint =
    sha256(source_label, entry_index, C_j ...)
```

The F.6.3 contract permits aggregate identity/count fingerprints but explicitly
keeps event-level `C_j` transient and non-persisted. The accepted F.4
correction/artifact fingerprints already own the correction-value authority.

### Required repair

Remove factor-value-derived provenance from the retained F.6.3 source.

Retain only what is needed for identity/population provenance, for example:

```text
transient_factor_population_count
transient_factor_identity_fingerprint
accepted F.4 source/correction/artifact fingerprints
live_cache_identity_count
live_cache_identity_fingerprint
```

Do not retain:

- event-level `C_j`;
- a list/table of event factors;
- a digest whose payload includes event-level `C_j`.

Add a regression proving that retained provenance contains no factor values or
factor-value-derived fingerprint field.

The transient in-memory identity-to-factor mapping remains allowed only for the
immediate Method-A fill and must be discarded afterward.

---

# Preserved accepted architecture

Everything below remains accepted and must not be redesigned:

- deterministic F.1/F.3/F.4 file loading;
- accepted F.4 source/fingerprint authority;
- shared F.4 transient correction reconstruction;
- exact selected-setting `(source_label, entry_index) -> C_j` pairing;
- all-or-nothing selected-setting Method-A branch;
- same authoritative post-proton/post-prune pion input for both branches;
- optional `method_a_event_multipliers` on the existing template filler;
- same multiplier for allcuts and nommcuts;
- no child normalization;
- no Method B;
- no empirical Fit1/Fit2;
- no second ROOT/tree traversal;
- no mutation of public baseline histograms;
- private `_f6_3_parallel_method_a_source`;
- E.8.4 remains a later consumer and is not implemented here.

Do not change F.4/F.5/F.6.1/F.6.2 source or accepted artifacts.

---

# Allowed files

Source:

```text
src/binning/calculate_yield.py
src/cuts/pion_component_subtraction.py
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
```

Tests:

```text
testing/test_f6_3_parallel_full_procedure_method_a.py
```

Memory/history:

```text
docs/memory/CURRENT.md
docs/memory/decisions/e8-full-analysis-procedure-roadmap.md
docs/memory/roadmap/STATUS.md
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/phases/f6-3-parallel-full-procedure-method-a.md
docs/memory/phases/f6-3-parallel-full-procedure-method-a-fix1-task-contract.md
docs/memory/manifest.json
```

The original F.6.3 task contract remains frozen history.

No other source/test file may change without a concrete blocker. Stop and
report rather than broadening scope.

---

# Frozen files / forbidden work

Do not edit:

```text
src/main.py
src/cuts/rand_sub.py
src/cuts/full_background_subtraction_plots.py
src/cuts/pion_hgcer_method_a_acceptance_contract.py
src/cuts/pion_hgcer_method_a_acceptance_map.py
src/cuts/pion_hgcer_method_a_parent_preserving_correction.py
src/cuts/pion_hgcer_method_a_tphi_propagation.py
src/cuts/pion_hgcer_method_a_reweighting_validation.py
src/cuts/pion_hgcer_method_a_acceptance_refinement_validation.py
src/utility/background_config.py
src/utility/root_histogram_ownership.py
```

Do not:

- change cuts, fits, windows, priors, amplitudes or normalizations;
- modify `w0`;
- modify F.4 `C_j`;
- use a factor fallback;
- persist event factors;
- normalize a child;
- activate Method B;
- activate empirical residual fits;
- change slow-proton handling;
- change pruning;
- change canonical binning;
- change SIMC;
- change cross-section/LT logic;
- implement E.8.4;
- commit/push;
- run the farm.

---

# Required focused tests after repair

The focused F.6.3 suite must explicitly cover at least:

1. branch `RuntimeError` at clone boundary -> unavailable, baseline survives;
2. branch `RuntimeError` at template-fill boundary -> unavailable, baseline
   survives;
3. branch `RuntimeError` at Method-A histogram arithmetic boundary ->
   unavailable, baseline survives;
4. no partial Method-A children survive any branch exception;
5. live parity uses exact effective `_component_cache_event_coefficient`;
6. cached/effective coefficient mismatch is detected;
7. accepted F.4 scaled tolerance semantics are used;
8. selected-setting identity inventory remains exact;
9. omitted multiplier preserves existing template output;
10. all-one multiplier preserves template output;
11. nontrivial multiplier changes only pion event weight;
12. allcuts and nommcuts use identical `C_j`;
13. missing/invalid multiplier fails;
14. public `calculate_yield_data(...)` baseline return is unchanged;
15. E.8.2/public baseline side products are unchanged;
16. unavailable F.6.3 cannot change public baseline;
17. branch exception cannot change public baseline;
18. all-one branch reproduces Y0/stat/total exactly;
19. no factor-value-derived provenance is retained;
20. helper remains ROOT-free and tree-traversal-free;
21. no Method-B / empirical-fit dependency is introduced.

Retain and update the existing focused tests rather than replacing them with
weaker assertions.

---

# Deterministic local validation

Run the same actual Python command/environment established by local
`AGENTS.md`.

At minimum:

```bash
<PYTHON> -m py_compile \
  src/cuts/pion_hgcer_method_a_parallel_full_procedure.py \
  src/cuts/pion_component_subtraction.py \
  src/binning/calculate_yield.py \
  testing/test_f6_3_parallel_full_procedure_method_a.py

<PYTHON> -m unittest testing.test_f6_3_parallel_full_procedure_method_a -v
<PYTHON> -m unittest testing.test_e8_2_baseline_stage_audit -v
<PYTHON> -m unittest testing.test_e8_3_detached_method_a_reweighting_audit -v
<PYTHON> -m unittest testing.test_pion_hgcer_method_a_parent_preserving_correction -v
<PYTHON> -m unittest testing.test_pion_hgcer_method_a_tphi_propagation -v
<PYTHON> -m unittest testing.test_pion_hgcer_phase_f_runtime_contract -v
<PYTHON> -m unittest testing.test_pion_component_dynamic_alignment -v

<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/check_memory_health.py --root .
<PYTHON> -B tools/memory_bootstrap.py --root . --json
<PYTHON> -B -m unittest testing.test_memory_health -v

git -c core.safecrlf=false diff --check
```

If the new public regression uses an existing E.8.2 support test module
directly, run that module as listed above. If another deterministic directly
affected module is evident, run it too.

Do not claim ROOT/PyROOT/full-analysis/farm validation.

---

# Memory/history update

Append Fix.1 to:

```text
docs/memory/phases/f6-3-parallel-full-procedure-method-a.md
```

Record:

- the independent review findings;
- exact repairs;
- exact tests Codex actually ran and results;
- explicit `NOT farm validated`;
- status remains `ACTIVE`;
- NEXT = independent ChatGPT review of the refreshed cumulative diff.

Update CURRENT/roadmap/status only as necessary to keep them internally
consistent. They must continue to say:

```text
E.8   — ACTIVE
E.8.2 — SOURCE REVIEWED
E.8.3 — SOURCE REVIEWED
F.6.3 — ACTIVE
E.8.4 — BLOCKED pending F.6.3 independent source review
```

Regenerate `docs/memory/manifest.json`.

Do not mark F.6.3 SOURCE REVIEWED.

---

# Diff audit

Before stopping:

```bash
git status --short
git diff --name-only
git diff --stat
git -c core.safecrlf=false diff --check
```

Every changed path must be within the cumulative F.6.3 allowlist plus this
Fix.1 contract.

Refresh root:

```text
kaonlt_review.diff
```

as the complete cumulative diff from committed HEAD:

```text
c9ed0d6b4d0013f7475eedbfca57a096730ea840
```

Include every intended untracked file with:

```bash
git diff --no-index /dev/null <path>
```

as required.

Do not stage or commit `kaonlt_review.diff`.

---

# Acceptance criteria for Fix.1

Return the candidate for independent review only when all are true:

1. committed HEAD remains exactly
   `c9ed0d6b4d0013f7475eedbfca57a096730ea840`;
2. branch-local ordinary exceptions cannot abort the already-completed public
   baseline;
3. no partial Method-A result survives branch failure;
4. live F.1 parity uses the exact effective coefficient applied by the template
   filler;
5. live signed baseline contribution is exactly that effective coefficient
   times current `w0`;
6. parity comparisons use accepted F.4 scaled tolerance semantics;
7. public Y0 and private YA use the same yield/error implementation;
8. public baseline outputs and E.8.2 side products are regression-identical;
9. all-one Method A reproduces Y0/stat/total exactly;
10. no event factor or factor-value-derived digest is retained;
11. no second tree traversal exists;
12. no child renormalization exists;
13. no Method B or empirical residual numerical path exists;
14. E.8.4 remains unimplemented and BLOCKED;
15. deterministic tests pass;
16. memory remains F.6.3 ACTIVE pending independent review;
17. no ROOT/PyROOT/farm/runtime claim is made.

---

# Hard stop

Stop and report a blocker instead of expanding scope if any repair requires:

- changing accepted baseline scientific behavior;
- changing F.4 mathematics or accepted artifacts;
- weakening F.1/F.4 authority;
- a second event traversal;
- changing pruning/proton subtraction;
- enabling empirical residual fits;
- Method B;
- modifying downstream cross-section/LT physics;
- implementing E.8.4;
- guessing farm-only behavior.

Do not commit, push, or run the farm.
