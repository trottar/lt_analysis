# KaonLT — E.8.4 Fix.5.4 current-lineage numerical identity and SIMC-normalization audit/fix

## 1. Purpose

Implement the next active KaonLT gate after the pushed post-farm memory checkpoint.

This task addresses the numerical/provenance blockers exposed by the fresh
`Q4p4W2p74 / Left / lowe` Fix.5 procedure PDF **before** any visualization
restyling and before any new farm run.

The task has four scientific-integrity objectives:

1. prove the exact current-F.6.3 per-child algebra
   `MM_A - MM_0 = -(B_pi_A - B_pi_0)` bin by bin;
2. prove that the displayed current-branch MM histograms close to their stored
   scalar yields under the producer-owned integration semantics;
3. construct and validate a current-F.6.3 aggregate view from the same current
   candidate lineage rather than treating historical E.8.3/F.6.1 aggregates as
   the current branch;
4. trace and prove, or explicitly fail closed on, the absolute-normalization
   comparability of the E.8.4 per-child SIMC MM support.

This task may repair only deterministic source-linkage/provenance defects
exposed by that audit. It must not tune or redesign the physics.

No plotting-style/visibility improvement is in scope here. That is the next,
separately reviewed source-changing stage after this task is reviewed, pushed,
and synchronized.

No farm run is authorized by this task.

---

## 2. Exact starting state

Required branch:

```text
test
```

Required starting HEAD:

```text
71a912bf14c59eee6d52956462aa0aa39c7f7c84
```

Commit subject:

```text
E8.4 Fix.5.4: record post-farm identity audit checkpoint
```

Required relationship to the prior source:

```text
parent scientific source:
f9d70732290ea461096374ca1270b47452644991

71a912bf1 changes memory only relative to f9d707322.
The analysis/runtime source is therefore unchanged from the already inspected
Fix.5 source.
```

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short --untracked-files=all
```

Requirements:

- branch must be `test`;
- HEAD must equal the exact required starting HEAD;
- inspect and preserve the worktree;
- root `AGENTS.md` and `.codex/` remain local-only/untracked;
- unrelated tracked changes are a blocker;
- do not reset, stash, clean, commit, push, or run the farm.

---

## 3. Required startup reading

Read, in order:

1. root `AGENTS.md`;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read only the records/source needed for this task:

### Active memory/evidence

- `docs/memory/phases/e8-4-fix5-4-current-lineage-identity-audit.md`
- `docs/memory/investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md`
- `docs/memory/phases/e8-4-fix5-shareable-method-a-impact-pages.md`
- `docs/memory/evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md`
- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/LEARNINGS.md`
- `docs/memory/CODEX.md`

### Relevant source

- `src/binning/calculate_yield.py`
- `src/cuts/full_background_subtraction_plots.py`
- `src/cuts/pion_hgcer_method_a_parallel_full_procedure.py`
- only the directly called SIMC normalization/helper source required to trace
  `find_yield_simc`, `normfac_simc`, and
  `hist["_xsect_support_simc"]["mm"]`
- only the directly called pion-template helper source required to prove the
  existing `B_pi_0/B_pi_A` fill semantics.

### Relevant tests

- `testing/test_e8_4_production_impact_audit.py`
- `testing/test_full_background_subtraction_plots.py`
- relevant existing F.6.3 tests
- relevant yield/SIMC tests
- any focused new Fix.5.4 test file created by this contract.

Do not reopen unrelated phases.

---

## 4. Established evidence and source facts

Treat these as the starting evidence, not as hypotheses to rediscover.

### 4.1 Fresh farm-render observation

The fresh Fix.5 `Q4p4W2p74 / Left / lowe` procedure PDF rendered 97 pages with:

```text
renderer_failures = []
```

and contained all new Fix.5 E.8.4 pages.

That proves rendering occurred. It does **not** prove the numerical linkage
between every displayed histogram, scalar yield, lineage, and SIMC object.

### 4.2 Confirmed lineage mismatch

Historical E.8.3 consumes the accepted historical F.6.1 lineage.

The current F.6.3/E.8.4 branch consumes the separately validated
current-baseline candidate F.4 lineage.

Therefore:

```text
historical E.8.3 aggregate
!= automatically
current F.6.3 aggregate
```

Do not use historical E.8.3/F.6.1 objects as current-F.6.3 authority.

Historical E.8.3 remains valid for its own recorded historical scope.

### 4.3 Current producer path already observed

The current F.6.3 branch in `calculate_yield.py` constructs, per canonical
`(t,phi)` child:

```text
pion_input
B_pi_0
B_pi_A
MM_0
MM_A
Y0
YA
```

The current source constructs:

```text
MM_A = pion_input - B_pi_A
```

and stores `MM_0` from the baseline pion-subtraction path.

`YA` is measured from the constructed `MM_A`.

`Y0` is taken from the existing baseline yield measurement.

This creates exact identities that must be enforced explicitly rather than
assumed.

### 4.4 SIMC render path already observed

Fix.5 consumes the existing normalized per-child SIMC support from:

```python
hist["_xsect_support_simc"]["mm"]
```

after the existing SIMC-yield producer runs.

The renderer does not intentionally perform a new E.8.4 normalization.

That fact alone does not prove that the plotted SIMC MM object is in the same
absolute normalization and units as the compared data histograms. This task
must trace that producer path exactly.

---

## 5. Scientific ownership and frozen boundaries

### 5.1 This task owns

Only:

- numerical identity auditing of the current F.6.3 sidecar;
- exact current-lineage aggregate construction from already-produced current
  F.6.3 child objects;
- provenance/normalization auditing of existing SIMC yield support;
- narrow object-linkage/provenance repairs where existing source semantics
  establish the intended identity unambiguously;
- fail-closed validation of those identities in the E.8.4 consumer;
- deterministic tests and warranted memory updates.

### 5.2 Frozen

Do not change:

- random subtraction;
- dummy subtraction;
- slow-proton subtraction;
- pion-component fit definitions;
- pion-control weights `w0`;
- Method-A map/scientific calculation;
- candidate F.3/F.4 inputs or hashes;
- parent-preserving correction factors `C`;
- Method-B code or numerical role;
- cuts;
- MM windows;
- fit windows;
- priors;
- templates;
- binning;
- efficiency;
- acceptance;
- normalization formulas;
- SIMC physics/model normalization;
- yield formulas;
- cross-section formulas;
- L/T separation;
- production baseline objects;
- production output semantics.

No child may be renormalized independently.

No data or SIMC object may be rescaled to improve agreement.

---

## 6. Source audit before editing

Before changing source, trace and document in Codex's task notes the exact path:

```text
baseline data child
-> baseline pion-subtraction payload
-> final production H_MM_DATA
-> baseline scalar Y0

current F.6.3 factor map
-> current child B_pi_A
-> current child MM_A
-> current scalar YA

SIMC source/root input
-> existing SIMC normalization
-> per-child SIMC MM support
-> hist["_xsect_support_simc"]["mm"]
-> E.8.4 payload
-> renderer
```

Also establish exactly whether:

```text
payload["H_MM_after_pion_subtraction"]
```

and the histogram whose integral produced the stored baseline `Y0` are exactly
the same final object/content under the active `no_empirical_residual` profile.

Do not infer this merely because legacy empirical residual scales are zero.

If the source path proves a different final-baseline object is used for Y0,
classify that as an object-linkage defect before editing.

---

## 7. Required current-lineage identity audit

Implement one private, non-production audit owned by the current F.6.3 branch.

It must operate on detached/current branch objects and must not mutate
production histograms.

Use one narrowly named schema/version for the audit.

### 7.1 Per-child histogram algebra

For every canonical child that is current-F.6.3 `available`, validate identical
MM binning for:

```text
B_pi_0
B_pi_A
MM_0
MM_A
```

Then prove, for every normal MM bin:

```text
(MM_A - MM_0) + (B_pi_A - B_pi_0) = 0
```

within a strict floating-identity tolerance only.

This is an algebraic identity check, not a physics cut or fit criterion.

Use a single source-owned scaled floating tolerance of order `1e-12` unless an
already-established repository identity helper provides the same or stricter
semantics. Do not introduce a configurable user threshold.

Persist/report at minimum:

```text
maximum_absolute_bin_residual
maximum_scaled_bin_residual
failing_bin_count
algebra_identity_passed
```

Do not hide negative bins, zero bins, overflow behavior, or sparse children.

Define explicitly whether underflow/overflow are excluded or included, and use
the same choice consistently for the scalar-yield audit. Match the
producer-owned yield semantics.

### 7.2 Histogram-to-scalar yield closure

For every available child, prove using the exact producer-owned integration
semantics:

```text
Integral(MM_0) = Y0
Integral(MM_A) = YA
Integral(MM_A - MM_0) = YA - Y0
```

Persist/report:

```text
MM_0_integral
MM_A_integral
delta_MM_integral
Y0
YA
delta_Y
MM_0_minus_Y0
MM_A_minus_YA
delta_integral_minus_delta_Y
yield_identity_passed
```

The audit must not silently recompute and replace Y0 or YA merely to make this
pass.

### 7.3 Signed-support decomposition

For both `MM_0` and `MM_A`, report:

```text
positive_support
negative_support
signed_integral
absolute_support
```

with:

```text
signed_integral = positive_support + negative_support
absolute_support = positive_support - negative_support
```

when `negative_support` is stored as the signed negative contribution.

Also report the corresponding decomposition for:

```text
MM_A - MM_0
```

or sufficient equivalent information to explain signed-cancellation behavior.

This diagnostic must never become an uncertainty or correction.

### 7.4 Baseline histogram ownership identity

Explicitly audit the relationship between:

```text
the MM_0 histogram placed into the F.6.3/E.8.4 sidecar
```

and:

```text
the final baseline histogram whose yield was stored as Y0
```

Require equal geometry and content/error identity using the repository's
existing histogram fingerprint helper if applicable.

If they are not the same content:

- do not overwrite the scalar Y0;
- do not silently substitute one histogram for another;
- trace the existing downstream baseline operations;
- repair the sidecar linkage only if the source semantics establish which
  existing final-baseline histogram is the proper E.8.4 comparison object and
  the Method-A branch can be put through the same already-existing downstream
  operation without changing physics;
- otherwise fail closed with a precise reason and hard-stop this task before
  inventing new production logic.

---

## 8. Required current-lineage aggregate

The current branch must expose an aggregate made only from the same current
F.6.3 child objects.

For each canonical t parent, construct detached aggregate summaries from the
current children:

```text
B_pi_0_current(t, MM) = sum_phi B_pi_0(t, phi, MM)
B_pi_A_current(t, MM) = sum_phi B_pi_A(t, phi, MM)
Delta_B_pi_current     = B_pi_A_current - B_pi_0_current
```

No child normalization is allowed.

Record:

```text
child_inventory
populated_child_inventory
aggregate_B_pi_0_integral
aggregate_B_pi_A_integral
aggregate_delta_integral
```

and exact current candidate-lineage provenance/fingerprints.

Where the current candidate F.4 parent sums are in the same defined units,
validate aggregate-integral closure to those current candidate parent sums.

If the source semantics show that the persisted F.4 parent sum and the plotted
pion-template integral are not the same quantity/units, do not force the
comparison. Record the distinction explicitly.

The aggregate must be clearly marked:

```text
current_f6_3_candidate_lineage
```

Historical E.8.3/F.6.1 aggregate objects must not satisfy or substitute for
this current-lineage record.

This task prepares the validated current aggregate for the later visualization
stage; it does not need to redesign the historical E.8.3 pages.

---

## 9. Required SIMC normalization/provenance audit

Trace the existing SIMC producer rather than reasoning from the renderer.

The audit must establish, from current source, the exact semantics of:

```python
hist["_xsect_support_simc"]["mm"][t][phi]
```

including at minimum:

- source SIMC file/object identity;
- setting identity;
- canonical t/phi geometry;
- normalization producer/function;
- normalization factor(s) actually applied;
- trial/luminosity/model-weight semantics as represented in existing source;
- resulting histogram units/meaning;
- whether those units are directly comparable in absolute amplitude to the
  E.8.4 `MM_0/MM_A` normalized-yield histograms.

### 9.1 If absolute comparability is proven

Persist explicit provenance/contract metadata sufficient for E.8.4 to validate
that fact without recomputing or rescaling SIMC.

Also record per child:

```text
SIMC_integral_in_the_same_lambda_window
```

using existing histogram contents only.

Do not add a new display scale.

### 9.2 If absolute comparability is NOT proven

Do not invent a conversion or rescale.

Instead:

- mark absolute SIMC comparison as unavailable with an explicit reason;
- preserve the existing SIMC yield producer unchanged;
- prevent E.8.4 from labeling the corresponding absolute overlay as a validated
  authoritative amplitude comparison;
- record the source-level blocker in CURRENT/phase/investigation memory;
- hard-stop any attempted physics conclusion about data/SIMC amplitude.

A source-proven lack of comparable units is an acceptable scientific result of
this audit.

---

## 10. E.8.4 consumer requirements

`full_background_subtraction_plots.py` may validate and clone the new audit
records, but it must not recompute the scientific branch.

The E.8.4 payload must carry sufficient detached structured information for the
later visualization stage to consume:

```text
current-lineage child identity audit
current-lineage aggregate audit
signed-support diagnostics
SIMC normalization/provenance audit
```

Fail closed on malformed, missing, stale, cross-lineage, or inconsistent audit
records.

Do not silently fall back to historical E.8.3/F.6.1 aggregates.

Do not silently fall back to setting-wide SIMC shapes.

Do not create a new renderer-side normalization.

### Lineage wording repair

Any E.8.4 explanatory text that implies historical E.8.3 is the aggregate
source for the current F.6.3 candidate branch must be corrected.

The repair must preserve historical E.8.3's own identity and status rather than
rewriting it as current-lineage science.

No line/marker/draw-order visibility changes are allowed in this task.

---

## 11. Required behavior on discovered defects

This task is an audit plus narrow repair, not an open-ended redesign.

Allowed repair classes:

1. **wrong sidecar object linkage** where current source unambiguously already
   owns the correct object;
2. **missing identity/provenance validation**;
3. **cross-lineage presentation/provenance linkage**;
4. **missing detached audit summary**;
5. **missing SIMC normalization metadata** that can be copied from existing
   producer semantics without changing the normalization.

Forbidden repair classes:

- changing Method-A `C` values;
- changing baseline pion weights;
- changing pion templates to make closure pass;
- changing SIMC scale;
- changing data normalization;
- altering yield arithmetic;
- changing MM integration window;
- clipping negative content;
- changing binning;
- smoothing/interpolating;
- introducing child normalization;
- changing physics to improve agreement.

If a required identity fails for reasons outside the allowed repair classes,
report the exact blocker and stop rather than broadening scope.

---

## 12. Allowed substantive files

Expected source scope:

```text
src/binning/calculate_yield.py
src/cuts/full_background_subtraction_plots.py
```

Focused tests:

```text
testing/test_e8_4_production_impact_audit.py
testing/test_full_background_subtraction_plots.py
```

One new focused test file is allowed if useful, preferably:

```text
testing/test_e8_4_fix5_4_identity_audit.py
```

A directly called SIMC helper source file may be modified **only** if the source
audit proves that metadata required to describe existing normalization cannot be
attached correctly from `calculate_yield.py` without duplicating or guessing
producer semantics. If that occurs:

- record the exact path and reason before editing it;
- do not change its numerical normalization behavior;
- include it in the final diff inventory.

Do not modify other scientific source without a concrete blocker and a new
contract.

---

## 13. Frozen substantive files

Unless the narrow metadata exception above applies, freeze:

```text
src/main.py
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
src/cuts/pion_hgcer_method_a_parent_preserving_correction.py
src/cuts/pion_hgcer_method_a_tphi_propagation.py
src/cuts/pion_component_subtraction.py
src/cuts/pion_component_fits.py
src/cuts/background_config.py
run_Prod_Analysis.sh
testing/run_e8_4_fix5_left_lowe_plot_gate.py
testing/pion_hgcer_validation_bundle_profile_e8_4_fix5.json
```

Freeze all accepted F.3/F.4/F.5/F.6.1/F.6.2 artifact identities.

Freeze owner/collector/profile behavior in this task.

The known owner post-render packaging failure is a separate operational issue
and must not be “fixed” inside this scientific identity task.

---

## 14. Required tests

### 14.1 Positive identity tests

Construct deterministic fake histograms/records proving:

- exact bin-by-bin pion/MM algebra passes;
- `Integral(MM_0) == Y0`;
- `Integral(MM_A) == YA`;
- delta integral equals delta scalar;
- signed positive/negative/absolute support is computed correctly;
- baseline sidecar histogram identity matches the baseline-yield histogram;
- current aggregate equals the sum of current children;
- historical-lineage aggregate is not used;
- comparable SIMC provenance is accepted without rescaling.

### 14.2 Negative tests

Require explicit failure/rejection for:

- one-bin `MM_A` perturbation;
- one-bin `B_pi_A` perturbation;
- Y0 mismatch;
- YA mismatch;
- MM binning mismatch;
- baseline final-histogram fingerprint mismatch;
- duplicated/missing current child;
- stale/cross-lineage aggregate provenance;
- historical F.6.1 aggregate substituted for current F.6.3 aggregate;
- missing SIMC normalization provenance;
- ambiguous/non-comparable SIMC units;
- any attempted display normalization fallback.

### 14.3 Regression tests

Prove that, for already-consistent synthetic inputs:

- stored Y0 is unchanged;
- stored YA is unchanged;
- `MM_0` contents are unchanged;
- `MM_A` contents are unchanged;
- `B_pi_0/B_pi_A` contents are unchanged;
- Method-A factors are unchanged;
- production objects are not mutated;
- Method B remains numerically absent;
- no child renormalization occurs;
- existing E.8.4 payload/page IDs remain available;
- no visualization style changes occur.

Run the relevant pre-existing F.6.3/E.8.4/yield tests as regressions.

---

## 15. Local validation

Use the repository-selected Python interpreter.

At minimum run syntax checks for every changed/new Python file and focused unit
tests covering:

```text
testing/test_e8_4_fix5_4_identity_audit.py
testing/test_e8_4_production_impact_audit.py
testing/test_full_background_subtraction_plots.py
```

plus relevant F.6.3/yield/SIMC regression tests identified by the source audit.

Do not claim ROOT/PyROOT/full-analysis behavior from local tests.

If PyROOT-dependent tests are unavailable locally, preserve explicit skips and
report them.

Also run:

```bash
git diff --check
```

and the ordinary memory-health command after warranted memory updates:

```bash
<PYTHON> -B tools/check_memory_health.py --root .
```

---

## 16. Warranted memory updates

Update memory in the same scoped change to reflect what the source audit
actually established.

Allowed memory files:

```text
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/LEARNINGS.md
docs/memory/roadmap/STATUS.md
docs/memory/phases/e8-4-fix5-4-current-lineage-identity-audit.md
docs/memory/investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md
docs/memory/manifest.json
```

Update only warranted durable knowledge.

Required status semantics after successful local implementation:

```text
E.8.4 Fix.5.4:
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

only if deterministic source/tests close the intended source-level work.

If the source audit instead exposes a blocker that cannot be repaired within
scope, keep Fix.5.4:

```text
ACTIVE
```

or:

```text
BLOCKED
```

as warranted, and record the one exact next repair.

Do not mark runtime closure.

`CURRENT.md` must preserve the next dependency:

```text
after Fix.5.4 actual-diff review -> user commit/push -> pushed-state review,
the next source-changing stage is visualization-only improvement.
```

Do not authorize farm execution yet.

Regenerate `docs/memory/manifest.json` whenever versioned memory changes.

---

## 17. Diff audit

Before stopping, report:

```bash
git status --short --untracked-files=all
git diff --stat
git diff --check
```

Audit the complete actual diff.

If the diff is too large for terminal review, create exactly one temporary
review file in the repository root:

```text
kaonlt_review.diff
```

It must contain:

- the complete tracked diff for every changed tracked file;
- complete `git diff --no-index /dev/null ...` additions for every new file.

Do not stage merely to make the diff reviewable.

The review bundle is temporary and must not be committed.

---

## 18. Acceptance criteria

PASS only if all of the following hold:

1. exact starting HEAD is respected;
2. the real source path was traced before editing;
3. current F.6.3 and historical F.6.1 lineages remain distinct;
4. per-bin pion/MM algebra is explicitly validated;
5. histogram/scalar yield closure is explicitly validated;
6. signed-support decomposition is available;
7. current-lineage aggregates are constructed only from current F.6.3 children;
8. SIMC absolute comparability is either source-proven with provenance or
   explicitly unavailable without rescaling;
9. no physics/normalization/binning/cut/factor changes are made merely to pass;
10. production objects remain unmutated;
11. Method B remains numerically absent;
12. deterministic local tests pass, apart from explicitly documented
    environment-only skips;
13. warranted memory is updated and manifest regenerated;
14. `git diff --check` passes;
15. no visualization-style work is mixed into this task;
16. no farm run, commit, or push is performed.

---

## 19. Farm boundary

This task does not establish:

- ROOT/PyROOT integration;
- full `main.py` runtime behavior;
- actual farm histogram closure;
- rendered numerical-audit behavior;
- SIMC/data agreement;
- Fix.5.4 runtime closure;
- Method-A promotion.

Those remain farm-only.

Do not provide or execute a farm command in this task.

---

## 20. Hard stop

After implementation, deterministic local checks, warranted memory updates, and
complete diff preparation:

**STOP.**

Do not implement the later visualization-style changes.

Do not modify the farm owner/profile/collector for the packaging issue.

Do not commit.

Do not push.

Do not run the farm.

Return the actual `kaonlt_review.diff` to ChatGPT for independent review.
