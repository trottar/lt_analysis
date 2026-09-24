# KaonLT E.8.3 Fix.1 — Gap-Safe Ratio Rendering and Contract-Coverage Repair

## 1. Objective

Repair the current local E.8.3 detached Method-A reweighting-audit candidate after
independent ChatGPT actual-diff review.

The E.8.3 scientific/runtime architecture is **accepted in principle and must not
be redesigned**:

- E.8.3 remains detached and presentation-only.
- It consumes only accepted persisted F.4/F.5/F.6.1/F.6.2 artifacts.
- F.4 owns the accepted parent-preserving correction.
- F.5 owns the accepted canonical `3 x 9` baseline/Method-A signed aggregates and
  parent closure.
- F.6.1 owns the accepted parent-level signed missing-mass baseline/Method-A/delta
  arrays.
- the existing frozen F.6.2 E.8 reader remains the accepted explanatory authority.
- Method B remains numerically absent.
- empirical residual Fit 1/Fit 2 remain absent.
- F.6.3 remains the only future owner of the parallel Method-A production branch.

Fix.1 owns only three defects/gaps found in source review:

1. make the ratio display **gap-safe** so undefined denominator bins cannot be
   visually connected/interpolated by a ROOT line;
2. preserve the persisted F.5 empty-child identity explicitly rather than
   validating `event_counts` and then discarding it;
3. complete the deterministic tests required by the original E.8.3 contract and
   record the exact Codex-run commands/results in the E.8.3 phase record.

No farm run is authorized.

---

## 2. Exact starting identity and dirty-worktree rule

Committed branch:

```text
test
```

Committed HEAD must remain exactly:

```text
91bb7809d27d84d7709a6600dbb0dc9ab514a458
```

Remote `test` was independently checked at that exact commit during Fix.1
contract preparation.

The worktree is intentionally dirty with the already reviewed/accepted
post-E.8.2 memory reconciliation plus the current local E.8.3 candidate.

**Do not reset, stash, clean, checkout away, discard, or rewrite that cumulative
work.**

Before editing, run:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
git diff --name-only
```

Hard stop if:

- branch is not `test`;
- committed HEAD is not exactly the SHA above;
- unrelated changes are present beyond the already accepted cumulative
  reconciliation/E.8.3 candidate and this Fix.1 contract.

Do not pull, commit, push, update remote refs, or run the Jefferson Lab farm.

---

## 3. Required startup reading

Read the normal repository-memory startup sequence in order:

1. root `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read:

- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/phases/phase-f6-method-a-production-promotion.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/e8-3-detached-method-a-reweighting-audit-task-contract.md`
- `docs/memory/phases/e8-3-detached-method-a-reweighting-audit.md`
- this Fix.1 contract

Inspect the current cumulative candidate in:

- `src/cuts/full_background_subtraction_plots.py`
- `src/cuts/rand_sub.py`
- `testing/test_e8_3_detached_method_a_reweighting_audit.py`

Do not restart the E.8.3 design.

---

## 4. Independent-review findings that Fix.1 must repair

### 4.1 Undefined ratio bins are masked numerically but visually bridged

The current helper correctly omits bins whose persisted baseline denominator is
zero or nonfinite.

However the current ROOT renderer builds one `TGraph` from all surviving points
and draws it with:

```python
graph.Draw("ALP")
```

The `L` option connects adjacent surviving graph points. If an undefined ratio
bin lies between them, the rendered line bridges that missing bin and visually
implies a ratio through a location where `B_pi^0 == 0` or is otherwise
undefined.

That conflicts with the E.8.3 contract:

```text
ratio only where the persisted baseline denominator is finite and nonzero
```

and:

```text
Do not ... interpolate ... or fabricate undefined ratio bins.
```

### Required repair

Make the ratio representation gap-safe.

The narrow preferred repair is a **points-only** ratio display:

```text
AP
```

or an equivalently explicit representation that cannot connect across undefined
bins.

Do not:

- fill undefined bins with zero;
- assign sentinel ratios;
- connect across missing bins;
- interpolate;
- smooth;
- cap or clip.

Add a deterministic test that proves the renderer/draw policy cannot introduce a
line across masked ratio bins.

---

### 4.2 Persisted F.5 empty-child identity is validated and then discarded

The current F.5 validator validates the persisted `event_counts` `3 x 9`
matrix, but does not copy it into the selected E.8.3 payload.

As a result, the renderer shows all nine numeric child values but cannot
distinguish:

```text
event_count == 0  -> accepted explicit empty child
```

from:

```text
event_count > 0 and signed aggregate happens to be zero
```

The original contract explicitly requires accepted empty children to remain
explicit rather than merely surviving as an unlabeled numerical zero.

### Required repair

Pass through the accepted persisted F.5 `event_counts` matrix without
recalculation or reinterpretation.

For the selected parent, retain the nine persisted child event counts and use
them only for presentation identity/status.

For each child row:

- if persisted `event_count == 0`, label it explicitly as `EMPTY` (or an
  equivalently unambiguous presentation label);
- otherwise display it as populated, without changing any signed aggregate.

Do not infer emptiness from `B_pi^0 == 0`, `B_pi^A == 0`, or delta values.

Add tests proving:

- a true persisted empty child stays explicitly empty;
- a populated child with zero signed content is **not** relabeled as empty;
- all nine canonical children remain present;
- no child renormalization occurs.

---

### 4.3 Required deterministic contract coverage is incomplete

The original E.8.3 contract says **“At minimum cover”** the listed reader,
geometry, presentation, and forbidden-call regression cases.

The current focused test module contains only eight test methods and does not
demonstrate several required negative/regression boundaries.

Fix.1 must add explicit deterministic coverage for the missing boundaries listed
below.

---

## 5. Allowed files

Fix.1 may edit only:

```text
src/cuts/full_background_subtraction_plots.py
testing/test_e8_3_detached_method_a_reweighting_audit.py
docs/memory/phases/e8-3-detached-method-a-reweighting-audit.md
docs/memory/manifest.json
```

and add this contract:

```text
docs/memory/phases/e8-3-detached-method-a-reweighting-audit-fix1-task-contract.md
```

No other file may change in Fix.1.

The complete review bundle will still contain the pre-existing cumulative
changes from the post-E.8.2 reconciliation and initial E.8.3 candidate. Fix.1
must not rewrite those unrelated patch sections.

---

## 6. Frozen files

From the current E.8.3 candidate onward, keep these byte-for-byte unchanged in
Fix.1:

```text
src/cuts/rand_sub.py
src/binning/calculate_yield.py
src/main.py
src/cuts/pion_component_subtraction.py
src/cuts/pion_component_fits.py
src/cuts/particle_subtraction.py
src/cuts/pion_hgcer_method_a_parent_preserving_correction.py
src/cuts/pion_hgcer_method_a_tphi_propagation.py
src/cuts/pion_hgcer_method_a_reweighting_validation.py
src/cuts/pion_hgcer_method_a_acceptance_refinement_validation.py
testing/analyze_pion_hgcer_method_a_parent_preserving_correction.py
testing/analyze_pion_hgcer_method_a_tphi_propagation.py
testing/analyze_pion_hgcer_method_a_reweighting_validation.py
testing/analyze_pion_hgcer_method_a_acceptance_refinement_validation.py
```

Also do not alter the already accepted cumulative post-E.8.2 reconciliation
records except for automatic `manifest.json` updates required by this Fix.1.

No changes to profiles, collectors, launchers, farm scripts, F.6.2 artifacts,
cuts, pion weights, proton factors, normalization, binning, templates, fits,
SIMC, yields, efficiencies, acceptance, L/T separation, or cross sections.

If any frozen file appears necessary, hard stop and report the blocker.

---

## 7. Required source behavior after Fix.1

### 7.1 Keep the existing authority path

Preserve the current byte-pinned reader sequence exactly in scientific meaning:

```text
accepted F.4 JSON
accepted F.5 JSON
accepted F.6.1 JSON
accepted F.6.2 payload from the existing frozen E.8 reader
    -> fail-closed E.8.3 presentation payload
    -> existing procedure-PDF lifecycle
```

Do not add any alternate file discovery, fallback, search, recomputation, or
producer.

### 7.2 Ratio

Keep the existing rule:

```text
ratio defined iff persisted B_pi^0 is finite and nonzero
```

A finite Method-A numerator is also required.

The renderer must make missing ratio bins visibly absent and must not connect
across them.

### 7.3 F.5 child identity

Retain, for every selected canonical t parent:

```text
f5_event_counts[0:9]
f5_baseline_signed_contents[0:9]
f5_method_a_signed_contents[0:9]
f5_signed_delta_contents[0:9]
```

all copied directly from persisted F.5.

The E.8.3 renderer may use `event_count == 0` only to print/display the accepted
empty-cell status. It must not alter any numerical aggregate.

### 7.4 Parent closure

Continue displaying the accepted persisted F.5 parent closure.

Do not invent a relative residual if the accepted F.5 schema does not persist
one.

---

## 8. Required focused tests

Retain the useful existing tests and add explicit coverage so the original
contract is actually demonstrated.

At minimum the focused E.8.3 test module must prove all of the following.

### Authority / fail-closed tests

1. accepted/synthetic-authorized F.4 input passes;
2. wrong F.4 raw SHA fails closed;
3. wrong F.4 correction fingerprint fails closed;
4. wrong F.4 artifact fingerprint fails closed;
5. accepted/synthetic-authorized F.5 input passes;
6. wrong F.5 raw SHA fails closed;
7. wrong F.5 propagation fingerprint fails closed;
8. wrong F.5 artifact fingerprint fails closed;
9. accepted/synthetic-authorized F.6.1 input passes;
10. wrong F.6.1 raw SHA fails closed;
11. wrong F.6.1 artifact fingerprint fails closed;
12. wrong F.6.1 validation fingerprint fails closed;
13. nonfinite JSON fails closed;
14. wrong kinematic token fails closed.

### Canonical inventory / geometry tests

15. missing canonical F.5 setting fails closed;
16. duplicate canonical F.5 setting fails closed;
17. missing canonical F.4 parent fails closed;
18. duplicate canonical F.4 parent fails closed;
19. missing canonical F.6.1 parent fails closed;
20. duplicate canonical F.6.1 parent fails closed;
21. malformed F.5 `3 x 9` matrix fails closed;
22. malformed F.6.1 MM edges fail closed;
23. malformed F.6.1 MM content length fails closed;
24. exact current-setting selection returns only its three t parents.

### Persisted-value / presentation tests

25. F.6.1 baseline/Method-A/delta arrays pass through unchanged;
26. F.5 baseline/Method-A/delta child aggregates pass through unchanged;
27. F.5 persisted event counts pass through unchanged;
28. accepted empty child remains explicitly `EMPTY`;
29. populated zero-valued child is not called empty;
30. all nine phi children remain present;
31. no child renormalization occurs;
32. persisted parent closure values are the ones rendered/displayed;
33. ratio helper includes only finite nonzero-denominator bins;
34. ratio renderer is gap-safe and does not use a connecting line across masked
    bins;
35. no clipping, capping, smoothing, or interpolation helper is introduced.

Use small fakes/mocks for ROOT-facing renderer tests. Do not require PyROOT for
these deterministic checks.

### Runtime/forbidden-call tests

36. existing `OUTPATH` route remains the only F.4/F.5/F.6.1 runtime input path;
37. E.8.3 remains in the existing E.8.2 pair-safe PDF/manifest lifecycle;
38. unavailable E.8.3 is presentation-only and nonfatal to baseline production;
39. prior E.8/E.8.1/E.8.2 page ordering remains unchanged except the intended
    E.8.3 insertion;
40. E.8.3 source/runtime path has no Method-B numerical builder call;
41. E.8.3 source/runtime path has no empirical residual Fit 1/Fit 2 builder call;
42. E.8.3 source/runtime path has no Method-A production application;
43. E.8.3 source/runtime path has no yield recalculation;
44. E.8.3 payload/rendering does not mutate supplied baseline production/yield
    sentinel objects.

Prefer executable mocks/spies for forbidden callable behavior where practical.
A source-text assertion may supplement but should not be the sole proof of a
callable boundary when an executable test is straightforward.

Do not duplicate scientific calculations in the tests.

---

## 9. Existing regression suites

Use the repository-appropriate interpreter discovered from repository
instructions/environment.

Run at minimum:

```bash
<PYTHON> -B -m py_compile \
  src/cuts/full_background_subtraction_plots.py \
  src/cuts/rand_sub.py \
  testing/test_e8_3_detached_method_a_reweighting_audit.py

<PYTHON> -B -m unittest testing.test_e8_3_detached_method_a_reweighting_audit -v
<PYTHON> -B -m unittest testing.test_full_background_subtraction_plots -v
<PYTHON> -B -m unittest testing.test_pion_hgcer_phase_e_runtime_contract -v

<PYTHON> -B -m unittest testing.test_pion_hgcer_method_a_parent_preserving_correction -v
<PYTHON> -B -m unittest testing.test_pion_hgcer_method_a_tphi_propagation -v
<PYTHON> -B -m unittest testing.test_pion_hgcer_method_a_reweighting_validation -v
<PYTHON> -B -m unittest testing.test_pion_hgcer_method_a_acceptance_refinement_validation -v
```

If an exact module name differs in the repository, use the actual existing
module and report the exact command.

Run the required memory maintenance/checks and:

```bash
git diff --check
```

No ROOT/PyROOT/full `main.py`/farm execution is required or authorized.

---

## 10. Durable phase record

Update:

```text
docs/memory/phases/e8-3-detached-method-a-reweighting-audit.md
```

Keep E.8.3:

```text
ACTIVE
```

pending independent ChatGPT Fix.1 actual-diff review.

Add a narrow Fix.1 chronology describing:

- the gap-safe ratio repair;
- persisted `event_counts` empty-child identity;
- completed focused negative/regression coverage;
- exact local interpreter used;
- every test/check command actually run;
- exact PASS/FAIL/SKIP counts or result for each command;
- any pre-existing unrelated failure, if present, with evidence it is unchanged.

State explicitly:

```text
Codex-reported checks; NOT RUN by ChatGPT.
```

Do not mark E.8.3 `SOURCE REVIEWED`.

Do not claim ROOT/PyROOT, full-analysis, farm, or runtime validation.

Regenerate:

```text
docs/memory/manifest.json
```

after versioned memory changes.

---

## 11. Diff audit

Before stopping:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
git diff --name-only
git diff --stat
git diff --check
```

Fix.1 incremental changes must be limited to the allowlist in Section 5.

Explicitly confirm that `src/cuts/rand_sub.py` is unchanged from the current
pre-Fix.1 E.8.3 candidate.

Audit for absence of changes to:

- production weights;
- pion subtraction;
- proton subtraction;
- random/dummy treatment;
- canonical binning;
- public yield/error calculation;
- Method B;
- empirical residual fits;
- F.6.3 production construction.

---

## 12. Complete cumulative review bundle

Refresh root:

```text
kaonlt_review.diff
```

It must be a **complete cumulative diff from committed HEAD
`91bb7809d27d84d7709a6600dbb0dc9ab514a458`**, not merely the Fix.1 delta.

Include:

1. complete tracked `git diff`;
2. complete `git diff --no-index /dev/null ...` content for every intended
   untracked/new file, including:
   - the accepted post-E.8.2 reconciliation contract;
   - the original E.8.3 contract;
   - this Fix.1 contract;
   - the E.8.3 phase record;
   - the focused E.8.3 test file.

Do not add `kaonlt_review.diff` to the intended commit.

---

## 13. Acceptance criteria

Fix.1 is ready for independent review only if:

- committed HEAD remains exact;
- no existing cumulative accepted work was discarded;
- ratio rendering cannot connect across undefined bins;
- no undefined ratio value is fabricated;
- F.5 event-count/empty identity is preserved and displayed;
- a populated zero-valued child is distinguishable from a true empty child;
- all original contract-required negative/regression boundaries are explicitly
  tested;
- exact commands/results are recorded in the phase record;
- no frozen source file changes;
- no production/scientific behavior changes;
- E.8.3 remains `ACTIVE`;
- F.6.3 remains `BLOCKED`;
- memory manifest/checks pass;
- `git diff --check` passes;
- complete cumulative `kaonlt_review.diff` is refreshed;
- Codex does not commit, push, update remote refs, or run the farm.

---

## 14. Hard stop

Stop and report instead of broadening scope if:

- a gap-safe ratio requires changing accepted scientific arrays;
- explicit empty-child status cannot be obtained directly from persisted F.5
  `event_counts`;
- required executable test coverage would require editing frozen production
  source;
- any accepted upstream artifact/schema is inconsistent with the current
  E.8.3 assumptions;
- the cumulative dirty worktree cannot be preserved;
- starting committed HEAD is no longer exact.

Do not redesign E.8.3 or begin F.6.3 inside Fix.1.
