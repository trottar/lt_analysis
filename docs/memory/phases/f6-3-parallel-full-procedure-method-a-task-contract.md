# KaonLT F.6.3 — Parallel Full Procedure Plus Method A — Task Contract

## Status and objective

Starting status:

```text
E.8   — ACTIVE
E.8.2 — SOURCE REVIEWED
E.8.3 — SOURCE REVIEWED
F.6.3 — NEXT
E.8.4 — BLOCKED pending F.6.3
```

Implement **F.6.3 only**: construct a shadow/parallel full-analysis Method-A
branch beside the unchanged authoritative baseline branch.

The governing roadmap is
`docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`.

The scientific operation is exactly:

```text
baseline pion event contribution:
    b_j^0 = s_j * w0_j

parallel Method-A pion event contribution:
    b_j^A = s_j * w0_j * C_j
```

where `C_j` is exactly the accepted F.4 parent-preserving Method-A correction.

F.6.3 is the only current step allowed to construct this actual parallel
Method-A analysis branch. It does **not** promote Method A to production.
The public/baseline analysis remains authoritative and unchanged.

Do not implement E.8.4 presentation in this task.

---

## Exact starting identity

Branch:

```text
test
```

Committed HEAD:

```text
c9ed0d6b4d0013f7475eedbfca57a096730ea840
```

Commit subject:

```text
Add E8.3 detached Method-A reweighting audit
```

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
git log -1 --oneline
```

Hard stop if the branch or committed HEAD differs, or if unrelated dirty work
is present.

`AGENTS.md` and `.codex/` are local-only/untracked repository control files.
Read the local root `AGENTS.md`; do not require it to exist in GitHub history.

---

## Required startup reading

Read, in order:

1. root `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read:

- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/phase-f4-method-a-parent-preserving-correction.md`
- `docs/memory/phases/phase-f5-method-a-tphi-propagation.md`
- `docs/memory/phases/phase-f6-method-a-production-promotion.md`
- `docs/memory/phases/e8-2-baseline-stage-audit.md`
- `docs/memory/phases/e8-3-detached-method-a-reweighting-audit.md`
- this contract

Do not reopen or reinterpret closed F.1–F.6.2 science.

---

## Governing scientific boundary

Preserve all of the following:

- the existing public baseline yield branch;
- prompt/random algebra;
- dummy subtraction and normalization;
- accepted slow-proton treatment;
- existing pruning order and semantics;
- pion component definitions;
- pion-control fits, amplitudes, fit windows, templates and priors;
- accepted baseline missing-mass pion weight `w0`;
- canonical t/phi binning;
- SIMC;
- active background profile `no_empirical_residual`;
- efficiencies, acceptance, L/T separation and cross-section formulas;
- every F.1–F.6.2 accepted artifact and scientific fingerprint.

Method B is forbidden numerically.

Legacy empirical residual Fit 1 / Fit 2 remain dormant. Do not activate, rerun,
tune, compare, or repurpose them.

No canonical `(t,phi)` child may be independently normalized.

Do not clip, cap, winsorize, smooth, interpolate, replace, or otherwise alter
accepted F.4 `C_j`.

The baseline branch must remain the public return value and downstream
production input.

---

## Source audit that this implementation must preserve

### 1. Existing production pion-template application

`src/cuts/pion_component_subtraction.py` currently computes baseline pion
weights through:

```python
simc_shape_pion_weight_from_value(...)
```

and fills the accepted pion template in:

```python
fill_simc_shape_pion_subtraction_templates(...)
```

using the established cached signed source coefficient and baseline pion weight.
The existing event algebra must remain unchanged when no Method-A multiplier is
supplied.

The accepted component path in `src/binning/calculate_yield.py` calls that
template filler through `_apply_component_pion_subtraction_for_bin(...)` and
subtracts the resulting pion template from the current proton-cleaned/pruned
kaon histograms.

### 2. Existing authoritative event identity is already available

The existing current-analysis caches contain stable event identity fields. The
F.1 Method-A application contract uses:

```text
(source_label, entry_index)
```

as the stable application identity, with canonical t/phi assignment,
`signed_source_coefficient`, `baseline_pion_weight_w0`,
`signed_baseline_event_contribution`, `analysis_t`, and `analysis_MM`.

The current yield-side cache must be joined to accepted F.1/F.4 authority by
that identity. Do not add a second ROOT/tree traversal.

### 3. Accepted F.4 transient-factor reconstruction already exists

F.4 intentionally does not persist an event-factor table. Its public shared
calculator:

```python
build_pion_hgcer_method_a_parent_preserving_correction_with_review_data(...)
```

reconstructs event-level `C_j` transiently.

F.5 already proves the accepted pairing pattern:

1. validate accepted F.1/F.3/F.4 authority;
2. reproduce F.4 with the shared calculator;
3. require reproduced F.4 aggregate output/fingerprint to equal the accepted
   persisted F.4 artifact;
4. obtain the validated F.1 application rows;
5. require exact row/factor length alignment;
6. pair rows and factors only after those checks;
7. never persist the event-level factors.

F.6.3 must reuse this scientific calculation and authority boundary. Do not
reimplement the F.3/F.4 correction mathematics.

### 4. Correct branch point

The current production sequence is:

```text
prompt/random
-> dummy
-> slow-proton treatment
-> established prune_hist
-> accepted baseline pion template B_pi^0
-> baseline pion subtraction
-> no empirical residual correction under no_empirical_residual
-> final MM_0
-> existing yield extraction Y0
```

The Method-A shadow branch must start from the **same authoritative
post-proton/post-prune pion input** and differ only in the pion event
contribution `w0 -> w0*C`.

---

## Allowed implementation files

Source:

```text
src/cuts/pion_component_subtraction.py
src/binning/calculate_yield.py
```

Add one narrowly owned pure-Python authority/application helper:

```text
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
docs/memory/phases/f6-3-parallel-full-procedure-method-a-task-contract.md
docs/memory/manifest.json
```

Only update additional durable memory if an actually new durable fact requires
it. Do not perform general memory cleanup.

---

## Frozen implementation files

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
```

Do not modify F.1–F.6.2 accepted JSON/PDF artifacts.

E.8.4 rendering/finalization is explicitly out of scope here.

---

# Required implementation

## A. Pure-Python F.6.3 authority helper

Add:

```text
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
```

This module owns only:

1. accepted-artifact loading/validation needed for F.6.3;
2. transient F.4 factor reconstruction;
3. exact selected-setting identity-to-factor pairing;
4. live-cache parity validation;
5. aggregate provenance/fingerprints safe to retain after event factors are
   discarded.

It must not import ROOT or traverse analysis trees.

### Accepted authority

For `Q4p4W2p74`, consume the already accepted F.1/F.3/F.4 authority. Preserve
the frozen accepted F.4 identity:

```text
F.4 JSON SHA-256:
adcc01900b8adc211e296aaedaa314fdb86207d51483ebfd4b3edf69f558f188

F.4 correction fingerprint:
362241005c02f2149e260c391b5c3d35793287573128b42cf5ed693419d9d2f3

F.4 artifact fingerprint:
c4b9f513d5918ca77179bbe5fab28c5d9d14e63a440338545e73961fa50b9d67
```

Preserve inherited accepted F.3 authority:

```text
F.3 JSON SHA-256:
04a9a576767ba3f81c8c03daab2461998894c655f0d875bd532ad30cb46c0d95

F.3 map fingerprint:
81b2a1e89ef9689b24c6dd145f53b26f7ac8fee9d5666cbc2ae8caa9da2830e6

F.3 algorithm fingerprint:
ba29630b2f40a87cbadbe751504ce48f23e2b17a08378a2c8a219f93131cb912

F.3 artifact fingerprint:
f8a12313bbd81aca48402a8c7eb4773c2dbe3ad0c7c6f64a2f212e98202d6ee2
```

Use the existing source-owned accepted-runtime-authority constants and
validators wherever they already exist. Do not create a competing scientific
definition of acceptance.

Load deterministic artifact names using existing filename helpers where
available. Do not add filename guessing/fallback search.

### Transient factor construction

Invoke the existing F.4 shared calculator with review data. Require:

- accepted F.1 inventory and hashes;
- accepted F.3 authority;
- accepted F.4 file hash/fingerprints;
- reproduced F.4 aggregate object exactly matches the accepted persisted F.4
  correction;
- exact factor-row length equality;
- all `C_j` finite and strictly positive;
- selected setting is one of the canonical accepted settings.

Construct a transient mapping:

```python
(source_label, entry_index) -> C_j
```

for the selected setting only.

The mapping must never be written to JSON, stored on `hist`, returned as a
public production product, or included in durable diagnostics.

### Live current-cache parity

Before applying any factor, validate the live authoritative pion-control cache
against the accepted F.1 row for the same identity.

At minimum require parity for:

- source label;
- entry index;
- canonical t index;
- canonical phi index;
- current analysis missing mass versus accepted `analysis_MM`;
- current analysis t versus accepted `analysis_t`;
- current signed source coefficient versus accepted
  `signed_source_coefficient`;
- current baseline pion weight versus accepted `baseline_pion_weight_w0`;
- current signed baseline contribution versus accepted
  `signed_baseline_event_contribution`.

Use the accepted existing floating-point tolerance; do not invent a looser
physics tolerance.

Missing, extra, duplicate, non-finite, mismatched, or unsupported identities
make the Method-A branch unavailable/fail closed.

Do **not** substitute `C=1` for an authority failure.

The ordinary baseline analysis must still be allowed to continue because this
is a non-promoted parallel branch.

---

## B. Minimal extension of the existing pion-template filler

Extend:

```python
fill_simc_shape_pion_subtraction_templates(...)
```

with one optional Method-A event multiplier input.

Requirements:

- default/omitted multiplier preserves the existing baseline path exactly;
- current baseline `w0` is still computed by
  `simc_shape_pion_weight_from_value(...)`;
- current signed source coefficient remains exactly the existing cached
  coefficient;
- Method-A changes only:

```python
event_weight_baseline
```

to:

```python
event_weight_method_a = event_weight_baseline * C_j
```

for the validated identity;

- use the same multiplier in the `nommcuts` and `allcuts` template fills;
- never modify `pion_mm_weights`;
- never alter source coefficients;
- never alter proton factors already encoded by the accepted cache;
- never apply Method A to data, dummy, random, SIMC, proton-cleaning, or
  empirical-background weights.

No event-level factor may be retained after the fill.

A lookup miss while Method A is requested is an error for the Method-A branch;
it is not a neutral factor.

---

## C. Construct the parallel branch in `calculate_yield.py`

Do not duplicate the data tree traversal.

Use the already available authoritative current pion-control child cache.

Construct the Method-A branch only for the supported kaon/component/t-bin
path and only under the frozen:

```text
BG_OPT_ACTIVE_PROFILE = "no_empirical_residual"
```

For unsupported kinematics/settings/modes, retain an explicit unavailable
F.6.3 sidecar and leave baseline behavior untouched.

### Per-cell construction

For each valid canonical `(t,phi)` child:

1. retain the exact same authoritative post-proton/post-prune pion input used by
   the baseline branch;
2. retain the existing baseline pion template `B_pi^0`;
3. fill a separate Method-A pion template `B_pi^A` using the same accepted
   source cache, same `w0`, and only the additional accepted `C_j`;
4. construct:

```text
MM_0 = pion_input - B_pi^0
MM_A = pion_input - B_pi^A
```

on separate ROOT objects;

5. do not mutate the existing baseline `hist_bin_dict`, existing
   `component_subtraction_payloads`, existing E.8.2 captures, or public return
   values;
6. preserve identical Lambda integration/signal window geometry.

The Method-A branch must not call the empirical Fit 1/Fit 2 machinery.

### Yield extraction

Carry `MM_A` through the **same current yield-extraction mathematics** used by
the baseline branch.

Use the existing `integral_with_stat_error(...)` authority and the same:

- data normalization uncertainty;
- dummy normalization uncertainty;
- active background-profile terms.

Under `no_empirical_residual`, the empirical residual contributions remain
zero for both branches.

Refactor common yield-measurement arithmetic only if necessary to guarantee
the baseline and Method-A branches invoke the same implementation. Any such
refactor must have an explicit regression proving the public baseline yield and
error outputs are unchanged.

The Method-A branch produces private:

```text
YA
YA_statistical_error
YA_total_error
```

It must not replace public `Y0`, public errors, or downstream yield dictionaries.

---

## D. F.6.3 runtime source object for later E.8.4

After yields exist, attach one private source-owned F.6.3 object to each
setting `hist`:

```python
hist["_f6_3_parallel_method_a_source"]
```

Use a versioned schema, e.g.:

```text
f6_3_parallel_method_a_source/v1
```

This object is the later E.8.4 consumer boundary.

It must clearly state:

```text
available / reason
branch_role = parallel_nonproduction_method_a_full_analysis
baseline_production_mutated = false
production_promotion_performed = false
method_b_numerical_dependency = false
empirical_residual_used = false
event_correction_persisted = false
```

Persist/retain only aggregate provenance and branch outputs. Do **not** retain
the identity-to-`C_j` map.

For each canonical child retain enough authoritative source material for E.8.4
to consume without recomputation:

```text
t/phi identity and edges
child validity/status
Lambda integration window
pion_input
B_pi^0
B_pi^A
MM_0
MM_A
Y0
Y0 statistical error
Y0 total error
YA
YA statistical error
YA total error
```

Also retain setting-level authority/provenance:

- selected canonical setting identity;
- accepted F.1 source hashes/fingerprints actually consumed;
- accepted F.3 source hash/fingerprints;
- accepted F.4 source hash/correction/artifact fingerprints;
- transient-factor population/count fingerprint(s), but not factor values;
- current live-cache identity/count fingerprint;
- explicit no-child-renormalization flag;
- explicit baseline-public-output-unchanged flag.

The source object may contain detached ROOT histograms because it is an
in-memory producer/consumer boundary for later E.8.4. Do not write a new
event-factor artifact.

F.6.3 itself does not add E.8.4 pages.

---

# Failure behavior

F.6.3 must fail closed **for the Method-A branch** when accepted authority is
missing or mismatched.

Examples:

```text
F.4 SHA/fingerprint mismatch
F.4 shared reproduction mismatch
F.1/F.3 authority mismatch
unsupported setting
missing/extra/duplicate application identity
live cache t/phi mismatch
live MM/t mismatch
current signed coefficient mismatch
current w0 mismatch
current signed baseline contribution mismatch
non-finite/non-positive C_j
active empirical residual profile
```

In these cases:

- do not fabricate Method-A output;
- do not use `C=1` as a fallback;
- do not partially apply Method A;
- mark `_f6_3_parallel_method_a_source` unavailable with a precise reason;
- leave the authoritative baseline run unchanged and usable.

No child-by-child or setting-by-setting partial Method-A success is allowed
inside one selected setting. The selected setting's Method-A branch is atomic.

---

# Required tests

Add:

```text
testing/test_f6_3_parallel_full_procedure_method_a.py
```

Tests must be deterministic and not require the Jefferson Lab farm.

At minimum cover all of the following.

## Authority / factor reconstruction

1. accepted-shaped F.1/F.3/F.4 fixtures reconstruct a transient identity factor
   map;
2. reproduced F.4 must equal the accepted persisted F.4 aggregate;
3. F.4 source SHA mismatch fails;
4. F.4 correction/artifact fingerprint mismatch fails;
5. F.1/F.3 mismatch fails;
6. row/factor length mismatch fails;
7. duplicate identity fails;
8. missing identity fails;
9. extra live application identity fails;
10. non-finite or non-positive factor fails;
11. no event-level factor values appear in the retained F.6.3 source.

## Live baseline parity

12. matching live cache identity/t/phi/MM/t/coefficient/w0/baseline
    contribution passes;
13. coefficient mismatch fails;
14. w0 mismatch fails;
15. signed baseline contribution mismatch fails;
16. t/phi mismatch fails;
17. MM/t mismatch fails.

## Template behavior

18. omitted multiplier preserves baseline template output exactly;
19. all-one multiplier reproduces baseline template exactly;
20. nontrivial synthetic factors change only pion event weights by `C_j`;
21. `nommcuts` and `allcuts` both receive the same event factor;
22. no child normalization is introduced;
23. no Method-B input exists.

## Full parallel branch

24. baseline public histograms remain unchanged after Method-A construction;
25. same pion input is used for both branches;
26. `MM_0 = input - B_pi^0`;
27. `MM_A = input - B_pi^A`;
28. all-one Method-A factor gives `B_pi^A == B_pi^0`,
    `MM_A == MM_0`, and `YA == Y0` including statistical/total error;
29. nontrivial factors give the expected signed `B_pi^A - B_pi^0` and
    `MM_A - MM_0`;
30. identical Lambda window is used for Y0 and YA;
31. active profile other than `no_empirical_residual` makes F.6.3 unavailable;
32. unsupported kinematic/mode leaves baseline unchanged and marks F.6.3
    unavailable.

## Regression

33. existing public `calculate_yield_data(...)` return structure and baseline
    values/errors are unchanged when the Method-A side branch is unavailable;
34. public `Y0`, E.8.2 source, stage-window yields, scale factors, component
    payloads and baseline ROOT objects are unchanged;
35. no new ROOT/tree traversal is introduced into the F.6.3 helper.

Also rerun the most relevant existing regressions:

```text
testing.test_e8_2_baseline_stage_audit
testing.test_e8_3_detached_method_a_reweighting_audit
testing.test_pion_hgcer_method_a_parent_preserving_correction
testing.test_pion_hgcer_method_a_tphi_propagation
testing.test_pion_hgcer_phase_f_runtime_contract
testing.test_pion_component_dynamic_alignment
```

If another directly affected existing deterministic test module is evident
from the diff, run it too.

---

# Local validation

Use the repository's actual Python command/environment established by
`AGENTS.md` / local tooling.

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

Do not claim ROOT/PyROOT/full-analysis/farm validation from these checks.

---

# Memory/history update

This implementation is not independently source reviewed merely because Codex
finishes it.

After implementation, current-state records should say:

```text
E.8   — ACTIVE
E.8.2 — SOURCE REVIEWED
E.8.3 — SOURCE REVIEWED
F.6.3 — ACTIVE
E.8.4 — BLOCKED pending independent F.6.3 source review
```

`CURRENT.md` must also remove the now-stale pre-push statement that F.6.3 waits
for the already-completed E.8.3 user commit/push.

Add:

```text
docs/memory/phases/f6-3-parallel-full-procedure-method-a.md
```

Record:

- exact starting HEAD;
- changed paths;
- actual runtime/source path;
- scientific ownership;
- frozen baseline behavior;
- factor authority/reconstruction;
- failure behavior;
- tests actually run by Codex and exact results;
- explicit farm boundary;
- status `ACTIVE`;
- NEXT = independent ChatGPT actual-diff/source/runtime-path review.

Regenerate `docs/memory/manifest.json`.

Do not mark F.6.3 SOURCE REVIEWED yourself.

Do not unblock or implement E.8.4 in this task.

---

# Forbidden shortcuts

Do not:

- modify the public baseline pion weight;
- overwrite baseline histograms with Method-A versions;
- alter public yields;
- run a second data/pion ROOT traversal;
- reconstruct F.3/F.4 mathematics independently;
- persist event-level `C_j`;
- derive `C_j` from F.5 aggregate children;
- infer `C_j` from E.8.3 plots/payloads;
- normalize an individual `(t,phi)` child;
- clip/cap/smooth/interpolate F.4 factors;
- use Method B numerically;
- use empirical Fit 1/Fit 2;
- add a silent `C=1` fallback;
- accept missing identities;
- alter pruning;
- alter cuts, component fits, amplitudes, windows, priors or normalization;
- change SIMC;
- change downstream xsec/LT machinery;
- implement E.8.4 presentation;
- commit, push, update remote refs, or run the farm.

---

# Diff audit

Before finishing:

```bash
git status --short
git diff --name-only
git diff --stat
git -c core.safecrlf=false diff --check
```

Confirm every changed path is allowlisted.

Create/refresh root:

```text
kaonlt_review.diff
```

as a complete cumulative diff from committed HEAD
`c9ed0d6b4d0013f7475eedbfca57a096730ea840`, including every intended
untracked file using `git diff --no-index /dev/null <path>` where needed.

Do not stage or commit `kaonlt_review.diff`.

---

# Acceptance criteria

The candidate is ready for independent ChatGPT review only if:

1. committed HEAD is still exactly
   `c9ed0d6b4d0013f7475eedbfca57a096730ea840`;
2. the ordinary baseline analysis remains byte/behavior-equivalent at its
   public outputs;
3. Method A exists only as a private parallel branch;
4. the only scientific difference in that branch is `w0_j -> w0_j*C_j`;
5. `C_j` comes from exact accepted F.4 shared reconstruction;
6. live current event identity and baseline contribution match accepted F.1
   before factor application;
7. event factors are transient and never persisted;
8. no child renormalization exists;
9. no Method B numerical dependency exists;
10. no empirical residual fit is activated;
11. Method-A authority failure leaves baseline production unchanged and marks
    the parallel branch unavailable;
12. all-one factor regression reproduces baseline branch including Y0/Y errors;
13. F.6.3 source object contains the authoritative two-branch material needed
    by later E.8.4 without E.8.4 recomputing the branch;
14. all required deterministic checks pass;
15. memory records F.6.3 as `ACTIVE`, not SOURCE REVIEWED;
16. no farm/runtime claim is made.

---

# Farm boundary

Do **not** run the Jefferson Lab farm for this task.

Per the governing E.8 roadmap, coherent F.6.3 and E.8.4 development proceeds
through local deterministic checks and independent source review first.

The later farm milestone remains:

```text
F.6.3 source reviewed
-> E.8.4 implemented and source reviewed
-> one narrow Q4p4W2p74 / Left / lowe end-to-end farm gate
-> inspect fresh runtime artifacts
-> PASS or one coherent repair
```

---

# Hard stop

Stop and report a blocker rather than broaden scope if:

- accepted F.4 cannot be reproduced from frozen F.1/F.3 authority;
- the live production cache cannot be matched exactly to accepted F.1
  application identities;
- current `w0`/signed source contribution differs from accepted F.1;
- the Method-A branch would require a second tree traversal;
- the baseline production branch would need scientific modification;
- empirical residual machinery would need activation;
- F.4/F.5/F.6.1/F.6.2 source would need modification;
- E.8.4 presentation would need to be implemented to make F.6.3 function;
- ROOT/farm behavior would need to be guessed locally.

Do not redesign adjacent architecture.
