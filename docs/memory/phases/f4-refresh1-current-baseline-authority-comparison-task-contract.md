# KaonLT F.4.Refresh.1 — Current-Baseline Method-A Authority Compatibility Audit

## 1. Purpose

Add one detached, deterministic diagnostic that determines the minimum scientific
authority refresh required after E.8.4.Fix.4 closed the persisted-alignment
semantic-version bug but the fresh `Q4p4W2p74 / Left / lowe` full-analysis gate
still failed F.6.3 with:

```text
f6_3_f4_shared_reproduction_failed:f3_fingerprint_input_content_mismatch
```

New direct read-only farm evidence establishes:

### Fix.4 runtime behavior

The current Left/lowe setting-wide alignment record carries:

```text
alignment_semantics_version = pion_component_dynamic_alignment_semantics/v2
persistence_status = rejected_stale_then_created
```

with rejection reasons including:

```text
alignment_semantics_version mismatch
```

Therefore the E.8.4.Fix.4 semantic-version gate executed as designed at runtime.
Fix.4 itself is now:

```text
CLOSED / RUNTIME VALIDATED
```

This closure is narrow: it validates stale semantic-cache rejection and
recomputation, not F.6.3/E.8.4.

### Exact frozen/current F.1 semantic comparison

A streaming record-by-record comparison between the accepted frozen
`Left / lowe` F.1 artifact and the current Fix.4-generated F.1 artifact found:

```text
F3 TRAINING INPUTS
records_compared = 52397
identity_mismatch_count = 0
projection_mismatch_count = 0
projection_match = True

F3 APPLICATION INPUTS
records_compared = 55380
identity_mismatch_count = 0
projection_mismatch_count = 0
projection_match = True
```

The compared F.3 training projection is exactly:

```text
source_label
entry_index
response_class
P_hgcer_npeSum
t_index
t_low
t_high
SHMS_delta
SHMS_xptar
SHMS_yptar
P_hgcer_xAtCer
P_hgcer_yAtCer
```

The compared F.3 application projection is exactly:

```text
source_label
entry_index
P_hgcer_npeSum
t_index
t_low
t_high
SHMS_delta
SHMS_xptar
SHMS_yptar
P_hgcer_xAtCer
P_hgcer_yAtCer
```

Those are the actual sanitized F.1 fields consumed by the F.2/F.3 acceptance
machinery.

By contrast, every current Left/lowe F.1 application record differs in the
F.4/F.6.3 baseline projection:

```text
F4/F6.3 BASELINE NUMERICS
records_compared = 55380
identity_mismatch_count = 0
projection_mismatch_count = 55380
projection_match = False
```

First mismatch:

```text
identity = ["prompt", 0]

changed fields:
analysis_MM
baseline_pion_weight_w0
signed_baseline_event_contribution

frozen:
analysis_MM = 0.9150867495356926
baseline_pion_weight_w0 = 0.0434198336219789
signed_baseline_event_contribution = 8.049822627351833e-06

current:
analysis_MM = 0.9150867495404144
baseline_pion_weight_w0 = 0.04235726703122994
signed_baseline_event_contribution = 7.852828031293563e-06
```

The current F.1 file SHA-256 remains:

```text
10086d7d4c42389c9fdd16471c49d59cd189a30980767bcdb87b85914169ef95
```

### Source-audit consequence

F.3 sanitizes the raw F.1 records before model construction. The exact
record-by-record comparison above proves that the actual F.2/F.3 numerical
training and application inputs for Left/lowe are unchanged.

F.4 is different. Its current source uses each
`signed_baseline_event_contribution` as `b`, constructs:

```text
B = sum(b)
U = sum(b * raw_factor)
parent_normalization = U / B
C_j = raw_factor_j / parent_normalization
```

Therefore the accepted F.4 correction cannot simply be assumed current when all
55,380 Left/lowe baseline application records changed.

The existing F.2/F.3/F.4 artifacts also intentionally bind to broad F.1/F.2/F.3
provenance and source-file fingerprints. Do not weaken those fail-closed gates
merely to make F.6.3 proceed.

This task adds a **diagnostic comparator only**. It rebuilds candidate F.2 and
F.3 from the current five F.1 artifacts, then rebuilds a candidate F.4 against
the current baseline using an explicitly detached, in-memory candidate-F.3
authority override. It compares those candidates to the accepted F.2/F.3/F.4
artifacts and emits structured evidence. It does not update any accepted
authority and does not change production physics.

---

## 2. Exact starting state

Required branch:

```text
test
```

Required committed starting HEAD:

```text
a8bd4dc4e0990dd45cfade0079dc1ce40bdd36ce
```

Required commit subject:

```text
E8.4 Fix.4: re-pin validation bundle profile
```

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
```

Hard requirements:

- branch must be `test`;
- committed HEAD must be exactly
  `a8bd4dc4e0990dd45cfade0079dc1ce40bdd36ce`;
- inspect the existing worktree before editing;
- temporary root `kaonlt_review*.diff` files may remain untracked and must not
  be removed, rewritten, staged, or treated as candidate files;
- root `AGENTS.md` and `.codex/` are local-only/untracked and must not be
  removed;
- if unrelated tracked edits exist, STOP and report the blocker;
- do not reset, stash, clean, discard, commit, push, run the full analysis, or
  mutate production artifacts.

---

## 3. Mandatory repository-memory startup

Read in this exact order:

1. root `AGENTS.md`;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read only task-relevant records/source, including:

- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/phase-f6-method-a-production-promotion.md`
- `docs/memory/phases/e8-4-fix4-alignment-cache-semantics.md`
- `docs/memory/evidence/e8-4-fix3-left-lowe-stale-alignment-cache-runtime-blocker.md`
- the current E.8.4 phase record
- `src/cuts/pion_hgcer_method_a_acceptance_representation.py`
- `src/cuts/pion_hgcer_method_a_acceptance_map.py`
- `src/cuts/pion_hgcer_method_a_parent_preserving_correction.py`
- `src/cuts/pion_hgcer_method_a_tphi_propagation.py`
- `src/cuts/pion_hgcer_method_a_parallel_full_procedure.py`
- task-relevant existing F.2/F.3/F.4 tests.

Do not redesign adjacent architecture.

---

## 4. Status and scientific ownership

Preserve historical closures:

- F.1 through F.6.2, including F.6.2.Fix.5 —
  **CLOSED / RUNTIME VALIDATED**

Do not downgrade those historical phases. Their accepted artifacts remain valid
evidence for the baseline/science state under which they were accepted.

Update:

- E.8.4.Fix.4 persisted-alignment semantic-version repair —
  **CLOSED / RUNTIME VALIDATED**

This is a narrow closure supported by the fresh alignment persistence evidence.
It does **not** close F.6.3, E.8.4, final E.8, or F.6.4.

Preserve:

- E.8 — **ACTIVE**
- E.8.2 — **SOURCE REVIEWED**
- E.8.3 — **SOURCE REVIEWED**
- F.6.3 — **SOURCE REVIEWED**
- E.8.4 — **SOURCE REVIEWED**
- final E.8 — **BLOCKED**
- F.6.4 — **BLOCKED**
- lifecycle-hook dispatch — **BLOCKED / DEFERRED**

Create:

- F.4.Refresh.1 current-baseline Method-A authority compatibility audit —
  **ACTIVE** pending independent ChatGPT actual-diff review.

The fresh Left/lowe E.8.4 gate remains **BLOCKED**. The new evidence changes
the diagnosis: Fix.4 worked, F.3 semantic inputs are unchanged, but current
F.4/F.6.3 baseline numerics differ from the accepted frozen baseline.

---

## 5. Scientific boundaries

This task is diagnostics/validation only.

Do not change or reinterpret:

- random/dummy subtraction;
- slow-proton subtraction;
- pion component definitions;
- pion-control fits;
- dynamic alignment configuration or scan;
- baseline pion-weight formula;
- accepted `w0` producer logic;
- Method-A response definition;
- F.2 candidate definitions;
- F.3 accepted basis `hgcer3`;
- F.3 model mathematics;
- F.4 parent-preserving mathematics;
- F.5 propagation mathematics;
- F.6.1/F.6.2 science;
- F.6.3 branch mathematics;
- E.8.4 consumer;
- Method B;
- SIMC normalization;
- canonical binning;
- yield extraction;
- efficiencies;
- acceptance;
- L/T separation;
- cross sections;
- active `no_empirical_residual` profile.

Do not update source-owned accepted F.3 or F.4 runtime-authority constants in
this task.

Do not promote Method A.

---

## 6. Allowed substantive changes

Add exactly one diagnostic implementation:

```text
testing/compare_method_a_current_baseline_authority.py
```

Add exactly one focused test:

```text
testing/test_compare_method_a_current_baseline_authority.py
```

No existing scientific source file may change.

In particular, do not modify:

```text
src/cuts/pion_hgcer_method_a_acceptance_representation.py
src/cuts/pion_hgcer_method_a_acceptance_map.py
src/cuts/pion_hgcer_method_a_parent_preserving_correction.py
src/cuts/pion_hgcer_method_a_tphi_propagation.py
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
src/cuts/pion_component_fits.py
src/cuts/pion_hgcer_event_contract.py
src/cuts/pion_hgcer_method_a_acceptance_contract.py
src/cuts/pion_hgcer_refinement_method_a.py
src/cuts/full_background_subtraction_plots.py
src/cuts/rand_sub.py
src/main.py
run_Prod_Analysis.sh
```

If an existing scientific source file must change to implement the comparator,
STOP and report the blocker.

---

## 7. Required diagnostic interface

Implement a pure-Python CLI.

It must accept explicit input paths for:

- all five current canonical F.1 artifacts:
  - Left/lowe
  - Left/highe
  - Center/lowe
  - Center/highe
  - Right/highe
- one accepted F.2 artifact;
- one accepted F.3 artifact;
- one accepted F.4 artifact;
- one output JSON path.

Do not embed farm paths in code.

Suggested interface:

```text
--f1 Left-lowe=/path/...
--f1 Left-highe=/path/...
--f1 Center-lowe=/path/...
--f1 Center-highe=/path/...
--f1 Right-highe=/path/...
--accepted-f2 /path/...
--accepted-f3 /path/...
--accepted-f4 /path/...
--output /path/diagnostic.json
```

Equivalent deterministic argument spelling is acceptable.

The tool must:

1. verify the exact canonical-five F.1 inventory;
2. load JSON with strict failures;
3. compute SHA-256 for every input file from raw bytes;
4. rebuild candidate F.2 using only the existing public F.2 builder;
5. derive deterministic candidate-F.2 serialized bytes using the repository's
   existing writer format and compute its source-file SHA-256;
6. rebuild candidate F.3 from the current F.1 artifacts plus candidate F.2
   using only the existing public F.3 builder;
7. derive deterministic candidate-F.3 serialized bytes and SHA-256 using the
   existing writer format;
8. rebuild candidate F.4 from the current F.1 artifacts plus candidate F.3
   using only the existing public F.4 builder;
9. for step 8 only, pass an explicit **in-memory diagnostic authority override**
   describing the candidate F.3 just built, because this tool is comparing a
   not-yet-accepted candidate. The override must never be written into
   production source, accepted authority constants, or ordinary runtime state;
10. compare candidate versus accepted scientific payloads;
11. write one JSON-safe structured diagnostic;
12. exit nonzero only for malformed/missing/internally inconsistent input or
    comparator execution failure, not merely because candidate and accepted
    scientific values differ.

No ROOT/PyROOT may be imported.

---

## 8. Required comparison semantics

### 8.1 F.1 input inventory

Persist:

- each setting ID;
- current F.1 raw file SHA-256;
- current F.1 stable-content fingerprint;
- current F.1 contract fingerprint;
- current training/application record counts;
- current training/application population fingerprints.

This section is provenance only.

### 8.2 F.2 comparison

The candidate F.2 will necessarily carry new F.1 provenance when Left/lowe
broad F.1 content differs. Therefore do **not** compare whole-artifact hashes as
a scientific equality test.

Compare and report exact equality for the F.2 numerical/scientific payload:

- candidate definitions;
- algorithm configuration;
- algorithm fingerprint;
- response-support result;
- all 15 groups':
  - geometry;
  - counts;
  - candidate metrics;
  - application-support metrics;
- candidate summaries;
- recommendation.

Exclude only:

- generated timestamp/git-status/path wrapper provenance;
- raw input-path strings;
- F.1 raw file hashes;
- broad stable-F.1 fingerprints;
- wrapper/artifact fingerprint;
- representation fingerprint portions whose only difference is inherited input
  provenance.

Report:

```text
f2_scientific_payload_match
f2_first_mismatch_path
```

and enough structured detail to inspect a mismatch.

### 8.3 F.3 comparison

Compare and report exact equality for the accepted-basis F.3 numerical/scientific
payload:

- accepted basis;
- ordered features;
- algorithm configuration;
- algorithm fingerprint;
- all 15 model identities/geometries;
- scaler median/divisor;
- logistic coefficients;
- intercept provenance value;
- optimizer/objective/convergence fields;
- support metrics;
- F.2 support-continuity fields;
- relative-response definition;
- any persisted review-grid/model numerical fields that affect or describe the
  accepted map.

Exclude only inherited source-file/input-provenance hashes/fingerprints and
wrapper provenance.

Report:

```text
f3_scientific_payload_match
f3_first_mismatch_path
```

The comparison must fail closed on missing fields; do not silently drop
unrecognized scientific fields.

### 8.4 F.4 comparison

The candidate F.4 is expected to use current F.1 baseline quantities. Compare
the accepted and candidate F.4 parent inventory by exact `(setting_id,t_index)`.

For every parent report accepted, candidate, absolute difference, and relative
difference where defined for at minimum:

- application event count;
- zero baseline-weight count;
- in-support count;
- OOD count/fraction;
- baseline parent sum;
- absolute baseline parent sum;
- baseline cancellation ratio;
- raw-shape parent sum;
- parent normalization;
- adjusted parent sum;
- closure residual;
- OOD final correction;
- raw-shape-factor summary;
- correction-factor summary;
- source diagnostics;
- canonical-phi diagnostics.

Also report setting/global maxima for:

```text
abs(parent_normalization_candidate - parent_normalization_accepted)
relative parent-normalization change
correction-factor-summary change
baseline-parent-sum change
```

Do not define an acceptance threshold in this task. This is evidence, not a
promotion decision.

### 8.5 Dependency decision summary

The output JSON must include a mechanical summary only:

```text
f2_scientific_payload_match
f3_scientific_payload_match
f4_scientific_payload_match
first_changed_stage
```

`first_changed_stage` is:

- `"F2"` if F.2 scientific payload differs;
- else `"F3"` if F.3 scientific payload differs;
- else `"F4"` if F.4 scientific payload differs;
- else `"none"`.

Do not convert this to a production/promotion verdict.

---

## 9. Forbidden comparator behavior

Do not:

- update accepted artifacts;
- overwrite any input;
- modify F.2/F.3/F.4 source-owned authority constants;
- weaken F.2/F.3/F.4 validation;
- bypass current F.1 validation;
- omit current F.1 raw SHA provenance;
- special-case Left/lowe numerical values;
- hard-code the expected output;
- use the old accepted F.1 in candidate construction;
- use Method B;
- import ROOT;
- run `main.py`;
- mutate production files;
- infer that matching F.2/F.3 numerics automatically approves F.4;
- infer that any F.4 difference is automatically a systematic uncertainty.

---

## 10. Output and write safety

The tool writes only the explicit `--output` path.

Refuse to overwrite an existing output unless an explicit `--overwrite` flag is
provided.

The output must be deterministic apart from an optional clearly separated
non-fingerprinted invocation/provenance block.

The diagnostic scientific payload itself must contain no timestamps.

No output is an accepted F.2/F.3/F.4 artifact.

Use a diagnostic schema name such as:

```text
method_a_current_baseline_authority_comparison/v1
```

---

## 11. Required local deterministic tests

The focused test must not require ROOT or farm artifacts.

Use small synthetic JSON fixtures and monkeypatching of the existing public
F.2/F.3/F.4 builders where appropriate to test orchestration separately from
the already-tested scientific builders.

At minimum test:

1. canonical-five F.1 inventory required;
2. duplicate/missing/unexpected setting rejected;
3. input raw SHA-256 recorded from exact bytes;
4. candidate F.2 SHA uses the deterministic writer representation;
5. candidate F.3 SHA uses the deterministic writer representation;
6. F.4 candidate authority override is in-memory only and exactly matches the
   just-built candidate F.3 identity;
7. F.2 scientific-equality projection ignores only declared provenance fields;
8. F.2 scientific numerical mismatch is detected;
9. F.3 scientific-equality projection ignores only declared provenance fields;
10. F.3 coefficient/scaler/support mismatch is detected;
11. F.4 parent comparison detects a parent-normalization change;
12. `first_changed_stage` resolves F2 -> F3 -> F4 -> none in that order;
13. existing output refuses overwrite without `--overwrite`;
14. malformed JSON fails nonzero;
15. no input file is modified.

Run at minimum:

```bash
python -B -m py_compile \
  testing/compare_method_a_current_baseline_authority.py \
  testing/test_compare_method_a_current_baseline_authority.py

python -B -m unittest \
  testing.test_compare_method_a_current_baseline_authority -v

python -B -m unittest \
  testing.test_pion_hgcer_method_a_acceptance_representation -v

python -B -m unittest \
  testing.test_pion_hgcer_method_a_acceptance_map -v

python -B -m unittest \
  testing.test_pion_hgcer_method_a_parent_preserving_correction -v

python -B -m unittest \
  testing.test_f6_3_parallel_full_procedure_method_a -v

python -B tools/update_memory_manifest.py --root . --write
python -B tools/update_memory_manifest.py --root . --check
python -B tools/check_memory_health.py --root .
python -B tools/memory_bootstrap.py --root . --json
python -B -m unittest testing.test_memory_health -v

git -c core.safecrlf=false diff --check
```

If the actual test module names differ, discover and use the repository's
existing focused modules rather than inventing substitutes.

Report skips exactly.

Tests do not establish farm/runtime validation.

---

## 12. Warranted durable memory/history

Allowed existing memory files:

```text
docs/memory/CURRENT.md
docs/memory/roadmap/STATUS.md
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/phases/e8-4-fix4-alignment-cache-semantics.md
docs/memory/manifest.json
```

Allowed new records:

```text
docs/memory/phases/f4-refresh1-current-baseline-authority-comparison-task-contract.md
docs/memory/phases/f4-refresh1-current-baseline-authority-comparison.md
docs/memory/evidence/e8-4-fix4-left-lowe-runtime-closure-and-f4-baseline-divergence.md
```

The new evidence record must distinguish direct farm evidence from source
interpretation.

Record direct evidence:

- Fix.4 current-semantics setting-wide record was
  `rejected_stale_then_created`;
- stale rejection included `alignment_semantics_version mismatch`;
- current Left/lowe F.1 SHA is
  `10086d7d4c42389c9fdd16471c49d59cd189a30980767bcdb87b85914169ef95`;
- F.3 training projection: 52,397 records, zero identity mismatch, zero
  projection mismatch;
- F.3 application projection: 55,380 records, zero identity mismatch, zero
  projection mismatch;
- F.4/F.6.3 baseline projection: 55,380 records, zero identity mismatch,
  55,380 projection mismatches;
- exact first mismatch fields and frozen/current values from Section 1;
- fresh procedure E.8.4 remained unavailable with
  `f6_3_f4_shared_reproduction_failed:f3_fingerprint_input_content_mismatch`.

Record source interpretation separately:

- F.3 numerically consumes sanitized acceptance-coordinate projections;
- F.4 parent normalization consumes `signed_baseline_event_contribution`;
- therefore the next evidence gate must determine whether the first
  scientifically changed detached stage is F.2, F.3, or F.4 rather than
  weakening F.6.3 provenance.

Do not claim F.4/F.5/F.6.1/F.6.2 current-baseline refresh is already required
until the comparator runs. Historical closures remain intact.

After implementation:

```text
F.4.Refresh.1 — ACTIVE
```

pending independent ChatGPT actual-diff review.

`CURRENT.md` exact NEXT:

```text
NEXT — independent ChatGPT actual-diff/source-runtime-path review of the F.4.Refresh.1 current-baseline Method-A authority comparator.
```

Do not add another NEXT.

Regenerate `docs/memory/manifest.json`.

Do not modify `MEMORY.md`, `USER.md`, `LEARNINGS.md`, decisions, or unrelated
phase/evidence records unless a concrete blocker requires it; STOP before
expanding scope.

---

## 13. Required diff audit

Before stopping:

```bash
git status --short
git -c core.safecrlf=false diff --stat
git -c core.safecrlf=false diff --no-ext-diff
```

Audit every intended new file with `git diff --no-index /dev/null ... || true`.

Expected substantive candidate paths:

```text
testing/compare_method_a_current_baseline_authority.py
testing/test_compare_method_a_current_baseline_authority.py
```

plus only the contract-allowed memory/history paths.

No unrelated source may change.

---

## 14. Required cumulative review bundle

Create one fresh temporary repository-root review bundle:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

It must contain byte-faithfully:

1. branch;
2. committed HEAD;
3. `git status --short`;
4. complete `git diff --stat`;
5. complete tracked `git -c core.safecrlf=false diff --no-ext-diff`;
6. complete `git diff --no-index /dev/null ...` for every intended new file;
7. exact local validation commands/results;
8. final changed-path inventory.

Do not normalize or rewrite raw diff text.

The review bundle is temporary and must not be staged or committed.

---

## 15. Acceptance criteria

The local candidate is acceptable for independent source review only if:

1. committed HEAD remains exactly
   `a8bd4dc4e0990dd45cfade0079dc1ce40bdd36ce`;
2. no scientific/runtime source file changes;
3. only the new detached comparator and focused test are substantive code;
4. comparator requires exact canonical-five current F.1 inputs;
5. comparator records exact raw input SHA-256 values;
6. candidate F.2 is built by the existing public F.2 builder;
7. candidate F.3 is built by the existing public F.3 builder;
8. candidate F.4 is built by the existing public F.4 builder;
9. the candidate-F.3 authority override exists only in memory inside this
   diagnostic path;
10. accepted source authority constants remain untouched;
11. F.2/F.3 scientific comparisons exclude only explicitly enumerated
    provenance fields;
12. F.4 parent comparison exposes current-baseline numerical differences;
13. output is deterministic and JSON-safe;
14. no input/production artifact is mutated;
15. Fix.4 is recorded as `CLOSED / RUNTIME VALIDATED` only for its narrow cache
    semantics gate;
16. historical F.1--F.6.2 closures are not downgraded;
17. F.6.3/E.8.4 remain unvalidated at runtime;
18. local deterministic tests pass subject to explicitly reported skips;
19. no commit, push, full farm analysis, authority update, or Method-A promotion
    occurs;
20. a fresh cumulative review bundle is produced.

---

## 16. Hard stop

After implementation, deterministic local checks, warranted memory/history
updates, manifest regeneration, diff audit, and creation of the fresh
timestamped cumulative review bundle:

**STOP.**

Do not commit, push, run the farm comparator, regenerate accepted F.2/F.3/F.4
artifacts, update runtime authorities, or rerun the full analysis.

NEXT is independent ChatGPT actual-diff review of the comparator candidate.
