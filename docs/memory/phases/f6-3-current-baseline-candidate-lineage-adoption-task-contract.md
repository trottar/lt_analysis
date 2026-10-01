# KaonLT — F.6.3 current-baseline candidate-lineage adoption for the detached Left/lowe reweighting gate

## 1. Purpose

Unblock the already-existing **private, non-production F.6.3 Method-A branch** so it can consume the now farm-validated **current-baseline candidate lineage** produced by F.4.Refresh.2.

This is the minimum source change required before the `Q4p4W2p74 / Left / lowe` debug full-analysis run can produce the actual baseline-vs-reweighted missing-mass plots and `(t,phi)` yield changes.

Do **not** redesign Method A, F.4, F.5, F.6.3, E.8.4, pion subtraction, or the production analysis.

The current source blocker is concrete:

- `src/cuts/pion_hgcer_method_a_parallel_full_procedure.py::accepted_f6_3_artifact_paths()` still resolves the historical canonical F.3/F.4 artifact basenames;
- F.6.3 calls F.5 `_validate_f4_artifact(..., authority_by_kinematic=None)`, which therefore enforces the historical accepted F.4 authority;
- F.6.3 recomputes F.4 through the public F.4 builder without an explicit current-baseline candidate-F.3 reconstruction authority;
- therefore the newly validated current-baseline candidate F.3/F.4 files cannot be consumed by F.6.3 as source currently stands.

This task fixes only that detached F.6.3 candidate-lineage boundary.

---

## 2. Exact starting state

Required branch:

```text
test
```

Required starting HEAD:

```text
b349967c0d4210a78b144ce6134d3c1f15970245
```

Commit subject:

```text
Adopt memory-bracketed scientific workflow
```

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
```

Requirements:

- branch must be `test`;
- HEAD must be exactly `b349967c0d4210a78b144ce6134d3c1f15970245`;
- inspect the current worktree before editing;
- root `AGENTS.md` and `.codex/` remain local-only/untracked;
- temporary `kaonlt_review*.diff` files may remain untracked and must not be altered;
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

Then read only the source/tests needed for this fix:

- `src/cuts/pion_hgcer_method_a_parallel_full_procedure.py`
- `src/cuts/pion_hgcer_method_a_parent_preserving_correction.py`
- `src/cuts/pion_hgcer_method_a_tphi_propagation.py`
- `testing/test_f6_3_parallel_full_procedure_method_a.py`
- relevant F.4/F.5 tests
- `testing/test_full_background_subtraction_plots.py`
- `testing/test_run_prod_analysis_debug_left_low.py`
- current E.8.4/F.6.3 phase/evidence records only as needed.

Do not reopen earlier science.

---

## 4. Direct farm evidence now accepted for this task

The user supplied:

```text
KaonLT_F4_Refresh2_materialization_Q4p4W2p74_20261001-120157.zip
```

Independent ChatGPT inspection established:

```text
ZIP SHA-256:
4cbdbdd2c0403b961614a8ccbb91104f13fef59194cd48e5e3570aded0954898

bundle git_head:
b349967c0d4210a78b144ce6134d3c1f15970245

required materializer source:
141a3d04f9e5d07be21dba14e0e63212c3990bf1

bundle complete:
true

bundle errors:
[]

unexpected committed files:
[]

scientific gate:
f2_scientific_payload_match = true
f3_scientific_payload_match = true
f4_scientific_payload_match = false
first_changed_stage = F4
```

The materialization manifest also proves:

```text
non_authoritative = true
accepted_authority_mutated = false
production_objects_mutated = false
production_application_performed = false
method_a_promoted = false
complete = true
errors = []
```

### Current F.1 raw SHA-256 identities

```text
Left-lowe:
10086d7d4c42389c9fdd16471c49d59cd189a30980767bcdb87b85914169ef95

Left-highe:
203f4c76f1a251e3e8f231fa3a5e50c9a4fa337420efffffd803205c6d7ea218

Center-lowe:
1593e22b55382b4a9e831d3a1114584e2e3057fcbc4aeeea4aa74c948edf39f8

Center-highe:
5de64b850735ebe70040a703bac999bdd1ac84dc821aba8e1a21ec86ae4b9db3

Right-highe:
77d98006ff81e772c466bf8d10ec83088509a0b8a4d430bded1172bc118bc15f
```

### Validated current-baseline candidate F.3 identity

Basename:

```text
Q4p4W2p74_kaon_pion-background_hgcer-method-a-acceptance-map-current-baseline-candidate.json
```

Identity:

```text
raw SHA-256:
eeb480d9b7ddffab1e97c7f9f27c3099d0780c2bbb4d42a0ee68c077a5752b3d

map fingerprint:
3a9787fc58d26cc0816012bd1b637ad0c8f201b448d54cb1625a841131154728

algorithm fingerprint:
ba29630b2f40a87cbadbe751504ce48f23e2b17a08378a2c8a219f93131cb912

artifact fingerprint:
8d2d12068dfa98922d01dfafedbaa4b994bce0bf1e39ce23b1333872313ef121
```

The materialized F.4 candidate was deliberately constructed with the reviewed diagnostic F.3 reconstruction sentinel:

```text
farm_source_head = 0000000000000000000000000000000000000000
```

That zero-head record is **not** an accepted farm authority. In this task it may be used only to reproduce the already-validated candidate F.4 byte/scientific lineage inside the detached F.6.3 branch. It must never be relabeled as production or general accepted F.3 authority.

### Validated current-baseline candidate F.4 identity

Basename:

```text
Q4p4W2p74_kaon_pion-background_hgcer-method-a-parent-preserving-correction-current-baseline-candidate.json
```

Identity:

```text
raw SHA-256:
1d545924eba89c7f9ffa28028e307aca9b434a89beec06863cf2893887b6b902

correction fingerprint:
bce0cc12ef283b0c361818b436c4e190fc4ea5165436f6c36912cfbb6de53368

artifact fingerprint:
4c935271a0b2723b58b02cc34d82a1e7d5757f7cd36b7707f400c2b28893668a
```

It inherits the candidate F.3 identity above.

The candidate numerical F.4 result reproduces F.4.Refresh.1 exactly. Only `Left-lowe` changes materially:

```text
t=0:
baseline parent sum accepted/current:
0.11020735426697968 -> 0.1081517597541286
relative change:
0.018652062981852853

parent normalization:
4.085998607893882 -> 4.114018603391793
relative change:
0.0068575636427722795

t=1:
difference only at floating-point scale

t=2:
baseline parent sum accepted/current:
0.1113882513418981 -> 0.10645779207740595
relative change:
0.044263728042183405

parent normalization:
1.988991016567433 -> 2.0555330576308606
relative change:
0.033455174261301984
```

---

## 5. Scientific ownership

Preserve all of the following:

- random/dummy subtraction;
- slow-proton subtraction;
- pion-component definitions and fits;
- dynamic alignment;
- baseline pion weight `w0`;
- F.1 event-contract science;
- F.2 representation science;
- F.3 `hgcer3` map mathematics;
- F.4 parent-preserving mathematics;
- F.5 `(t,phi)` propagation mathematics;
- frozen F.6.2 science;
- F.6.3 weighting formula;
- E.8.4 consumer mathematics;
- Method B diagnostic-only status;
- active `no_empirical_residual` profile;
- SIMC normalization;
- cuts, windows, priors, templates, binning;
- efficiencies, acceptance, yields, L/T separation, cross sections.

The only scientific operation F.6.3 may continue to perform is the existing private:

```text
w0 -> w0 * C
```

with parent-t normalization preserved and no child `(t,phi)` renormalization.

This task must not promote Method A.

---

## 6. Required design

### 6.1 Keep historical F.3/F.4 accepted authorities untouched

Do **not** modify:

```text
src/cuts/pion_hgcer_method_a_parent_preserving_correction.py
    ACCEPTED_F3_RUNTIME_AUTHORITY_BY_KINEMATIC

src/cuts/pion_hgcer_method_a_tphi_propagation.py
    ACCEPTED_F4_RUNTIME_AUTHORITY_BY_KINEMATIC
```

Those remain the historical detached F.3/F.4 accepted-authority records and preserve their earlier closures.

Do not replace historical accepted F.3/F.4 artifact files.

### 6.2 F.6.3 owns one explicit current-baseline candidate lineage

In:

```text
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
```

add narrowly named source-owned constants/records for **F.6.3 candidate-lineage reconstruction only**.

They must pin exactly:

- candidate F.3 basename and identities above;
- candidate F.4 basename and identities above;
- all five current F.1 raw SHA values above;
- the F.4.Refresh.2 validation source HEAD `b349967c0d4210a78b144ce6134d3c1f15970245`;
- the exact zero-head F.3 diagnostic reconstruction sentinel required to reproduce the already-materialized F.4 candidate.

Names and comments must make clear that:

- this is a detached F.6.3 current-baseline candidate lineage;
- it is not the general accepted F.3 or F.4 authority;
- the all-zero F.3 farm head remains a candidate-construction sentinel only;
- it does not authorize production promotion.

Do not add dynamic discovery or a fallback search.

### 6.3 Resolve exact candidate F.3/F.4 paths

For `Q4p4W2p74`, `accepted_f6_3_artifact_paths(...)` must resolve the F.3 and F.4 inputs to the two exact validated candidate basenames listed above.

The five F.1 paths remain the existing canonical current F.1 files.

No globbing, newest-file selection, environment fallback, archive extraction, or historical fallback is permitted.

Unsupported kinematics must continue to fail closed.

### 6.4 Validate the candidate F.4 through F.5's existing validator

In `reconstruct_transient_factor_map(...)`, call the existing:

```text
_f5._validate_f4_artifact(...)
```

with an explicit F.6.3-only candidate-F.4 authority record rather than `None`.

That record must pin:

- candidate F.4 raw SHA;
- candidate F.4 correction fingerprint;
- candidate F.4 artifact fingerprint;
- candidate F.3 raw SHA/map/algorithm/artifact fingerprints;
- all five current F.1 raw SHA values;
- farm validation source HEAD `b349967c0d4210a78b144ce6134d3c1f15970245`.

Do not weaken any F.5 validation.

### 6.5 Reproduce the validated candidate F.4 exactly

When F.6.3 calls:

```text
_f4.build_pion_hgcer_method_a_parent_preserving_correction_with_review_data(...)
```

pass the exact F.6.3-only candidate-F.3 reconstruction record that reproduces the materialized F.4 candidate, including the all-zero sentinel farm head.

This override exists only to reproduce the already-validated candidate artifact inside the detached F.6.3 branch.

Do not:

- install it into F.4's global accepted-authority constant;
- call it an accepted farm F.3 authority;
- expose it to production;
- use it outside F.6.3 candidate-lineage reconstruction.

The existing exact `recomputed == persisted` and fingerprint equality gate must remain unchanged. The task succeeds only if the validated candidate F.4 reproduces exactly.

### 6.6 Preserve live-cache parity

Do not weaken `validate_live_cache_parity(...)`.

It must still require exact identity inventory and current baseline parity for:

```text
t_index
phi_index
analysis_MM
analysis_t
signed_source_coefficient
baseline_pion_weight_w0
signed_baseline_event_contribution
```

before any Method-A multiplier is used.

### 6.7 Provenance

Extend F.6.3 aggregate provenance only as needed to make the boundary explicit.

At minimum the returned provenance must make it unambiguous that:

```text
branch_role = parallel_nonproduction_method_a_full_analysis
current_baseline_candidate_lineage = true
production_promotion_performed = false
event_correction_persisted = false
```

and identify the candidate F.3/F.4 raw hashes.

Do not change downstream science or public baseline outputs.

---

## 7. Allowed files

Substantive source:

```text
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
```

Focused tests:

```text
testing/test_f6_3_parallel_full_procedure_method_a.py
```

Integration/regression tests may change only if necessary to assert the unchanged consumer/debug behavior:

```text
testing/test_full_background_subtraction_plots.py
testing/test_run_prod_analysis_debug_left_low.py
```

Minimal durable memory after implementation only:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
```

One concise new phase record is allowed:

```text
docs/memory/phases/f6-3-current-baseline-candidate-lineage-adoption.md
```

No new workflow/maintenance/decision/investigation documents. No new generic infrastructure.

If another scientific/runtime source file appears necessary, STOP and report the blocker before expanding scope.

---

## 8. Required before/after behavior

### Before

At current HEAD, F.6.3 resolves historical F.3/F.4 authority and fails against the current-baseline lineage.

### After

For `Q4p4W2p74`, F.6.3:

1. loads the exact five current canonical F.1 files;
2. loads the exact validated candidate F.3 file;
3. loads the exact validated candidate F.4 file;
4. validates candidate F.4 against the F.6.3-only source-owned candidate record;
5. rebuilds F.4 from current F.1 + candidate F.3 using the exact diagnostic reconstruction sentinel;
6. requires byte/scientific persisted/recomputed equality;
7. joins correction factors to the current live cache;
8. requires exact live-cache parity;
9. exposes transient factors only to the existing private non-production branch;
10. leaves baseline public production output unchanged.

No fallback to historical F.4 is allowed if candidate validation fails.

---

## 9. Required tests

Add or update focused tests for at least:

1. exact current-baseline candidate F.3/F.4 path resolution;
2. exact candidate F.3 raw/map/algorithm/artifact identity pin;
3. exact candidate F.4 raw/correction/artifact identity pin;
4. exact five-current-F.1 raw hash inventory;
5. F.5 validator receives the explicit F.6.3 candidate-F.4 authority;
6. F.4 builder receives the explicit candidate-F.3 reconstruction sentinel;
7. sentinel remains all-zero farm head and is labelled candidate reconstruction only;
8. historical F.4/F.5 global accepted-authority constants are not modified by F.6.3;
9. wrong candidate F.4 raw hash fails closed;
10. wrong candidate F.3 fingerprint fails closed;
11. recomputed/persisted F.4 mismatch still fails closed;
12. live-cache identity mismatch still fails closed;
13. live-cache baseline numerical mismatch still fails closed;
14. transient correction factors remain positive/finite and use the unchanged `w0*C` branch;
15. no event correction is persisted;
16. no Method-B numerical dependency appears;
17. baseline public output remains unchanged;
18. unsupported kinematic fails closed.

Run existing regression tests as well.

---

## 10. Local validation

Run at minimum:

```bash
python -B -m py_compile \
  src/cuts/pion_hgcer_method_a_parallel_full_procedure.py \
  testing/test_f6_3_parallel_full_procedure_method_a.py \
  testing/test_full_background_subtraction_plots.py \
  testing/test_run_prod_analysis_debug_left_low.py

python -B -m unittest \
  testing.test_f6_3_parallel_full_procedure_method_a -v

python -B -m unittest \
  testing.test_pion_hgcer_method_a_parent_preserving_correction -v

python -B -m unittest \
  testing.test_pion_hgcer_method_a_tphi_propagation -v

python -B -m unittest \
  testing.test_run_prod_analysis_debug_left_low -v

python -B -m unittest \
  testing.test_full_background_subtraction_plots -v

python -B tools/update_memory_manifest.py --root . --write
python -B tools/update_memory_manifest.py --root . --check

git -c core.safecrlf=false diff --check
```

If ROOT-dependent tests are skipped locally, report exact skips. Do not treat local tests as farm validation.

Do not run `main.py` locally.

---

## 11. Diff audit

Before stopping:

```bash
git status --short
git -c core.safecrlf=false diff --stat
git -c core.safecrlf=false diff --no-ext-diff
```

Audit every changed file.

No change is allowed to:

- F.4 mathematics;
- F.5 mathematics;
- F.6.3 weighting mathematics;
- E.8.4 consumer mathematics;
- production subtraction/correction;
- `run_Prod_Analysis.sh`;
- Method B;
- normalizations/cuts/binning/yields/efficiencies/SIMC.

---

## 12. Status / memory

This source change does **not** close F.6.3 or E.8.4.

After implementation only:

```text
F.4.Refresh.2 candidate materialization — CLOSED / RUNTIME VALIDATED
```

for the detached materialize/verify/package gate represented by the supplied ZIP.

The new F.6.3 candidate-lineage adoption is:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

after independent ChatGPT actual-diff review passes.

`CURRENT.md` should contain exactly one scientific NEXT:

```text
NEXT — run Q4p4W2p74 / Left / lowe through the tracked -d debug full-analysis path and inspect baseline-versus-reweighted missing-mass spectra and per-(t,phi) yield changes.
```

Do not add maintenance NEXTs.

---

## 13. Farm boundary

Codex must not run the farm.

After source review and user commit/push, the next farm command is the existing tracked debug path:

```text
./run_Prod_Analysis.sh -d 4p4 2p74
```

That path must preserve paired low/high canonical preflight, execute only the actual `Left / lowe` full analysis, and stop before high-epsilon full processing.

Farm evidence must then show:

- baseline vs Method-A-reweighted missing-mass spectrum;
- baseline vs reweighted spectrum per canonical `t`;
- baseline and reweighted yields per `(t,phi)`;
- absolute yield change per `(t,phi)`;
- fractional yield change per `(t,phi)`;
- parent-t preservation;
- procedure-PDF representation of the effect;
- no baseline production mutation.

---

## 14. Forbidden shortcuts

Do not:

- overwrite the historical accepted F.3/F.4 JSONs;
- copy candidate bytes over historical canonical files;
- modify F.4/F.5 global accepted-authority constants;
- remove or weaken F.3/F.4/F.5 fingerprint checks;
- remove the recomputed-vs-persisted F.4 equality check;
- relax live-cache parity;
- silently fall back to historical F.4;
- introduce environment-selected candidate files;
- search directories for newest artifacts;
- make Method A production-authoritative;
- involve Method B numerically;
- alter pion weights outside the F.6.3 private branch;
- alter any yield extraction or normalization;
- create more workflow infrastructure;
- commit or push.

---

## 15. Acceptance criteria

The candidate passes source review only if:

1. starting HEAD is exactly `b349967c0d4210a78b144ce6134d3c1f15970245`;
2. only the allowed files change;
3. historical F.3/F.4 authority constants remain byte-identical;
4. F.6.3 resolves the exact validated current-baseline candidate F.3/F.4 files;
5. all exact candidate identities from the supplied farm bundle are pinned;
6. the zero-head record is used only for candidate-F.3 reconstruction needed to reproduce the validated candidate F.4;
7. candidate F.4 is independently validated by F.5's existing validator using an F.6.3-only authority record;
8. persisted/recomputed candidate F.4 equality remains exact;
9. live-cache parity remains exact;
10. `w0*C` mathematics is unchanged;
11. no child renormalization is introduced;
12. no production output is mutated;
13. Method B remains excluded;
14. focused and integration tests pass subject to exact reported skips;
15. diff audit finds no unrelated changes;
16. no farm run, commit, or push occurs.

---

## 16. Hard stop

After implementation, deterministic local checks, minimal memory update, and diff audit:

**STOP.**

Do not run the farm, create another workflow phase, broaden to canonical five settings, or attempt production promotion.

The immediate next operation after independent ChatGPT review + user commit/push is the one `Q4p4W2p74 / Left / lowe` debug full-analysis run that must finally produce the before/after plots and yields.
