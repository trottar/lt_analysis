# KaonLT F.4.Refresh.2 — Current-Baseline Method-A Candidate Materialization

## 1. Purpose

Implement one detached, fail-closed materializer for the **current-baseline
Method-A authority lineage** established by F.4.Refresh.1.

The new farm comparator result is direct runtime evidence:

```text
KaonLT_F4_Refresh1_authority_comparison_Q4p4W2p74_20260930-134654.json
SHA-256:
c16287e9192288ed5f116b5eebf95f49261d7fd09755ee923070154e047cf7b5
```

Its accepted artifact identities are:

```text
F.2 accepted JSON:
87e0a16870ea34a6611702e1bd4fc02c4dc6a8daa132e75b725a1b782a9a31da

F.3 accepted JSON:
04a9a576767ba3f81c8c03daab2461998894c655f0d875bd532ad30cb46c0d95

F.4 accepted JSON:
adcc01900b8adc211e296aaedaa314fdb86207d51483ebfd4b3edf69f558f188
```

The comparator rebuilt the current lineage and found:

```text
f2_scientific_payload_match = true
f3_scientific_payload_match = true
f4_scientific_payload_match = false
first_changed_stage = F4
```

The comparator-generated deterministic candidate byte identities are:

```text
candidate F.2 serialized SHA-256:
182433ccd2d13d3f7af80b963e69190ad8f2eb4e295c8216b36585e8ee7409b2

candidate F.3 serialized SHA-256:
eeb480d9b7ddffab1e97c7f9f27c3099d0780c2bbb4d42a0ee68c077a5752b3d
```

The current Left/lowe F.1 raw SHA-256 is:

```text
10086d7d4c42389c9fdd16471c49d59cd189a30980767bcdb87b85914169ef95
```

while the other four F.1 raw identities remain the accepted F.2 inputs.

The first F.4 scientific mismatch is:

```text
$.parents[0].absolute_baseline_parent_sum
accepted = 0.12294577672157099
candidate = 0.12072578661296547
```

The comparator global maxima are:

```text
baseline_parent_sum:
  absolute_difference = 0.0049304592644921486
  relative_difference = 0.044263728042183405

parent_normalization:
  absolute_difference = 0.06654204106342765
  relative_difference = 0.033455174261301984

correction_factor_summary:
  absolute_difference = 1.1986283985416364
  relative_difference = 0.03237215807179576
```

Only the three `Left-lowe` canonical-t parents differ numerically in the F.4
detail. `t=1` differs only at floating-point scale; `t=0` and `t=2` carry the
substantive baseline/current-normalization changes.

This result closes the **diagnostic question** posed by F.4.Refresh.1:

- F.2 accepted science remains numerically unchanged;
- F.3 accepted `hgcer3` map science remains numerically unchanged;
- F.4 is the first scientifically changed stage because its parent-preserving
  normalization consumes the current signed baseline contributions.

This task must **not** revise accepted source-owned authorities. It creates only
detached, explicitly non-authoritative **candidate artifacts** plus a
materialization manifest suitable for later validation/profile packaging.

---

## 2. Exact starting state

Required branch:

```text
test
```

Required committed starting HEAD:

```text
08f8c5be84ab7a54278e8c623eae5d9d3c43d938
```

Required commit subject:

```text
F4 Refresh.1: add current-baseline authority comparator
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
  `08f8c5be84ab7a54278e8c623eae5d9d3c43d938`;
- inspect the current worktree before editing;
- pre-existing temporary `kaonlt_review*.diff` files may remain untracked and
  must not be removed, rewritten, staged, or treated as candidate files;
- root `AGENTS.md` and `.codex/` remain local-only/untracked;
- if unrelated tracked changes exist, STOP and report the blocker;
- do not reset, stash, clean, discard, commit, push, run the farm, mutate
  accepted artifacts, or run the full analysis.

---

## 3. Mandatory repository-memory startup

Read in this exact order:

1. root `AGENTS.md`;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read only task-relevant records/source:

- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/phase-f6-method-a-production-promotion.md`
- `docs/memory/phases/f4-refresh1-current-baseline-authority-comparison.md`
- `docs/memory/evidence/e8-4-fix4-left-lowe-runtime-closure-and-f4-baseline-divergence.md`
- `docs/memory/evidence/f2-fix1-runtime-closure.md`
- `docs/memory/evidence/f3-runtime-closure.md`
- `docs/memory/evidence/f4-runtime-closure.md`
- `testing/compare_method_a_current_baseline_authority.py`
- `testing/test_compare_method_a_current_baseline_authority.py`
- `src/cuts/pion_hgcer_method_a_acceptance_representation.py`
- `src/cuts/pion_hgcer_method_a_acceptance_map.py`
- `src/cuts/pion_hgcer_method_a_parent_preserving_correction.py`
- `src/cuts/pion_hgcer_method_a_tphi_propagation.py`
- task-relevant existing F.2/F.3/F.4 tests.

Do not redesign adjacent architecture.

---

## 4. Status reconciliation owned by this task

Update durable memory to record the supplied F.4.Refresh.1 farm-comparator
evidence.

Advance:

```text
F.4.Refresh.1 — SOURCE REVIEWED
F.4.Refresh.1.Fix.1 — SOURCE REVIEWED
```

to:

```text
F.4.Refresh.1 — CLOSED / RUNTIME VALIDATED
F.4.Refresh.1.Fix.1 — CLOSED / RUNTIME VALIDATED
```

The closure is **only** for the detached current-baseline authority-comparison
gate. It proves:

```text
F.2 scientific payload unchanged
F.3 scientific payload unchanged
F.4 scientific payload changed
first changed stage = F4
```

It does not:

- replace accepted F.2/F.3/F.4 artifacts;
- reopen historical F.1-F.6.2 closures;
- validate a refreshed F.4 authority;
- validate F.5 against a refreshed F.4;
- validate F.6.3 or E.8.4 runtime;
- promote Method A.

Create:

```text
F.4.Refresh.2 — ACTIVE
```

pending independent ChatGPT actual-diff/source-runtime-path review of the
materializer candidate.

Preserve:

- historical F.1 through F.6.2 — **CLOSED / RUNTIME VALIDATED**
- E.8.4.Fix.4 — **CLOSED / RUNTIME VALIDATED** for its narrow cache-semantics
  gate only
- E.8 — **ACTIVE**
- E.8.2 — **SOURCE REVIEWED**
- E.8.3 — **SOURCE REVIEWED**
- F.6.3 — **SOURCE REVIEWED**
- E.8.4 — **SOURCE REVIEWED**
- final E.8 — **BLOCKED**
- F.6.4 — **BLOCKED**
- lifecycle-hook dispatch — **BLOCKED / DEFERRED**

The fresh F.6.3/E.8.4 Left/lowe runtime gate remains blocked until the
current-baseline authority lineage is resolved.

---

## 5. Scientific ownership and frozen boundaries

This task is detached validation infrastructure only.

Do not modify:

- random/dummy subtraction;
- slow-proton subtraction;
- pion component fits or definitions;
- dynamic alignment semantics;
- baseline pion-weight formula;
- F.1 event-contract producer;
- F.2 scientific representation calculation;
- F.3 `hgcer3` basis or map mathematics;
- F.4 parent-normalization mathematics;
- F.5 propagation mathematics;
- F.6.1/F.6.2 accepted science;
- F.6.3 branch mathematics;
- E.8.4 consumer;
- Method B;
- SIMC;
- cuts/windows/priors/templates/binning;
- yield extraction;
- efficiencies/acceptance;
- L/T separation;
- cross sections;
- active `no_empirical_residual` profile.

Do not update:

```text
ACCEPTED_F3_RUNTIME_AUTHORITY_BY_KINEMATIC
ACCEPTED_F4_RUNTIME_AUTHORITY_BY_KINEMATIC
```

or any equivalent source-owned acceptance constant.

Do not replace or overwrite accepted F.2/F.3/F.4 JSON files.

Do not promote Method A.

---

## 6. Allowed substantive implementation

Add exactly:

```text
testing/materialize_method_a_current_baseline_authority.py
testing/test_materialize_method_a_current_baseline_authority.py
```

Do not modify:

```text
testing/compare_method_a_current_baseline_authority.py
testing/test_compare_method_a_current_baseline_authority.py
```

No scientific/runtime source file may change.

The materializer may import and reuse the already-reviewed comparator's:

- canonical aliases;
- strict JSON reader semantics;
- deterministic writer-byte helper;
- scientific-projection exclusions;
- body validation;
- first-mismatch logic;
- F.4 detailed comparison logic.

It must invoke only the existing public F.2/F.3/F.4 artifact builders and public
JSON writers for candidate construction/materialization.

If existing scientific source must change to implement this task, STOP and
report the blocker.

---

## 7. Required CLI

Implement a pure-Python detached CLI with explicit inputs:

```text
--f1 Left-lowe=/path/current_f1.json
--f1 Left-highe=/path/current_f1.json
--f1 Center-lowe=/path/current_f1.json
--f1 Center-highe=/path/current_f1.json
--f1 Right-highe=/path/current_f1.json

--accepted-f2 /path/accepted_f2.json
--accepted-f3 /path/accepted_f3.json
--accepted-f4 /path/accepted_f4.json

--comparison /path/reviewed_f4_refresh1_comparison.json
--expected-comparison-sha256 <64-hex>
--source-head <40-hex>
--output-dir /path/existing/output_directory
```

Optional:

```text
--overwrite
```

Requirements:

- exact canonical-five aliases only;
- duplicate/missing/unexpected alias fails;
- every input must exist and be distinct;
- output directory must already exist and be a directory;
- no output may resolve to an input;
- no output may overwrite an existing file without `--overwrite`;
- `--expected-comparison-sha256` must be valid lowercase-normalized 64-hex;
- `--source-head` must be valid lowercase-normalized 40-hex and is recorded as
  supplied source provenance, not silently inferred;
- no ROOT/PyROOT import;
- scientific mismatch or malformed provenance fails closed before final output
  publication.

---

## 8. Deterministic output set

For kinematic token `{kinematic}`, derive exactly these distinct basenames:

```text
{kinematic}_kaon_pion-background_hgcer-method-a-current-baseline-authority-comparison-input.json

{kinematic}_kaon_pion-background_hgcer-method-a-acceptance-representation-current-baseline-candidate.json

{kinematic}_kaon_pion-background_hgcer-method-a-acceptance-map-current-baseline-candidate.json

{kinematic}_kaon_pion-background_hgcer-method-a-parent-preserving-correction-current-baseline-candidate.json

{kinematic}_kaon_pion-background_hgcer-method-a-current-baseline-authority-materialization-manifest.json
```

The first file must be a **byte-for-byte copy** of the validated comparator
input. Its SHA-256 must equal `--expected-comparison-sha256`.

The candidate F.2/F.3/F.4 artifacts must remain explicitly non-authoritative and
must use the existing artifact schemas.

The materialization manifest must be a new diagnostic schema:

```text
method_a_current_baseline_authority_materialization/v1
```

No output name may collide with accepted F.2/F.3/F.4 basenames.

---

## 9. Required preflight and scientific gates

All gates below must complete before final outputs are published.

### 9.1 Strict comparator identity

Read comparator raw bytes and require:

```text
SHA-256 == --expected-comparison-sha256
schema_version == method_a_current_baseline_authority_comparison/v1
non_authoritative == true
```

For the supplied F.4.Refresh.1 evidence, the expected SHA is:

```text
c16287e9192288ed5f116b5eebf95f49261d7fd09755ee923070154e047cf7b5
```

Do not hard-code that hash as a universal scientific constant. It is supplied
through the CLI and recorded in the manifest.

### 9.2 Accepted artifact identity continuity

Compute raw SHA-256 of the supplied accepted F.2/F.3/F.4 files and require exact
equality to the comparator's:

```text
accepted_file_sha256.f2
accepted_file_sha256.f3
accepted_file_sha256.f4
```

The current reviewed comparator records:

```text
f2 = 87e0a16870ea34a6611702e1bd4fc02c4dc6a8daa132e75b725a1b782a9a31da
f3 = 04a9a576767ba3f81c8c03daab2461998894c655f0d875bd532ad30cb46c0d95
f4 = adcc01900b8adc211e296aaedaa314fdb86207d51483ebfd4b3edf69f558f188
```

### 9.3 Current F.1 identity continuity

Require the current five F.1 raw SHA values to match the comparator's exact
`f1_inputs` inventory by setting ID.

Do not hard-code one setting independently from the comparison input.

Require one common kinematic token across all five F.1 artifacts and exact
agreement with the comparator F.1 inventory.

### 9.4 Rebuild deterministic candidate F.2

Invoke only:

```text
build_pion_hgcer_method_a_acceptance_representation_artifact(...)
```

with the same deterministic candidate-wrapper policy as F.4.Refresh.1:

```text
input_paths = {}
generated_at_utc = None
git_head = None
git_status_short = None
```

Compute writer-format bytes and require their SHA-256 to equal:

```text
comparison.candidate_serialized_sha256.f2
```

Then independently repeat the exact provenance-excluded scientific comparison
against accepted F.2 using the reviewed comparator projection.

Require:

```text
scientific_payload_match == true
first_mismatch == null
```

If not, fail closed.

### 9.5 Rebuild deterministic candidate F.3

Invoke only:

```text
build_pion_hgcer_method_a_acceptance_map_artifact(...)
```

using the current F.1 artifacts, candidate F.2 artifact, current F.1 raw hashes,
and candidate F.2 serialized SHA.

Use the same deterministic candidate-wrapper policy:

```text
input_paths = {}
generated_at_utc = None
git_head = None
git_status_short = None
```

Require writer-format candidate F.3 SHA-256 to equal:

```text
comparison.candidate_serialized_sha256.f3
```

Then independently repeat the exact reviewed F.3 scientific comparison against
accepted F.3.

Require:

```text
scientific_payload_match == true
first_mismatch == null
```

If not, fail closed.

### 9.6 Rebuild deterministic candidate F.4

Construct the candidate-F.3 runtime-authority override exactly as
F.4.Refresh.1 does:

- candidate F.3 serialized SHA;
- candidate F.3 map fingerprint;
- candidate F.3 algorithm fingerprint;
- candidate F.3 artifact fingerprint;
- `farm_source_head = "0" * 40`.

This is an explicit **diagnostic non-authority sentinel** only.

Invoke only:

```text
build_pion_hgcer_method_a_parent_preserving_correction_artifact(...)
```

with current F.1 artifacts and current F.1 raw hashes.

Use:

```text
input_paths = {}
generated_at_utc = None
git_head = None
git_status_short = None
```

Repeat the reviewed F.4 scientific comparison and detailed parent comparison.

Require exact reproduction of the supplied comparator's:

```text
f4.scientific_payload_match
f4.first_mismatch_path
f4.first_mismatch
f4.parents
f4.global_maxima
f4.setting_maxima
summary
```

For this accepted gate the required summary is:

```text
f2_scientific_payload_match = true
f3_scientific_payload_match = true
f4_scientific_payload_match = false
first_changed_stage = F4
```

Do not impose a numerical acceptance threshold on the F.4 change. The
materializer proves reproducibility and lineage only.

---

## 10. Candidate artifact byte requirements

Candidate F.2 and F.3 final written bytes must exactly match the deterministic
writer bytes used to satisfy the comparator SHA gate.

For F.2/F.3/F.4:

1. write via the module's existing public JSON writer into temporary files;
2. read the temporary raw bytes back;
3. verify those bytes match the deterministic expected writer bytes;
4. compute final candidate SHA-256 from those exact bytes.

The candidate F.2 written SHA must therefore remain:

```text
182433ccd2d13d3f7af80b963e69190ad8f2eb4e295c8216b36585e8ee7409b2
```

for the supplied current input set.

The candidate F.3 written SHA must therefore remain:

```text
eeb480d9b7ddffab1e97c7f9f27c3099d0780c2bbb4d42a0ee68c077a5752b3d
```

for the supplied current input set.

Do not hard-code either candidate SHA as a general rule; verify against the
comparison input.

Candidate F.4 SHA is newly computed and recorded in the materialization
manifest.

---

## 11. Materialization manifest

The manifest must contain at minimum:

```text
schema_version
non_authoritative = true
accepted_authority_mutated = false
production_objects_mutated = false
production_application_performed = false
method_a_promoted = false
source_head
kinematic_token
comparison_input
current_f1_inputs
accepted_inputs
candidate_outputs
scientific_gate
diagnostic_f3_authority_override
complete
errors
```

### `comparison_input`

Include:

- source basename;
- raw SHA-256;
- copied-output basename;
- copied-output raw SHA-256;
- schema version.

### `current_f1_inputs`

For each canonical setting include:

- setting ID;
- raw source-file SHA-256;
- stable F.1 content fingerprint;
- F.1 contract fingerprint;
- training/application record counts;
- training/application population fingerprints.

### `accepted_inputs`

Include accepted F.2/F.3/F.4:

- input basename;
- raw SHA-256;
- scientific fingerprint where available;
- artifact fingerprint.

### `candidate_outputs`

For candidate F.2/F.3/F.4 include:

- deterministic candidate basename;
- raw final SHA-256;
- artifact fingerprint;
- scientific payload fingerprint;
- explicit `non_authoritative = true`.

For F.3 also include:

- algorithm fingerprint;
- accepted basis.

For F.4 also include:

- F.3 source SHA inherited by the candidate;
- F.3 map/algorithm/artifact fingerprints;
- correction fingerprint.

### `scientific_gate`

Persist:

```text
f2_scientific_payload_match
f3_scientific_payload_match
f4_scientific_payload_match
first_changed_stage
f4_first_mismatch_path
f4_global_maxima
```

Do not copy all 9,000+ lines of detailed parent deltas into the manifest; the
byte-faithful comparison-input copy remains the detailed source.

### `diagnostic_f3_authority_override`

Persist the exact in-memory zero-head sentinel record and label it explicitly:

```text
accepted_authority = false
purpose = diagnostic_candidate_construction_only
```

---

## 12. Write safety and publication semantics

Build and validate **everything in memory first**.

Then:

1. create one temporary directory under the explicit output directory;
2. write the comparison-input copy and candidate artifacts only inside that
   temporary directory;
3. verify all raw bytes and hashes there;
4. build/write the manifest inside the temporary directory;
5. publish final candidate files into the output directory only after every
   scientific/provenance gate has passed;
6. publish the manifest **last**.

On any pre-publication failure:

- remove only the task-created temporary directory;
- publish no final output;
- never delete or rewrite inputs;
- never clean the repository or output directory.

Without `--overwrite`, any pre-existing target causes a fail-closed error before
temporary writes.

With `--overwrite`, replacement is permitted only for this exact derived five
file target set. No wildcard deletion.

The manifest's `complete = true` means only that all five detached candidate
outputs were written and verified. It is not scientific acceptance or
production promotion.

---

## 13. Forbidden shortcuts

Do not:

- write accepted artifact basenames;
- update source-owned F.3/F.4 runtime-authority constants;
- treat candidate F.2/F.3 byte hashes as accepted authority;
- treat the zero farm-head sentinel as an accepted farm source;
- weaken F.2/F.3 scientific equality to a tolerance;
- silently ignore unknown scientific fields;
- use a numerical threshold to approve the F.4 difference;
- skip the supplied comparison-input SHA gate;
- trust comparison summary without independently rebuilding F.2/F.3/F.4;
- use the old accepted Left/lowe F.1 artifact for candidate construction;
- mutate any current F.1 input;
- invoke F.5, F.6.1, F.6.2, F.6.3, E.8.4, Method B, ROOT, or `main.py`;
- modify production physics.

---

## 14. Required focused tests

Tests must use temporary directories and small synthetic fixtures with
monkeypatched public builders/writers where appropriate.

At minimum test:

1. exact canonical-five F.1 inventory required;
2. malformed/duplicate/unexpected aliases rejected;
3. invalid comparator SHA rejected before output;
4. accepted F.2/F.3/F.4 raw-SHA mismatch rejected;
5. current F.1 raw-SHA mismatch against comparison input rejected;
6. candidate F.2 writer SHA must reproduce comparator candidate SHA;
7. candidate F.3 writer SHA must reproduce comparator candidate SHA;
8. independent F.2 scientific mismatch rejects materialization;
9. independent F.3 scientific mismatch rejects materialization;
10. F.4 detail/summary mismatch against comparison input rejects
    materialization;
11. F.3 authority override uses the candidate identities and all-zero farm-head
    sentinel only;
12. public F.2/F.3/F.4 builders are each called once in dependency order;
13. public F.2/F.3/F.4 writers are used for candidate artifact bytes;
14. comparison-input final copy is byte-identical to the validated input;
15. deterministic output basenames cannot collide with accepted basenames;
16. existing outputs reject without `--overwrite`;
17. a pre-publication failure leaves no final outputs;
18. materialization manifest is published last and records all required hashes,
    flags, and scientific-gate fields;
19. no input is modified;
20. no ROOT/PyROOT import.

Preserve and rerun the existing comparator suite.

---

## 15. Required local deterministic checks

Run:

```bash
python -B -m py_compile \
  testing/materialize_method_a_current_baseline_authority.py \
  testing/test_materialize_method_a_current_baseline_authority.py

python -B -m unittest \
  testing.test_materialize_method_a_current_baseline_authority -v

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

If actual existing test module names differ, discover and use the exact
repository names rather than inventing substitutes.

Report skips exactly.

These checks do not establish farm/runtime validation.

---

## 16. Warranted durable memory/history

Allowed existing memory updates:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/f4-refresh1-current-baseline-authority-comparison.md
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/roadmap/STATUS.md
```

Create:

```text
docs/memory/evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md
docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization-task-contract.md
docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization.md
```

The new F.4.Refresh.1 evidence record must distinguish direct evidence from
interpretation.

### Direct evidence to record

- comparison filename;
- comparison SHA-256
  `c16287e9192288ed5f116b5eebf95f49261d7fd09755ee923070154e047cf7b5`;
- accepted F.2/F.3/F.4 raw SHA values;
- candidate F.2/F.3 serialized SHA values;
- exact comparator summary:
  - F.2 match true;
  - F.3 match true;
  - F.4 match false;
  - first changed stage F4;
- first F.4 mismatch path and values;
- global maxima listed in Section 1;
- current Left/lowe F.1 raw SHA;
- only three Left/lowe parents carry nonzero F.4 differences in the detailed
  comparison; `t=1` is floating-point scale while `t=0` and `t=2` carry the
  substantive current-baseline changes.

### Interpretation to record separately

- accepted F.2 science does not require scientific refresh;
- accepted F.3 map science does not require scientific refresh;
- current-baseline provenance lineage must still be rebuilt because F.3/F.4
  fail-closed input identities bind to current F.1 artifacts;
- F.4 is the first scientifically changed stage;
- accepted historical F.4 remains valid for its accepted historical baseline;
- a current-baseline F.4 authority candidate must be materialized and reviewed
  before any downstream source-owned authority re-pin.

Do not claim F.4.Refresh.2 runtime acceptance.

---

## 17. CURRENT.md exact NEXT after implementation

Set exactly one ordinary NEXT:

```text
NEXT — independent ChatGPT actual-diff/source-runtime-path review of the F.4.Refresh.2 current-baseline Method-A candidate materializer.
```

Do not add another NEXT.

---

## 18. Required diff audit

Before stopping:

```bash
git status --short
git -c core.safecrlf=false diff --stat
git -c core.safecrlf=false diff --no-ext-diff
```

Audit every intended new/untracked file with:

```bash
git diff --no-index /dev/null <path> || true
```

Expected substantive code paths:

```text
testing/materialize_method_a_current_baseline_authority.py
testing/test_materialize_method_a_current_baseline_authority.py
```

plus only the contract-allowed memory/history paths.

No scientific/runtime source may change.

---

## 19. Required cumulative review bundle

Create one fresh temporary repository-root review bundle:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

It must contain byte-faithfully:

1. branch;
2. committed HEAD;
3. `git status --short`;
4. complete cumulative `git diff --stat`;
5. complete tracked `git -c core.safecrlf=false diff --no-ext-diff`;
6. complete `git diff --no-index /dev/null ...` for every intended untracked
   candidate file;
7. exact local validation commands/results;
8. final changed-path inventory.

The review bundle is temporary and must remain untracked.

---

## 20. Acceptance criteria

The local candidate is acceptable for independent source review only if:

1. committed HEAD remains exactly
   `08f8c5be84ab7a54278e8c623eae5d9d3c43d938`;
2. no scientific/runtime source file changes;
3. F.4.Refresh.1 direct farm-comparator evidence is durably recorded;
4. F.4.Refresh.1 and Fix.1 advance only to
   **CLOSED / RUNTIME VALIDATED** for their detached comparator gate;
5. historical F.1-F.6.2 closures remain unchanged;
6. F.6.3/E.8.4 remain **SOURCE REVIEWED** and runtime blocked;
7. materializer consumes exact canonical-five current F.1 inputs;
8. accepted F.2/F.3/F.4 raw identities are checked against the reviewed
   comparison input;
9. comparison raw SHA is explicitly checked;
10. candidate F.2 and F.3 are rebuilt through public builders;
11. candidate F.2/F.3 serialized SHA values reproduce the reviewed comparator;
12. candidate F.2/F.3 scientific payloads are independently proved identical
    to accepted F.2/F.3 before publication;
13. candidate F.4 is rebuilt through the public F.4 builder with only the
    explicit in-memory zero-head candidate-F.3 override;
14. candidate F.4 detailed comparison and summary exactly reproduce the reviewed
    F.4.Refresh.1 comparison;
15. accepted authority constants remain untouched;
16. candidate artifacts use distinct non-authoritative filenames;
17. public writers generate candidate artifact bytes;
18. final outputs are published only after all gates pass and manifest is
    published last;
19. no input or accepted artifact is modified;
20. local deterministic tests pass subject to explicitly reported skips;
21. no farm command, full analysis, ROOT/PyROOT, commit, push, authority update,
    Method-A promotion, or downstream re-pin occurs;
22. one fresh cumulative review bundle is produced.

---

## 21. Hard stop

After implementation, deterministic local checks, warranted memory/history
updates, manifest regeneration, diff audit, and creation of the fresh
timestamped cumulative review bundle:

**STOP.**

Do not commit, push, create the F.4.Refresh.2 farm validation profile, run the
materializer on the farm, package candidate artifacts, update accepted
F.2/F.3/F.4 authority, modify F.5/F.6.3 authority, rerun the full analysis, or
promote Method A.

NEXT is independent ChatGPT actual-diff/source-runtime-path review of the
F.4.Refresh.2 materializer candidate.
