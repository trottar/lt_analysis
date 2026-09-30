# KaonLT E.8.4.Fix.4 — Persisted Alignment Semantic-Version Gate

## 1. Purpose

Repair one newly demonstrated runtime/provenance defect in the persisted pion-component dynamic-alignment cache before any further E.8.4 / F.6.3 farm validation.

The fresh `Q4p4W2p74 / Left / lowe` Fix.3 farm gate completed the ordinary low-epsilon analysis but E.8.4 failed closed while F.6.3 attempted to reproduce the accepted F.4 authority:

```text
f6_3_f4_shared_reproduction_failed:f3_fingerprint_input_content_mismatch
```

Read-only farm comparison established:

- accepted F.3 file SHA-256 still exactly:
  `04a9a576767ba3f81c8c03daab2461998894c655f0d875bd532ad30cb46c0d95`;
- accepted F.4 file SHA-256 still exactly:
  `adcc01900b8adc211e296aaedaa314fdb86207d51483ebfd4b3edf69f558f188`;
- four canonical settings still reproduce the accepted F.1 raw/stable identity;
- only `Left / lowe` differs;
- frozen `Left / lowe` F.1 SHA-256:
  `f16c5e89f62e848ec221335bd4916265757040af9970b005d8f67cc800e0077d`;
- fresh Fix.3 `Left / lowe` F.1 SHA-256:
  `10086d7d4c42389c9fdd16471c49d59cd189a30980767bcdb87b85914169ef95`;
- frozen stable F.1 fingerprint:
  `76bcf74ac1de6c19d402e078bc8b8ac7f1c0f821e6d7fda1ef97338541d7c988`;
- fresh stable F.1 fingerprint:
  `250d95f3f210e3bca7315439777a71adbdb1c21c5c74bfeab36490e432f7dfee`.

The F.1 partition comparison further established that the following remain unchanged:

```text
part1_config_fingerprint
method_a_event_population_fingerprint
acceptance_feature_metadata_fingerprint
application_child_assignment_projection_fingerprint
```

while these change for `Left / lowe`:

```text
phase_a_contract_fingerprint
phase_a_pion_event_population_fingerprint
coordinate_fingerprint
method_a_fingerprint
method_a_training_population_fingerprint
application_population_fingerprint
```

The frozen accepted `Left / lowe` coordinate fingerprint is:

```text
207456612e2f4b88daa8de7d89b55e0a56d52e01096608fb3ef7d9316d19f2b4
```

The fresh Fix.3 `Left / lowe` coordinate fingerprint is:

```text
e41d4d004e047076c213cab7bedf5a51fe56f6f7b83c3c9e0fd5e2a2e52a4066
```

That fresh fingerprint is the same coordinate fingerprint recorded by the September 29 pre-Fix.3 alignment-blocker artifact.

The source audit identifies the missing cache boundary:

1. E.8.4.Fix.3 changed the **alignment candidate-acceptance semantics** by making the configured `minimum_template_integral` comparison machine-scale floating-point safe.
2. Fix.3 deliberately retained alignment schema v2 and the same resolved scientific configuration.
3. `load_or_resolve_pion_component_alignment(...)` builds only cache metadata, attempts to load a compatible persisted record, and **returns the cache hit before running `resolve_pion_component_alignment(...)`**.
4. `_alignment_compatibility_reasons(...)` currently validates schema/config/scope/physical identity/axis/parent/template identity and stable pion-control checksum+axis, but contains **no implementation/semantic revision for the candidate-scan acceptance algorithm**.
5. Therefore a pre-Fix.3 schema-v2 record produced under the old strict floating-point boundary can still be declared compatible after Fix.3. On such a cache hit, the repaired minimum-integral predicate is never evaluated.

This new direct farm evidence invalidates the earlier Fix.3 assumption that all valid schema-v2 stores could remain reusable across the changed scan semantics.

This task adds one explicit persisted **alignment-resolution semantic-version gate** so records produced under pre-Fix.3 scan semantics fail closed and are deterministically recomputed by the already-existing resolver. It does not change the scientific configuration, scan grid, scoring, templates, windows, component ordering, subtraction formulas, Method A, Method B, or accepted F.1--F.6.2 artifacts.

---

## 2. Exact starting state

Required branch:

```text
test
```

Required committed starting HEAD:

```text
8d62dcfdf2fd8d08298c075940b4ab1b28d7f079
```

Commit subject:

```text
E8.4 Fix.3: re-pin validation bundle profile
```

The Fix.4 implementation is based on source commit ancestry containing the pushed Fix.3 analysis source:

```text
29d7b7f9635db899939efeb3508e941e994e8928
```

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
```

Hard requirements:

- branch must be `test`;
- committed HEAD must be exactly `8d62dcfdf2fd8d08298c075940b4ab1b28d7f079`;
- inspect the existing worktree before editing;
- do not reset, stash, clean, discard, overwrite, or stage unrelated user work;
- root `AGENTS.md`, `.codex/`, the task contract, and temporary `kaonlt_review*.diff` files may be local/untracked and must not be removed;
- if unrelated tracked edits already exist, STOP and report the blocker;
- do not commit, push, run the Jefferson Lab farm, package a farm bundle, or mutate farm artifacts.

---

## 3. Mandatory repository-memory startup

Read in this exact order:

1. root `AGENTS.md`;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read only task-relevant records/source, including at minimum:

- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/e8-4-fix3-pion-alignment-determinism.md`
- `docs/memory/phases/e8-4-fix3-pion-alignment-determinism-task-contract.md`
- `docs/memory/phases/phase-f6-method-a-production-promotion.md`
- the current E.8.4 production-impact phase record
- directly referenced Fix.3 blocker evidence
- `src/cuts/pion_component_fits.py`
- `testing/test_pion_component_dynamic_alignment.py`
- `src/cuts/pion_hgcer_method_a_acceptance_contract.py`
- `src/cuts/pion_hgcer_method_a_acceptance_map.py`
- `src/cuts/pion_hgcer_method_a_parent_preserving_correction.py`
- `src/cuts/pion_hgcer_method_a_parallel_full_procedure.py`

Read additional source only when required to trace the existing runtime path.

---

## 4. Status and scientific ownership

Preserve:

- F.1 through F.6.2, including F.6.2.Fix.5 — **CLOSED / RUNTIME VALIDATED**
- E.8.2 — **SOURCE REVIEWED**
- E.8.3 — **SOURCE REVIEWED**
- F.6.3 — **SOURCE REVIEWED**
- E.8.4 — **SOURCE REVIEWED**
- E.8 — **ACTIVE**
- final E.8 — **BLOCKED**
- F.6.4 — **BLOCKED**
- lifecycle-hook dispatch — **BLOCKED / DEFERRED**

Preserve E.8.4.Fix.3 as **SOURCE REVIEWED**. Its implementation/source review is not erased by this task, but its earlier assumption that all valid schema-v2 alignment stores remained reusable has been invalidated by the new runtime evidence.

Create:

- E.8.4.Fix.4 persisted-alignment semantic-version repair — **ACTIVE** pending independent ChatGPT actual-diff/source-runtime-path review.

Do not reopen F.1--F.6.2 science. The accepted F.3/F.4 files and frozen F.6.2 scientific JSON remain unchanged authorities until a later explicit gate says otherwise.

The fresh `Left / lowe` E.8.4 runtime gate remains **BLOCKED** pending this repair, source review, user-controlled push, pushed-state review, any required validation-profile re-pin, and a new narrow farm run.

---

## 5. Source-audit facts that define this repair

### 5.1 Cache load occurs before the scan

At the starting source, `load_or_resolve_pion_component_alignment(...)`:

1. calls `build_expected_pion_alignment_metadata(...)`;
2. calls `load_pion_component_alignment(...)`;
3. returns immediately when a compatible persisted record is found;
4. calls `resolve_pion_component_alignment(...)` only after a cache miss/rejection.

Therefore an old cache hit bypasses the Fix.3 scan predicate entirely.

### 5.2 Existing compatibility lacks algorithm semantics

At the starting source, `_alignment_compatibility_reasons(...)` checks:

- alignment schema version;
- resolved configuration hash;
- analysis scope;
- complete physical-bin identity;
- histogram axis;
- parent alignment hash when required;
- immutable source-template identifier/checksum;
- stable pion-control content checksum plus axis.

It does **not** distinguish records produced under:

- pre-Fix.3 strict minimum-template-integral comparison; versus
- Fix.3 floating-point-safe minimum-template-integral comparison.

### 5.3 Fix.3 changed executable acceptance semantics without changing config

Fix.3 correctly preserves the configured scientific threshold:

```text
minimum_template_integral = 1.0
```

and evaluates the boundary with the production helper equivalent to:

```python
value < threshold and not math.isclose(
    value,
    threshold,
    rel_tol=1e-12,
    abs_tol=1e-12,
)
```

This is a source-level algorithm semantic change, not a user configuration change. Consequently `resolved_configuration_hash` cannot distinguish pre-/post-Fix.3 records.

### 5.4 Downstream fail-closed behavior is correct

F.3 deliberately records broad frozen F.1 input identity.

F.4 verifies the current F.1 raw SHA and stable F.1 fingerprint against the accepted F.3 input identity and raises:

```text
f3_fingerprint_input_content_mismatch
```

when they differ.

F.6.3 deliberately reconstructs accepted F.4 from current F.1 + accepted F.3 and requires exact persisted F.4 reproduction before any Method-A factor is exposed.

Do not weaken these gates.

---

## 6. Allowed substantive source changes

Only:

- `src/cuts/pion_component_fits.py`
- `testing/test_pion_component_dynamic_alignment.py`

No other substantive production/test file may change.

In particular, do **not** modify:

- `src/utility/background_config.py`
- `src/cuts/pion_hgcer_method_a_acceptance_contract.py`
- `src/cuts/pion_hgcer_method_a_acceptance_map.py`
- `src/cuts/pion_hgcer_method_a_parent_preserving_correction.py`
- `src/cuts/pion_hgcer_method_a_tphi_propagation.py`
- `src/cuts/pion_hgcer_method_a_parallel_full_procedure.py`
- `src/cuts/full_background_subtraction_plots.py`
- `src/cuts/rand_sub.py`
- `src/cuts/pion_component_subtraction.py`
- `src/cuts/particle_subtraction.py`
- `src/binning/calculate_yield.py`
- `src/binning/ave_per_bin.py`
- `src/main.py`
- `run_Prod_Analysis.sh`
- collector/wrapper/profile files
- accepted scientific artifacts.

If implementation requires another substantive source file, STOP and report the blocker.

---

## 7. Required implementation

### 7.1 Add one explicit alignment-resolution semantics version

Define a source-owned, deterministic constant in `src/cuts/pion_component_fits.py`, with a clear name such as:

```python
PION_COMPONENT_ALIGNMENT_SEMANTICS_VERSION = (
    "pion_component_dynamic_alignment_semantics/v2"
)
```

Equivalent naming is acceptable, but the purpose must be explicit: this is the semantic revision of the executable alignment resolver/cache result, distinct from user scientific configuration and distinct from the existing serialized schema version.

Do not derive it from timestamps, Git HEAD, environment, hostname, ROOT name, or runtime state.

### 7.2 Persist the semantics version in every alignment result/metadata path

The same current semantics version must be present in:

- `build_expected_pion_alignment_metadata(...)`;
- normal successful `resolve_pion_component_alignment(...)` output;
- disabled/fallback `resolve_pion_component_alignment(...)` output;
- persisted JSON records;
- CSV diagnostic rows, where practical for direct provenance inspection.

No scientific quantity may depend numerically on the text value.

### 7.3 Make it an exact fail-closed cache compatibility boundary

`_alignment_compatibility_reasons(stored, expected)` must require exact semantics-version equality.

Required behavior:

- record has current semantics version + every existing compatibility field matches => eligible for reuse;
- record lacks semantics version => stale/incompatible;
- record has different semantics version => stale/incompatible;
- malformed/non-string semantics value => stale/incompatible.

Use a clear rejection reason such as:

```text
alignment_semantics_version mismatch
```

Do not silently map a missing field to the current version.

The existing old schema-v2 record must therefore be rejected and the ordinary resolver must run once, producing a new current-semantics record.

### 7.4 Preserve schema v2

Do **not** change the existing alignment schema version in `background_config.py`.

This repair intentionally separates:

- serialized payload schema; from
- executable alignment-resolution semantics.

The new evidence shows that schema-v2 compatibility alone is insufficient, but it does not require redefining the payload schema.

### 7.5 Parent-to-fine boundary must also be semantics-safe

The direct `parent_valid` check in `resolve_pion_component_alignment(...)` must require that a parent alignment carries the current semantics version.

A pre-Fix.3 parent passed directly to a fine-bin resolver must not be accepted merely because schema/config/parent hash match.

Do not invent a fallback that strips or substitutes the version.

Changing `_alignment_parent_hash(...)` is not required merely to force invalidation; prefer the explicit semantics boundary. If Codex finds a concrete correctness reason that the parent hash itself must include this version, STOP and report that reasoning before expanding behavior.

### 7.6 Preserve Fix.3 semantic pion-control identity

Do not regress the Fix.3 behavior:

- same pion-control checksum + same complete axis + different generated ROOT histogram name => compatible;
- checksum mismatch => incompatible;
- axis mismatch => incompatible;
- malformed semantic identity => incompatible/fail closed.

The new semantics-version gate is additive to those boundaries.

### 7.7 Preserve Fix.3 floating-point boundary

Do not change:

```python
_template_integral_meets_minimum(...)
```

except for comments/tests if necessary.

Its current scientific threshold and `1e-12` machine-scale comparison remain frozen for this task.

### 7.8 No forced physics result

Do not hard-code or force:

- `pi_n = -0.003` or `-0.004 GeV`;
- `pi_delta = +0.008` or `+0.016 GeV`;
- any specific parent score;
- any specific `w0`;
- any coordinate fingerprint;
- a `Left / lowe` special case;
- cache deletion;
- unconditional rescan on every run.

A current-semantics compatible cache must still be reusable normally.

---

## 8. Before / after behavior

### Before

A persisted schema-v2 alignment record produced under the pre-Fix.3 strict integral boundary can satisfy all current compatibility fields. `load_or_resolve_pion_component_alignment(...)` then returns it before the corrected candidate scan runs.

This permits stale pre-Fix.3 alignment semantics to propagate into Phase-A/F.1 provenance, after which F.6.3 correctly fails against accepted F.3/F.4 authority.

### After

The same pre-Fix.3 record lacks the current alignment semantics version and is rejected as stale.

The existing resolver then executes with:

- the same scientific configuration;
- the same component order;
- the same templates;
- the same scan grids/windows;
- the Fix.3 semantic histogram identity;
- the Fix.3 machine-scale minimum-integral boundary.

The newly persisted record carries the current semantics version. A later run with the same current semantics and all existing provenance equal may reuse it.

No farm file needs to be manually deleted to obtain the corrected behavior.

---

## 9. Required deterministic tests

Extend only:

```text
testing/test_pion_component_dynamic_alignment.py
```

### 9.1 Current-semantics reuse

Retain and strengthen the existing test proving:

- same checksum;
- same axis;
- different generated histogram name;
- same current semantics version;

=> cache is reused and resolver is not called.

### 9.2 Pre-Fix/missing-semantics cache rejection

Construct a valid persisted cache record equivalent to the current metadata, then remove the semantics-version field to represent a pre-Fix record.

Call the actual `load_or_resolve_pion_component_alignment(...)` path.

The test must prove:

- the old record is not returned as a cache hit;
- rejection contains the semantics-version mismatch reason;
- the resolver is invoked exactly once.

It is acceptable to monkeypatch only `resolve_pion_component_alignment` with a narrow deterministic sentinel payload so the cache-control path is executable without PyROOT. Do not duplicate scan physics in the test.

### 9.3 Explicit wrong-version rejection

A record with a non-current semantics version must fail compatibility.

### 9.4 Malformed semantics rejection

At minimum cover missing and non-string/malformed values. Fail closed.

### 9.5 Parent semantics rejection

A fine-bin resolver must not treat a parent lacking/currently mismatching semantics version as valid even when its schema/config/parent hash otherwise match.

Use the existing focused parent helper/tests; do not alter the scientific parent/fine algorithm.

### 9.6 Preserve Fix.3 tests

Keep passing coverage for:

- renamed histogram with same checksum+axis;
- checksum mismatch;
- axis mismatch;
- malformed pion-control semantic identity;
- exact `minimum_template_integral` threshold;
- machine-scale below threshold;
- materially below threshold.

### 9.7 Regression boundaries

Preserve existing coverage for:

- config/source/axis/physical-identity cache mismatch;
- immutable raw-template single-shift provenance;
- parent-to-fine fallback;
- independent boundary policies;
- support metric;
- mixed accepted/fallback component behavior;
- no kaon-side data entering alignment calibration;
- F.6.3 and E.8.4 fail-closed contracts.

---

## 10. Frozen scientific/runtime interfaces

Do not change:

- random/dummy subtraction;
- slow-proton subtraction;
- baseline pion subtraction formula;
- component definitions/order;
- common-setting shift;
- scan min/max/step;
- fit/evaluation windows;
- window expansions;
- interpolation;
- template renormalization;
- scientific acceptance thresholds;
- scoring/ranking;
- component amplitudes;
- canonical t/phi binning;
- baseline yields;
- Method-A training/application definitions;
- F.3 map mathematics;
- F.4 parent normalization;
- F.5/F.6.1/F.6.2 accepted results;
- F.6.3 private branch mathematics;
- E.8.4 consumer;
- Method B;
- SIMC normalization;
- efficiencies;
- cross sections;
- L/T separation;
- active `no_empirical_residual` profile.

Do not change production physics to make downstream authority checks pass.

---

## 11. Local deterministic validation

Use the repository's actual local Python interpreter.

Run at minimum:

```bash
python -B -m py_compile \
  src/cuts/pion_component_fits.py \
  testing/test_pion_component_dynamic_alignment.py

python -B -m unittest testing.test_pion_component_dynamic_alignment -v

python -B -m unittest testing.test_t_bin_pion_parent_integrity -v
python -B -m unittest testing.test_binning_pre_particle_subtraction -v

python -B -m unittest testing.test_f6_3_parallel_full_procedure_method_a -v
python -B -m unittest testing.test_e8_4_production_impact_audit -v

python -B tools/update_memory_manifest.py --root . --write
python -B tools/update_memory_manifest.py --root . --check
python -B tools/check_memory_health.py --root .
python -B tools/memory_bootstrap.py --root . --json
python -B -m unittest testing.test_memory_health -v

git -c core.safecrlf=false diff --check
```

If PyROOT-dependent tests skip, report the exact skip count/reason. Do not claim ROOT/PyROOT or farm validation.

Do not run the farm.

---

## 12. Warranted durable memory/history

Update only what this task establishes.

Allowed existing memory/history files:

- `docs/memory/CURRENT.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/manifest.json`
- `docs/memory/phases/e8-4-fix3-pion-alignment-determinism.md`

Allowed new records:

- `docs/memory/phases/e8-4-fix4-alignment-cache-semantics-task-contract.md`
- `docs/memory/phases/e8-4-fix4-alignment-cache-semantics.md`
- `docs/memory/evidence/e8-4-fix3-left-lowe-stale-alignment-cache-runtime-blocker.md`

The new evidence record should capture only direct established evidence:

- fresh Left/lowe analysis completed;
- E.8.4 unavailable reason is
  `f6_3_f4_shared_reproduction_failed:f3_fingerprint_input_content_mismatch`;
- accepted F.3/F.4 SHA values still match;
- only Left/lowe F.1 differs across canonical five;
- exact frozen/current Left/lowe F.1 raw SHA and stable fingerprints;
- unchanged versus changed partition-fingerprint inventory listed in Section 1;
- frozen/current coordinate fingerprints;
- fresh coordinate fingerprint equals the September 29 blocker artifact's coordinate fingerprint;
- source audit proves cache load precedes scan and cache compatibility lacks scan-semantics identity.

Clearly distinguish the last item as source proof and the farm values as runtime evidence.

Do not claim the fresh farm run directly printed `persistence_status=reused`; it did not. The defect is established by the combination of current source behavior and the supplied artifact/fingerprint evidence.

`CURRENT.md` after implementation must make Fix.4 **ACTIVE pending independent ChatGPT actual-diff review** and set the exact NEXT to that review. Do not authorize a farm rerun yet.

Regenerate `docs/memory/manifest.json` whenever versioned memory changes.

Do not modify `MEMORY.md`, `LEARNINGS.md`, `USER.md`, `CODEX.md`, `TOOLS.md`, decisions, or unrelated phase/evidence files unless a concrete blocker requires it; if so STOP before expanding scope.

---

## 13. Required diff audit

Before stopping:

```bash
git status --short
git diff --stat
git diff -- \
  src/cuts/pion_component_fits.py \
  testing/test_pion_component_dynamic_alignment.py \
  docs/memory/CURRENT.md \
  docs/memory/roadmap/STATUS.md \
  docs/memory/manifest.json \
  docs/memory/phases/e8-4-fix3-pion-alignment-determinism.md
```

Audit all intended new files with:

```bash
git diff --no-index /dev/null docs/memory/phases/e8-4-fix4-alignment-cache-semantics-task-contract.md || true
git diff --no-index /dev/null docs/memory/phases/e8-4-fix4-alignment-cache-semantics.md || true
git diff --no-index /dev/null docs/memory/evidence/e8-4-fix3-left-lowe-stale-alignment-cache-runtime-blocker.md || true
```

No unrelated source may change.

---

## 14. Required review bundle

Create one fresh temporary cumulative review bundle at repository root:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

It must contain, byte-faithfully:

1. branch;
2. committed HEAD;
3. `git status --short`;
4. complete `git diff --stat`;
5. complete tracked `git -c core.safecrlf=false diff --no-ext-diff`;
6. complete `git diff --no-index /dev/null ...` output for every intended untracked candidate file;
7. exact local test commands and results;
8. final changed-path inventory.

Do not Markdown-normalize, sanitize, replace quotes, or rewrite the raw diff text.

The review bundle is temporary and must not be staged or committed.

Stop for independent ChatGPT review after producing it.

---

## 15. Acceptance criteria

The implementation candidate is acceptable for independent source review only if all are true:

1. committed HEAD remained exactly `8d62dcfdf2fd8d08298c075940b4ab1b28d7f079`;
2. substantive edits are limited to the two allowed files;
3. a source-owned alignment semantics version exists;
4. current expected metadata and resolved payloads carry it;
5. persisted pre-Fix/missing-version cache records fail closed;
6. current-semantics compatible records still reuse;
7. a histogram-name-only change still reuses when checksum+axis match;
8. current checksum/axis/malformed identity checks remain;
9. the Fix.3 minimum-integral predicate is unchanged;
10. parent-to-fine validity requires current semantics;
11. no scan grid/window/config/score/physics logic changed;
12. F.3/F.4/F.6.3/E.8.4 source remains untouched;
13. accepted F.1--F.6.2 artifacts remain untouched;
14. focused/local regression tests pass subject to explicitly reported PyROOT skips;
15. memory records distinguish source proof from farm proof;
16. Fix.4 remains **ACTIVE** pending independent ChatGPT review;
17. no commit, push, farm run, bundle packaging, or production promotion occurs.

---

## 16. Forbidden shortcuts

Do not:

- delete or rename the farm cache as the repair;
- force an unconditional rescan forever;
- bump the scientific configuration merely to invalidate the cache;
- change `minimum_template_integral`;
- increase the Fix.3 floating tolerance;
- hard-code the observed Left/lowe shifts or fingerprints;
- weaken F.3/F.4/F.6.3 provenance checks;
- refresh accepted F.3/F.4 authority in this task;
- regenerate F.1--F.6.2 accepted science;
- special-case a kinematic/setting;
- alter Method A or Method B mathematics;
- alter baseline production physics;
- change the validation profile in this task;
- run the Jefferson Lab farm;
- commit or push.

---

## 17. Hard stop

After implementation, deterministic local checks, warranted memory updates, manifest regeneration, diff audit, and creation of the fresh cumulative review bundle:

**STOP.**

NEXT is independent ChatGPT review of the actual cumulative candidate diff/source-runtime path.

No profile re-pin, Git handoff, or farm command is authorized by this task.
