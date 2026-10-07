# E.8.4 / F.6.3 canonical-five provenance/identity repair — task contract

## 1. Task class and starting authority

This is a **narrow source-changing provenance/identity repair task** under
`docs/memory/CODEX.md`.

Authoritative repository:

`https://github.com/trottar/lt_analysis/tree/test`

Exact pushed starting identity:

```text
branch: test
HEAD: ab3f29c9268a80cd872903003da0422195d118e8
origin/test: ab3f29c9268a80cd872903003da0422195d118e8
```

The user has explicitly chosen to **reopen** the canonical-five
provenance/identity repair that CURRENT previously marked `DEFERRED`.

This task does **not** reopen the scientific design of Method A. It owns only
the provenance/lineage boundary that prevents the current canonical-five
runtime from consuming scientifically equivalent freshly regenerated F.1
lineage.

Do not run the Jefferson Lab farm.
Do not commit, push, update refs, stage merely for review, reset, clean, stash,
or overwrite unrelated local work.

---

## 2. Required startup and hard start gate

Before editing, read in this exact order:

1. `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read:

- `docs/memory/MAINTENANCE.md`
- `docs/memory/CODEX.md`
- `docs/memory/COMMUNICATION.md`
- `docs/memory/TOOLS.md`
- `docs/memory/LEARNINGS.md`
- `docs/memory/evidence/e8-4-canonical-five-second-full-run-provenance-blocker-2026-10-05.md`
- `docs/memory/evidence/e8-4-canonical-five-fresh-f1-lineage-materialization-2026-10-05.md`
- `docs/memory/evidence/e8-4-canonical-five-lineage-preflight-runtime-closure-2026-10-05.md`
- `docs/memory/phases/e8-4-method-a-current-lineage-canonical-five-validation-orchestration-task-contract.md`
- `docs/memory/evidence/e8-3-left-lowe-runtime-closure-2026-10-07.md`
- `docs/memory/roadmap/STATUS.md`

Then inspect the task-relevant source and tests listed below before changing
anything.

Report:

```bash
git branch --show-current
git rev-parse HEAD
git rev-parse origin/test
git status --short --untracked-files=all
```

Required branch/HEAD/origin are exactly the identity above.

The worktree may contain unrelated pre-existing local state. Preserve it
byte-for-byte. Hard stop if its ownership is ambiguous or this task cannot be
completed without resetting, cleaning, stashing, checking out over, deleting,
or staging unrelated paths.

---

## 3. Established blocker — do not rediscover from scratch

The second full canonical-five farm attempt at source

`ace8688a71431d13b40ed19713a27746f3da6a8e`

completed the child analysis with return code zero but failed the owner artifact
gate. All five page manifests contained only
`full_background.e8_4.unavailable`, with the same reason:

```text
f6_3_f4_shared_reproduction_failed:f3_fingerprint_input_content_mismatch
```

The supplied post-run diagnosis established:

```text
f2_scientific_payload_match = true
f3_scientific_payload_match = true
f4_scientific_payload_match = true
first_changed_stage = none
```

while all five regenerated F.1 raw SHA-256 and stable F.1 content fingerprints
differed from the pre-analysis candidate lineage.

This is accepted diagnostic runtime evidence. It does **not** establish final
canonical-five closure.

### Source-level failure point already audited

Current
`src/cuts/pion_hgcer_method_a_parent_preserving_correction.py`
strictly validates the persisted F.3 artifact against the **current** F.1
lineage.

In `_validate_f3_artifact(...)`, it reconstructs:

```text
expected_input_content =
  setting_id
  stable_f1_content_fingerprint
  source_file_sha256
```

from the current F.1 artifacts and rejects when the stored F.3
`fingerprint_inputs["input_content"]` differs:

```text
f3_fingerprint_input_content_mismatch
```

It then separately requires each stored F.3 input fingerprint to match the
current F.1 setting, raw SHA, stable fingerprint and F.1 contract fingerprint.

This strict check is correct for ordinary persisted F.3 lineage validation.
The blocker arises because F.6.3 currently presents a staged candidate F.3/F.4
from one F.1 provenance lineage to that validator after the full run has
regenerated a different F.1 provenance lineage.

### Comparator behavior already audited

`testing/compare_method_a_current_baseline_authority.py` rebuilds current
F.2/F.3/F.4 from the supplied F.1 artifacts through the public scientific
builders and compares **scientific projections** after excluding explicitly
identified provenance/fingerprint fields. The accepted post-run invocation
reported exact scientific equality through F.4.

The repair must consume that scientific fact safely. It must not simply suppress
the F.3 input-content check.

---

## 4. Objective

Implement the smallest fail-closed provenance/identity bridge that allows the
F.6.3 canonical-five runtime path to proceed when, and **only when**, a newly
regenerated current F.1 lineage is deterministically proven scientifically
equivalent to the reviewed staged current-baseline F.2/F.3/F.4 candidate
lineage.

The desired behavior is:

```text
current F.1 lineage exact-match candidate lineage
    -> existing strict path succeeds unchanged

otherwise
    -> reconstruct current F.2/F.3/F.4 using the existing public scientific
       builders and frozen algorithms
    -> compare the reconstructed scientific payload to the staged reviewed
       candidate F.2/F.3/F.4 scientific payload under an explicit, source-owned,
       deterministic equivalence contract
    -> PASS only if every required scientific field matches exactly under that
       contract
    -> expose/consume current-lineage factors only after the equivalence gate
       passes
    -> record both candidate authority identity and current regenerated
       provenance identity
```

Any scientific mismatch must remain a hard failure.

Do not make the stored candidate appear to have been generated from the new
F.1 lineage. Preserve both provenance roles explicitly.

---

## 5. Design constraints

### 5.1 Preserve the strict F.4 persisted-artifact validator by default

Preferred architecture:

- do **not** weaken `_validate_f3_artifact(...)` globally;
- preserve its exact ordinary behavior for callers validating a persisted F.3
  against its own F.1 lineage;
- solve the mismatch at the F.6.3 candidate/current-lineage boundary.

If the source audit proves a change inside
`pion_hgcer_method_a_parent_preserving_correction.py` is unavoidable, the
default strict behavior must remain unchanged and any new equivalence path must
be explicit, opt-in, narrow to F.6.3 candidate/current-lineage reconstruction,
and covered by negative tests.

A global “ignore raw hashes”, “ignore stable fingerprints”, or
“accept matching map fingerprint” change is forbidden.

### 5.2 Reconstruct rather than relabel

Do not mutate or rewrite the staged candidate F.3/F.4 provenance fields to
pretend they came from the current F.1 artifacts.

Do not refresh a static F.1 SHA pin to the latest run. The runtime evidence
already demonstrates that run-to-run F.1 provenance identities can change.

Do not bypass validation by copying candidate fingerprints into regenerated
objects.

### 5.3 Scientific equivalence must be explicit and complete

A valid equivalence gate must be based on the actual current F.1 artifacts and
the existing public F.2/F.3/F.4 builders.

At minimum it must prove:

- exact canonical-five setting inventory;
- current F.1 artifacts pass their existing full contracts;
- current reconstructed F.2 scientific payload equals the staged candidate F.2
  scientific payload under one explicit projection;
- current reconstructed F.3 scientific payload equals the staged candidate F.3
  scientific payload under one explicit projection;
- current reconstructed F.4 scientific payload equals the staged candidate F.4
  scientific payload under one explicit projection;
- exact 15-parent F.4 inventory;
- exact parent correction factors/normalizations and existing parent-closure
  semantics;
- no Method-B numerical dependency;
- no production promotion;
- no child renormalization;
- no change to thresholds/tolerances/scientific algorithms.

The scientific projection may remove only fields proven to be provenance or
serialized identity metadata. It must be centrally defined, named, documented
and tested. Do not use an open-ended recursive “ignore anything different”
rule.

If a field’s classification as provenance versus science is ambiguous, stop
`BLOCKED` rather than excluding it.

### 5.4 Runtime consumption after equivalence

After a successful equivalence gate, F.6.3 must consume factors that are
scientifically tied to the **current regenerated F.1 population**.

Preferred behavior is to consume the transient current-lineage reconstruction
whose scientific payload has just been proven equal to the reviewed candidate,
rather than forcing the stale-provenance candidate through the strict validator.

If a different design is necessary, document why it is equally fail-closed and
preserves current-population identity.

### 5.5 Provenance record

The F.6.3 provenance returned downstream must clearly distinguish:

- reviewed candidate authority:
  - candidate F.2/F.3/F.4 file identities/fingerprints;
- current runtime lineage:
  - current five F.1 raw SHA-256 values;
  - current stable F.1 content fingerprints;
  - transient current F.2/F.3/F.4 identities/fingerprints as applicable;
- equivalence decision:
  - schema/version;
  - exact fields/projections compared;
  - pass/fail;
  - first mismatch path on failure when deterministically available.

Never collapse candidate authority and current runtime provenance into one
identity.

---

## 6. Scientific source frozen

Do not change the scientific mathematics or accepted definitions in:

- F.1 population construction;
- F.2 candidate definitions or representation algorithm;
- F.3 accepted basis, features, fitting, support/OOD logic or model parameters;
- F.4 correction mathematics, signed parent preservation, support or closure
  tolerance;
- F.5 propagation;
- F.6.2 validation;
- baseline pion subtraction;
- random/dummy/slow-proton treatment;
- SIMC;
- yields, efficiencies, acceptance, L/T separation or cross sections;
- `no_empirical_residual`.

Method A remains detached/non-production.
Method B remains diagnostic/cross-check only and numerically excluded.

No production promotion follows from this task.

---

## 7. Mandatory source audit before implementation

Inspect at minimum:

```text
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
src/cuts/pion_hgcer_method_a_parent_preserving_correction.py
src/cuts/pion_hgcer_method_a_acceptance_map.py
src/cuts/pion_hgcer_method_a_acceptance_representation.py

testing/compare_method_a_current_baseline_authority.py
testing/test_compare_method_a_current_baseline_authority.py
testing/test_f6_3_parallel_full_procedure_method_a.py
testing/test_pion_hgcer_method_a_parent_preserving_correction.py
testing/run_e8_4_fix5_canonical_five_plot_gate.py
testing/test_run_e8_4_fix5_canonical_five_plot_gate.py
```

Trace:

```text
fresh F.1 generation
-> F.6.3 path discovery/load
-> staged candidate F.3/F.4 load
-> reconstruct_transient_factor_map()
-> shared F.4 builder
-> _validate_f3_artifact()
-> exact mismatch
-> downstream E.8.4 availability
```

Also trace the existing comparator’s scientific projection and verify every
excluded field is genuinely provenance/identity metadata for the active
candidate.

Before changing code, state in the task notes:

1. the exact current mismatch path;
2. why current post-run scientific equality can coexist with F.1/F.3 provenance
   inequality;
3. the narrowest source boundary that can own the repair;
4. why the proposed repair cannot silently accept a genuinely changed F.2/F.3/F.4
   result.

Hard stop if that argument cannot be made from source and accepted evidence.

---

## 8. Allowed versioned paths

The task may change only the following **if source audit demonstrates each is
needed**:

```text
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
src/cuts/pion_hgcer_method_a_parent_preserving_correction.py

testing/compare_method_a_current_baseline_authority.py
testing/test_compare_method_a_current_baseline_authority.py
testing/test_f6_3_parallel_full_procedure_method_a.py
testing/test_pion_hgcer_method_a_parent_preserving_correction.py
testing/test_run_e8_4_fix5_canonical_five_plot_gate.py

docs/memory/phases/e8-4-canonical-five-provenance-identity-repair-task-contract.md
docs/memory/phases/e8-4-canonical-five-provenance-identity-repair.md
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/roadmap/STATUS.md
docs/memory/manifest.json
```

Prefer fewer paths.

`testing/run_e8_4_fix5_canonical_five_plot_gate.py` is **frozen by default**.
If audit proves the tracked owner itself must change to expose or enforce the
new equivalence gate, stop and report `BLOCKED` before editing it; ChatGPT will
decide whether to broaden the contract.

No other versioned path may change.

In particular, do not modify:

```text
run_Prod_Analysis.sh
set_SymLinks.sh
src/cuts/pion_hgcer_method_a_acceptance_map.py
src/cuts/pion_hgcer_method_a_acceptance_representation.py
src/cuts/pion_hgcer_method_a_tphi_propagation.py
src/cuts/full_background_subtraction_plots.py
testing/collect_pion_hgcer_validation_bundle.py
testing/pion_hgcer_validation_bundle_profile_e8_4_canonical_five.json
farm_env/
background_samples/
```

The F.2/F.3 builders are audit references and must remain unchanged.

---

## 9. Required positive tests

Add focused deterministic tests proving at least:

### Exact-lineage path

- an exact candidate/current F.1 lineage still uses the existing strict path;
- results are unchanged from the pre-repair path.

### Equivalent-new-lineage path

Construct synthetic or existing fixtures where:

- current F.1 raw identity differs;
- current stable F.1 identity differs where the fixture supports it;
- transient F.2/F.3/F.4 rebuilt from current F.1 has the same defined scientific
  payload as the reviewed candidate;
- the equivalence gate passes;
- the current-lineage factors are returned;
- candidate and current provenance remain separately recorded.

### Fail-closed scientific differences

Independently perturb representative scientific content at each stage and
require failure:

- F.2 scientific field;
- F.3 model/scaler/coefficient/support field;
- F.4 parent correction/normalization field;
- parent inventory;
- setting inventory;
- Method-B flag;
- child-renormalization or production-promotion flag.

A change that only modifies an explicitly allowlisted provenance field may be
accepted only if the reconstructed scientific payload still matches exactly.

### No accidental weakening

Direct tests of the ordinary F.4 persisted-artifact validator must continue to
reject:

- mismatched F.1 raw SHA;
- mismatched stable F.1 content fingerprint;
- mismatched F.1 contract fingerprint;
- malformed F.3 fingerprint inputs;
- altered F.3 model content.

If the F.4 module is not changed, retain these existing regressions and add only
what is necessary to prove the new F.6.3 path does not bypass them globally.

---

## 10. Regression tests

Run at minimum:

```text
testing/test_f6_3_parallel_full_procedure_method_a.py
testing/test_pion_hgcer_method_a_parent_preserving_correction.py
testing/test_compare_method_a_current_baseline_authority.py
testing/test_run_e8_4_fix5_canonical_five_plot_gate.py
testing/test_e8_4_production_impact_audit.py
testing/test_e8_4_fix5_main_order.py
testing/test_collect_pion_hgcer_validation_bundle.py
```

Also run any directly affected F.2/F.3 tests even though those sources are
frozen.

Do not run the real analysis, ROOT/PyROOT, or the farm.

---

## 11. Memory/status update after successful local implementation

If the source audit and deterministic repair both pass, update CURRENT so that:

- canonical-five provenance/identity repair is `ACTIVE` during implementation;
- after successful local implementation/source checks it is
  `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`;
- the prior second full runtime gate remains `BLOCKED` historical evidence;
- E.8.2 Left/lowe remains `CLOSED / RUNTIME VALIDATED`;
- E.8.3 Left/lowe remains `CLOSED / RUNTIME VALIDATED`;
- F.6.3/E.8.4 existing Left/lowe closures retain their exact scopes;
- final E.8 remains `BLOCKED` pending new canonical-five runtime evidence;
- F.6.4 remains `BLOCKED`;
- Method A remains detached/non-production;
- Method B remains diagnostic/cross-check only and numerically excluded;
- absolute-SIMC interpretation remains separately `BLOCKED`.

Push-stable sole NEXT must be:

```text
NEXT — after ChatGPT actual-diff review, user commit/push, and pushed-state
synchronization, perform the farm-readiness audit for one fresh isolated
Q4p4W2p74 canonical-five owner gate using the reviewed provenance/identity
repair. Only Farm readiness: PASS may authorize that run.
```

Do not mark the canonical-five runtime itself validated locally.

Record a concise phase note explaining:

- exact low-level mismatch;
- repair boundary;
- scientific-equivalence contract;
- provenance separation;
- deterministic checks;
- farm boundary.

Update durable MEMORY only with a general lesson that exact serialized lineage
identity and scientific-equivalence identity are distinct and must never be
interchanged silently; any bridge requires deterministic reconstruction and a
fail-closed explicit scientific projection.

Update roadmap only as needed to reflect the reopened repair and its source
state. CURRENT remains the only exact NEXT owner.

---

## 12. Memory health and deterministic completion

Discover `<PYTHON>` according to `docs/memory/TOOLS.md`.

After all allowed changes:

```bash
<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/check_memory_health.py --root .
<PYTHON> -B tools/memory_bootstrap.py --root . --json
git diff --check
```

Report:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
warning classification: <none or each warning: blocking | nonblocking>
manifest check: PASS | FAIL
bootstrap: PASS | FAIL
git diff --check: PASS | FAIL
```

Local checks are `SOURCE VERIFIED` only. They do not establish ROOT/PyROOT,
full `main.py`, procedure-PDF, farm, canonical-five runtime, or production
validation.

---

## 13. Review bundle

Create one temporary repository-root review file:

```text
kaonlt_review.diff
```

It must contain:

- complete diffs for every modified tracked file;
- complete `git diff --no-index /dev/null ...` additions for every intended new
  file;
- no unrelated paths.

Do not stage merely to make the review bundle.

Stop before commit/push.

ChatGPT will review the actual diff.

---

## 14. Acceptance criteria

Return:

`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`

only if all of the following hold:

- exact start identity passed;
- unrelated local state was preserved;
- the low-level mismatch was traced from current source;
- no scientific algorithm, threshold, basis, cut, correction formula or parent
  normalization was changed;
- ordinary strict F.3/F.4 lineage validation remains fail-closed;
- current-lineage reconstruction uses existing public scientific builders;
- the equivalence projection is explicit and closed, not heuristic;
- exact scientific equality through F.2/F.3/F.4 is required before use;
- current and candidate provenance remain separately visible;
- any scientific mismatch fails before factors can be consumed;
- all focused and regression tests pass;
- memory/manifest/bootstrap/diff checks pass;
- only allowlisted paths changed;
- complete review bundle exists;
- no ROOT/farm/stage/commit/push/ref update occurred.

---

## 15. Hard stops

Return `BLOCKED` without inventing a workaround if:

- branch/HEAD/origin differ from the starting identity;
- unrelated local state cannot be preserved;
- the exact runtime mismatch cannot be explained from current source;
- safe repair requires changing F.1/F.2/F.3 scientific algorithms;
- safe repair requires changing F.4 correction mathematics or tolerances;
- safe repair requires owner changes outside this contract;
- scientific equality requires a tolerance or field exclusion not already
  source-justified;
- the repair would trust only raw hashes, only one high-level fingerprint, or
  only a candidate pin refresh;
- the repair would relabel stale candidate provenance as current provenance;
- the repair imports test code into runtime source;
- any focused negative test demonstrates a scientifically changed F.2/F.3/F.4
  payload can pass;
- a path outside the allowlist is required;
- memory/manifest hard checks fail.

Do not rerun the farm from this task.
