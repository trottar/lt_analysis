# E.8.4 / F.6.3 canonical-five reviewed-F.2 handoff and scientific-equivalence repair — task contract

## 1. Task class and starting authority

This is a **narrow source-changing provenance/runtime-wiring repair task** under
`docs/memory/CODEX.md`.

Authoritative repository:

`https://github.com/trottar/lt_analysis/tree/test`

Exact pushed starting identity:

```text
branch: test
HEAD: 5c8736c20e94e2f07776cfccc47bbbfd790fd057
origin/test: 5c8736c20e94e2f07776cfccc47bbbfd790fd057
```

This contract follows the accepted blocked audit in:

```text
docs/memory/phases/e8-4-canonical-five-provenance-identity-repair.md
```

The audit established that the previous repair contract could not be implemented
safely because the canonical-five owner verifies a reviewed F.2/F.3/F.4
materialization but stages only F.3/F.4, while the F.6.3 runtime path discovers
only fresh F.1 plus staged F.3/F.4.

ChatGPT actual-diff review explicitly approved one narrow follow-up:

> hand off the already reviewed/hash-pinned candidate F.2 alongside F.3/F.4,
> then implement a fail-closed current-lineage F.2/F.3/F.4 scientific-equivalence
> gate before current-lineage factors may be consumed.

This task owns exactly that follow-up.

It does **not** reopen Method-A science, F.1/F.2/F.3/F.4 algorithms, parent
normalization, baseline pion subtraction, or production physics.

Do not run the Jefferson Lab farm.
Do not stage, commit, push, update refs, reset, clean, stash, or overwrite
unrelated local work.

---

## 2. Required startup and hard start gate

Before editing, read in exact order:

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
- `docs/memory/phases/e8-4-canonical-five-provenance-identity-repair.md`
- `docs/memory/phases/e8-4-canonical-five-provenance-identity-repair-task-contract.md`
- `docs/memory/evidence/e8-4-canonical-five-second-full-run-provenance-blocker-2026-10-05.md`
- `docs/memory/evidence/e8-4-canonical-five-fresh-f1-lineage-materialization-2026-10-05.md`
- `docs/memory/evidence/e8-4-canonical-five-lineage-preflight-runtime-closure-2026-10-05.md`
- `docs/memory/phases/e8-4-method-a-current-lineage-canonical-five-validation-orchestration-task-contract.md`
- `docs/memory/roadmap/STATUS.md`

Then inspect all source/test paths listed in Sections 6 and 7.

Report:

```bash
git branch --show-current
git rev-parse HEAD
git rev-parse origin/test
git status --short --untracked-files=all
```

Required:

```text
branch = test
HEAD = 5c8736c20e94e2f07776cfccc47bbbfd790fd057
origin/test = 5c8736c20e94e2f07776cfccc47bbbfd790fd057
```

Unrelated pre-existing local state may exist. Preserve it byte-for-byte.

Hard stop if:

- branch/HEAD/origin differ;
- unrelated local state cannot be preserved;
- implementation requires reset/clean/stash/checkout/delete of unrelated work;
- a scientific algorithm or frozen production path must change;
- a path outside this contract's allowlist is required.

---

## 3. Established failure and exact repair boundary

The second canonical-five full runtime attempt completed the child analysis but
failed downstream with:

```text
f6_3_f4_shared_reproduction_failed:f3_fingerprint_input_content_mismatch
```

Accepted diagnostic evidence established that the post-run reconstructed
scientific payload matched through all three stages:

```text
f2_scientific_payload_match = true
f3_scientific_payload_match = true
f4_scientific_payload_match = true
first_changed_stage = none
```

while the freshly regenerated F.1 raw and stable identities differed from the
pre-analysis candidate lineage.

The strict F.4 persisted-F.3 validator is correct and must remain strict. It
rejects when a persisted F.3 claims an F.1 provenance lineage different from
the current F.1 lineage.

The repair therefore belongs at the **private F.6.3 candidate/current-lineage
boundary**, not inside the general scientific validators.

### SOURCE VERIFIED current wiring

Current owner:

```text
testing/run_e8_4_fix5_canonical_five_plot_gate.py
```

already verifies all materialized candidate files:

```text
F.2 raw SHA-256 = 2fa715b2b5c2d5077e38416fb9eceb814e013e72f308c43526fddff36ddd962e
F.3 raw SHA-256 = c5b86452b790ecbaf5b8f0df05da67efa2fa92aab157b12b153ed7e491a38228
F.4 raw SHA-256 = 79e7ceda7221cbeeead4ed5bc306b0e0e670741a27beaa980e22349c555e96d7
```

with reviewed fingerprints:

```text
F.2 representation:
e9df6e09c31dda4d7d724792ce3c40125e27530c4c76e0dbe7f5b9861d0e6216

F.2 artifact:
87d9a8780cb2de9f2830afd9e388064151f796654ddb5c4a8e1938443806a0c6

F.3 map:
6042e9485827d9884d2a841d5c6cb5e23401e1026c8e2d4bd2022ede15843548

F.3 algorithm:
ba29630b2f40a87cbadbe751504ce48f23e2b17a08378a2c8a219f93131cb912

F.3 artifact:
e1cf0dd14644034f5d21df27d29f48355cde7c8e85c6123b709a6975c152a302

F.4 correction:
71527eecd5828766ebcb3b9ef5e4e1931eb242fabfc1059f04e4eec104f6cc98

F.4 artifact:
0b3899e5a3f24a34e04aa550ef162414c1d95af9a786e748612cbcb018dbcad4
```

But current owner `CANDIDATES` stages only F.3/F.4.

Current:

```text
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
```

discovers/loads fresh F.1 plus staged F.3/F.4 only.

Current:

```text
src/binning/calculate_yield.py
```

passes only F.1/F.3/F.4 into `reconstruct_transient_factor_map()`.

The repair must make reviewed F.2 available through this path and use it only
for the explicit scientific-equivalence gate described below.

---

## 4. Required after behavior

### 4.1 Owner handoff

The canonical-five owner must:

- continue to verify the complete reviewed materialization;
- stage **F.2, F.3 and F.4**, not only F.3/F.4;
- require exact reviewed raw SHA-256 for all three;
- never synthesize or rewrite candidate provenance;
- record all three staged identities in owner status/summary;
- preserve all existing isolation, preservation, freshness, artifact, page,
  collection and ZIP gates;
- leave the ordinary checkout and installed ltsep untouched as before.

The F.2 candidate is:

```text
Q4p4W2p74_kaon_pion-background_hgcer-method-a-acceptance-representation-current-baseline-candidate.json
```

with raw SHA-256:

```text
2fa715b2b5c2d5077e38416fb9eceb814e013e72f308c43526fddff36ddd962e
```

### 4.2 Runtime authority discovery

The private F.6.3 source path must discover/load:

```text
fresh current F.1 canonical five
reviewed candidate F.2
reviewed candidate F.3
reviewed candidate F.4
```

No fallback to historical F.2/F.3/F.4 is allowed.

No directory search, newest-file heuristic, or guessed path is allowed.

### 4.3 Exact-lineage path

When the current F.1 raw lineage exactly equals the reviewed candidate F.1
lineage, preserve the existing strict reconstruction behavior.

The reviewed candidate F.2 identity must still be verified, but the scientific
result and transient factors must remain unchanged from the current exact-lineage
path.

### 4.4 New-lineage scientific-equivalence path

When current F.1 identity differs from the reviewed candidate F.1 lineage:

1. Validate all current F.1 artifacts with the existing public contracts.
2. Rebuild transient current-lineage F.2 from current F.1 using the existing
   public F.2 builder and frozen configuration.
3. Compare current F.2 to reviewed candidate F.2 using the explicit scientific
   projection in Section 5.
4. Only if F.2 passes, rebuild current-lineage F.3 from current F.1 + rebuilt
   current F.2 using the existing public F.3 builder.
5. Compare current F.3 to reviewed candidate F.3 under the explicit projection.
6. Only if F.3 passes, rebuild current-lineage F.4 from current F.1 + rebuilt
   current F.3 using the existing public F.4 builder with a temporary
   current-lineage F.3 authority record.
7. Compare current F.4 to reviewed candidate F.4 under the explicit projection.
8. Require exact 15-parent inventory and unchanged parent-closure semantics.
9. Only after all three stages pass may F.6.3 consume **current-lineage**
   transient event factors from the rebuilt F.4 review data.

Do not force the stale-lineage reviewed F.3 through the strict validator against
new F.1.

Do not consume stale-lineage event factors.

### 4.5 Fail closed

Any scientific mismatch at F.2, F.3 or F.4 must make the private branch
unavailable before factor application.

The failure must identify:

```text
stage
first mismatch path
```

where deterministically available.

No fallback to candidate factors, historical authorities, Method B, or baseline
production is permitted as a substitute for the private Method-A branch.

Baseline production itself remains unchanged and available under its existing
rules.

---

## 5. Scientific-equivalence contract

### 5.1 Single source-owned projection

Move or define the exact scientific-projection semantics in source-owned code,
preferably:

```text
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
```

The testing comparator may consume that source-owned definition.

Do **not** import `testing/` code into runtime source.

The projection must be closed and explicit. Unknown fields remain compared.

The accepted exclusions are only the source-audited provenance/serialized
identity fields already used by the tracked comparator:

```text
F.2:
  input_fingerprints
  fingerprint_inputs
  fingerprint

F.3:
  input_fingerprints
  fingerprint_inputs
  fingerprint
  f2_representation_fingerprint
  f2_source_file_sha256

F.4:
  input_fingerprints
  fingerprint_inputs
  fingerprint
  f3_source_file_sha256
  f3_map_fingerprint
  f3_artifact_fingerprint
  f3_runtime_authority
```

Do not expand these exclusion sets without returning `BLOCKED`.

Retain and compare all other fields, including:

- F.2 candidate definitions/configuration/support/groups/summaries/recommendation;
- F.3 accepted basis, ordered features, algorithm configuration and fingerprint,
  models, scalers, coefficients, support/OOD fields and continuity content;
- F.4 accepted basis, algorithm configuration, all 15 parents and all parent
  scientific metrics/diagnostics;
- all unknown/new fields not explicitly excluded.

Comparisons are exact recursive equality. Do not introduce a numerical tolerance
for the equivalence projection.

Existing scientific builders may retain their already-defined internal
tolerances; this repair must not change them.

### 5.2 Candidate wrapper/identity validation

Before scientific projection, the runtime path must fail closed unless reviewed
candidate wrappers have their expected:

- schema versions;
- non-authoritative/validation flags;
- no-production/no-Method-B flags;
- raw SHA-256;
- primary scientific fingerprint;
- artifact fingerprint.

The owner hash gate is not a reason to omit runtime fail-closed identity checks.

### 5.3 Current-lineage provenance

The returned F.6.3 authority/provenance must clearly separate:

```text
reviewed_candidate
  F.2 raw SHA + representation/artifact fingerprints
  F.3 raw SHA + map/algorithm/artifact fingerprints
  F.4 raw SHA + correction/artifact fingerprints
  reviewed candidate validation/materialization source identity

current_runtime_lineage
  five current F.1 raw SHA-256
  five current stable F.1 content fingerprints
  rebuilt F.2 serialized/scientific fingerprint identities
  rebuilt F.3 serialized/scientific fingerprint identities
  rebuilt F.4 serialized/scientific fingerprint identities

scientific_equivalence
  schema/version
  exact-lineage | equivalent-new-lineage mode
  F.2 pass/fail + first mismatch path
  F.3 pass/fail + first mismatch path
  F.4 pass/fail + first mismatch path
  all_stages_passed
```

Do not relabel reviewed candidate provenance as current runtime provenance.

Do not persist event-level correction-factor arrays in the public authority
record if the current design intentionally keeps them transient.

---

## 6. Allowed source/test paths

Only the following source/test paths may change, and only as required by this
contract:

```text
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
src/binning/calculate_yield.py

testing/run_e8_4_fix5_canonical_five_plot_gate.py
testing/test_run_e8_4_fix5_canonical_five_plot_gate.py

testing/pion_hgcer_validation_bundle_profile_e8_4_canonical_five.json
testing/test_pion_hgcer_validation_bundle_profile_e8_4_canonical_five.py

testing/compare_method_a_current_baseline_authority.py
testing/test_compare_method_a_current_baseline_authority.py

testing/test_f6_3_parallel_full_procedure_method_a.py
```

Prefer fewer paths where possible.

If another source/test path is required, return `BLOCKED`.

---

## 7. Explicitly frozen scientific/runtime source

Do not modify:

```text
src/cuts/pion_hgcer_method_a_acceptance_contract.py
src/cuts/pion_hgcer_method_a_acceptance_representation.py
src/cuts/pion_hgcer_method_a_acceptance_map.py
src/cuts/pion_hgcer_method_a_parent_preserving_correction.py
src/cuts/pion_hgcer_method_a_tphi_propagation.py
src/cuts/full_background_subtraction_plots.py

run_Prod_Analysis.sh
set_SymLinks.sh
farm_env/
background_samples/

testing/collect_pion_hgcer_validation_bundle.py
```

Also freeze all unrelated `src/` and `testing/` paths.

In `src/binning/calculate_yield.py`, only the private F.6.3 authority
load/pass-through wiring may change. Do not alter baseline yield extraction,
histogram arithmetic, uncertainty propagation, cuts, SIMC, or final public
yield behavior.

---

## 8. Canonical-five bundle/profile behavior

Because reviewed F.2 becomes a real canonical-five runtime input, the dedicated
canonical-five validation profile must collect it as a required global artifact.

Expected global artifact keys after this repair:

```text
candidate_f2
candidate_f3
candidate_f4
run_summary
```

Per-setting artifact inventory remains unchanged.

The generic collector remains unchanged.

Update profile tests for the new exact inventory and total artifact count.

Owner ZIP verification must require the same four global keys.

Do not mutate the historical Left/lowe profile.

---

## 9. Required focused tests

### 9.1 Owner handoff tests

Prove:

- materialization verification still checks F.2/F.3/F.4 exact hashes;
- staging now installs exactly reviewed F.2/F.3/F.4;
- wrong F.2 raw hash fails before analysis;
- unknown existing candidate-F.2 target fails closed;
- staged F.2 is not rewritten;
- owner summary records F.2/F.3/F.4 candidate identities;
- final artifact/ZIP inventory contains required candidate F.2;
- existing ordinary-checkout and ltsep preservation behavior is unchanged.

Replace the old regression that explicitly required “stages only F.3/F.4” with
the new exact three-candidate contract.

### 9.2 Exact-lineage F.6.3 tests

Prove:

- exact reviewed F.1 lineage + reviewed F.2/F.3/F.4 succeeds;
- exact-lineage factors and parent closure match the pre-repair behavior;
- strict persisted F.3/F.4 validators are still exercised;
- candidate and current provenance fields are not conflated.

### 9.3 Equivalent-new-lineage tests

Using deterministic synthetic fixtures, create current F.1 identity changes that
produce an exactly equal projected F.2/F.3/F.4 scientific result.

Require:

- F.2 is rebuilt from current F.1;
- F.3 is rebuilt from current F.1 + rebuilt current F.2;
- F.4 is rebuilt from current F.1 + rebuilt current F.3;
- all three exact scientific projections pass;
- returned event factors come from current-lineage F.4 review data;
- current F.1 raw/stable identities appear under current provenance;
- reviewed candidate identities remain under candidate provenance;
- no candidate event-factor array is substituted.

### 9.4 Scientific mismatch negative tests

Independently perturb and require failure before factor use for:

- F.2 candidate definition/config/support/group/summary/recommendation content;
- F.3 model coefficient;
- F.3 scaler;
- F.3 support/OOD content;
- F.3 accepted basis/feature ordering/algorithm content;
- F.4 parent normalization;
- F.4 correction-factor summary;
- F.4 source diagnostics;
- F.4 canonical-phi diagnostics;
- missing/extra parent;
- wrong setting inventory;
- Method-B numerical-dependency flag;
- production-promotion flag;
- child-renormalization flag.

A difference only in one explicitly excluded provenance field may pass **only**
when the full reconstructed scientific projection remains exactly equal and all
wrapper/identity gates pass.

### 9.5 No weakening of default validators

Existing tests must continue to prove ordinary strict validation rejects:

- F.1 raw SHA mismatch;
- stable F.1 fingerprint mismatch;
- F.1 contract fingerprint mismatch;
- malformed F.3 fingerprint input content;
- altered F.3 model science;
- F.4 authority mismatch.

No change to the default strict validator is expected.

---

## 10. Deterministic test-fixture clock issue

The blocked audit recorded a local `/tmp` filesystem mtime discrepancy that made
some owner fixture freshness tests fail even though the owner behavior was
unchanged.

Do not weaken production freshness semantics:

```text
artifact st_mtime_ns >= started_ns
```

If that local discrepancy recurs, this task may stabilize **tests only** by
mocking the owner start-time sample or explicitly setting synthetic fixture
mtime under test control.

The production owner freshness rule must remain byte/semantics unchanged except
for unrelated F.2 inventory wiring.

Do not use filesystem-clock behavior as justification to skip owner tests.

---

## 11. Regression tests

Run at minimum:

```text
testing/test_f6_3_parallel_full_procedure_method_a.py
testing/test_compare_method_a_current_baseline_authority.py
testing/test_run_e8_4_fix5_canonical_five_plot_gate.py
testing/test_pion_hgcer_validation_bundle_profile_e8_4_canonical_five.py

testing/test_pion_hgcer_method_a_acceptance_representation.py
testing/test_pion_hgcer_method_a_acceptance_map.py
testing/test_pion_hgcer_method_a_parent_preserving_correction.py
testing/test_pion_hgcer_method_a_tphi_propagation.py

testing/test_e8_4_production_impact_audit.py
testing/test_e8_4_fix5_main_order.py
testing/test_collect_pion_hgcer_validation_bundle.py
```

Also run any existing test module that directly loads
`src/binning/calculate_yield.py` if not already covered above.

No test may run the real farm analysis.

---

## 12. Memory/status paths allowed

After successful implementation and deterministic checks, only these memory
paths may change:

```text
docs/memory/phases/e8-4-canonical-five-f2-handoff-equivalence-repair-task-contract.md
docs/memory/phases/e8-4-canonical-five-f2-handoff-equivalence-repair.md
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/roadmap/STATUS.md
docs/memory/manifest.json
```

Do not edit historical evidence records.

### CURRENT after successful local implementation

Record:

- the previous audit remains the accepted explanation of the missing F.2 handoff;
- reviewed candidate F.2 is now handed off with F.3/F.4;
- private F.6.3 has a fail-closed exact scientific-equivalence path for
  regenerated F.1 lineage;
- ordinary strict F.3/F.4 validation remains unchanged;
- candidate and current provenance remain separate;
- local deterministic checks passed;
- canonical-five runtime is **not** validated locally;
- repair state becomes:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

Keep:

```text
E.8.2 Left/lowe: CLOSED / RUNTIME VALIDATED
E.8.3 Left/lowe: CLOSED / RUNTIME VALIDATED
existing F.6.3/E.8.4 Left/lowe scopes: unchanged
canonical-five second full runtime attempt: BLOCKED historical evidence
final E.8: BLOCKED
F.6.4: BLOCKED
absolute-SIMC interpretation: BLOCKED
Method A: detached/non-production
Method B: diagnostic/cross-check only, numerically excluded
```

Push-stable sole NEXT:

```text
NEXT — after ChatGPT actual-diff review, user commit/push, and pushed-state
synchronization, perform the farm-readiness audit for one fresh isolated
Q4p4W2p74 canonical-five owner gate using the reviewed F.2/F.3/F.4
provenance-equivalence repair. Only Farm readiness: PASS may authorize that run.
```

### MEMORY

Add only the durable rule:

- exact serialized lineage identity and scientific-equivalence identity are
  distinct;
- a permitted bridge must reconstruct the current lineage through the same
  public scientific builders, compare an explicit closed scientific projection
  exactly, preserve both provenance roles, and consume only current-lineage
  transient factors.

Do not add mutable HEAD/NEXT detail to MEMORY.

---

## 13. Local health and integrity checks

Discover `<PYTHON>` using repository convention.

After all changes:

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
warning classification: <none or blocking/nonblocking>
manifest check: PASS | FAIL
bootstrap: PASS | FAIL
git diff --check: PASS | FAIL
```

Local tests are `SOURCE VERIFIED` only. They do not establish farm, ROOT/PyROOT,
full `main.py`, procedure-PDF, canonical-five runtime, or production validation.

---

## 14. Review bundle

Create repository-root temporary:

```text
kaonlt_review.diff
```

It must include:

- complete tracked diffs for every modified path;
- complete additions for every new file using `git diff --no-index /dev/null`;
- no unrelated paths.

Do not stage merely for review.

Stop before commit/push.

---

## 15. Acceptance criteria

Return:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

only if:

- exact start identity passes;
- unrelated local state is preserved;
- reviewed F.2 is staged/collected alongside F.3/F.4;
- runtime discovers exact reviewed F.2/F.3/F.4;
- exact-lineage behavior is preserved;
- new-lineage path rebuilds current F.2/F.3/F.4 through unchanged public builders;
- one source-owned explicit scientific projection is used;
- all three projections must match exactly;
- unknown scientific fields remain compared;
- candidate wrapper/hash/fingerprint identity is fail-closed;
- current factors come only from current-lineage F.4 review data;
- candidate and current provenance remain separate;
- default strict F.3/F.4 validators remain unchanged;
- no scientific algorithms/tolerances/cuts/normalizations/yields change;
- all focused negative tests and regressions pass;
- memory/manifest/bootstrap/diff checks pass;
- only allowlisted paths change;
- complete review bundle exists;
- no farm/stage/commit/push/ref update occurred.

---

## 16. Hard stops

Return `BLOCKED` without workaround if:

- branch/HEAD/origin differ;
- unrelated state cannot be preserved;
- candidate F.2 cannot be handed off through the tracked owner without changing
  another unallowlisted path;
- exact F.2/F.3/F.4 scientific equivalence cannot be demonstrated with the
  closed projection above;
- safe implementation requires changing F.1/F.2/F.3/F.4 scientific builders;
- safe implementation requires changing F.4 parent-preservation mathematics;
- safe implementation requires weakening default persisted-artifact validation;
- a new tolerance is required;
- a scientific field would need to be newly excluded from comparison;
- candidate provenance would need to be relabeled as current provenance;
- runtime source would import `testing/`;
- factors would come from stale candidate F.1 lineage after the current F.1
  identity changed;
- a path outside the allowlist is required;
- hard memory/manifest checks fail.

Do not run the farm from this task.
