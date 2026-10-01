# F.4.Refresh.2.Validation.2 — tracked materialize → verify → package execution owner Codex contract

## Objective and exact starting HEAD

Implement the missing tracked execution owner for the **F.4.Refresh.2 current-baseline candidate farm review**.

The exact starting pushed `test` HEAD is:

```text
c406c138285727503b12115177dfb8bc7efcb7fe
Harden post-push CURRENT continuity
```

Before making any change, Codex must independently verify:

```bash
git branch --show-current
git rev-parse HEAD
git status --short --untracked-files=all
```

Required start state:

- branch: `test`;
- HEAD: exactly `c406c138285727503b12115177dfb8bc7efcb7fe`;
- the only permitted pre-existing worktree change is this newly placed task contract:
  `?? docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner-task-contract.md`.

If the branch, HEAD, or worktree differs from that exact state, **STOP** and report the exact state. Do not adapt this contract to a different base without a new ChatGPT audit/contract.

The task is orchestration only. It closes the currently documented ownership gap:

```text
validated current F.1 canonical-five inputs
+
accepted F.2/F.3/F.4 authorities
+
validated F.4.Refresh.1 comparison
        ↓
reviewed F.4.Refresh.2 materializer
        ↓
explicit post-materialization verification
        ↓
existing F.4.Refresh.2 generic collector/profile via existing bundle-only wrapper
        ↓
exactly one fresh review ZIP
```

No physics, production, accepted authority, Method-A promotion, or Method-B behavior is in scope.

---

## Repository-memory startup

Before implementation, read in full and in the repository-required order:

1. `AGENTS.md` if present in the local checkout;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read only the CURRENT-direct/task-relevant records needed for this task, including at minimum:

- `docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization.md`
- `docs/memory/phases/f4-refresh2-validation1-farm-materialization-bundle-profile.md`
- `docs/memory/investigations/f4-refresh2-validation1-operational-readiness-failure.md`
- `docs/memory/phases/memory-health-operational-completeness-hardening.md`
- `docs/memory/phases/post-push-current-continuity.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/templates/CODEX_CONTRACT.md`

Preserve CURRENT as the active-state/NEXT authority.

---

## Established source facts that this task must preserve

The reviewed/pushed F.4.Refresh.2 materializer source is:

```text
141a3d04f9e5d07be21dba14e0e63212c3990bf1
```

The validated F.4.Refresh.1 comparison artifact SHA-256 is:

```text
c16287e9192288ed5f116b5eebf95f49261d7fd09755ee923070154e047cf7b5
```

The materializer is:

```text
testing/materialize_method_a_current_baseline_authority.py
```

Its reviewed CLI owns scientific candidate construction and already accepts:

- five repeated `--f1 ALIAS=PATH` inputs;
- `--accepted-f2 PATH`;
- `--accepted-f3 PATH`;
- `--accepted-f4 PATH`;
- `--comparison PATH`;
- `--expected-comparison-sha256 SHA256`;
- `--source-head SHA40`;
- `--output-dir PATH`;
- optional `--overwrite`.

This new execution owner must invoke that CLI. It must not import or duplicate the scientific builders to reconstruct F.2/F.3/F.4 itself.

The F.4.Refresh.2 bundle profile is:

```text
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
```

It declares exactly five required global JSON artifacts:

1. reviewed F.4.Refresh.1 comparison-input copy;
2. current-baseline candidate F.2;
3. current-baseline candidate F.3;
4. current-baseline candidate F.4;
5. F.4.Refresh.2 materialization manifest.

Its `required_analysis_commit` remains exactly `141a3d04f9e5d07be21dba14e0e63212c3990bf1`.

The generic collector is:

```text
testing/collect_pion_hgcer_validation_bundle.py
```

The existing package-only farm wrapper is:

```text
testing/package_pion_hgcer_validation_bundle.tcsh
```

That wrapper deliberately packages **existing** artifacts only. It owns the canonical artifact-root restriction, detached bundle worktree, generic collector invocation, immutable-file SHA checks, canonical Globus destination, ZIP integrity check, and returned ZIP printout. It must remain bundle-only.

The existing canonical artifact root used by that reviewed wrapper is:

```text
/group/c-kaonlt/USERS/trottar/lt_analysis/OUTPUT/Analysis/KaonLT
```

The existing canonical Globus transfer directory is:

```text
/volatile/hallc/c-kaonlt/trottar/globus
```

Do not invent alternate farm paths.

---

## Allowed files

Substantive changes are limited to:

```text
NEW  docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner-task-contract.md

NEW  testing/run_f4_refresh2_materialize_verify_package.py
NEW  testing/test_run_f4_refresh2_materialize_verify_package.py

MOD  testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
MOD  testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py

NEW  docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner.md
MOD  docs/memory/CURRENT.md
MOD  docs/memory/roadmap/STATUS.md
MOD  docs/memory/manifest.json
```

The profile change is limited to extending `source_identity.allowed_committed_files` so the future pushed execution-owner source and its focused test are permitted after the frozen materializer source. Do not change:

- profile schema;
- validation-profile identity;
- canonical-five setting inventory;
- required artifact declarations;
- `required_analysis_commit`;
- `allowed_non_analysis_path_prefixes`.

No other file may change substantively.

The task contract itself belongs directly in tracked `docs/memory/phases/`; do not stage it through the repository root.

The temporary repository-root review bundle required by repository workflow is the root-level untracked artifact:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

Do not stage or commit that review bundle.

---

## Frozen files and interfaces

The following are frozen for this task:

```text
testing/materialize_method_a_current_baseline_authority.py
testing/test_materialize_method_a_current_baseline_authority.py

testing/compare_method_a_current_baseline_authority.py
testing/test_compare_method_a_current_baseline_authority.py

testing/collect_pion_hgcer_validation_bundle.py
testing/test_collect_pion_hgcer_validation_bundle.py

testing/package_pion_hgcer_validation_bundle.tcsh
testing/test_package_pion_hgcer_validation_bundle_tcsh.py
```

Also frozen:

- all `src/` scientific/runtime files;
- `main.py` and analysis entry points;
- `run_Prod_Analysis.sh`;
- random subtraction;
- slow-proton subtraction;
- pion-component subtraction;
- all F.2/F.3/F.4 scientific builders;
- accepted F.2/F.3/F.4 authority artifacts and their identities;
- frozen F.6.2 science/fingerprints;
- F.5/F.6.3/E.8.4 scientific/runtime interfaces;
- Method-B calculations and diagnostics;
- production pion weights and yields;
- cuts, normalizations, priors, templates, binning, fit windows, subtraction formulas, and component definitions.

If implementation appears to require any frozen file to change, **STOP** and report the concrete blocker.

---

## Scientific ownership

This phase owns **orchestration and deterministic verification only**.

It must not:

- calculate a new scientific quantity;
- change any F.2/F.3/F.4 formula;
- change accepted authority;
- reinterpret F.4.Refresh.1 results;
- alter Method-A correction values;
- apply Method A to production;
- involve Method B numerically;
- alter baseline production;
- alter F.6.3 or E.8.4 runtime behavior.

The scientific result expected from the existing materializer remains exactly:

```text
f2_scientific_payload_match = true
f3_scientific_payload_match = true
f4_scientific_payload_match = false
first_changed_stage = F4
```

The new verifier may confirm that the materializer manifest records this reviewed result. It must not recompute the scientific comparison independently.

---

## Required implementation architecture

Implement a narrow **pure-Python CLI orchestration driver**:

```text
testing/run_f4_refresh2_materialize_verify_package.py
```

The Python driver is preferred because the missing middle step is structured JSON/hash verification and deterministic process orchestration. Do not create a second scientific implementation and do not replace the existing `tcsh` package wrapper.

The driver must make the complete operation explicit and fail closed.

### 1. Preflight before any candidate output is written

Before invoking the materializer, verify at minimum:

- the driver is running from the KaonLT repository;
- current Git HEAD is exactly the user-supplied full 40-character `--bundle-commit`;
- the worktree is clean;
- the frozen materializer source commit `141a3d04f9e5d07be21dba14e0e63212c3990bf1` is an ancestor of HEAD;
- changes after that source commit are allowed by the F.4.Refresh.2 profile source-identity rule;
- the exact F.4.Refresh.2 profile loads successfully;
- its `required_analysis_commit` remains the frozen materializer source;
- the requested kinematic is exactly `Q4p4W2p74` for this narrow gate;
- the artifact directory is the canonical KaonLT artifact root already enforced by the reviewed package wrapper;
- all explicitly supplied five F.1 paths, accepted F.2/F.3/F.4 paths, and comparison path exist as regular files;
- the five F.1 aliases are exactly the canonical-five aliases expected by the existing materializer, with no duplicate/missing alias;
- the supplied comparison file SHA-256 equals the frozen reviewed comparison SHA above;
- the intended ZIP does not already exist;
- none of the five expected F.4.Refresh.2 materialization outputs already exists.

The last requirement is deliberate: this execution owner is for a **fresh** materialization gate. Do **not** use `--overwrite`, delete old outputs automatically, rename stale outputs automatically, or silently reuse stale/partial outputs.

If preflight fails, return nonzero before materialization and print a specific STOP/failure reason.

### 2. Invoke the existing reviewed materializer

Invoke:

```text
testing/materialize_method_a_current_baseline_authority.py
```

through its real CLI as a subprocess.

Pass:

- the five explicit F.1 inputs;
- the explicit accepted F.2/F.3/F.4 inputs;
- the explicit reviewed comparison file;
- `--expected-comparison-sha256 c16287e9192288ed5f116b5eebf95f49261d7fd09755ee923070154e047cf7b5`;
- `--source-head 141a3d04f9e5d07be21dba14e0e63212c3990bf1`;
- the canonical artifact directory as `--output-dir`.

Do not pass `--overwrite`.

A nonzero materializer return code is fatal. Do not package anything after it.

### 3. Explicitly verify the completed materialization

After the materializer exits successfully and before the package wrapper is invoked, independently verify the completed output set.

At minimum verify:

- exactly the five expected task output files now exist as regular files;
- the copied comparison-output SHA-256 equals the frozen reviewed comparison SHA;
- the materialization manifest parses as strict JSON;
- manifest schema is exactly:

```text
method_a_current_baseline_authority_materialization/v1
```

- manifest flags are exactly consistent with detached/non-production ownership:
  - `non_authoritative == true`;
  - `accepted_authority_mutated == false`;
  - `production_objects_mutated == false`;
  - `production_application_performed == false`;
  - `method_a_promoted == false`;
- manifest `source_head` equals the frozen reviewed materializer source `141a3d04f9e5d07be21dba14e0e63212c3990bf1`;
- manifest `kinematic_token == "Q4p4W2p74"`;
- manifest comparison-input source/copy hashes equal the frozen reviewed comparison SHA;
- the materializer scientific-gate summary is exactly F.2 match / F.3 match / F.4 mismatch / first changed stage F4;
- candidate F.2/F.3/F.4 basenames in the manifest are exactly the expected current-baseline candidate basenames;
- each candidate F.2/F.3/F.4 raw SHA-256 recorded by the manifest equals the actual published file SHA-256;
- all five output files have their SHA-256 calculated for immutable packaging.

This verifier is a **completion/provenance check**, not a second scientific calculator.

Any mismatch is fatal. Do not invoke packaging after a verification failure.

### 4. Delegate packaging to the existing reviewed wrapper

Invoke the existing:

```text
testing/package_pion_hgcer_validation_bundle.tcsh
```

Do not duplicate its detached-worktree or collector logic.

Pass the exact current `--bundle-commit`, the existing F.4.Refresh.2 profile, canonical artifact root, `Q4p4W2p74`, the explicitly requested fresh ZIP destination/basename, and all five newly materialized files as repeated:

```text
--immutable PATH SHA256
```

using the hashes captured by the post-materialization verifier.

The wrapper must remain the owner of:

- detached bundle worktree creation/removal;
- generic collector invocation;
- canonical Globus path enforcement;
- before/after immutable SHA checks;
- ZIP integrity testing.

A nonzero wrapper return code is fatal.

### 5. Verify the returned review ZIP

After wrapper success, the driver must verify the exact requested ZIP, without scanning for or accepting some other stale ZIP.

At minimum:

- exact requested ZIP exists as a regular file;
- ZIP opens and is structurally readable;
- `manifest.json` exists and parses;
- collector manifest has `complete == true`;
- collector manifest validation-profile identity is exactly the F.4.Refresh.2 profile identity;
- requested kinematic is `Q4p4W2p74`;
- all five required global artifact records are present;
- each packaged global artifact SHA-256 equals the immutable SHA captured before packaging;
- there is no silent substitution of an alternate artifact;
- print the exact returned ZIP path and its SHA-256 on success.

Do not interpret `complete=true` as scientific acceptance. It establishes successful packaging/source checks only; ChatGPT will inspect the returned materialization/candidate evidence in the later farm-evidence review.

---

## CLI requirements

Keep the driver input-driven for artifact inputs. Do not guess accepted-authority or F.1 farm paths.

The CLI must take enough explicit information to run the operation reproducibly, including at least:

```text
--python PATH
--bundle-commit SHA40
--f1 ALIAS=PATH          # exactly five, repeat
--accepted-f2 PATH
--accepted-f3 PATH
--accepted-f4 PATH
--comparison PATH
--artifact-dir PATH
--kinematic Q4p4W2p74
--output ZIP_OR_CANONICAL_GLOBUS_PATH
```

The F.4.Refresh.2 profile path and frozen comparison/materializer identities may be fixed task constants rather than user-overridable values.

Reject unknown/duplicate aliases, unsafe/malformed commit identities, and unsupported kinematics.

Do not add arbitrary command execution, `eval`, shell-string interpolation, scheduler submission, analysis execution, or a general-purpose hook.

Use argument-vector subprocess execution, never `shell=True`.

---

## Failure and stale-state behavior

The driver must fail closed.

Required behavior:

- preflight failure: no candidate outputs created;
- materializer failure: no package invocation;
- verification failure: no package invocation;
- package failure: no success claim;
- returned-ZIP verification failure: no success claim;
- stale/partial candidate outputs at next invocation: preflight fails rather than reusing or overwriting them;
- existing ZIP: fail rather than overwrite;
- no automatic deletion of accepted or candidate artifacts;
- no automatic cleanup that could conceal failed candidate evidence;
- no fallback to an older comparison, older authority, alternate kinematic, alternate profile, or alternate ZIP.

Do not use `git reset`, `git clean`, `git stash`, commits, pushes, scheduler commands, or farm analysis commands.

---

## Positive tests

Create focused deterministic tests in:

```text
testing/test_run_f4_refresh2_materialize_verify_package.py
```

Use temporary directories, synthetic JSON/ZIP fixtures, and an injectable/mock command runner so the local tests require neither ROOT nor the JLab farm.

Cover at minimum:

1. exact CLI/preflight accepts a clean expected source state and canonical-five inputs;
2. subprocess order is materializer → verification → package wrapper → ZIP verification;
3. successful synthetic materializer outputs produce five immutable SHA arguments for the wrapper;
4. valid materialization manifest fields/hashes pass;
5. valid synthetic collector ZIP with exact five packaged hashes passes;
6. final success returns/prints the exact requested ZIP and ZIP SHA-256;
7. F.4.Refresh.2 profile accepts the new execution-owner source/test as reviewed post-materializer paths.

---

## Negative tests

Cover at minimum:

- wrong/short/malformed bundle commit;
- current HEAD differs from supplied bundle commit;
- dirty worktree;
- required materializer commit not an ancestor;
- unexpected changed scientific/materializer path after the frozen source;
- wrong profile identity or changed required-analysis commit;
- unsupported kinematic;
- missing F.1/accepted/comparison input;
- duplicate/missing/unknown F.1 alias;
- comparison SHA mismatch;
- stale existing comparison-copy/candidate/manifest output;
- existing ZIP;
- materializer nonzero return;
- missing materialization manifest;
- malformed manifest JSON;
- wrong materialization schema;
- any detached/non-production flag mismatch;
- wrong materializer source head;
- wrong kinematic token;
- wrong scientific-gate summary;
- candidate basename mismatch;
- candidate raw-hash mismatch;
- package wrapper nonzero return;
- package wrapper not called after any earlier failure;
- missing returned ZIP;
- malformed ZIP;
- missing/malformed ZIP manifest;
- ZIP manifest `complete != true`;
- wrong validation-profile identity;
- missing global artifact record;
- packaged artifact SHA mismatch.

---

## Regression tests

The scoped change must preserve existing reviewed behavior.

At minimum run:

```bash
python -m py_compile \
  testing/run_f4_refresh2_materialize_verify_package.py \
  testing/test_run_f4_refresh2_materialize_verify_package.py

python -m unittest testing.test_run_f4_refresh2_materialize_verify_package
python -m unittest testing.test_pion_hgcer_validation_bundle_profile_f4_refresh2
python -m unittest testing.test_materialize_method_a_current_baseline_authority
python -m unittest testing.test_collect_pion_hgcer_validation_bundle
python -m unittest testing.test_package_pion_hgcer_validation_bundle_tcsh
```

If local `tcsh` is unavailable, the existing wrapper test may report its existing syntax-test skip; record that exact skip. Do not treat an unavailable local `tcsh` as farm validation.

Also run the repository memory-health test suite warranted by the memory changes.

---

## Forbidden shortcuts

Do not:

- change the materializer to make orchestration easier;
- change the generic collector;
- change the existing package wrapper;
- duplicate F.2/F.3/F.4 builders or comparison logic;
- weaken source-identity checks;
- move `required_analysis_commit` forward from `141a3d04...` merely to permit the new owner;
- broaden `allowed_committed_files` beyond the new owner, its focused test, and the already-reviewed profile/test;
- add a blanket `testing/` allowlist;
- treat JSON parseability as sufficient materialization verification;
- treat collector `complete=true` as scientific acceptance;
- use `--overwrite`;
- delete stale candidate outputs automatically;
- silently reuse stale outputs;
- mutate accepted authorities;
- alter production physics;
- apply Method A to production;
- make Method B numerical;
- run ROOT, `main.py`, `run_Prod_Analysis.sh`, or a farm analysis;
- commit, push, or update remote refs.

---

## Local validation and diff audit

After implementation and warranted memory updates:

1. Run the focused/regression tests above.
2. Regenerate the memory manifest:

```bash
python -B tools/update_memory_manifest.py --root . --write
```

3. Verify the memory manifest:

```bash
python -B tools/update_memory_manifest.py --root . --check
```

4. Run deterministic memory/bootstrap checks:

```bash
python -B tools/memory_bootstrap.py --root .
python -m unittest testing.test_memory_health
```

5. Run strict task-final health:

```bash
python -B tools/check_memory_health.py --root . --fail-on-warning
```

6. Run:

```bash
git diff --check
git status --short --untracked-files=all
```

7. Audit the actual changed-path inventory. No substantive path outside the allowlist in this contract is permitted.

8. Create one fresh temporary repository-root cumulative review bundle:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

It must contain:

- branch and committed starting HEAD;
- `git status --short --untracked-files=all`;
- cumulative tracked diff/stat;
- no-index diffs for intended untracked source/memory files;
- exact test/check outputs;
- memory-health report;
- candidate path inventory.

The review bundle must remain untracked.

---

## Memory updates

Create:

```text
docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner.md
```

Record:

- exact base HEAD;
- execution-owner scope;
- delegated existing components;
- scientific/farm boundary;
- exact deterministic checks run and results;
- no farm execution claim;
- no accepted-authority mutation;
- no production or Method-A-promotion claim.

Update `CURRENT.md` and `roadmap/STATUS.md` only as warranted by the completed local implementation.

If all implementation checks pass, the new local implementation may be recorded as:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

but the **farm execution gate remains BLOCKED** until:

1. ChatGPT independently reviews the actual cumulative diff/runtime path;
2. the separate required final-pre-push memory/status reconciliation is contracted, implemented, and reviewed;
3. the user commits/pushes;
4. ChatGPT reviews the pushed state.

Do not mark this phase `SOURCE REVIEWED`; Codex cannot confer independent source review on its own work.

CURRENT's immediate substantive NEXT after this implementation must be the independent ChatGPT actual-diff/source-runtime-path review, followed—only if that passes—by the separately contracted final-pre-push reconciliation. Do not make commit/push the sole NEXT.

Regenerate `docs/memory/manifest.json` after all versioned memory changes.

Do not modify `CURRENT_HANDOFF.md` unless an actual exceptional transfer state arises.

---

## Required memory-health report

At task completion report exactly:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
manifest check: PASS | FAIL
```

Task-final health must use:

```bash
python -B tools/check_memory_health.py --root . --fail-on-warning
```

Any warning blocks handoff unless the user explicitly approves a recorded maintenance exception.

---

## Farm-validation boundary

This task does **not** run the farm and does not validate:

- ROOT/PyROOT;
- `main.py`;
- production analysis;
- full kinematics;
- procedure PDFs;
- farm filesystem integration;
- actual candidate materialization from farm artifacts;
- actual packaging on farm;
- candidate authority acceptance;
- F.6.3/E.8.4 runtime closure;
- Method-A production promotion.

The farm gate remains `BLOCKED` throughout Codex implementation and until the full source-review/push/pushed-state-review chain above is complete.

Do not provide or execute a farm command in this task.

---

## Operational readiness chain after successful future push review

The intended later farm operation, if and only if this source task passes all later review/push gates, is:

```text
explicit validated F.1/F.2/F.3/F.4/comparison inputs
        ↓
NEW tracked execution owner
        ↓
existing reviewed materializer
        ↓
NEW explicit manifest/hash completion verifier
        ↓
existing reviewed package wrapper
        ↓
existing generic collector + F.4.Refresh.2 profile
        ↓
exactly one fresh Globus review ZIP
        ↓
ChatGPT evidence review
```

Every scientific calculation remains delegated to already-reviewed components. The new owner supplies missing invocation ownership and completion verification only.

---

## Acceptance criteria

This task passes implementation review only if:

- starting HEAD is exact;
- only allowed files changed substantively;
- new driver is orchestration-only;
- materializer, collector, and package wrapper remain unchanged;
- profile source identity remains anchored to `141a3d04...` and is extended only for the new owner/test;
- fresh-output preflight prevents stale/partial reuse;
- materializer is invoked through its real CLI without `--overwrite`;
- post-materialization schema/flags/provenance/scientific-summary/file-hash verification is explicit;
- package wrapper receives immutable hashes for all five task outputs;
- returned ZIP is explicitly verified against those hashes;
- all required positive/negative/regression tests pass locally, except only an explicitly recorded environment-dependent existing `tcsh` syntax-test skip if `tcsh` is unavailable;
- strict memory health has zero warnings;
- manifest verification passes;
- `git diff --check` passes;
- cumulative review bundle is complete and untracked;
- memory status does not overclaim source review, farm validation, accepted authority, or production promotion.

---

## Hard stop

Stop after:

1. implementation;
2. deterministic local checks;
3. warranted memory/history updates;
4. manifest regeneration/check;
5. strict memory-health report;
6. cumulative review-bundle creation.

Do **not** commit.
Do **not** push.
Do **not** run the farm.
Do **not** mutate accepted authority.
Do **not** begin F.6.3/E.8.4.
Do **not** perform the later final-pre-push reconciliation without a separate ChatGPT contract.

Return the exact changed-path inventory, exact test/check results, memory-health report, and review-bundle path to ChatGPT for independent actual-diff review.
