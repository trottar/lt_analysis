# KaonLT — E.8.4 Fix.5.6 farm-owner readiness and failure provenance

## 1. Purpose

Repair one narrow operational weakness in the tracked `Q4p4W2p74 / Left / lowe`
E.8.4 Fix.5 owner before another expensive farm run.

The previous tracked owner completed the full analysis/render stage but returned
no final ZIP. The existing per-run `.log` is only the captured
`run_Prod_Analysis.sh` subprocess output. Failures occurring afterward in
artifact verification, page-manifest checks, detached collector/source checks,
ZIP collection, or ZIP verification are emitted by the owner to stderr only and
are not persisted in that analysis log.

Therefore the prior grep of the analysis log could not have shown a post-analysis
`E.8.4.Fix.5 gate failed: ...` message even if such a failure occurred.

Fix.5.6 must:

1. preserve a durable, per-attempt owner status/failure record with exact stage
   and reason;
2. run the collector's deterministic source checks in the exact detached pushed
   worktree **before** the expensive analysis;
3. preserve the existing post-analysis artifact verification and final detached
   collection/ZIP verification;
4. change no analysis physics, presentation science, profile scope, collector
   semantics, launcher behavior, or artifact definitions.

No farm run is authorized by this contract.

---

## 2. Exact starting state

Required branch:

```text
test
```

Required starting HEAD:

```text
df957a6414fc9c515d1f82228517cb801dc90350
```

Commit subject:

```text
E8.4 Fix.5.5: improve current-lineage visualization clarity
```

Required parent:

```text
761fbb6c03d2d7a10bb911cf84e9ba898496fab6
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
- preserve the worktree;
- `docs/memory/phases/workflow-continuity-hardening-task-contract.md` remains
  unrelated user-owned untracked work;
- `kaonlt_review.diff` remains temporary review material;
- known farm-generated model output dirt remains preserved/reported, never
  cleaned;
- root `AGENTS.md` and `.codex/` remain local-only/untracked;
- unrelated tracked changes are a blocker;
- do not reset, clean, stash, commit, push, or run the farm.

---

## 3. Required startup reading

Read in order:

1. root `AGENTS.md`;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read:

- `docs/memory/phases/e8-4-fix5-5-current-lineage-visualization-clarity.md`
- `docs/memory/phases/e8-4-fix5-5-current-lineage-visualization-clarity-task-contract.md`
- `docs/memory/phases/e8-4-fix5-shareable-method-a-impact-pages.md`
- `docs/memory/investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/TOOLS.md`
- `docs/memory/CODEX.md`

Relevant owner/collector source and tests:

- `testing/run_e8_4_fix5_left_lowe_plot_gate.py`
- `testing/test_run_e8_4_fix5_left_lowe_plot_gate.py`
- `testing/collect_pion_hgcer_validation_bundle.py`
- `testing/test_collect_pion_hgcer_validation_bundle.py`
- `testing/pion_hgcer_validation_bundle_profile_e8_4_fix5.json`
- `run_Prod_Analysis.sh`

Do not reopen analysis-source implementation.

---

## 4. Pushed-state facts to record

Independent pushed-state review established:

```text
remote test:
df957a6414fc9c515d1f82228517cb801dc90350
```

with exactly one commit over:

```text
761fbb6c03d2d7a10bb911cf84e9ba898496fab6
```

and exactly the 12 previously reviewed Fix.5.5 paths.

Fix.5.5 remains:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

Actual ROOT/PDF visual quality remains farm-only.

Fix.5.4 numerical/provenance source remains:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

with source review complete and no new farm closure.

Method A remains detached/non-production.
Method B remains diagnostic-only/numerically absent.
Absolute SIMC amplitude interpretation remains blocked by the existing
source-provenance limitation.

---

## 5. Existing owner path — preserve its scientific behavior

The tracked owner currently performs:

```text
exact source/branch/origin preflight
-> dynamic in-memory profile binding to --source-commit
-> frozen candidate-hash checks
-> ./run_Prod_Analysis.sh -d 4p4 2p74
-> completion-marker checks
-> fresh artifact checks
-> page-manifest checks
-> run-summary creation
-> exact detached-source collector
-> exact Left/lowe generic bundle collection
-> ZIP provenance/hash/inventory verification
-> exactly one ZIP path on stdout
```

Preserve this behavior.

The profile remains canonical-five declaratively while the owner requests only:

```text
Left / lowe
```

The collector must continue to run from a clean detached worktree at the exact
pushed source commit while reading artifacts from the ordinary canonical farm
output directory.

Do not change the frozen candidate hashes.

---

## 6. Confirmed diagnostic gap from the previous failed owner attempt

The existing:

```text
KaonLT_E8_4_Fix5_Left_lowe_Q4p4W2p74_<timestamp>.log
```

is created inside `run_analysis()` and contains only the child
`run_Prod_Analysis.sh` stdout/stderr stream.

The following failures occur outside that child log:

- completion marker failure;
- artifact missing/stale/not-refreshed failure;
- strict JSON failure;
- page-manifest contract failure;
- head change after run;
- run-summary write failure;
- detached worktree/source-check failure;
- collector exception/failure;
- ZIP integrity/provenance/inventory failure.

`main()` catches the owner exception and prints:

```text
E.8.4.Fix.5 gate failed: <reason>
```

to owner stderr only.

Therefore absence of that string from the analysis subprocess `.log` is not
evidence that no post-analysis owner failure occurred.

Record this distinction durably. Do not speculate which exact post-analysis
stage failed in the prior attempt without direct evidence.

---

## 7. Required repair A — persistent owner-status sidecar

Add one deterministic per-attempt owner-status JSON sidecar.

Recommended basename:

```text
KaonLT_E8_4_Fix5_Left_lowe_Q4p4W2p74_<same-output-timestamp>-gate-status.json
```

Use the same unique output stem as the intended ZIP.

Place it in the ordinary canonical analysis artifact directory:

```text
/group/c-kaonlt/USERS/trottar/lt_analysis/OUTPUT/Analysis/KaonLT
```

It is an owner diagnostic/provenance artifact, not a production-analysis
artifact and not a new scientific input.

### Required schema content

At minimum persist:

```text
schema_version
source_commit
kinematic
phi
epsilon
expected_zip_path
analysis_log_path
started_at_utc
updated_at_utc
status
stage
failure_reason
analysis_started
analysis_completed
artifact_verification_completed
collector_source_preflight_completed
collection_completed
zip_verification_completed
```

Allowed status values:

```text
running
failed
success
```

Use stable explicit stage names, for example:

```text
preflight
profile
candidate_identity
collector_source_preflight
analysis
completion_markers
verify_artifacts
write_run_summary
collection
verify_zip
complete
```

The exact names may differ if equally deterministic and tested.

### Persistence requirements

- Create/update the sidecar atomically.
- Write it before the expensive analysis begins.
- Update it at meaningful stage transitions.
- On any caught owner failure, set:
  - `status = "failed"`
  - exact `stage`
  - literal failure reason/string
  - final timestamp.
- On success, set:
  - `status = "success"`
  - `stage = "complete"`
  - `collection_completed = true`
  - `zip_verification_completed = true`.
- Do not put Python traceback noise into the structured reason unless the
  exception itself is unexpected; preserve deterministic known reasons.
- Do not overwrite an unrelated existing status sidecar.
- Do not change or append to the child analysis log after its hash is recorded
  in the run summary.

The owner should print the status-sidecar path to stderr on failure so the user
can inspect one durable file.

The success stdout contract remains exactly one POSIX ZIP path and nothing else.

### Bundle scope

Do **not** add the gate-status sidecar to the existing declarative validation
profile in this task.

It exists to diagnose owner orchestration. The accepted ZIP inventory remains
unchanged.

---

## 8. Required repair B — pre-run detached collector source checks

Before calling:

```text
./run_Prod_Analysis.sh -d 4p4 2p74
```

run the collector's existing deterministic source checks from a clean detached
worktree at the exact supplied/pushed source commit.

Reuse the existing:

```text
clean_collection_worktree(...)
collection_module(...)
```

and the collector's existing source-check implementation.

The effective profile must already be dynamically bound in memory to:

```text
required_analysis_commit = source_commit
```

with no tracked repin.

The pre-run source check must verify at minimum:

- detached worktree HEAD equals source commit;
- source commit is its own required analysis commit for this effective profile;
- source checks all return zero;
- committed source identity produces no unexpected paths under that effective
  profile;
- the collector/profile can be imported/read in the farm Python environment.

If any deterministic collector/source check fails:

```text
STOP BEFORE ANALYSIS
```

and persist the exact failure in the owner-status sidecar.

Do not create the final ZIP during this preflight.

The ordinary final collector invocation after fresh artifact verification must
remain in place and must still re-check the source and build/verify the actual
bundle.

---

## 9. Artifact and page verification — preserve current contract

Do not weaken:

- source HEAD/origin identity;
- branch identity;
- known-dirt allowlist;
- candidate hash identity;
- artifact freshness;
- strict JSON parsing;
- PDF magic check;
- `renderer_failures == []`;
- critical E.8.2/E.8.3/E.8.4 page presence;
- all Fix.5 new page IDs;
- canonical nine-phi represented inventory;
- final ZIP completeness;
- exact source identity;
- exact Left/lowe requested setting;
- profile identity;
- ZIP hashes/byte sizes/inventory.

The Fix.5.4 absolute-SIMC page families may remain page-manifest records with
their existing `available=false` provenance state. Do not reinterpret that as a
renderer failure.

No page ID or expected bundle artifact is changed here.

---

## 10. Required failure-path tests

Extend `testing/test_run_e8_4_fix5_left_lowe_plot_gate.py`.

At minimum test:

### Pre-run source-check failure

Simulate a detached collector source check failure.

Require:

- analysis is not called;
- output ZIP does not exist;
- owner-status sidecar exists;
- `status == "failed"`;
- stage identifies collector/source preflight;
- literal reason is present.

### Analysis failure

Require:

- owner-status sidecar identifies analysis stage;
- analysis log path is recorded;
- no ZIP success stdout semantics.

### Completion-marker failure

Require exact stage/reason.

### Artifact freshness/verification failure after a successful analysis

Require:

- analysis completed flag is true;
- artifact verification flag remains false;
- exact failing stage/reason is persisted;
- no ZIP is reported.

### Collector failure after successful analysis/artifact checks

Require:

- analysis completed true;
- artifact verification completed true;
- collection completed false;
- stage/reason persisted.

### ZIP verification failure

Require corresponding stage and flag behavior.

### Full success

Require:

- existing ZIP contents/provenance tests still pass;
- status sidecar ends `success / complete`;
- all completion booleans are true;
- stdout remains exactly one ZIP path;
- status sidecar is not added to the bundle profile/inventory.

### Existing safety regressions

Preserve tests proving:

- ordinary checkout is never reset/cleaned/stashed;
- temporary detached worktree cleanup is path-bounded;
- known farm-generated model-output dirt is preserved;
- wrong branch/head/origin/unexpected source dirt fail before analysis;
- existing ZIP is refused;
- exact `-d 4p4 2p74` launcher invocation remains.

---

## 11. Allowed substantive files

Expected:

```text
testing/run_e8_4_fix5_left_lowe_plot_gate.py
testing/test_run_e8_4_fix5_left_lowe_plot_gate.py
```

Read-only/frozen unless a concrete blocker requires a new contract:

```text
testing/collect_pion_hgcer_validation_bundle.py
testing/test_collect_pion_hgcer_validation_bundle.py
testing/pion_hgcer_validation_bundle_profile_e8_4_fix5.json
run_Prod_Analysis.sh
src/
```

Do not modify scientific source.

---

## 12. Memory updates

Allowed:

```text
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/LEARNINGS.md
docs/memory/TOOLS.md
docs/memory/roadmap/STATUS.md
docs/memory/investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md
docs/memory/phases/e8-4-fix5-5-current-lineage-visualization-clarity.md
docs/memory/phases/e8-4-fix5-6-owner-farm-readiness-and-failure-provenance-task-contract.md
docs/memory/phases/e8-4-fix5-6-owner-farm-readiness-and-failure-provenance.md
docs/memory/manifest.json
```

Update only warranted durable knowledge.

### Required reconciliation

Record pushed Fix.5.5 source:

```text
df957a6414fc9c515d1f82228517cb801dc90350
```

and that independent actual-diff and pushed-state review passed.

Remove stale wording that calls Fix.5.5 local/unreviewed.

### Fix.5.6 phase

Create:

```text
docs/memory/phases/e8-4-fix5-6-owner-farm-readiness-and-failure-provenance.md
```

During implementation:

```text
ACTIVE
```

If local deterministic checks pass:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

Do not call the farm owner runtime-validated before a real farm attempt.

### CURRENT / NEXT after local completion

Use:

```text
Fix.5.6 actual-diff review
-> user commit/push
-> pushed-state review
-> one narrow Q4p4W2p74 / Left / lowe tracked-owner farm gate
-> fresh ZIP/status/log/PDF/manifest evidence review
```

No farm command is authorized inside this task.

Regenerate `docs/memory/manifest.json` for the intended versioned candidate,
excluding unrelated user-owned untracked files.

---

## 13. Local validation

Use the repository-selected Python interpreter.

Run syntax checks for changed/new Python files.

Run at minimum:

```text
testing.test_run_e8_4_fix5_left_lowe_plot_gate
testing.test_collect_pion_hgcer_validation_bundle
testing.test_e8_4_fix5_5_visualization
testing.test_e8_4_fix5_4_identity_audit
```

plus any directly affected owner/profile regression suites.

Run:

```bash
git diff --check
```

Run manifest regeneration/check and:

```bash
<PYTHON> -B tools/check_memory_health.py --root .
```

Report exact pass/skip counts and warnings.

Local tests do not establish farm filesystem permissions, ROOT/PyROOT, full
analysis execution, final ZIP delivery, or PDF rendering.

---

## 14. Diff audit

Before stopping:

```bash
git status --short --untracked-files=all
git diff --stat
git diff --check
```

Preserve and exclude from the intended candidate:

```text
docs/memory/phases/workflow-continuity-hardening-task-contract.md
kaonlt_review.diff
```

Produce one complete replacement:

```text
kaonlt_review.diff
```

relative to:

```text
df957a6414fc9c515d1f82228517cb801dc90350
```

Include:

- complete tracked diff;
- complete additions for every intended new file;
- no unrelated user-owned untracked file.

Do not stage merely to construct the review bundle.

---

## 15. Acceptance criteria

PASS only if:

1. exact starting HEAD is respected;
2. no scientific or presentation source changes;
3. the owner-status sidecar persists exact stage/failure reason;
4. a collector/source deterministic failure blocks before expensive analysis;
5. existing post-analysis artifact/page/ZIP validation is not weakened;
6. final collector still runs from exact detached pushed source;
7. successful stdout remains exactly one ZIP path;
8. failure stdout remains empty;
9. no destructive ordinary-checkout operation is introduced;
10. profile/bundle inventory remains unchanged;
11. deterministic owner/collector tests pass;
12. warranted memory is accurate and manifest regenerated;
13. no commit, push, or farm run occurs.

---

## 16. Farm boundary

This task does not establish:

- JLab farm execution;
- ROOT/PyROOT behavior;
- procedure-PDF visual quality;
- actual Fix.5.4 numerical closure on farm objects;
- actual signed-cancellation explanation;
- final ZIP delivery;
- filesystem/runtime behavior of the new status sidecar;
- Method-A promotion.

Those require the next narrow farm gate.

---

## 17. Hard stop

After implementation, deterministic local checks, warranted memory updates,
manifest regeneration, and complete review-diff preparation:

**STOP.**

Do not commit.

Do not push.

Do not run the farm.

Return the complete `kaonlt_review.diff` to ChatGPT for independent actual-diff
review.
