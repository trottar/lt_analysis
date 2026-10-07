# KaonLT E.8.2 Left/lowe runtime-owner OUTPUT-isolation repair — task contract

## 1. Task class and authority

This is a **narrow source-changing validation-infrastructure repair** for the
active E.8.2 Left/lowe runtime scientific-audit owner.

It does **not** change KaonLT production physics, the E.8.2 scientific chain,
the frozen collector, the validation profile, the production launcher, the
Diamond-cut implementation, Method A, Method B, or canonical-five provenance
work.

Follow repository-root `AGENTS.md` and the exact five-file startup core before
editing:

1. `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Authority remains:

```text
current source/diff
-> newest applicable farm/runtime evidence
-> tracked repository memory
-> older chat/history
```

The failed farm attempt is diagnostic evidence only. Do not treat it as E.8.2
runtime acceptance.

---

## 2. Exact starting identity and start gate

Required repository:

```text
https://github.com/trottar/lt_analysis/tree/test
```

Required branch:

```text
test
```

Required starting HEAD and local `origin/test`:

```text
f860f832f3e4b05f468cb6bb636f82f4a26d726b
```

Commit:

```text
Add isolated E8.2 Left lowe scientific audit owner
```

Before editing, establish:

```bash
git branch --show-current
git rev-parse HEAD
git rev-parse refs/remotes/origin/test
git status --short --untracked-files=all
git diff --check
```

The only expected new task file before implementation is:

```text
docs/memory/phases/e8-2-left-lowe-runtime-owner-output-isolation-repair-task-contract.md
```

Temporary `kaonlt_review*.diff` files are permitted only as review artifacts.

If branch, HEAD, local `origin/test`, or unrelated worktree state differs,
**STOP as `BLOCKED`**. Do not reset, stash, clean, overwrite, or repair unrelated
user state.

The Jefferson Lab farm attempt already failed and must **not** be rerun during
this implementation task.

---

## 3. Failed runtime gate and established root cause

The supplied farm attempt reached the ordinary Center Diamond-cut stage and
failed with:

```text
Kinematics:  Q4p4W2p74
Phi Setting:  Center
!!!!! ERROR !!!!!
 No valid file found!
!!!!! ERROR !!!!!

E.8.2 gate failed: analysis_failed
```

`RUNTIME VERIFIED` from the supplied terminal output:

- the E.8.2 owner started its child analysis;
- the child reached `DiamondPlot(...)` for the Center setting;
- no matching input ROOT file was visible through the child `OUTPATH`;
- the child exited nonzero and the owner correctly stopped as
  `analysis_failed`;
- no E.8.2 runtime acceptance follows.

`SOURCE VERIFIED` at starting HEAD:

1. `testing/run_e8_2_left_lowe_scientific_audit_gate.py` creates an owner-owned
   detached analysis worktree.
2. The reused copied-ltsep overlay changes only copied `LTANAPATH` to that
   detached worktree.
3. The resulting child `OUTPATH` is therefore the worktree-local lexical path:

   ```text
   <owned-worktree>/OUTPUT/Analysis/KaonLT
   ```

4. `run_Prod_Analysis.sh` executes:

   ```bash
   mkdir -p "${OUTPATH}"
   ```

   before its cleanup/symlink setup.
5. In `-d` mode the launcher intentionally skips `git clean -fdx`.
6. Therefore a fresh detached worktree acquires a real local `OUTPUT/`
   directory before `set_SymLinks.sh`.
7. `set_SymLinks.sh` expects `${LTANAPATH}/OUTPUT` to be absent or an existing
   symlink. It does not replace a pre-existing real directory with the intended
   external-output symlink.
8. `src/cuts/diamond.py` then searches the isolated worktree-local `OUTPATH`
   for the normal Center kaon ROOT inputs and finds none.

The ordinary farm checkout itself was not the cause. The current repair is the
new owner's filesystem topology, not the production analysis.

---

## 4. Scientific and production boundaries

The following remain frozen and **must not change**:

```text
run_Prod_Analysis.sh
set_SymLinks.sh
src/
tools/
farm_env/
background_samples/
testing/collect_pion_hgcer_validation_bundle.py
testing/pion_hgcer_validation_bundle_profile_e8_2_left_lowe.json
testing/run_e8_4_fix5_canonical_five_plot_gate.py
```

Do not change:

- random subtraction;
- dummy subtraction;
- slow-proton treatment;
- `prune_hist`;
- baseline pion subtraction;
- final kaon missing-mass/yield extraction;
- cuts;
- Diamond physics;
- SIMC;
- normalization;
- weights;
- uncertainties;
- binning;
- efficiency;
- acceptance;
- L/T separation;
- cross sections;
- active `no_empirical_residual` profile;
- Method-A/Method-B scientific ownership.

Method A remains detached/non-production and is not required by this gate.
Method B remains diagnostic/cross-check only and numerically excluded.

Canonical-five provenance/identity repair remains `DEFERRED`.

---

## 5. Exact allowed versioned paths

Only these paths may change:

```text
testing/run_e8_2_left_lowe_scientific_audit_gate.py
testing/test_run_e8_2_left_lowe_scientific_audit_gate.py
docs/memory/phases/e8-2-left-lowe-runtime-owner-output-isolation-repair-task-contract.md
docs/memory/phases/e8-2-left-lowe-runtime-owner-output-isolation-repair.md
docs/memory/CURRENT.md
docs/memory/manifest.json
```

No other versioned path may change.

The final diff must prove that all frozen production/scientific files remain
byte-identical.

---

## 6. Repair objective

Repair only the detached E.8.2 debug owner's repository-local OUTPUT topology.

Before the owner launches:

```text
./run_Prod_Analysis.sh -d 4p4 2p74
```

the owner must establish, inside the **owner-created disposable analysis
worktree only**, the same required OUTPUT indirection that the ordinary farm
checkout already has:

```text
<owned-worktree>/OUTPUT
    -> <validated VOLATILEPATH>/OUTPUT
```

The helper must use only paths already returned by the reviewed ltsep path
probe. Do not hard-code or invent a Jefferson Lab filesystem path.

After this link exists, the child lexical `OUTPATH` remains:

```text
<owned-worktree>/OUTPUT/Analysis/KaonLT
```

while resolving to the configured external artifact root already checked by
the owner.

The ordinary checkout and installed ltsep package must remain unchanged.

---

## 7. Required implementation

In:

```text
testing/run_e8_2_left_lowe_scientific_audit_gate.py
```

add one small owner-local helper with an explicit purpose such as:

```python
prepare_debug_output_link(worktree, paths, outdir)
```

Naming may differ if clearer, but behavior must satisfy all requirements below.

### 7.1 Inputs and authority

The helper must consume:

- the owner-created analysis worktree path;
- the path dictionary returned by the existing `probe_paths(...)`;
- the already resolved/validated `outdir`.

It must not call ltsep again to discover alternative paths and must not guess a
path.

### 7.2 Required preconditions

Fail closed unless all are true:

1. `worktree` is an absolute owner-created worktree path supplied by
   `owned_worktree(...)`;
2. probed `LTANAPATH` resolves/identifies that exact worktree;
3. probed `VOLATILEPATH` is absolute;
4. probed lexical `OUTPATH` is exactly the worktree-local expected form derived
   from the probed fields, i.e. the worktree `OUTPUT/Analysis/<ANATYPE>LT`
   location;
5. the configured external target `<VOLATILEPATH>/OUTPUT` exists and is a
   directory;
6. the configured selected artifact root under that target resolves exactly to
   the `outdir` supplied to this owner;
7. `<worktree>/OUTPUT` does not already exist in any form.

Use `os.path.lexists(...)` or equivalent when checking the worktree link path so
a dangling symlink cannot be silently replaced.

If any precondition fails, stop the gate with a specific fail-closed reason.
Do not delete or replace an unexpected pre-existing path.

### 7.3 Link creation and postconditions

Create exactly one symlink:

```text
<worktree>/OUTPUT -> <VOLATILEPATH>/OUTPUT
```

inside the disposable worktree.

Immediately verify:

- the worktree `OUTPUT` path is a symlink;
- it is not dangling;
- it resolves to the exact configured external OUTPUT directory;
- the probed lexical child `OUTPATH` now resolves to the exact configured
  `outdir`.

Return a small structured record suitable for status/run-summary provenance,
including at minimum:

```text
link_path
literal_target
resolved_target
resolved_child_outpath
passed
```

Do not create, delete, or modify anything in the ordinary checkout.

### 7.4 Execution order

The owner order must be:

```text
owned analysis worktree
-> copied-ltsep runtime overlay
-> existing path probes
-> existing outdir/VOLATILEPATH authority check
-> NEW debug OUTPUT-link preparation/verification
-> existing external symlink preflight
-> exactly one analysis command
```

The new link must exist **before** `run_analysis(...)`.

No second analysis command is allowed.

### 7.5 Runtime record

Persist the returned OUTPUT-link record in:

- gate status before `analysis_started=True` or in the same transition;
- the E.8.2 run summary.

The record is validation provenance only; it creates no scientific claim.

---

## 8. Do not change the launcher to solve this

Do **not** modify:

```text
run_Prod_Analysis.sh
set_SymLinks.sh
```

Do not add a production special case for this owner.

The ordinary launcher and ordinary farm checkout are established mature
interfaces. This bug exists because the new owner changed execution topology.

The repair belongs in the owner that created that topology.

---

## 9. Deterministic regression tests

Extend:

```text
testing/test_run_e8_2_left_lowe_scientific_audit_gate.py
```

with deterministic pure-Python/filesystem coverage.

At minimum cover the following.

### 9.1 Fresh-worktree positive case

Using a temporary directory:

1. create a fake owner worktree;
2. create a fake absolute `VOLATILEPATH/OUTPUT/Analysis/KaonLT`;
3. provide a path dictionary whose:
   - `LTANAPATH` is the fake worktree;
   - `VOLATILEPATH` is the fake volatile root;
   - `ANATYPE` is `Kaon`;
   - lexical `OUTPATH` is
     `<worktree>/OUTPUT/Analysis/KaonLT`;
4. call the new helper;
5. prove:
   - `<worktree>/OUTPUT` is a symlink;
   - its target is the configured external OUTPUT;
   - lexical child `OUTPATH` resolves to the external KaonLT artifact root.

### 9.2 Launcher-prefix regression

After the helper succeeds, reproduce the launcher's relevant first filesystem
operation deterministically:

```text
mkdir -p <lexical OUTPATH>
```

using Python or a local subprocess.

Prove that:

- the worktree `OUTPUT` path remains a symlink;
- the mkdir traverses the link rather than replacing it with a real local
  directory;
- the lexical child OUTPATH still resolves to the configured external artifact
  root.

This test must require no ROOT, PyROOT, farm filesystem, or production run.

### 9.3 Input visibility regression

Place a harmless fixture file under the fake external artifact root with a
Center/kaon/Q4p4W2p74-like basename.

Prove it is visible through the lexical worktree `OUTPATH`.

This is only a filesystem regression; do not fake or interpret ROOT contents.

### 9.4 Fail-closed cases

Reject at minimum:

- pre-existing real `<worktree>/OUTPUT` directory;
- pre-existing wrong symlink;
- dangling pre-existing symlink;
- wrong `LTANAPATH`;
- wrong lexical `OUTPATH`;
- relative/missing `VOLATILEPATH`;
- missing external `VOLATILEPATH/OUTPUT`;
- external artifact root not equal to supplied `outdir`.

The helper must never remove or replace these unexpected paths.

### 9.5 Owner-flow ordering

Update the owner flow/mocking test so it explicitly proves:

```text
prepare_runtime_overlay
< probe_paths
< prepare_debug_output_link
< external_symlink_preflight
< run_analysis
```

The test must no longer mock away the existence of this required new
filesystem-topology step without asserting it occurred.

### 9.6 Existing guarantees retained

All existing E.8.2 owner tests must continue to prove:

- canonical-five profile declaration remains unchanged;
- requested/packaged scope remains Left/lowe only;
- 18 E.8.2 pages are required;
- Method-A artifacts are not required;
- E.8.4 unavailable does not fail E.8.2 acceptance;
- one analysis command only;
- no canonical-five execution;
- no F.6.3/candidate materialization;
- ordinary checkout preservation is fail closed;
- installed-ltsep preservation is fail closed;
- collection and ZIP verification remain unchanged.

---

## 10. Local checks

Do not run ROOT, PyROOT, the farm, or the production analysis.

Discover the local Python interpreter according to tracked `docs/memory/TOOLS.md`.

Run at minimum:

```text
py_compile:
  testing/run_e8_2_left_lowe_scientific_audit_gate.py
  testing/test_run_e8_2_left_lowe_scientific_audit_gate.py

new/updated E.8.2 owner unit suite
existing E.8.2 source regression suite
directly relevant isolation/collector regressions
git diff --check
```

If a directly relevant deterministic isolation test already exists, run it
rather than creating redundant broad test coverage.

Report exact commands and results.

---

## 11. Memory/status update

Create:

```text
docs/memory/phases/e8-2-left-lowe-runtime-owner-output-isolation-repair.md
```

Record concisely:

- starting source identity;
- supplied failed-gate symptom;
- source-verified root cause;
- exact repair boundary;
- files changed;
- deterministic test results;
- no farm execution;
- runtime validation still pending.

Update `docs/memory/CURRENT.md` so the active state no longer claims that the
pre-repair owner is simply ready for farm validation.

Required meaning:

```text
E.8 remains ACTIVE.
The first E.8.2 isolated Left/lowe farm attempt is BLOCKED at owner runtime
isolation because the detached debug worktree exposed an empty local OUTPATH
to DiamondPlot.
Production/scientific source is unchanged.
A narrow owner-only OUTPUT-link repair is under source review / development.
Canonical-five provenance repair remains DEFERRED.
Final E.8 and F.6.4 remain BLOCKED.
```

After successful local implementation but before farm validation, use the
appropriate work-state:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

for this repaired owner only.

The sole ordinary NEXT after successful implementation must be equivalent to:

```text
NEXT — after ChatGPT actual-diff review, user commit/push and pushed-state
synchronization/farm-readiness review, rerun exactly one isolated
Q4p4W2p74 Left/lowe E.8.2 scientific-audit owner gate and return its fresh
evidence for review.
```

Do not make commit/push itself the scientific NEXT.

`CURRENT.md` is already near its 8 KiB soft threshold. Consolidate stale wording
as needed so it remains concise and preferably below 8192 bytes. Do not perform
a broad memory rewrite.

Regenerate/check `docs/memory/manifest.json` according to the tracked memory
procedure.

---

## 12. Final health/integrity gate

At the end:

1. regenerate/check the memory manifest;
2. run the ordinary memory-health command from tracked `TOOLS.md`;
3. run bootstrap JSON if required by the tracked task-final procedure;
4. run `git diff --check`;
5. confirm only the six allowed versioned paths changed;
6. confirm every frozen production/scientific path is byte-identical;
7. produce one complete `kaonlt_review.diff` in repository root containing:
   - tracked diffs;
   - complete no-index additions for new/untracked files.

Do not stage files merely to create the review bundle.

Report:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
warning classification: <none or classification>
manifest check: PASS | FAIL
bootstrap: PASS | FAIL
git diff --check: PASS | FAIL
```

---

## 13. Hard stops

Stop as `BLOCKED` rather than broadening scope if any of the following appears
necessary:

- changing `run_Prod_Analysis.sh`;
- changing `set_SymLinks.sh`;
- changing any `src/` file;
- changing the frozen collector;
- changing the E.8.2 profile;
- changing canonical-five owner/source;
- deleting/replacing an unexpected worktree path;
- mutating ordinary farm checkout state;
- mutating installed ltsep;
- introducing a second analysis execution;
- constructing a new scientific correction;
- using Method A or Method B numerically;
- running the farm before actual-diff/push/pushed-state review.

---

## 14. Acceptance criteria for this implementation

The source-changing repair is ready for ChatGPT actual-diff review only when all
of the following hold:

1. only the six allowed paths changed;
2. production/scientific source is byte-identical;
3. the E.8.2 profile and frozen collector are byte-identical;
4. the owner creates exactly one disposable-worktree OUTPUT symlink from
   validated path authority before child launch;
5. the owner proves lexical child OUTPATH resolves to configured external
   `outdir`;
6. unexpected pre-existing paths fail closed and are never removed/replaced;
7. deterministic regression reproduces the launcher's early mkdir behavior and
   proves the symlink survives;
8. owner flow proves link preparation precedes child execution;
9. existing E.8.2 scientific/page/collector boundaries remain unchanged;
10. all required local checks pass;
11. memory manifest/health/bootstrap/diff checks pass;
12. no farm execution occurred.

Status after implementation:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

This status applies to the repaired E.8.2 owner infrastructure only. It is not
E.8.2 runtime acceptance.

---

## 15. Codex stop point

Codex must stop after:

- implementation;
- deterministic local checks;
- memory update;
- manifest/health/bootstrap checks;
- complete `kaonlt_review.diff` creation.

Codex must **not**:

- stage merely for review;
- commit;
- push;
- update refs;
- run the farm;
- run ROOT/PyROOT production analysis.

Return:

- starting/ending branch;
- starting/ending HEAD and local `origin/test`;
- exact changed paths;
- concise implementation summary;
- exact local checks and results;
- memory-health report;
- `kaonlt_review.diff` size/path;
- explicit confirmation that no farm command was run.
