# E.8.4 / F.6.3 current-lineage canonical-five validation — isolated-runtime orchestration task contract

## 1. Task class and objective

This is a **narrow source-changing operational-readiness task** under
`docs/memory/CODEX.md`.

It supersedes the prior canonical-five orchestration drafts that stopped at
hard gates. Those blocks were correct and exposed an incomplete runtime-mutation
audit.

The exact live starting identity is:

```text
branch: test
HEAD: ed378e0f30a357c6293da64f8f8f3bfbe187d1fe
origin/test: ed378e0f30a357c6293da64f8f8f3bfbe187d1fe
```

The scientific Method-A implementation is already established. This task must
not redesign Method A, F.3, F.4, F.6.2, F.6.3, E.8.4, baseline pion subtraction,
or production physics.

### Objective

Add one tracked canonical-five validation owner/profile that runs the existing
ordinary full Q4p4W2p74 kaon analysis **inside an owner-created detached
temporary worktree**, so all repository-local runtime mutations are isolated
from the user's ordinary Jefferson Lab checkout.

Do **not** add a `-v` launcher mode.
Do **not** modify `run_Prod_Analysis.sh`.
Do **not** modify `set_SymLinks.sh`.
Do **not** modify `src/setup/set_sig_fortran.py`.

The existing full launcher remains the analysis producer. The owner supplies
the isolation boundary.

---

## 2. Why the previous contract was blocked

The prior source audit stopped too early at `git clean -fdx`.

Current source proves the ordinary full launcher is mutation-bearing through at
least these paths:

1. `run_Prod_Analysis.sh` executes `git clean -fdx` on the current worktree for
   ordinary non-iteration/non-automated analysis.
2. It invokes `./set_SymLinks.sh`, which creates/replaces repository symlinks and
   can reach through `simc_gfortran` / background-SIMC links.
3. For low epsilon it invokes:

   ```text
   python3 set_sig_fortran.py ${Q2} ${W} ${ParticleType} ${POL}
   ```

4. `src/setup/set_sig_fortran.py` opens and rewrites:

   ```text
   ${LTANAPATH}/src/models/xmodel_${ParticleType}_${pol_str}.f
   ```

   so Q4p4W2p74 kaon/+1 rewrites:

   ```text
   src/models/xmodel_kaon_pl.f
   ```

5. Accepted farm provenance has also observed:

   ```text
    M src/models/xmodel_kaon_pl.f
   ?? src/kaon/functions/Q4p4W2p74.model
   ```

Therefore a primary-checkout "skip git clean" exception is insufficient.

### Durable lesson

For farm-bound orchestration, recursively audit the runtime mutation graph,
including called scripts and symlink-target side effects. When the established
analysis legitimately mutates repository-local runtime files, isolate the
entire execution in a bounded disposable worktree rather than trying to preserve
the primary checkout command-by-command.

---

## 3. Exact start gate

Before editing, report:

```bash
git branch --show-current
git rev-parse HEAD
git rev-parse origin/test
git status --short --untracked-files=all
```

Required:

```text
branch = test
HEAD = ed378e0f30a357c6293da64f8f8f3bfbe187d1fe
origin/test = ed378e0f30a357c6293da64f8f8f3bfbe187d1fe
```

The local repository may already contain this uncommitted contract path because
the user moved it into place for Codex. That is expected.

Hard stop if:

- branch/HEAD/origin differ;
- unrelated local state cannot be preserved;
- implementation requires scientific-source changes;
- implementation requires changing the ordinary launcher or its physics;
- implementation requires cleanup/reset/stash/checkout of the user's ordinary
  worktree.

Do not reset, clean, stash, stage, commit, push, update refs, or run the farm.

---

## 4. Mandatory startup reads

Read in exact order:

1. `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read:

- `docs/memory/MAINTENANCE.md`
- `docs/memory/CODEX.md`
- `docs/memory/COMMUNICATION.md`
- `docs/memory/LEARNINGS.md`
- `docs/memory/TOOLS.md`
- `docs/memory/templates/CODEX_CONTRACT.md`

Then task-relevant source/records:

- `docs/memory/evidence/e8-4-left-lowe-method-a-current-lineage-diagnostic-runtime-closure-2026-10-04.md`
- `docs/memory/evidence/f6-2-scientific-runtime-closure.md`
- `docs/memory/evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md`
- `docs/memory/evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md`
- `docs/memory/roadmap/STATUS.md`
- `run_Prod_Analysis.sh`
- `set_SymLinks.sh`
- `src/setup/set_sig_fortran.py`
- `farm_env/print_ltsep_path_fields.py`
- `farm_env/ltsep_paths.py`
- `background_samples/background_samples.conf`
- `testing/run_e8_4_fix5_left_lowe_plot_gate.py`
- `testing/test_run_e8_4_fix5_left_lowe_plot_gate.py`
- `testing/pion_hgcer_validation_bundle_profile_e8_4_fix5.json`
- `testing/test_pion_hgcer_validation_bundle_profile_e8_4_fix5.py`
- `testing/collect_pion_hgcer_validation_bundle.py`
- `testing/test_collect_pion_hgcer_validation_bundle.py`

Do not reopen closed F stages.

---

## 5. Scientific/authority state to preserve

`SOURCE VERIFIED`:

- accepted F.1 scientific/runtime authority already spans the canonical five:
  - Left / lowe
  - Center / lowe
  - Left / highe
  - Center / highe
  - Right / highe
- there is no canonical Right / lowe setting;
- current-baseline F.2/F.3 science is unchanged from accepted authority;
- F.4 is the first scientifically changed stage;
- Method A remains detached/non-production;
- Method B remains diagnostic/cross-check only and numerically excluded;
- baseline `w0` owns the accepted pion-control -> kaon-background transfer;
- F.4 owns signed canonical-t parent preservation;
- no new response map, correction, normalization, cut, PID rule, MM rule,
  acceptance model, threshold, or physics quantity is required.

This is orchestration/validation work only.

---

## 6. Correct runtime isolation architecture

### 6.1 Primary checkout

The user's ordinary farm checkout is an **input/source-control anchor only**.

The canonical-five analysis subprocess must never run with the ordinary checkout
as its cwd/LTANAPATH.

Before creating the temporary analysis worktree, the owner must record:

- branch;
- HEAD;
- local `origin/test`;
- complete `git status --porcelain=v1 --untracked-files=all`;
- exact byte SHA-256 for each allowed dirty farm-local source/output file that
  exists, including at minimum:
  - `src/models/xmodel_kaon_pl.f`
  - `src/kaon/functions/Q4p4W2p74.model`
- `git diff --check HEAD --` over gate-owned source, excluding only the already
  accepted farm-output paths.

The owner must never call destructive Git operations on the ordinary checkout.

### 6.2 Dedicated analysis worktree

Create one owner-bounded temporary directory and add:

```text
<temporary-root>/analysis
```

with:

```text
git worktree add --detach <temporary-root>/analysis <source_commit>
```

Require before execution:

- exact `HEAD == source_commit`;
- detached branch;
- clean worktree;
- path is not the ordinary checkout;
- parent temporary directory is owner-created.

All full analysis execution occurs there.

The existing ordinary command is:

```text
./run_Prod_Analysis.sh 4p4 2p74
```

No `-d`.
No new `-v`.
No launcher modification.

The existing launcher therefore retains its normal paired prepass, full low
analysis and full high analysis.

### 6.3 LTANAPATH isolation proof

Before analysis, the owner must source-prove/runtime-probe path resolution from
the detached worktree.

Using the tracked `farm_env/print_ltsep_path_fields.py`, resolve path fields for
at least:

```text
<analysis-worktree>
<analysis-worktree>/src/setup/set_sig_fortran.py
<analysis-worktree>/set_SymLinks.sh
```

For each probe, require the resolved `LTANAPATH` to normalize exactly to the
detached analysis worktree.

If any probe resolves to the ordinary checkout or another repository tree,
return `BLOCKED` before running analysis.

This gate exists specifically because `set_sig_fortran.py` writes through its
own `LTANAPATH`.

### 6.4 External-symlink no-mutation preflight

`set_SymLinks.sh` is allowed to create/replace links **inside the disposable
analysis worktree**.

It is not allowed to repair/mutate external SIMC-tree links during this gate.

Before launching, inspect the resolved external paths and prove that
`set_SymLinks.sh` will be a no-op for link objects outside the disposable
analysis worktree.

At minimum:

#### Primary SIMC tree

For:

```text
${SIMCPATH}/OUTPUTS
${SIMCPATH}/input
${SIMCPATH}/worksim
```

require each path to:

- already be a symlink;
- resolve to an existing target.

Because current `set_SymLinks.sh` changes those only when the symlink is
missing/broken, this guarantees no external repair there.

#### Background SIMC tree

Resolve `BACKGROUND_SIMC_PATH` from the tracked
`background_samples/background_samples.conf` in the detached worktree.

If configured and the directory exists, require:

```text
${BACKGROUND_SIMC_PATH}/worksim
```

to already be a symlink whose **literal readlink target** is exactly the target
current `set_SymLinks.sh` would request:

```text
${VOLATILEPATH}/worksim/
```

and require that target to exist.

If this condition is not satisfied, return `BLOCKED` before analysis rather
than allowing `sync_symlink` to replace it.

No owner code may "fix" these external links.

### 6.5 Runtime mutations inside the temporary worktree

Once all isolation gates pass, the following are allowed **only inside the
owner-created detached analysis worktree**:

- `git clean -fdx` from the unchanged ordinary launcher;
- temporary repository-local symlink creation by `set_SymLinks.sh`;
- `src/models/xmodel_kaon_pl.f` rewrite by `set_sig_fortran.py`;
- any generated Q4p4W2p74 repository-local model/cache file produced by the
  established analysis;
- ignored Python/cache files.

These are runtime working-state mutations, not source changes.

Do not copy them back into the ordinary checkout.
Do not stage them.
Do not use them as new source authority.

### 6.6 Ordinary checkout preservation after analysis

After the subprocess returns—whether success or failure—the owner must inspect
the ordinary checkout again.

Require:

- branch unchanged;
- HEAD unchanged;
- local `origin/test` unchanged;
- complete porcelain status exactly equal to the pre-run snapshot;
- every pre-snapshotted allowed dirty farm-local file still exists with the
  exact same byte SHA-256;
- no new gate-owned source path appeared.

A mismatch is a failed preservation gate. Stop before collection/acceptance.

Do not restore the ordinary checkout automatically; preservation must be proven,
not repaired after the fact.

---

## 7. Bounded cleanup rule

The owner must remove only worktrees/directories it created.

For the detached analysis worktree, cleanup may use:

```text
git worktree remove --force <exact-owner-created-analysis-worktree>
```

because the worktree is disposable and expected to contain runtime mutations.

Before force removal require:

- the path is under the unique owner-created temporary parent;
- the path is not equal to or inside the ordinary checkout;
- the temporary parent is not the ordinary checkout.

The same bounded policy applies to collector source-check worktrees.

No `git clean`, `git reset`, `git checkout`, `git stash`, or worktree removal
may target the ordinary checkout.

---

## 8. Canonical-five owner

Create:

```text
testing/run_e8_4_fix5_canonical_five_plot_gate.py
```

Use the accepted Left/lowe owner architecture where applicable, but make this
owner responsible for the complete canonical-five operation:

```text
ordinary-checkout preflight/snapshot
-> collector source preflight
-> create detached analysis worktree
-> LTANAPATH isolation probes
-> external-symlink no-mutation preflight
-> run ordinary full Q4p4W2p74 launcher in detached analysis worktree
-> verify low+high completion markers
-> verify ordinary checkout preservation
-> verify exact five-setting artifacts/pages
-> write neutral canonical-five owner summary
-> generic five-setting collection from clean detached collector source
-> final source/ordinary-checkout recheck
-> ZIP verification
-> bounded temporary-worktree cleanup
```

Success stdout must be one ZIP path.
Diagnostics/logging go to stderr / owner status artifacts as in the accepted
owner pattern.

### Exact settings

Exactly:

```text
Left / lowe
Center / lowe
Left / highe
Center / highe
Right / highe
```

Reject:

- Right / lowe;
- duplicates;
- missing settings;
- extra settings;
- unknown tokens.

---

## 9. Analysis completion gate

The subprocess command must be exactly:

```text
./run_Prod_Analysis.sh 4p4 2p74
```

executed with:

```text
cwd = detached analysis worktree
```

Require return code zero.

Require the log to contain the ordinary full-analysis completion markers:

```text
Low Epsilon Completed!
High Epsilon Completed!
```

Do not accept the Left/lowe debug markers as canonical-five success.

Do not use `|| true`.

A child analysis failure stops before collection.

---

## 10. Dedicated canonical-five profile

Create:

```text
testing/pion_hgcer_validation_bundle_profile_e8_4_canonical_five.json
```

Use the existing compatible generic bundle-profile schema and
`collection_mode = generic_artifacts`.

Exact setting inventory is the canonical five above.

Global artifacts:

- accepted current-lineage candidate F.3;
- accepted current-lineage candidate F.4;
- new setting-neutral canonical-five run summary.

Per-setting artifacts preserve the accepted E.8.4 Fix.5 inventory:

- procedure PDF;
- page manifest;
- full-analysis JSON;
- correction-ledger JSON;
- correction-ledger CSV.

Do not reuse the Left/lowe owner-summary name.

Do not mutate the historical Left/lowe profile.

Use the existing effective-profile source-pin pattern so the future farm bundle
is pinned to the exact reviewed/pushed SHA.

---

## 11. Five-setting artifact/page verification

For each setting require freshness and exact setting identity.

At minimum preserve existing accepted E.8.4 page semantics:

- procedure PDF valid/nonempty;
- page manifest valid/nonempty;
- renderer failures empty;
- required E.8.2/E.8.3/E.8.4 page families present;
- new E.8.4 per-t pages have exact t identity;
- represented phi child inventory is complete under the existing accepted rule;
- setting-scope parent-closure page exists;
- full-analysis JSON fresh;
- correction ledgers fresh;
- candidate F.3/F.4 identities exact;
- Method B numerical dependency absent;
- detached branch does not mutate baseline production objects;
- parent-preserving closure uses the existing accepted tolerance/semantics.

Do not invent a new threshold.

A single failed setting fails the whole canonical-five gate.

---

## 12. Canonical-five owner summary

Create a neutral summary artifact:

```text
Q4p4W2p74_e8_4_fix5_canonical_five_run-summary.json
```

with a new explicit schema.

Record:

- source SHA;
- exact five-setting inventory;
- analysis command;
- detached analysis-worktree identity;
- LTANAPATH isolation-probe results;
- external-symlink no-mutation preflight;
- ordinary-checkout pre-run status/hashes;
- ordinary-checkout post-run preservation result;
- analysis return code;
- log SHA-256;
- per-setting artifact/page verification;
- candidate F.3/F.4 identities;
- collector source checks;
- ZIP identity on success;
- explicit boundaries:
  - Method A detached/non-production;
  - Method B numerically excluded;
  - no production promotion;
  - no absolute pion mis-ID claim;
  - no absolute-SIMC amplitude claim.

Do not claim the ordinary farm checkout was clean if it was not.

---

## 13. Deterministic tests

Create:

```text
testing/test_run_e8_4_fix5_canonical_five_plot_gate.py
testing/test_pion_hgcer_validation_bundle_profile_e8_4_canonical_five.py
```

### Owner tests must cover

- exact ordinary launcher command with **no flag**;
- analysis cwd is the detached analysis worktree, never the ordinary checkout;
- exact five-setting inventory;
- Right/lowe rejection;
- wrong branch/HEAD/origin blocks before worktree creation/analysis;
- unexpected dirty gate source blocks before analysis;
- accepted farm-local dirty paths are recorded, not cleaned;
- analysis worktree is detached/exact/clean before execution;
- LTANAPATH resolving to ordinary checkout blocks before analysis;
- set-sig caller LTANAPATH resolving outside analysis worktree blocks;
- missing/broken `${SIMCPATH}/OUTPUTS|input|worksim` blocks before analysis;
- mismatched background-SIMC `worksim` target blocks before analysis;
- no external symlink repair is performed by owner;
- synthetic runtime mutation of analysis-worktree
  `src/models/xmodel_kaon_pl.f` does not change the primary fixture;
- primary status/hash mismatch after analysis blocks before collection;
- analysis failure blocks before collection;
- missing low/high completion markers block;
- stale/missing setting artifact blocks;
- page-manifest setting mismatch blocks;
- renderer failure blocks;
- page/child inventory failure blocks;
- candidate authority mismatch blocks;
- collector failure blocks;
- ZIP verification failure blocks;
- worktree cleanup is path-bounded to owner-created worktrees;
- owner never calls destructive Git operations on the ordinary checkout;
- successful synthetic flow has exact ordering:
  `source preflight -> isolated analysis -> preservation check ->
  five-setting verify -> collection -> source recheck -> ZIP verify`.

### Profile tests must cover

- exact five settings;
- no Right/lowe;
- neutral canonical-five summary;
- expected global/per-setting artifact keys;
- compatibility with unchanged generic collector;
- no Left/lowe summary dependency.

### Existing regressions

Run unchanged:

```text
testing/test_run_e8_4_fix5_left_lowe_plot_gate.py
testing/test_pion_hgcer_validation_bundle_profile_e8_4_fix5.py
testing/test_collect_pion_hgcer_validation_bundle.py
testing/test_f6_3_parallel_full_procedure_method_a.py
testing/test_e8_4_production_impact_audit.py
testing/test_e8_4_fix5_main_order.py
testing/test_run_prod_analysis_debug_left_low.py
```

No test may run the real farm analysis.

---

## 14. Explicitly frozen source

Do not modify:

```text
run_Prod_Analysis.sh
set_SymLinks.sh
src/
farm_env/
testing/run_e8_4_fix5_left_lowe_plot_gate.py
testing/pion_hgcer_validation_bundle_profile_e8_4_fix5.json
testing/collect_pion_hgcer_validation_bundle.py
```

Do not change any scientific or production logic.

The existing ordinary full launcher is the producer being isolated and tested,
not redesigned.

---

## 15. Allowed versioned paths

Only:

```text
testing/run_e8_4_fix5_canonical_five_plot_gate.py
testing/test_run_e8_4_fix5_canonical_five_plot_gate.py
testing/pion_hgcer_validation_bundle_profile_e8_4_canonical_five.json
testing/test_pion_hgcer_validation_bundle_profile_e8_4_canonical_five.py

docs/memory/phases/e8-4-method-a-current-lineage-canonical-five-validation-orchestration-task-contract.md
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/LEARNINGS.md
docs/memory/roadmap/STATUS.md
docs/memory/manifest.json
```

If any other versioned path is required, return `BLOCKED`.

---

## 16. Memory updates

### CURRENT.md

After successful local implementation/checks, record:

- user authorization lifted canonical-five `DEFERRED`;
- scientific authority required no Method-A redesign;
- audit found ordinary full runtime mutates repository-local state through
  cleanup, symlink setup and `set_sig_fortran`;
- canonical-five orchestration now isolates the unchanged full analysis in a
  detached temporary worktree;
- ordinary checkout preservation is a mandatory gate;
- canonical-five is not runtime validated yet;
- Method A remains detached/non-production;
- Method B remains diagnostic only;
- final E.8 and F.6.4 remain `BLOCKED`;
- absolute-SIMC interpretation remains separately `BLOCKED`.

Push-stable NEXT:

```text
NEXT — after ChatGPT actual-diff review, user commit/push, and pushed-state
synchronization, perform the farm-readiness audit of the tracked isolated
canonical-five owner/profile. Only a PASS may authorize one user-run
Q4p4W2p74 canonical-five farm gate.
```

### MEMORY.md

Add durable operational knowledge only:

- the ordinary full launcher is intentionally mutation-bearing;
- canonical validation must not run it in the ordinary farm checkout;
- detached-worktree isolation preserves production semantics while containing
  `git clean`, model rewrites, generated repository-local runtime files and
  local symlink creation;
- path isolation and external-symlink no-mutation checks precede execution.

### LEARNINGS.md

Record the user-identified workflow failure:

- do not stop a farm-orchestration audit at top-level destructive commands;
  recursively trace invoked scripts and symlink targets;
- a preservation contract must audit the full runtime mutation graph before
  implementation;
- when legitimate runtime mutation is broad, isolate the runtime as a unit
  instead of adding successive primary-checkout exceptions.

### roadmap/STATUS.md

Canonical-five moves from user-deferred to active implementation/source state
only as warranted.

Do not mark runtime closure or promotion.

---

## 17. Local validation and memory health

Discover `<PYTHON>` by repository convention.

Run the two new test files plus all regressions in Section 13.

Then:

```bash
<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/check_memory_health.py --root .
<PYTHON> -B tools/memory_bootstrap.py --root . --json
git diff --check
```

Report exact:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
warning classification: <none or each warning>
manifest check: PASS | FAIL
bootstrap: PASS | FAIL
git diff --check: PASS | FAIL
```

Local checks do not establish ROOT/farm validation.

---

## 18. Future farm chain — not authorized by this task

After successful actual-diff review, user push and pushed-state review:

```text
input authority:
  accepted current-lineage F.1/F.3/F.4

producer:
  unchanged run_Prod_Analysis.sh 4p4 2p74
  executed only in owner-created detached analysis worktree

path/preservation checker:
  canonical-five owner LTANAPATH + external-symlink + primary-preservation gates

artifact/page checker:
  canonical-five owner exact five-setting verification

collector:
  unchanged collect_pion_hgcer_validation_bundle.py
  with dedicated canonical-five effective profile

invocation owner:
  testing/run_e8_4_fix5_canonical_five_plot_gate.py

returned artifact:
  one fresh canonical-five ZIP plus owner status/summary/log
```

No interactive multi-step farm sequence is allowed.

No farm command is authorized by this contract.

---

## 19. Review bundle

Create repository-root temporary:

```text
kaonlt_review.diff
```

Include complete tracked diffs and complete additions for every new allowlisted
file.

Do not stage merely for review.

Stop before commit/push.

---

## 20. Acceptance criteria

PASS only if:

- start identity exact;
- only allowlisted files change;
- no launcher/source/farm_env changes;
- ordinary launcher command remains unchanged;
- analysis runs only in detached temporary worktree in the owner design;
- LTANAPATH isolation is explicitly checked for launcher/set-sig contexts;
- external SIMC symlink mutation is fail-closed before analysis;
- primary checkout exact status/byte preservation is checked after analysis;
- owner-created worktree cleanup is path-bounded;
- exact canonical-five setting/profile/artifact verification exists;
- all deterministic/regression/memory checks pass;
- complete review bundle exists;
- no farm, stage, commit or push occurred.

If all pass, report only for the orchestration candidate:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

Canonical-five runtime itself remains unvalidated.

---

## 21. Hard stops

Return `BLOCKED` if:

- exact identity fails;
- LTANAPATH cannot be proven isolated to the detached analysis worktree;
- any external symlink would be repaired/mutated by `set_SymLinks.sh`;
- ordinary checkout preservation cannot be guaranteed and checked;
- analysis requires copying the user's dirty model state into the detached
  worktree;
- analysis requires a source path outside the allowlist to change;
- generic collector cannot support the dedicated profile unchanged;
- a new physics threshold/calculation is required;
- memory/manifest hard checks fail.

Do not invent another primary-checkout exception.
Do not repair the ordinary checkout after a failed preservation check.
Do not run the farm.
