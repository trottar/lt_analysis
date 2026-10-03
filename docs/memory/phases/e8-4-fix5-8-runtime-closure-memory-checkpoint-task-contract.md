# KaonLT — Fix.5.8 runtime-closure memory checkpoint and consistency cleanup

## 1. Task class and objective

This is a **closure task** under `docs/memory/CODEX.md`.

It is memory/evidence-only. Do **not** modify scientific, analysis, renderer,
owner, collector, launcher, profile, test, or farm-execution source.

Objective:

1. record the accepted `Q4p4W2p74 / Left / lowe` Fix.5.8 farm evidence;
2. close Fix.5.8 narrowly as `CLOSED / RUNTIME VALIDATED`;
3. reconcile CURRENT, durable MEMORY, roadmap, validation-history, phase, and
   investigation records with that accepted evidence;
4. remove consumed/stale "pending actual-diff / pending push / pending farm /
   repaired rendering NOT VERIFIED" wording where it is now false;
5. preserve all remaining scientific and scope boundaries exactly;
6. regenerate/check the memory manifest and run the full ordinary memory-health
   gate;
7. stop for independent ChatGPT actual-diff review before user commit/push.

No farm run is authorized by this task.

---

## 2. Exact starting source

Required branch:

```text
test
```

Required starting local HEAD and local `origin/test`:

```text
2ddeab47d55edb57d2f313022a948c4376730c19
```

Commit subject:

```text
E8.4 Fix.5.8: repair presentation legibility
```

Before editing, establish and report:

```bash
git branch --show-current
git rev-parse HEAD
git rev-parse refs/remotes/origin/test
git status --short --untracked-files=all
```

Hard requirements:

- branch exactly `test`;
- HEAD exactly the full SHA above;
- local `origin/test` exactly the same SHA;
- preserve unrelated worktree state;
- this task-contract file may be the expected untracked input;
- any unrelated tracked modification is a hard stop;
- do not reset, clean, stash, checkout over, commit, push, update remote refs,
  or run the farm.

---

## 3. Mandatory startup order

Read in full and in exact order:

1. repository-root `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read only task-relevant records:

- `docs/memory/MAINTENANCE.md`
- `docs/memory/CODEX.md`
- `docs/memory/COMMUNICATION.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/evidence/VALIDATION_HISTORY.md`
- `docs/memory/evidence/e8-4-fix5-7-left-lowe-runtime-and-fix5-visual-gate-2026-10-03.md`
- `docs/memory/phases/e8-4-fix5-8-presentation-legibility-repair.md`
- `docs/memory/investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md`
- the tracked memory manifest.

Do not reopen unrelated historical records unless one of the above has a direct
reference necessary to resolve a concrete inconsistency.

---

## 4. Accepted fresh farm evidence to record

The user supplied and ChatGPT independently inspected:

```text
KaonLT_E8_4_Fix5_Left_lowe_Q4p4W2p74_20261003-032934-987132.zip
```

Exact ZIP properties:

```text
bytes: 30851374
sha256: 966e667b36b2626b5c16f00fc684a099a5d24e5d2015bb547612a1e63eb2a25e
```

Scope:

```text
Q4p4W2p74 / Left / lowe
```

Farm-evaluated source:

```text
2ddeab47d55edb57d2f313022a948c4376730c19
```

The independently reviewed returned bundle establishes:

- ZIP integrity PASS;
- bundle schema `pion_hgcer_validation_bundle/v4`;
- `complete=true`;
- `errors=[]`;
- `git_head == required_analysis_commit ==
  2ddeab47d55edb57d2f313022a948c4376730c19`;
- requested setting exactly `Left / lowe`;
- validation profile
  `phase_e8_4_fix5_left_lowe_shareable_pages/v1`;
- all eight declared artifact entries match the manifest byte counts and
  SHA-256 values;
- analysis return code zero;
- procedure PDF page count exactly 97;
- page-manifest `renderer_failures=[]`;
- full seven-field page-manifest producer provenance is present;
- source/checker tests embedded in the returned validation record returned zero.

### Targeted rendered-page review

The contracted Fix.5.8 visual defects are repaired in the fresh farm PDF:

- pages 67, 69, 71:
  E.8.3 t-phi title is ASCII-safe; no `â€”` mojibake;
- page 73:
  all seven non-production flag name/value pairs are visibly present without
  right-edge clipping;
- pages 80, 87, 94:
  all nine phi children and all six stored support quantities per child are
  visibly present, with no right-edge clipping and no overlap with the upper
  yield panels;
- page 96:
  E.8.4 setting-summary title is ASCII-safe; no mojibake;
- full 97-page scan found no new black-square renderer failure, broken target
  glyph, or obvious collateral clipping;
- previously repaired baseline-versus-Method-A visibility remains intact on
  representative t-bin comparison pages.

This establishes repaired ROOT/PDF legibility for this one accepted setting.

### Preserved scientific boundaries in the fresh bundle

The unresolved absolute-SIMC blocker remains explicit:

```text
SIMC_normfac_luminosity_and_charge_units_not_source_proven
```

The six absolute-SIMC pages remain unavailable by design. Do not recast this as
a Fix.5.8 failure or as evidence of incorrect SIMC normalization.

Current F.4 candidate SHA-256 remains:

```text
1d545924eba89c7f9ffa28028e307aca9b434a89beec06863cf2893887b6b902
```

Stored parent closure remains passed for the three canonical t parents:

```text
t1:
  baseline = 0.1081517597541286
  adjusted = 0.10815175975412862
  residual = 1.3877787807814457e-17

t2:
  baseline = 0.3374726818838243
  adjusted = 0.3374726818838243
  residual = 0

t3:
  baseline = 0.10645779207740595
  adjusted = 0.10645779207740595
  residual = 0
```

The candidate still records:

```text
production_application_performed = false
production_objects_mutated = false
correction_applied_to_production = false
method_b_numerical_dependency = false
```

Therefore:

- Method A remains detached/non-production;
- Method B remains diagnostic/cross-check only and numerically excluded;
- no production promotion follows;
- no canonical-five closure follows;
- no new cancellation or SIMC amplitude interpretation follows.

The raw owner-status sidecar and raw analysis log are not required for this
closure because the tracked owner reached a complete, internally verified ZIP
and the accepted question here is the returned structural/visual gate. Do not
invent facts about sidecar/log contents.

---

## 5. Required work-state reconciliation

Use only repository-defined work-state labels.

Record exactly these scope-limited consequences:

- E.8 remains `ACTIVE`.
- Fix.5.7 remains `CLOSED / RUNTIME VALIDATED` only for its narrow
  `Q4p4W2p74 / Left / lowe` owner/checker setting-provenance repair.
- Fix.5.8 becomes `CLOSED / RUNTIME VALIDATED` only for its narrow
  `Q4p4W2p74 / Left / lowe` presentation-legibility repair.
- The narrow Fix.5 Left/lowe owner + visual presentation gate is accepted.
- F.6.3 and the existing E.8.4 Left/lowe production-impact scope retain their
  prior `CLOSED / RUNTIME VALIDATED` status; do not broaden their scope.
- E.8.1 canonical-five expansion remains `DEFERRED` by the existing user
  decision.
- Final canonical-five E.8 remains `BLOCKED`.
- F.6.4 remains `BLOCKED`; no automatic Method-A promotion.
- Absolute-SIMC interpretation remains `BLOCKED` for the existing unit-
  provenance reason only.

Do not mark all E.8 closed. Do not mark Method A production-ready. Do not
promote Method B. Do not infer runtime acceptance for the other four canonical
settings.

---

## 6. Exact stale-memory findings to clean up

The current repository memory is intentionally stale because the accepted farm
artifact arrived after commit `2ddeab47d...`.

Repair these concrete inconsistencies only.

### `docs/memory/CURRENT.md`

Remove consumed wording that still says:

- Fix.5.8 is farm-validation pending;
- repaired rendering is NOT VERIFIED;
- actual-diff review is pending;
- commit/push/pushed-state synchronization is pending;
- the narrow Left/lowe farm run is still NEXT.

Replace it with the accepted runtime closure and the new push-stable NEXT
defined below.

Keep CURRENT concise. It is already below its 8 KiB soft limit; do not expand
it back above the warning threshold unless unavoidable.

### `docs/memory/MEMORY.md`

Target only the durable statements now made false by the accepted bundle:

- the statement that actual Fix.5.5/Fix.5 visual legibility still awaits a farm
  PDF;
- the statement that real PDF/farm acceptance remains pending.

Replace them with a durable, narrow statement that the later Fix.5.8 farm gate
runtime-validates the Left/lowe presentation/owner result while canonical-five
closure and absolute-SIMC interpretation remain unresolved.

Do not add ephemeral CURRENT/NEXT prose to MEMORY.

### `docs/memory/phases/e8-4-fix5-8-presentation-legibility-repair.md`

Update status from:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

to the narrow runtime closure.

Add a concise closure section with:

- exact farm source SHA;
- exact returned ZIP name/bytes/SHA-256;
- structural PASS;
- targeted page PASS;
- no baseline-versus-Method-A regression;
- preserved SIMC blocker and production boundaries;
- no canonical-five or F.6.4 closure.

Replace its consumed pre-farm NEXT/hard-stop wording with historical closure
language. This phase record must no longer present the farm run as future work.

### `docs/memory/investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md`

At the top checkpoint only, mark the listed presentation defects as resolved by
the accepted Fix.5.8 bundle.

Preserve the historical investigation body and original hypotheses/findings.
Do not rewrite old observations as if they were made after Fix.5.8.

### `docs/memory/roadmap/STATUS.md`

The current Fix.5 region is stale. Reconcile narrowly:

- preserve historical Fix.5/Fix.5.4/Fix.5.5/Fix.5.6 chronology;
- add/record Fix.5.7 and Fix.5.8 narrow runtime closure;
- remove the obsolete dependency chain that says Fix.5.6 review/push/run is
  still pending;
- replace the stale statement that the Fix.5 shareable interpretation still
  awaits fresh Fix.5.4/Fix.5.5 visual evidence;
- keep E.8 overall `ACTIVE`;
- keep E.8.1 canonical-five expansion `DEFERRED`;
- keep final E.8 and F.6.4 `BLOCKED`.

Do not convert the roadmap into CURRENT or give it sole NEXT ownership.

### `docs/memory/evidence/VALIDATION_HISTORY.md`

This file contains old "current downstream gate" and non-closure statements
that predate the accepted F-stage/E.8 evidence.

Clean it up narrowly so it functions as a validation/evidence history rather
than a competing stale live-state record:

- add the accepted Fix.5.7/Fix.5.8 Left/lowe runtime closure entries;
- remove or explicitly historicalize obsolete "ACTIVE" claims that contradict
  current accepted evidence;
- replace the stale current-gate section with a reference that CURRENT owns the
  live NEXT;
- preserve historical evidence entries and their original scope.

Do not make this file a second CURRENT.

---

## 7. New canonical evidence record

Create:

```text
docs/memory/evidence/e8-4-fix5-8-left-lowe-runtime-closure-2026-10-03.md
```

It must include:

1. exact source identity;
2. exact ZIP filename, byte size, and SHA-256;
3. structural/bundle checks;
4. 97-page / renderer-failure result;
5. exact targeted visual-page acceptance;
6. no-regression statement for baseline-versus-Method-A overlays;
7. absolute-SIMC blocker remains;
8. current F.4 candidate identity and the three parent-closure values;
9. non-production/Method-B boundaries;
10. exact narrow status consequences;
11. explicit non-closures:
    - other four canonical settings;
    - canonical-five E.8;
    - F.6.4;
    - Method-A production promotion;
    - absolute-SIMC amplitude interpretation.

Use evidence labels exactly where helpful:
`SOURCE VERIFIED`, `RUNTIME VERIFIED`, `MEMORY ONLY`, `INFERENCE`,
`NOT VERIFIED`.

Do not embed workstation or farm filesystem paths beyond artifact names already
present in accepted evidence.

---

## 8. New authoritative NEXT

After this closing memory checkpoint is independently actual-diff reviewed,
user-committed/pushed, and pushed-state synchronized, the substantive NEXT is a
**decision gate**, not another automatic farm run:

```text
NEXT — resolve the existing user-deferred canonical-five E.8/F.6.3 expansion
decision before any additional E.8 farm work. If the deferral is explicitly
lifted, create a dedicated contract for the remaining four canonical settings
and re-establish the exact validation scope before issuing any farm command.
If the deferral remains, keep it DEFERRED and select the next approved roadmap
item. This checkpoint authorizes no farm run by itself.
```

CURRENT may shorten the wording, but it must preserve those semantics and remain
push-stable.

Do not silently lift the canonical-five deferral.

---

## 9. Allowed files

Expected tracked modifications:

```text
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/roadmap/STATUS.md
docs/memory/evidence/VALIDATION_HISTORY.md
docs/memory/phases/e8-4-fix5-8-presentation-legibility-repair.md
docs/memory/investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md
docs/memory/manifest.json
```

Expected new files:

```text
docs/memory/evidence/e8-4-fix5-8-left-lowe-runtime-closure-2026-10-03.md
docs/memory/phases/e8-4-fix5-8-runtime-closure-memory-checkpoint-task-contract.md
```

No other file is authorized without a hard stop and explicit review.

---

## 10. Frozen files

All non-memory source is frozen, including but not limited to:

```text
src/
testing/
tools/
farm_env/
run_Prod_Analysis.sh
```

Also freeze:

```text
AGENTS.md
docs/memory/USER.md
docs/memory/handoffs/CURRENT_HANDOFF.md
docs/memory/CODEX.md
docs/memory/MAINTENANCE.md
docs/memory/COMMUNICATION.md
```

unless the memory-health checker exposes a hard integrity failure that cannot
be resolved inside the allowlisted files. In that case, stop rather than
broadening scope.

---

## 11. Scientific ownership

This task changes no physics and no executable behavior.

Preserve exactly:

- random/dummy subtraction;
- slow-proton subtraction;
- pion subtraction;
- active `no_empirical_residual` profile;
- F.6.3 private-branch semantics;
- Method-A candidate values/fingerprints;
- Method-B diagnostic-only role;
- SIMC production/normalization;
- yields and uncertainties;
- cuts, templates, priors, binning, efficiencies, acceptance;
- L/T separation and cross sections;
- all current runtime artifacts and hashes.

This task records evidence; it does not reinterpret or recompute it.

---

## 12. Required consistency audit

Before finalizing edits, search the allowlisted memory files for stale phrases
and semantic equivalents of:

```text
Fix.5.8 ... FARM VALIDATION PENDING
repaired rendering remains NOT VERIFIED
actual-diff review remains pending
pushed-state synchronization ... pending
Fix.5.6 ... next farm gate
Fix.5 shareable interpretation awaits ... visual farm evidence
current downstream gate ... E.8 S1
fresh F.1.Fix.5 farm gate
```

Do not mechanically replace text outside its historical context.

Classify every match:

- current-state stale -> repair;
- historical chronology -> preserve but, if needed, label historical;
- unrelated scope -> leave unchanged.

Then audit the final memory set for contradictory current claims about:

- Fix.5.7 status;
- Fix.5.8 status;
- E.8 overall status;
- canonical-five deferral;
- F.6.4;
- Method A;
- Method B;
- absolute SIMC;
- sole NEXT.

There must be one and only one ordinary live NEXT, in CURRENT.

---

## 13. Memory manifest and health checks

Discover the repository-selected local Python interpreter using the existing
TOOLS guidance.

After all versionable memory changes:

1. regenerate `docs/memory/manifest.json` using the repository's existing
   manifest procedure;
2. verify it;
3. run:

```text
<PYTHON> -B tools/check_memory_health.py --root .
```

This checkpoint explicitly owns consistency cleanup, so run the zero-warning
variant as well if the tool supports the existing documented
`--fail-on-warning` mode.

Report:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
warning classification: <none or each warning: blocking | nonblocking>
manifest check: PASS | FAIL
```

Hard failures block completion.

A warning caused by newly contradictory active status/NEXT is blocking.
A purely historical warning may be classified nonblocking only if the ordinary
health gate passes and the exact reason is recorded.

Run:

```bash
git diff --check
```

No executable tests are required for a memory-only closure task unless a
repository memory tool itself requires them.

---

## 14. Diff audit

Before stopping:

```bash
git status --short --untracked-files=all
git diff --stat
git diff --check
```

Confirm:

- no file outside the allowlist changed;
- no analysis/test/owner/collector/launcher/profile source changed;
- the new evidence file contains only accepted facts;
- CURRENT is concise and push-stable;
- roadmap and validation history do not compete with CURRENT;
- manifest is regenerated after the final memory contents.

Create one temporary complete review bundle at repository root:

```text
kaonlt_review.diff
```

It must contain:

- complete tracked diff for every changed tracked path;
- complete `git diff --no-index /dev/null ...` addition for both intended new
  files;
- no unrelated user-owned content;
- no staging merely for review.

---

## 15. Acceptance criteria

PASS only if all are true:

1. exact starting branch/HEAD/origin identity is verified;
2. no executable/scientific source changes;
3. exact Fix.5.8 accepted farm identity and ZIP hash/size are recorded;
4. Fix.5.8 is narrowly `CLOSED / RUNTIME VALIDATED`;
5. Fix.5.7 remains narrowly closed;
6. E.8 remains `ACTIVE`;
7. canonical-five E.8.1 remains `DEFERRED`;
8. final E.8 and F.6.4 remain `BLOCKED`;
9. Method A remains detached/non-production;
10. Method B remains diagnostic-only/numerically excluded;
11. absolute-SIMC provenance blocker remains explicit;
12. consumed pre-farm/pre-push wording is removed from current-state records;
13. historical chronology is preserved where appropriate;
14. CURRENT is the sole ordinary live-state/NEXT authority;
15. validation history no longer advertises an obsolete current gate;
16. roadmap reflects the accepted narrow gate without broadening scope;
17. manifest regeneration/check passes;
18. ordinary memory health passes;
19. zero-warning health passes if supported by the existing documented mode;
20. `git diff --check` passes;
21. complete review bundle is prepared;
22. no commit, push, remote-ref update, or farm run occurs.

---

## 16. Hard stop

After memory/evidence edits, consistency audit, manifest regeneration, health
checks, and complete `kaonlt_review.diff` preparation:

**STOP.**

Do not commit.

Do not push.

Do not run the farm.

Return the exact changed paths, memory-health report, stale-reference audit
summary, and complete actual-diff review bundle to ChatGPT.
