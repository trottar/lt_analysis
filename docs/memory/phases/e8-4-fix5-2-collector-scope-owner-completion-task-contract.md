# KaonLT E.8.4 Fix.5.2 — collector-scope and farm-owner completion

## Objective

Continue the existing local E.8.4 Fix.5/Fix.5.1 candidate from exact committed `test` HEAD

```text
da38444e7aa60efd62d6638780776344daf40276
```

without discarding, restarting, or redesigning the already-implemented E.8 presentation work.

The immediate objective is to remove the single operational blocker preventing the new tracked
`Q4p4W2p74 / Left / lowe` farm plot gate from running and returning a complete validation ZIP.

The scientific deliverable that must remain intact is the already-implemented E.8.4 shareable
Method-A evaluation:

1. final Method-A clean-kaon MM versus authoritative SIMC for every canonical `(t,phi)` child;
2. baseline MM0 versus Method-A MMA versus authoritative SIMC for every canonical `(t,phi)` child;
3. `Y0` and `YA` versus phi for each t parent;
4. `DeltaY = YA - Y0` versus phi for each t parent;
5. `DeltaY / Y0` versus phi for each t parent where defined;
6. explicit parent-t normalization-closure sanity page;
7. all existing E.8/E.8.2/E.8.3/E.8.4 pages preserved.

All 27 canonical `(t,phi)` cells must remain represented.

No production physics change is allowed.

---

## Starting authority and evidence

Required branch:

```text
test
```

Required committed HEAD:

```text
da38444e7aa60efd62d6638780776344daf40276
```

Commit subject:

```text
F6.3: adopt detached current-baseline candidate lineage
```

Newest cumulative local review candidate already supplied before this task:

```text
kaonlt_review(20261001-214401).diff
```

That candidate already contains the substantive Fix.5/Fix.5.1 implementation. Preserve it.

Previously accepted runtime evidence:

```text
KaonLT_F6_3_Left_lowe_runtime_evidence_Q4p4W2p74_20261001-222142.zip
SHA-256:
200fda66fe410274df1c9a8252b9e87114d8fd10a520fb7d71691ec6b3772874
```

Farm/source HEAD:

```text
da38444e7aa60efd62d6638780776344daf40276
```

Current-baseline candidate identities:

```text
F.3 SHA-256:
eeb480d9b7ddffab1e97c7f9f27c3099d0780c2bbb4d42a0ee68c077a5752b3d

F.4 SHA-256:
1d545924eba89c7f9ffa28028e307aca9b434a89beec06863cf2893887b6b902
```

Established status from that runtime evidence:

```text
F.4.Refresh.2:
CLOSED / RUNTIME VALIDATED

F.6.3 Q4p4W2p74 / Left / lowe:
CLOSED / RUNTIME VALIDATED

Existing E.8.4 Left/lowe Method-A impact/runtime gate:
CLOSED / RUNTIME VALIDATED

E.8.4 Fix.5 shareable SIMC/yield presentation:
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING after this local source task passes
```

Do not close canonical-five E.8.
Do not close F.6.4.
Do not promote Method A.

Representative already-observed Left/lowe yield effects are evidence only, not source constants:

```text
t1 phi [-180,-140): DeltaY/Y0 ≈ -25.8%
t1 phi [ 140, 180): DeltaY/Y0 ≈ -45.1%

t2 phi [-180,-140): ≈ -4.55%
t2 phi [ 100, 140): ≈ +8.65%

t3 phi [-180,-140): ≈ -1.83%
t3 phi [ 100, 140): ≈ +2.84%
```

Parent-t preservation closed to floating-point precision.

---

## Mandatory startup audit

Before editing:

1. establish the actual local branch, committed HEAD, and complete worktree;
2. do not reset, stash, clean, checkout-over, discard, or reconstruct the existing Fix.5/Fix.5.1 candidate;
3. read repository memory in the standard order.

Read in order:

```text
AGENTS.md
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/handoffs/CURRENT_HANDOFF.md
docs/memory/USER.md
```

Then read:

```text
docs/memory/CODEX.md
docs/memory/COMMUNICATION.md
docs/memory/MAINTENANCE.md
docs/memory/LEARNINGS.md
docs/memory/TOOLS.md
docs/memory/roadmap/STATUS.md
docs/memory/decisions/e8-full-analysis-procedure-roadmap.md
docs/memory/decisions/farm-validation-bundle-procedure.md
docs/memory/phases/e8-4-fix5-shareable-method-a-impact-pages-task-contract.md
docs/memory/phases/e8-4-fix5-1-simc-runtime-handoff-unblock-task-contract.md
docs/memory/phases/e8-4-fix5-2-collector-scope-owner-completion-task-contract.md
```

Inspect:

```bash
git branch --show-current
git rev-parse HEAD
git status --short --untracked-files=all
git log --oneline -1
git -c core.safecrlf=false diff --check
```

Required branch is `test`.

Required committed HEAD is exactly:

```text
da38444e7aa60efd62d6638780776344daf40276
```

The pre-existing Fix.5/Fix.5.1 worktree changes are expected. Unrelated source changes are a blocker.

---

## Scientific ownership and frozen physics

Preserve exactly:

- random subtraction;
- dummy subtraction;
- slow-proton subtraction;
- pion-component models/fits/windows/amplitudes;
- active `no_empirical_residual` profile;
- baseline pion weight `w0`;
- Method-A correction factors `C`;
- F.4 parent-preserving mathematics;
- F.5 propagation mathematics;
- F.6.3 `w0 -> w0*C` mathematics;
- live-cache parity;
- canonical t/phi binning;
- existing yield extraction;
- Lambda integration window;
- data normalization;
- SIMC event weighting and normalization;
- efficiencies;
- acceptance;
- L/T separation;
- cross-section formulas;
- Method-B diagnostic-only status.

Method A remains detached/non-production.

Method B remains diagnostic/cross-check only and must never enter the numerical Method-A path.

Do not modify physics, normalization, cuts, fit windows, templates, priors, binning,
subtraction formulas, component definitions, yields, efficiencies, acceptance, L/T separation,
or cross-section logic to improve agreement or presentation.

E.8 remains a consumer/presentation layer.

---

## Existing Fix.5/Fix.5.1 candidate that must be preserved

The existing local candidate already established the authoritative SIMC path:

```text
find_yield_data
-> find_yield_simc
-> existing normfac_simc
-> hist["_xsect_support_simc"]["mm"][t][phi]
-> E.8 finalizer
-> validate/clone existing SIMC child MM
-> render E.8 pages
```

No second SIMC producer is allowed.

Do not reopen SIMC files in E.8.
Do not refill SIMC.
Do not recalculate normalization.
Do not shape-normalize data or SIMC.
Do not add a new comparison normalization.

The existing candidate already contains:

- the narrow `src/main.py` caller-order repair;
- authoritative E.8 SIMC handoff;
- `full_background.e8_4.method_a_vs_simc.t1/t2/t3`;
- `full_background.e8_4.baseline_method_a_simc.t1/t2/t3`;
- `full_background.e8_4.yield_summary.t1/t2/t3`;
- `full_background.e8_4.parent_closure`;
- all 27 canonical `(t,phi)` cells;
- preservation of existing E.8 pages;
- focused deterministic tests;
- the tracked Left/lowe farm owner/profile candidate.

These files may already differ from committed HEAD as part of the existing candidate:

```text
src/main.py
src/cuts/full_background_subtraction_plots.py
testing/test_e8_4_fix5_main_order.py
testing/test_e8_4_production_impact_audit.py
```

Preserve their current candidate behavior. Do not redesign them. Edit them only if a deterministic
test reveals a concrete defect required to satisfy this contract, and report that explicitly.

---

## Concrete blocker

The generic collector:

```text
testing/collect_pion_hgcer_validation_bundle.py
```

requires every validation profile to declare exactly the canonical five settings:

```text
Left / lowe
Left / highe
Center / lowe
Center / highe
Right / highe
```

The current Fix.5 profile declares only `Left / lowe`, causing:

```text
validation_bundle_profile_invalid
```

The existing generic collector already supports selecting one authorized setting with explicit
`phi` and `epsilon`.

Therefore the correct repair is:

```text
canonical-five profile declaration
+
explicit collection selection phi="Left", epsilon="lowe"
```

Do not modify the generic collector.

---

## Required Fix.5.2 source changes

### 1. Validation profile

Edit:

```text
testing/pion_hgcer_validation_bundle_profile_e8_4_fix5.json
```

Change only the settings declaration to:

```json
[
  {"phi": "Left", "epsilon": "lowe"},
  {"phi": "Left", "epsilon": "highe"},
  {"phi": "Center", "epsilon": "lowe"},
  {"phi": "Center", "epsilon": "highe"},
  {"phi": "Right", "epsilon": "highe"}
]
```

Keep the existing artifact inventory unchanged.

Global artifacts remain:

```text
candidate F.3
candidate F.4
Fix.5 run summary
```

Setting-scoped artifacts remain:

```text
procedure PDF
page manifest
full-analysis JSON
no_empirical_residual correction-ledger JSON
no_empirical_residual correction-ledger CSV
```

Do not add the large yield-data PDF.

Keep the tracked source-base placeholder design. The farm owner binds only its in-memory effective
profile to the exact reviewed/pushed source SHA at runtime.

### 2. Tracked farm owner

Edit:

```text
testing/run_e8_4_fix5_left_lowe_plot_gate.py
```

`resolved_profile(...)` must:

- require exact canonical-five profile declaration;
- require `("Left", "lowe")` to be authorized;
- require `allowed_committed_files == []`;
- require the tracked base placeholder to remain exact `BASE_HEAD`;
- change only the in-memory effective profile's
  `source_identity.required_analysis_commit` to the supplied exact pushed `source_commit`;
- never rewrite the tracked profile on disk.

The owner must call the existing generic collector with explicit:

```python
phi="Left",
epsilon="lowe",
```

The resulting bundle manifest must contain exactly:

```json
"requested_settings": [
  {"phi": "Left", "epsilon": "lowe"}
]
```

and exactly one setting entry.

### 3. Deterministic success/failure stdout

Success stdout must contain exactly one ZIP path and it must be POSIX formatted:

```python
print(result.as_posix())
```

Failure stdout must contain no ZIP path.

No extra success chatter may be printed to stdout.

### 4. Farm dirt and source checking

The ordinary farm checkout may contain known farm-generated model outputs:

```text
M src/models/xmodel_kaon_pl.f
?? src/kaon/functions/Q4p4W2p74.model
```

Preserve and report them. Do not clean, reset, stash, delete, or rewrite them.

Before the expensive analysis run, fail closed on any unexpected gate-relevant modified/untracked source.

Do not impose a contradictory repository-wide pre-run `git diff --check` on the ordinary farm
checkout if it would fail solely because of known allowed farm-generated outputs.

The generic collector itself performs repository-wide source checks. Therefore, after the analysis
and artifact verification, collection must run from a temporary clean detached Git worktree at the
exact reviewed/pushed `source_commit` while reading scientific artifacts from the ordinary canonical
farm output directory.

Required structure:

```text
ordinary farm checkout
  -> exact source/provenance preflight
  -> expensive Left/lowe analysis
  -> artifact/page verification
  -> write run summary
  -> create temporary clean detached worktree at exact source_commit
  -> run existing generic collector with:
       repo_root = temporary clean worktree
       outdir = ordinary canonical artifact directory
       effective canonical-five profile
       phi = Left
       epsilon = lowe
  -> verify fresh ZIP
  -> remove/prune only that temporary worktree
```

Never run analysis in the temporary worktree.

Never clean/reset/stash the ordinary farm checkout.

---

## Allowed files for Fix.5.2 edits

Operational/profile/test files:

```text
testing/pion_hgcer_validation_bundle_profile_e8_4_fix5.json
testing/run_e8_4_fix5_left_lowe_plot_gate.py
testing/test_pion_hgcer_validation_bundle_profile_e8_4_fix5.py
testing/test_run_e8_4_fix5_left_lowe_plot_gate.py
```

Existing scientific/presentation candidate files may remain changed from HEAD but should not receive
new edits unless a deterministic defect requires a narrowly scoped repair:

```text
src/main.py
src/cuts/full_background_subtraction_plots.py
testing/test_e8_4_fix5_main_order.py
testing/test_e8_4_production_impact_audit.py
```

Memory/status files authorized for the already-warranted alignment:

```text
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/USER.md
docs/memory/LEARNINGS.md
docs/memory/TOOLS.md
docs/memory/roadmap/STATUS.md
docs/memory/decisions/e8-full-analysis-procedure-roadmap.md
docs/memory/manifest.json
docs/memory/phases/f6-3-current-baseline-candidate-lineage-adoption.md
docs/memory/evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md
docs/memory/phases/e8-4-fix5-shareable-method-a-impact-pages.md
docs/memory/phases/e8-4-fix5-2-collector-scope-owner-completion-task-contract.md
```

Retain the earlier Fix.5 and Fix.5.1 contracts as history.

---

## Frozen files / forbidden changes

Do not modify:

```text
testing/collect_pion_hgcer_validation_bundle.py
testing/package_pion_hgcer_validation_bundle.tcsh

src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
src/cuts/pion_hgcer_method_a_parent_preserving_correction.py
src/cuts/pion_hgcer_method_a_tphi_propagation.py
src/cuts/pion_component_subtraction.py
src/cuts/pion_component_fits.py
src/cuts/particle_subtraction.py
src/cuts/calculate_yield.py
binning/calculate_yield.py
run_Prod_Analysis.sh
```

Do not alter any Method-B source.

Do not alter F.1-F.6.2 accepted authority artifacts/identities.

Do not create another SIMC producer, normalization, or scaling path.

Do not create ad-hoc plots as substitutes for the tracked E.8 procedure-PDF pages.

Do not commit, push, or run the farm.

---

## Owner/profile deterministic tests

Update the two focused owner/profile test modules to prove at minimum:

1. profile loader accepts the canonical-five profile;
2. all five canonical settings occur in exact required order;
3. effective profile changes only the required analysis source commit;
4. tracked profile on disk remains unchanged;
5. owner explicitly selects `Left / lowe` only during collection;
6. synthetic bundle has `requested_settings == [Left/lowe]`;
7. synthetic bundle contains exactly one setting;
8. candidate F.3/F.4 and all required Left/lowe artifacts are present;
9. wrong branch/head/origin source fails before analysis;
10. unexpected gate-relevant source dirt fails before analysis;
11. known farm model-output dirt is recorded and never cleaned;
12. allowed farm outputs do not spuriously block the gate;
13. failed analysis stops before collection;
14. stale required artifact fails before collection;
15. collection failure fails closed;
16. missing/duplicate new E.8 page IDs fail;
17. nonempty renderer failures fail;
18. missing old critical E.8 pages fail;
19. pre-existing output ZIP is refused;
20. success prints exactly one POSIX ZIP path;
21. failure prints no ZIP path;
22. temporary detached worktree is used only for collector/source checks;
23. analysis is never invoked from the temporary worktree;
24. only the created temporary worktree is removed/pruned afterward.

Use the real generic collector in the synthetic success-path test where practical.

---

## Required new E.8 page inventory to verify

The owner must require each new page ID exactly once:

```text
full_background.e8_4.method_a_vs_simc.t1
full_background.e8_4.method_a_vs_simc.t2
full_background.e8_4.method_a_vs_simc.t3

full_background.e8_4.baseline_method_a_simc.t1
full_background.e8_4.baseline_method_a_simc.t2
full_background.e8_4.baseline_method_a_simc.t3

full_background.e8_4.yield_summary.t1
full_background.e8_4.yield_summary.t2
full_background.e8_4.yield_summary.t3

full_background.e8_4.parent_closure
```

Each t-scoped MM page must represent all nine canonical phi children.

Preserve existing critical E.8.4 page IDs, including:

```text
full_background.e8_4.authority
full_background.e8_4.pion_consequence.t1/t2/t3
full_background.e8_4.final_mm.t1/t2/t3
full_background.e8_4.signed_difference.t1/t2/t3
full_background.e8_4.yield_impact.t1/t2/t3
full_background.e8_4.setting_summary
```

Require:

```text
renderer_failures == []
```

---

## Farm-owner acceptance path

After user commit/push and independent pushed-state review, the tracked owner must own the complete
farm operation:

```text
preflight exact pushed HEAD / origin/test / gate-relevant worktree
-> verify exact candidate F.3/F.4 SHA-256
-> run existing ./run_Prod_Analysis.sh -d 4p4 2p74
-> confirm Left/lowe debug completion
-> verify fresh procedure PDF + page manifest
-> require renderer_failures == []
-> require every new E.8.4 page ID exactly once
-> require all existing critical E.8 pages
-> require all nine phi children on each new t-scoped page
-> create run-summary provenance JSON
-> create temporary clean detached worktree at exact pushed commit
-> invoke existing generic collector from that clean worktree with phi=Left, epsilon=lowe
-> verify ZIP integrity/source identity/artifact hashes
-> remove/prune only the created temporary worktree
-> print exactly one fresh POSIX ZIP path
```

No second manual packaging command is permitted.

---

## Memory alignment

Complete the already-authorized milestone alignment in this same task.

Record narrowly:

```text
F.4.Refresh.2:
CLOSED / RUNTIME VALIDATED

F.6.3 Q4p4W2p74 / Left / lowe:
CLOSED / RUNTIME VALIDATED

Existing E.8.4 Left/lowe Method-A impact/runtime gate:
CLOSED / RUNTIME VALIDATED

E.8.4 Fix.5 shareable SIMC/yield presentation:
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

Do not close canonical-five E.8.
Do not close F.6.4.
Do not promote Method A.

Preserve the post-run worktree provenance caveat:

```text
M src/models/xmodel_kaon_pl.f
?? src/kaon/functions/Q4p4W2p74.model
```

without reclassifying it as a demonstrated scientific failure.

`TOOLS.md` ordinary memory-health policy must be:

```text
<PYTHON> -B tools/check_memory_health.py --root .
```

`--fail-on-warning` is reserved for an explicit milestone/zero-warning audit or a warning already
classified as materially blocking.

Leave exactly one ordinary CURRENT/NEXT:

```text
NEXT — after user commit/push and pushed-state synchronization, run the tracked Q4p4W2p74 / Left / lowe E.8.4.Fix.5 farm owner to regenerate the full-analysis procedure PDF with the authoritative per-(t,phi) SIMC comparisons and shareable Method-A yield-impact pages, verify them, and return its single fresh validation ZIP.
```

Commit/push is a synchronization condition, not the scientific NEXT.

Regenerate and check:

```text
docs/memory/manifest.json
```

---

## Local deterministic validation

Use the repository-authoritative Python identified by the repository.

Run at minimum:

```bash
<PYTHON> -B -m py_compile \
  src/cuts/full_background_subtraction_plots.py \
  testing/run_e8_4_fix5_left_lowe_plot_gate.py \
  testing/test_e8_4_fix5_main_order.py \
  testing/test_e8_4_production_impact_audit.py \
  testing/test_pion_hgcer_validation_bundle_profile_e8_4_fix5.py \
  testing/test_run_e8_4_fix5_left_lowe_plot_gate.py

<PYTHON> -B -m unittest testing.test_e8_4_fix5_main_order -v
<PYTHON> -B -m unittest testing.test_e8_4_production_impact_audit -v
<PYTHON> -B -m unittest testing.test_full_background_subtraction_plots -v
<PYTHON> -B -m unittest testing.test_f6_3_parallel_full_procedure_method_a -v
<PYTHON> -B -m unittest testing.test_run_prod_analysis_debug_left_low -v
<PYTHON> -B -m unittest testing.test_pion_hgcer_validation_bundle_profile_e8_4_fix5 -v
<PYTHON> -B -m unittest testing.test_run_e8_4_fix5_left_lowe_plot_gate -v
<PYTHON> -B -m unittest testing.test_collect_pion_hgcer_validation_bundle -v
```

Then:

```bash
<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/memory_bootstrap.py --root . --json
<PYTHON> -m unittest testing.test_memory_health
<PYTHON> -B tools/check_memory_health.py --root .
git -c core.safecrlf=false diff --check
```

Report exact test counts and skips.

Do not run `main.py` locally.
Do not claim ROOT/PyROOT/farm rendering validation.

---

## Diff audit

Before stopping:

```bash
git status --short --untracked-files=all
git -c core.safecrlf=false diff --stat
git -c core.safecrlf=false diff --check
git -c core.safecrlf=false diff --no-ext-diff
```

Explicitly verify:

- generic collector unchanged;
- SIMC producer unchanged;
- SIMC normalization unchanged;
- Method-A mathematics unchanged;
- Method-B source unchanged;
- random/proton/pion production subtraction unchanged;
- only the already-reviewed caller-order move remains in `src/main.py`;
- existing shareable E.8 page implementation remains intact;
- farm owner/profile now form a complete one-command Left/lowe gate;
- memory manifest is fresh;
- ordinary memory health passes.

Create one fresh complete repository-root review bundle:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

It must include:

- branch/HEAD/status;
- complete cumulative tracked diff from the committed base;
- complete `git diff --no-index /dev/null ...` for every intended new/untracked file;
- changed-path inventory;
- all local validation outputs;
- memory manifest/health output;
- CURRENT/MEMORY/CURRENT_HANDOFF byte counts.

Do not stage merely for review.

---

## Acceptance criteria

PASS only if:

1. branch is `test`;
2. committed HEAD remains exact `da38444e7aa60efd62d6638780776344daf40276`;
3. existing Fix.5/Fix.5.1 scientific/presentation candidate is preserved;
4. all 27 canonical `(t,phi)` cells remain represented in the new tracked E.8 pages;
5. authoritative per-cell SIMC path remains the existing `find_yield_simc -> normfac_simc -> hist["_xsect_support_simc"]["mm"]` path;
6. no SIMC normalization/shape normalization/reload/refill is added;
7. generic collector is unchanged;
8. profile validates with canonical five;
9. owner explicitly collects only Left/lowe;
10. output manifest contains exactly one requested setting: Left/lowe;
11. effective profile is pinned in memory to the exact supplied source commit;
12. tracked profile on disk remains unchanged at runtime;
13. success stdout is exactly one POSIX ZIP path;
14. failure stdout contains no ZIP path;
15. known farm-generated model outputs are preserved and recorded;
16. unexpected gate-relevant source dirt fails before the expensive run;
17. generic collector source checks run from a temporary clean detached worktree at exact source commit;
18. scientific analysis never runs in that temporary worktree;
19. only that temporary worktree is removed/pruned;
20. all E.8/SIMC/yield tests pass;
21. all owner/profile/collector tests pass;
22. memory/status alignment is accurate and narrow;
23. manifest is fresh;
24. ordinary memory health passes;
25. one fresh complete `kaonlt_review(...).diff` is returned;
26. Codex does not commit, push, or run the farm.

---

## Hard stop

After deterministic local completion and creation of the fresh review bundle:

**STOP.**

Do not:

- commit;
- push;
- run the farm;
- broaden to Left/highe or canonical five;
- modify production physics;
- promote Method A;
- close F.6.4;
- create another memory-only phase.

Return the fresh `kaonlt_review(...).diff` for independent ChatGPT actual-diff review.
