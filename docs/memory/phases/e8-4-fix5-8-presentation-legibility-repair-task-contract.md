# KaonLT — E.8.4 Fix.5.8 presentation legibility repair

## 1. Objective and exact starting source

Repair only the minor presentation defects found during independent review of the fresh
`Q4p4W2p74 / Left / lowe` Fix.5.7 farm bundle.

Required branch:

```text
test
```

Required starting remote/source identity:

```text
4da101ef2a0f633d3b98f95cc5dbe8f903b94a35
```

Commit subject at that identity:

```text
E8.4 Fix.5.7: repair page-manifest setting validation
```

This is a presentation-only repair. It must not alter any scientific calculation,
production subtraction/correction, Method-A numerical value, Method-B behavior, SIMC
normalization, yield, uncertainty, cut, template, prior, binning, efficiency,
acceptance, L/T separation, cross section, candidate artifact, page identity, or
bundle inventory.

No farm run is authorized by this contract.

---

## 2. Mandatory startup and re-anchor

Read in exact order:

1. repository-root `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read only task-relevant records/source:

- `docs/memory/MAINTENANCE.md`
- `docs/memory/CODEX.md`
- `docs/memory/phases/e8-4-fix5-7-page-manifest-setting-provenance-owner-repair.md`
- `docs/memory/investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md`
- `src/cuts/full_background_subtraction_plots.py`
- `testing/test_e8_4_fix5_5_visualization.py`
- directly relevant fake-ROOT helpers/tests only if needed

Before editing, establish and report:

```bash
git branch --show-current
git rev-parse HEAD
git rev-parse refs/remotes/origin/test
git status --short --untracked-files=all
```

Hard requirements:

- branch is exactly `test`;
- local HEAD is exactly `4da101ef2a0f633d3b98f95cc5dbe8f903b94a35`;
- local `origin/test` is exactly the same identity;
- preserve all unrelated worktree state;
- this task-contract file may be present as the expected untracked input;
- any unrelated tracked modification is a hard stop;
- do not reset, clean, stash, checkout over, commit, push, update remote refs, or run the farm.

`CURRENT.md` is stale with respect to the newly reviewed farm evidence. Perform the
opening memory checkpoint first: reconcile that fresh evidence and the now-completed
Fix.5.7 runtime gate before making the presentation edit. Do not create a separate
maintenance-only phase or commit.

---

## 3. Fresh farm evidence to record

The user supplied and ChatGPT independently inspected:

```text
KaonLT_E8_4_Fix5_Left_lowe_Q4p4W2p74_20261003-002936-622785.zip
```

ZIP SHA-256:

```text
d39b81544b483f332ae45cb0342aa7d95b2c1931d866cfcb04a52e46449fe38d
```

ZIP size:

```text
30851202 bytes
```

The reviewed bundle establishes, for this one `Q4p4W2p74 / Left / lowe` attempt:

- ZIP integrity passes;
- bundle `complete=true` and `errors=[]`;
- `git_head` and `required_analysis_commit` are both exactly
  `4da101ef2a0f633d3b98f95cc5dbe8f903b94a35`;
- requested scope is exactly `Left / lowe`;
- validation profile is `phase_e8_4_fix5_left_lowe_shareable_pages/v1`;
- analysis return code is zero;
- page count is 97;
- renderer failures are empty;
- the page manifest contains the full seven-field producer setting:
  `Q2`, `W`, `epsilon_filename_token`, `epsilon_setting`, `kinematic_token`,
  `particle_type`, and `phi_setting`;
- therefore the Fix.5.7 owner/checker setting-provenance repair succeeds on the real
  farm path through collection and ZIP verification;
- the unresolved absolute-SIMC provenance blocker remains unchanged;
- Method A remains detached/non-production and Method B remains numerically excluded.

Record Fix.5.7 as `CLOSED / RUNTIME VALIDATED` for its narrow owner/checker repair at
`Q4p4W2p74 / Left / lowe`. Do not generalize that closure to canonical-five E.8,
Fix.5.5 visual acceptance, Method-A promotion, or F.6.4.

The same fresh PDF fails the final rendered-page legibility gate only for the minor
presentation defects listed below. Those defects are runtime evidence and are the
sole reason for Fix.5.8.

---

## 4. Exact visual defects from the fresh 97-page PDF

The page-manifest ordering maps the affected pages as follows:

```text
page 67  full_background.e8_3.tphi                t1
page 69  full_background.e8_3.tphi                t2
page 71  full_background.e8_3.tphi                t3
page 73  full_background.e8_4.authority            setting
page 80  full_background.e8_4.yield_summary.t1     t1
page 87  full_background.e8_4.yield_summary.t2     t2
page 94  full_background.e8_4.yield_summary.t3     t3
page 96  full_background.e8_4.setting_summary      setting
```

Observed defects:

1. Pages 67/69/71 render the literal Unicode em dash in the E.8.3 title as mojibake
   (`â€”`).
2. Page 96 has the same problem in the E.8.4 setting-summary title.
3. Page 73 clips the long single-line `Non-production flags:` text at the right edge.
4. Pages 80/87/94 clip the long per-phi stored signed-support lines at the right edge.

These are presentation-legibility failures only. The bundle/page checker, numerical
payload, parent preservation, page inventory, and source provenance passed.

---

## 5. SOURCE VERIFIED presentation mechanism

Current source at the required starting identity shows:

### E.8.3 mojibake source

`_e8_3_render_tphi_page(...)` sends this literal ROOT text:

```python
"E.8.3 persisted canonical (t,phi) pion redistribution — {} t{}"
```

### E.8.4 setting-summary mojibake source

`_e8_4_render_setting_summary_page(...)` sends:

```python
"E.8.4 stored canonical t-by-phi impact summary — {}"
```

### E.8.4 authority clipping source

`_e8_4_render_authority_page(...)` concatenates all seven flags into one
`Non-production flags:` line and passes it to `_e8_text_page(...)`.
`_e8_text_page(...)` is intentionally non-wrapping.

### E.8.4 yield-summary clipping source

`_e8_4_stored_support_lines(...)` currently formats, for every physical phi child,
all `MM_0`, `MM_A`, and `delta_MM` signed/absolute support values into one long line.
`_e8_4_render_yield_summary_page(...)` passes those lines to one fixed lower
`TPaveText` box.

Repair these display strings/line breaks only. Do not change how any value is obtained.

---

## 6. Scientific ownership and frozen interfaces

This task owns presentation formatting only.

Preserve exactly:

- random/dummy subtraction;
- slow-proton subtraction;
- pion subtraction;
- `no_empirical_residual` profile;
- F.6.3 branch construction and parent preservation;
- Method-A values, fingerprints, candidates, and event weights;
- Method-B diagnostic-only/numerically-absent role;
- E.8.3 historical F.6.1 lineage identity;
- current F.6.3/E.8.4 lineage identity;
- SIMC support and unresolved absolute-unit blocker;
- all yields and statistical errors;
- all stored signed/absolute support values;
- all aggregate values;
- all page IDs, page scopes, page order, and page-manifest schemas;
- 97-page expected procedure-PDF inventory for this Left/lowe gate;
- owner/collector/profile/launcher behavior;
- eight-artifact bundle inventory.

Forbidden:

- recomputing any diagnostic or yield for presentation;
- integrating/re-summing histograms;
- changing payload schemas;
- shortening text by deleting provenance, flags, values, or diagnostic meaning;
- changing numerical precision solely to make text fit;
- changing scientific source outside the single presentation module;
- changing owner/checker thresholds to hide presentation failures.

---

## 7. Allowed substantive files

Expected implementation/test changes:

```text
src/cuts/full_background_subtraction_plots.py
testing/test_e8_4_fix5_5_visualization.py
```

Allowed memory/evidence files when warranted:

```text
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md
docs/memory/evidence/e8-4-fix5-7-left-lowe-runtime-and-fix5-visual-gate-2026-10-03.md
docs/memory/phases/e8-4-fix5-8-presentation-legibility-repair-task-contract.md
docs/memory/phases/e8-4-fix5-8-presentation-legibility-repair.md
docs/memory/manifest.json
```

Do not modify `MEMORY.md` unless a genuinely durable cross-phase rule changes; the
fresh runtime result and this repair normally belong in CURRENT/evidence/phase/investigation.

Frozen unless an unexpected contradiction forces a hard stop and a new contract:

```text
src/cuts/rand_sub.py
src/cuts/pion_hgcer_refinement_checkpoint.py
testing/run_e8_4_fix5_left_lowe_plot_gate.py
testing/test_run_e8_4_fix5_left_lowe_plot_gate.py
testing/collect_pion_hgcer_validation_bundle.py
testing/test_collect_pion_hgcer_validation_bundle.py
testing/pion_hgcer_validation_bundle_profile_e8_4_fix5.json
run_Prod_Analysis.sh
all other src/
```

---

## 8. Required repair A — ASCII-safe E.8.3/E.8.4 titles

For the two observed ROOT-facing title strings only:

- replace the Unicode em dash with an ASCII-safe separator such as ` - `;
- preserve all surrounding text and setting/t identity;
- do not alter E.8.3 lineage semantics or E.8.4 current-lineage semantics.

Required target renderers:

```text
_e8_3_render_tphi_page
_e8_4_render_setting_summary_page
```

Do not broadly rewrite unrelated E.7.2 presentation text merely because other Unicode em
dashes exist elsewhere in the module; they were not part of this failed visual gate.

---

## 9. Required repair B — E.8.4 authority-page bounded text

Keep `_e8_text_page(...)` non-wrapping unless a concrete deterministic need proves that
changing it is safer. Prefer a narrow authority-specific formatting repair.

The authority page must still display all seven exact flag/value pairs:

```text
baseline_production_mutated
production_promotion_performed
method_b_numerical_dependency
empirical_residual_used
event_correction_persisted
canonical_child_renormalization_performed
baseline_public_output_unchanged
```

Requirements:

- every flag appears exactly once with its stored Boolean value;
- `Non-production flags` remains visibly identified;
- split the content across deterministic bounded lines rather than relying on ROOT wrapping;
- target every authority-page text line to be <= 96 characters where practical;
- preserve the existing authority statements and their meaning;
- no payload mutation.

If another existing authority line exceeds the same bound, split it semantically rather
than lowering the whole page to unreadably small text.

---

## 10. Required repair C — E.8.4 yield-summary bounded stored-support text

`_e8_4_stored_support_lines(...)` remains a formatter of already stored producer-owned
audit values only.

For each physical phi child, retain exact formatted values for:

```text
MM_0 signed_integral
MM_0 absolute_support
MM_A signed_integral
MM_A absolute_support
delta_MM signed_integral
delta_MM absolute_support
```

Do not change the existing `.8g` precision unless a concrete rendering blocker forces a
hard stop for review.

Required display behavior:

- replace each overlong per-child line with deterministic bounded line(s);
- retain the phi child identity and edges;
- retain all six numerical support quantities;
- retain the stored aggregate `B_pi_0`, `B_pi_A`, and `delta` line;
- retain the current-lineage label and the F.4 non-comparability statement;
- no histogram integration/re-summing or child reaggregation;
- no numerical inference;
- no payload mutation.

The lower text region may be adjusted narrowly to fit the additional line count, but it
must not overlap the three upper yield panels or the page header. Keep a legible text
size; do not "fix" clipping by making the text microscopically small.

As a practical acceptance target, all support lines supplied to ROOT should be <= 96
characters. If preserving an indivisible numerical token makes one line slightly longer,
record the exact reason and use the smallest necessary bound; do not silently drop content.

---

## 11. Required deterministic regressions

Extend `testing/test_e8_4_fix5_5_visualization.py` rather than creating a new broad test
module unless a direct dependency makes that impossible.

At minimum add/extend tests proving:

### ASCII-safe titles

- E.8.3 t-phi title contains no `\u2014`;
- E.8.4 setting-summary title contains no `\u2014`;
- neither target title contains the observed mojibake sequence `â€”`;
- expected setting/t identity remains present.

### Authority page

- capture the exact lines passed to `_e8_text_page(...)`;
- every seven flag names appears exactly once with its stored value;
- `Non-production flags` label remains present;
- bounded line-length assertion covers the affected authority page;
- payload is unchanged before/after rendering/formatting.

### Yield-summary support block

Using deliberately distinctive stored support values:

- no `GetBinContent`/`GetBinError` integration is allowed, preserving the existing
  presentation-ownership test;
- every child still exposes `MM_0`, `MM_A`, and `delta` stored signed/absolute values;
- the aggregate line remains the exact stored aggregate values;
- current F.6.3 lineage and F.4 non-comparability text remain;
- all generated support lines meet the chosen deterministic character bound;
- payload remains unchanged.

### Page inventory and science isolation

Preserve existing regression proving:

- all 24 E.8.4 page records remain present;
- no page ID/scope changes;
- blocked SIMC placeholders remain six and no SIMC histogram is drawn when absolute
  comparison is unavailable;
- source histograms/payloads remain unchanged.

Add a source-level assertion that this Fix.5.8 change introduces no numerical producer,
histogram integration, or normalization call into the affected formatter/renderer region.

---

## 12. Memory checkpoint and status rules

Create:

```text
docs/memory/evidence/e8-4-fix5-7-left-lowe-runtime-and-fix5-visual-gate-2026-10-03.md
docs/memory/phases/e8-4-fix5-8-presentation-legibility-repair.md
```

The evidence record must distinguish:

- `RUNTIME VERIFIED`: Fix.5.7 owner/checker/page-setting path through final ZIP;
- `RUNTIME VERIFIED`: 97-page render and exact visual defects listed in this contract;
- `BLOCKED`: final Fix.5 visual acceptance due only to those presentation defects;
- `NOT VERIFIED`: repaired Fix.5.8 rendering until a later farm run.

Update `CURRENT.md` so that:

- Fix.5.7 is `CLOSED / RUNTIME VALIDATED` only for its narrow Left/lowe owner/checker
  repair;
- Fix.5.8 is the current work item;
- after local implementation/tests, Fix.5.8 becomes
  `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`;
- canonical-five E.8 remains `BLOCKED`/not expanded;
- Method A remains detached/non-production;
- Method B remains diagnostic-only/numerically absent;
- the SIMC absolute-unit blocker remains intact;
- stale wording that actual-diff/push/farm review of Fix.5.7 is still pending is removed.

`CURRENT.md` is currently above its 8 KiB soft warning and this is the scheduled next
memory checkpoint. Compact CURRENT narrowly while preserving authoritative active state,
accepted evidence, blockers, source identities, and exact NEXT. Target `< 8192` bytes.
Do not turn this into broad memory refactoring.

Push-stable NEXT after local completion must be:

```text
independent ChatGPT actual-diff review
-> user commit/push
-> ChatGPT pushed-state synchronization
-> one narrow Q4p4W2p74 / Left / lowe tracked-owner farm validation
-> fresh structural + targeted visual review of pages 67,69,71,73,80,87,94,96
```

No canonical-five run is authorized by this task.

Regenerate `docs/memory/manifest.json` after versionable memory changes.

---

## 13. Local validation

Use the repository-selected local Python interpreter.

Run syntax checks for changed Python files.

Run at minimum:

```text
python -B -m unittest testing.test_e8_4_fix5_5_visualization
python -B -m unittest testing.test_full_background_subtraction_plots
python -B -m unittest testing.test_run_e8_4_fix5_left_lowe_plot_gate
```

plus directly affected focused suites if source inspection shows another dependency.

Run:

```bash
git diff --check
```

After memory updates:

```text
<PYTHON> -B tools/check_memory_health.py --root .
```

Report exact pass/skip counts and warnings.

The ordinary health gate must have no hard failure. The existing CURRENT soft-size
warning should be eliminated in this checkpoint if possible without losing material state.
If it cannot be eliminated narrowly, stop and explain why rather than deleting required
authority/evidence.

Local/fake-ROOT tests do not establish actual ROOT font rendering, procedure-PDF
legibility, farm integration, or final visual acceptance.

---

## 14. Diff audit and review bundle

Before stopping:

```bash
git status --short --untracked-files=all
git diff --stat
git diff --check
```

Expected substantive implementation paths are only:

```text
src/cuts/full_background_subtraction_plots.py
testing/test_e8_4_fix5_5_visualization.py
```

plus the allowlisted memory/evidence/task-contract files above.

No owner, collector, launcher, profile, candidate, or unrelated scientific source may change.

Create one complete temporary repository-root review bundle named:

```text
kaonlt_review.diff
```

It must include:

- complete tracked diff for every changed tracked path;
- complete `git diff --no-index /dev/null ...` addition for every intended new file;
- no unrelated user-owned file;
- no staging merely for review.

Do not commit or push.

---

## 15. Acceptance criteria

PASS local implementation only if all are true:

1. exact starting branch/HEAD/origin identity is respected;
2. fresh Fix.5.7 runtime evidence is reconciled before substantive editing;
3. Fix.5.7 is recorded as narrowly `CLOSED / RUNTIME VALIDATED` without broader promotion;
4. only presentation source changes;
5. the two observed Unicode em dashes are replaced by ASCII-safe separators;
6. E.8.4 authority text retains every flag/value and no affected line clips by construction;
7. E.8.4 stored-support text retains every stored quantity and is bounded by construction;
8. no numerical precision/value, payload schema, histogram, page ID/scope/order, or bundle inventory changes;
9. existing baseline-versus-Method-A visibility styling is untouched;
10. six blocked absolute-SIMC pages remain explicitly unavailable with the same reason;
11. deterministic regressions prove immutability and presentation-only ownership;
12. required local suites, syntax, diff, manifest, and memory-health checks pass;
13. CURRENT/NEXT is push-stable and materially accurate;
14. no commit, push, remote-ref update, or farm run occurs.

---

## 16. Farm-validation boundary

This task does **not** establish repaired PDF legibility.

After independent actual-diff review, user-controlled commit/push, and pushed-state
synchronization, the next substantive gate is one narrow tracked-owner
`Q4p4W2p74 / Left / lowe` farm validation.

The later evidence review must first verify owner/ZIP/provenance/page-manifest success,
then visually inspect at minimum:

```text
67, 69, 71, 73, 80, 87, 94, 96
```

and confirm no collateral regression in the previously improved baseline-versus-Method-A
overlays.

Do not expand to canonical-five settings until that visual gate passes.

---

## 17. Hard stop

After implementation, deterministic local checks, memory/evidence updates, manifest
regeneration, and complete `kaonlt_review.diff` preparation:

**STOP.**

Do not commit.

Do not push.

Do not run the farm.

Return the complete actual-diff review material and KaonLT memory-health report to
ChatGPT.
