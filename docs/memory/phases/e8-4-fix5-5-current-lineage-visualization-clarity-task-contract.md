# KaonLT — E.8.4 Fix.5.5 current-lineage visualization clarity

## 1. Purpose

Implement the visualization-only stage that follows the pushed and
pushed-state-reviewed E.8.4 Fix.5.4 numerical/provenance audit.

This task exists because the fresh `Q4p4W2p74 / Left / lowe` procedure PDF made
the baseline blue curve difficult or impossible to distinguish beneath the
Method-A magenta curve in many overlays. It also exposed a narrative risk:
historical E.8.3/F.6.1 aggregate pages sit immediately before the current
F.6.3/E.8.4 branch and can be mistaken for the same lineage.

Fix.5.5 must improve visual separability and lineage clarity **without changing
any scientific value**.

It must also expose the already-validated Fix.5.4 signed-support/current-lineage
audit records in the existing E.8.4 presentation so that the next farm run can
show whether large fractional signed-yield shifts arise in cancellation-dominated
cells.

No new physics, normalization, fitting, correction, yield calculation, or SIMC
scale is allowed.

No farm run is authorized by this contract.

---

## 2. Exact starting state

Required branch:

```text
test
```

Required starting HEAD:

```text
761fbb6c03d2d7a10bb911cf84e9ba898496fab6
```

Commit subject:

```text
E8.4 Fix.5.4: audit current lineage identity and SIMC provenance
```

Required parent:

```text
71a912bf14c59eee6d52956462aa0aa39c7f7c84
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
- preserve the current worktree;
- `docs/memory/phases/workflow-continuity-hardening-task-contract.md` is
  unrelated user-owned untracked work and must remain untouched/untracked;
- `kaonlt_review.diff` is temporary review material and must not be committed;
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

- `docs/memory/phases/e8-4-fix5-4-current-lineage-identity-audit.md`
- `docs/memory/phases/e8-4-fix5-4-current-lineage-identity-audit-and-fix-task-contract.md`
- `docs/memory/phases/e8-4-fix5-4-fix1-simc-availability-and-manifest-scope-task-contract.md`
- `docs/memory/investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/LEARNINGS.md`
- `docs/memory/CODEX.md`

Relevant source/tests:

- `src/cuts/full_background_subtraction_plots.py`
- `testing/test_full_background_subtraction_plots.py`
- `testing/test_e8_4_production_impact_audit.py`
- `testing/test_e8_4_fix5_4_identity_audit.py`

Read only additional source directly required to verify that no producer/runtime
ownership is being changed.

---

## 4. Established source/evidence facts

Treat these as fixed inputs.

### 4.1 Fix.5.4 numerical/provenance implementation

At the starting HEAD:

- the private current-F.6.3 branch owns and validates the exact current
  `B_pi_0`, `B_pi_A`, `MM_0`, `MM_A`, `Y0`, and `YA`;
- producer-side identity checks validate:
  - common MM geometry;
  - exact final-baseline histogram fingerprint ownership;
  - binwise
    `MM_A - MM_0 = -(B_pi_A - B_pi_0)`;
  - `Integral(MM_0) = Y0`;
  - `Integral(MM_A) = YA`;
  - delta-integral closure;
  - signed positive/negative/absolute support;
  - current-F.6.3 child aggregation;
- these audits do not alter production objects.

Do not recompute those scientific audits in the renderer.

### 4.2 Historical/current lineage distinction

Historical E.8.3 consumes accepted persisted F.6.1/F.5/F.4 lineage.

Current E.8.4 consumes current-baseline F.6.3 candidate lineage.

They are both valid in their own scopes but are not the same lineage.

Historical E.8.3 must remain historical context, not the current-F.6.3 aggregate
explanation.

### 4.3 SIMC blocker

Fix.5.4 source tracing established that absolute SIMC luminosity/effective-charge
units are not source-proven.

Therefore:

```text
simc_absolute_comparison_available = False
```

for the actual current source until a separate authority establishes otherwise.

Only the absolute-SIMC comparison page families are provenance-blocked.
The rest of the current-F.6.3 E.8.4 payload remains available.

Fix.5.5 must preserve this exactly.

### 4.4 Confirmed visual problem

In the prior farm PDF, the baseline blue and Method-A magenta curves frequently
overlap so closely that the blue curve is hidden beneath magenta.

That is a presentation defect. It is not a numerical defect.

---

## 5. Scientific ownership and frozen interfaces

### 5.1 Fix.5.5 owns only

- E.8.3 lineage labeling/visual context;
- E.8.4 plot draw order;
- E.8.4 comparison line/marker styles;
- legends/notes needed to identify the curves independently;
- presentation of **already-stored** Fix.5.4 signed-support/current-lineage
  aggregate audit values;
- deterministic presentation tests;
- warranted memory updates.

### 5.2 Frozen

Do not change:

- `src/binning/calculate_yield.py`;
- `src/main.py`;
- `run_Prod_Analysis.sh`;
- Method-A factors, factor reconstruction, or event multipliers;
- F.3/F.4/F.5/F.6 artifact identities;
- `B_pi_0`, `B_pi_A`, `MM_0`, `MM_A`, `Y0`, `YA`;
- random/dummy subtraction;
- slow-proton subtraction;
- pion fits/templates/weights;
- Method B;
- cuts;
- MM/Lambda windows;
- binning;
- priors;
- efficiencies;
- acceptance;
- yield arithmetic;
- SIMC normalization/weighting/model;
- cross sections;
- L/T separation;
- production objects;
- owner/profile/collector behavior.

No independent child normalization.

No data or SIMC display scale.

No clipping/smoothing/interpolation intended to make curves look different.

No histogram value may be changed merely to improve visibility.

---

## 6. Required visual comparison style

Create one E.8.4-specific comparison-style helper or equivalent narrowly scoped
logic. Do **not** change the generic repository-wide histogram style helper.

Preserve the established semantic colors:

```text
baseline/current w0: blue
Method A:           magenta
SIMC:               black
common input:       black
```

Add redundant non-color encoding so the curves remain distinguishable in print
and under near-overlap.

Required roles:

### Baseline

Use:

```text
blue
solid line
thin-to-medium line width
open marker or another clearly visible non-filled marker
```

Draw the baseline **after** Method A in overlapping baseline/Method-A panels so
blue cannot be completely hidden beneath magenta.

### Method A

Use:

```text
magenta
dashed line
line width visibly wider than the baseline line
```

Draw Method A before the baseline.

The combination of wider dashed magenta beneath a narrower solid/marked blue
baseline should leave both identities visible when the numerical curves are
nearly identical.

### SIMC

When and only when absolute comparison is source-authorized:

```text
black
visually distinct solid/step line
```

SIMC must not obscure both data curves.

When absolute SIMC comparison is unavailable, retain the existing explicit
unavailable page and do not draw a SIMC histogram.

### Common pion input

Retain a neutral black style distinct from both `B_pi_0` and `B_pi_A`.

---

## 7. Required page-specific behavior

### 7.1 `full_background.e8_4.pion_consequence`

Current page content is:

```text
common pion input
B_pi_0
B_pi_A
```

Required draw order:

```text
common input
Method A B_pi_A
baseline B_pi_0 last
```

Use the E.8.4 comparison styles above.

Add a compact legend or equivalent unambiguous style note that identifies:

- common input;
- baseline `B_pi_0`;
- Method A `B_pi_A`.

Preserve:

- exact histogram contents;
- exact common y range;
- exact child inventory;
- exact page ID;
- exact canonical t/phi geometry.

### 7.2 `full_background.e8_4.final_mm`

Current page compares:

```text
MM_0
MM_A
```

Required draw order:

```text
Method A MM_A first
baseline MM_0 last
```

Use:

```text
MM_0 = blue solid + visible open markers
MM_A = magenta dashed + wider line
```

Add a compact per-pad legend or an equally unambiguous presentation element.

Retain:

- Lambda-window markers;
- stored `Y0` and `YA` text;
- common y range;
- page ID;
- child inventory;
- histogram values/errors.

Do not change the y-axis merely to magnify differences.

### 7.3 Absolute-SIMC pages

For:

```text
full_background.e8_4.method_a_vs_simc.t1/t2/t3
full_background.e8_4.baseline_method_a_simc.t1/t2/t3
```

When:

```text
simc_absolute_comparison_available == False
```

retain the Fix.5.4 explicit unavailable pages unchanged in scientific meaning.

Do not draw any black SIMC spectrum.

When a deterministic synthetic test provides a source-authorized comparable
SIMC case, the triple overlay order must be:

```text
SIMC
Method A
baseline last
```

and all three must have independently recognizable style.

No display renormalization is permitted.

### 7.4 Signed-difference page

Do not change the scientific construction:

```text
Delta B_pi = B_pi_A - B_pi_0
Delta K    = MM_A - MM_0
```

Only presentation labels may be improved if required for consistency.

No new difference arithmetic is needed.

### 7.5 Yield-summary page: expose Fix.5.4 signed-support audit

Keep the same page ID:

```text
full_background.e8_4.yield_summary.tN
```

Do not add a new page ID solely for this diagnostic.

Re-layout the existing yield-summary page if needed so that one panel or text
region presents the **already-stored** current-F.6.3 Fix.5.4 audit values for
that t parent.

For every canonical child, expose enough stored information to distinguish a
small signed yield from a large positive/negative cancellation.

Use the stored audit's:

```text
MM_0:
  positive_support
  negative_support
  signed_integral
  absolute_support

MM_A:
  positive_support
  negative_support
  signed_integral
  absolute_support

delta_MM:
  positive_support
  negative_support
  signed_integral
  absolute_support
```

At minimum the rendered page must show, per child:

```text
MM_0 signed integral and absolute support
MM_A signed integral and absolute support
delta signed integral and absolute support
```

Positive and negative supports should also be visible if they fit legibly.

The renderer must consume the stored Fix.5.4 audit values. It must **not**
re-integrate the ROOT histograms for this page.

Also show the corresponding current-F.6.3 aggregate values for that t parent:

```text
aggregate_B_pi_0_integral
aggregate_B_pi_A_integral
aggregate_delta_integral
```

with an explicit label:

```text
current F.6.3 candidate lineage; Lambda/allcut MM-template aggregate
```

and a reminder that it is not the broader F.4 application-population parent
sum when the stored audit marks those quantities non-comparable.

This is diagnostic presentation only.

### 7.6 Setting summary

Keep the existing setting-summary page ID.

Add a concise source-owned statement that:

- all displayed current-branch numerical identities come from the stored
  Fix.5.4 current-lineage audit;
- historical E.8.3 is separate;
- signed-support values are diagnostics, not uncertainties or corrections;
- absolute SIMC amplitude interpretation remains blocked unless explicitly
  source-authorized.

Do not turn this page into a new scientific calculation.

---

## 8. E.8.3 historical-lineage visual labeling

The existing E.8.3 pages are valid historical presentation but can be mistaken
for the current F.6.3 branch.

Improve labels only.

### Authority page

The title or first visible line must explicitly state:

```text
historical accepted F.6.1 lineage
```

and:

```text
not the current F.6.3/E.8.4 candidate lineage
```

### E.8.3 MM pages

The visible page title or note must state that the plotted aggregate arrays are:

```text
historical accepted F.6.1 aggregate context
```

and are not the current E.8.4 branch.

Do not alter the persisted arrays.

### E.8.3 t/phi pages

Add the same lineage distinction to the visible heading.

### E.8.3 cross-reference

Clarify that E.8.3 remains historical accepted-lineage context and must not be
used as the aggregate explanation for current F.6.3/E.8.4.

Do not rename or replace page IDs.

Do not change E.8.3 status or authority.

---

## 9. Page and manifest stability

Preserve all existing procedure page IDs.

Do not add or remove E.8.4 page families in this task.

Preserve the existing ordering unless a minimal local order change is required
inside one page's draw sequence.

In particular, retain:

```text
full_background.e8_4.authority
full_background.e8_4.pion_consequence
full_background.e8_4.final_mm
full_background.e8_4.signed_difference
full_background.e8_4.yield_impact
full_background.e8_4.method_a_vs_simc.tN
full_background.e8_4.baseline_method_a_simc.tN
full_background.e8_4.yield_summary.tN
full_background.e8_4.parent_closure
full_background.e8_4.setting_summary
```

SIMC page manifest records must retain their Fix.5.4
`simc_absolute_comparison_available` / literal reason semantics.

No page count should change solely because of Fix.5.5.

---

## 10. Allowed substantive files

Expected source:

```text
src/cuts/full_background_subtraction_plots.py
```

Expected tests:

```text
testing/test_full_background_subtraction_plots.py
testing/test_e8_4_production_impact_audit.py
```

The existing:

```text
testing/test_e8_4_fix5_4_identity_audit.py
```

may be modified only if required to assert that presentation consumes the
stored audit without mutation/recomputation.

One new focused visualization test is allowed if that is cleaner:

```text
testing/test_e8_4_fix5_5_visualization.py
```

Do not modify scientific producers.

---

## 11. Required deterministic tests

### 11.1 Draw-order/style tests

Using the repository's fake ROOT/test infrastructure, prove for `final_mm`:

```text
Method A is drawn before baseline
baseline is drawn last
baseline style != Method-A style beyond color alone
```

Specifically verify line-style/width/marker distinctions.

Prove the same baseline-last principle for `pion_consequence`.

For a synthetic comparable-SIMC payload, prove:

```text
SIMC -> Method A -> baseline
```

draw order and distinct styles.

For source-unproven SIMC, prove no SIMC histogram draw occurs and the explicit
unavailable page/status remains.

### 11.2 No numerical mutation

Snapshot before and after rendering:

- histogram contents/errors;
- `Y0`, `YA`, `delta_y`;
- current-lineage identity audit;
- aggregate audit;
- SIMC normalization audit.

Require exact equality.

### 11.3 Signed-support presentation

Provide synthetic stored Fix.5.4 signed-support audit values with strong
positive/negative cancellation.

Prove the yield-summary renderer displays values from the stored audit.

Do not permit the renderer to obtain those displayed values by re-integrating
the histograms.

A useful negative test is to make the stored audit support values intentionally
distinct from what a simple histogram re-sum would produce and confirm the
display follows the stored audit record. This is a presentation ownership test,
not an accepted scientific mismatch.

### 11.4 Lineage labeling

Require visible E.8.3 text to contain both concepts:

```text
historical F.6.1
not current F.6.3/E.8.4
```

or an equivalent unambiguous formulation.

Require E.8.4 current-lineage text to identify its aggregate as current F.6.3.

### 11.5 Page/manifest regression

Prove:

- no E.8.4 page ID is removed;
- no new page ID is introduced by Fix.5.5;
- SIMC blocked-page records remain available=false with the literal reason;
- existing non-SIMC E.8.4 page records remain available;
- renderer failure behavior remains fail-closed.

Run the existing Fix.5.4 identity tests unchanged as regressions.

---

## 12. Forbidden shortcuts

Do not:

- multiply one curve by an arbitrary factor;
- normalize curves to common area/maximum;
- change axes to hide/show agreement;
- clip negative bins;
- smooth/interpolate;
- alter binning;
- alter marker positions;
- recompute Method-A factors;
- recompute yields;
- derive a new SIMC scale;
- use historical E.8.3 arrays in current E.8.4 rendering;
- sum current children again in the renderer when the stored aggregate audit
  already owns the result;
- convert signed-support diagnostics into uncertainties;
- infer physics agreement/disagreement from a presentation choice.

---

## 13. Memory updates

Update only warranted durable memory.

Allowed:

```text
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/LEARNINGS.md
docs/memory/roadmap/STATUS.md
docs/memory/phases/e8-4-fix5-4-current-lineage-identity-audit.md
docs/memory/investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md
docs/memory/phases/e8-4-fix5-5-current-lineage-visualization-clarity-task-contract.md
docs/memory/phases/e8-4-fix5-5-current-lineage-visualization-clarity.md
docs/memory/manifest.json
```

Do not change `USER.md` or `CURRENT_HANDOFF.md` unless their owned durable
meaning actually changed. It should not change for this task.

### Required memory reconciliation

Record the now-pushed numerical source identity:

```text
761fbb6c03d2d7a10bb911cf84e9ba898496fab6
```

and that independent ChatGPT actual-diff and pushed-state review passed.

Remove stale wording that calls the Fix.5.4 source edits uncommitted or still
awaiting numerical actual-diff review.

Do not claim farm validation for Fix.5.4.

### New phase record

Create:

```text
docs/memory/phases/e8-4-fix5-5-current-lineage-visualization-clarity.md
```

During implementation it is `ACTIVE`.

If local deterministic implementation and tests pass, it may become:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

with explicit wording that actual PDF/ROOT rendering remains farm-only.

### CURRENT / NEXT

At local completion, CURRENT should state the dependency chain without
authorizing the farm:

```text
Fix.5.5 actual-diff review
-> user commit/push
-> pushed-state review
-> farm-readiness review of the tracked owner/packaging path
-> one narrow Q4p4W2p74 / Left / lowe farm gate only after that path is ready
```

The known prior missing-owner-ZIP issue remains separate. Do not repair it in
this contract.

Regenerate `docs/memory/manifest.json` for the intended versioned candidate.
Exclude unrelated user-owned untracked memory files.

---

## 14. Local validation

Use the repository-selected Python interpreter.

Run syntax checks for every changed/new Python file.

Run at minimum:

```text
testing/test_e8_4_fix5_5_visualization.py          # if created
testing/test_full_background_subtraction_plots.py
testing/test_e8_4_production_impact_audit.py
testing/test_e8_4_fix5_4_identity_audit.py
```

plus directly relevant E.8.3/E.8.4 regressions discovered by source audit.

Run:

```bash
git diff --check
```

Regenerate/check memory manifest for the intended candidate and run:

```bash
<PYTHON> -B tools/check_memory_health.py --root .
```

Report exact pass/skip counts.

Do not claim ROOT/PyROOT/full `main.py` or procedure-PDF validation from local
tests.

---

## 15. Diff audit

Before stopping:

```bash
git status --short --untracked-files=all
git diff --stat
git diff --check
```

The unrelated user-owned:

```text
docs/memory/phases/workflow-continuity-hardening-task-contract.md
```

may remain untracked and must not be included in this candidate.

The temporary:

```text
kaonlt_review.diff
```

may exist locally but must not be intended for commit.

Produce one complete replacement `kaonlt_review.diff` relative to committed
starting HEAD:

```text
761fbb6c03d2d7a10bb911cf84e9ba898496fab6
```

It must contain:

- complete tracked diff;
- complete additions for every intended new file;
- no unrelated user-owned untracked file.

Do not stage merely to make the diff reviewable.

---

## 16. Acceptance criteria

PASS only if:

1. starting HEAD is exact;
2. only presentation/test/warranted-memory scope changes;
3. blue baseline is independently visible under near-overlap conditions;
4. Method A remains independently visible;
5. visual distinction does not rely only on color;
6. actual histogram values/errors are unchanged;
7. no normalization or physics changes;
8. Fix.5.4 signed-support values are visibly presented from the stored audit;
9. current-F.6.3 aggregate values are presented from the stored audit, not
   recomputed;
10. historical E.8.3 is visibly labeled as historical/not-current lineage;
11. absolute SIMC blocker remains intact;
12. no page IDs/page count change solely from this task;
13. deterministic tests pass with documented environment-only skips;
14. memory/manifest are scoped correctly;
15. no owner/profile/collector repair occurs;
16. no farm run, commit, or push occurs.

---

## 17. Farm boundary

This task does not establish:

- ROOT/PyROOT rendering quality;
- whether blue/magenta styles are visually satisfactory in the real PDF;
- actual farm signed-support values;
- whether cancellation explains the previously observed large fractional yield
  shifts;
- Fix.5.4 or Fix.5.5 runtime closure;
- absolute SIMC amplitude comparability;
- owner ZIP packaging repair;
- Method-A promotion.

Those require fresh farm evidence or a separate source authority.

---

## 18. Hard stop

After implementation, deterministic local checks, warranted memory updates,
manifest regeneration, and complete diff preparation:

**STOP.**

Do not commit.

Do not push.

Do not repair the owner/profile/collector.

Do not run the farm.

Return the complete `kaonlt_review.diff` to ChatGPT for independent actual-diff
review.
