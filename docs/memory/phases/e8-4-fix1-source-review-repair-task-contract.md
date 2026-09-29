# KaonLT E.8.4.Fix.1 — Source-Review Runtime-Identity and Display Repair

## Purpose

Repair three narrow defects found by independent ChatGPT review of
`kaonlt_review(20260928-175335).diff` for the local E.8.4 candidate.

This is **not** a redesign of E.8.4 and does not reopen F.6.3 or any accepted
upstream science. Preserve the presentation-only consumer architecture exactly:

```text
completed baseline yield + existing F.6.3 sidecar
    -> E.8.4 consumer/validator
    -> detached display objects
    -> existing post-yield pair-safe procedure rerender
```

The repair is required because the current local candidate cannot become
available on the real runtime path as written and one required comparison is
not actually visible on its page.

Do not run the Jefferson Lab farm. Local checks are not ROOT/PyROOT/full-runtime
validation.

---

# 1. Exact committed base and expected dirty candidate

Required branch:

```text
test
```

Required committed HEAD:

```text
d84f9610427cf38ccaaa8d6d8b14a8a995ecbac9
Add F6.3 parallel Method-A full procedure
```

This repair starts from the existing uncommitted E.8.4 candidate reviewed in:

```text
kaonlt_review(20260928-175335).diff
```

Before editing, run:

```bash
git branch --show-current
git rev-parse HEAD
git log -1 --oneline
git status --short
```

The intended dirty candidate paths are limited to:

```text
docs/memory/CURRENT.md
docs/memory/decisions/e8-full-analysis-procedure-roadmap.md
docs/memory/manifest.json
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/roadmap/STATUS.md
src/cuts/full_background_subtraction_plots.py
docs/memory/phases/e8-4-production-impact-audit-task-contract.md
docs/memory/phases/e8-4-production-impact-audit.md
testing/test_e8_4_production_impact_audit.py
```

plus this new repair contract:

```text
docs/memory/phases/e8-4-fix1-source-review-repair-task-contract.md
```

Local-only `AGENTS.md` / `.codex/` may also exist. The prior temporary review
bundle `kaonlt_review(20260928-175335).diff` may remain untracked and read-only.
Do not stage it.

Any other tracked modification or unrelated untracked file is a hard stop.
Do not reset, stash, clean, commit, push, or overwrite unrelated user work.

---

# 2. Mandatory startup reading

Read in this exact order:

1. local root `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read only:

- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/phase-f6-method-a-production-promotion.md`
- `docs/memory/phases/e8-2-baseline-stage-audit.md`
- `docs/memory/phases/e8-3-detached-method-a-reweighting-audit.md`
- `docs/memory/phases/f6-3-parallel-full-procedure-method-a.md`
- `docs/memory/phases/e8-4-production-impact-audit-task-contract.md`
- `docs/memory/phases/e8-4-production-impact-audit.md`
- this Fix.1 contract

Inspect the actual current candidate source and test before editing.

---

# 3. Independent-review findings that must be repaired

## Finding 1 — E.8.2 and F.6.3 use different epsilon token semantics

The current candidate derives the F.6.3 setting identity from the E.8.2 display
payload as:

```text
<setting>-<epsilon>
```

but the actual producers intentionally use different epsilon representations:

- `main.py` / `inpDict["EPSSET"]` and therefore the E.8.2 source use semantic
  runtime values `low` or `high`;
- `_f6_3_setting_id(...)` converts those to canonical filename/setting tokens
  `lowe` or `highe`.

Therefore a real Left/low run presents:

```text
E.8.2 identity: Left + low
F.6.3 identity: Left-lowe
```

and the current local E.8.4 candidate incorrectly constructs `Left-low`, causing
`e8_4_setting_identity_mismatch` before an available E.8.4 page can exist.

### Required repair

Keep the E.8.2 payload's semantic epsilon unchanged for display/provenance.
Introduce one narrow explicit conversion for comparison with the frozen F.6.3
setting token:

```text
low  -> lowe
high -> highe
```

Do not modify `main.py`, E.8.2 producer state, F.6.3 producer state, artifact
filenames, or any upstream epsilon semantics.

Unexpected E.8.2 epsilon values must fail closed with a literal E.8.4 identity
reason; do not guess or use substring heuristics.

The available E.8.4 payload should retain:

- semantic display epsilon (`low` / `high`); and
- canonical F.6.3 selected setting ID (`Left-lowe`, etc.).

Tests must use the real runtime semantic token, not `lowe` as the E.8.2 fixture.
Cover both low->lowe and high->highe mapping.

---

## Finding 2 — E.8.2 wide-MM geometry is not the F.6.3 analysis-MM geometry

The current candidate requires all F.6.3 branch histograms to have the same
`mm_edges` as the E.8.2 baseline presentation.

That is incorrect for the real source path:

- E.8.2 intentionally audits the **wide** missing-mass chain built from
  `H_MM_nosub_*`, whose production histogram range is the diagnostic
  `BG_OPT_MM_PLOT_MIN .. BG_OPT_MM_PLOT_MAX` (currently 0.7..1.30 GeV);
- F.6.3 intentionally stores the actual **analysis/yield** branch objects built
  from `H_MM_DATA`, `H_pion_subtraction_template_MM`, and the corresponding
  Method-A branch. `H_MM_DATA` uses `mm_min .. mm_max` (currently 1.10..1.16
  GeV for the production runtime defaults).

These are different authoritative views of the same analysis and must not be
forced to share histogram binning.

The current test accidentally hides this defect by giving E.8.2 and F.6.3 the
same fake MM edges.

### Required repair

E.8.4 must continue to use E.8.2 only for:

- setting/epsilon identity;
- canonical t/phi inventory and geometry;
- Lambda integration-window identity; and
- baseline `Y0`, statistical-error, and total-error identity.

Do **not** compare F.6.3 MM histogram edges to E.8.2 wide-MM edges.

Instead, validate the F.6.3 branch histogram geometry internally:

1. each required F.6.3 histogram
   `pion_input`, `B_pi_0`, `B_pi_A`, `MM_0`, `MM_A` must have a valid finite,
   strictly ordered 1D axis;
2. those five histograms must have identical bin edges within one child;
3. all available F.6.3 canonical children must use one identical analysis-MM
   edge vector;
4. expose that F.6.3 analysis edge vector as the E.8.4 payload's `mm_edges`
   (or an equivalently explicit `analysis_mm_edges` field), not the E.8.2
   wide-MM edge vector.

The Lambda integration window must still match exactly between the E.8.2 and
F.6.3 sidecars. Do not alter the window or rebin either source.

A real wide-E.8.2 / narrow-F.6.3 fixture must produce an **available** E.8.4
payload when all other identities match.

A binning mismatch among the F.6.3 branch outputs themselves must fail closed.

No normalization, scaling, rebinning, interpolation, smoothing, clipping, or
child renormalization is authorized.

---

## Finding 3 — the required common pion input is overwritten on the page

The current pion-consequence renderer draws:

1. `pion_input` with `Draw("hist e")`;
2. `B_pi_0` through the same helper, again with `Draw("hist e")`;
3. `B_pi_A` with `Draw("hist e same")`.

The second non-`same` draw replaces the first display, so the required common
proton-cleaned pion input is not actually visible even though the note says it
is.

### Required repair

On each pion-consequence child panel, visibly retain all three already-produced
objects:

```text
common pion_input
B_pi_0
B_pi_A
```

Use the first object to establish the axes and draw the latter two with an
explicit `same` option. A legend or existing color/note convention may identify
the three curves, but the input must genuinely remain on the pad.

Do not alter histogram values to improve visibility.

Also replace identity checks of the form:

```python
if None in (...):
```

with explicit `is None` logic in the E.8.4 renderer so PyROOT object equality is
never invoked merely to test absence.

---

# 4. Preserved architecture

Do not alter the accepted E.8.4 architecture beyond the three repairs above.

Preserve:

- `full_background_subtraction_e8_4/v1` unless a schema change is genuinely
  required by a field semantic change; if the existing field can simply be
  corrected before any runtime acceptance, keep v1;
- exact F.6.3 source schema/branch-role validation;
- live-cache parity requirement;
- all non-production flags;
- exact t/phi child inventory checks;
- exact E.8.2 baseline `Y0`, statistical-error, and total-error identity checks;
- detached clones of all retained histograms;
- `DeltaY = YA - Y0` and defined `DeltaY/Y0` only;
- explicit undefined fractional shift for `Y0 == 0`;
- no `DeltaY` or fractional-shift uncertainty;
- display-only `Delta B_pi` and `Delta K` from clones;
- Method B numerical absence;
- dormant Fit 1/Fit 2;
- existing E.8.3 -> E.8.4 -> terminal-handoff ordering;
- existing post-yield pair-safe PDF/page-manifest transaction;
- unavailable E.8.4 as an optional presentation state that does not invalidate
  the public baseline result;
- final E.8 and F.6.4 still blocked.

---

# 5. Allowed source/test scope

Source edit remains limited to:

```text
src/cuts/full_background_subtraction_plots.py
```

Test edit remains limited to:

```text
testing/test_e8_4_production_impact_audit.py
```

Warranted memory/history edits are limited to the existing E.8.4 candidate
records plus this repair contract and manifest:

```text
docs/memory/CURRENT.md
docs/memory/decisions/e8-full-analysis-procedure-roadmap.md
docs/memory/manifest.json
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/phases/e8-4-production-impact-audit.md
docs/memory/phases/e8-4-production-impact-audit-task-contract.md
docs/memory/phases/e8-4-fix1-source-review-repair-task-contract.md
docs/memory/roadmap/STATUS.md
```

Do not edit any other source/test/memory path.

---

# 6. Explicitly frozen files

Remain byte-unchanged:

```text
src/main.py
src/binning/calculate_yield.py
src/cuts/rand_sub.py
src/cuts/pion_component_subtraction.py
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
src/utility/background_config.py
```

Also do not modify Method-B files, accepted F-stage calculators/validators,
launchers, collectors, profiles, wrappers, cross-section code, or production
subtraction logic.

---

# 7. Required deterministic test repairs/additions

Update the focused E.8.4 suite so it tests the **actual runtime contracts**, not
the accidental fake identities that masked Findings 1 and 2.

At minimum add/repair coverage for:

1. E.8.2 `epsilon="low"` + F.6.3 `selected_setting_id="Left-lowe"` -> available.
2. E.8.2 `epsilon="high"` + an appropriate `*-highe` F.6.3 setting -> available.
3. Unsupported semantic epsilon -> fail closed.
4. E.8.2 wide MM edges and F.6.3 narrow analysis MM edges -> available when all
   other identities match.
5. A bin-edge mismatch among `pion_input/B_pi_0/B_pi_A/MM_0/MM_A` -> fail closed.
6. A bin-edge mismatch in one later canonical child relative to the first
   F.6.3 child -> fail closed.
7. E.8.4 output MM edges are the F.6.3 analysis edges, not E.8.2 wide edges.
8. Baseline statistical-error mismatch -> fail closed.
9. Baseline total-error mismatch -> fail closed.
10. Builder and all E.8.4 renderers leave source histogram **contents and bin
    errors** unchanged.
11. Pion-consequence draw order/options prove the common input remains drawn and
    both B_pi curves are overlaid with `same` semantics.
12. E.8.4 unavailable + successful unavailable-page rendering still allows the
    existing finalizer transaction to return baseline finalization success while
    exposing `e8_4_available=False` and the literal reason.
13. Existing renderer-failure rollback/recovery test remains intact.
14. Existing zero-Y0 gap-safe point behavior remains intact.
15. Existing static producer/Method-B/uncertainty guards remain intact.

Prefer at least one fixture with multiple t parents and multiple phi children;
a full 3x9 geometry fixture is useful if it stays lightweight, but do not
recreate production calculations in the test.

---

# 8. Required regression checks

Rerun the original E.8.4 required checks after Fix.1:

```bash
python -B -m py_compile \
  src/cuts/full_background_subtraction_plots.py \
  testing/test_e8_4_production_impact_audit.py

python -B -m unittest testing.test_e8_4_production_impact_audit -v
python -B -m unittest testing.test_e8_2_baseline_stage_audit -v
python -B -m unittest testing.test_e8_3_detached_method_a_reweighting_audit -v
python -B -m unittest testing.test_f6_3_parallel_full_procedure_method_a -v
python -B -m unittest testing.test_full_background_subtraction_plots -v
python -B -m unittest testing.test_pion_hgcer_phase_e_runtime_contract -v
python -B -m unittest testing.test_pion_hgcer_phase_f_runtime_contract -v
```

Then rerun memory integrity:

```bash
python -B tools/update_memory_manifest.py --root . --write
python -B tools/update_memory_manifest.py --root . --check
python -B tools/check_memory_health.py --root .
python -B tools/memory_bootstrap.py --root . --json
python -B -m unittest testing.test_memory_health -v

git -c core.safecrlf=false diff --check
```

Record the actual interpreter/version. PyROOT-unavailable local skips remain
local limitations, not runtime validation.

Do not repair unrelated pre-existing failures.

---

# 9. Memory/history requirements

Append a clearly labeled **E.8.4.Fix.1 — independent source-review repair**
section to:

```text
docs/memory/phases/e8-4-production-impact-audit.md
```

Record the three findings and exact repairs:

- semantic `low/high` -> F.6.3 `lowe/highe` identity normalization;
- wide E.8.2 MM geometry kept distinct from narrow F.6.3 analysis-MM geometry;
- actual common pion input retained visibly on the consequence page.

Also record the strengthened deterministic tests and results.

Keep E.8.4 status:

```text
ACTIVE
```

until independent ChatGPT review of the refreshed complete diff passes.

Do not self-promote to `SOURCE REVIEWED`, `DEVELOPMENT COMPLETE, FARM VALIDATION
PENDING`, or `CLOSED / RUNTIME VALIDATED`.

Update CURRENT/roadmaps only as needed to say Fix.1 is ACTIVE and the exact NEXT
is refreshed independent source review. Final E.8 and F.6.4 remain BLOCKED.

Regenerate `docs/memory/manifest.json` after all versioned memory edits.

---

# 10. Diff audit and refreshed review bundle

Before stopping, inspect:

```bash
git status --short
git diff --stat
git diff -- src/cuts/full_background_subtraction_plots.py
git diff -- testing/test_e8_4_production_impact_audit.py
git diff -- docs/memory
git -c core.safecrlf=false diff --check
```

Confirm all frozen files remain byte-unchanged and no scope expansion occurred.

Create a **new** complete temporary review bundle with a fresh timestamp that
contains:

- the entire tracked diff from committed HEAD
  `d84f9610427cf38ccaaa8d6d8b14a8a995ecbac9`; and
- complete `git diff --no-index /dev/null ...` sections for every intended
  untracked/new file not otherwise present in ordinary `git diff` output.

Do not stage merely to make the review bundle.

Do not modify or stage the older `kaonlt_review(20260928-175335).diff`.

---

# 11. Acceptance criteria

Fix.1 is locally implementation-complete only when:

1. committed base remains exactly `d84f961...`;
2. only the allowlisted E.8.4 candidate/repair paths changed;
3. E.8.2 semantic epsilon `low/high` maps explicitly to F.6.3
   `lowe/highe` for identity comparison;
4. real `Left / low` source identity can pass without becoming `Left-low`;
5. E.8.2 wide-MM edges are not required to equal F.6.3 analysis-MM edges;
6. F.6.3 branch histograms have internally consistent analysis-MM geometry and
   E.8.4 exposes that geometry;
7. no source histogram is rebinned/scaled/normalized or mutated;
8. the common pion input is genuinely visible together with `B_pi_0` and
   `B_pi_A`;
9. source histogram contents and errors are regression-protected;
10. baseline Y0/statistical/total identity remains fail-closed;
11. Method B remains numerically absent and Fit 1/Fit 2 dormant;
12. pair-safe finalization/recovery semantics remain unchanged;
13. unavailable E.8.4 remains optional and does not invalidate baseline success;
14. required focused/regression/memory checks pass;
15. `git diff --check` passes;
16. memory truthfully says E.8.4.Fix.1 is ACTIVE pending independent review;
17. no commit, push, farm run, final-E.8 closure, F.6.4 work, or Method-A
    promotion occurred.

---

# 12. Hard stop

Stop after the narrow Fix.1 implementation, deterministic local checks, memory
update, diff audit, and refreshed review-bundle creation.

Do not commit or push.
Do not run the farm.
Do not begin final E.8 or F.6.4.
Do not change F.6.3 production-side source.

The exact NEXT is:

```text
independent ChatGPT review of the refreshed complete cumulative E.8.4 + Fix.1 diff/source/runtime path
```
