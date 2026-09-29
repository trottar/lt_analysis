# KaonLT E.8.4 — Baseline-versus-Method-A Production-Impact Audit

## Purpose

Implement **E.8.4**, the presentation-only audit of the actual F.6.3 parallel
full-analysis Method-A branch against the unchanged baseline branch.

This is the next approved source phase after pushed/source-reviewed F.6.3. It
must **consume the F.6.3 two-branch runtime sidecar and never construct, rerun,
refit, renormalize, or otherwise reproduce either branch**.

The implementation must make the actual production-impact comparison visible in
the existing full-background-subtraction procedure PDF while preserving all
accepted production/scientific ownership boundaries.

This task is **source development only**. Do not run the Jefferson Lab farm.
Do not claim ROOT/PyROOT, full `main.py`, procedure-PDF runtime, farm, or
production validation from local checks.

---

# 1. Exact starting state

Required branch:

```text
test
```

Required committed starting HEAD:

```text
d84f9610427cf38ccaaa8d6d8b14a8a995ecbac9
Add F6.3 parallel Method-A full procedure
```

Before editing, run:

```bash
git branch --show-current
git rev-parse HEAD
git log -1 --oneline
git status --short
```

Hard stop if the branch or committed HEAD differs.

The worktree may contain only:

- this task-contract file as the intentional new repository file; and
- local-only/untracked `AGENTS.md` and/or `.codex/` if they already exist.

Any other tracked modification or unrelated untracked file is a hard stop.
Do not reset, stash, clean, checkout over, commit, push, or delete unrelated
user work.

`AGENTS.md` and `.codex/` are local-only and must never be staged merely because
they exist.

---

# 2. Mandatory startup reading

Read in this exact order:

1. local root `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read only the task-relevant records:

- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/phase-f6-method-a-production-promotion.md`
- `docs/memory/phases/e8-2-baseline-stage-audit.md`
- `docs/memory/phases/e8-3-detached-method-a-reweighting-audit.md`
- `docs/memory/phases/f6-3-parallel-full-procedure-method-a.md`
- this contract

Before editing, inspect the current source/runtime path in:

- `src/binning/calculate_yield.py`
- `src/cuts/pion_hgcer_method_a_parallel_full_procedure.py`
- `src/cuts/pion_component_subtraction.py`
- `src/cuts/full_background_subtraction_plots.py`
- `src/cuts/rand_sub.py`
- `src/main.py`

The first three producer/scientific files are inspection-only for this task.

---

# 3. Authoritative status and dependency boundary

Preserve the current status exactly at task start:

```text
F.6.2 / F.6.2.Fix.5 — CLOSED / RUNTIME VALIDATED
E.8                 — ACTIVE
E.8.2               — SOURCE REVIEWED
E.8.3               — SOURCE REVIEWED
F.6.3               — SOURCE REVIEWED
E.8.4               — NEXT
final E.8            — BLOCKED
F.6.4               — BLOCKED
```

The governing sequence remains:

```text
E.8.2 baseline full-analysis audit
-> E.8.3 detached Method-A reweighting audit
-> F.6.3 parallel full procedure plus Method A
-> E.8.4 baseline-versus-Method-A production-impact audit
-> final E.8 closure
-> F.6.4 explicit production-promotion decision
```

E.8.4 is presentation-only. It is not a production-promotion decision.

---

# 4. Frozen scientific/runtime ownership

The following are frozen and must not change:

- baseline production branch;
- random subtraction;
- dummy subtraction;
- slow-proton subtraction/weights/gates;
- pion components, fits, windows, amplitudes, priors, templates and baseline
  `w0`;
- F.3/F.4/F.5/F.6.1/F.6.2 scientific algorithms, artifacts, fingerprints and
  accepted authorities;
- F.6.3 event-factor authority and private branch construction;
- canonical `(t,phi)` binning;
- pruning;
- normalizations;
- SIMC;
- efficiencies and acceptance;
- yield extraction/error arithmetic;
- L/T separation and cross-section formulas;
- Method B, which remains diagnostic/cross-check only and numerically absent;
- active background profile `no_empirical_residual`;
- dormant legacy empirical residual Fit 1 / Fit 2.

Do **not** modify Method-A factors, event weights, F.6.3 authority, or public
production results to make a comparison/checker/plot look better.

The Method-A minus baseline shift is a **correction effect**, not automatically
a systematic uncertainty. E.8.4 must not invent an uncertainty for that shift.

---

# 5. Current source/runtime path that E.8.4 must consume

The reviewed F.6.3 producer already provides everything E.8.4 needs.

After the ordinary baseline yield extraction, `calculate_yield_data(...)`
attaches:

```text
hist["_f6_3_parallel_method_a_source"]
```

with source schema:

```text
f6_3_parallel_method_a_source/v1
```

For an available setting it contains, per canonical `(t,phi)` child:

```text
pion_input
B_pi_0
B_pi_A
MM_0
MM_A
Y0
Y0_statistical_error
Y0_total_error
YA
YA_statistical_error
YA_total_error
```

plus the exact Lambda integration window, selected setting identity, authority,
and explicit non-production flags.

The F.6.3 source is built from the already-existing authoritative cache after the
public baseline result exists. E.8.4 must not trigger another tree traversal,
F.4 reconstruction, factor lookup, pion-template fill, or yield calculation.

`main.py` already calls the established post-yield procedure finalizer
immediately after `find_yield_data(...)`. That finalizer already rerenders the
procedure PDF/page-manifest pair through temporary paths with recovery of the
preliminary pair on handled failure.

**Preserve this runtime path.** E.8.4 should join that existing post-yield
rerender rather than introducing a new PDF lifecycle or a second main-level
finalization stage.

---

# 6. Allowed implementation files

Scientific/presentation source change is limited to:

```text
src/cuts/full_background_subtraction_plots.py
```

Add one focused test module:

```text
testing/test_e8_4_production_impact_audit.py
```

Warranted memory/history files may include only:

```text
docs/memory/CURRENT.md
docs/memory/decisions/e8-full-analysis-procedure-roadmap.md
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/phases/e8-4-production-impact-audit.md
docs/memory/phases/e8-4-production-impact-audit-task-contract.md
docs/memory/roadmap/STATUS.md
docs/memory/manifest.json
```

Do not edit any other source/test/memory file unless a concrete blocker makes
that impossible. If a blocker requires scope expansion, stop and report it
rather than expanding scope autonomously.

---

# 7. Explicitly frozen implementation files

Do not modify:

```text
src/main.py
src/binning/calculate_yield.py
src/cuts/rand_sub.py
src/cuts/pion_component_subtraction.py
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
src/utility/background_config.py
```

Do not modify Method-B files, accepted F-stage analyzers/validators, production
subtraction code, collector/profile infrastructure, bundle wrappers, or launcher
scripts.

The existing finalizer name may remain historically E.8.2-named. Do not refactor
or rename it merely for aesthetics.

---

# 8. E.8.4 presentation source contract

Add a presentation schema such as:

```text
full_background_subtraction_e8_4/v1
```

and a narrow builder in `full_background_subtraction_plots.py` that consumes:

1. the existing F.6.3 sidecar; and
2. the already-built E.8.2 baseline presentation payload for geometry/baseline
   identity cross-checks.

The builder is a **consumer/validator**, not a producer.

## 8.1 Fail-closed source identity

For an available E.8.4 payload require:

- F.6.3 source schema exactly `f6_3_parallel_method_a_source/v1`;
- `available is True`;
- selected setting exactly matches the current procedure setting;
- branch role exactly `parallel_nonproduction_method_a_full_analysis`;
- the source authority is present;
- live-cache parity is explicitly passed;
- exact required non-production flags:
  - `baseline_production_mutated is False`;
  - `production_promotion_performed is False`;
  - `method_b_numerical_dependency is False`;
  - `empirical_residual_used is False`;
  - `event_correction_persisted is False`;
  - `canonical_child_renormalization_performed is False`;
  - `baseline_public_output_unchanged is True`.

Do not reopen or re-run the F.6.3 authority calculation. E.8.4 validates the
published sidecar contract and copies the authority/provenance needed for the
presentation.

If the F.6.3 source is missing, unavailable, malformed, wrong-setting, stale in
geometry, or violates its flags, return an explicit unavailable E.8.4 payload.
Do not fall back to E.8.3, F.5, F.6.1, unity factors, reconstructed Method-A
histograms, or the baseline branch.

## 8.2 Geometry and baseline identity

Use the E.8.2 presentation payload as the current-baseline geometry/yield
cross-check only. Do not use it to construct Method-A output.

Require:

- identical setting and epsilon identity;
- identical Lambda integration window;
- exact one-to-one canonical `(t_index, phi_index)` inventory;
- matching `t_low/t_high` and `phi_low/phi_high` for each child;
- no duplicate, missing, or unexpected child;
- for every available F.6.3 child, `Y0`, `Y0_statistical_error`, and
  `Y0_total_error` match the corresponding authoritative E.8.2 baseline values.

An inventory/geometry/baseline mismatch makes E.8.4 unavailable. Do not repair
it in presentation code.

## 8.3 Histogram ownership

Required source histograms are:

```text
pion_input
B_pi_0
B_pi_A
MM_0
MM_A
```

Every histogram retained by the E.8.4 presentation payload must be a detached
clone. Never keep a mutable alias to the F.6.3 source object.

Reuse the existing display-clone/ROOT-ownership helpers in
`full_background_subtraction_plots.py`; do not create a new ownership system.

The E.8.4 builder must not mutate source histograms and must not call production
fit/weight/yield machinery.

## 8.4 Scalar impact quantities

For each canonical child copy the authoritative `Y0` and `YA` values and their
producer-owned statistical/total errors.

Presentation may derive only:

```text
DeltaY = YA - Y0
DeltaY_over_Y0 = DeltaY / Y0
```

The fractional quantity is defined only for a finite, nonzero `Y0`.

If `Y0 == 0`, retain an explicit undefined status/reason. Do not use an epsilon
denominator, cap, clip, winsorization, or replacement value.

Do not derive an uncertainty for `DeltaY` or `DeltaY/Y0`. In particular, do not
quadrature-combine the baseline and Method-A errors as though the two branches
were independent.

Nonfinite required authoritative values make the E.8.4 payload unavailable.

---

# 9. Display-only signed histogram differences

E.8.4 must show:

```text
Delta B_pi = B_pi_A - B_pi_0
Delta K    = MM_A - MM_0
```

These are display comparisons of already-produced branch outputs, not new
production corrections.

Construct them only from detached display clones. Do not write them back to the
F.6.3 source or use them downstream.

No normalization, scaling to common area, clipping, smoothing, interpolation,
rebinning, or child renormalization is allowed.

The builder itself should remain an authority/geometry/detachment reader. Keep
histogram subtraction in the renderer/display helper where practical so the
payload remains a detached representation of the actual branch outputs.

---

# 10. Required E.8.4 procedure pages

Append E.8.4 to the **existing** full-background-subtraction procedure PDF.
Do not create another PDF file or lifecycle.

Placement must be:

```text
... E.8 / E.8.2 baseline pages
-> E.8.3 detached Method-A reweighting pages
-> E.8.4 actual F.6.3 production-impact pages
-> terminal E.8 handoff page
```

## 10.1 Authority/context page — one per setting

Show at minimum:

- `E.8.4 baseline-versus-Method-A production-impact audit`;
- current setting identity;
- F.6.3 source schema/branch role;
- Lambda integration window;
- canonical child count;
- live-cache parity status;
- non-production flags;
- explicit statements that Method B is numerically absent, empirical residual
  Fit 1/Fit 2 are inactive, and no production promotion is performed;
- statement that the Method-A minus baseline shift is a correction effect, not
  automatically a systematic uncertainty.

## 10.2 Per-canonical-t pion consequence page

For all canonical phi children of each t parent, show the common
proton-cleaned pion input and overlay:

```text
B_pi_0
B_pi_A
```

The page must make clear that both branches use the same input and differ only
through the already-produced F.6.3 Method-A pion template.

## 10.3 Per-canonical-t clean-kaon/final-MM page

For every populated canonical cell overlay:

```text
MM_0 == K_0
MM_A == K_A
```

Draw the same Lambda signal/integration window on the comparison.

Each child title/note must identify the canonical t parent, phi index/range,
`Y0`, and `YA` sufficiently to audit the cell.

## 10.4 Per-canonical-t signed-difference page

For every child show both signed display differences:

```text
Delta B_pi = B_pi_A - B_pi_0
Delta K    = MM_A - MM_0
```

Do not force them to agree with an expected algebraic relation and do not alter
source content to close a checker. They are direct visual consequences of the
stored branch outputs.

## 10.5 Per-canonical-t yield-impact page

For the nine phi children show:

- `Y0(phi)` and `YA(phi)` using the stored producer-owned values;
- `DeltaY(phi)`;
- `DeltaY/Y0(phi)` only where defined.

Use gap-safe point rendering for undefined fractional values; never connect a
line across an undefined child.

Do not display an invented uncertainty on `DeltaY` or `DeltaY/Y0`.

## 10.6 Setting-level t/phi impact summary

Add one setting-level summary that exposes the t dependence without constructing
new integrated yields. A canonical t-by-phi table/map of `DeltaY` and/or
`DeltaY/Y0` is acceptable.

Do **not** sum child yields to invent a new t-integrated observable merely for
this page.

The roadmap's all-five-setting summary remains later work after the narrow farm
gate and subsequent broadening. Do not build a multi-setting collector in this
task.

---

# 11. Unavailable behavior

If the F.6.3 source is unavailable or invalid, render an explicit E.8.4
unavailable page containing the literal reason.

An unavailable optional F.6.3/E.8.4 branch must **not** abort, alter, or
invalidate the public baseline yield calculation.

The existing post-yield procedure rerender may still complete successfully with
an E.8.4-unavailable page. This lets farm evidence distinguish:

- baseline analysis success; from
- E.8.4 branch availability/acceptance.

Never replace an unavailable Method-A comparison with E.8.3 detached evidence
or baseline-equals-Method-A placeholders.

A true renderer/filesystem failure remains subject to the existing pair-safe
E.8.2 finalization/recovery semantics.

---

# 12. Existing post-yield finalizer integration

Extend the existing post-yield full-background-subtraction finalization in
`full_background_subtraction_plots.py` only.

The finalizer already:

1. consumes retained detached Step-3 render state;
2. builds the E.8.2 post-yield payload;
3. rerenders the complete PDF to temporary paths;
4. writes the page manifest;
5. preserves/restores the preliminary PDF/manifest pair on handled failures.

Add E.8.4 by:

- reading `hist.get("_f6_3_parallel_method_a_source")` only after baseline
  yields are complete;
- building the E.8.4 presentation against the already-built E.8.2 payload;
- passing that payload into the existing renderer invocation;
- preserving the existing temporary/recovery/install transaction unchanged.

Do not add a second finalizer call in `main.py` and do not alter the ordinary
Step-3 `rand_sub.py` route.

The returned finalization status should expose E.8.4 availability/reason so farm
logs can distinguish a rendered unavailable page from a successful available
comparison, while retaining existing E.8.2 status fields used by current tests.

---

# 13. Renderer/API compatibility

Extend `render_full_background_subtraction_procedure_pages(...)` with one
optional keyword-only E.8.4 payload argument.

Existing callers that omit it must retain their previous page inventory and
behavior. This is important for deterministic regressions and non-E.8.4 helper
usage.

On the production post-yield route, the finalizer supplies an explicit available
or unavailable E.8.4 payload.

Add stable page-manifest IDs, for example:

```text
full_background.e8_4.authority
full_background.e8_4.pion_consequence
full_background.e8_4.final_mm
full_background.e8_4.signed_difference
full_background.e8_4.yield_impact
full_background.e8_4.setting_summary
full_background.e8_4.unavailable
```

Use canonical t scope where applicable. These pages are presentation-only and
must not be marked as production-authoritative.

---

# 14. Narrow stale-text reconciliation

Within `full_background_subtraction_plots.py`, update only E.8/E.8.3/handoff
text that would become factually stale after F.6.3/E.8.4 is present.

In particular, wording that says the actual production-impact comparison
"belongs to later F.6.3 work" must no longer imply F.6.3 is unimplemented.

Replacement wording must preserve these facts:

- E.8.3 remains detached evidence;
- F.6.3 alone constructs the private parallel branch;
- E.8.4 consumes that branch for the actual impact audit;
- Method A is still not production-promoted;
- F.6.4 remains the only promotion decision.

Do not otherwise rewrite accepted E.8.1/E.8.2/E.8.3 presentation content.

---

# 15. Focused deterministic tests

Create:

```text
testing/test_e8_4_production_impact_audit.py
```

Use lightweight fake histogram/ROOT fixtures where needed. Do not require the
Jefferson Lab farm for deterministic source contracts.

At minimum cover:

## Positive contracts

1. A valid F.6.3 source plus matching E.8.2 payload produces an available E.8.4
   payload.
2. Exact selected-setting, Lambda-window, t/phi geometry, child inventory and
   baseline `Y0` identity are retained.
3. All canonical children are retained in canonical order.
4. `DeltaY` and defined `DeltaY/Y0` are numerically correct.
5. An all-one/no-change branch gives zero `DeltaY`, zero defined fractional
   shift, and zero display differences without mutating baseline objects.
6. Source histogram identity/content/errors are unchanged by builder and
   renderer.
7. E.8.4 page ordering is after E.8.3 and before the terminal E.8 handoff.
8. Available E.8.4 pages append stable manifest IDs/scopes.
9. The existing post-yield finalizer passes the available E.8.4 payload through
   the existing one-rerender transaction.

## Negative/fail-closed contracts

10. Missing F.6.3 source -> explicit E.8.4 unavailable payload/page.
11. F.6.3 `available=False` -> preserve literal reason; no fallback.
12. Wrong schema, setting, branch role, or non-production flag -> unavailable.
13. Missing/false live-cache parity authority -> unavailable.
14. Lambda-window mismatch -> unavailable.
15. Missing, duplicate, extra, or geometrically mismatched child -> unavailable.
16. Baseline F.6.3 `Y0`/error mismatch versus E.8.2 -> unavailable.
17. Missing required histogram -> unavailable.
18. Nonfinite required yield/error -> unavailable.
19. `Y0 == 0` -> fractional shift explicitly undefined; no epsilon fallback.
20. Undefined fractional children are rendered gap-safe; no line bridges the
    missing point.
21. E.8.4 unavailable does not make the public baseline result unavailable and
    does not replace E.8.2/E.8.3 evidence.
22. Renderer failure still exercises the existing preliminary-pair recovery
    behavior; E.8.4 does not weaken that transaction.

## Static ownership guards

23. E.8.4 builder/renderer contains no Method-B numerical dependency.
24. It does not call F.4/F.5/F.6.1/F.6.2 calculators, `bin_data`,
    `calculate_yield_data`, pion-template fill, `bg_fit`, or a tree traversal.
25. It does not normalize/scale/rebin authoritative histograms.
26. It does not calculate a `DeltaY` uncertainty.

Do not duplicate production calculations in the tests as the implementation.
Tests should assert the consumer contract around already-produced F.6.3 objects.

---

# 16. Required regression tests

Run, with the locally discovered interpreter:

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

If the local interpreter is not literally invoked as `python`, use the actual
interpreter discovered on the host and record it. Do not pretend PyROOT is
available if it is not.

Do not repair unrelated pre-existing failures. Establish whether a failure is
introduced by this change; stop/report if the distinction cannot be made
cleanly.

---

# 17. Memory/history update

After the implementation and focused local checks are complete, add:

```text
docs/memory/phases/e8-4-production-impact-audit.md
```

Record:

- starting HEAD `d84f9610427cf38ccaaa8d6d8b14a8a995ecbac9`;
- exact changed source/test paths;
- exact runtime source -> consumer -> renderer path;
- E.8.4's presentation-only ownership;
- unavailable/fail-closed behavior;
- local deterministic test results;
- explicit statement that these tests are not farm/ROOT/full-runtime
  validation.

For the local implementation candidate, keep E.8.4 status **ACTIVE** pending
independent ChatGPT actual-diff/source-runtime-path review. Do not self-promote
it to `SOURCE REVIEWED`, `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`, or
`CLOSED / RUNTIME VALIDATED`.

Update, as warranted:

```text
docs/memory/CURRENT.md
docs/memory/decisions/e8-full-analysis-procedure-roadmap.md
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/roadmap/STATUS.md
```

The exact post-Codex `NEXT` must be:

```text
independent ChatGPT review of the complete E.8.4 actual diff/source/runtime path
```

Do not request a farm run until that source review passes.

Do not modify `MEMORY.md`, `USER.md`, `CURRENT_HANDOFF.md`, `LEARNINGS.md`,
`TOOLS.md`, or `CODEX.md` unless a genuinely new fact owned by one of those
records emerges. If none does, leave them unchanged.

Regenerate the versioned memory manifest:

```bash
python -B tools/update_memory_manifest.py --root . --write
python -B tools/update_memory_manifest.py --root . --check
python -B tools/check_memory_health.py --root .
python -B -m unittest testing.test_memory_health -v
```

---

# 18. Diff audit

Before stopping, inspect the actual diff, not a summary:

```bash
git status --short
git diff --stat
git diff -- src/cuts/full_background_subtraction_plots.py
git diff -- testing/test_e8_4_production_impact_audit.py
git diff -- docs/memory
git -c core.safecrlf=false diff --check
```

Confirm that none of the frozen producer/runtime files changed.

If the new test/phase files are still untracked, include them in review using
`git diff --no-index /dev/null <path>`; ordinary `git diff` alone is not enough.

If the complete review is too large for terminal copy/paste, create one
temporary review bundle in the repository root containing:

- all tracked diffs; and
- complete `git diff --no-index /dev/null ...` output for every intended
  untracked file.

Use a clearly temporary filename such as:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

Do not stage merely to make the diff reviewable. The review bundle itself must
remain untracked and must be removed before the eventual user-controlled
commit.

---

# 19. Acceptance criteria

The local E.8.4 candidate is implementation-complete only if all of the
following are true:

1. starting branch/HEAD contract was satisfied;
2. only allowlisted source/test/memory paths changed;
3. baseline production and F.6.3 producer source are byte-unchanged;
4. E.8.4 consumes `_f6_3_parallel_method_a_source` only after baseline yields;
5. no second tree traversal, Method-A reconstruction, fit, factor calculation,
   pion-template refill, or yield integration was added;
6. exact setting/window/geometry/baseline identity is fail-closed;
7. Method B is numerically absent;
8. empirical residual Fit 1/Fit 2 remain dormant;
9. no child normalization, clipping, smoothing, interpolation, or rebinning is
   introduced;
10. source histograms/public outputs are unmodified;
11. `DeltaY` and defined `DeltaY/Y0` come only from stored `Y0`/`YA` scalars;
12. no uncertainty is invented for the Method-A shift;
13. actual pion-background, clean-kaon/final-MM, signed-difference, and yield
    impacts are all visible in the existing procedure PDF route;
14. E.8.4 unavailable state is explicit and does not alter public baseline
    success;
15. existing pair-safe PDF/page-manifest replacement/recovery remains intact;
16. focused and required regression checks pass, apart from clearly established
    unchanged/pre-existing out-of-scope failures;
17. memory/status/manifest accurately describe **ACTIVE pending independent
    source review**, not runtime acceptance;
18. `git diff --check` passes;
19. Codex has not committed, pushed, or run the farm.

---

# 20. Hard stop

Stop after:

- implementation;
- deterministic local checks;
- warranted memory/history updates;
- manifest regeneration/check;
- actual diff audit; and
- preparation of a complete review bundle if needed.

Do **not**:

- commit;
- push;
- run the Jefferson Lab farm;
- claim ROOT/PyROOT/full-analysis/procedure-PDF runtime validation;
- begin final E.8 closure;
- begin F.6.4;
- promote Method A;
- broaden to the all-five-setting production-impact summary.

The next human/ChatGPT gate is independent review of the actual E.8.4 diff.
