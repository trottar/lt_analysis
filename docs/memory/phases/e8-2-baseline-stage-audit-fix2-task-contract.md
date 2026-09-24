# KaonLT E.8.2.Fix.2 — Production-pruning boundary and residual contract closure

## Objective

Repair the remaining E.8.2 source-review issues found during independent review
of the complete cumulative Fix.1 staging.

This is a **narrow E.8.2 source-review repair**, not a redesign of the KaonLT
analysis and not a new scientific subtraction.

The newly established source fact is that the authoritative yield path applies
the existing production `prune_hist(...)` treatment to the canonical wide
missing-mass histograms **after random/dummy subtraction and before pion
subtraction**. Because slow-proton cleaning is already carried by the event
factor during the authoritative fill, the production wide histogram immediately
before this pruning treatment is the true post-proton/pre-prune state, while the
histogram entering the accepted pion subtraction is the post-prune production
state.

E.8.2 must expose that existing boundary faithfully. It must not hide the
`prune_hist` treatment inside the slow-proton subtraction, and it must not rerun,
change, bypass, tune, or reinterpret `prune_hist`.

The repaired presentation chain is therefore:

```text
prompt pre-proton
    -> random subtraction
    -> dummy subtraction
    -> slow-proton event-factor removal
    -> authoritative post-proton / pre-prune production snapshot
    -> existing production prune_hist treatment
    -> authoritative post-prune pion input
    -> exact baseline pion background B_pi^0 using w0
    -> final baseline clean kaon MM_0(t,phi)
    -> existing Y_0(t,phi)
```

This repair also closes the remaining deterministic negative/regression-test
requirements and makes the valid-child E.8.2 stage inventory fail closed to
forbidden historical/Method-A/Method-B stage injection.

No production scientific calculation may change.

---

## Starting repository identity

Required branch:

`test`

Required repository HEAD:

`5b2ba7f92dc51b303257dbf209df687415364602`

Commit message:

`Correct E8 no-empirical-residual chain`

This repair starts from the **existing local cumulative E.8.2 + Fix.1
worktree**. Do not recreate E.8.2 from the remote checkout.

Before editing, run:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
git diff --name-only
git diff --check
```

Required conditions:

- branch is exactly `test`;
- HEAD is exactly
  `5b2ba7f92dc51b303257dbf209df687415364602`;
- the dirty worktree contains only the intended cumulative E.8.2/Fix.1 work and
  temporary review artifact(s);
- there are no unrelated user changes.

Do **not** reset, stash, clean, discard, checkout over, or recreate the existing
E.8.2/Fix.1 implementation.

Hard stop on a different branch, different HEAD, or unrelated worktree change.

Do not pull, commit, push, switch branches, update remote refs, or run the
Jefferson Lab farm.

---

## Required startup reading

Read repository memory in the normal required order:

1. `docs/memory/AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read:

- `docs/memory/CODEX.md`
- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/decisions/e8-no-empirical-residual-chain-correction.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/e8-2-baseline-stage-audit-task-contract.md`
- `docs/memory/phases/e8-2-baseline-stage-audit-fix1-task-contract.md`
- `docs/memory/phases/e8-2-baseline-stage-audit.md`
- this Fix.2 contract.

Then inspect the complete current cumulative diff and surrounding live source.

The original E.8.2 and Fix.1 contracts remain authoritative except where this
Fix.2 contract explicitly corrects the newly discovered production-pruning
boundary.

---

## Independently established source fact

At starting HEAD `5b2ba7f...`, the normal `process_hist_data(...)` path performs:

```text
event filling with accepted proton factor
    -> random normalization/subtraction
    -> data normalization
    -> dummy normalization/subtraction
    -> prune_hist(H_MM_nosub_DATA, event_threshold)
    -> baseline pion subtraction
```

The existing `prune_hist(...)` implementation mutates the histogram. It:

1. zeroes bins with unusually large fractional uncertainty; and
2. resets the histogram if entries are at/below threshold or the integral is
   non-positive.

Therefore the current Fix.1 E.8.2 relation:

```text
after_dummy_pre_proton
    - proton_component_removed
    = after_proton
```

cannot in general close if `after_proton` is defined as the already-pruned
production histogram while the left-hand diagnostic terms are built from the
unpruned same-traversal source contributions.

This is a presentation-authority mismatch only. It does **not** authorize any
change to `prune_hist` or production physics.

---

# Scientific ownership and frozen interfaces

## Frozen production behavior

Do not change:

- event selection or cuts;
- shifted missing mass;
- shifted t;
- canonical t/phi binning;
- random-window definition or coefficients;
- data normalization;
- dummy normalization;
- slow-proton PID model;
- slow-proton event-factor calculation;
- setting-wide proton application gate;
- accepted proton factor values;
- the location, arguments, thresholds, factor, or behavior of any production
  `prune_hist(...)` call;
- baseline pion-control population;
- baseline pion event weight `w0`;
- pion components;
- pion fits;
- pion fit/control windows;
- pion amplitudes/scales;
- pion application/fallback policy;
- public `stage_window_yields`;
- final baseline production yield formula;
- final baseline production yield error formula;
- SIMC;
- efficiencies;
- acceptance;
- L/T separation;
- cross sections;
- all accepted F.6.2 artifacts/fingerprints;
- all existing E.8.1/F.6.2 numerical payloads.

The accepted active profile remains:

`no_empirical_residual`

Historical empirical residual Fit 1/Fit 2 remain dormant and must not be
activated, called by E.8.2, displayed as E.8.2 stages, tuned, or reintroduced.

Method A remains numerically absent from E.8.2.

Method B remains numerically absent from E.8.2.

## Presentation ownership

E.8.2 may:

- clone an authoritative production histogram immediately before the existing
  production pruning call;
- clone the authoritative production histogram after that call;
- retain the already-built event-derived proton-removed diagnostic;
- validate identities/closures;
- compute producer-side diagnostic Lambda-window integrals from detached
  snapshots;
- render the existing production pruning as an explicit treatment boundary.

E.8.2 must **not** call `prune_hist` on a diagnostic clone to manufacture an
expected result.

E.8.2 must **not** infer the pre-prune state by undoing the post-prune state.

E.8.2 must **not** define pruning as a new physics correction or uncertainty.

---

# Allowed files

Repair source edits are limited to:

- `src/binning/calculate_yield.py`
- `src/cuts/full_background_subtraction_plots.py`

Focused test edits are limited to:

- `testing/test_e8_2_baseline_stage_audit.py`

Warranted memory edits are limited to:

- `docs/memory/phases/e8-2-baseline-stage-audit.md`
- new
  `docs/memory/phases/e8-2-baseline-stage-audit-fix2-task-contract.md`
- `docs/memory/manifest.json`

`docs/memory/CURRENT.md` must remain `ACTIVE` and should remain byte-unchanged
unless the new source fact makes one sentence factually incorrect. If that
occurs, make only the minimum wording correction while preserving `ACTIVE`,
E.8.3 `BLOCKED`, and F.6.3 `BLOCKED`.

Do not edit the original E.8.2 or Fix.1 contracts.

---

# Explicitly frozen files for Fix.2

Do not modify:

- `src/cuts/rand_sub.py`
- `src/main.py`
- `src/utility/utility.py`
- `src/utility/background_config.py`
- proton-cleaning scientific source;
- pion-subtraction scientific source outside the allowlisted producer;
- validation profiles;
- collectors;
- wrappers;
- launchers;
- farm scripts;
- accepted E.8.1/F.6.2 artifacts;
- existing regression-test modules outside
  `testing/test_e8_2_baseline_stage_audit.py`.

If another source file appears necessary, hard stop and report the blocker.

---

# Fix.2 producer contract

## A. Capture the actual post-proton / pre-prune production snapshot

Within the existing `process_hist_data(...)` canonical `(t,phi)` loop, capture a
detached clone of the **actual production**
`H_MM_nosub_DATA_{j}_{k}` at the exact point:

- after the established random subtraction;
- after data normalization;
- after established dummy normalization/subtraction;
- before the first production `prune_hist(...)` that mutates that histogram;
- before pion subtraction.

Recommended semantic name:

`after_proton_pre_prune`

or an equally explicit name.

This clone is observational only.

Do not move the production `prune_hist` call.

Do not add a second event loop.

Do not rerun proton cleaning.

Do not reconstruct this state from a later histogram.

## B. Retain the actual post-prune production state

The existing production histogram after its established `prune_hist(...)`
treatment and immediately before pion subtraction remains authoritative.

Recommended semantic name:

`after_proton_post_prune`

or:

`pion_input_production`

The exact accepted pion payload's:

`H_MM_nosub_before_pion_subtraction`

must be validated against this post-prune production state for consistent
binning and bin contents within the existing strict numerical tolerance.

Do not replace the exact accepted pion payload with a reconstructed clone.

## C. Slow-proton closure

The proton closure must be defined against the **pre-prune** post-proton
production snapshot:

```text
after_dummy_pre_proton
    - proton_component_removed
    = after_proton_pre_prune
```

Require bin-by-bin closure.

The right-hand side must be the actual detached production snapshot captured
before `prune_hist`, not the auxiliary post-proton calculation.

The existing auxiliary post-proton diagnostic may remain only as an internal
consistency check if useful, but it must not replace the authoritative production
snapshot in the E.8.2 payload.

## D. Existing production-pruning treatment

Represent the frozen production transition explicitly:

```text
after_proton_pre_prune
    -> existing production prune_hist treatment
    -> after_proton_post_prune
```

Do not require equality across this transition.

Do not manufacture a fitted or modeled pruning component.

If a diagnostic signed difference is retained producer-side, it must be labeled
only as the observed effect:

```text
prune_effect = after_proton_pre_prune - after_proton_post_prune
```

and marked presentation-only/non-authoritative. It is optional; a before/after
treatment presentation is sufficient.

## E. Pion handoff identity

Require the actual production post-prune state to match the exact accepted pion
input:

```text
after_proton_post_prune == pion_input
```

bin-for-bin within deterministic numerical tolerance.

If this identity fails for an otherwise accepted/current-chain child, fail the
entire E.8.2 source closed with an explicit reason.

Then retain the existing exact pion closure:

```text
pion_input
    - pion_component_removed
    = after_pion_final
```

No pion reconstruction is allowed.

---

# Revised E.8.2 stage inventory

For every valid child, the authoritative E.8.2 stage record must now distinguish
the two post-proton production states.

At minimum retain:

1. `prompt_pre_proton`
2. `random_component_pre_proton`
3. `after_random_pre_proton`
4. `dummy_component_pre_proton`
5. `after_dummy_pre_proton`
6. `proton_component_removed`
7. `after_proton_pre_prune`
8. `after_proton_post_prune`
9. `pion_input`
10. `pion_component_removed`
11. `after_pion_final`

`after_proton_post_prune` and `pion_input` are distinct authoritative snapshots
of the same scientific handoff and must validate equal; retaining both makes
their ownership explicit.

Do not include any Fit 1/Fit 2, Method-A, or Method-B numerical stage.

---

# Revised diagnostic stage-window sequence

Do not modify legacy public `stage_window_yields`.

The private E.8.2 producer-side diagnostic Lambda-window progression must expose
the pruning transition rather than hide it.

Use this ordered sequence:

1. `prompt_pre_proton`
2. `after_random_pre_proton`
3. `after_dummy_pre_proton`
4. `after_proton_pre_prune`
5. `after_proton_post_prune`
6. `after_pion_final`

Continue to show authoritative final `Y_0(t,phi)` separately.

These are diagnostic stage-window integrals, not replacements for final yield.

The renderer must not calculate them.

---

# Revised page contract

Retain the same total E.8.2 placement:

- after the existing E.8.1/E.8 parent pages;
- before the terminal handoff;
- handoff remains last.

Preserve the existing E.8.2 page count if practical.

## Random page

Unchanged.

## Dummy page

Unchanged.

## Slow-proton / production-handoff page

Update the existing proton page so that it does not silently equate a pruned
histogram with the unpruned proton closure.

It must show, for every canonical phi child:

1. `after_dummy_pre_proton`
2. `proton_component_removed`
3. `after_proton_pre_prune`
4. `after_proton_post_prune`

The first three columns own the exact slow-proton closure.

The fourth column is explicitly labeled as the result after the **existing
production pruning treatment** and as the input state passed to the pion stage.

A clear page note/header must state that production `prune_hist` occurs between
columns 3 and 4.

Do not label pruning as proton subtraction.

## Baseline pion page

Continue to show:

1. exact post-prune `pion_input`
2. exact applied baseline `B_pi^0`
3. exact `after_pion_final = K_0`

No other pion behavior changes.

## Final MM page

Unchanged except for any provenance text required to identify the post-prune
pion input correctly.

## Stage-yield page

Show the revised six producer-computed diagnostic stage integrals, with
pre-prune and post-prune proton states separately labeled, plus final
authoritative `Y_0` separately.

---

# Fail-closed authority rules

The Fix.1 fail-closed authority behavior remains in force.

Additionally fail the entire E.8.2 source/payload closed for an otherwise
accepted/current-chain child if:

- actual pre-prune production snapshot is missing;
- actual post-prune production snapshot is missing;
- either snapshot cannot be detached;
- their MM binning is inconsistent;
- post-prune production state and exact pion input have inconsistent binning;
- post-prune production state and exact pion input differ in bin content beyond
  the strict deterministic tolerance;
- the pre-prune proton closure fails;
- any required valid-child stage is absent.

Legitimate authoritative `zero`, `skip_bin`, rejected, or otherwise defined
non-applicable children remain explicit child-level unavailable states.

Do not invent substitute data.

---

# Exact valid-child stage schema

For a valid child, make the accepted E.8.2 stage-key inventory explicit.

The presentation builder must reject:

- missing required stage keys; and
- forbidden/unknown numerical stage keys that would silently introduce a
  scientific stage outside the E.8.2 contract.

At minimum explicitly reject attempted valid-child inclusion of:

- Fit 1;
- Fit 2;
- empirical residual stage aliases;
- Method A numerical stages/factors/templates;
- Method B numerical stages/factors/templates.

Do not reject private metadata that is stored outside the `stages` mapping.

Prefer exact stage-key equality for the valid-child `stages` mapping if that can
be implemented without hiding required provenance.

This is a schema/presentation guard only. Do not delete or edit dormant
historical production source.

---

# Required deterministic tests

Extend only:

`testing/test_e8_2_baseline_stage_audit.py`

Retain all existing Fix.1 tests.

## A. Production-pruning boundary

Add a deterministic producer/path test that proves the exact ordering in the
edited source path:

```text
dummy subtraction
-> E.8.2 actual pre-prune production clone
-> existing production prune_hist
-> E.8.2 post-prune production capture/use
-> pion subtraction
```

The test may use narrow mocks/fake histograms, but it must execute the relevant
producer helper/path rather than merely search text when practical.

Use a fixture where the pruning treatment changes the production histogram.

Verify:

- pre-prune production snapshot retains the original content;
- post-prune production snapshot reflects the existing production mutation;
- E.8.2 does not call `prune_hist` itself to manufacture the transition;
- production `prune_hist` call count and arguments are unchanged from the
  existing path.

## B. Proton closure versus prune transition

Verify:

```text
after_dummy_pre_proton
- proton_component_removed
= after_proton_pre_prune
```

for a nontrivial proton factor.

Then verify that `after_proton_post_prune` is allowed to differ when the frozen
production pruning treatment changes bins.

The slow-proton closure checker must not compare against the post-prune state.

## C. Pion handoff

Verify:

- exact post-prune production state equals exact accepted `pion_input`;
- mismatch in binning fails closed;
- mismatch in bin content fails closed;
- exact `pion_input - B_pi^0 = after_pion_final` closure remains unchanged.

## D. Stage-integral sequence

Verify the producer-side private sequence is exactly:

```text
prompt_pre_proton
after_random_pre_proton
after_dummy_pre_proton
after_proton_pre_prune
after_proton_post_prune
after_pion_final
```

Verify the renderer does not recompute these integrals.

## E. Forbidden-stage negative tests

Inject otherwise-valid source children containing:

- a Fit-1/Fit-2 stage;
- a Method-A numerical stage;
- a Method-B numerical stage.

Each must make the E.8.2 presentation unavailable/fail closed.

## F. Original negative-contract coverage still missing

Add focused deterministic checks for:

- malformed canonical t/phi geometry;
- inconsistent MM binning;
- malformed E.8.2 source schema;
- missing Step-3 render state at finalization.

## G. Public production-regression coverage

Add the narrowest practical deterministic fixture proving E.8.2 observation
does not change:

- existing public `stage_window_yields`;
- existing returned public `Y_0`;
- existing returned public yield error.

Also prove the private E.8.2 statistical error is sourced from the existing
`integral_with_stat_error(final_hist)` result rather than independently
recomputed by the renderer.

Do not create a duplicate fake analysis.

If directly calling `calculate_yield_data(...)` requires excessive unrelated
stubbing, extract only the smallest pure assembly helper needed to test the
unchanged public-return/private-sidecar split. Production scientific functions
remain authoritative.

## H. `bg_fit` non-use

Strengthen the existing spy so the tested E.8.2 producer/builder/render helper
path cannot invoke `bg_fit`.

Do not alter dormant historical Fit 1/Fit 2 source to satisfy the test.

## I. Finalizer tests

Retain all Fix.1 pair-safe finalization tests unchanged unless a genuine bug is
found.

---

# Existing regression suites

Run unchanged:

```bash
python -B -m unittest testing.test_full_background_subtraction_plots -v
python -B -m unittest testing.test_proton_cleaning_runtime_status -v
python -B -m unittest testing.test_pion_t_bin_particle_stage -v
python -B -m unittest testing.test_pion_hgcer_phase_f_runtime_contract -v
```

Do not edit those test modules merely to make E.8.2 pass.

The already-established
`testing.test_proton_cleaning_runtime_status` 34-pass/one-error result is
pre-existing at exact starting HEAD and outside E.8.2 ownership. Reconfirm it
did not change; do not repair proton source here.

If that regression differs from the exact starting-HEAD behavior, hard stop.

---

# Local deterministic validation

At minimum run:

```bash
python -B -m py_compile \
  src/binning/calculate_yield.py \
  src/cuts/full_background_subtraction_plots.py \
  testing/test_e8_2_baseline_stage_audit.py

python -B -m unittest testing.test_e8_2_baseline_stage_audit -v
python -B -m unittest testing.test_full_background_subtraction_plots -v
python -B -m unittest testing.test_proton_cleaning_runtime_status -v
python -B -m unittest testing.test_pion_t_bin_particle_stage -v
python -B -m unittest testing.test_pion_hgcer_phase_f_runtime_contract -v
```

Also run:

- established memory-manifest regeneration;
- manifest `--check`;
- memory health/integrity;
- memory bootstrap check;
- `git diff --check`.

Report the exact interpreter used if it differs from `python`.

Local checks do not establish ROOT/PyROOT, full `main.py`, actual procedure-PDF
rendering on the farm, or runtime validation.

---

# Memory update

Update `docs/memory/phases/e8-2-baseline-stage-audit.md` narrowly.

Record:

- Fix.2;
- the newly confirmed production `prune_hist` boundary;
- the authoritative pre-prune and post-prune snapshot ownership;
- slow-proton closure against the pre-prune production snapshot;
- explicit production-pruning treatment before pion input;
- post-prune/pion-input identity;
- exact valid-child stage inventory guard;
- expanded deterministic tests actually run;
- unchanged production/scientific behavior;
- unchanged pre-existing proton regression result;
- farm boundary.

Correct any earlier phase-record wording that says the old three-column proton
page directly closed onto the already-pruned production histogram.

Do not rewrite historical runtime evidence.

E.8.2 remains:

`ACTIVE`

pending independent ChatGPT review of the complete repaired cumulative diff.

Do not mark E.8.2 `SOURCE REVIEWED`.

E.8.3 remains:

`BLOCKED`

F.6.3 remains:

`BLOCKED`

Regenerate `docs/memory/manifest.json`.

---

# Diff audit

Before finishing, run:

```bash
git status --short
git diff --stat
git diff --check
git diff --name-only
```

Audit that Fix.2 changed only its allowlisted repair files beyond the prior
E.8.2/Fix.1 staging.

Explicitly confirm no new changes to:

- `src/cuts/rand_sub.py`;
- `src/main.py`;
- `src/utility/utility.py`;
- `src/utility/background_config.py`;
- proton-cleaning science;
- pion science outside `calculate_yield.py`;
- random/dummy coefficients;
- proton-factor calculation;
- production `prune_hist` calls/arguments/thresholds;
- `w0`;
- pion components/fits/windows;
- canonical binning;
- yield formulas;
- public stage-yield schema;
- Method A;
- Method B numerical use;
- dormant empirical residual configuration.

---

# Required complete review bundle

Regenerate the root temporary:

`kaonlt_review.diff`

It must be the **complete cumulative E.8.2 + Fix.1 + Fix.2 diff from starting
HEAD `5b2ba7f...`**, not an incremental Fix.2-only diff.

Include:

1. complete tracked `git diff`;
2. complete `git diff --no-index -- /dev/null ...` additions for every intended
   untracked file.

At minimum include the complete additions for:

- `docs/memory/phases/e8-2-baseline-stage-audit-task-contract.md`
- `docs/memory/phases/e8-2-baseline-stage-audit-fix1-task-contract.md`
- `docs/memory/phases/e8-2-baseline-stage-audit-fix2-task-contract.md`
- `docs/memory/phases/e8-2-baseline-stage-audit.md`
- `testing/test_e8_2_baseline_stage_audit.py`

Do not add `kaonlt_review.diff` to Git.

---

# Acceptance criteria

Fix.2 is ready for independent review only if:

1. branch remains `test`;
2. HEAD remains exactly `5b2ba7f92dc51b303257dbf209df687415364602`;
3. only allowlisted Fix.2 paths changed beyond the cumulative staging;
4. no production scientific behavior changed;
5. no second scientific tree traversal was added;
6. production `prune_hist` behavior/order/arguments are unchanged;
7. actual production post-proton/pre-prune state is captured before pruning;
8. actual production post-prune state is retained separately;
9. slow-proton closure is against the authoritative pre-prune production state;
10. the pruning transition is explicit and is not mislabeled as proton
    subtraction;
11. post-prune production state equals exact accepted pion input or E.8.2 fails
    closed;
12. exact baseline pion before/component/after authority remains unchanged;
13. revised diagnostic stage integrals expose both pre-prune and post-prune
    states;
14. public `stage_window_yields` is unchanged;
15. public `Y_0` is unchanged;
16. public yield error is unchanged;
17. private statistical error remains separately sourced from existing yield
    integration;
18. valid-child stage inventory fails closed to Fit 1/Fit 2 and Method-A/Method-B
    numerical injection;
19. malformed geometry/binning/schema and missing Step-3 state are tested
    fail-closed;
20. Fix.1 renderer-failure and pair-recovery behavior remains intact;
21. existing E.8.1/F.6.2 payloads remain unchanged;
22. `no_empirical_residual` remains literal;
23. Method A remains numerically absent;
24. Method B remains numerically absent;
25. required deterministic and memory checks are run and reported;
26. the pre-existing proton-regression result is unchanged and remains
    out-of-scope;
27. memory remains `ACTIVE` pending ChatGPT review;
28. no farm/runtime claim is made;
29. a complete cumulative review bundle is generated.

---

# Hard stop

Stop and report rather than broadening scope if:

- the authoritative pre-prune production snapshot cannot be captured without
  moving or changing production pruning;
- `prune_hist` would need to be changed or rerun diagnostically;
- repair requires changing proton physics;
- repair requires changing pion physics;
- repair requires changing `w0`;
- repair requires changing random/dummy normalization;
- repair requires changing canonical binning;
- repair requires changing yield formulas;
- repair requires introducing Method A or Method B numerically;
- repair requires activating empirical Fit 1/Fit 2;
- repair requires another source file;
- the existing proton regression changes relative to exact starting HEAD;
- unrelated worktree changes are discovered.

Codex must not commit, push, update remote refs, or run the Jefferson Lab farm.
