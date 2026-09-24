# KaonLT E.8.2.Fix.3 — Contract-closure test and memory reconciliation

## Objective

Close the remaining E.8.2 source-review gap after Fix.2.

Independent review of the complete cumulative E.8.2 + Fix.1 + Fix.2 diff finds
the scientific/runtime-path implementation itself consistent with the accepted
E.8.2 architecture. Fix.3 does **not** redesign or change the E.8.2 scientific
path.

The remaining blocker is narrower:

1. the Fix.2 contract explicitly required deterministic regression coverage
   proving that E.8.2 observation does not change the existing public
   `stage_window_yields`, returned `Y_0`, or returned yield error;
2. the current focused test module does not exercise or assert those public
   outputs;
3. the existing `bg_fit` spy only wraps the presentation-reader call and should
   be strengthened so the E.8.2 producer/builder path is explicitly covered;
4. the E.8.2 phase record still describes the PDF/manifest installation as
   "atomic", while Fix.1 correctly implemented a pair-safe two-file transaction
   with recovery rather than impossible literal multi-file filesystem
   atomicity.

This task is therefore a **deterministic contract-closure and memory wording
repair only**.

No production scientific behavior may change.

---

## Starting repository identity

Required branch:

`test`

Required repository HEAD:

`5b2ba7f92dc51b303257dbf209df687415364602`

Commit message:

`Correct E8 no-empirical-residual chain`

This task starts from the existing local cumulative E.8.2 + Fix.1 + Fix.2
worktree.

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
- the worktree contains only the intended cumulative E.8.2/Fix.1/Fix.2 changes
  and temporary review artifact(s);
- there are no unrelated user changes.

Do not reset, stash, clean, discard, checkout over, or recreate the existing
worktree.

Do not commit, push, switch branches, update remote refs, or run the Jefferson
Lab farm.

Hard stop on a different branch/HEAD or unrelated worktree changes.

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
- `docs/memory/phases/e8-2-baseline-stage-audit-fix2-task-contract.md`
- `docs/memory/phases/e8-2-baseline-stage-audit.md`
- this Fix.3 contract.

Then inspect the complete cumulative diff before editing.

All earlier E.8.2/Fix.1/Fix.2 contracts remain authoritative.

---

# Review conclusion carried into Fix.3

The current cumulative source implementation is accepted **for the purposes of
this repair contract** as the architecture to preserve:

```text
same authoritative event traversal
    -> true pre-proton auxiliary capture
    -> random/dummy normalization and subtraction
    -> actual post-proton/pre-prune production snapshot
    -> unchanged production prune_hist treatment
    -> actual post-prune production snapshot
    -> exact accepted pion input / B_pi^0 / K_0
    -> existing Y_0
```

The Fix.1 pair-safe procedure-PDF/manifest finalization behavior is also frozen.

Fix.3 must not reopen these implementations absent a deterministic test that
reveals an actual defect.

---

# Scientific ownership and frozen behavior

Do not change:

- event selection/cuts;
- shifted MM or t;
- canonical t/phi binning;
- event traversal;
- random windows or coefficients;
- data/dummy normalization;
- proton PID or proton factors;
- production `prune_hist` calls, ordering, arguments, thresholds, or behavior;
- baseline pion control;
- pion `w0`;
- pion fits/components/windows/scales/application policy;
- public `stage_window_yields` implementation;
- public yield formula;
- public yield-error formula;
- SIMC;
- efficiencies;
- acceptance;
- L/T separation;
- cross sections;
- E.8.1/F.6.2 frozen payloads/fingerprints;
- Fix.1 finalizer semantics;
- E.8.2 pre-prune/post-prune stage ownership;
- E.8.2 page ordering.

Accepted active profile remains:

`no_empirical_residual`

Fit 1/Fit 2 remain dormant and excluded from E.8.2.

Method A remains numerically absent.

Method B remains numerically absent.

---

# Allowed files

Expected edits are limited to:

- `testing/test_e8_2_baseline_stage_audit.py`
- `docs/memory/phases/e8-2-baseline-stage-audit.md`
- new
  `docs/memory/phases/e8-2-baseline-stage-audit-fix3-task-contract.md`
- `docs/memory/manifest.json`

`docs/memory/CURRENT.md` should remain byte-unchanged and status must remain
`ACTIVE` pending independent ChatGPT review.

## Source files are frozen by default

Do not modify:

- `src/binning/calculate_yield.py`
- `src/cuts/full_background_subtraction_plots.py`
- `src/cuts/rand_sub.py`
- `src/main.py`
- any other scientific/runtime source.

If a deterministic public-output fixture cannot be written without a tiny
testability extraction from `src/binning/calculate_yield.py`, hard stop and
report the exact blocker rather than editing source automatically. ChatGPT will
decide whether a separate source-allowlisted repair is warranted.

---

# Required deterministic coverage

## A. Public yield/output non-regression

Add a focused deterministic test that exercises the existing
`calculate_yield_data(...)` public-return path with the narrowest practical fake
histogram/processed-data fixture.

The test must establish that adding the private E.8.2 sidecar does not change
the public production outputs.

At minimum verify:

1. the returned public `Y_0` value is exactly the value produced by the existing
   yield integration path;
2. the returned public yield error is exactly the existing production error;
3. the existing public `stage_window_yields` object/value remains unchanged by
   E.8.2 sidecar assembly;
4. the private E.8.2 record copies the same `Y_0`;
5. the private `statistical_error` is the result from the existing
   `integral_with_stat_error(final_hist)` call;
6. the private `total_error` is the same existing public production yield error;
7. building the private E.8.2 source does not mutate the public returned yield
   container or the existing processed-entry `stage_window_yields`.

Prefer mocks/stubs around unrelated I/O and plotting dependencies.

Do not reproduce the scientific yield calculation in test code.

Do not merely assert source-text strings for these public-output invariants:
the fixture must execute the relevant production function/path.

## B. Public `stage_window_yields` preservation

The fixture must contain a sentinel/known `stage_window_yields` mapping on the
processed entry and verify it is byte/value-equivalent after E.8.2 source
assembly.

E.8.2-specific `stage_window_integrals` must remain separate.

## C. Strengthened dormant-`bg_fit` exclusion

Strengthen the existing test so a `bg_fit` spy covers the executable E.8.2
producer-side path as well as the presentation builder.

At minimum exercise, as applicable:

- `_build_e8_2_early_stage_capture`
- `_build_e8_2_pion_stage_capture`
- `_build_e8_2_baseline_stage_source`
- `build_full_background_subtraction_e8_2_payload`

with a spy/mock that raises if `bg_fit` is invoked.

The test must prove the E.8.2 helpers under test do not invoke the dormant fit.

Do not alter the historical Fit 1/Fit 2 source.

## D. Preserve all existing focused tests

All current Fix.1/Fix.2 tests remain and must continue to pass, including:

- one-pass proton capture;
- random/dummy/proton closure;
- pre-prune/post-prune boundary;
- post-prune/pion-input identity;
- exact pion authority;
- yield-authority failures;
- forbidden stage injection;
- malformed geometry/binning/schema;
- missing Step-3 state;
- render-state detachment;
- page ordering;
- renderer-failure behavior;
- pair-safe installation/recovery.

---

# Memory wording repair

Update only `docs/memory/phases/e8-2-baseline-stage-audit.md`.

Correct the earlier sentence that says the preliminary PDF/manifest pair is
"atomically replaced".

The durable wording must match the actual Fix.1 implementation:

- final temporary PDF and manifest are completed first;
- recoverable copies of the preliminary pair are created;
- the two final names are installed sequentially;
- handled partial installation failure restores the preliminary pair;
- preservation is claimed only when recovery succeeds;
- transient render state is released only after successful installation of both.

Do not claim literal multi-file filesystem atomicity.

Add Fix.3 chronology recording:

- the missing public-output regression coverage;
- the strengthened `bg_fit` non-use coverage;
- exact tests actually run;
- no source/scientific change;
- farm boundary.

E.8.2 remains:

`ACTIVE`

pending independent ChatGPT actual-diff review.

E.8.3 remains:

`BLOCKED`

F.6.3 remains:

`BLOCKED`

Do not alter historical runtime evidence.

Regenerate `docs/memory/manifest.json`.

---

# Existing regression suites

Run unchanged:

```bash
python -B -m unittest testing.test_full_background_subtraction_plots -v
python -B -m unittest testing.test_proton_cleaning_runtime_status -v
python -B -m unittest testing.test_pion_t_bin_particle_stage -v
python -B -m unittest testing.test_pion_hgcer_phase_f_runtime_contract -v
```

Do not edit those suites.

The known `testing.test_proton_cleaning_runtime_status` 34-pass/one-error result
is pre-existing at exact starting HEAD and outside this task. Reconfirm that the
same result persists. Do not repair or suppress it.

Hard stop if that result changes relative to the exact starting-HEAD behavior.

---

# Local deterministic validation

At minimum run:

```bash
python -B -m py_compile \
  testing/test_e8_2_baseline_stage_audit.py

python -B -m unittest testing.test_e8_2_baseline_stage_audit -v
python -B -m unittest testing.test_full_background_subtraction_plots -v
python -B -m unittest testing.test_proton_cleaning_runtime_status -v
python -B -m unittest testing.test_pion_t_bin_particle_stage -v
python -B -m unittest testing.test_pion_hgcer_phase_f_runtime_contract -v
```

Also run:

- memory-manifest regeneration;
- manifest `--check`;
- memory health/integrity;
- memory bootstrap check;
- `git diff --check`.

Local checks do not establish ROOT/PyROOT, full `main.py`, procedure-PDF farm
rendering, or runtime validation.

---

# Diff audit

Before finishing, run:

```bash
git status --short
git diff --stat
git diff --check
git diff --name-only
```

Confirm Fix.3 itself changes no production/scientific source.

Explicitly confirm no Fix.3 changes to:

- `src/binning/calculate_yield.py`;
- `src/cuts/full_background_subtraction_plots.py`;
- `src/cuts/rand_sub.py`;
- `src/main.py`;
- utility/background configuration;
- proton science;
- pion science;
- `prune_hist`;
- random/dummy normalization;
- yield formulas;
- public output schema;
- Method A;
- Method B;
- empirical residual configuration.

---

# Required complete review bundle

Regenerate root temporary:

`kaonlt_review.diff`

It must contain the **complete cumulative E.8.2 + Fix.1 + Fix.2 + Fix.3 diff
from starting HEAD `5b2ba7f...`**, not an incremental Fix.3 delta.

Include:

1. complete tracked `git diff`;
2. complete `git diff --no-index -- /dev/null ...` for every intended untracked
   file.

At minimum include complete additions for:

- `docs/memory/phases/e8-2-baseline-stage-audit-task-contract.md`
- `docs/memory/phases/e8-2-baseline-stage-audit-fix1-task-contract.md`
- `docs/memory/phases/e8-2-baseline-stage-audit-fix2-task-contract.md`
- `docs/memory/phases/e8-2-baseline-stage-audit-fix3-task-contract.md`
- `docs/memory/phases/e8-2-baseline-stage-audit.md`
- `testing/test_e8_2_baseline_stage_audit.py`

Do not add `kaonlt_review.diff` to Git.

---

# Acceptance criteria

Fix.3 is ready for independent review only if:

1. branch remains `test`;
2. HEAD remains exactly `5b2ba7f92dc51b303257dbf209df687415364602`;
3. no production/scientific source changed in Fix.3;
4. all prior E.8.2/Fix.1/Fix.2 scientific behavior is preserved;
5. a deterministic executable fixture proves public `Y_0` is unchanged;
6. a deterministic executable fixture proves public yield error is unchanged;
7. a deterministic executable fixture proves public `stage_window_yields` is
   unchanged;
8. private statistical error is demonstrably sourced from the existing
   `integral_with_stat_error(final_hist)` result;
9. private total error matches the existing public error;
10. private sidecar construction does not mutate public containers;
11. strengthened executable coverage proves the E.8.2 producer/builder helpers
    do not invoke `bg_fit`;
12. all existing Fix.1/Fix.2 focused tests remain green;
13. existing regression suites retain their prior outcomes;
14. the known proton regression remains identical and out-of-scope;
15. phase memory no longer claims literal multi-file atomic replacement;
16. memory remains `ACTIVE` pending ChatGPT review;
17. no farm/runtime claim is made;
18. manifest/memory checks and `git diff --check` pass;
19. complete cumulative `kaonlt_review.diff` is regenerated.

---

# Hard stop

Stop and report rather than broadening scope if:

- the public-output regression fixture reveals any actual change to `Y_0`,
  yield error, or `stage_window_yields`;
- a source change appears necessary;
- `bg_fit` is actually reached from an E.8.2 helper;
- an existing Fix.1/Fix.2 focused test regresses;
- the proton regression differs from the exact starting-HEAD behavior;
- unrelated worktree changes are discovered.

Codex must not commit, push, update remote refs, or run the Jefferson Lab farm.
