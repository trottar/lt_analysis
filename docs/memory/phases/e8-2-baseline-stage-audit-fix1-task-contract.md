# KaonLT E.8.2.Fix.1 — Source-review repair task contract

## Objective

Repair only the narrow E.8.2 source-review defects identified during ChatGPT's
independent review of the complete `kaonlt_review.diff`.

This is **not** a redesign of E.8.2 and must not broaden into E.8.3, F.6.3,
Method-A production application, Method-B numerical use, proton-physics changes,
pion-model changes, or farm validation.

The existing E.8.2 scientific architecture is retained:

```text
prompt/random
    -> dummy
    -> accepted slow proton cleaning
    -> accepted baseline pion background using w0
    -> final baseline clean kaon MM_0(t,phi)
    -> existing Y_0(t,phi)
```

The repair owns only:

1. fail-closed finalization when any renderer failure is reported;
2. pair-safe PDF/page-manifest replacement with recovery on partial replacement
   failure;
3. fail-closed source authority for populated/current-chain children that lose
   required authoritative stage or yield inputs;
4. the deterministic E.8.2 producer/finalizer tests required by the original
   E.8.2 contract;
5. warranted E.8.2 phase-memory reconciliation and manifest regeneration.

No production scientific calculation may change.

---

## Starting repository identity

Required branch:

`test`

Required repository HEAD:

`5b2ba7f92dc51b303257dbf209df687415364602`

Commit message:

`Correct E8 no-empirical-residual chain`

This repair starts from the **existing local E.8.2 implementation worktree**,
not from a clean checkout.

The expected worktree already contains the intended E.8.2 changes represented
by the review bundle generated after the original E.8.2 implementation. Those
changes are local/proposed and must be preserved.

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
- the existing dirty worktree consists only of the intended E.8.2 implementation
  paths and temporary review artifact(s);
- there are no unrelated user changes.

Do **not** reset, stash, clean, checkout over, discard, or recreate the existing
E.8.2 implementation.

Hard stop if the branch/HEAD differs or unrelated worktree changes exist.

Do not pull, commit, push, switch branches, or run the Jefferson Lab farm.

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
- `docs/memory/phases/e8-2-baseline-stage-audit.md`
- this repair contract.

Then inspect the complete current local diff and surrounding source before
editing.

The original E.8.2 task contract remains authoritative for all requirements not
explicitly narrowed or clarified here.

---

## Independent review findings to repair

### Finding 1 — renderer failures currently still permit artifact replacement

Current local behavior collects:

```python
failures.extend(rendered.get("failures") or ())
```

but continues through temporary manifest creation and replacement even when
`failures` is non-empty.

That violates the original E.8.2 failure contract.

### Required behavior

A non-empty renderer-failure list is a finalization failure.

Before any final artifact replacement:

- close the temporary PDF cleanly;
- detect any non-empty renderer failure list;
- do not install the temporary PDF;
- do not install the temporary manifest;
- remove incomplete temporary replacement files;
- preserve the existing preliminary PDF and manifest;
- retain the transient Step-3 render state for diagnosis/retry;
- record explicit private finalization status;
- report `status = "unavailable"` or the existing equivalent fail-closed status;
- never return E.8.2 finalization as available with non-empty renderer failures.

Do not change the ordinary production-analysis failure policy.

---

### Finding 2 — sequential `os.replace` calls do not preserve the PDF/manifest pair

Current local behavior performs the equivalent of:

```python
os.replace(temporary_pdf, pdf_path)
os.replace(temporary_manifest, manifest_path)
```

If the first replacement succeeds and the second fails, a mixed new-PDF /
old-manifest pair can remain while status claims the preliminary pair was
preserved.

### Required behavior

Implement a **pair-safe transaction with recovery**.

The contract does not require impossible multi-file filesystem atomicity. It
requires the observable pair to be recoverable and not left mixed after a
handled failure.

Required semantics:

1. final temporary PDF exists and is closed;
2. final temporary manifest exists and is complete;
3. renderer failure list is empty;
4. preserve recoverable copies/backups of the existing preliminary PDF and
   manifest before installing either new final artifact;
5. install the new PDF and new manifest;
6. only after both installations succeed:
   - remove recovery copies;
   - update in-memory final manifest/failure state;
   - remove the transient Step-3 render state;
   - return available/success;
7. if either installation fails:
   - restore the original preliminary PDF and manifest pair;
   - remove temporary/new/recovery files as appropriate;
   - return unavailable/failure;
   - retain Step-3 render state;
   - set `preliminary_artifact_preserved = True` only if recovery actually
     succeeded;
   - if recovery itself fails, report that explicitly and do not falsely claim
     preservation.

Use ordinary Python/POSIX file operations only.

Do not use external PDF merge/concatenation tools.

Do not rerun any scientific calculation.

---

### Finding 3 — missing authority can be downgraded to an ordinary invalid child

The original E.8.2 contract distinguishes:

- legitimate authoritative child states that may be represented as unavailable,
  such as `zero`, `skip_bin`, rejected application, or an explicitly invalid
  child state;

from:

- missing required authority for a child that otherwise claims the current
  populated/accepted subtraction path.

The latter must fail E.8.2 closed.

Current local source can record
`final_yield_authority_missing` or an accepted pion application with missing
exact objects as an invalid child and still later mark the overall E.8.2 source
available.

### Required behavior

Preserve legitimate authoritative invalid/fallback states exactly and continue
to represent those children explicitly.

However, fail the **entire E.8.2 source/payload closed** when a current-chain
child that otherwise claims authoritative applicability loses required
authority.

At minimum this includes:

- true pre-proton/early-stage capture unexpectedly missing for a canonical
  current-chain child;
- accepted pion application missing any exact required authoritative
  before/template/after object;
- accepted pion application clone/detachment failure;
- required stage histogram missing after a child has entered the authoritative
  current-chain path;
- missing final `Y_0` authority for an otherwise valid populated child;
- missing statistical uncertainty for such a child;
- missing existing total uncertainty for such a child;
- any of those yield quantities non-finite.

Use explicit literal source-level failure reasons.

Do not manufacture substitute spectra.

Do not reconstruct `B_pi^0` from displayed input/output.

Do not silently convert an authority failure into a normal empty/invalid child.

It is acceptable and required to keep explicit child-level unavailable records
for genuine authoritative `zero`, `skip_bin`, rejected, or otherwise defined
non-applicable states.

---

## Scientific ownership and frozen behavior

The following remain frozen and must not change:

- event selection and cuts;
- shifted missing mass;
- shifted t;
- canonical t/phi binning;
- random-window definition;
- random subtraction coefficients;
- data normalization;
- dummy normalization;
- slow-proton PID model;
- slow-proton event-factor calculation;
- setting-wide proton application gate;
- production proton-cleaning factors;
- baseline pion-control definition;
- baseline pion event weight `w0`;
- pion components;
- pion fits;
- pion fit/control windows;
- pion amplitudes/scales;
- pion fallback/application policy;
- public `stage_window_yields`;
- final baseline production yield formula;
- final baseline production yield error formula;
- SIMC;
- efficiencies;
- acceptance;
- L/T separation;
- cross sections;
- F.6.2 frozen artifact and fingerprints;
- all existing E.8.1/F.6.2 numerical payloads.

The accepted profile remains:

`no_empirical_residual`

Historical empirical residual Fit 1/Fit 2 remain dormant and must not be
activated, called for E.8.2, presented, tuned, or reintroduced as stages.

Method A remains numerically absent from E.8.2.

Method B remains numerically absent from E.8.2.

---

## Allowed source files

Repair source edits are limited to:

- `src/binning/calculate_yield.py`
- `src/cuts/full_background_subtraction_plots.py`

Focused test edits are limited to:

- `testing/test_e8_2_baseline_stage_audit.py`

Warranted memory edits are limited to:

- `docs/memory/phases/e8-2-baseline-stage-audit.md`
- new
  `docs/memory/phases/e8-2-baseline-stage-audit-fix1-task-contract.md`
- `docs/memory/manifest.json`

`docs/memory/CURRENT.md` should remain unchanged unless the repair makes an
existing proposed E.8.2 sentence factually false. Its status must remain
`ACTIVE` and still await independent actual-diff review.

Do not edit the original E.8.2 task contract.

---

## Explicitly frozen files for this repair

Do not modify:

- `src/cuts/rand_sub.py`
- `src/main.py`
- `src/utility/background_config.py`
- proton-cleaning scientific source;
- pion-subtraction scientific source outside the allowlisted producer file;
- validation profiles;
- collectors;
- wrappers;
- launchers;
- farm scripts;
- E.8.1/F.6.2 accepted artifacts;
- `testing/test_proton_cleaning_runtime_status.py`
- `testing/test_pion_t_bin_particle_stage.py`
- `testing/test_pion_hgcer_phase_f_runtime_contract.py`

If another source file appears necessary, hard stop and report the blocker.

---

## Producer-level deterministic test requirements

The original focused test module was insufficient because it exercised mostly
presentation helpers.

Extend `testing/test_e8_2_baseline_stage_audit.py` so that it directly tests the
new producer behavior without creating a duplicate scientific analysis.

Use fake histograms, fake trees/events, and mocks/stubs where practical.

### A. Same-traversal proton capture

Test accepted factor `1.0`:

- pre-proton fill equals post-proton fill;
- proton-removed fill is zero;
- ordinary production fill remains the existing factor-weighted fill.

Test one nontrivial accepted factor, for example `0.4`:

- pre-proton contribution is `1.0`;
- production/post-proton contribution is `0.4`;
- removed contribution is `0.6`;
- ordinary production fill is still `0.4`;
- canonical t/phi assignment is identical;
- the input tree is traversed exactly once.

The test must detect accidental introduction of a second E.8.2 tree loop.

### B. Early-stage source-role algebra

Test exact deterministic closure for:

```text
prompt - random_component = after_random
after_random - dummy_component = after_dummy
after_dummy - proton_component_removed = after_proton
```

Verify the established normalization ownership:

- prompt uses `normfac_data`;
- random uses `normfac_data / nWindows`;
- dummy contribution is
  `(dummy_prompt - dummy_random / nWindows) * normfac_dummy`.

Do not change the production normalization formulas to satisfy the test.

### C. Baseline pion authority

Test:

- an accepted child uses the exact supplied production pion template;
- exact before - template = after closure;
- exact post-pion wide object is retained;
- genuine `zero`/`skip_bin`/rejected/unavailable status remains an explicit
  child-level state with its authoritative reason;
- accepted child with missing exact template/input/output causes whole-source
  fail-closed authority failure;
- accepted child clone failure causes fail-closed authority failure;
- no renderer reconstruction is introduced.

### D. Final yield authority

Test:

- existing final `Y_0` value is copied into the private E.8.2 record;
- existing histogram statistical error is stored separately;
- existing total error is stored separately;
- missing measurement for an otherwise valid child fails the whole source
  closed;
- malformed/non-finite yield/statistical/total values fail the whole source
  closed;
- legitimate authoritative invalid child does not require fabricated yield.

### E. Current-chain stage inventory

Retain and strengthen the existing checks that E.8.2 contains only:

- prompt;
- random;
- dummy;
- proton;
- baseline pion/final;

and contains no numerical:

- Fit 1;
- Fit 2;
- Method A;
- Method B.

Add a spy/mock or equivalent deterministic call-path test proving E.8.2
producer/builder/renderer helpers do not call `bg_fit`.

Do not achieve this by deleting or altering dormant historical production
source.

### F. Production regression

Within deterministic fixture coverage, verify E.8.2 observation does not alter:

- ordinary production post-proton histogram contents;
- proton factor;
- baseline pion template object/content;
- baseline post-pion object/content;
- existing public `stage_window_yields`;
- existing public `Y_0`;
- existing public yield error.

Use the narrowest practical fixtures. Do not reimplement the full analysis in
the test.

---

## Finalizer deterministic test requirements

Retain the existing render-state-detachment and ordering tests, and add:

### A. Renderer failures

Mock the complete renderer so it returns a non-empty failure list without
raising.

Verify:

- finalization returns unavailable/failure;
- neither final PDF nor final manifest is replaced;
- preliminary PDF content remains unchanged;
- preliminary manifest content remains unchanged;
- transient render state is retained;
- no success status is recorded.

### B. Replacement failure after PDF installation begins

Inject a failure on the manifest installation/replacement after the PDF
installation step has started.

Verify:

- finalization returns unavailable/failure;
- original preliminary PDF is restored;
- original preliminary manifest is restored;
- no mixed old/new pair remains;
- temporary/recovery artifacts are cleaned as appropriate;
- preservation status is true only when restoration succeeded;
- render state is retained.

### C. Recovery failure

Inject a restoration failure.

Verify the status does **not** claim
`preliminary_artifact_preserved = True`.

Record an explicit recovery failure reason.

### D. Successful transaction

Verify:

- both final PDF and final manifest are the new pair;
- page manifest in memory corresponds to the installed manifest;
- renderer failure list is empty;
- transient Step-3 render state is removed only after both installs succeed.

### E. Page order

Strengthen the existing ordering fixture so it actually includes at least one
existing E.8 parent page and proves:

```text
E.8 context
-> existing E.8 parent page(s)
-> E.8.2 block
-> terminal handoff
```

The handoff remains last.

---

## Existing regression suites

Run unchanged:

```bash
python -B -m unittest testing.test_full_background_subtraction_plots -v
python -B -m unittest testing.test_proton_cleaning_runtime_status -v
python -B -m unittest testing.test_pion_t_bin_particle_stage -v
python -B -m unittest testing.test_pion_hgcer_phase_f_runtime_contract -v
```

Do not edit those tests merely to make E.8.2 pass.

### Existing proton-cleaning regression error

The initial E.8.2 implementation report recorded one error in
`testing.test_proton_cleaning_runtime_status`, involving an unchanged
proton-cleaning fixture/path.

This repair does **not** own proton-cleaning science.

If that exact error persists:

1. verify the failing test file and implicated proton-cleaning source are
   unchanged relative to starting HEAD;
2. reproduce the same failing test against an isolated checkout/worktree of
   exact starting HEAD `5b2ba7f...` if the local environment permits doing so
   without altering the current worktree;
3. record the exact baseline and repaired-worktree results in the E.8.2 phase
   record;
4. do not repair proton source in this contract.

If the error does not reproduce at the exact starting HEAD, or if its behavior
changes after the E.8.2 repair, hard stop and report it as a possible regression.

Do not suppress the failure with `|| true`.

---

## Local deterministic validation

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

Also run the established repository-memory operations:

- regenerate `docs/memory/manifest.json`;
- manifest `--check`;
- memory health/integrity;
- memory bootstrap check;
- `git diff --check`.

Report the exact interpreter used if it differs from `python`.

Local checks do not establish ROOT/PyROOT, full `main.py`, procedure-PDF farm
rendering, or production runtime validation.

---

## Memory update

Update `docs/memory/phases/e8-2-baseline-stage-audit.md` narrowly to record:

- this Fix.1 source-review repair;
- the repaired finalization failure semantics;
- the pair-safe recovery behavior;
- the fail-closed authority distinction;
- the expanded deterministic producer/finalizer tests;
- the exact test results actually run;
- the unchanged scientific boundaries;
- any reproduced pre-existing proton regression result;
- the farm boundary.

Do not rewrite the record into a new phase design.

E.8.2 remains:

`ACTIVE`

until ChatGPT independently reviews the repaired complete actual diff.

Do not mark E.8.2 `SOURCE REVIEWED`.

E.8.3 remains:

`BLOCKED`

F.6.3 remains:

`BLOCKED`

Do not alter historical runtime evidence.

Regenerate `docs/memory/manifest.json` because versioned memory changes.

---

## Diff audit

Before finishing, run:

```bash
git status --short
git diff --stat
git diff --check
git diff --name-only
```

Audit explicitly that this Fix.1 repair changed only its allowlisted repair
paths beyond the already-intended E.8.2 staging.

Confirm absence of new changes to:

- `src/cuts/rand_sub.py`
- `src/main.py`
- `src/utility/background_config.py`
- proton scientific source;
- pion scientific source outside the allowlist;
- random/dummy coefficients;
- proton factor calculation;
- pion weights;
- pion components/fits/windows;
- canonical binning;
- yield formulas;
- public stage-yield schema;
- Method A;
- Method B numerical use;
- dormant residual-fit configuration.

---

## Required complete review bundle

Overwrite/regenerate the root temporary review bundle:

`kaonlt_review.diff`

It must contain the **complete cumulative E.8.2 worktree diff from starting HEAD
`5b2ba7f...`**, not merely the incremental Fix.1 delta.

Include:

1. complete tracked `git diff`;
2. complete `git diff --no-index -- /dev/null ...` additions for every intended
   untracked file.

At minimum include complete additions for:

- `docs/memory/phases/e8-2-baseline-stage-audit-task-contract.md`
- `docs/memory/phases/e8-2-baseline-stage-audit-fix1-task-contract.md`
- `docs/memory/phases/e8-2-baseline-stage-audit.md`
- `testing/test_e8_2_baseline_stage_audit.py`

Do not add `kaonlt_review.diff` to Git.

The next handoff is the regenerated complete review bundle for independent
ChatGPT actual-diff review.

---

## Acceptance criteria

This repair is ready for independent review only if:

1. branch remains `test`;
2. HEAD remains exactly `5b2ba7f92dc51b303257dbf209df687415364602`;
3. the original E.8.2 scientific architecture is unchanged;
4. no production physics is changed;
5. no second scientific tree traversal is added;
6. renderer failures prevent artifact replacement;
7. finalization never reports success with non-empty renderer failures;
8. partial PDF/manifest installation cannot leave a handled mixed pair;
9. rollback/recovery restores the preliminary pair when possible;
10. preservation status is truthful when recovery fails;
11. transient render state is retained on failure and released only after
    complete success;
12. accepted/current-chain missing authority fails the whole E.8.2 source
    closed;
13. legitimate authoritative zero/skip/rejected child states remain explicit
    child-level unavailable states;
14. missing/non-finite final yield authority for an otherwise valid child fails
    closed;
15. exact baseline pion before/template/after objects remain the authority;
16. public `Y_0` and existing yield error are unchanged;
17. public `stage_window_yields` is unchanged;
18. Fit 1/Fit 2 remain dormant and absent from E.8.2;
19. Method A and Method B remain numerically absent;
20. deterministic producer tests cover the true pre-proton factor algebra and
    one-pass traversal;
21. deterministic tests cover random/dummy/proton closure;
22. deterministic tests cover pion and yield authority failures;
23. deterministic finalizer tests cover renderer failure, second-install
    failure, recovery failure, and success;
24. page-order test includes existing E.8 parent content and terminal handoff;
25. required local and memory checks are run and reported;
26. any unchanged proton-cleaning regression error is handled exactly as the
    baseline-comparison rule above, without an out-of-scope repair;
27. memory remains `ACTIVE` pending ChatGPT review;
28. no farm/runtime validation claim is made;
29. a complete cumulative `kaonlt_review.diff` is generated.

---

## Hard stop

Stop and report rather than broadening scope if:

- repair requires changing proton physics;
- repair requires changing pion physics;
- repair requires changing `w0`;
- repair requires changing random/dummy normalization;
- repair requires changing canonical binning;
- repair requires changing yield formulas;
- repair requires introducing Method A or Method B numerically;
- repair requires activating empirical Fit 1/Fit 2;
- repair requires another source file;
- pair-safe recovery cannot be implemented without changing external workflow;
- the proton-cleaning regression differs from the exact starting-HEAD behavior;
- unrelated worktree changes are discovered.

Codex must not commit, push, update remote refs, or run the Jefferson Lab farm.
