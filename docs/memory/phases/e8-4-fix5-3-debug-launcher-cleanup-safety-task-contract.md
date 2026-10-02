# KaonLT E.8.4 Fix.5.3 — debug-launcher ordinary-checkout cleanup safety

## Objective

Continue the existing local E.8.4 Fix.5/Fix.5.1/Fix.5.2 candidate from exact committed `test` HEAD

```text
da38444e7aa60efd62d6638780776344daf40276
```

without resetting, stashing, cleaning, discarding, reconstructing, or redesigning the current worktree.

Independent ChatGPT actual-diff review of
`kaonlt_review(20261001-222330).diff` found one narrow farm-readiness blocker:

```text
testing/run_e8_4_fix5_left_lowe_plot_gate.py
    -> ./run_Prod_Analysis.sh -d 4p4 2p74
```

correctly preflights and preserves known farm-local model-output dirt, but the tracked
`run_Prod_Analysis.sh` still executes:

```bash
git clean -fdx
```

for the `-d` debug invocation before it calls `set_SymLinks.sh`.

That means the one-command Fix.5 farm owner would indirectly delete untracked/ignored files
from the ordinary farm checkout, contradicting the Fix.5/Fix.5.2 contract and the repository
farm-safety policy that the ordinary checkout must not be cleaned/reset/stashed by this gate.

Repair only that debug-mode cleanup behavior. Preserve the ordinary non-debug launcher path.

No scientific calculation, normalization, cut, fit, template, binning, subtraction, yield,
SIMC normalization, Method-A mathematics, Method-B behavior, efficiency, acceptance,
L/T separation, or cross-section logic may change.

---

## Required starting state

Required branch:

```text
test
```

Required committed HEAD:

```text
da38444e7aa60efd62d6638780776344daf40276
```

The existing cumulative Fix.5/Fix.5.1/Fix.5.2 worktree is expected and must be preserved.

Before editing run:

```bash
git branch --show-current
git rev-parse HEAD
git status --short --untracked-files=all
git -c core.safecrlf=false diff --check
```

Do not reset, stash, clean, checkout-over, discard, or reconstruct any existing candidate file.

---

## Mandatory startup read

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
docs/memory/TOOLS.md
docs/memory/phases/e8-4-fix5-shareable-method-a-impact-pages-task-contract.md
docs/memory/phases/e8-4-fix5-1-simc-runtime-handoff-unblock-task-contract.md
docs/memory/phases/e8-4-fix5-2-collector-scope-owner-completion-task-contract.md
```

Inspect the current cumulative candidate and the live source of:

```text
run_Prod_Analysis.sh
testing/test_run_prod_analysis_debug_left_low.py
testing/run_e8_4_fix5_left_lowe_plot_gate.py
testing/test_run_e8_4_fix5_left_lowe_plot_gate.py
```

---

## Concrete source-level blocker

The current launcher contains the ordinary cleanup block:

```bash
if [[ $i_flag != "true" && $a_flag != "true" ]]; then
    if ! validate_external_sigma0_paths_before_cleanup; then
        exit 1
    fi
    git clean -fdx
    ./set_SymLinks.sh $ParticleType
fi
```

Because `-d` does not set `i_flag` or `a_flag`, the reviewed Fix.5 owner would call this
`git clean -fdx` during the ordinary farm checkout.

This violates the gate contract even though the owner itself contains no cleanup call.

The repair must make `-d` non-destructive while leaving the ordinary non-debug path unchanged.

---

## Allowed source changes

Edit only:

```text
run_Prod_Analysis.sh
testing/test_run_prod_analysis_debug_left_low.py
testing/test_run_e8_4_fix5_left_lowe_plot_gate.py
```

The owner source itself:

```text
testing/run_e8_4_fix5_left_lowe_plot_gate.py
```

should remain byte-identical unless a deterministic test demonstrates that a tiny source change is
strictly required to prove/enforce the repaired launcher contract. Prefer no owner edit.

The existing Fix.5 scientific/presentation candidate must remain byte-identical.

Create this task contract:

```text
docs/memory/phases/e8-4-fix5-3-debug-launcher-cleanup-safety-task-contract.md
```

If versioned memory changes solely because this contract is added, regenerate
`docs/memory/manifest.json`. Do not otherwise rewrite CURRENT/status: the repair does not change
the intended NEXT or scientific state.

---

## Frozen scientific and presentation files

Do not modify:

```text
src/main.py
src/cuts/full_background_subtraction_plots.py
testing/test_e8_4_production_impact_audit.py
testing/test_e8_4_fix5_main_order.py

testing/collect_pion_hgcer_validation_bundle.py
testing/package_pion_hgcer_validation_bundle.tcsh
testing/pion_hgcer_validation_bundle_profile_e8_4_fix5.json

src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
src/cuts/pion_hgcer_method_a_parent_preserving_correction.py
src/cuts/pion_hgcer_method_a_tphi_propagation.py
src/cuts/pion_component_subtraction.py
src/cuts/pion_component_fits.py
src/cuts/particle_subtraction.py
src/cuts/calculate_yield.py
binning/calculate_yield.py
```

Do not modify Method-B source or accepted F.1-F.6.2 authorities.

Preserve all existing memory/evidence/status content except the manifest update required by adding
this contract.

---

## Required launcher behavior

Change only the cleanup behavior for `-d`.

Required behavior:

### Ordinary non-debug path

For ordinary non-debug execution where the existing code historically performs cleanup:

```text
validate_external_sigma0_paths_before_cleanup
-> git clean -fdx
-> ./set_SymLinks.sh $ParticleType
```

must remain unchanged.

### `-d` debug path

For `./run_Prod_Analysis.sh -d 4p4 2p74`:

```text
validate_external_sigma0_paths_before_cleanup
-> DO NOT run git clean -fdx
-> still run ./set_SymLinks.sh $ParticleType
-> preserve all existing paired low/high canonical preflight behavior
-> run full analysis only Left / lowe
-> stop before full high-epsilon processing
```

A narrow implementation is expected, for example guarding only the existing `git clean -fdx`
line with `d_flag != "true"` while leaving `set_SymLinks.sh` outside that new guard.

Do not bypass the external-Sigma0 safety check.
Do not bypass symlink setup.
Do not change any other flag semantics.
Do not change ordinary non-debug cleanup behavior.

---

## Required deterministic tests

Extend `testing/test_run_prod_analysis_debug_left_low.py` to prove:

1. the launcher still contains exactly one `git clean -fdx`;
2. ordinary non-debug execution still reaches that cleanup;
3. `-d` explicitly bypasses that cleanup;
4. `-d` still runs `set_SymLinks.sh $ParticleType`;
5. the external-Sigma0 validation still occurs before cleanup/symlink setup;
6. paired low/high canonical preflight remains unchanged;
7. only Left/lowe gets full debug analysis;
8. high-epsilon full analysis remains skipped in `-d`;
9. ordinary full low/high paths remain present;
10. shell syntax remains valid.

Extend `testing/test_run_e8_4_fix5_left_lowe_plot_gate.py` with a static/runtime-path assertion that
the tracked farm owner invokes the `-d` launcher and that the launcher source used by that owner has
the non-destructive debug cleanup guard. Do not mock away this specific contract.

Run all existing owner/profile/E.8 regressions as listed below.

---

## Positive / regression validation

Use the repository-authoritative Python.

Run at minimum:

```bash
<PYTHON> -B -m py_compile \
  testing/test_run_prod_analysis_debug_left_low.py \
  testing/run_e8_4_fix5_left_lowe_plot_gate.py \
  testing/test_run_e8_4_fix5_left_lowe_plot_gate.py

<PYTHON> -B -m unittest testing.test_run_prod_analysis_debug_left_low -v
<PYTHON> -B -m unittest testing.test_run_e8_4_fix5_left_lowe_plot_gate -v
<PYTHON> -B -m unittest testing.test_pion_hgcer_validation_bundle_profile_e8_4_fix5 -v
<PYTHON> -B -m unittest testing.test_e8_4_fix5_main_order -v
<PYTHON> -B -m unittest testing.test_e8_4_production_impact_audit -v
<PYTHON> -B -m unittest testing.test_full_background_subtraction_plots -v
<PYTHON> -B -m unittest testing.test_f6_3_parallel_full_procedure_method_a -v
<PYTHON> -B -m unittest testing.test_collect_pion_hgcer_validation_bundle -v
```

Then:

```bash
<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/memory_bootstrap.py --root . --json
<PYTHON> -m unittest testing.test_memory_health -v
<PYTHON> -B tools/check_memory_health.py --root .
git -c core.safecrlf=false diff --check
```

Report exact test counts/skips.

Do not run `main.py`.
Do not run the farm.
Do not claim ROOT/PyROOT/farm validation.

---

## Diff audit

Before stopping:

```bash
git status --short --untracked-files=all
git -c core.safecrlf=false diff --stat
git -c core.safecrlf=false diff --check
git -c core.safecrlf=false diff --no-ext-diff
```

Explicitly prove:

- the existing Fix.5/Fix.5.1/Fix.5.2 scientific/presentation candidate is byte-identical to the
  supplied `kaonlt_review(20261001-222330).diff` candidate;
- generic collector unchanged;
- validation profile unchanged;
- farm owner unchanged unless strictly required by this contract;
- Method-A mathematics unchanged;
- Method-B unchanged;
- production physics unchanged;
- ordinary non-debug launcher cleanup behavior unchanged;
- debug launcher no longer executes `git clean -fdx`;
- debug launcher still executes symlink setup and the existing canonical preflight/full-analysis path;
- manifest fresh;
- ordinary memory health passes.

Create one fresh complete cumulative review bundle:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

It must include:

- branch/HEAD/status;
- complete cumulative tracked diff from committed HEAD;
- complete no-index diffs for every intended new/untracked file;
- changed-path inventory;
- all deterministic test output;
- manifest/memory-health output;
- explicit source-preservation audit.

---

## Acceptance criteria

PASS only if:

1. committed HEAD remains exactly `da38444e7aa60efd62d6638780776344daf40276`;
2. existing Fix.5/Fix.5.1/Fix.5.2 candidate is preserved;
3. `run_Prod_Analysis.sh` ordinary non-debug cleanup path is unchanged;
4. `-d` no longer executes `git clean -fdx`;
5. `-d` still runs external-Sigma0 validation and `set_SymLinks.sh`;
6. paired low/high canonical preflight remains unchanged;
7. full analysis remains Left/lowe only in `-d`;
8. high-epsilon full processing remains skipped in `-d`;
9. owner still invokes exactly `./run_Prod_Analysis.sh -d 4p4 2p74`;
10. owner/profile/collector behavior from Fix.5.2 is unchanged;
11. scientific/presentation source is unchanged from the reviewed candidate;
12. all focused/regression tests pass;
13. manifest and ordinary memory health pass;
14. one complete fresh cumulative review bundle is produced;
15. Codex does not commit, push, or run the farm.

---

## Hard stop

After the deterministic repair and fresh cumulative review bundle:

**STOP.**

Do not commit.
Do not push.
Do not run the farm.
Do not broaden scope.
Do not change physics.
Do not promote Method A.
Do not close F.6.4.

Return the fresh review bundle for independent ChatGPT actual-diff review.
