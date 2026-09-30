# KaonLT E.8.4.Fix.4 — Final Pre-Push Source-Review Reconciliation

## 1. Purpose

Perform the repository-memory/status reconciliation required after independent ChatGPT actual-diff/source-runtime-path review **PASSED** the local E.8.4.Fix.4 persisted-alignment semantic-version implementation candidate.

This task is memory/status reconciliation only. It must not alter the independently reviewed source/test candidate.

Reviewed cumulative bundle:

```text
kaonlt_review(20260929-215203).diff
```

Independent review result:

```text
PASS
```

The review established that the candidate:

- starts from committed `test` HEAD
  `8d62dcfdf2fd8d08298c075940b4ab1b28d7f079`;
- changes substantive source only in:
  - `src/cuts/pion_component_fits.py`
  - `testing/test_pion_component_dynamic_alignment.py`;
- adds source-owned alignment resolver semantics
  `pion_component_dynamic_alignment_semantics/v2`;
- carries that semantics value in expected metadata, normal resolver output,
  disabled/fallback resolver output, persisted JSON, and CSV diagnostics;
- rejects missing, stale, wrong, or malformed persisted semantics;
- preserves current-semantics cache reuse;
- preserves Fix.3 checksum+axis pion-control identity semantics;
- preserves the Fix.3 minimum-template-integral predicate;
- requires current semantics for direct fine-bin parent validity;
- does not alter alignment schema v2, scientific configuration, scan grids,
  thresholds, templates, component ordering, F.3/F.4/F.6.3/E.8.4 source,
  validation profile, baseline physics, Method A mathematics, or Method B;
- retains the fresh `Q4p4W2p74 / Left / lowe` farm observation as a blocker,
  not runtime closure.

Codex-reported deterministic checks passed, with 11 PyROOT-dependent focused
alignment tests skipped because PyROOT was unavailable. Those checks were
**NOT RUN by ChatGPT** and do not constitute ROOT/PyROOT, farm, or runtime
validation.

The purpose of this reconciliation is to record that source review result
durably, advance Fix.4 to `SOURCE REVIEWED`, and prepare the exact candidate for
a final independent pre-push reconciliation review.

---

## 2. Exact starting state

Required branch:

```text
test
```

Required committed HEAD:

```text
8d62dcfdf2fd8d08298c075940b4ab1b28d7f079
```

The remote `test` branch was independently observed at that same HEAD during
the ChatGPT review.

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
```

The expected cumulative candidate worktree before this reconciliation contains:

Tracked modifications:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/e8-4-fix3-pion-alignment-determinism.md
docs/memory/roadmap/STATUS.md
src/cuts/pion_component_fits.py
testing/test_pion_component_dynamic_alignment.py
```

Intended untracked candidate files:

```text
docs/memory/evidence/e8-4-fix3-left-lowe-stale-alignment-cache-runtime-blocker.md
docs/memory/phases/e8-4-fix4-alignment-cache-semantics-task-contract.md
docs/memory/phases/e8-4-fix4-alignment-cache-semantics.md
```

Root-local `AGENTS.md`, `.codex/`, and temporary review bundles may also exist
and must not be removed, staged, rewritten, or treated as candidate files.

If the committed HEAD differs, or unrelated tracked changes are present,
**STOP** and report the blocker.

Do not reset, stash, clean, discard, commit, push, package, rerender, or run the
Jefferson Lab farm.

---

## 3. Mandatory repository-memory startup

Read in this exact order:

1. root `AGENTS.md`;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read only task-relevant records:

- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/e8-4-fix3-pion-alignment-determinism.md`
- `docs/memory/phases/e8-4-fix4-alignment-cache-semantics-task-contract.md`
- `docs/memory/phases/e8-4-fix4-alignment-cache-semantics.md`
- `docs/memory/evidence/e8-4-fix3-left-lowe-stale-alignment-cache-runtime-blocker.md`
- the applicable source-changing workflow record.

Do not reopen adjacent architecture or science.

---

## 4. Files that must remain byte-identical

The following independently reviewed candidate files must not be edited in this
reconciliation:

```text
src/cuts/pion_component_fits.py
testing/test_pion_component_dynamic_alignment.py
docs/memory/phases/e8-4-fix3-pion-alignment-determinism.md
docs/memory/evidence/e8-4-fix3-left-lowe-stale-alignment-cache-runtime-blocker.md
docs/memory/phases/e8-4-fix4-alignment-cache-semantics-task-contract.md
```

Also frozen:

```text
src/utility/background_config.py
src/cuts/pion_hgcer_method_a_acceptance_contract.py
src/cuts/pion_hgcer_method_a_acceptance_map.py
src/cuts/pion_hgcer_method_a_parent_preserving_correction.py
src/cuts/pion_hgcer_method_a_tphi_propagation.py
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
src/cuts/full_background_subtraction_plots.py
src/cuts/rand_sub.py
src/cuts/pion_component_subtraction.py
src/cuts/particle_subtraction.py
src/binning/calculate_yield.py
src/binning/ave_per_bin.py
src/main.py
run_Prod_Analysis.sh
testing/pion_hgcer_validation_bundle_profile_e8_1.json
testing/test_pion_hgcer_validation_bundle_profile_e8_1.py
collector/wrapper source
accepted F.1--F.6.2 artifacts
```

No scientific/runtime source change is allowed.

If any byte in the reviewed source/test candidate needs changing, **STOP**.
That would require a new implementation repair, not reconciliation.

---

## 5. Allowed reconciliation changes

Only these existing files may be edited:

```text
docs/memory/CURRENT.md
docs/memory/roadmap/STATUS.md
docs/memory/phases/e8-4-fix4-alignment-cache-semantics.md
docs/memory/manifest.json
```

Create exactly one new reconciliation contract:

```text
docs/memory/phases/e8-4-fix4-final-pre-push-source-review-reconciliation-task-contract.md
```

No other file may change.

---

## 6. Required durable status reconciliation

### 6.1 E.8.4.Fix.4

Change E.8.4.Fix.4 from:

```text
ACTIVE
```

to:

```text
SOURCE REVIEWED
```

Record that independent ChatGPT actual-diff/source-runtime-path review of:

```text
kaonlt_review(20260929-215203).diff
```

passed.

Record the source-review boundaries:

- source/test diff matches the narrow contract;
- no production physics changed;
- no F.3/F.4/F.6.3/E.8.4 consumer source changed;
- no alignment schema/config/scan/threshold change occurred;
- persisted pre-Fix semantics fail closed and current semantics remain reusable;
- parent-to-fine semantics boundary is fail closed;
- Fix.3 pion-control checksum+axis identity and minimum-integral predicate remain
  unchanged.

### 6.2 Validation boundary

Explicitly state:

- Codex-reported local deterministic checks passed;
- `testing.test_pion_component_dynamic_alignment` reported 19 tests OK with 11
  PyROOT histogram-path tests skipped because PyROOT was unavailable;
- the F.6.3 and E.8.4 local deterministic suites passed according to the
  uploaded review bundle;
- tests were **NOT RUN by ChatGPT**;
- no ROOT/PyROOT, full `main.py`, procedure-PDF, farm, runtime, production, or
  Method-A promotion acceptance is claimed.

### 6.3 Preserve project statuses

Keep:

- F.1 through F.6.2, including F.6.2.Fix.5 —
  **CLOSED / RUNTIME VALIDATED**
- E.8 — **ACTIVE**
- E.8.2 — **SOURCE REVIEWED**
- E.8.3 — **SOURCE REVIEWED**
- F.6.3 — **SOURCE REVIEWED**
- E.8.4 — **SOURCE REVIEWED**
- E.8.4.Fix.3 — **SOURCE REVIEWED**
- final E.8 — **BLOCKED**
- F.6.4 — **BLOCKED**
- lifecycle-hook dispatch — **BLOCKED / DEFERRED**

Do not upgrade E.8.4 or F.6.3 to runtime validated.

### 6.4 Preserve the fresh runtime blocker

The fresh `Q4p4W2p74 / Left / lowe` E.8.4 gate remains `BLOCKED` by:

```text
f6_3_f4_shared_reproduction_failed:f3_fingerprint_input_content_mismatch
```

until the reviewed Fix.4 source is user-committed/pushed, pushed-state reviewed,
any required validation-profile re-pin is separately reviewed, and a fresh
narrow farm gate supplies applicable runtime evidence.

Do not claim that Fix.4 has resolved the blocker on the farm.

---

## 7. CURRENT.md exact NEXT

After reconciliation, `CURRENT.md` must contain exactly one ordinary NEXT:

```text
NEXT — user-controlled commit/push of the independently reviewed E.8.4.Fix.4 persisted-alignment semantic-version repair.
```

This repository-memory NEXT does not authorize Codex to commit or push.

The current Codex task still hard-stops for independent ChatGPT review of the
final pre-push reconciliation bundle before the user performs that Git action.

Do not add another `NEXT` elsewhere.

---

## 8. Manifest and deterministic reconciliation checks

Regenerate and verify memory:

```bash
python -B tools/update_memory_manifest.py --root . --write
python -B tools/update_memory_manifest.py --root . --check
python -B tools/check_memory_health.py --root .
python -B tools/memory_bootstrap.py --root . --json
python -B -m unittest testing.test_memory_health -v
git -c core.safecrlf=false diff --check
```

The existing `CURRENT.md` soft-size warning is non-fatal only if
`check_memory_health.py` exits 0.

Do not rerun the farm.

A source test rerun is not required for this memory-only reconciliation. If
Codex chooses to rerun any deterministic source suite, report it accurately,
but do not claim new runtime evidence.

---

## 9. Required cumulative diff audit

Before stopping, inspect:

```bash
git status --short
git -c core.safecrlf=false diff --stat
git -c core.safecrlf=false diff --no-ext-diff
```

The cumulative candidate must contain exactly these intended paths:

```text
src/cuts/pion_component_fits.py
testing/test_pion_component_dynamic_alignment.py
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/e8-4-fix3-pion-alignment-determinism.md
docs/memory/roadmap/STATUS.md
docs/memory/evidence/e8-4-fix3-left-lowe-stale-alignment-cache-runtime-blocker.md
docs/memory/phases/e8-4-fix4-alignment-cache-semantics-task-contract.md
docs/memory/phases/e8-4-fix4-alignment-cache-semantics.md
docs/memory/phases/e8-4-fix4-final-pre-push-source-review-reconciliation-task-contract.md
```

Temporary `kaonlt_review*.diff` artifacts may exist in the repository root but
must remain untracked and must not be part of the intended candidate inventory.

No `git add -A`.

---

## 10. Required final pre-push review bundle

Create one **fresh timestamped cumulative** review bundle at repository root:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

It must contain byte-faithful command output and no Markdown rewriting.

Include:

1. branch;
2. committed HEAD;
3. `git status --short`;
4. cumulative `git diff --stat`;
5. complete cumulative tracked
   `git -c core.safecrlf=false diff --no-ext-diff`;
6. complete `git diff --no-index /dev/null ...` sections for every intended
   untracked candidate file:
   - `docs/memory/evidence/e8-4-fix3-left-lowe-stale-alignment-cache-runtime-blocker.md`
   - `docs/memory/phases/e8-4-fix4-alignment-cache-semantics-task-contract.md`
   - `docs/memory/phases/e8-4-fix4-alignment-cache-semantics.md`
   - `docs/memory/phases/e8-4-fix4-final-pre-push-source-review-reconciliation-task-contract.md`
7. exact reconciliation check commands/results;
8. final cumulative changed-path inventory.

The bundle must preserve literal quote characters and raw Git diff text.

Do not stage the bundle.

---

## 11. Acceptance criteria

The reconciliation candidate is acceptable only if:

1. committed HEAD remains exactly
   `8d62dcfdf2fd8d08298c075940b4ab1b28d7f079`;
2. the previously reviewed source/test diff is byte-identical;
3. the original Fix.4 task contract is byte-identical;
4. the blocker evidence is byte-identical;
5. the Fix.3 phase update is byte-identical;
6. only the four allowed existing memory files plus this new reconciliation
   contract changed during reconciliation;
7. Fix.4 is `SOURCE REVIEWED`;
8. the ChatGPT source-review PASS is recorded without upgrading runtime status;
9. tests are described as Codex-reported and `NOT RUN by ChatGPT`;
10. the fresh Left/lowe farm gate remains `BLOCKED`;
11. E.8/F.6 ownership and statuses remain correct;
12. `CURRENT.md` has exactly the required single NEXT;
13. the manifest is regenerated and passes;
14. memory-health/bootstrap/tests and `diff --check` pass, allowing only the
    known non-fatal CURRENT soft-size warning;
15. no source, physics, validation-profile, package, farm, commit, or push action
    occurs;
16. a fresh faithful cumulative review bundle is created for independent
    ChatGPT review.

---

## 12. Hard stop

After memory/status reconciliation, manifest regeneration, deterministic memory
checks, cumulative diff audit, and creation of the fresh timestamped review
bundle:

**STOP.**

Do not commit, push, package, rerender, or run the farm.

Return the fresh review bundle for independent ChatGPT final pre-push
reconciliation review.
