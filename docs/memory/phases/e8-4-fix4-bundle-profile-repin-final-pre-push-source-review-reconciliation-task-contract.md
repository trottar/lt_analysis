# KaonLT E.8.4.Fix.4 Bundle Profile Re-Pin — Final Pre-Push Source-Review Reconciliation

## 1. Purpose

Perform the repository-memory/status reconciliation required after independent
ChatGPT actual-diff/source-provenance review **PASSED** the local
E.8.4.Fix.4 validation-bundle profile re-pin candidate.

Reviewed cumulative bundle:

```text
kaonlt_review(20260929-222602).diff
```

Independent review result:

```text
PASS
```

The review established that the candidate:

- starts from committed `test` HEAD
  `6e2adf7a37ac9e79cad99242686804cf51701644`;
- changes substantive files only in:
  - `testing/pion_hgcer_validation_bundle_profile_e8_1.json`
  - `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`;
- changes the profile's `source_identity.required_analysis_commit` exactly from
  `29d7b7f9635db899939efeb3508e941e994e8928`
  to
  `6e2adf7a37ac9e79cad99242686804cf51701644`;
- changes the focused test's `REVIEWED_SOURCE` to the same pushed Fix.4 source;
- leaves canonical-five settings unchanged;
- leaves global/setting artifact declarations unchanged;
- leaves the exact profile/test-only `allowed_committed_files` unchanged;
- leaves the sole `docs/memory/` non-analysis prefix unchanged;
- leaves collector, wrapper, Fix.4 analysis/test source, scientific/runtime
  source, and accepted F.6.2 science unchanged;
- records pushed-state provenance for Fix.4 without upgrading runtime status.

Codex-reported deterministic checks passed:

- profile/collector `py_compile`;
- profile `json.tool`;
- focused profile suite: 5 tests OK;
- generic collector suite: 28 tests OK;
- dynamic alignment suite: 19 tests OK with 11 PyROOT-dependent tests skipped
  because PyROOT was unavailable;
- manifest check;
- memory health/bootstrap;
- 35 memory-health tests;
- `git diff --check`.

These checks were **NOT RUN by ChatGPT**. They establish no ROOT/PyROOT,
full-analysis, procedure-PDF, farm, runtime, production, or Method-A-promotion
acceptance.

This task is memory/status reconciliation only. It must not alter the reviewed
profile/test candidate.

---

## 2. Exact starting state

Required branch:

```text
test
```

Required committed HEAD:

```text
6e2adf7a37ac9e79cad99242686804cf51701644
```

Remote `test` was independently observed at that same HEAD during ChatGPT
review.

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
```

The expected cumulative candidate worktree before reconciliation contains:

Tracked modifications:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/e8-4-fix4-alignment-cache-semantics.md
docs/memory/roadmap/STATUS.md
testing/pion_hgcer_validation_bundle_profile_e8_1.json
testing/test_pion_hgcer_validation_bundle_profile_e8_1.py
```

Intended untracked candidate files:

```text
docs/memory/phases/e8-4-fix4-bundle-profile-repin-task-contract.md
docs/memory/phases/e8-4-fix4-bundle-profile-repin.md
```

Temporary root `kaonlt_review*.diff` files may remain untracked. Root
`AGENTS.md` and `.codex/` are local-only/untracked. None of those temporary/local
files may be staged, removed, rewritten, or treated as candidate files.

If committed HEAD differs or unrelated tracked changes exist, **STOP** and
report the blocker.

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
- `docs/memory/phases/e8-4-fix4-alignment-cache-semantics.md`
- `docs/memory/phases/e8-4-fix4-bundle-profile-repin-task-contract.md`
- `docs/memory/phases/e8-4-fix4-bundle-profile-repin.md`
- the applicable source-changing workflow record.

Do not reopen adjacent architecture or science.

---

## 4. Files that must remain byte-identical

The following independently reviewed candidate files must not be edited in this
reconciliation:

```text
testing/pion_hgcer_validation_bundle_profile_e8_1.json
testing/test_pion_hgcer_validation_bundle_profile_e8_1.py
docs/memory/phases/e8-4-fix4-alignment-cache-semantics.md
docs/memory/phases/e8-4-fix4-bundle-profile-repin-task-contract.md
```

Also frozen:

```text
testing/collect_pion_hgcer_validation_bundle.py
testing/package_pion_hgcer_validation_bundle.tcsh
src/cuts/pion_component_fits.py
testing/test_pion_component_dynamic_alignment.py
src/cuts/full_background_subtraction_plots.py
src/cuts/rand_sub.py
src/cuts/pion_hgcer_method_a_acceptance_contract.py
src/cuts/pion_hgcer_method_a_acceptance_map.py
src/cuts/pion_hgcer_method_a_parent_preserving_correction.py
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
src/cuts/pion_component_subtraction.py
src/cuts/particle_subtraction.py
src/binning/calculate_yield.py
src/binning/ave_per_bin.py
src/main.py
run_Prod_Analysis.sh
src/utility/background_config.py
```

No scientific/runtime/profile behavior change is allowed.

If any byte in the reviewed profile/test candidate needs changing, **STOP**.
That would require a new repair, not reconciliation.

---

## 5. Allowed reconciliation changes

Only these existing files may be edited:

```text
docs/memory/CURRENT.md
docs/memory/roadmap/STATUS.md
docs/memory/phases/e8-4-fix4-bundle-profile-repin.md
docs/memory/manifest.json
```

Create exactly one new reconciliation contract:

```text
docs/memory/phases/e8-4-fix4-bundle-profile-repin-final-pre-push-source-review-reconciliation-task-contract.md
```

No other file may change.

---

## 6. Required durable status reconciliation

### 6.1 Profile re-pin

Change:

```text
E.8.4.Fix.4 validation-bundle profile re-pin — ACTIVE
```

to:

```text
E.8.4.Fix.4 validation-bundle profile re-pin — SOURCE REVIEWED
```

Record that independent ChatGPT actual-diff/source-provenance review of:

```text
kaonlt_review(20260929-222602).diff
```

passed.

Record the source-review boundaries:

- exact required-analysis source pin is pushed Fix.4
  `6e2adf7a37ac9e79cad99242686804cf51701644`;
- `REVIEWED_SOURCE` matches that exact source;
- only the two approved profile/test substantive lines changed;
- canonical-five setting inventory is unchanged;
- global/setting artifact inventory is unchanged;
- exact committed-file allowlist is unchanged;
- the sole non-analysis prefix remains `docs/memory/`;
- collector and wrapper are unchanged;
- Fix.4 analysis/test source is unchanged;
- no scientific/runtime behavior changed.

### 6.2 Validation boundary

Explicitly state:

- Codex-reported local deterministic checks passed;
- focused profile suite reported 5 tests OK;
- collector suite reported 28 tests OK;
- dynamic alignment suite reported 19 tests OK with 11 PyROOT-dependent skips
  because PyROOT was unavailable;
- tests were **NOT RUN by ChatGPT**;
- no ROOT/PyROOT, full `main.py`, procedure-PDF, farm, runtime, production, or
  Method-A-promotion acceptance is claimed.

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
- E.8.4.Fix.4 persisted-alignment repair — **SOURCE REVIEWED**
- final E.8 — **BLOCKED**
- F.6.4 — **BLOCKED**
- lifecycle-hook dispatch — **BLOCKED / DEFERRED**

Do not upgrade E.8.4 or F.6.3 to runtime validated.

### 6.4 Preserve the current runtime blocker

The fresh `Q4p4W2p74 / Left / lowe` E.8.4 gate remains **BLOCKED** by:

```text
f6_3_f4_shared_reproduction_failed:f3_fingerprint_input_content_mismatch
```

until the reviewed profile re-pin is user-committed/pushed, pushed-state
reviewed, and a fresh narrow farm gate supplies applicable runtime evidence.

Do not claim that the profile re-pin or Fix.4 has resolved the blocker on the
farm.

---

## 7. CURRENT.md exact NEXT

After reconciliation, `CURRENT.md` must contain exactly one ordinary NEXT:

```text
NEXT — user-controlled commit/push of the independently reviewed E.8.4.Fix.4 validation-bundle profile re-pin.
```

This repository-memory NEXT does not authorize Codex to commit or push.

The current Codex task still hard-stops for independent ChatGPT review of the
final pre-push reconciliation bundle before the user performs that Git action.

Do not add another NEXT elsewhere.

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

Do not run the farm.

A profile/source test rerun is not required for this memory-only
reconciliation. If Codex chooses to rerun deterministic tests, report them
accurately without claiming new runtime evidence.

---

## 9. Required cumulative diff audit

Before stopping:

```bash
git status --short
git -c core.safecrlf=false diff --stat
git -c core.safecrlf=false diff --no-ext-diff
```

The cumulative candidate must contain exactly these intended paths:

```text
testing/pion_hgcer_validation_bundle_profile_e8_1.json
testing/test_pion_hgcer_validation_bundle_profile_e8_1.py
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/e8-4-fix4-alignment-cache-semantics.md
docs/memory/roadmap/STATUS.md
docs/memory/phases/e8-4-fix4-bundle-profile-repin-task-contract.md
docs/memory/phases/e8-4-fix4-bundle-profile-repin.md
docs/memory/phases/e8-4-fix4-bundle-profile-repin-final-pre-push-source-review-reconciliation-task-contract.md
```

Temporary `kaonlt_review*.diff` files may exist at repository root but must
remain untracked and must not be part of the intended candidate inventory.

No `git add -A`.

---

## 10. Required final pre-push review bundle

Create one fresh timestamped cumulative review bundle at repository root:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

It must contain byte-faithful command output with no Markdown rewriting.

Include:

1. branch;
2. committed HEAD;
3. `git status --short`;
4. cumulative `git diff --stat`;
5. complete cumulative tracked
   `git -c core.safecrlf=false diff --no-ext-diff`;
6. complete `git diff --no-index /dev/null ...` sections for every intended
   untracked candidate file:
   - `docs/memory/phases/e8-4-fix4-bundle-profile-repin-task-contract.md`
   - `docs/memory/phases/e8-4-fix4-bundle-profile-repin.md`
   - `docs/memory/phases/e8-4-fix4-bundle-profile-repin-final-pre-push-source-review-reconciliation-task-contract.md`
7. exact reconciliation check commands/results;
8. final cumulative changed-path inventory.

The review bundle must preserve raw Git diff text and literal quote characters.

Do not stage the review bundle.

---

## 11. Acceptance criteria

The reconciliation candidate is acceptable only if:

1. committed HEAD remains exactly
   `6e2adf7a37ac9e79cad99242686804cf51701644`;
2. the reviewed profile JSON diff is byte-identical;
3. the reviewed focused-test diff is byte-identical;
4. the original profile-repin task contract is byte-identical;
5. the Fix.4 alignment-semantic phase pushed-state update is byte-identical;
6. only the four allowed existing memory files plus this new reconciliation
   contract change during reconciliation;
7. profile re-pin is `SOURCE REVIEWED`;
8. the ChatGPT source/provenance PASS is recorded without upgrading runtime
   status;
9. tests are described as Codex-reported and `NOT RUN by ChatGPT`;
10. the fresh Left/lowe farm gate remains `BLOCKED`;
11. E.8/F.6 ownership and statuses remain correct;
12. `CURRENT.md` has exactly the required single NEXT;
13. manifest regeneration/check passes;
14. memory-health/bootstrap/tests and `diff --check` pass, allowing only the
    known non-fatal CURRENT soft-size warning;
15. no source, scientific, collector, wrapper, package, farm, commit, or push
    action occurs;
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
