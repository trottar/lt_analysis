# KaonLT E.8.4.Fix.3 — Pion-Alignment Determinism Repair

## 1. Purpose

Repair one runtime-exposed determinism defect in the existing pion-component dynamic-alignment machinery before any further E.8.4 / F.6.3 farm validation.

Fresh `Q4p4W2p74 / Left / lowe` farm evidence from
`KaonLT_E8_4_Left_lowe_alignment_20260929-102321.json` established two coupled source defects:

1. **Persistence compatibility is sensitive to an ephemeral ROOT histogram name.**
   The September 24 and September 29 setting-wide pion-control histograms have the same axis and the same content checksum

   `0dc9c642ade833a0b8fd727be1ba55dd79a09392082c07e12afafa82cbdb0ab1`

   but different generated ROOT names. The September 29 store therefore reports
   `persistence_status="rejected_stale_then_created"` with
   `pion control histogram identifier mismatch` instead of reusing the semantically identical stored alignment.

2. **The renormalized-template minimum-integral acceptance boundary is not floating-point safe.**
   The active configuration has
   `minimum_template_integral=1.0` and `renormalize_shifted_templates=true`.
   Candidate validity changes across numerically equivalent runs because the current strict
   `_hist_integral(template) < minimum_template_integral` comparison can classify a unit-normalized shifted template differently for machine-scale roundoff around `1.0`.

The farm artifact demonstrates the consequence directly:

- September 24 accepts the `pi_n` residual-shift candidate `-0.003 GeV` and rejects the neighboring `-0.004 GeV` candidate for `insufficient template integral`.
- September 29 rejects the effectively identical `-0.003 GeV` candidate for `insufficient template integral` and accepts `-0.004 GeV`.
- Since the alignment resolver applies components sequentially to the residual, the `pi_n` flip changes the residual seen by `pi_delta`; the parent `pi_delta` solution then moves from `+0.008 GeV` to `+0.016 GeV`.
- The parent baseline score remains effectively unchanged (`73.77888464083372` versus `73.77888471147615`) and the common setting shift differs by only about `4.7e-12 GeV`, so this is a determinism defect, not evidence for a changed pion-control population.

This task repairs only those two defects. It must not alter the accepted scientific model, scan definition, component ordering, fit windows, normalization, production subtraction formulas, or Method-A/Method-B ownership.

---

## 2. Exact starting state

Required branch:

`test`

Required starting HEAD:

`9daca79051a4eb451a158dd7c34b7028c6b96456`

Commit message at that HEAD:

`Re-pin E8.4 validation bundle profile`

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
```

Hard requirements:

- branch must be `test`;
- committed HEAD must be exactly `9daca79051a4eb451a158dd7c34b7028c6b96456`;
- inspect the existing worktree before editing;
- do **not** reset, stash, clean, discard, or overwrite user work;
- the local root `AGENTS.md` and task contract may be untracked/local-only and must not be removed;
- if the committed HEAD differs, or unrelated tracked changes are already present, **STOP** and report the blocker.

Do not commit, push, run the Jefferson Lab farm, package a farm bundle, or mutate farm artifacts.

---

## 3. Mandatory repository-memory startup

Read in this order before implementation:

1. root `AGENTS.md`;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read only the task-relevant records/source:

- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/e8-4-production-impact-audit.md`
- the current E.8.4 bundle-profile re-pin phase record
- the relevant pion-alignment history/decision/evidence records referenced by CURRENT or source comments
- `src/cuts/pion_component_fits.py`
- `src/utility/background_config.py`
- `testing/test_pion_component_dynamic_alignment.py`
- only additional directly relevant producer/consumer source needed to verify the runtime path.

Reconcile the already-pushed bundle-profile re-pin in CURRENT: HEAD `9daca790...` means the old `NEXT — user-controlled commit/push of the reviewed bundle-profile re-pin` is stale.

---

## 4. Current status and scientific ownership

Preserve these statuses:

- F.6.2 / F.6.2.Fix.5 — **CLOSED / RUNTIME VALIDATED**
- E.8.2 — **SOURCE REVIEWED**
- E.8.3 — **SOURCE REVIEWED**
- F.6.3 — **SOURCE REVIEWED**
- E.8.4 — **SOURCE REVIEWED**
- E.8 — **ACTIVE**
- final E.8 — **BLOCKED**
- F.6.4 — **BLOCKED**
- lifecycle-hook dispatch — **BLOCKED / DEFERRED**

Create the new narrow work item:

- E.8.4.Fix.3 pion-alignment determinism repair — **ACTIVE** pending independent ChatGPT actual-diff review.

This fix is a **runtime determinism/provenance repair inside the existing baseline pion-alignment machinery**. It does not reopen or redefine accepted F.1–F.6.2 physics, and it does not promote Method A.

The September 29 F.6.3/E.8.4 runtime gate remains **BLOCKED** until this repair is implemented, independently source-reviewed, committed/pushed by the user, pushed-state reviewed, and the exact validation profile is ready for that later gate.

---

## 5. Source audit facts that must remain true

### 5.1 Current pion-control identifier

At the starting source, `_pion_control_histogram_identifier(hist)` serializes:

- `hist_name`;
- content `checksum`;
- histogram `axis`.

The generated ROOT name is useful diagnostic provenance, but it is not stable semantic content identity.

### 5.2 Current persistence compatibility

`_alignment_compatibility_reasons(stored, expected)` currently validates, among other things:

- alignment schema version;
- resolved configuration hash;
- analysis scope;
- complete physical-bin identity;
- histogram axis specification;
- parent alignment hash when required;
- immutable source-template identifiers/checksums;
- full equality of `pion_control_histogram_identifier`.

Only the final pion-control comparison is in scope for the identity repair. All other fail-closed compatibility checks remain authoritative.

### 5.3 Current template construction and acceptance

`_build_alignment_template(...)` delegates to `build_shifted_template_histogram(...)` with the configured interpolation mode and `renormalize_shifted_templates`.

The active farm configuration has:

```text
minimum_template_integral = 1.0
renormalize_shifted_templates = true
interpolation_mode = linear
```

The candidate scan currently rejects a candidate when:

```python
_hist_integral(template) < minimum_template_integral
```

That strict boundary is in scope for a machine-scale tolerance repair only.

### 5.4 Sequential scientific behavior

The resolver scans/applies pion components sequentially to the residual. Do not alter that ordering or residual construction. The observed `pi_delta` movement is a downstream consequence of the unstable `pi_n` acceptance, not authorization to redesign component fitting.

---

## 6. Allowed substantive source changes

Only these existing substantive files are allowed:

- `src/cuts/pion_component_fits.py`
- `testing/test_pion_component_dynamic_alignment.py`

A new narrowly scoped pure-Python test helper/module is allowed **only if absolutely necessary** to test the numeric boundary deterministically on a host without PyROOT. Prefer extending the existing focused test file. Do not introduce a new production architecture solely for testability.

If implementation requires another production source file, **STOP** and report the blocker rather than expanding scope.

---

## 7. Required repair A — stable pion-control persistence compatibility

### 7.1 Required behavior

Keep `hist_name` in serialized diagnostic/provenance output if useful. Do **not** require the generated ROOT histogram name to match for persistence reuse.

For compatibility, pion-control identity must be based on stable semantic information:

- content checksum;
- histogram axis specification.

Required behavior:

- same checksum + same axis + different `hist_name` => compatible;
- changed checksum => incompatible;
- changed axis => incompatible;
- missing/malformed required semantic checksum or axis => incompatible/fail closed.

Do not remove pion-control identity checking.

### 7.2 Preserve every other compatibility boundary

Do not weaken or bypass:

- alignment schema version;
- resolved configuration hash;
- analysis scope;
- complete physical-bin identity;
- parent alignment hash;
- immutable source-template identifier/checksum checks;
- setting/phi/epsilon store identity;
- any existing exact source/provenance checks outside the ephemeral pion-control ROOT name.

Do not special-case `Q4p4W2p74`, `Left`, `lowe`, the September 29 run, or a specific histogram name/checksum.

### 7.3 Schema policy

Do not bump the alignment schema merely to make the repair work. The existing serialized payload already contains checksum and axis, so a compatibility-level repair should remain backward compatible with valid existing schema-v2 stores.

If a schema bump proves technically unavoidable, **STOP** and report why before implementing it. A schema bump that simply invalidates all prior stores is not an acceptable substitute for the requested repair.

---

## 8. Required repair B — floating-point-safe minimum-template-integral boundary

### 8.1 Required behavior

Keep the configured scientific threshold exactly:

`minimum_template_integral = 1.0`

Do not lower it in configuration.

The candidate should be rejected only when the measured template integral is **materially below** the configured threshold, not when it differs from the threshold only by machine-scale floating-point roundoff.

Implement a narrow, explicit, documented floating-point comparison at **this minimum-template-integral boundary only**. A suitable pattern is equivalent to:

```python
value < threshold and not math.isclose(
    value,
    threshold,
    rel_tol=1e-12,
    abs_tol=1e-12,
)
```

The exact local helper form may differ, but acceptance must satisfy the tests below and must remain far tighter than any physically meaningful template normalization change.

Do not use a large epsilon, percentage tolerance, bin-width-dependent tolerance, or data-dependent tuning.

### 8.2 Required boundary behavior

At threshold `1.0`:

- `1.0` => pass;
- a value below `1.0` only at approximately machine-scale / `1e-13` level => pass;
- a clearly sub-threshold value such as `1.0 - 1e-6` or `0.99` => reject.

The repair must not change score ordering, candidate scoring, relative-improvement logic, data-support thresholds, evaluation-bin thresholds, boundary penalties, localization requirements, or lost-integral limits.

### 8.3 No forced result

Do not hard-code or force:

- `pi_n = -0.003 GeV`;
- `pi_delta = +0.008 GeV`;
- any particular window expansion;
- a particular parent score;
- a particular `w0`;
- cache reuse regardless of semantic provenance.

The repaired code must deterministically evaluate the existing scientific algorithm.

---

## 9. Frozen scientific/runtime interfaces

Unless a concrete blocker is found, the following are frozen and must remain unchanged:

- `src/main.py`
- `src/cuts/rand_sub.py`
- `src/cuts/pion_component_subtraction.py`
- `src/cuts/particle_subtraction.py`
- `src/cuts/full_background_subtraction_plots.py`
- `src/cuts/pion_hgcer_method_a_parallel_full_procedure.py`
- `src/binning/calculate_yield.py`
- `src/binning/ave_per_bin.py`
- `src/utility/background_config.py`
- `run_Prod_Analysis.sh`
- bundle collector/profile/wrapper files
- F.6.3 producer and E.8.4 consumer
- frozen F.6.2 artifacts/fingerprints.

Also frozen:

- component definitions;
- component ordering;
- common-setting shift;
- scan ranges and scan step sizes;
- fit/evaluation windows;
- window-expansion candidates;
- interpolation mode;
- shifted-template renormalization behavior;
- acceptance thresholds as configured;
- score metric and penalties;
- templates and priors;
- random/dummy and slow-proton subtraction;
- pion subtraction formulas and public baseline production branch;
- baseline-yield calculation;
- canonical binning;
- Method-A factors/correction ownership;
- Method B;
- SIMC normalization;
- efficiencies, cross sections, and L/T separation;
- active `no_empirical_residual` profile and dormant Fit 1/Fit 2.

Do not alter physics to make the checker or presentation agree.

---

## 10. Required tests

Extend the existing focused coverage in:

`testing/test_pion_component_dynamic_alignment.py`

### 10.1 Persistence positive test

Construct/store a valid alignment, then build the expected payload from a pion-control histogram with:

- the same binning/axis;
- the same bin contents/errors/checksum;
- a different ROOT histogram name.

The cache must load as compatible/reused without calling the resolver a second time.

This test must specifically prove that `hist_name` alone does not determine compatibility.

### 10.2 Persistence negative tests

Prove fail-closed rejection when:

1. pion-control content checksum changes while axis remains the same;
2. pion-control axis changes;
3. required semantic pion-control identity is malformed/missing, if that can occur through the helper boundary.

Retain the existing parent/config/source/axis/physical-identity rejection tests.

### 10.3 Numeric-boundary tests

Prove the actual acceptance predicate used by the candidate scan:

- exact threshold passes;
- machine-scale below-threshold value passes;
- materially below-threshold value rejects.

If PyROOT is available locally, include at least one actual scan/template regression that exercises the production path.

If PyROOT is unavailable, do not fake a runtime PASS. Add the narrowest deterministic pure-Python/source-contract coverage needed to execute the boundary logic locally while keeping the production implementation single-sourced. Do not duplicate the scientific algorithm into test code as an alternate implementation.

### 10.4 Regression tests

Preserve existing focused behaviors, including:

- compatible cache reuse;
- parent/config/source/axis/identity mismatch rejection;
- immutable raw-template single-shift provenance;
- parent-to-fine fallback behavior;
- independent boundary policies;
- support metric behavior;
- mixed accepted/fallback component behavior;
- no kaon-side data entering the alignment API.

The fix must not require changes to F.6.3/E.8.4 tests to make them pass.

---

## 11. Local deterministic validation

Use the repository's actual local Python interpreter. Do not claim ROOT/PyROOT coverage when ROOT is unavailable.

Run at minimum:

```bash
python -B -m py_compile \
  src/cuts/pion_component_fits.py \
  testing/test_pion_component_dynamic_alignment.py

python -B -m unittest testing.test_pion_component_dynamic_alignment -v

python -B -m unittest testing.test_t_bin_pion_parent_integrity -v
python -B -m unittest testing.test_binning_pre_particle_subtraction -v

python -B -m unittest testing.test_f6_3_parallel_full_procedure_method_a -v
python -B -m unittest testing.test_e8_4_production_impact_audit -v

python -B tools/update_memory_manifest.py --root . --write
python -B tools/update_memory_manifest.py --root . --check
python -B tools/check_memory_health.py --root .
python -B tools/memory_bootstrap.py --root . --json
python -B -m unittest testing.test_memory_health -v

git -c core.safecrlf=false diff --check
```

If the focused alignment suite reports PyROOT skips, report the exact skipped tests and do not convert that into runtime evidence.

Do not run the farm.

---

## 12. Warranted durable memory/history

Update only what this task establishes.

Allowed memory/history files:

- `docs/memory/CURRENT.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/manifest.json`
- `docs/memory/phases/e8-4-fix3-pion-alignment-determinism-task-contract.md`
- `docs/memory/phases/e8-4-fix3-pion-alignment-determinism.md`
- `docs/memory/evidence/e8-4-left-lowe-pion-alignment-determinism-blocker.md`

No `USER.md` change is needed.

Do not change `MEMORY.md`, `LEARNINGS.md`, `CODEX.md`, `TOOLS.md`, decisions, or handoff unless a concrete owned durable fact changes; if so, stop and explain before expanding scope.

### 12.1 Blocker evidence record

Create a concise evidence record for the fresh
`KaonLT_E8_4_Left_lowe_alignment_20260929-102321.json` farm artifact.

Record, at minimum:

- setting `Q4p4W2p74 / Left / lowe`;
- the identical pion-control checksum across the compared September 24 and September 29 setting-wide records;
- differing generated pion-control ROOT names;
- September 29 `rejected_stale_then_created` due to `pion control histogram identifier mismatch`;
- active `minimum_template_integral=1.0`,
  `renormalize_shifted_templates=true`,
  `interpolation_mode=linear`;
- the reciprocal `pi_n -0.003/-0.004` candidate acceptance flip;
- the downstream `pi_delta +0.008 -> +0.016` change;
- effectively identical baseline score/common-setting shift;
- conclusion: direct runtime evidence of a determinism blocker, **not** runtime closure.

Do not invent an artifact checksum that was not actually measured.

### 12.2 Current/roadmap status

`CURRENT.md` must make this exact fix the active work item and remove the stale pre-`9daca` commit/push NEXT.

After local implementation, status remains:

- E.8.4.Fix.3 — **ACTIVE** pending independent ChatGPT actual-diff review;
- E.8.4 — **SOURCE REVIEWED**;
- F.6.3 — **SOURCE REVIEWED**;
- final E.8 — **BLOCKED**;
- F.6.4 — **BLOCKED**.

Do not claim runtime validation from local tests.

Regenerate `docs/memory/manifest.json` after all versioned memory changes.

---

## 13. Forbidden shortcuts

Do not:

- delete the persisted alignment store as the "fix";
- disable persistence;
- disable pion-control identity checking;
- accept checksum or axis mismatches;
- ignore parent/config/template/source incompatibility;
- bump schema solely to force a fresh store;
- lower `minimum_template_integral`;
- disable `renormalize_shifted_templates`;
- change interpolation mode;
- widen/narrow scan ranges;
- change scan steps or expansion candidates;
- change scoring/penalties;
- hard-code the September 24 alignment values;
- hard-code a `w0`;
- reorder pion components;
- alter pion subtraction or public yields;
- touch Method A/B scientific ownership;
- alter F.6.3 or E.8.4 to hide the blocker;
- modify the validation profile in this task;
- add a fallback that silently treats malformed provenance as compatible;
- perform unrelated cleanup/refactors;
- commit, push, run the farm, or package farm artifacts.

---

## 14. Required actual-diff audit

Before stopping, inspect the real worktree:

```bash
git status --short
git diff -- src/cuts/pion_component_fits.py
git diff -- testing/test_pion_component_dynamic_alignment.py
git diff -- docs/memory/CURRENT.md
git diff -- docs/memory/roadmap/STATUS.md
git diff -- docs/memory/manifest.json
git diff -- docs/memory/phases/e8-4-fix3-pion-alignment-determinism-task-contract.md
git diff -- docs/memory/phases/e8-4-fix3-pion-alignment-determinism.md
git diff -- docs/memory/evidence/e8-4-left-lowe-pion-alignment-determinism-blocker.md
git -c core.safecrlf=false diff --check
```

Confirm there are no substantive changes outside the allowlist.

---

## 15. Review bundle

Create one complete temporary review bundle in the repository root, for example:

`kaonlt_review.diff`

It must contain:

1. starting branch and HEAD;
2. current HEAD;
3. `git status --short`;
4. `git diff --stat`;
5. the complete tracked diff;
6. complete `git diff --no-index /dev/null <file>` representations for every intended new/untracked repository file, including this task contract, phase record, and blocker-evidence record as applicable;
7. the exact commands and results of local deterministic checks;
8. a final changed-path inventory.

The review bundle itself is temporary and must **not** be staged or committed.

Do not rely on a prose Codex summary as proof.

---

## 16. Acceptance criteria

This task is ready for independent source review only when all are true:

1. implementation started from `test` HEAD
   `9daca79051a4eb451a158dd7c34b7028c6b96456`;
2. same pion-control checksum+axis with a different generated ROOT name is persistence-compatible;
3. checksum mismatch still rejects;
4. axis mismatch still rejects;
5. malformed/missing semantic pion-control identity fails closed;
6. every non-pion-control compatibility boundary remains unchanged;
7. configured `minimum_template_integral=1.0` remains unchanged;
8. machine-scale roundoff below that threshold does not reject an otherwise valid renormalized template;
9. materially sub-threshold template integrals still reject;
10. no scan range/window/step/scoring/interpolation/renormalization/component-order change occurred;
11. no baseline production, F.6.3, E.8.4, Method-A, Method-B, yield, or cross-section logic changed;
12. focused and regression tests were run locally with skips reported honestly;
13. blocker evidence and active-state memory are updated without claiming runtime closure;
14. memory manifest is regenerated and checks pass, apart from any already-known explicitly non-fatal health warning that must be reported verbatim;
15. actual diff is confined to the allowlist;
16. complete review bundle is produced;
17. no commit, push, farm run, farm packaging, or later-gate action was performed.

---

## 17. Hard stop / NEXT

**STOP after producing the implementation, warranted memory/evidence updates, local checks, actual-diff audit, and complete review bundle.**

Do not commit or push.

Do not re-pin the validation profile in this task.

Do not provide or authorize a farm-run command.

**NEXT — independent ChatGPT review of the complete E.8.4.Fix.3 actual diff/review bundle only.**
