# KaonLT F.6.3 Fix.2 — Public Regression Contract Closure

## Purpose

Close one remaining **test/contract coverage gap only** in the current local
F.6.3 candidate.

Independent review of the refreshed cumulative candidate found that the Fix.1
source repairs are correct:

- the private F.6.3 side branch is fail-closed for ordinary branch-local
  exceptions;
- live F.1 parity uses the exact effective
  `_component_cache_event_coefficient(...)`;
- parity uses the accepted F.4 scaled tolerance;
- public `Y0` and private `YA` share the unchanged final yield/error helper;
- factor-value-derived retained provenance is removed.

Do **not** change those source implementations.

The remaining issue is that the mandatory public
`calculate_yield_data(...)` regression does not yet prove every output required
by the original F.6.3 contract and Fix.1 contract. In particular, it reduces
the public `groups` result to a numeric tuple and does not snapshot/compare the
processed component payload. The repository contract explicitly requires the
public return **structure and values** plus component payloads to remain
unchanged.

This is a **test-only narrow repair** unless a concrete blocker is found.

---

## Exact committed base and current worktree

Authoritative committed base:

```text
branch: test
HEAD: c9ed0d6b4d0013f7475eedbfca57a096730ea840
subject: Add E8.3 detached Method-A reweighting audit
```

The worktree already contains the cumulative F.6.3 + Fix.1 candidate reviewed
in:

```text
kaonlt_review(20260924-193936).diff
```

Do not reset, stash, clean, checkout over, or discard the current candidate.

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
git log -1 --oneline
```

Hard stop if committed HEAD differs from
`c9ed0d6b4d0013f7475eedbfca57a096730ea840` or unrelated dirty work is present.

---

## Required startup reading

Read, in order:

1. local root `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read:

- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/phase-f6-method-a-production-promotion.md`
- `docs/memory/phases/e8-2-baseline-stage-audit.md`
- `docs/memory/phases/e8-3-detached-method-a-reweighting-audit.md`
- `docs/memory/phases/f6-3-parallel-full-procedure-method-a-task-contract.md`
- `docs/memory/phases/f6-3-parallel-full-procedure-method-a-fix1-task-contract.md`
- `docs/memory/phases/f6-3-parallel-full-procedure-method-a.md`
- this Fix.2 contract

Status remains:

```text
E.8   — ACTIVE
E.8.2 — SOURCE REVIEWED
E.8.3 — SOURCE REVIEWED
F.6.3 — ACTIVE
E.8.4 — BLOCKED
```

Do not mark F.6.3 `SOURCE REVIEWED`.

---

# Independent-review result

The source changes themselves pass this review. Do not modify:

```text
src/binning/calculate_yield.py
src/cuts/pion_component_subtraction.py
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
```

unless the required test exposes a concrete source defect. If that occurs,
stop and report the blocker rather than silently broadening this task.

The only remaining acceptance gap is deterministic public regression coverage.

The original F.6.3 contract requires:

```text
existing public calculate_yield_data(...) return structure and baseline
values/errors are unchanged when the Method-A side branch is unavailable

public Y0, E.8.2 source, stage-window yields, scale factors, component
payloads and baseline ROOT objects are unchanged
```

Fix.1 additionally requires the public regression to snapshot the historical
baseline/public result and prove that the returned `groups` **structure and
values** and component payloads remain unchanged.

The current test instead stores:

```python
"groups": (
    float(groups[(0, 0)]["kaon"]),
    float(groups[(0, 0)]["kaon_err"]),
)
```

and does not snapshot or compare a component payload.

That is insufficient for the explicit regression contract even though the
source path itself is currently correct.

---

# Required repair

Modify only:

```text
testing/test_f6_3_parallel_full_procedure_method_a.py
```

plus the warranted memory/history files listed below.

## 1. Exercise the actual public yield-mode return shape

The production caller uses:

```python
calculate_yield_data("yield", ...)
```

The F.6.3 public regression should exercise that public mode rather than
reducing the return to a generic `"kaon"` tuple.

Using the existing detached `_public_yield_fixture(...)`, adapt the local
F.6.3 test fixture without changing the shared E.8.2 fixture file. It is
acceptable to alias/move the one fixture payload from the `"kaon"` key to the
`"yield"` key inside this F.6.3 test.

Call:

```python
calculate_yield_data("yield", ...)
```

and retain a serializable snapshot of the complete returned structure.

For the one-child fixture, assert the complete public shape and values,
including:

```text
top-level child-key inventory
per-child field-key inventory
yield value
yield_err value
```

Do not merely extract two floats and compare those.

Preserve the existing expected historical numerical result.

## 2. Add a representative component-payload sentinel

Before the public call, attach a representative nested serializable component
payload to the processed child, for example under the existing
`particle_subtraction_component_payload` key.

The payload should be intentionally nontrivial enough to detect mutation:
nested dict/list/scalar content is sufficient. It need not be a ROOT object
because this regression is checking public state mutation, not ROOT behavior.

Before calling `calculate_yield_data(...)`:

- take a deep-value snapshot;
- retain the original object identity where useful.

After both:

1. ordinary unavailable F.6.3; and
2. injected branch-local `RuntimeError`;

assert:

- payload contents are byte/value-equivalent to the pre-call snapshot;
- the payload was not replaced;
- no child field inside it was added, removed, or modified.

This is a regression for preservation only. Do not make the unavailable branch
consume the payload.

## 3. Preserve the existing public checks

Retain the current public assertions for:

- historical public yield/error value;
- unavailable vs branch-runtime-failure equality;
- E.8.2 final yield/statistical/total values;
- `stage_window_yields`;
- scale factor;
- final histogram contents;
- final histogram errors;
- F.6.3 unavailable state and precise runtime-failure reason.

The regression should now cover the complete contract:

```text
public groups structure + values
public Y0 / total error
E.8.2 source values
stage_window_yields
scale factors
component payloads
baseline histogram contents/errors
```

## 4. Do not change the accepted private all-one test

The existing private branch test already proves:

```text
B_pi^A == B_pi^0
MM_A == MM_0
YA == Y0
YA_statistical_error == Y0_statistical_error
YA_total_error == Y0_total_error
```

for all-one Method A.

Do not redesign that fixture in this task.

---

# Allowed files

Tests:

```text
testing/test_f6_3_parallel_full_procedure_method_a.py
```

Memory/history:

```text
docs/memory/phases/f6-3-parallel-full-procedure-method-a-fix2-task-contract.md
docs/memory/phases/f6-3-parallel-full-procedure-method-a.md
docs/memory/manifest.json
```

No other file should change.

If CURRENT/roadmap/status genuinely require a correction because this test
exposes a contradiction, stop and report rather than changing them
preemptively. Their present status is already correct.

---

# Frozen files

Do not edit any source file, including:

```text
src/binning/calculate_yield.py
src/cuts/pion_component_subtraction.py
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
src/cuts/rand_sub.py
src/cuts/full_background_subtraction_plots.py
src/utility/background_config.py
src/utility/root_histogram_ownership.py
```

Do not edit accepted F.1-F.6.2 scientific modules/artifacts.

Do not:

- change cuts;
- change fits/windows/templates/priors;
- change normalization;
- change `w0`;
- change `C_j`;
- change proton subtraction;
- change pruning;
- change canonical binning;
- activate Method B;
- activate empirical residual Fits 1/2;
- implement E.8.4;
- commit/push;
- run the farm.

---

# Deterministic validation

Use the repository's actual Python command/environment established by
`AGENTS.md`.

At minimum:

```bash
<PYTHON> -m py_compile \
  testing/test_f6_3_parallel_full_procedure_method_a.py

<PYTHON> -m unittest testing.test_f6_3_parallel_full_procedure_method_a -v
<PYTHON> -m unittest testing.test_e8_2_baseline_stage_audit -v
<PYTHON> -m unittest testing.test_e8_3_detached_method_a_reweighting_audit -v
<PYTHON> -m unittest testing.test_pion_hgcer_method_a_parent_preserving_correction -v
<PYTHON> -m unittest testing.test_pion_hgcer_method_a_tphi_propagation -v
<PYTHON> -m unittest testing.test_pion_hgcer_phase_f_runtime_contract -v
<PYTHON> -m unittest testing.test_pion_component_dynamic_alignment -v

<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/check_memory_health.py --root .
<PYTHON> -B tools/memory_bootstrap.py --root . --json
<PYTHON> -B -m unittest testing.test_memory_health -v

git -c core.safecrlf=false diff --check
```

Do not claim ROOT/PyROOT/full-analysis/farm/runtime validation.

---

# Memory/history

Append a short Fix.2 subsection to:

```text
docs/memory/phases/f6-3-parallel-full-procedure-method-a.md
```

Record only:

- independent source review found the Fix.1 source repairs correct;
- one deterministic public-regression coverage gap remained;
- Fix.2 closed the exact public `groups` structure/value and component-payload
  preservation assertions;
- exact tests Codex actually ran and their results;
- `NOT farm validated`;
- F.6.3 remains `ACTIVE`;
- NEXT = independent ChatGPT review of the refreshed cumulative diff.

Do not rewrite prior F.6.3/Fix.1 history.

Regenerate `docs/memory/manifest.json`.

---

# Diff audit

Before stopping:

```bash
git status --short
git diff --name-only
git diff --stat
git -c core.safecrlf=false diff --check
```

The new Fix.2 delta itself should touch only:

```text
testing/test_f6_3_parallel_full_procedure_method_a.py
docs/memory/phases/f6-3-parallel-full-procedure-method-a-fix2-task-contract.md
docs/memory/phases/f6-3-parallel-full-procedure-method-a.md
docs/memory/manifest.json
```

Refresh root:

```text
kaonlt_review.diff
```

as the complete cumulative diff from committed HEAD:

```text
c9ed0d6b4d0013f7475eedbfca57a096730ea840
```

Include every intended untracked file with `git diff --no-index /dev/null ...`
as required.

Do not stage or commit `kaonlt_review.diff`.

---

# Acceptance criteria

Return for independent review only when all are true:

1. committed HEAD remains exactly
   `c9ed0d6b4d0013f7475eedbfca57a096730ea840`;
2. no source file changed in Fix.2;
3. the public regression invokes `calculate_yield_data("yield", ...)`;
4. the complete one-child `groups` structure and values are asserted;
5. unavailable and runtime-failed F.6.3 produce identical public structures and
   values;
6. a representative component payload is proven content-identical and
   unreplaced across both public calls;
7. E.8.2 values, stage-window yields, scale factor, baseline histogram
   contents/errors remain asserted unchanged;
8. the existing all-one private branch test remains intact;
9. deterministic tests pass;
10. memory remains `F.6.3 — ACTIVE`;
11. E.8.4 remains `BLOCKED`;
12. no ROOT/PyROOT/farm/runtime claim is made.

---

# Hard stop

Stop and report a blocker rather than changing source if the strengthened
public regression exposes an actual F.6.3 source mutation or runtime-path
defect.

Do not broaden scope, commit, push, or run the farm.
