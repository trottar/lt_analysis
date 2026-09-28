# KaonLT F.6.3 — Final Pre-Push Source-Review Memory Reconciliation

## Purpose

Reconcile durable repository memory after independent ChatGPT review of the
complete local F.6.3 + Fix.1 + Fix.2 candidate.

This task is **memory/history only**.

Independent ChatGPT actual-diff review of:

```text
kaonlt_review(20260928-161516).diff
```

passed.

The accepted local candidate is based on committed:

```text
branch: test
HEAD: c9ed0d6b4d0013f7475eedbfca57a096730ea840
subject: Add E8.3 detached Method-A reweighting audit
```

The review established:

- the F.6.3 scientific/source architecture passes;
- Fix.1 source repairs pass;
- Fix.2 closes the remaining public-regression contract gap;
- Fix.2 made no substantive source change;
- the public regression now exercises `calculate_yield_data("yield", ...)`,
  freezes the complete one-child public return structure/values, and proves
  processed component-payload identity/content preservation for both an
  unavailable F.6.3 branch and an injected branch-local runtime failure;
- the private all-one Method-A regression remains intact;
- no Method-B numerical path, empirical residual Fit 1/Fit 2 activation,
  child renormalization, second tree traversal, or production promotion was
  introduced.

This review establishes **SOURCE REVIEWED only**. It does not establish
ROOT/PyROOT, full-analysis, farm, rendering, or runtime validation.

The next implementation phase is E.8.4, but do not implement E.8.4 in this
task.

---

## Exact starting state

Committed branch/HEAD must remain:

```text
test
c9ed0d6b4d0013f7475eedbfca57a096730ea840
```

The worktree must contain the already-reviewed cumulative F.6.3 candidate.

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
git log -1 --oneline
```

Hard stop if:

- committed HEAD differs;
- the cumulative F.6.3 source/test candidate has changed since
  `kaonlt_review(20260928-161516).diff`;
- unrelated dirty work is present.

Do not reset, stash, clean, checkout over, commit, push, or run the farm.

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
- `docs/memory/phases/f6-3-parallel-full-procedure-method-a.md`
- the original F.6.3 task contract;
- the F.6.3 Fix.1 contract;
- the F.6.3 Fix.2 contract;
- this reconciliation contract.

---

# Required status reconciliation

Update durable memory to the post-review state:

```text
E.8   — ACTIVE
E.8.2 — SOURCE REVIEWED
E.8.3 — SOURCE REVIEWED
F.6.3 — SOURCE REVIEWED
E.8.4 — NEXT
final E.8 — BLOCKED pending E.8.4 and later runtime/visual validation
F.6.4 — BLOCKED pending completed F.6.3/E.8.4 production-impact evidence
```

Do not claim F.6.3 runtime validation.

Do not mark E.8.4 ACTIVE; implementation has not started.

Do not alter the existing E.8.1/F.6.2 closures or deferred canonical-five
status.

---

## CURRENT.md

Update only the F.6.3/E.8.4 active-state text.

Record that:

- independent ChatGPT actual-diff/source-runtime-path review of the cumulative
  F.6.3 + Fix.1 + Fix.2 candidate passed;
- review artifact:
  `kaonlt_review(20260928-161516).diff`;
- Codex-reported deterministic checks were not run independently by ChatGPT;
- F.6.3 is `SOURCE REVIEWED`, not ROOT/PyROOT/farm/runtime validated;
- E.8.4 is `NEXT`;
- no farm run occurs between F.6.3 and E.8.4 source development under the
  governing roadmap;
- immediate NEXT is user-controlled commit/push of the reviewed cumulative
  F.6.3 set -> ChatGPT pushed-state review -> E.8.4 source/runtime-path audit
  and standalone implementation contract.

Remove stale statements that F.6.3 is `ACTIVE` pending independent review or
that E.8.4 remains blocked on that review.

---

## E.8 roadmap and STATUS

In:

```text
docs/memory/decisions/e8-full-analysis-procedure-roadmap.md
docs/memory/roadmap/STATUS.md
```

reconcile only the status text:

```text
F.6.3 — SOURCE REVIEWED
E.8.4 — NEXT
```

Preserve the governing architecture verbatim in substance:

- baseline production remains unchanged;
- Method-A branch changes only `w0_j -> w0_j * C_j`;
- same proton-cleaned input;
- no Method B numerical dependency;
- no empirical residual Fit 1/Fit 2;
- E.8.4 consumes F.6.3 branch outputs and never constructs them;
- final E.8 and F.6.4 remain downstream;
- farm milestone remains after coherent F.6.3 + E.8.4 source review.

Do not rewrite the roadmap.

---

## Phase-F record

In:

```text
docs/memory/phases/phase-f6-method-a-production-promotion.md
```

reconcile F.6.3 from `ACTIVE` to `SOURCE REVIEWED` and E.8.4 from `BLOCKED`
to `NEXT`.

Preserve all scientific ownership and promotion boundaries.

F.6.4 remains blocked.

---

## F.6.3 phase record

In:

```text
docs/memory/phases/f6-3-parallel-full-procedure-method-a.md
```

change the top-level status from `ACTIVE` to `SOURCE REVIEWED`.

Append a concise independent-review closure subsection recording:

- review artifact:
  `kaonlt_review(20260928-161516).diff`;
- remote committed base remained
  `c9ed0d6b4d0013f7475eedbfca57a096730ea840`;
- Fix.1 source implementation passed independent review;
- comparison against the prior candidate showed no substantive source change
  in Fix.2;
- Fix.2 changed only the focused test plus warranted memory/contract/manifest
  state;
- Fix.2 public regression now verifies:
  - actual `"yield"` mode;
  - complete public one-child `groups` structure and values;
  - unavailable/runtime-failed branch equality;
  - E.8.2 final yield/statistical/total values;
  - stage-window yields;
  - scale factor;
  - baseline histogram contents/errors;
  - processed component-payload object identity and deep/serialized value
    preservation;
- private all-one Method-A baseline reproduction remains covered;
- Codex deterministic checks are recorded as Codex-reported and
  **NOT RUN by ChatGPT**;
- no ROOT/PyROOT/farm/runtime acceptance claim;
- NEXT = user commit/push -> pushed-state review -> E.8.4 audit/contract.

Do not rewrite the implementation chronology.

---

# Allowed files

```text
docs/memory/CURRENT.md
docs/memory/decisions/e8-full-analysis-procedure-roadmap.md
docs/memory/roadmap/STATUS.md
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/phases/f6-3-parallel-full-procedure-method-a.md
docs/memory/phases/f6-3-final-pre-push-source-review-reconciliation-task-contract.md
docs/memory/manifest.json
```

No source file or test file may change.

Do not modify the original F.6.3, Fix.1, or Fix.2 contracts.

---

# Frozen files

All analysis and test code is frozen in this task, including:

```text
src/binning/calculate_yield.py
src/cuts/pion_component_subtraction.py
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
testing/test_f6_3_parallel_full_procedure_method_a.py
```

Also freeze:

```text
src/main.py
src/cuts/rand_sub.py
src/cuts/full_background_subtraction_plots.py
src/utility/background_config.py
```

and all accepted F.1-F.6.2 science/artifacts.

Do not implement E.8.4.

---

# Validation

Because this is memory-only reconciliation, do not rerun the entire source test
suite unless repository policy requires it.

Run at minimum:

```bash
<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/check_memory_health.py --root .
<PYTHON> -B tools/memory_bootstrap.py --root . --json
<PYTHON> -B -m unittest testing.test_memory_health -v

git -c core.safecrlf=false diff --check
```

Also verify that the reviewed source/test sections have not changed:

```bash
git diff -- src/binning/calculate_yield.py
git diff -- src/cuts/pion_component_subtraction.py
git diff -- src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
git diff -- testing/test_f6_3_parallel_full_procedure_method_a.py
```

Those diffs should remain exactly the already-reviewed cumulative F.6.3
candidate, with no new reconciliation delta.

---

# Diff audit

Before stopping:

```bash
git status --short
git diff --name-only
git diff --stat
git -c core.safecrlf=false diff --check
```

Confirm the reconciliation delta itself touches only the allowlisted memory
paths.

Refresh the root review bundle:

```text
kaonlt_review.diff
```

as the complete cumulative diff from committed:

```text
c9ed0d6b4d0013f7475eedbfca57a096730ea840
```

including every intended untracked file with complete
`git diff --no-index /dev/null ...` representations.

Do not stage or commit `kaonlt_review.diff`.

---

# Acceptance criteria

Return for final ChatGPT pre-push review only if:

1. committed HEAD is still
   `c9ed0d6b4d0013f7475eedbfca57a096730ea840`;
2. no source/test implementation changed;
3. F.6.3 is consistently `SOURCE REVIEWED`;
4. E.8.4 is consistently `NEXT`;
5. final E.8 remains blocked pending E.8.4 and later runtime/visual validation;
6. F.6.4 remains blocked pending production-impact evidence;
7. CURRENT NEXT is:
   user commit/push -> pushed-state review -> E.8.4 audit/contract;
8. farm cadence still defers the next farm gate until coherent F.6.3 + E.8.4
   source review;
9. no ROOT/PyROOT/farm/runtime claim is added;
10. memory manifest/health checks pass;
11. cumulative review bundle is refreshed;
12. no commit or push is performed.

---

# Hard stop

Stop and report rather than altering source/test code, scientific ownership,
accepted artifacts, runtime behavior, or E.8.4 implementation.

Do not run the farm.

Do not commit or push.
