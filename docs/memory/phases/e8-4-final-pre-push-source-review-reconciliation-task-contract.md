# KaonLT E.8.4 — Final Pre-Push Source-Review Memory Reconciliation

## Purpose

Reconcile durable repository memory after independent ChatGPT review of the
complete local cumulative E.8.4 + Fix.1 + Fix.2 candidate.

Independent review of:

```text
kaonlt_review(20260928-233040).diff
```

**passed** for source/runtime-path ownership.

This task is **memory/history only**. It must not change E.8.4 source, tests,
F.6.3, production physics, or runtime logic.

The review establishes **SOURCE REVIEWED only** for the local E.8.4 candidate.
It does not establish ROOT/PyROOT rendering, full `main.py`, procedure-PDF
runtime behavior, Jefferson Lab farm integration, or production validation.

---

# 1. Exact starting state

Required branch:

```text
test
```

Required committed HEAD:

```text
d84f9610427cf38ccaaa8d6d8b14a8a995ecbac9
Add F6.3 parallel Method-A full procedure
```

The worktree must contain the already-reviewed cumulative E.8.4 + Fix.1 + Fix.2
candidate corresponding exactly to:

```text
kaonlt_review(20260928-233040).diff
```

plus this reconciliation contract after the user copies it into repository
memory.

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git log -1 --oneline
git status --short
```

Hard stop if:

- committed HEAD differs from `d84f9610427cf38ccaaa8d6d8b14a8a995ecbac9`;
- `src/cuts/full_background_subtraction_plots.py` or
  `testing/test_e8_4_production_impact_audit.py` differs from the reviewed
  `kaonlt_review(20260928-233040).diff`;
- any frozen producer/runtime file changed;
- unrelated tracked or untracked work is present.

Do not reset, stash, clean, checkout over, commit, push, or run the farm.

Local-only `AGENTS.md` / `.codex/` may remain untracked. Temporary
`kaonlt_review(...).diff` bundles may remain untracked and must not be staged.

---

# 2. Required startup reading

Read in this exact order:

1. local root `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read only:

- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/phase-f6-method-a-production-promotion.md`
- `docs/memory/phases/e8-2-baseline-stage-audit.md`
- `docs/memory/phases/e8-3-detached-method-a-reweighting-audit.md`
- `docs/memory/phases/f6-3-parallel-full-procedure-method-a.md`
- `docs/memory/phases/e8-4-production-impact-audit-task-contract.md`
- `docs/memory/phases/e8-4-fix1-source-review-repair-task-contract.md`
- `docs/memory/phases/e8-4-fix2-malformed-sidecar-fail-closed-task-contract.md`
- `docs/memory/phases/e8-4-production-impact-audit.md`
- this reconciliation contract

Do not reopen settled scientific architecture.

---

# 3. Independent source-review verdict to record

Record that independent ChatGPT actual-diff/source-runtime-path review of the
complete cumulative candidate:

```text
kaonlt_review(20260928-233040).diff
```

**passed**.

The review established, at source level:

- E.8.4 remains a presentation-only consumer of the already-produced
  F.6.3 `_f6_3_parallel_method_a_source`;
- the public baseline branch is unchanged;
- F.6.3 remains the sole owner of the private parallel Method-A branch;
- Method B has no numerical input;
- empirical residual Fit 1 / Fit 2 remain dormant under
  `no_empirical_residual`;
- semantic E.8.2 epsilon `low/high` is explicitly mapped only for comparison
  with F.6.3 `lowe/highe`;
- E.8.2 wide diagnostic-MM geometry remains separate from F.6.3 analysis-MM
  geometry;
- F.6.3 analysis-MM histograms are internally fail-closed for geometry;
- the common pion input, `B_pi_0`, and `B_pi_A` remain simultaneously visible;
- malformed optional F.6.3 Lambda-window and child-inventory containers fail
  closed into explicit E.8.4-unavailable payloads;
- an unavailable or malformed optional E.8.4 branch does not invalidate an
  otherwise successful baseline finalization;
- retained histograms are detached;
- `DeltaY` and defined `DeltaY/Y0` are display comparisons from stored
  producer-owned yields only;
- no uncertainty is invented for the Method-A shift;
- no child renormalization, fit, factor reconstruction, tree traversal,
  template refill, or second yield calculation was added;
- the existing single post-yield pair-safe PDF/page-manifest transaction is
  preserved;
- E.8.4 pages remain after E.8.3 and before the terminal E.8 handoff;
- no production promotion is performed.

Also record explicitly:

```text
Codex-reported local deterministic tests were NOT RUN by ChatGPT.
```

The review does not claim ROOT/PyROOT, farm, rendered-PDF, full-analysis, or
runtime acceptance.

---

# 4. Required status reconciliation

Update durable memory to:

```text
F.6.2 / F.6.2.Fix.5 — CLOSED / RUNTIME VALIDATED
E.8                 — ACTIVE
E.8.2               — SOURCE REVIEWED
E.8.3               — SOURCE REVIEWED
F.6.3               — SOURCE REVIEWED
E.8.4               — SOURCE REVIEWED
final E.8            — BLOCKED
F.6.4               — BLOCKED
```

Do not mark E.8.4:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

or:

```text
CLOSED / RUNTIME VALIDATED
```

because no E.8.4 farm/runtime evidence exists yet.

Fix.1 and Fix.2 are incorporated into the source-reviewed E.8.4 cumulative
candidate; do not present them as still ACTIVE after this reconciliation.

---

# 5. Exact NEXT

The exact next workflow is:

```text
user-controlled commit/push of the reviewed cumulative E.8.4 + Fix.1 + Fix.2
+ final source-review memory reconciliation
-> ChatGPT pushed-state review
-> one narrow Jefferson Lab farm gate for Q4p4W2p74 / Left / lowe
-> fresh artifacts
-> ChatGPT evidence review
```

Do not run the farm before the user-controlled commit/push and ChatGPT
pushed-state review.

The narrow farm gate is the governing post-source-review milestone. Do not
broaden immediately to the canonical-five-setting campaign.

Final E.8 closure remains blocked pending the required runtime/visual evidence.
F.6.4 remains blocked pending final production-impact evidence and an explicit
promotion decision.

---

# 6. Allowed changes

This reconciliation is memory/history only.

Allowed:

```text
docs/memory/CURRENT.md
docs/memory/decisions/e8-full-analysis-procedure-roadmap.md
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/phases/e8-4-production-impact-audit.md
docs/memory/phases/e8-4-final-pre-push-source-review-reconciliation-task-contract.md
docs/memory/roadmap/STATUS.md
docs/memory/manifest.json
```

Do not edit the original E.8.4, Fix.1, or Fix.2 task contracts merely to update
their historical status language. They are historical records.

Do not modify `MEMORY.md`, `USER.md`, `CURRENT_HANDOFF.md`, `LEARNINGS.md`,
`TOOLS.md`, or `CODEX.md` unless a genuinely new fact owned by one of those
records emerged. The independent review itself does not require such a change.

---

# 7. Explicitly frozen source/test files

These must remain byte-identical to the independently reviewed cumulative
candidate:

```text
src/cuts/full_background_subtraction_plots.py
testing/test_e8_4_production_impact_audit.py
```

These pushed producer/runtime files must also remain byte-unchanged:

```text
src/main.py
src/binning/calculate_yield.py
src/cuts/rand_sub.py
src/cuts/pion_component_subtraction.py
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
src/utility/background_config.py
```

Do not modify any Method-B file, F-stage calculator/validator, launcher,
collector, wrapper, yield/cross-section implementation, or production
subtraction code.

---

# 8. Phase-record reconciliation

In:

```text
docs/memory/phases/e8-4-production-impact-audit.md
```

append a concise final source-review reconciliation section recording:

- reviewed artifact:
  `kaonlt_review(20260928-233040).diff`;
- independent source/runtime-path verdict: PASS;
- E.8.4 cumulative status: `SOURCE REVIEWED`;
- Fix.1 and Fix.2 findings are resolved in the reviewed cumulative candidate;
- frozen source ownership was preserved;
- Codex-reported checks were `NOT RUN by ChatGPT`;
- no farm/ROOT/PyROOT/full-runtime claim;
- exact next workflow is commit/push -> pushed-state review -> narrow
  `Q4p4W2p74 / Left / lowe` farm gate.

Do not erase the historical ACTIVE sections for the implementation and fixes;
append the reconciliation so the record preserves chronology.

---

# 9. Memory integrity checks

After memory edits:

```bash
python -B tools/update_memory_manifest.py --root . --write
python -B tools/update_memory_manifest.py --root . --check
python -B tools/check_memory_health.py --root .
python -B tools/memory_bootstrap.py --root . --json
python -B -m unittest testing.test_memory_health -v
```

The bootstrap output must still report committed HEAD
`d84f9610427cf38ccaaa8d6d8b14a8a995ecbac9` and only the intended cumulative
E.8.4/reconciliation paths as dirty.

Because source/test code must not change in this reconciliation, do not rerun
or alter scientific/source tests merely to manufacture new evidence. Existing
Codex-reported E.8.4 test results remain historical local evidence.

---

# 10. Diff audit

Inspect the actual cumulative diff:

```bash
git status --short
git diff --stat
git diff -- docs/memory
git diff -- src/cuts/full_background_subtraction_plots.py
git diff -- testing/test_e8_4_production_impact_audit.py
git -c core.safecrlf=false diff --check
```

Verify:

- source/test portions are identical to
  `kaonlt_review(20260928-233040).diff`;
- only warranted memory/history changed after the reviewed candidate;
- all frozen files remain unchanged;
- no temporary review bundle is staged.

Create a fresh complete timestamped review bundle if needed for final
pre-commit inspection. It must include the complete cumulative diff and all
intended untracked/new files via `git diff --no-index /dev/null ...` as needed.

Do not stage merely to create the review bundle.

---

# 11. Acceptance criteria

This reconciliation is complete only if:

1. branch remains `test`;
2. committed HEAD remains exactly `d84f961...`;
3. the independently reviewed E.8.4 source/test candidate is unchanged;
4. E.8.4 is recorded as `SOURCE REVIEWED`;
5. E.8 remains `ACTIVE`;
6. E.8.2, E.8.3, and F.6.3 remain `SOURCE REVIEWED`;
7. final E.8 and F.6.4 remain `BLOCKED`;
8. memory names `kaonlt_review(20260928-233040).diff` as the independent review
   artifact;
9. local Codex tests remain clearly distinguished from ChatGPT source review
   and farm/runtime evidence;
10. exact NEXT is user commit/push -> ChatGPT pushed-state review -> narrow
    `Q4p4W2p74 / Left / lowe` farm gate;
11. manifest and memory-health checks pass;
12. `git diff --check` passes;
13. no source/test/physics change was made;
14. Codex did not commit, push, or run the farm.

---

# 12. Hard stop

Stop after:

- memory/status reconciliation;
- manifest regeneration/check;
- memory-health checks;
- actual cumulative diff audit;
- optional refreshed review bundle.

Do not commit.
Do not push.
Do not run the farm.
Do not begin final E.8 closure.
Do not begin F.6.4.
Do not promote Method A.

The user alone controls commit/push.
