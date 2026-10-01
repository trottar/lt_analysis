# KaonLT Codex workflow

This record owns source-changing Codex workflow. Codex may inspect live source
and targeted repository memory, work one allowlisted phase or fix contract at a
time, run deterministic local checks, trace relevant source/runtime paths,
audit the actual diff, and update source-review memory when the contract allows.

Codex must not commit, push, update remote refs, initiate Jefferson Lab farm
execution, or claim farm execution, ROOT/PyROOT integration, full `main.py`
validation, rendered procedure-PDF farm behavior, or phase closure without
direct applicable farm evidence. Passing local tests is not farm integration.
ChatGPT audits the actual diff. The user alone commits/pushes accepted changes and runs farm validation.
Workflow: Codex local changes -> ChatGPT audit -> user commit/push -> ChatGPT pushed-state synchronization review -> user farm run -> ChatGPT evidence review.
Default cycle: contract -> Codex implementation + local tests -> ChatGPT
actual-diff review -> user commit/push -> ChatGPT pushed-state synchronization
review -> user narrow farm gate if required -> fresh evidence review.
Opening and closing memory checkpoints bracket this substantive cycle; see
the [workflow decision](decisions/memory-bracketed-scientific-throughput.md).
Commit/push and pushed-state review are required synchronization stages,
not independent scientific phases.

## Task classes

1. **Implementation task:** make the approved source change, run local checks,
   trace the preserved path, audit the diff, and update source-review memory
   only when allowlisted.
2. **Closure task:** reconcile documentation, evidence, and status only after
   accepted farm evidence; it must not change scientific source.

## Actual-diff review before user handoff

1. Codex implements locally and stops before commit/push.
2. ChatGPT reviews the actual diff, not the Codex summary.
3. For a small diff, terminal output may be pasted directly.
4. For a large diff, create one review bundle directly in the repository root
   with a clearly temporary name such as `kaonlt_review.diff`, and have the
   user upload it to ChatGPT. Remove that review file before commit/push unless
   it is explicitly intended to be tracked.
5. Review material must include every tracked file changed by the task contract.
6. It must also include complete proposed additions for new/untracked
   implementation or memory files; ordinary `git diff` does not show them.
7. ChatGPT gives PASS, one narrow repair, or blocker before the user-controlled
   commit/push.
8. After push, ChatGPT reviews the pushed repository before farm validation.

CURRENT must already have a materially correct, push-stable NEXT: the sole
ordinary NEXT cannot be the commit/push itself. State the substantive next
gate conditional on user-controlled commit/push and pushed-state review, so it
remains accurate across that transition. During pushed-state review, explicitly
check source identity, exact gate-relevant changed paths/blobs, and CURRENT/NEXT
continuity. When the pushed state matches the source-reviewed candidate and
memory remains materially accurate, proceed directly to the substantive gate.
No new source-change contract or broad memory reconciliation is required merely
because a push occurred. A separate reconciliation is warranted only for
concrete ambiguity about source identity, accepted evidence/status, active
ownership/blocker, frozen scientific interface, or exact NEXT. Record cosmetic
or historical wording drift and batch it at the next checkpoint or milestone
audit unless it changes active meaning.

Do not require staging merely to review a diff. For each task, enumerate its
actual scoped paths. A safe generic POSIX/Git-Bash pattern is:

```bash
git diff -- <tracked paths> > kaonlt_review.diff
git diff --no-index -- /dev/null <new-file> >> kaonlt_review.diff || true
```

Here `|| true` is permitted only because `git diff --no-index` returns status 1
when it displays a difference; it must never suppress an analysis or test
failure. Repeat the second command for every intended new/untracked file. The
root-level review bundle is temporary and must be removed before commit/push
unless it is explicitly intended to be tracked.

## Contract policy

Before a farm handoff, audit input authority -> producer/materializer/analyzer/
renderer (if any) -> verification/checker (if any) -> collector/packager ->
invocation owner -> expected returned artifact. Every required executable step
must have tracked source, be locally deterministic where possible, independently
source reviewed, pushed, and pushed-state reviewed. If multi-step orchestration
lacks a tracked, reviewed driver, mark the farm gate `BLOCKED` and return to a
source-changing contract. A single reviewed CLI may directly own a complete
single-command operation. A reviewed profile alone does not establish readiness.

Each task contract must state:

- exact starting HEAD;
- objective;
- allowed files;
- frozen files;
- scientific ownership;
- exact before/after behavior;
- preserved source/runtime path;
- positive checks;
- negative checks;
- regression checks;
- forbidden shortcuts or fallbacks;
- local validation;
- farm-validation boundary;
- diff audit;
- acceptance criteria; and
- hard stop.

For a farm-bound task, the contract also names the complete execution chain,
the source owner of each executable step, whether a driver is needed, and the
missing-orchestration blocking rule. Its memory-health section requires byte
counts for CURRENT/MEMORY/CURRENT_HANDOFF, manifest write/check after memory
changes, and the ordinary task-final check:
`<PYTHON> -B tools/check_memory_health.py --root .`.
This command returns nonzero for hard health failures. Require
`--fail-on-warning` only for an explicit milestone/zero-warning audit, a
memory-hardening task that owns warning elimination, or an active warning
already classified as materially blocking. A normal source/science cycle may
complete with all hard checks passing and a recorded nonblocking warning
explicitly scheduled for the next checkpoint or milestone audit.
After substantial implementation, source review reconciliation, closure, or
pushed-state handoff, Codex reports local health and ChatGPT independently
verifies available evidence using the [maintenance format](MAINTENANCE.md).
Hard failures remain blocking. Classify each warning as blocking or nonblocking.
Warnings block only for material ambiguity about current source identity,
accepted evidence/status, frozen scientific interfaces, active scientific
ownership/blocker, or exact NEXT, or a hard integrity violation; record
nonblocking warnings and batch correction at a checkpoint or milestone audit.
Each substantive loop
must deliver a plot, yield table, validated numerical comparison, accepted
runtime artifact, or one directly evidenced runtime/scientific blocker with
one coherent repair. Documentation alone is not scientific completion.

Codex's summary is not proof: actual source and diff must be reviewed. If the
branch and HEAD have not moved, work remains local/proposed rather than
implemented on the branch. The frozen [CODEX contract template](templates/CODEX_CONTRACT.md)
instantiates this policy.

This record does not own farm paths, user preferences, memory maintenance, or
detailed farm packaging commands.
