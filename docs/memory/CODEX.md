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
Workflow: Codex local changes -> ChatGPT audit -> user commit/push -> user farm run -> ChatGPT evidence review.

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

Codex's summary is not proof: actual source and diff must be reviewed. If the
branch and HEAD have not moved, work remains local/proposed rather than
implemented on the branch. The frozen [CODEX contract template](templates/CODEX_CONTRACT.md)
instantiates this policy.

This record does not own farm paths, user preferences, memory maintenance, or
detailed farm packaging commands.
