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
