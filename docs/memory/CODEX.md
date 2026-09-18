# KaonLT Codex workflow

Codex may inspect live source, read targeted repository memory, execute a
tracked phase or fix contract, change only allowlisted files, run deterministic
local checks, trace source/runtime paths, audit its own diff, and update
source-review memory when that contract permits it.

Codex must not commit, push, update remote refs, initiate Jefferson Lab farm
execution, or claim Jefferson Lab farm execution, ROOT/PyROOT integration,
full `main.py` validation, rendered procedure-PDF farm behavior, final
kinematic production validation, or phase closure without accepted farm
evidence where farm evidence is required. ChatGPT independently audits Codex's
actual diff. The user alone commits/pushes accepted changes and runs farm
validation. Workflow: Codex local changes -> ChatGPT audit -> user commit/push
-> user farm run -> ChatGPT evidence review.

## Task classes

1. **Implementation task:** make the approved source change, run local tests,
   audit the diff, and update source-review memory only when allowlisted.
2. **Closure task:** reconcile documentation, evidence, and status only after
   accepted farm evidence. It must not change scientific source.

Neither class may infer runtime acceptance from local checks. A contract must
identify its exact starting source, allowed and frozen files, preserved path,
local validation, farm-validation boundary, and hard-stop condition.
