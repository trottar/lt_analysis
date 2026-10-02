# KaonLT workflow-continuity/chat-health hardening closure — task contract

## Status

`ACTIVE`

This is the final narrow closure task for the workflow-continuity and in-chat
health hardening completed on 2026-10-02.

Scientific analysis remains paused by user instruction until this closure is
finished. This task owns no scientific/runtime behavior.

## Required starting identity

Required branch:

`test`

Required starting HEAD and locally recorded `origin/test`:

`8a6cad9c4c85bbe8fe9b32c15cc31402d4758af9`

Expected commit subject:

`Track KaonLT root workflow authority`

Codex must establish and report, before editing:

```text
git status --short --branch --untracked-files=all
git rev-parse HEAD
git branch --show-current
git rev-parse origin/test
git log --oneline -1
```

Required:

- branch is exactly `test`;
- local HEAD is exactly the required starting HEAD;
- `origin/test` is exactly the required starting HEAD;
- there are no tracked modifications.

Do not pull, reset, clean, stash, commit, push, move, rename, delete, or discard
pre-existing local state.

## Expected pre-existing untracked state

The latest user-supplied worktree status contained exactly these three
pre-existing untracked files:

```text
docs/memory/phases/workflow-continuity-hardening-task-contract.md
kaonlt_hardening_review.diff
kaonlt_review.diff
```

At task start they are known, intentional local artifacts.

The two root-level review bundles remain preserved throughout implementation:

```text
kaonlt_hardening_review.diff
kaonlt_review.diff
```

Do not edit, move, delete, stage, include them in the closure review bundle, or
include them in the user Git handoff before pushed-state synchronization.

The superseded untracked memory contract

```text
docs/memory/phases/workflow-continuity-hardening-task-contract.md
```

has a narrower rule because the manifest intentionally inventories every
tracked plus nonignored versionable file under `docs/memory`. Keeping that
superseded contract in place makes an actual-worktree manifest check
deterministically fail.

Therefore:

- preserve it through implementation and the first ChatGPT actual-diff review;
- after that actual-diff review passes, the **user may delete only that exact
  superseded untracked contract** before the final actual-worktree manifest and
  memory-health validation;
- Codex must not delete it;
- its deletion is local cleanup of an untracked superseded artifact, not a
  repository source change and not a scientific/runtime change;
- do not modify the manifest writer/checker to ignore it;
- do not add it to the manifest;
- do not stage or commit it.

After the new closure contract is placed, that contract is also an intended
untracked file for this task:

```text
docs/memory/phases/workflow-continuity-and-chat-health-hardening-closure-task-contract.md
```

Any other tracked modification or unrelated untracked file at preflight is a
hard stop unless it is proven to be an intentional artifact of this task.

## Mandatory startup sequence

Read in full, in this exact order:

1. repository-root `AGENTS.md`;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read only the task-relevant records:

- `docs/memory/MAINTENANCE.md`;
- `docs/memory/phases/workflow-continuity-and-chat-health-hardening.md`;
- this contract.

Use current source/diff over stale chat summaries. Do not reopen scientific
architecture or historical closures.

## Established closure facts

The following facts are established before this task:

- repository-side workflow/chat-health hardening received independent
  ChatGPT actual-diff review;
- the reviewed hardening was published in two user-controlled commits;
- pushed-state synchronization subsequently passed at
  `8a6cad9c4c85bbe8fe9b32c15cc31402d4758af9`;
- root `AGENTS.md` is tracked publicly and is the universal startup authority;
- the matching ChatGPT Project instructions were updated outside the
  repository;
- the Project file `KaonLT_Environment.md` was also updated outside the
  repository;
- those Project-side changes contain stable behavior/environment configuration,
  not mutable scientific HEAD/phase/NEXT state;
- no farm/runtime/scientific acceptance was created by the hardening;
- scientific analysis remains paused until the user explicitly resumes it.

This closure task does not re-review or redesign those already accepted
hardening behaviors. It fixes the final omissions identified by the post-
hardening audit and makes repository memory push-stable.

## Objective

Perform exactly four closure actions:

1. restore three explicitly discussed synchronization-receipt fields to the
   full in-chat health-check template;
2. deterministically enforce those fields in repository memory-health checks;
3. update CURRENT and the hardening phase record so they accurately state that
   repository pushed-state synchronization and external Project configuration
   are complete while preserving the scientific state and substantive NEXT;
4. regenerate integrity metadata and provide one complete actual-diff review
   bundle.

No other change is authorized.

## Allowed paths

Only these repository paths may change:

- `docs/memory/MAINTENANCE.md`
- `docs/memory/CURRENT.md`
- `docs/memory/phases/workflow-continuity-and-chat-health-hardening.md`
- `docs/memory/phases/workflow-continuity-and-chat-health-hardening-closure-task-contract.md`
- `docs/memory/manifest.json`
- `tools/check_memory_health.py`
- `testing/test_memory_health.py`

Everything else is frozen.

The three pre-existing untracked files named above are frozen-in-place and are
not task outputs.

## Frozen scientific/runtime scope

Do not modify:

- `src/**`;
- any scientific/runtime test other than the explicitly allowlisted
  `testing/test_memory_health.py`;
- any tool other than the explicitly allowlisted
  `tools/check_memory_health.py`;
- farm profiles, owners, wrappers, collectors or scientific/runtime checkers;
- `docs/memory/MEMORY.md`;
- `docs/memory/USER.md`;
- `docs/memory/CODEX.md`;
- `docs/memory/COMMUNICATION.md`;
- `docs/memory/TOOLS.md`;
- `docs/memory/LEARNINGS.md`;
- `docs/memory/README.md`;
- `docs/memory/AGENTS.md`;
- root `AGENTS.md`;
- `docs/memory/handoffs/**`;
- `docs/memory/roadmap/**`;
- `docs/memory/evidence/**`;
- `docs/memory/decisions/**`;
- `docs/memory/investigations/**`;
- any existing scientific phase/fix record;
- any ChatGPT Project instruction or Project file.

Do not modify `.gitignore` or local/global Git ignore configuration.

## Scientific ownership

This task owns no physics quantity.

Preserve exactly:

- all accepted historical closures;
- Fix.5.4/Fix.5.5/Fix.5.6 scientific/runtime implementation;
- current background profile and subtraction ordering;
- SIMC;
- yield extraction;
- cuts/templates/priors/binning/efficiencies/acceptance;
- Method A detached/non-production role;
- Method B diagnostic-only/numerically excluded role;
- F.6.3/E.8.4 accepted narrow runtime closure already recorded;
- canonical-five/F.6.4 blockers;
- current failed Fix.5.6 farm-gate interpretation.

No farm command, rerun, artifact packaging, scientific interpretation, Method-A
promotion, or production change is authorized.

## 1. Complete the full health-check receipt

In `docs/memory/MAINTENANCE.md`, keep the existing full health-check structure
and add the three fields explicitly discussed during hardening:

```text
startup core: READ
CURRENT/source consistency: PASS | DRIFT
task class: <source review | contract | pushed-state review | farm evidence | other>
```

They must appear in the full `KaonLT health check` receipt.

Requirements:

- `startup core: READ` asserts that the universal five-file core was actually
  read for the substantial task; do not permit a fabricated `READ` value before
  the startup sequence is completed.
- `CURRENT/source consistency: PASS | DRIFT` explicitly exposes material
  alignment between live source identity and CURRENT's active/source claims.
- `task class` makes the workflow state visible and prevents a chat from
  silently switching from one gate class to another.
- retain the existing `current gate` and `gate status` fields;
- retain the existing evidence labels and ownership section;
- retain the existing periodic health pulse;
- do not create a second health-check format.

The final full receipt must visibly expose both workflow identity and evidence
identity.

## 2. Deterministically enforce the receipt markers

Update only the existing memory-health checker and its dedicated deterministic
test:

```text
tools/check_memory_health.py
testing/test_memory_health.py
```

Add the smallest deterministic invariant sufficient to fail if
`docs/memory/MAINTENANCE.md` loses any of these required full-health receipt
markers:

```text
startup core: READ
CURRENT/source consistency: PASS | DRIFT
task class:
current gate:
gate status:
SOURCE VERIFIED:
RUNTIME VERIFIED:
MEMORY ONLY:
INFERENCE:
NOT VERIFIED:
```

Requirements:

- validation must inspect the canonical MAINTENANCE document;
- it must not parse or validate chat transcripts;
- it must not create mutable active state;
- it must preserve all existing schema-3/startup/manifest/authority checks;
- add focused regression tests for at least the three newly restored fields;
- do not broaden into unrelated linting or formatting policy.

## 3. Reconcile `CURRENT.md`

Make only the minimum push-stable edits necessary.

The Current Work Item must state that:

- workflow/chat-health hardening is `SOURCE REVIEWED`;
- repository-side hardening passed pushed-state synchronization through
  `8a6cad9c4c85bbe8fe9b32c15cc31402d4758af9`;
- the matching ChatGPT Project instructions and `KaonLT_Environment.md` are now
  updated externally;
- scientific work remains paused until explicit user resumption;
- this closure creates no farm/runtime acceptance.

Preserve Fix.5.6 as:

`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`

Preserve the failed owner-gate facts:

- analysis child completed;
- failure at `verify_artifacts/page_manifest_setting_invalid`;
- collection/ZIP did not begin;
- generated PDF/manifest are failed-gate diagnostic artifacts only.

Preserve all existing scientific blockers and accepted closures.

### Authoritative NEXT

Remove the now-consumed condition "after workflow hardening is complete" from
the sole active NEXT, because the hardening is now complete.

The substantive scientific NEXT must remain:

when the user explicitly resumes scientific work, diagnose the exact
`page_manifest_setting_invalid` failure and corresponding page/payload
provenance before any rerun, packaging, or scientific interpretation, while
preserving the existing SIMC absolute-unit blocker.

Do not authorize a farm command.

Keep CURRENT compact. Avoid introducing a soft-size warning if the same facts
can be expressed by replacing stale hardening wording rather than appending new
history.

## 4. Reconcile the hardening phase record

Update:

`docs/memory/phases/workflow-continuity-and-chat-health-hardening.md`

Keep status:

`SOURCE REVIEWED`

Replace consumed pre-push wording so the record now states:

- independent actual-diff review passed;
- user-controlled repository publication occurred;
- pushed-state synchronization passed through
  `8a6cad9c4c85bbe8fe9b32c15cc31402d4758af9`;
- root `AGENTS.md` is publicly tracked;
- matching ChatGPT Project instructions and `KaonLT_Environment.md` are
  configured externally;
- the final closure restores and enforces the three omitted synchronization
  receipt fields;
- no farm/runtime/scientific validation follows;
- science remains paused until explicit user resumption.

Do not imply that pushed-state synchronization is still pending.

Do not claim the external Project configuration is repository evidence or
scientific evidence; describe it only as external workflow/environment
configuration confirmed for continuity.

Do not rewrite historical scientific phase records.

## 5. Manifest and deterministic validation

Discover a working Python interpreter according to repository policy.

After edits:

```text
<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/check_memory_health.py --root .
<PYTHON> -B tools/memory_bootstrap.py --self-test
<PYTHON> -B testing/test_memory_health.py
git diff --check
```

Run targeted textual checks confirming:

- the full MAINTENANCE health receipt contains all three newly restored fields;
- the checker enforces them;
- dedicated tests fail when each newly restored field is removed or changed;
- CURRENT contains exactly one active `NEXT —`;
- CURRENT no longer says workflow hardening still needs to be completed;
- the hardening phase no longer says pushed-state synchronization is unclaimed
  or pending;
- no scientific/runtime source changed;
- changed paths are exactly within this task's allowlist.

Report:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
warning classification: <none | blocking | nonblocking>
manifest check: PASS | FAIL
bootstrap self-test: PASS | FAIL
memory-health tests: <count passed> | FAIL
```

## 6. Diff audit and review bundle

Report:

```text
git status --short --branch --untracked-files=all
git diff --stat
git diff --name-only
git diff --check
```

Do not stage files merely for review.

Create exactly one new temporary closure review bundle in the repository root:

`kaonlt_hardening_closure_review.diff`

It must contain:

1. the complete tracked diff for every modified allowlisted tracked file;
2. the complete `git diff --no-index /dev/null ...` addition for the new closure
   contract while it is untracked.

It must exclude:

- `docs/memory/phases/workflow-continuity-hardening-task-contract.md`;
- `kaonlt_hardening_review.diff`;
- `kaonlt_review.diff`.

Do not overwrite any pre-existing review bundle.

Use `|| true` only for `git diff --no-index`, where exit status 1 means a
difference was emitted.

## 7. User-controlled Git handoff

After all deterministic checks pass, provide exact non-executed commands for the
user to stage only:

```text
docs/memory/MAINTENANCE.md
docs/memory/CURRENT.md
docs/memory/phases/workflow-continuity-and-chat-health-hardening.md
docs/memory/phases/workflow-continuity-and-chat-health-hardening-closure-task-contract.md
docs/memory/manifest.json
tools/check_memory_health.py
testing/test_memory_health.py
```

Use commit message:

```text
Close KaonLT workflow hardening
```

Do not stage or commit any review bundle or the superseded untracked contract.

Codex must not commit or push.

## 8. Local cleanup sequencing

Codex must not delete known temporary untracked files.

### Pre-validation cleanup after first actual-diff PASS

The first ChatGPT actual-diff review has passed for the closure candidate.

Because the repository manifest intentionally inventories nonignored
versionable files under `docs/memory`, the user may now delete only:

```text
docs/memory/phases/workflow-continuity-hardening-task-contract.md
```

before the final actual-worktree manifest and memory-health validation.

This is required to remove the already-diagnosed manifest-inventory conflict.
Do not replace this cleanup with a manifest exception, ignore rule, or tracked
entry for the superseded contract.

After that single deletion, Codex may regenerate the allowlisted manifest and
rerun the full deterministic validation. Because this contract itself is
revised, the closure review bundle must then be regenerated and returned for a
final ChatGPT actual-diff review.

### Final cleanup after pushed-state synchronization

Only after:

1. final ChatGPT actual-diff review PASS;
2. user commit/push;
3. final ChatGPT pushed-state synchronization PASS;

the user may remove the remaining temporary local review bundles:

```text
kaonlt_review.diff
kaonlt_hardening_review.diff
kaonlt_hardening_closure_review.diff
```

Then verify:

```text
git status --short --branch --untracked-files=all
```

The intended final local baseline is a clean tracked/untracked worktree on
`test`.

## Before/after behavior

### Before

The repository hardening is substantively complete, but:

- the full health receipt does not explicitly expose `startup core: READ`;
- it does not explicitly expose `CURRENT/source consistency: PASS | DRIFT`;
- it does not expose `task class`;
- CURRENT still phrases hardening/Project configuration as unfinished;
- the hardening phase record still says pushed-state synchronization is not
  claimed;
- the superseded untracked memory contract must be removed after the first
  actual-diff PASS so the actual-worktree manifest can validate;
- root-level temporary review bundles remain preserved until pushed-state PASS.

### After

- the full health receipt exposes startup-read state, source/CURRENT consistency,
  task class, current gate, gate status, evidence labels and ownership;
- deterministic memory-health checks prevent regression of the restored fields;
- CURRENT accurately records hardening and external Project configuration as
  complete while keeping science paused and preserving the substantive
  scientific NEXT;
- the hardening phase is push-stable and no longer contains consumed
  pre-synchronization wording;
- scientific/runtime source is unchanged;
- after the first actual-diff PASS, the superseded untracked memory contract
  is removed so actual-worktree manifest validation can pass;
- after pushed-state PASS, the remaining root-level review bundles may be
  removed to restore a clean local baseline.

## Acceptance criteria

PASS only if:

- branch/HEAD/origin preflight matches the required starting identity;
- only allowlisted files change;
- all three omitted full-health receipt fields are restored;
- deterministic checker/tests enforce the restored fields;
- existing health/pulse/evidence-label/re-anchor semantics remain intact;
- CURRENT accurately records completed hardening and external Project setup;
- CURRENT retains exactly one substantive scientific NEXT;
- Fix.5.6 remains `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`;
- failed-gate artifact admissibility remains unchanged;
- hardening phase remains `SOURCE REVIEWED`;
- no farm/runtime/scientific source or scientific test changes;
- manifest check passes;
- memory health has no hard failure;
- bootstrap and dedicated memory tests pass;
- diff audit is allowlist-clean;
- one complete `kaonlt_hardening_closure_review.diff` exists.

## Hard stop

Stop and report `BLOCKED` without broadening scope if:

- branch is not `test`;
- local HEAD or `origin/test` is not
  `8a6cad9c4c85bbe8fe9b32c15cc31402d4758af9`;
- a tracked modification exists before this task;
- pre-existing untracked state differs materially from the explicitly permitted
  files plus the newly placed closure contract, except that the superseded
  untracked memory contract may be absent after the explicitly authorized
  post-review user cleanup;
- a required change would touch a frozen file;
- scientific/runtime source changes appear;
- an existing scientific status/NEXT would need reinterpretation;
- a hard memory-health failure cannot be repaired within the allowlist.

Do not reset, stash, clean, delete, commit, push, run the farm, or invent a
workaround.
