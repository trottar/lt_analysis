# E.8.1.Debug.1.Fix.1 — source-review reconciliation and workflow-memory update

## Starting state

Branch: `test`

Committed HEAD remains:

`4c4a18a098d184fe0b13dd1fc57010ac7f5b30c0`

The worktree intentionally contains the uncommitted E.8.1.Debug.1.Fix.1
implementation and its existing task/history records.

Independent ChatGPT review has now inspected the actual local Fix.1 diff.

The reviewed source repair is a PASS:

- early `LT_ANALYSIS_DEBUG_LEFT_LOW` validation remains;
- `phisetlist = ["Center", "Left", "Right"]` remains active through the
  ordinary Step-2 Diamond stage;
- Diamond artifact registration completes before the debug restriction;
- only afterward does debug mode set `phisetlist = ["Left"]`;
- that restriction occurs before `shift_prep`;
- `src/cuts/diamond.py` was not modified;
- the launcher/canonical-preflight implementation from the pushed Debug.1
  source remains unchanged.

Therefore E.8.1.Debug.1.Fix.1 may now be recorded as `SOURCE REVIEWED`.

No ROOT/PyROOT, full `main.py`, Jefferson Lab farm, or runtime validation has
been performed by ChatGPT.

Codex-reported deterministic tests passed previously, but ChatGPT did not run
them. Preserve that distinction explicitly as `NOT RUN by ChatGPT`.

## Purpose of this pass

This is a memory/history reconciliation only.

Do not modify analysis source or tests in this pass.

Correct the status language for the original flawed Debug.1 implementation,
record Fix.1 as SOURCE REVIEWED, and add the newly requested durable
collaboration workflow for large diff review.

## Required memory correction

### `docs/memory/CURRENT.md`

The pushed original Debug.1 source:

`4c4a18a098d184fe0b13dd1fc57010ac7f5b30c0`

must NOT itself be labeled `SOURCE REVIEWED`.

Independent review of that implementation found a concrete Diamond-cut
blocker: its pre-Step-2 Left restriction bypassed the ordinary
Center-produced common Diamond polygon and allowed the Left-low fallback cut.

Describe that source as reviewed with a blocker / superseded by Fix.1 without
assigning it the `SOURCE REVIEWED` status.

Record:

- E.8.1.Debug.1.Fix.1 — `SOURCE REVIEWED`
- source review was based on ChatGPT inspection of the actual local diff;
- tests reported by Codex were not rerun by ChatGPT;
- no ROOT/PyROOT/farm/runtime validation exists;
- user commit/push remains pending;
- after push, the next narrow operational step is the targeted `-d`
  Jefferson Lab debug run;
- scientific E.8.1 remains `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`;
- after the debug prerequisite is resolved, the E.8.1 canonical-five-setting
  gate remains unchanged, with `Q4p4W2p74 / Left / highe` inspected first.

Do not distort closed or unrelated phase status.

### Existing Fix.1 phase record

Update:

`docs/memory/phases/e8-1-debug-1-fix1-preserve-diamond-cut.md`

to record the independent actual-diff review outcome.

It must state:

- starting committed HEAD was `4c4a18a098d184fe0b13dd1fc57010ac7f5b30c0`;
- original Debug.1 review found the Diamond-cut ownership regression;
- Fix.1 moves only the downstream debug restriction;
- ordinary Center/Left/Right Diamond preparation and artifact registration
  precede the Left-only restriction;
- downstream `shift_prep` and later processing are Left-only in debug mode;
- analysis source outside `src/main.py` is unchanged by Fix.1;
- Fix.1 status is `SOURCE REVIEWED`;
- Codex-reported deterministic checks passed;
- `NOT RUN by ChatGPT`;
- no farm/runtime claim;
- user commit/push and targeted farm debug run remain pending.

Preserve historical chronology. Do not rewrite history to imply the first
Debug.1 implementation passed review.

## Durable workflow-memory update

The user explicitly requested the following workflow to become durable
repository memory.

### `docs/memory/USER.md`

Add the collaboration preference:

After Codex completes a local source-changing task, ChatGPT must inspect the
actual diff before the user commits/pushes.

If the diff is too large to paste comfortably, do not require manual terminal
copy/paste. Give an exact command that writes the review material to a diff
file and have the user upload that file for review.

Prefer creating transient review artifacts outside the repository worktree,
for example under `/tmp`, so they do not pollute `git status` or risk being
committed.

### `docs/memory/CODEX.md`

Add the corresponding source-change workflow rule:

1. Codex implements locally and stops before commit/push.
2. ChatGPT reviews the actual diff, not the Codex summary.
3. For a small diff, terminal output may be pasted directly.
4. For a large diff, create a review bundle/diff file outside the worktree and
   upload it to ChatGPT.
5. The review material must contain all changed tracked files relevant to the
   contract.
6. New/untracked implementation or memory files must also be included in the
   review material; ordinary `git diff` alone does not show their contents.
7. ChatGPT gives PASS / narrow repair / blocker before user-controlled
   commit/push.
8. After push, ChatGPT reviews the pushed repository before farm validation.

Document a safe example pattern equivalent to:

`git diff -- <tracked paths> > /tmp/kaonlt_review.diff`

and, for new untracked files, append their complete proposed additions using
an appropriate non-destructive diff representation such as:

`git diff --no-index -- /dev/null <new-file> >> /tmp/kaonlt_review.diff || true`

The `|| true` here is allowed only because `git diff --no-index` returns 1
when differences exist; it must not be used to suppress analysis/test
failures.

The command example must remain generic enough to list the actual scoped
paths for each task.

Do not require users to stage files merely to make them reviewable.

## Files allowed in this reconciliation pass

Memory changes only:

- `docs/memory/CURRENT.md`
- `docs/memory/USER.md`
- `docs/memory/CODEX.md`
- `docs/memory/phases/e8-1-debug-1-fix1-preserve-diamond-cut.md`
- this task contract
- `docs/memory/manifest.json`

The already-local source/test changes:

- `src/main.py`
- `testing/test_run_prod_analysis_debug_left_low.py`

must remain byte-for-byte unchanged during this reconciliation pass.

The original Fix.1 implementation task contract must remain unchanged.

Do not modify:

- `run_Prod_Analysis.sh`
- `src/cuts/diamond.py`
- canonical-binning source
- subtraction/background/SIMC/yield source
- E.8.1 collectors/profiles/wrappers
- unrelated repository memory.

## Required checks

Before finishing:

1. verify `src/main.py` and the focused test have no new diff relative to the
   worktree state at the beginning of this reconciliation pass;
2. regenerate `docs/memory/manifest.json`;
3. run applicable manifest/integrity checks;
4. run `tools/check_memory_health.py --root .`;
5. run `tools/memory_bootstrap.py --root . --json`;
6. run applicable memory-health tests;
7. run `git diff --check`;
8. inspect the final changed-path list.

No farm command.
No ROOT/PyROOT execution.
No commit.
No push.

## Final review bundle

Because this worktree contains both tracked modifications and new untracked
memory files, provide exact commands for the user to produce one complete
review file under `/tmp`.

The review file must include:

- tracked diffs;
- full diff representations of all new/untracked task/phase records that are
  intended for the eventual commit.

Do not place the review artifact in the repository.

## Acceptance

This reconciliation passes locally only if:

- original `4c4a18a...` Debug.1 is not mislabeled SOURCE REVIEWED;
- Fix.1 is recorded as SOURCE REVIEWED;
- local-vs-farm evidence boundaries remain explicit;
- actual-diff review-before-push is durable workflow memory;
- large diffs use uploaded review files rather than manual pasting;
- untracked files are explicitly included in future review bundles;
- analysis source/test files are untouched by this memory-only pass;
- manifest and memory checks pass.

## Hard stop

After memory updates and deterministic checks:

STOP.

Do not commit.
Do not push.
Do not run the farm.

Return:

- exact files modified;
- checks and results;
- confirmation that source/test bytes were untouched in this pass;
- exact `/tmp` review-bundle commands;
- no Git commit/push handoff yet, because ChatGPT will inspect the final
  memory reconciliation diff first.
