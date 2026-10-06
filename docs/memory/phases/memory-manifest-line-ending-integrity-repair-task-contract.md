# KaonLT pushed-state memory-manifest line-ending integrity repair task contract

## 1. Authority and task class

This is a **narrow repository-memory integrity repair** discovered during the
post-push synchronization audit after the science-first E.8 re-anchor.

It is not a scientific phase, not a Method-A/Method-B change, not a farm task,
and not permission to reopen any accepted closure.

Repository authority remains:

```text
current source/diff
-> applicable runtime evidence
-> tracked repository memory
-> older chat/history
```

Follow repository-root `AGENTS.md`, then the five-file startup core in its exact
required order before doing any work.

The authoritative repository is:

```text
https://github.com/trottar/lt_analysis/tree/test
```

## 2. Exact starting identity and hard start gate

Expected pushed source:

```text
branch: test
HEAD: 9bf56d2cb7b43be23ac0419df9985a2c9518b369
origin/test: 9bf56d2cb7b43be23ac0419df9985a2c9518b369
```

Before changing anything, run the ordinary branch/HEAD/origin/worktree audit.

The only pre-existing untracked path permitted by this task is this contract
after the user moves it into the repository:

```text
docs/memory/phases/memory-manifest-line-ending-integrity-repair-task-contract.md
```

If branch, HEAD, or `origin/test` differs; if any unrelated tracked or untracked
state exists; or if the contract path is not the only allowed pre-existing
untracked path, **STOP as `BLOCKED`**. Do not reset, stash, clean, overwrite, or
repair unrelated user state.

## 3. Evidence that triggered this repair

The pushed `test` branch is at:

```text
9bf56d2cb7b43be23ac0419df9985a2c9518b369
```

The push from the prior farm source changed only repository-memory paths; no
scientific/runtime source changed.

The live Git blobs and the committed `docs/memory/manifest.json` disagree.

### 3.1 `docs/memory/CURRENT.md`

Manifest entry:

```text
bytes: 8180
sha256: c17155508da55054beac8ad3e3ebb896dbe6bfd3b1716be60852481521fc5e66
```

Live Git blob:

```text
bytes: 8041
sha256: 3d7cefff483ddc87993a98295a02821cf4a8f5737edf599d004a9ec394356ebc
```

The file has 139 LF newlines, and the manifest byte count exceeds the live blob
by exactly 139 bytes. Replacing every LF in the live blob by CRLF produces
exactly the manifest byte count and manifest SHA-256.

### 3.2 `docs/memory/roadmap/STATUS.md`

Manifest entry:

```text
bytes: 23367
sha256: 8e73663b034b6ae2160ebc6f2c660027caea0d54d159e03f7b5bcc5d9fb85b23
```

Live Git blob:

```text
bytes: 23026
sha256: a5f8e1ff6e2596e38f3b50a50cb780153e16bb1e5dfe50ee73c282588d16190f
```

The file has 341 LF newlines, and the manifest byte count exceeds the live blob
by exactly 341 bytes. Replacing every LF in the live blob by CRLF produces
exactly the manifest byte count and manifest SHA-256.

### 3.3 `docs/memory/MEMORY.md`

Manifest entry:

```text
bytes: 20386
sha256: 7b01e4c2831315cebb4fe130aba3691b25f13a1446796f82382bbe68611b59c8
```

Live Git blob:

```text
bytes: 20085
sha256: ac8bac01e02214e9f90438e635e661f76bed874b7681d26f6f6c93fe09e69bc2
```

This file also has a pushed manifest mismatch. The exact pre-push mixed-EOL
layout is not independently reconstructed here.

### 3.4 Evidence classification

```text
SOURCE VERIFIED:
  - live test HEAD above;
  - exact changed-path set of the preceding push;
  - live Git blob identities and sizes above;
  - committed manifest entries above;
  - exact CRLF reconstruction for CURRENT.md and roadmap/STATUS.md;
  - update_memory_manifest.py hashes raw working-tree bytes;
  - check_memory_health.py invokes update_memory_manifest.py --check and treats
    manifest mismatch as a hard health failure.

INFERENCE:
  - the prior local manifest was generated from CRLF or mixed-EOL working-tree
    representations while Git stored normalized LF blobs.

RUNTIME VERIFIED:
  - NONE for this repair.

NOT VERIFIED:
  - the user's present local WSL worktree state until Codex performs the start gate.
```

This is a real pushed-state manifest-integrity failure. It blocks completion of
the pushed-state synchronization gate. It does **not** establish any scientific
defect.

## 4. Objective

Make the memory manifest stable across the user's Windows/WSL Git line-ending
configuration **without changing scientific or active-state content**.

The smallest accepted repair is:

1. add a repository-root `.gitattributes` rule forcing LF working-tree and Git
   representation for all `docs/memory/**` files;
2. normalize the local working-tree EOL representation of every versionable
   `docs/memory/**` text file to LF without changing semantic content;
3. regenerate `docs/memory/manifest.json` only after that normalization;
4. verify that the manifest, ordinary memory health, bootstrap, and diff gates
   pass;
5. prepare one complete review bundle and stop before commit/push.

Do **not** redesign the manifest tool unless this exact LF-discipline repair
cannot satisfy the deterministic checks. If it cannot, stop `BLOCKED` and
report the concrete reason rather than inventing a second design.

## 5. Exact allowed versioned paths

Only these paths may appear in the final Git diff:

```text
.gitattributes
docs/memory/manifest.json
docs/memory/phases/memory-manifest-line-ending-integrity-repair-task-contract.md
```

The contract is already user-supplied task authority and may appear as a new
tracked file.

Local-only EOL normalization of existing `docs/memory/**` files is permitted
only when it produces **zero semantic/content diff** against the index.

If any other path appears in `git diff`, `git diff --cached`, or the intended
review bundle, **STOP as `BLOCKED`**.

## 6. Frozen files and scientific ownership

Everything not listed in Section 5 is frozen as versioned content.

In particular, do not modify:

```text
src/
testing/
farm_env/
background_samples/
run_Prod_Analysis.sh
set_SymLinks.sh
tools/
AGENTS.md
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/USER.md
docs/memory/handoffs/CURRENT_HANDOFF.md
docs/memory/roadmap/STATUS.md
docs/memory/evidence/
docs/memory/decisions/
docs/memory/investigations/
```

No scientific ownership changes.

Preserve exactly:

```text
random subtraction
dummy subtraction
slow-proton subtraction
baseline pion subtraction
no_empirical_residual
SIMC
yield extraction
uncertainties
cuts
templates
priors
binning
efficiencies
acceptance
L/T separation
cross sections
```

Method A remains detached/non-production.

Method B remains diagnostic/cross-check only and numerically excluded.

Canonical-five full runtime remains `BLOCKED`.

Canonical-five provenance/identity repair remains `DEFERRED`.

E.8 remains `ACTIVE`.

Final E.8 remains `BLOCKED`.

F.6.4 remains `BLOCKED`.

The substantive science NEXT remains the E.8.2 baseline missing-mass/stage-yield
audit **after** this pushed-state synchronization defect is repaired and the
pushed repair is independently synchronized.

## 7. Required `.gitattributes` behavior

Create repository-root:

```text
.gitattributes
```

with exactly the narrow memory rule:

```gitattributes
docs/memory/** text eol=lf
```

A final newline is required.

Do not add broad repository-wide EOL rules. Do not change source-file attributes.

## 8. EOL audit before normalization

Before normalizing any bytes, record:

```text
git ls-files --eol -- docs/memory
```

Audit the **index** representation.

Expected safe case:

```text
tracked docs/memory text files are index LF
while one or more working-tree files may be CRLF or mixed
```

Hard stop if any tracked `docs/memory/**` file reports an indexed representation
whose committed content would be changed by applying the LF attribute rule,
including an `i/crlf`, `i/mixed`, or otherwise anomalous indexed text state.

Do not use `git add --renormalize`.

Do not stage files merely to perform this audit.

## 9. Local-only EOL normalization

After adding `.gitattributes`, normalize line endings to LF for all
**versionable text files under `docs/memory/**`**, including this untracked task
contract, before regenerating the manifest.

Requirements:

- preserve decoded text exactly apart from CRLF/CR -> LF conversion;
- preserve the required final newline;
- do not alter whitespace other than line-ending bytes;
- do not change encoding;
- do not touch files outside `docs/memory/**`;
- do not use checkout/reset/stash/clean as a shortcut.

After normalization, existing tracked memory files must have zero semantic Git
diff. If any existing tracked memory file other than `manifest.json` shows a
content diff, stop `BLOCKED` and report it.

Re-run:

```text
git ls-files --eol -- docs/memory
```

All existing tracked versionable memory text files must now have LF working-tree
representation consistent with the new attribute rule.

## 10. Manifest regeneration

Discover `<PYTHON>` using `docs/memory/TOOLS.md`.

Only after Section 9 passes, run:

```text
<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/update_memory_manifest.py --root . --check
```

The regenerated manifest must include the new task contract and must reflect LF
bytes for all memory files.

Do not hand-edit manifest hashes or byte counts.

## 11. Positive checks

All of the following must pass:

1. branch/HEAD/origin start gate;
2. worktree allowlist gate;
3. pre-normalization `git ls-files --eol` index audit;
4. exact `.gitattributes` rule;
5. LF normalization with no semantic changes to existing tracked memory files;
6. post-normalization EOL audit;
7. manifest write;
8. manifest check;
9. ordinary memory health;
10. bootstrap;
11. `git diff --check`;
12. final diff allowlist;
13. complete review bundle.

The expected live core sizes after canonical LF normalization should be based on
the repository's LF blobs, not the stale pre-push manifest entries. At minimum,
the current pushed values established in the audit are:

```text
CURRENT.md: 8041 bytes
MEMORY.md: 20085 bytes
CURRENT_HANDOFF.md: 323 bytes
```

If the task itself legitimately changes none of those three contents, those
values must remain unchanged.

## 12. Negative and regression checks

The final candidate must not:

- change any scientific/runtime source;
- change CURRENT's objective, blockers, status, or NEXT;
- change MEMORY's durable scientific content;
- change roadmap status prose;
- modify accepted evidence;
- reopen canonical-five provenance repair;
- promote Method A;
- numerically include Method B;
- authorize a farm run;
- hide the manifest mismatch by weakening or skipping manifest checks;
- modify `tools/update_memory_manifest.py`;
- use a broad repository-wide EOL policy;
- stage unrelated files;
- rely on a machine-specific Git configuration as the fix.

Search the final diff for unexpected changes in:

```text
src/
testing/
tools/
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/roadmap/STATUS.md
docs/memory/evidence/
```

There must be none.

## 13. Local validation and memory health

Run:

```text
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/check_memory_health.py --root .
<PYTHON> -B tools/memory_bootstrap.py --root . --json
git diff --check
```

Do not use `--fail-on-warning` unless the ordinary health command reports a
warning that is materially blocking under `docs/memory/MAINTENANCE.md`.

Report exactly:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
warning classification: <none or each warning: blocking | nonblocking>
manifest check: PASS | FAIL
bootstrap: PASS | FAIL
git diff --check: PASS | FAIL
```

Also report the pre- and post-normalization `git ls-files --eol` summary,
especially any file that was `w/crlf` or `w/mixed` before normalization.

## 14. Farm-validation boundary

There is **no farm validation** for this task.

Do not run Jefferson Lab analysis.

Do not run ROOT/PyROOT.

Do not run `main.py`.

Do not generate or inspect procedure PDFs as part of this repair.

This task restores pushed repository-memory integrity only.

## 15. Diff audit and review bundle

Before stopping, run:

```text
git status --short --untracked-files=all
git diff --stat
git diff --check
```

The final intended paths must be exactly:

```text
.gitattributes
docs/memory/manifest.json
docs/memory/phases/memory-manifest-line-ending-integrity-repair-task-contract.md
```

Create one complete temporary review bundle in repository root:

```text
kaonlt_review.diff
```

It must contain:

- complete tracked diff for `.gitattributes` if Git treats it as tracked in the
  review method used;
- complete `git diff --no-index /dev/null ...` addition for the new
  `.gitattributes` when still untracked;
- complete tracked diff for `docs/memory/manifest.json`;
- complete `git diff --no-index /dev/null ...` addition for this task contract
  while it remains untracked;
- no unrelated user-owned content.

Do not stage merely for review.

If both `.gitattributes` and the contract are untracked, include both complete
`--no-index` additions.

## 16. Acceptance criteria

PASS only if all are true:

1. exact starting branch/HEAD/origin gate passes;
2. no unrelated local state exists;
3. pre-normalization index EOL audit shows the LF attribute will not alter
   committed semantic content;
4. `.gitattributes` contains only the narrow `docs/memory/** text eol=lf` rule;
5. all versionable memory text working-tree bytes are LF before manifest write;
6. existing tracked memory files other than `manifest.json` have no semantic
   content diff;
7. the manifest is regenerated from the LF working tree;
8. manifest check passes;
9. ordinary memory health passes;
10. bootstrap passes;
11. `git diff --check` passes;
12. only the exact allowed versioned paths remain changed/new;
13. CURRENT/MEMORY/roadmap scientific/status text is unchanged;
14. no farm work occurs;
15. `kaonlt_review.diff` is complete;
16. no commit, push, remote-ref update, staging-for-review, reset, stash, or
   cleanup of unrelated state occurs.

## 17. Hard stop

After implementation, deterministic checks, memory-health report, final
allowlist audit, and `kaonlt_review.diff` creation:

**STOP.**

Do not commit.

Do not push.

Do not run the farm.

Return:

- exact starting and ending branch/HEAD/origin observations;
- exact changed/untracked paths;
- concise root-cause confirmation from the EOL audit;
- pre/post EOL summary;
- exact memory-health report;
- manifest/bootstrap/diff-check results;
- path to `kaonlt_review.diff`.

ChatGPT must review the actual diff before the user performs any commit/push.
