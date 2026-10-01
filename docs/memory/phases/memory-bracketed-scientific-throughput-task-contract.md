# KaonLT workflow reset — memory-bracketed scientific throughput

## Objective

Perform one **memory-only workflow checkpoint** before the next substantive KaonLT scientific/runtime gate.

The repository-memory system remains mandatory. The correction is that memory maintenance must **bracket** substantive work rather than recursively become the work.

The normal loop is:

```text
memory checkpoint
  -> substantive science / implementation
  -> deterministic local test
  -> ChatGPT actual-diff review
  -> user commit/push
  -> ChatGPT pushed-state synchronization review
  -> narrow farm/runtime test
  -> evidence review
  -> memory checkpoint
  -> repeat
```

The commit/push and pushed-state review in the middle are **unavoidable synchronization stages**. They exist so ChatGPT is reviewing the same pushed source that Codex implemented and the farm will run. They are not separate scientific phases and must not trigger recursive memory/reconciliation work unless they expose a concrete source-identity or memory-authority blocker.

This task changes workflow policy and active-state wording only. It must not change scientific source, production behavior, accepted authority, Method-A numerical content, Method-B behavior, plots, yields, cuts, normalizations, or runtime logic.

---

## Exact starting state

Required branch:

```text
test
```

Required pushed HEAD:

```text
712ba32b062772d87fe44efa6865346e5c438827
F4 Refresh.2 Validation.2: add tracked execution owner
```

Before editing, Codex must verify:

```bash
git branch --show-current
git rev-parse HEAD
git status --short --untracked-files=all
git log -1 --oneline
```

Required committed HEAD is exactly the SHA above.

The worktree may contain only this newly placed task contract as an intended pre-existing versionable change. If unrelated versionable changes exist, STOP and report them.

Independent ChatGPT pushed-state review already established that remote `test` points to the required HEAD and that the pushed F.4.Refresh.2 execution owner/profile match the previously source-reviewed candidate. Do not reopen that implementation here.

---

## Mandatory startup

Read, in order:

1. root `AGENTS.md` if present locally;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read only task-relevant records:

```text
docs/memory/CODEX.md
docs/memory/MAINTENANCE.md
docs/memory/LEARNINGS.md
docs/memory/roadmap/STATUS.md
docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner.md
docs/memory/phases/f4-refresh2-validation2-fix1-post-hardening-source-allowlist.md
docs/memory/decisions/farm-validation-bundle-procedure.md
```

Do not reopen scientific architecture or closed phases.

---

## Allowed files

Modify only:

```text
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/USER.md
docs/memory/CODEX.md
docs/memory/MAINTENANCE.md
docs/memory/LEARNINGS.md
docs/memory/roadmap/STATUS.md
docs/memory/manifest.json
```

Create exactly:

```text
docs/memory/decisions/memory-bracketed-scientific-throughput.md
docs/memory/phases/memory-bracketed-scientific-throughput-task-contract.md
```

No other versionable file may change.

---

## Frozen files and scientific boundary

Freeze all analysis/runtime/testing/operational source:

```text
src/
testing/
tools/
run_Prod_Analysis.sh
main.py
```

Also freeze all accepted evidence and authority artifacts.

Do not change:

- F.1-F.6.2 accepted science;
- current F.4.Refresh.2 owner/profile implementation;
- E.8 scientific implementation;
- Method-A values, response construction, application populations, or parent-t normalization;
- Method-B diagnostic-only status;
- random, proton, or pion subtraction;
- cuts, templates, priors, normalizations, binning, efficiencies, acceptance, yield formulas, or cross-section logic.

This is a workflow-policy checkpoint only.

---

## Required policy changes

### 1. Define the normal scientific loop

Repository memory must explicitly state that the default KaonLT development loop is:

```text
checkpoint memory
-> implement one meaningful scientific milestone
-> run deterministic local checks
-> ChatGPT actual-diff review
-> user commit/push
-> ChatGPT pushed-state synchronization review
-> run one narrow farm/runtime gate when required
-> inspect fresh evidence
-> update memory once
-> continue
```

The push and pushed-state review are required synchronization stages, not optional shortcuts and not separate scientific milestones.

### 2. Prevent recursive memory work

Add a durable rule:

> A memory inconsistency blocks scientific progress only when it creates concrete ambiguity about current source identity, accepted evidence, frozen scientific interfaces, active scientific ownership, or the exact next operation.

Otherwise:

- record the issue;
- continue the active scientific milestone;
- repair it at the next memory checkpoint or milestone audit.

Do not start a standalone Codex repair cycle for cosmetic, historical, or wording-only memory drift unless it makes the active state materially ambiguous.

### 3. Batch nonblocking memory corrections

Add a durable rule that accumulated nonblocking memory defects are corrected in one batch at meaningful milestones.

Suitable milestones include:

- accepted farm/runtime evidence;
- a completed scientific implementation;
- a real blocker changing NEXT;
- production-promotion decisions;
- completion of a major E.8/F.6 gate.

Memory maintenance must not become an indefinitely recursive sequence between scientific gates.

### 4. Keep push synchronization lightweight

Update `CODEX.md` / `USER.md` / `MAINTENANCE.md` so the required middle synchronization remains:

```text
Codex local implementation
-> ChatGPT actual-diff review
-> user commit/push
-> ChatGPT pushed-state review
```

But state explicitly:

- no new source-change contract is required merely because a push occurred;
- pushed-state review should verify source identity, exact changed paths/blobs relevant to the gate, and CURRENT/NEXT continuity;
- if the pushed state matches the source-reviewed candidate and memory remains materially accurate, proceed directly to the farm/scientific gate;
- do not perform another broad memory reconciliation merely to restate that the push happened.

### 5. Require a scientific deliverable from each loop

Add a durable throughput rule:

Every substantive scientific loop must end with at least one of:

- a new plot;
- a new yield table;
- a new validated numerical comparison;
- a new accepted runtime artifact;
- or one concrete runtime/scientific blocker with direct evidence and one coherent repair.

Documentation-only output does not count as completion of a scientific loop.

### 6. Current immediate objective

Update CURRENT so that after this workflow checkpoint and its normal source-review/push synchronization, the next substantive gate is:

```text
Q4p4W2p74 F.4.Refresh.2 farm materialize -> verify -> package
```

If that passes, proceed directly to:

```text
Q4p4W2p74 / Left / lowe Method-A reweighting runtime demonstration
```

whose concrete deliverables are:

```text
baseline vs reweighted missing-mass spectrum
per-t baseline vs reweighted comparison
per-(t,phi) baseline yield
per-(t,phi) reweighted yield
absolute yield change
fractional yield change
parent-t preservation check
procedure-PDF presentation of the effect
```

Do not expand to canonical-five Method-A yield presentation before Left/lowe has produced actual before/after evidence.

Do not begin unrelated memory hardening or presentation cleanup while this chain is active unless a concrete blocker requires it.

---

## CURRENT requirements

`docs/memory/CURRENT.md` must remain concise and contain exactly one ordinary NEXT.

It must:

- retain E.8 as `ACTIVE`;
- retain all closed/runtime-validated statuses already established;
- retain F.4.Refresh.2.Validation.2/Fix.1 as `SOURCE REVIEWED`;
- record pushed source identity `712ba32b062772d87fe44efa6865346e5c438827`;
- remove stale wording that says the execution owner still awaits commit/push;
- state that pushed-state synchronization has passed if and only if the repository source directly supports that conclusion from the supplied ChatGPT review;
- identify the next substantive gate as F.4.Refresh.2 farm materialization;
- preserve the downstream Left/lowe Method-A reweighting/yield goal;
- contain no competing NEXT.

The NEXT must be push-stable. It must not make this memory checkpoint's own commit/push the scientific next action.

---

## Decision record

Create:

```text
docs/memory/decisions/memory-bracketed-scientific-throughput.md
```

It must durably record:

1. why the workflow was changed;
2. the normal memory -> science -> validation -> memory loop;
3. the unavoidable diff-review / commit-push / pushed-state synchronization stages;
4. the rule that those synchronization stages do not become independent phases unless they expose a blocker;
5. blocking vs nonblocking memory-error criteria;
6. milestone batching of nonblocking maintenance;
7. the requirement for visible scientific output or a concrete evidence-backed blocker from each substantive loop;
8. preservation of all scientific ownership and farm-validation boundaries.

Keep it concise enough to serve as an operational decision, not a retrospective essay.

---

## CODEX workflow update

Update `docs/memory/CODEX.md` narrowly.

Preserve actual-diff review and user-controlled push.

Clarify that the default source-changing cycle is:

```text
contract
-> Codex implementation + local tests
-> ChatGPT actual-diff review
-> user commit/push
-> ChatGPT pushed-state synchronization
-> farm gate if required
```

Do not require a separate final-pre-push reconciliation task when CURRENT is already materially correct and push-stable.

A separate reconciliation is warranted only if the candidate's memory is materially wrong about source identity, accepted status, blocker, frozen interface, or exact NEXT.

Cosmetic/historical wording drift is deferred to the next checkpoint unless it changes active meaning.

---

## MAINTENANCE update

Update `docs/memory/MAINTENANCE.md` so maintenance occurs primarily:

- at the opening checkpoint for a new substantive milestone;
- after accepted substantive implementation/runtime evidence;
- when a concrete blocker changes active state;
- at major milestone audits.

Retain memory-health tooling and manifest integrity.

Modify the existing rule that every warning automatically blocks every next gate. A warning may block only when it represents a material active-state/provenance ambiguity or violates a hard integrity requirement. Nonblocking warnings must be recorded and batched for the next checkpoint/milestone audit.

Hard failures remain blocking.

---

## USER collaboration update

Update `docs/memory/USER.md` to preserve the user's existing conservative source/farm workflow while adding the stable preference that:

- repository memory exists to improve scientific throughput and prevent regression;
- memory maintenance must not displace substantive physics work;
- push/pushed-state synchronization remains required so ChatGPT and Codex operate on the same source;
- once synchronization passes, proceed directly to the substantive gate;
- favor early narrow runtime tests over prolonged speculative hardening;
- batch nonblocking memory cleanup at milestones.

---

## LEARNINGS update

Add concise generalized lessons:

- test earlier; process cannot substitute for runtime evidence;
- memory is a checkpoint mechanism, not a parallel deliverable stream;
- push synchronization is necessary but should be lightweight;
- nonblocking memory drift should be batched;
- every scientific loop should produce visible evidence or a concrete blocker.

Do not duplicate phase chronology.

---

## Roadmap update

Update only enough to reflect the workflow decision and current substantive next gate.

Do not redesign phase dependencies.

Do not create a second active NEXT in the roadmap.

---

## Local validation

Because this is memory-only, no scientific/unit regression suite is required.

Run:

```bash
python -B tools/update_memory_manifest.py --root . --write
python -B tools/update_memory_manifest.py --root . --check
python -B tools/memory_bootstrap.py --root .
python -m unittest testing.test_memory_health
python -B tools/check_memory_health.py --root . --fail-on-warning
git -c core.safecrlf=false diff --check
```

If the repository uses a different discovered Python interpreter, use the repository-authoritative interpreter from `TOOLS.md`.

Also run targeted text checks proving:

- exactly one ordinary NEXT in CURRENT;
- CURRENT contains the pushed HEAD `712ba32b062772d87fe44efa6865346e5c438827`;
- no current wording says the F.4.Refresh.2 owner still awaits user commit/push;
- the memory -> science -> validation -> memory loop is recorded;
- push/pushed-state synchronization is retained;
- milestone batching of nonblocking memory repair is recorded;
- no scientific/runtime source changed.

---

## Diff audit

Before stopping, report:

```bash
git status --short --untracked-files=all
git diff --stat
git diff --check
git diff -- docs/memory/CURRENT.md docs/memory/MEMORY.md docs/memory/USER.md \
  docs/memory/CODEX.md docs/memory/MAINTENANCE.md docs/memory/LEARNINGS.md \
  docs/memory/roadmap/STATUS.md docs/memory/manifest.json
git diff --no-index -- /dev/null docs/memory/decisions/memory-bracketed-scientific-throughput.md || true
```

The task-contract file itself is already supplied and must remain byte-identical.

Produce one temporary repository-root cumulative review bundle containing the complete tracked diff and complete no-index diff for the new decision record.

Do not stage merely for review.

---

## Acceptance criteria

PASS only if all are true:

1. committed base remains exactly `712ba32b062772d87fe44efa6865346e5c438827`;
2. no scientific/runtime/testing/operational source changed;
3. repository memory clearly defines the memory -> science -> validation -> memory loop;
4. required ChatGPT review -> user push -> pushed-state synchronization remains intact;
5. synchronization stages are explicitly not independent scientific phases;
6. blocking vs nonblocking memory defects are distinguished;
7. nonblocking memory cleanup is batched at checkpoints/milestones;
8. every substantive loop requires visible scientific evidence or a concrete blocker;
9. CURRENT has exactly one push-stable NEXT;
10. CURRENT points to the F.4.Refresh.2 farm materialization gate, followed by Left/lowe Method-A reweighting/yield evidence if that gate passes;
11. all prior scientific statuses and frozen interfaces remain intact;
12. manifest and strict memory-health checks pass;
13. diff check passes;
14. one fresh review bundle is returned;
15. no commit, push, or farm run occurs.

---

## Hard stop

Stop after the memory-policy update, deterministic memory checks, manifest regeneration,
diff audit, and review-bundle creation.

Do not:

- modify scientific source;
- modify testing/runtime tools;
- commit;
- push;
- run the farm;
- begin F.4.Refresh.2 materialization;
- begin F.6.3/E.8.4;
- begin canonical-five expansion;
- create another reconciliation phase.

Return the fresh review bundle for one ChatGPT actual-diff review. If it passes, the user
performs the normal commit/push, ChatGPT performs the lightweight pushed-state
synchronization review, and then the project proceeds directly to the F.4.Refresh.2 farm
materialization gate.
