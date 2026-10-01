# Memory-bracketed scientific throughput

## Decision and reason

Repository memory remains mandatory and authoritative. It brackets substantive
KaonLT work to prevent regression and improve scientific throughput. Repeated
wording/reconciliation cycles must not displace the scientific/runtime gate.
This memory-only checkpoint is based on committed `test` HEAD
`712ba32b062772d87fe44efa6865346e5c438827` under the
[task contract](../phases/memory-bracketed-scientific-throughput-task-contract.md).
It creates no new scientific phase or runtime acceptance.

## Normal loop

```text
memory checkpoint
  -> one meaningful scientific milestone / implementation
  -> deterministic local checks
  -> ChatGPT actual-diff review
  -> user commit/push
  -> ChatGPT pushed-state synchronization review
  -> one narrow user-run farm/runtime gate when required
  -> fresh evidence review
  -> memory checkpoint
  -> continue
```

Actual-diff review, user commit/push, and pushed-state review are required.
They synchronize the source Codex implemented, ChatGPT reviewed, and the farm
will run; they are not independent scientific milestones. Pushed-state review
checks identity, exact gate-relevant changed paths/blobs, and CURRENT/NEXT
continuity. If pushed source matches the reviewed candidate and memory remains
materially accurate, proceed directly to the substantive gate. A push alone
requires no new source-change contract, separate final-pre-push reconciliation,
or broad memory update merely to restate the push.

## Material blockers and batching

A memory inconsistency blocks scientific progress only when it creates
concrete ambiguity about current source identity, accepted evidence, frozen
scientific interfaces, active scientific ownership, or the exact next operation.
Hard integrity failures remain blocking. A warning blocks only for that material
ambiguity or a hard integrity violation. Otherwise record the issue, continue
the active milestone, and correct it at the next checkpoint or milestone audit.
Cosmetic, historical, or wording-only drift cannot start a standalone repair
cycle unless it changes active meaning.

Batch nonblocking corrections at completed scientific implementations, accepted
farm/runtime evidence, real blocker changes to NEXT, production-promotion
decisions, and major E.8/F.6 gate audits. Maintenance must not become an
indefinitely recursive sequence between scientific gates. Retain manifest
integrity, memory-health reporting, and push-stable CURRENT with one NEXT.

## Scientific deliverable and boundaries

Each substantive scientific loop must produce at least one new plot, yield
table, validated numerical comparison, accepted runtime artifact, or one
directly evidenced runtime/scientific blocker with one coherent repair.
Documentation alone is not scientific completion. Test early with a narrow
runtime gate; process cannot substitute for fresh evidence.

The immediate substantive chain is `Q4p4W2p74` F.4.Refresh.2 farm materialize
-> verify -> package, followed only after accepted evidence by
`Q4p4W2p74 / Left / lowe` Method-A reweighting/yield demonstration. That
demonstration must show baseline/reweighted missing-mass spectra and per-t
comparisons, per-(t,phi) baseline and reweighted yields, absolute/fractional
changes, parent-t preservation, and procedure-PDF presentation. Do not expand
to canonical-five or begin unrelated hardening/presentation cleanup before
Left/lowe before/after evidence unless a concrete blocker requires it.

Preserve all accepted scientific statuses, authorities, frozen interfaces,
separate subtraction/diagnostic/presentation owners, and uncertainty handling.
Method B remains diagnostic-only; Method A remains detached until an explicit
validated F.6.4 promotion decision. Source review and synchronization are not
farm, ROOT/PyROOT, full-analysis, or production validation. The user alone
commits/pushes and runs farm gates. This checkpoint runs none and changes no
scientific, runtime, testing, plotting, operational, or authority source.
