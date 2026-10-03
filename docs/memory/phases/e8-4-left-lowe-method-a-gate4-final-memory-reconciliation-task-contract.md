# KaonLT — Gate 4 final memory audit/reconciliation after Method-A scientific audit

## 1. Task class and objective

This is a **closure/reconciliation task** under `docs/memory/CODEX.md`.

It is the fourth and final gate in the approved pre-implementation sequence:

```text
Gate 1 — initial memory re-anchor
-> Gate 2 — complete repository-memory consistency audit
-> Gate 3 — complete code/science audit
-> Gate 4 — final memory audit/reconciliation
-> only then — next scientific-direction decision before any implementation contract
```

This task is **Gate 4 only**.

The objective is to reconcile repository memory with the completed Gate-2
repository-memory consistency audit and completed Gate-3 code/science audit,
preserve the Gate-3 scientific findings with correct evidence labels, remove
consumed live-gate wording, compact `CURRENT.md`, and establish one accurate,
push-stable substantive NEXT.

Gate 4 does **not** decide the eventual pion-response implementation and does
**not** perform external literature research. After Gate 4 is reviewed, pushed,
and synchronized, the next step is a scientific-direction checkpoint to decide
whether to perform an external analogous-experiment/methods review before any
implementation contract.

Gate 4 authorizes:

- memory-only edits within the explicit allowlist;
- one new Gate-3 source/science investigation record;
- deterministic repository-memory checks;
- one complete review bundle for independent ChatGPT actual-diff review.

Gate 4 does **not** authorize:

- any Jefferson Lab farm run;
- any scientific-source edit;
- any change to pion/slow-proton/random/dummy/SIMC/yield/cross-section behavior;
- any Method-A production promotion;
- any numerical Method-B dependency;
- any pion probability implementation;
- any external literature/web research;
- canonical-five expansion;
- final E.8 closure;
- F.6.4 closure;
- any implementation contract.

Stop after the local candidate and review bundle are ready.

---

## 2. Exact starting source and hard precondition

Required branch:

```text
test
```

Required starting local HEAD and local `origin/test`:

```text
3a3a4217c78dde07123e71cd36e43ee54666d1e4
```

Before editing, run and report:

```text
git branch --show-current
git rev-parse HEAD
git rev-parse origin/test
git status --short --untracked-files=all
```

Hard stop if any of the following is true:

- branch is not `test`;
- `HEAD` differs from the required SHA;
- local `origin/test` differs from the required SHA;
- the worktree contains unexpected unrelated state that cannot be preserved safely;
- this contract file is absent from the exact path below;
- current tracked source/evidence materially contradicts the Gate-2 or Gate-3
  findings recorded in this contract.

Do not reset, clean, stash, overwrite, commit, push, update refs, or run the farm.

Preserve unrelated local state exactly.

---

## 3. Required startup reads

Read in the repository-mandated order:

1. `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read:

- `docs/memory/MAINTENANCE.md`
- `docs/memory/CODEX.md`

Then read only the records needed for this Gate-4 reconciliation, including at minimum:

- `docs/memory/roadmap/STATUS.md`
- `docs/memory/evidence/VALIDATION_HISTORY.md`
- `docs/memory/phases/e8-4-left-lowe-method-a-scientific-validity-memory-reanchor-task-contract.md`
- `docs/memory/evidence/f1-fix5-v4-bundle-inspection-2026-09-14.md`
- `docs/memory/evidence/f2-fix1-runtime-closure.md`
- `docs/memory/evidence/f3-runtime-closure.md`
- `docs/memory/evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md`
- `docs/memory/evidence/f6-1-runtime-closure.md`
- `docs/memory/evidence/f6-2-scientific-runtime-closure.md`
- `docs/memory/evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md`
- `docs/memory/evidence/e8-4-fix5-8-left-lowe-runtime-closure-2026-10-03.md`
- `docs/memory/investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md`

Inspect current source only as needed to confirm that the Gate-3 findings still
match the required HEAD. Do not repeat the complete Gate-3 code/science audit.

Use the repository authority order:

```text
current source/diff
-> newest applicable farm/runtime evidence
-> newest applicable validation/handoff artifact
-> tracked durable memory
-> older chat/history
```

---

## 4. Gate-2 repository-memory audit findings to reconcile

Gate 2 completed as a read-only repository-memory consistency audit and found
**no blocking inconsistency** in current scientific ownership, accepted runtime
evidence, Method-A/Method-B boundaries, absolute-SIMC separation, or the approved
pre-implementation sequence.

Gate 2 identified the following nonblocking reconciliation set:

1. `CURRENT.md` still contains consumed Gate-1/Gate-2 live wording.
2. `roadmap/STATUS.md` still contains consumed Gate-1/Gate-2 transition wording.
3. `docs/memory/investigations/KNOWN_GAPS.md` contains an obsolete historical
   `NEXT — Begin detached F.2 ...`.
4. `docs/memory/investigations/2026-09-11-phase-f1-source-identity-reconciliation.md`
   contains an obsolete historical F.1 farm NEXT.
5. `docs/memory/phases/phase-f6-method-a-production-promotion.md` contains stale
   downstream chronology for F.4.Refresh.2 / F.6.3 / E.8.4.
6. `docs/memory/phases/PHASE_HISTORY.md` contains an older F.6.3/E.8.4 frontier.
7. `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md` contains stale
   historical Fix.5 presentation status.
8. Older F.6.2 records contain obsolete then-current downstream “E.8 is NEXT”
   wording.
9. `CURRENT.md` was 8152 bytes at the Gate-2 audit, close to the 8-KiB soft
   warning threshold and should be compacted rather than expanded.

Gate 2 also established that not every stale historical sentence should be
rewritten. Historical records that are clearly chronology, dated investigation,
accepted evidence, or immutable imported provenance do not acquire active-state
authority merely because their historical NEXT is old.

Gate 4 must therefore:

- repair **live/canonical ambiguity**;
- compact `CURRENT.md`;
- update the active dependency/status surface;
- preserve historical records as historical unless a minimal explicit
  superseded/consumed marker is needed to avoid current ambiguity;
- avoid broad history rewriting.

---

## 5. Gate-3 SOURCE VERIFIED scientific findings to preserve

Gate 3 completed as a read-only code/science audit. It introduced no new farm
runtime evidence and changed no source.

### 5.1 Method-A population ownership

Current source and accepted F.1/F.2/F.3 structure preserve distinct populations:

```text
training:
  prompt / noRF / nommcuts / P_hgcer_npeSum > 0

training low-response class:
  0 < P_hgcer_npeSum <= 2

training control-response class:
  P_hgcer_npeSum > 2

physical application:
  authoritative physical-pion-control population with P_hgcer_npeSum > 2
```

Do not merge these populations or describe the application population as the
training population.

### 5.2 Accepted F.3 scientific meaning

The accepted F.3 basis is:

```text
hgcer3:
  SHMS_delta
  P_hgcer_xAtCer
  P_hgcer_yAtCer
```

F.3 constructs a support-aware **relative detector-response representation**.
It explicitly does **not** construct:

- an absolute pion leakage probability;
- an absolute zero-photoelectron probability;
- a production correction;
- a production-side weight.

Do not rename F.3 as a probability map.

### 5.3 F.4 parent-preserving mathematics

For physical application events the relevant conceptual structure is:

```text
b_j = signed_source_coefficient_j * w0_j
r_j = relative Method-A response
B   = sum_j b_j
U   = sum_j b_j * r_j
C_j = r_j / (U / B)
```

within each canonical-`t` parent.

Therefore:

```text
sum_j b_j * C_j = sum_j b_j
```

by construction.

This signed parent-preserving equality is a constraint. It does **not** require:

- a `(t,phi)` child to be preserved;
- a missing-mass subregion to be preserved;
- a Lambda-window yield to be preserved;
- an individual event contribution to remain near baseline.

No independent child renormalization is permitted.

### 5.4 Baseline `w0` scientific ownership

The existing baseline pion weight `w0` already represents the established
pion-control-to-kaon-background transfer derived from the accepted pion
component model.

A future HGCer response quantity cannot automatically be treated as another
independent contamination normalization or multiplied into `w0` without first
establishing distinct scientific ownership and showing that no normalization
or transfer physics is double-counted.

This is a durable design constraint, not an approved replacement algorithm.

### 5.5 Current-baseline lineage

Applicable accepted comparator evidence established:

```text
F.2 scientific payload: exact match
F.3 scientific payload: exact match
first scientifically changed stage: F.4
```

for the current-baseline refresh comparison.

Therefore the current Left/lowe scientific question is concentrated at the
interaction between:

```text
accepted relative HGCer response
x
current physical application population
x
baseline w0 pion transfer
x
signed F.4 parent normalization
```

Do not claim that F.2/F.3 themselves changed scientifically on this evidence.

### 5.6 HGCer-hole interpretation

The configured HGCer geometric hole is already excluded from the relevant
diagnostic/application populations.

Therefore the current low-`t` redistribution must not be described simply as
“correcting events inside the known HGCer hole.”

Detector-response dependence elsewhere in the HGCer plane remains possible.

No particular mirror, PMT, optical alignment, track-geometry, or hardware cause
is established.

### 5.7 Signed-normalization amplification is a hypothesis

Because F.4 uses signed baseline contributions:

```text
B = sum_j b_j
U = sum_j b_j * r_j
normalization = U / B
```

substantial positive/negative cancellation could make the normalized correction
more sensitive than a positive-measure average.

This is `INFERENCE`, not an accepted cause.

The current-lineage `t1` decomposition needed to test it remains `NOT VERIFIED`,
including at minimum:

```text
sum b_j
sum |b_j|
|sum b_j| / sum |b_j|
sum b_j r_j
source-separated signed support
raw r_j distribution
final C_j distribution
(t,phi) redistribution
missing-mass-region redistribution
response-coordinate dependence
support/OOD dependence
```

Do not write that signed cancellation explains the observed `t1` effect.

### 5.8 Existing zero-photoelectron pion diagnostic

Current source already contains a detached zero-photoelectron HGCer transfer
diagnostic using zero-truncated detector-response models to infer a control-to-
zero-photoelectron response transfer.

It remains:

```text
non-authoritative
diagnostic-only
production-side-effect-free
```

It predates the later accepted `hgcer3` program and is not automatically a
replacement for Method A.

Do not promote, copy, or apply it in Gate 4.

### 5.9 Slow-proton analogue

The slow-proton architecture is useful only as an **architectural analogue**:

```text
event probability
-> proposed effect
-> support/preservation diagnostic
-> applied-effect decision
```

The pion problem is statistically different because the target
zero-photoelectron pion population is absent from the positive-response pion
control sample.

Do not transplant the slow-proton numerical model into pion subtraction.

---

## 6. Applicable RUNTIME VERIFIED boundaries to preserve

Gate 4 must preserve, without widening:

- accepted F.1/F.2/F.3/F.4 detached scientific/runtime closures at their recorded
  scopes;
- accepted F.6.1/F.6.2 detached validation scopes;
- `F.4.Refresh.2` as `CLOSED / RUNTIME VALIDATED` at its accepted scope;
- F.6.3 / E.8.4 `CLOSED / RUNTIME VALIDATED` only for
  `Q4p4W2p74 / Left / lowe` branch execution, current-lineage application,
  real child changes, and signed parent preservation;
- Fix.5.7 and Fix.5.8 `CLOSED / RUNTIME VALIDATED` only for their accepted
  Q4p4W2p74 / Left / lowe owner/checker/presentation scopes.

No Gate-3 source conclusion upgrades any scientific/runtime scope.

No new runtime evidence exists from Gates 2–4.

---

## 7. Scientific blocker after Gate 4

The live unresolved scientific issue remains:

```text
Why does the current Q4p4W2p74 / Left / lowe Method-A candidate produce a large
low-t / t1 redistribution and downstream pion-background/yield impact, and is
that redistribution supported by a physically appropriate detector-response
model rather than by an inappropriate normalization or transfer construction?
```

Bookkeeping consistency is not the blocker.

No accepted cause currently exists.

The exact current-lineage diagnostics still needed include:

```text
signed vs absolute parent support
source-sign/source-class decomposition
raw relative-response tails
final correction-factor tails
(t,phi) redistribution
missing-mass-region redistribution
hgcer3 coordinate dependence
support/OOD dependence
comparison to a physically interpretable zero-photoelectron response transfer
```

Gate 4 does not implement these diagnostics.

---

## 8. Global frozen boundaries

Preserve:

```text
no_empirical_residual
zero legacy empirical residual scales
baseline production authoritative
random/dummy subtraction frozen
slow-proton subtraction frozen
baseline pion subtraction frozen
pion templates/components/fits/windows/amplitudes frozen
Method A detached/non-production
Method B diagnostic/cross-check only and numerically excluded
SIMC production/normalization frozen
yield formulas frozen
uncertainty propagation frozen
binning frozen
efficiencies frozen
acceptance frozen
L/T separation frozen
cross sections frozen
```

Preserve:

```text
canonical-five expansion: DEFERRED
final E.8: BLOCKED
F.6.4: BLOCKED
absolute-SIMC interpretation: BLOCKED as a separate provenance/units issue
```

No farm run is authorized by Gate 4.

---

## 9. Allowed versioned files

The intended versioned changes are limited to:

```text
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/roadmap/STATUS.md
docs/memory/evidence/VALIDATION_HISTORY.md
docs/memory/manifest.json
docs/memory/phases/e8-4-left-lowe-method-a-gate4-final-memory-reconciliation-task-contract.md
docs/memory/investigations/e8-4-left-lowe-method-a-detector-response-source-science-audit-2026-10-03.md
```

The last two paths are intended new tracked files.

`MEMORY.md` is allowlisted only for compact durable scientific lessons from
Gate 3. Do not duplicate the full Gate-3 investigation there.

Do not edit the Gate-1 contract.

---

## 10. Frozen files and directories

Everything outside the allowlist is frozen.

In particular, do not modify:

```text
src/
testing/
tools/
farm_env/
run_Prod_Analysis.sh
AGENTS.md
docs/memory/USER.md
docs/memory/handoffs/CURRENT_HANDOFF.md
docs/memory/CODEX.md
docs/memory/MAINTENANCE.md
docs/memory/COMMUNICATION.md
docs/memory/LEARNINGS.md
docs/memory/decisions/
docs/memory/import/
all existing runtime evidence files
all existing task contracts and historical phase records
```

The Gate-2 stale-history findings outside the allowlist must be preserved as
reported historical drift unless they create active-state ambiguity after the
canonical/live records are corrected.

Do not broaden Gate 4 into a history-rewrite exercise.

---

## 11. Exact required changes

### 11.1 `docs/memory/CURRENT.md`

Compact rather than append.

Keep E.8 `ACTIVE`.

Replace consumed Gate-1/Gate-2 wording with the current state:

- Gates 1, 2, and 3 are consumed;
- Gate 2 found no blocking repository-memory contradiction and recorded
  nonblocking historical/canonical drift for this closing checkpoint;
- Gate 3 completed a source/science audit and supplied no new farm evidence;
- stale scalar / wrong histogram is not the primary scientific issue;
- the current scientific blocker is the physical validity and detector-response
  origin of the large low-`t` / `t1` Method-A redistribution;
- F.2/F.3 science remained unchanged in the current-baseline comparison and F.4
  was the first changed scientific stage;
- signed parent-normalization amplification is a hypothesis, not an accepted
  explanation;
- no particular PMT/mirror/track/hardware cause is established;
- the existing zero-photoelectron pion transfer machinery is diagnostic-only;
- the slow-proton approach is analogue-only;
- Method A remains detached/non-production;
- Method B remains diagnostic/cross-check only and numerically excluded;
- canonical-five expansion remains `DEFERRED`;
- final E.8/F.6.4 remain `BLOCKED`;
- absolute-SIMC interpretation remains a separate blocker.

Keep `CURRENT.md` below the 8-KiB soft-warning threshold by consolidating stale
source-history detail and linking canonical records rather than repeating them.

The sole authoritative NEXT must be push-stable and must be:

```text
NEXT — after Gate-4 actual-diff review, user commit/push, and pushed-state
synchronization: scientific-direction checkpoint to decide whether an external
analogous-experiment / PID-background methods review should be performed before
any implementation contract.
```

Do not make commit/push itself the NEXT.

Do not authorize the external research automatically in this contract.

### 11.2 `docs/memory/MEMORY.md`

Add only compact durable scientific lessons from Gate 3:

1. A relative detector-response representation and an absolute leakage
   probability are different scientific objects.
2. When an established control-to-background transfer weight already owns
   normalization, a second response/transfer object cannot be multiplied in
   without proving distinct ownership and avoiding double counting.
3. Signed parent-preserving normalization is not statistically equivalent to a
   positive probability-measure normalization and must be decomposed explicitly
   when cancellation can matter.
4. Zero-truncated/censored detector-response inference is relevant when the
   population of interest is absent from the selected control sample.
5. The slow-proton proposed/applied architecture is a useful structural analogue
   but not a numerical pion prescription.
6. For the current Left/lowe case, F.2/F.3 science remained unchanged in the
   accepted current-baseline comparison and the first scientific change appeared
   at F.4.

Preserve all accepted scientific ownership.

Do not perform broad MEMORY consolidation except where needed to remove direct
duplication created by these additions.

### 11.3 `docs/memory/roadmap/STATUS.md`

Preserve all accepted closures and scopes.

Replace the consumed live dependency entry with:

```text
accepted Left/lowe branch/runtime evidence
-> Gate 1 initial memory re-anchor: consumed
-> Gate 2 repository-memory consistency audit: consumed
-> Gate 3 code/science audit: consumed
-> Gate 4 final memory reconciliation: current closing checkpoint
-> scientific-direction checkpoint
-> implementation contract only if later explicitly warranted
-> later reconsider canonical-five expansion
-> final E.8
-> F.6.4 explicit production-promotion decision
```

Do not declare external literature research mandatory yet.

Keep:

```text
canonical-five expansion: DEFERRED
final E.8: BLOCKED
F.6.4: BLOCKED
```

Do not create a new production phase or promotion status.

### 11.4 `docs/memory/evidence/VALIDATION_HISTORY.md`

Preserve accepted farm evidence exactly.

Do not add a new `RUNTIME VERIFIED` claim from Gate 3.

Make only the minimum update needed so the live limitation section reflects:

- Gates 2 and 3 have been completed;
- accepted narrow runtime closures remain unchanged;
- Gate 3 did not establish a detector cause;
- Gate 3 did not validate Method A for production;
- absolute-SIMC provenance/units remain separate;
- current scientific interpretation remains blocked pending further evidence.

The detailed Gate-3 source/science findings belong in the new investigation
record, not in runtime evidence history.

### 11.5 New Gate-3 investigation record

Create:

```text
docs/memory/investigations/e8-4-left-lowe-method-a-detector-response-source-science-audit-2026-10-03.md
```

This record must preserve the detailed Gate-3 audit and explicitly separate:

```text
SOURCE VERIFIED
RUNTIME VERIFIED
MEMORY ONLY
INFERENCE
NOT VERIFIED
```

At minimum it must include:

- Method-A training/application population distinction;
- F.3 `hgcer3` relative-response meaning;
- F.4 signed parent-preserving algebra;
- baseline `w0` ownership;
- accepted current-baseline F.2/F.3 equality and F.4 first-change result;
- configured HGCer-hole exclusion;
- signed-normalization amplification as a hypothesis only;
- existing detached zero-photoelectron pion transfer diagnostic;
- slow-proton architectural analogy;
- current unresolved `t1` scientific questions;
- exact current-lineage diagnostic quantities still required.

It must not prescribe an implementation.

It may end by stating that after Gate 4 the project will hold a
scientific-direction checkpoint before any implementation contract.

### 11.6 `docs/memory/manifest.json`

Regenerate using the repository's existing manifest procedure after the final
intended memory tree is complete.

Do not hand-edit hashes or byte counts.

### 11.7 Track this Gate-4 contract

Track this exact contract at:

```text
docs/memory/phases/e8-4-left-lowe-method-a-gate4-final-memory-reconciliation-task-contract.md
```

and include it in the regenerated manifest.

---

## 12. Before/after behavior

### Before

At starting HEAD `3a3a4217...`:

- Gate 1 is still described as `ACTIVE` in `CURRENT.md`;
- Gate 2 is still described as the sole NEXT;
- roadmap status still says later audits have not begun;
- Gate-3 findings exist only in the completed source/science review, not in a
  durable detailed investigation record;
- `CURRENT.md` is near its soft size threshold.

### After

After Gate 4:

- repository memory accurately records Gates 1–3 as consumed and Gate 4 as the
  closing checkpoint;
- Gate-3 source/science findings are durable and evidence-labeled;
- stale scalar/wrong histogram is not the live blocker;
- the unresolved scientific problem is the physical validity/detector-response
  origin of the large Left/lowe low-`t`/`t1` redistribution;
- signed cancellation remains a hypothesis;
- Method A remains detached;
- Method B remains numerically excluded;
- accepted runtime closures remain unchanged;
- absolute-SIMC remains a separate blocker;
- canonical-five remains `DEFERRED`;
- final E.8/F.6.4 remain `BLOCKED`;
- `CURRENT.md` is compact and below its warning threshold;
- the sole substantive NEXT is a post-Gate-4 scientific-direction checkpoint,
  not an implementation contract and not an automatically authorized literature
  review.

---

## 13. Positive checks

Confirm all of the following:

```text
CURRENT has one active objective
CURRENT has one push-stable substantive NEXT
CURRENT is below the 8-KiB soft-warning threshold
Gate 1 consumed
Gate 2 consumed
Gate 3 consumed
Gate 4 is the current closing checkpoint
accepted runtime statuses preserved
Method A remains detached/non-production
Method B remains diagnostic-only and numerically excluded
absolute-SIMC blocker remains separate
canonical-five remains DEFERRED
final E.8 remains BLOCKED
F.6.4 remains BLOCKED
Gate-3 investigation has explicit evidence labels
signed-normalization explanation remains INFERENCE
specific detector-hardware cause remains NOT VERIFIED
no external-research decision is silently promoted
```

---

## 14. Negative checks

Confirm all of the following:

```text
no src changes
no testing changes
no tool/runtime changes
no farm command
no new probability implementation
no production correction
no Method-B numerical dependency
no new child renormalization
no claimed PMT/mirror/hardware cause
no claimed signed-cancellation explanation
no new runtime evidence
no historical evidence rewrite
no implementation contract
no external literature/web research
```

Search the intended live/canonical memory surfaces for stale active wording such
as:

```text
Gate 1 is ACTIVE
Gate 2 is the only next executable gate
Gate 2 has not begun
later audits have not begun
```

Historical task contracts, dated investigations, and chronology may retain their
original wording when clearly historical and non-authoritative.

---

## 15. Local validation and memory health

Discover `<PYTHON>` as described in `docs/memory/TOOLS.md`.

After all intended memory edits:

1. regenerate the memory manifest using the existing repository procedure;
2. check the manifest;
3. run:

```text
<PYTHON> -B tools/check_memory_health.py --root .
```

Because Gate 4 is the final memory/milestone audit and owns warning elimination,
also run:

```text
<PYTHON> -B tools/check_memory_health.py --root . --fail-on-warning
```

Gate 4 requires zero hard failures and zero warnings.

Report:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
warning classification: <none or classifications>
manifest check: PASS | FAIL
ordinary health: PASS | FAIL
--fail-on-warning: PASS | FAIL
```

Run:

```text
git diff --check
```

No ROOT/PyROOT/full-`main.py`/farm/PDF claim follows from these checks.

---

## 16. Complete actual-diff review bundle

Create one temporary repository-root file:

```text
kaonlt_review.diff
```

It must contain:

- complete tracked diffs for every modified tracked file;
- complete `git diff --no-index /dev/null ...` addition for this new Gate-4
  contract;
- complete `git diff --no-index /dev/null ...` addition for the new Gate-3
  investigation record.

Do not stage merely for review.

Do not include unrelated local changes.

Stop before commit/push and return the review bundle for independent ChatGPT
actual-diff review.

---

## 17. Acceptance criteria

Gate-4 local candidate passes only if:

- exact starting branch/HEAD/origin/worktree preconditions pass;
- only allowlisted paths change;
- consumed Gate-1/Gate-2/Gate-3 live wording is reconciled;
- Gate-3 findings are preserved with correct evidence labels;
- `CURRENT.md` is compact and below 8 KiB;
- accepted runtime evidence and scientific ownership remain unchanged;
- no inference is promoted to a verified cause;
- the next step remains a scientific-direction checkpoint before implementation;
- manifest check passes;
- ordinary health passes;
- `--fail-on-warning` passes;
- `git diff --check` passes;
- the complete review bundle is produced;
- no farm, commit, push, scientific-source change, literature research, or
  implementation occurs.

---

## 18. Hard stop

Return `BLOCKED` and stop if:

- repository identity differs;
- unexpected worktree state cannot be preserved safely;
- current source contradicts the Gate-3 findings;
- accepted runtime evidence would need to be downgraded;
- a scientific-source modification appears necessary;
- `CURRENT.md` cannot be brought below the milestone warning threshold within
  the allowed scope without losing necessary active-state information;
- memory health has a hard failure or warning that cannot be resolved within the
  allowlisted memory scope;
- exact live NEXT cannot be stated without inventing a scientific decision.

Do not invent a workaround.

Do not proceed automatically to literature research or implementation.
