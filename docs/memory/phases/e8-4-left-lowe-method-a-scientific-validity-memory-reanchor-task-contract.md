# KaonLT — Gate 1 initial memory re-anchor for Method-A low-t scientific validity

## 1. Task class and objective

This is a **closure/reconciliation task** under `docs/memory/CODEX.md`.

It is the first gate in the user-approved pre-implementation sequence:

```text
Gate 1 — initial memory re-anchor
-> Gate 2 — complete repository-memory consistency audit
-> Gate 3 — complete code/science audit
-> Gate 4 — final memory audit/reconciliation
-> only then — implementation contract
```

This task is **Gate 1 only**.

The objective is to re-anchor authoritative repository memory to the current scientific question without performing the later audits and without changing analysis behavior.

The current scientific focus is no longer the earlier apparent `MM_0 -> Y0` / `MM_A -> YA` stale-scalar or wrong-histogram bookkeeping hypothesis as the main issue. The re-anchor must instead make the live unresolved issue:

```text
the scientific validity and detector-response origin of the large low-t
Method-A redistribution, including the strong global-weight differences and
large t1 pion-background / yield impact
```

The slow-proton event-probability treatment may be recorded as a useful **scientific/architectural analogue** for thinking about a future rigorous pion HGCer response / contamination-probability framework. It is not an approved pion probability implementation, production correction, or promotion decision.

This task must preserve all accepted upstream/narrow runtime closures and all frozen production ownership.

This task authorizes:

- memory-only edits within the explicit allowlist;
- deterministic repository-memory checks;
- creation of one complete review bundle for ChatGPT actual-diff review.

This task does **not** authorize:

- any Jefferson Lab farm run;
- any scientific-source edit;
- any change to pion/slow-proton/random/dummy/SIMC/yield/cross-section behavior;
- any Method-A production promotion;
- any numerical Method-B dependency;
- any pion probability implementation;
- canonical-five expansion;
- any E.8/F.6.4 closure.

Stop after the local candidate and review bundle are ready for independent ChatGPT actual-diff review.

---

## 2. Exact starting source and hard precondition

Required branch:

```text
test
```

Required starting local HEAD and local `origin/test`:

```text
222af982c567078b8ddf804b5762cb42cc991e32
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
- the contract file itself is not present at the exact path given below;
- current tracked source/evidence materially contradicts the scientific re-anchor described in this contract.

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

Then read only the records needed for this Gate 1 re-anchor, including at minimum:

- `docs/memory/roadmap/STATUS.md`
- `docs/memory/evidence/VALIDATION_HISTORY.md`
- `docs/memory/evidence/e8-4-fix5-8-left-lowe-runtime-closure-2026-10-03.md`
- `docs/memory/evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md`
- `docs/memory/investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md`
- the immediately preceding memory re-anchor contract:
  `docs/memory/phases/e8-4-left-lowe-mm-yield-closure-memory-reanchor-task-contract.md`

Do **not** perform the complete repository-memory audit in Gate 1. Read additional records only when directly necessary to resolve a contradiction in the active-state re-anchor.

Use the repository authority order from `AGENTS.md`:

```text
current source/diff
-> newest applicable farm/runtime evidence
-> newest applicable validation/handoff artifact
-> tracked durable memory
-> older chat/history
```

---

## 4. Re-anchor facts and boundaries to preserve

### 4.1 Accepted narrow closures remain accepted

Do not downgrade or reopen the accepted narrow scopes already recorded for:

- Fix.5.7 owner/checker setting-provenance repair;
- Fix.5.8 Left/lowe presentation-legibility repair;
- the existing Q4p4W2p74 / Left / lowe F.6.3/E.8.4 branch-execution evidence;
- parent-preserving Method-A behavior already accepted for that narrow branch.

Preserve their exact work-state labels and runtime scope unless current tracked evidence directly contradicts them.

### 4.2 Earlier bookkeeping concern is not the main scientific issue

The previous live memory at starting HEAD `222af982...` restores the apparent
`MM_0 -> Y0` / `MM_A -> YA` numerical-closure question as the active gate.

Gate 1 must re-examine only enough tracked evidence/current source context to decide whether that wording is now stale.

The intended re-anchor is:

- the source/runtime path is internally fail-closed;
- the previously suspected stale-scalar / wrong-histogram identity failure is not the main unresolved scientific problem;
- accepted evidence and prior arithmetic checks support histogram-to-scalar bookkeeping consistency strongly enough that the active question should no longer be framed as “which object is wrong?”;
- do **not** overstate this into a new production-physics validation;
- do **not** claim that every detector-response interpretation is established;
- do **not** conflate arithmetic consistency with scientific validity of the redistribution.

If current source/evidence cannot support that distinction, stop and report the exact contradiction instead of forcing the new wording.

### 4.3 The unresolved scientific question

The active unresolved question after Gate 1 must be the **scientific validity and detector origin of the large low-t Method-A redistribution**, especially:

- strong global-weight differences;
- large t1 pion-background redistribution;
- large t1 downstream yield impact;
- whether the behavior is physically attributable to HGCer response / pion contamination structure rather than an analysis-bookkeeping defect;
- what detector-response variables/populations would be scientifically appropriate for diagnosing that behavior.

No cause is accepted in Gate 1.

Do not conclude:

- that HGCer inefficiency is the cause;
- that a particular PMT/mirror/track geometry is the cause;
- that a contamination probability is already defined;
- that Method A is correct or incorrect as a production correction.

Those belong to later audit/implementation gates.

### 4.4 Slow-proton treatment is an analogue only

Record, where appropriate, that the established slow-proton event-probability treatment is a useful scientific/architectural analogue for thinking about a rigorous pion HGCer-response / contamination-probability framework.

This means only that a future framework may need to distinguish:

- detector-response probability from parent normalization;
- event-level response from child-bin redistribution;
- training/control populations from the physical application population;
- detector geometry/kinematics from purely global reweighting.

It does **not** authorize copying the slow-proton implementation, changing proton subtraction, or creating a pion probability map in this task.

### 4.5 Frozen global boundaries

Preserve:

- `no_empirical_residual`;
- zero legacy residual scales;
- baseline production as authoritative;
- random/dummy subtraction;
- slow-proton subtraction;
- baseline pion subtraction;
- pion templates/components/fits/windows/amplitudes;
- Method-A detached/non-production status;
- Method-B diagnostic/cross-check-only status and numerical exclusion;
- SIMC production/normalization;
- yield formulas and uncertainty propagation;
- binning;
- efficiencies;
- acceptance;
- L/T separation;
- cross sections.

Preserve the **absolute-SIMC issue as a separate blocker**.

Preserve:

```text
canonical-five expansion: DEFERRED
final E.8: BLOCKED
F.6.4: BLOCKED
```

No farm run is authorized by this re-anchor.

---

## 5. Exact pre-implementation sequence to record

Repository memory after Gate 1 must preserve this exact sequence:

```text
Gate 1 — initial memory re-anchor
-> Gate 2 — complete repository-memory consistency audit
-> Gate 3 — complete code/science audit
-> Gate 4 — final memory audit/reconciliation
-> only then — implementation contract
```

The sequence is intentionally conservative.

During **all four pre-implementation gates**:

- no Jefferson Lab farm run;
- no scientific-source modification;
- no production behavior change;
- no Method-A promotion;
- no numerical Method-B use;
- no canonical-five expansion;
- no pion-probability implementation.

Gate 2 is the only next executable gate after Gate 1 passes.

Do not skip directly to the code/science audit or an implementation contract.

---

## 6. Allowed versioned files

The intended versioned changes are limited to:

```text
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/roadmap/STATUS.md
docs/memory/evidence/VALIDATION_HISTORY.md
docs/memory/manifest.json
docs/memory/phases/e8-4-left-lowe-method-a-scientific-validity-memory-reanchor-task-contract.md
```

`MEMORY.md` is allowlisted only for the **smallest direct correction needed to remove a material contradiction** with the new active-state re-anchor.

Do not use Gate 1 to perform broad consolidation or the complete repository-memory consistency audit. If other stale records are discovered, report them as Gate-2 audit targets unless they create a hard active-state ambiguity that must be corrected now.

---

## 7. Frozen files and directories

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
docs/memory/investigations/
docs/memory/evidence/e8-4-fix5-8-left-lowe-runtime-closure-2026-10-03.md
docs/memory/evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md
```

Do not rewrite historical evidence merely to fit the new active interpretation.

If historical wording is stale but correctly historical, leave it and let Gate 2 classify it.

---

## 8. Exact required changes

### 8.1 `docs/memory/CURRENT.md`

Keep E.8 as `ACTIVE`.

Replace the current live work item that frames the primary issue as unresolved
`MM_0 -> Y0` / `MM_A -> YA` object/scalar closure.

The new **Current Work Item** must state that:

- the previous apparent bookkeeping inconsistency has been re-examined;
- no stale-scalar / wrong-histogram identity defect is currently established as the main issue;
- the source/runtime path is fail-closed;
- accepted arithmetic/evidence supports treating histogram-to-scalar bookkeeping as internally consistent for the present purpose;
- this does **not** validate the scientific origin or production correctness of the Method-A redistribution;
- the unresolved issue is the large low-t Method-A redistribution, especially the strong global-weight differences and t1 pion-background/yield impact;
- detector-response science/HGCer origin is the next substantive subject of investigation;
- the slow-proton event-probability framework is only an analogue for a possible future rigorous pion-response framework;
- no pion probability implementation is approved;
- Method A remains detached/non-production;
- Method B remains diagnostic/cross-check only and numerically excluded;
- canonical-five expansion remains `DEFERRED`;
- final E.8/F.6.4 remain `BLOCKED`;
- absolute-SIMC interpretation remains a separate blocker.

The **Blockers** section must distinguish:

1. scientific validity / detector-response origin of the large low-t Method-A redistribution;
2. absolute-SIMC provenance/units as a separate blocker.

Do not retain stale wording that makes stale scalar/wrong histogram identity the primary active blocker unless the current evidence forces a hard stop.

The sole authoritative **NEXT** must be:

```text
Gate 2 — complete repository-memory consistency audit
```

and it must state that Gate 2 occurs before the code/science audit and before any implementation contract.

The NEXT must also state that Gates 1–4 authorize no farm run or scientific-source modification.

Do not make commit/push itself the NEXT.

### 8.2 `docs/memory/MEMORY.md`

Perform only the minimum direct correction needed so durable memory does not materially contradict CURRENT.

Specifically audit the paragraph(s) that currently require histogram-scalar closure as the live unresolved issue.

Update them only as needed to preserve the distinction:

```text
arithmetic/bookkeeping consistency
!=
scientific validity / detector-response explanation
```

Retain the reusable rule that histogram/scalar closure must be checked when interpreting large signed shifts.

Do not erase that general lesson.

Instead, make clear that for the present Left/lowe work item, the stale-scalar/wrong-histogram hypothesis is no longer the main unresolved question; the next science issue is detector-response validity of the redistribution.

Do not perform broad MEMORY consolidation. Any other inconsistency belongs to Gate 2 unless it creates immediate active-state ambiguity.

### 8.3 `docs/memory/roadmap/STATUS.md`

Preserve all accepted phase closures.

Replace or supersede the live dependency entry that makes the Left/lowe histogram-to-yield numerical-closure gate the immediate active dependency.

The roadmap dependency should become:

```text
accepted Left/lowe branch execution/presentation and arithmetic bookkeeping
-> Gate 1 memory re-anchor
-> Gate 2 complete repository-memory consistency audit
-> Gate 3 complete code/science audit of Method-A redistribution and detector-response origin
-> Gate 4 final memory audit/reconciliation
-> only then implementation contract if warranted
-> later reconsider canonical-five expansion
-> final E.8
-> F.6.4 explicit production-promotion decision
```

Keep canonical-five expansion `DEFERRED`.

Keep final E.8 and F.6.4 `BLOCKED`.

Do not create a new production phase or promotion status.

### 8.4 `docs/memory/evidence/VALIDATION_HISTORY.md`

Preserve all accepted evidence and narrow runtime scopes.

Remove the live-history wording that says the accepted Fix.5.7/Fix.5.8 state necessarily leaves the stale-scalar/wrong-histogram numerical closure as the current active gate.

Replace it with a narrow current limitation:

- the accepted evidence supports arithmetic/bookkeeping integrity strongly enough that stale scalar / wrong histogram is not the live primary scientific issue;
- those runtime closures still do **not** establish the detector-response origin or scientific validity of the large low-t Method-A redistribution;
- absolute-SIMC units remain separate;
- no Method-A promotion follows.

Do not invent new farm evidence.

Do not change the accepted scope of Fix.5.7 or Fix.5.8.

### 8.5 `docs/memory/manifest.json`

Regenerate using the existing repository procedure after the final intended memory tree is complete.

Do not hand-edit hashes or byte counts.

### 8.6 Track this contract

Place this exact contract at:

```text
docs/memory/phases/e8-4-left-lowe-method-a-scientific-validity-memory-reanchor-task-contract.md
```

and include it in the regenerated manifest.

---

## 9. Before/after behavior

### Before

At starting HEAD `222af982...`:

- CURRENT makes the Left/lowe `MM_0 -> Y0` / `MM_A -> YA` numerical-closure question the active gate;
- roadmap/validation history carry the same live framing;
- the next substantive step is a source/runtime-path audit of histogram-to-scalar provenance.

### After

After Gate 1:

- repository memory no longer treats stale scalar / wrong histogram identity as the main unresolved scientific problem;
- the active unresolved issue is the scientific validity and detector-response origin of the large low-t Method-A redistribution;
- the slow-proton probability treatment is recorded only as an analogue;
- the four pre-implementation gates are explicit and ordered;
- Gate 2, the complete repository-memory consistency audit, is the sole NEXT;
- no code/science audit has yet been performed;
- no implementation is authorized;
- no farm run is authorized;
- production behavior is unchanged.

---

## 10. Positive checks

Confirm all of the following:

1. `CURRENT.md` names scientific validity / detector-response origin of the large low-t Method-A redistribution as the active unresolved issue.
2. `CURRENT.md` does not present stale scalar / wrong histogram identity as the main active blocker.
3. The distinction between arithmetic consistency and scientific validity is explicit.
4. The strong low-t/t1 redistribution and global-weight behavior are preserved as unresolved scientific concerns.
5. The slow-proton event-probability treatment appears only as an analogue.
6. No pion probability implementation is authorized.
7. Method A remains detached/non-production.
8. Method B remains diagnostic/cross-check only and numerically excluded.
9. Canonical-five expansion remains `DEFERRED`.
10. Final E.8 and F.6.4 remain `BLOCKED`.
11. Absolute-SIMC remains a separate blocker.
12. The exact four-gate pre-implementation sequence is recorded.
13. Gate 2 is the sole authoritative NEXT.
14. No farm run is authorized.
15. No scientific-source change is authorized.
16. Accepted Fix.5.7/Fix.5.8 and F.6.3/E.8.4 narrow runtime closures are preserved.

---

## 11. Negative and regression checks

Search the final candidate memory state for live wording that incorrectly keeps the stale bookkeeping hypothesis primary.

At minimum inspect live/current occurrences of:

```text
MM_0 -> Y0
MM_A -> YA
histogram-to-scalar
numerical closure
which object is wrong
signed cancellation remains unproven
```

Historical evidence may retain such wording when clearly historical, but it must not compete with the new CURRENT/NEXT.

Also verify that no new text claims:

```text
Method A is production
Method A is scientifically validated
Method B contributes numerically
pion probability is implemented
slow-proton code is reused for pion correction
HGCer is proven to cause the redistribution
canonical-five is active
F.6.4 is unblocked
absolute SIMC is resolved
a farm run is authorized
code/science audit has already passed
implementation may begin before Gate 4
```

Run:

```text
git diff --check
```

No executable/scientific path may appear in the diff.

---

## 12. Preserved scientific/runtime path

This Gate 1 task changes memory only.

The scientific/runtime path must remain byte-for-byte untouched.

Preserve:

```text
random/dummy
-> frozen binning
-> slow-proton treatment
-> baseline pion subtraction
-> detached Method-A diagnostics/private branch where already implemented
-> authoritative yields / E.8 consumers
```

No producer, serializer, checkpoint, cache, histogram, scalar yield, detector cut, response map, weight, template, SIMC normalization, or renderer behavior may change.

---

## 13. Local validation

Discover the repository-selected Python interpreter using existing `TOOLS.md` guidance.

After all versioned memory changes:

1. regenerate `docs/memory/manifest.json`;
2. verify the manifest;
3. run the ordinary memory-health gate:

```text
<PYTHON> -B tools/check_memory_health.py --root .
```

Because this is an explicit memory re-anchor correcting material active-state/NEXT drift, also run:

```text
<PYTHON> -B tools/check_memory_health.py --root . --fail-on-warning
```

Report exactly:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
warning classification: <none or each warning: blocking | nonblocking>
manifest check: PASS | FAIL
```

A hard failure or a warning that leaves current source identity, accepted evidence, scientific ownership, blockers, or NEXT ambiguous is a hard stop.

Do not run ROOT/PyROOT, `main.py`, the farm, or any production analysis.

---

## 14. Actual-diff review bundle

Do not stage merely for review.

Produce one complete temporary review bundle in the repository root:

```text
kaonlt_review.diff
```

It must include:

- the complete tracked diff for every modified tracked allowlisted file;
- a complete `git diff --no-index /dev/null ...` addition for this new contract file if it remains untracked at review time.

The review bundle itself is temporary and must not be tracked.

Return:

- starting branch/HEAD/origin/test;
- exact changed paths;
- exact unchanged/frozen-path confirmation;
- memory-health report;
- stale-live-wording audit summary;
- the complete `kaonlt_review.diff`.

Stop before commit/push.

ChatGPT must review the **actual diff**, not only the Codex summary.

---

## 15. Farm-validation boundary

There is no farm-validation step in Gate 1.

The user-approved four pre-implementation gates explicitly prohibit a farm run.

No farm command may be supplied or executed during this task.

A future farm gate can only be reconsidered after:

```text
Gate 1 PASS
-> user commit/push
-> ChatGPT pushed-state synchronization
-> Gate 2 complete
-> Gate 3 complete
-> Gate 4 complete
-> reviewed implementation contract, if implementation is actually warranted
```

---

## 16. Acceptance criteria

Gate 1 passes only if:

- exact starting branch/HEAD/origin/test match the contract;
- the worktree is safe to modify without disturbing unrelated state;
- only allowlisted memory files and this contract change;
- current memory no longer makes stale scalar/wrong histogram identity the main unresolved scientific issue;
- the active issue becomes the scientific validity and detector origin of the large low-t Method-A redistribution;
- arithmetic/bookkeeping consistency is distinguished from scientific validity without inventing new farm evidence;
- slow-proton probability treatment is analogue-only;
- Method A remains detached/non-production;
- Method B remains diagnostic-only and numerically excluded;
- canonical-five remains `DEFERRED`;
- final E.8/F.6.4 remain `BLOCKED`;
- absolute-SIMC remains separate;
- the four-gate pre-implementation sequence is explicit;
- Gate 2 is the sole substantive NEXT;
- no farm run or scientific-source change is authorized;
- manifest verification passes;
- ordinary and zero-warning memory health pass;
- `git diff --check` passes;
- independent ChatGPT actual-diff review has not yet been bypassed.

---

## 17. Hard stop

Return `BLOCKED` and make no further edits if:

- branch/HEAD/origin/test do not match;
- the worktree has unexpected unrelated changes that cannot be safely preserved;
- current tracked source/evidence materially contradicts the intended bookkeeping/scientific-validity distinction;
- accepted runtime status would have to be downgraded or invented to perform the re-anchor;
- any scientific/executable source change appears necessary;
- the task would require implementing pion probability or changing Method-A/Method-B roles;
- manifest/health hard checks fail;
- the diff escapes the allowlist.

Do not invent a workaround.

Do not proceed to Gate 2 automatically.

Stop for independent ChatGPT actual-diff review.
