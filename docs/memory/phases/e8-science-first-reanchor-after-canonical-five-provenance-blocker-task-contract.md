# KaonLT — science-first memory re-anchor after canonical-five provenance blocker

## 1. Task class and objective

This is a **closure/reconciliation task** under `docs/memory/CODEX.md`.

It is a memory-only checkpoint after a failed canonical-five runtime gate and
subsequent read-only farm diagnosis.

The objective is to:

1. record the second canonical-five full-run failure and the post-run detached
   comparator result as diagnostic runtime evidence;
2. reconcile the live active state so the failed canonical-five provenance gate
   is no longer represented as awaiting an ordinary rerun;
3. mark the narrow canonical-five provenance/identity repair `DEFERRED` by the
   current user decision;
4. return the sole authoritative NEXT to the E.8 scientific audit of the
   baseline kaon missing-mass and stage-yield chain;
5. preserve every accepted narrow runtime closure and all frozen scientific
   ownership.

This task changes no analysis behavior.

It does **not** authorize:

- any Jefferson Lab farm run;
- any analysis/scientific-source edit;
- any owner/checker/profile/collector modification;
- any F.1/F.3/F.4 pin refresh;
- any provenance/identity implementation repair;
- any Method-A production promotion;
- any numerical Method-B dependency;
- any SIMC normalization change;
- final E.8 closure;
- F.6.4 promotion.

Stop after the local memory candidate, deterministic memory checks, and one
complete review bundle are ready for independent ChatGPT actual-diff review.

---

## 2. Exact starting identity and hard precondition

Required branch:

```text
test
```

Required starting local HEAD and local `origin/test`:

```text
ace8688a71431d13b40ed19713a27746f3da6a8e
```

Before editing, establish and report:

```text
git branch --show-current
git rev-parse HEAD
git rev-parse origin/test
git status --short --untracked-files=all
```

Hard stop if any of the following is true:

- branch is not `test`;
- HEAD differs from the required SHA;
- local `origin/test` differs from the required SHA;
- the worktree contains unexpected unrelated state that cannot be preserved;
- this contract is not present at its exact repository path;
- current tracked source/evidence materially contradicts the runtime facts in
  this contract.

Do not reset, clean, stash, overwrite unrelated files, commit, push, update
refs, or run the farm.

Preserve unrelated local state exactly.

---

## 3. Mandatory startup reads

Read in full and in exact repository-required order:

1. `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read:

- `docs/memory/MAINTENANCE.md`
- `docs/memory/CODEX.md`
- `docs/memory/decisions/memory-bracketed-scientific-throughput.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/evidence/VALIDATION_HISTORY.md`
- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/evidence/e8-4-canonical-five-lineage-preflight-runtime-closure-2026-10-05.md`
- `docs/memory/evidence/e8-4-canonical-five-fresh-f1-lineage-materialization-2026-10-05.md`
- `docs/memory/phases/e8-4-canonical-five-fresh-f1-candidate-lineage-refresh-task-contract.md`
- `docs/memory/phases/e8-4-canonical-five-lineage-preflight-runtime-closure-task-contract.md`

Read additional memory only if one of those records directly references it and
it is necessary to resolve a concrete contradiction.

Do not reopen closed phases.

Use repository authority order:

```text
current source/diff
-> newest applicable farm/runtime evidence
-> newest applicable validation/handoff artifact
-> tracked durable memory
-> older chat/history
```

---

## 4. Runtime evidence to record

Create one new canonical diagnostic-evidence record:

```text
docs/memory/evidence/e8-4-canonical-five-second-full-run-provenance-blocker-2026-10-05.md
```

Record only the facts below and their explicit interpretation limits.

### 4.1 Second full canonical-five owner attempt

Farm-evaluated source:

```text
ace8688a71431d13b40ed19713a27746f3da6a8e
```

Owner gate-status basename:

```text
KaonLT_E8_4_Fix5_CanonicalFive_Q4p4W2p74_20261005-202554_1227600-gate-status.json
```

Gate-status SHA-256:

```text
ddfadaf9946e262c35795689ab8d2d52a306ede9e775d98fb05df4e15924e91f
```

Directly supplied/inspected gate facts:

```text
status = failed
stage = verify_artifacts
failure_reason = new_page_missing_or_duplicate
analysis_started = true
analysis_returncode = 0
analysis_completed = true
artifact_verification_completed = false
collection_completed = false
zip_verification_completed = false
```

The normal owner lineage preflight passed before analysis.

The pre-analysis F.4 shared reproduction passed for all five settings.

This does **not** establish canonical-five runtime closure.

### 4.2 Fresh page-manifest diagnosis after the failed owner gate

For all five canonical settings:

```text
Left-lowe
Left-highe
Center-lowe
Center-highe
Right-highe
```

the fresh page manifests showed:

```text
page count = 74
renderer_failures = []
exactly one E.8.4 page = full_background.e8_4.unavailable
```

and the identical unavailable reason:

```text
f6_3_f4_shared_reproduction_failed:f3_fingerprint_input_content_mismatch
```

Therefore the generic owner failure `new_page_missing_or_duplicate` is a
downstream page-inventory consequence of F.6.3 becoming unavailable during the
real analysis. It is not evidence of a duplicate-page renderer defect.

### 4.3 Post-run F.1 identities

The full analysis rewrote all five F.1 artifacts after the successful
pre-analysis lineage preflight.

Record these post-run raw SHA-256 values:

```text
Left-lowe
34e98696fab7807e9c37580d145b762fdfee3e53025c32dba586dabbef6d4f36

Left-highe
3648346330ab01cd52e7d58fb3002d12f758e4f9f14634017ec1e5a02e5b4a9c

Center-lowe
4fd5cc5c691d2f181ffa5134fcf803328e824b7bb94ae04fb8e91a7c8144c724

Center-highe
fc5f340db724a0bd4dd9b2e527dd67003b088a8751bc53cce71956f28652c251

Right-highe
81a09c6dca80cbe5c5dc477ebcece19f297aadced33c7cb30cbf8492e3d7775a
```

Record these corresponding post-run stable F.1 content fingerprints:

```text
Left-lowe
69c1b044abb143a187dd4471efba4ca430b61122204c41d6ad4225ba70d152ba

Left-highe
e492e9ffb50b868855fd548ab99a97ff366ddb1a06fbac81a8be1535c3f8c3fd

Center-lowe
3e3fb1be62074a4e44c408b23674493a2c90f26ff93f6e1cb45bb7f5862ff7cc

Center-highe
4cd7924d2a136e89ca508d3b4d093263fdace20b9cb62a7f540471152a4312e7

Right-highe
2b41d1c72a2b0f3b62f4d4a64de0cc7731e465e3c405d085952d067a2ba6fdb8
```

They differ from the F.1 lineage that passed the accepted pre-analysis
lineage-preflight gate.

Do not infer from the identity change alone that the Method-A scientific payload
changed.

### 4.4 Tracked post-run scientific-equivalence comparator

The existing tracked read-only comparator:

```text
testing/compare_method_a_current_baseline_authority.py
```

was run on the five post-run F.1 artifacts against the reviewed current-baseline
materialization.

The supplied farm output established:

```text
f2_scientific_payload_match = true
f3_scientific_payload_match = true
f4_scientific_payload_match = true
first_changed_stage = none
```

Record this as **RUNTIME VERIFIED diagnostic evidence**.

Interpretation boundaries:

**RUNTIME VERIFIED**

The full run regenerated a different F.1 raw/stable identity set, while the
tracked detached comparator reconstructed F.2, F.3, and F.4 scientific payloads
that match the reviewed current-baseline materialization.

**INFERENCE**

The immediate canonical-five failure is a provenance/identity coupling in the
validation/runtime path rather than evidence, from this gate, of a newly changed
Method-A F.2/F.3/F.4 scientific payload.

**NOT VERIFIED**

- exact low-level cause of the run-to-run F.1 identity change;
- whether a future provenance repair should use scientific-equivalence matching,
  producer canonicalization, or another reviewed design;
- final canonical-five runtime closure;
- canonical-five PDF/visual closure;
- production promotion.

Do not claim the provenance issue is repaired.

---

## 5. Work-state reconciliation

Use only repository-defined work-state labels.

Preserve without downgrade or reopening:

- accepted fresh-F1 separate lineage preflight:
  `CLOSED / RUNTIME VALIDATED` for its exact pre-analysis scope;
- accepted Left/lowe F.6.3/E.8.4 branch execution:
  `CLOSED / RUNTIME VALIDATED` for its existing narrow scope;
- accepted Fix.5.7/Fix.5.8 Left/lowe owner/presentation closures;
- all earlier accepted F-stage closures at their existing scopes.

Reconcile the broader canonical-five item as follows:

```text
canonical-five full runtime gate: BLOCKED
```

Reason:

```text
the second full canonical-five owner attempt completed the analysis child but
failed artifact/page verification because F.6.3 became unavailable after the
analysis regenerated a different F.1 identity set
```

Also record:

```text
canonical-five provenance/identity repair: DEFERRED
```

This `DEFERRED` status is the current user decision to stop another
minor-repair loop and return to the E.8 science objective. It is not a claim
that the blocker is solved or unimportant for later final closure.

Preserve:

```text
E.8: ACTIVE
final E.8 closure: BLOCKED
F.6.4: BLOCKED
absolute-SIMC interpretation: BLOCKED
```

Method A remains detached/non-production.

Method B remains diagnostic/cross-check only and numerically excluded.

---

## 6. Authoritative NEXT

`docs/memory/CURRENT.md` must have one sole ordinary NEXT.

Replace the stale canonical-five farm-readiness/full-run NEXT with a
science-first E.8 NEXT.

Required semantics:

```text
NEXT — resume the E.8.2 scientific audit of the authoritative Q4p4W2p74
baseline kaon missing-mass and stage-yield chain. Trace the existing production
sequence from prompt/random through dummy, slow-proton cleaning, baseline pion
subtraction, final canonical clean-kaon missing mass, and Y0(t,phi), and identify
where the actual spectrum/yield changes occur before deciding whether any new
implementation or farm gate is required.
```

The NEXT must also make clear that:

- the canonical-five provenance repair is `DEFERRED`;
- no farm command is authorized merely by this memory checkpoint;
- no new Method-A design is authorized;
- E.8.2 is science/audit work, not a new maintenance phase;
- after E.8.2, subsequent E.8.3/E.8.4 interpretation follows the already
  approved roadmap unless fresh evidence establishes a blocker.

Commit/push must not be the sole NEXT.

---

## 7. Exact allowed paths

Expected tracked modifications:

```text
docs/memory/CURRENT.md
docs/memory/roadmap/STATUS.md
docs/memory/evidence/VALIDATION_HISTORY.md
docs/memory/manifest.json
```

Expected new files:

```text
docs/memory/evidence/e8-4-canonical-five-second-full-run-provenance-blocker-2026-10-05.md
docs/memory/phases/e8-science-first-reanchor-after-canonical-five-provenance-blocker-task-contract.md
```

No other file is authorized.

In particular, `docs/memory/MEMORY.md` is **not** allowlisted. Its durable
cross-phase rules remain applicable and this task must not restate ephemeral
active status there.

If a material contradiction is discovered that cannot be corrected inside this
allowlist, stop and report it. Do not broaden scope automatically.

---

## 8. Frozen files and directories

Everything outside the allowlist is frozen.

In particular, do not modify:

```text
src/
testing/
tools/
farm_env/
background_samples/
run_Prod_Analysis.sh
set_SymLinks.sh
AGENTS.md
docs/memory/MEMORY.md
docs/memory/USER.md
docs/memory/handoffs/CURRENT_HANDOFF.md
docs/memory/CODEX.md
docs/memory/MAINTENANCE.md
docs/memory/COMMUNICATION.md
docs/memory/LEARNINGS.md
docs/memory/decisions/
docs/memory/investigations/
```

Do not modify any existing accepted evidence record.

Do not change candidate identities, authority maps, profiles, owner source,
collectors, renderers, tests, or launchers.

---

## 9. Scientific ownership and preserved runtime path

This task changes no physics and no executable behavior.

Preserve exactly:

- random subtraction;
- dummy subtraction;
- slow-proton subtraction;
- baseline pion subtraction;
- active `no_empirical_residual` profile;
- pion templates, components, fits, windows, amplitudes, and priors;
- Method-A F.1/F.2/F.3/F.4/F.5 mathematics;
- F.6.3 private-branch weighting semantics;
- Method-B diagnostic/cross-check-only role and numerical exclusion;
- SIMC generation/normalization;
- yields and uncertainty propagation;
- cuts and canonical binning;
- efficiencies and acceptance;
- L/T separation and cross sections.

Preserve the normal scientific/runtime production ordering:

```text
prompt/random
-> dummy
-> slow-proton cleaning
-> baseline pion subtraction
-> final baseline canonical missing mass / yield extraction
```

Method A remains detached unless and until a later explicit validated F.6.4
decision promotes it.

This memory checkpoint does not alter or validate that runtime path.

---

## 10. Required file-specific edits

### 10.1 `docs/memory/CURRENT.md`

Keep E.8 `ACTIVE`.

Replace stale wording that says canonical-five full validation remains merely
farm-validation pending or that another full run is the ordinary NEXT.

Record concisely:

- accepted separate lineage-preflight closure and its narrow scope;
- second full canonical-five attempt completed analysis but failed artifact/page
  verification;
- all five E.8.4 pages were unavailable because of
  `f3_fingerprint_input_content_mismatch`;
- post-run tracked comparator reconstructed matching F.2/F.3/F.4 scientific
  payloads with `first_changed_stage = none`;
- no new Method-A scientific defect is established by this failed gate;
- canonical-five full runtime gate is `BLOCKED`;
- canonical-five provenance/identity repair is `DEFERRED`;
- final E.8/F.6.4 remain `BLOCKED`;
- absolute-SIMC blocker remains separate;
- sole NEXT is the E.8.2 baseline missing-mass/stage-yield scientific audit.

Keep CURRENT concise and below its soft size warning if practical. Consolidate
obsolete pre-run prose rather than appending indefinitely.

### 10.2 `docs/memory/roadmap/STATUS.md`

Preserve historical chronology and all accepted closures.

Update only the current canonical-five dependency/status surface.

The `Current-lineage canonical-five isolated orchestration` section must no
longer say only:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

as if the second full gate had not happened.

Record that:

- implementation/source readiness remains historical context;
- accepted preflight remains narrowly closed;
- the second full runtime gate is now `BLOCKED` by post-analysis F.1
  provenance/identity mismatch;
- tracked post-run comparator found exact F.2/F.3/F.4 scientific equality;
- the provenance repair is `DEFERRED` while science-first E.8 work resumes;
- final E.8 and F.6.4 remain `BLOCKED`.

Update the dependency prose so it does not require immediate provenance repair
before the active E.8.2 scientific audit.

Do not give the roadmap a competing sole NEXT.

### 10.3 `docs/memory/evidence/VALIDATION_HISTORY.md`

Keep this file evidence/history-oriented.

Add a concise entry for the failed second canonical-five full run and link the
new evidence record.

Correct current-state wording that still says canonical-five expansion is simply
`DEFERRED` from the older user decision or otherwise contradicts the later
user-authorized orchestration history.

The correct current distinction is:

```text
canonical-five expansion/orchestration was user-authorized;
the second full runtime gate is BLOCKED;
the narrow provenance/identity repair is now DEFERRED by the current user
decision while E.8 science resumes.
```

Do not turn VALIDATION_HISTORY into a second CURRENT.

End current-state prose by deferring live objective/NEXT ownership to
`CURRENT.md`.

### 10.4 New evidence record

Create the file specified in Section 4.

It must contain:

- source identity;
- gate-status basename and SHA;
- analysis-completed/owner-failed distinction;
- five-setting page-manifest diagnosis;
- post-run five F.1 raw hashes;
- post-run five stable F.1 fingerprints;
- tracked comparator result;
- exact evidence labels;
- explicit non-closures;
- Method-A/Method-B boundaries;
- no personal filesystem paths.

Do not claim accepted canonical-five runtime validation.

### 10.5 `docs/memory/manifest.json`

Regenerate using the repository tool after the final intended memory contents
are complete.

Do not hand-edit hashes or byte counts.

---

## 11. Positive checks

Confirm all of the following:

1. `CURRENT.md` contains exactly one ordinary live NEXT.
2. That NEXT is the E.8.2 baseline missing-mass/stage-yield scientific audit.
3. Another canonical-five farm run is not the NEXT.
4. The second full canonical-five owner attempt is recorded as `BLOCKED`, not
   accepted runtime closure.
5. Analysis child success is distinguished from owner/artifact-gate failure.
6. The accepted separate lineage preflight remains
   `CLOSED / RUNTIME VALIDATED` only for its narrow scope.
7. The post-run comparator result is recorded exactly:
   F.2=true, F.3=true, F.4=true, first changed stage=none.
8. No new Method-A scientific defect is claimed from the failed gate.
9. Canonical-five provenance/identity repair is `DEFERRED`.
10. E.8 remains `ACTIVE`.
11. Final E.8 remains `BLOCKED`.
12. F.6.4 remains `BLOCKED`.
13. Absolute-SIMC interpretation remains separately `BLOCKED`.
14. Method A remains detached/non-production.
15. Method B remains diagnostic/cross-check only and numerically excluded.
16. All accepted Left/lowe and earlier F-stage closures retain their exact
    scopes.
17. No farm run is authorized by this checkpoint.
18. No scientific/executable source changes.

---

## 12. Negative and regression checks

Audit the final candidate for contradictory live wording equivalent to:

```text
full canonical-five validation remains FARM VALIDATION PENDING
one full Q4p4W2p74 canonical-five farm run is NEXT
canonical-five runtime is validated
post-run F.1 change proves Method-A science changed
F.4 is the first changed stage in the second full run
Method A is production
Method B contributes numerically
absolute SIMC is resolved
F.6.4 is unblocked
provenance repair is complete
```

Historical statements may remain where clearly historical and correctly scoped.

Also search for stale current wording that says canonical-five expansion itself
is still deferred under the older decision without acknowledging that the user
later authorized the isolated canonical-five orchestration.

The final current-state distinction must be unambiguous:

```text
orchestration attempted
-> second full runtime gate BLOCKED
-> provenance/identity repair DEFERRED
-> science-first E.8.2 audit NEXT
```

Run:

```text
git diff --check
```

No executable/scientific path may appear in the diff.

---

## 13. Forbidden shortcuts and fallbacks

Do not:

- change F.1/F.3/F.4 pins;
- weaken or bypass lineage validation;
- add a scientific-equivalence fallback to executable source;
- canonicalize F.1 producers;
- rerun or redesign Method A;
- change accepted materialization;
- create a new farm owner;
- issue or execute a farm command;
- reinterpret the failed owner artifact as accepted closure;
- downgrade accepted narrow closures;
- use chat memory in place of the supplied/runtime facts;
- broaden the allowlist because another memory file contains merely historical
  wording;
- modify `MEMORY.md` merely to repeat CURRENT.

Any executable repair belongs to a later standalone implementation contract if
and when the user returns to the deferred provenance blocker.

---

## 14. Local validation and memory health

Discover the working local Python interpreter using `docs/memory/TOOLS.md`.

After all versioned memory changes:

1. regenerate the manifest:

```text
<PYTHON> -B tools/update_memory_manifest.py --root . --write
```

2. verify it:

```text
<PYTHON> -B tools/update_memory_manifest.py --root . --check
```

3. run ordinary memory health:

```text
<PYTHON> -B tools/check_memory_health.py --root .
```

4. run bootstrap:

```text
<PYTHON> -B tools/memory_bootstrap.py --root . --json
```

5. run:

```text
git diff --check
```

Do not use `--fail-on-warning` unless an existing warning is found to be
materially blocking under `MAINTENANCE.md`. Record nonblocking warnings and do
not create another maintenance cycle for them.

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

No executable tests are required for this memory-only closure unless a memory
tool itself requires them.

---

## 15. Diff audit and review bundle

Before stopping, run:

```text
git status --short --untracked-files=all
git diff --stat
git diff --check
```

Confirm:

- no path outside the allowlist changed;
- no scientific/runtime/test/owner/profile/collector source changed;
- new evidence contains only supplied/established diagnostic facts;
- CURRENT remains concise and push-stable;
- roadmap and validation history do not compete with CURRENT;
- manifest was regenerated after final contents.

Create one complete temporary review bundle in repository root:

```text
kaonlt_review.diff
```

It must contain:

- complete tracked diff for every changed tracked path;
- complete `git diff --no-index /dev/null ...` addition for both intended new
  files;
- no unrelated user-owned content;
- no staging merely for review.

Stop before commit/push.

---

## 16. Acceptance criteria

PASS only if all are true:

1. starting branch/HEAD/origin/worktree gate passes;
2. only the exact allowlisted memory files change;
3. second full-run failed-gate evidence is recorded accurately;
4. the accepted separate preflight remains narrowly closed;
5. post-run comparator result is recorded exactly;
6. canonical-five full runtime gate is `BLOCKED`;
7. provenance/identity repair is `DEFERRED`;
8. no executable repair is made;
9. E.8 remains `ACTIVE`;
10. final E.8 and F.6.4 remain `BLOCKED`;
11. Method A remains detached/non-production;
12. Method B remains diagnostic-only/numerically excluded;
13. absolute-SIMC blocker remains explicit;
14. sole CURRENT NEXT is the E.8.2 scientific audit;
15. roadmap and validation history reflect that state without taking NEXT
    ownership;
16. manifest regeneration/check passes;
17. ordinary memory health passes;
18. bootstrap passes;
19. `git diff --check` passes;
20. complete `kaonlt_review.diff` is prepared;
21. no commit, push, remote-ref update, or farm run occurs.

---

## 17. Hard stop

After memory edits, consistency audit, manifest regeneration, health/bootstrap
checks, and `kaonlt_review.diff` creation:

**STOP.**

Do not commit.

Do not push.

Do not run the farm.

Return:

- exact changed paths;
- concise reconciliation summary;
- exact memory-health report;
- exact warnings/classification;
- manifest/bootstrap results;
- `git diff --check` result;
- path to `kaonlt_review.diff`.

ChatGPT must review the actual diff before any user-controlled commit/push.
