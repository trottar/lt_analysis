# KaonLT — E.8.4 Left/lowe MM-versus-yield numerical-closure memory re-anchor

## 1. Task class and objective

This is a **closure/reconciliation task** under `docs/memory/CODEX.md`.

It is memory-only. Do **not** modify scientific, analysis, renderer, owner,
collector, launcher, profile, test, farm-execution, or production source.

The objective is to correct the live KaonLT memory state after a user-identified
workflow/state mismatch:

1. preserve the accepted narrow Fix.5.7/Fix.5.8 runtime closures exactly;
2. preserve the accepted narrow F.6.3/E.8.4 Left/lowe branch-execution and
   presentation evidence exactly;
3. restore the unresolved `Q4p4W2p74 / Left / lowe` histogram-to-scalar-yield
   consistency question as the **active** E.8.4 gate;
4. make clear that presentation/legibility acceptance did **not** validate the
   numerical closure between the displayed `MM_0`/`MM_A` spectra and stored
   `Y0`/`YA` values;
5. keep canonical-five expansion `DEFERRED` and prevent any broadening from
   becoming the next executable gate until the Left/lowe numerical closure is
   resolved;
6. keep final E.8 and F.6.4 `BLOCKED`;
7. preserve Method A as detached/non-production and Method B as diagnostic-only
   and numerically excluded;
8. preserve the absolute-SIMC normalization/unit blocker as a separate issue;
9. regenerate/check the memory manifest and run memory health;
10. stop for independent ChatGPT actual-diff review before user commit/push.

This task authorizes no farm run and no scientific-source change.

---

## 2. Exact starting source

Required branch:

```text
test
```

Required starting local HEAD and local `origin/test`:

```text
5addf32433d61e94391fd31dc454b3628b026735
```

Before editing, verify:

```text
git branch --show-current
git rev-parse HEAD
git rev-parse origin/test
git status --short --untracked-files=all
```

Hard stop if branch, HEAD, or local `origin/test` differs.

Preserve unrelated local state. Do not reset, clean, stash, overwrite, commit,
push, or run the farm.

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
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/evidence/VALIDATION_HISTORY.md`
- `docs/memory/investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md`
- `docs/memory/evidence/e8-4-fix5-8-left-lowe-runtime-closure-2026-10-03.md`
- `docs/memory/evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md`

Use current source/diff and applicable runtime evidence over stale summaries.

---

## 4. Evidence that must be preserved

The following accepted facts are not being reopened.

### Fix.5.8 presentation closure

`E.8.4 Fix.5.8` remains:

```text
CLOSED / RUNTIME VALIDATED
```

only for `Q4p4W2p74 / Left / lowe` presentation legibility.

Accepted farm source:

```text
2ddeab47d55edb57d2f313022a948c4376730c19
```

Accepted ZIP:

```text
KaonLT_E8_4_Fix5_Left_lowe_Q4p4W2p74_20261003-032934-987132.zip
```

Accepted ZIP SHA-256:

```text
966e667b36b2626b5c16f00fc684a099a5d24e5d2015bb547612a1e63eb2a25e
```

The 97-page structural/visual result, no-renderer-failure result, page-legibility
repairs, parent-preservation values, and non-production flags remain unchanged.

### Fix.5.7 owner/checker closure

Fix.5.7 remains narrowly:

```text
CLOSED / RUNTIME VALIDATED
```

for the Left/lowe owner/checker setting-provenance repair.

### Existing F.6.3/E.8.4 narrow runtime evidence

Do not erase or downgrade the accepted runtime evidence that the Left/lowe
current-baseline Method-A branch executes, changes real canonical child yields,
and preserves the signed parent constraint.

That accepted execution evidence does **not** by itself establish the unresolved
histogram-to-scalar-yield numerical closure described below.

---

## 5. Unresolved numerical inconsistency that must become live again

The existing post-farm investigation records a direct farm observation:

- baseline and Method-A missing-mass spectra are visually near-overlapping in
  several populated children;
- the stored scalar yields differ substantially in the same cells.

Two recorded examples are:

```text
t1, phi [-180,-140):
  Y0 = 0.033906
  YA = 0.025165

t1, phi [140,180):
  Y0 = 0.029873
  YA = 0.016411
```

These values are observations, not proof that either the histogram or scalar is
wrong.

The investigation explicitly leaves the following numerical gate unresolved:

```text
producer-owned integral(MM_0) == stored Y0
producer-owned integral(MM_A) == stored YA
```

for every populated canonical child, using the exact producer-owned integration
window, bin semantics, flow treatment, sign convention, and runtime object
identity.

It also requires, where relevant, exact reporting of:

- signed integral;
- positive-bin support;
- negative-bin support;
- absolute support.

Signed cancellation remains a hypothesis until those exact runtime quantities
close the discrepancy. Do not promote cancellation to an explanation merely
because stored support diagnostics exist.

The confirmed historical E.8.3/F.6.1 versus current E.8.4/F.6.3 lineage
difference remains real but is not sufficient to explain the within-E.8.4
`MM_0/MM_A` versus `Y0/YA` inconsistency.

The absolute-SIMC normalization/unit problem is a separate blocker and must not
be conflated with this data-side histogram-to-yield closure.

---

## 6. Allowed files

The intended versioned changes are limited to:

```text
docs/memory/CURRENT.md
docs/memory/roadmap/STATUS.md
docs/memory/evidence/VALIDATION_HISTORY.md
docs/memory/manifest.json
docs/memory/phases/e8-4-left-lowe-mm-yield-closure-memory-reanchor-task-contract.md
```

`docs/memory/MEMORY.md` already contains the durable rule requiring exact
histogram-scalar closure and treating cancellation as a hypothesis. Do not edit
it unless an actual contradiction is found during the audit. If such a
contradiction is found, stop and report it rather than silently broadening scope.

The existing post-farm investigation is the evidence-bearing record for the
observation and required audit. Do not rewrite its historical chronology merely
to restate current status.

---

## 7. Frozen files

Everything outside the allowlist is frozen, including all of:

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
docs/memory/MEMORY.md
docs/memory/investigations/
docs/memory/evidence/e8-4-fix5-8-left-lowe-runtime-closure-2026-10-03.md
docs/memory/evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md
```

unless this contract explicitly says otherwise.

No numerical, scientific, production, renderer, checker, collector, owner,
profile, launcher, or test behavior may change.

---

## 8. Exact required memory correction

### 8.1 `docs/memory/CURRENT.md`

Keep the active objective as E.8.

Replace the consumed canonical-five decision as the immediate live work item
with an explicit active numerical-closure gate for:

```text
Q4p4W2p74 / Left / lowe
```

The current work item must state, using exact repository work-state labels, that:

- the Left/lowe MM-versus-yield numerical-closure investigation is `ACTIVE`;
- Fix.5.7 and Fix.5.8 remain `CLOSED / RUNTIME VALIDATED` in their narrow scopes;
- accepted F.6.3/E.8.4 branch-execution evidence remains accepted;
- that execution/presentation evidence did not establish
  `MM_0 -> Y0` and `MM_A -> YA` numerical closure;
- the two recorded t1 examples are preserved as observations, not conclusions;
- signed cancellation is not yet an accepted explanation;
- the absolute-SIMC blocker remains separate;
- canonical-five expansion remains `DEFERRED`;
- final E.8/F.6.4 remain `BLOCKED`;
- Method A remains detached/non-production;
- Method B remains diagnostic/cross-check only and numerically excluded.

The sole authoritative `NEXT` must be push-stable and substantive. It should
direct the next work to a dedicated source/runtime-path audit of the exact
Left/lowe histogram-to-scalar-yield provenance and producer semantics before any
canonical-five broadening.

Do **not** make commit/push itself the NEXT.

No farm command is authorized by CURRENT after this memory repair.

### 8.2 `docs/memory/roadmap/STATUS.md`

Preserve all accepted closures.

Add a narrow live dependency/status entry under the existing E.8.4/F.6.3 area
for the Left/lowe histogram-to-yield numerical-closure gate.

Use an exact work-state label:

```text
ACTIVE
```

Do not invent a production promotion or claim a scientific explanation.

Make the dependency explicit:

```text
accepted Left/lowe branch execution/presentation
    -> Left/lowe MM-versus-yield numerical closure
    -> only then may canonical-five expansion be reconsidered
    -> final E.8 closure
    -> F.6.4 explicit production-promotion decision
```

Canonical-five expansion itself remains `DEFERRED` by user decision.

Final E.8 and F.6.4 remain `BLOCKED`.

Do not renumber or reopen historical phases unless the existing roadmap
structure requires a label for this gate. Prefer a descriptive gate title over
inventing a new Fix number.

### 8.3 `docs/memory/evidence/VALIDATION_HISTORY.md`

Preserve all accepted evidence status.

Add the narrow limitation that the accepted Fix.5.8 Left/lowe
presentation/owner closure does **not** establish the unresolved
histogram-to-scalar-yield numerical closure.

State that the post-farm investigation remains the supporting record for the
direct visual/scalar observation.

Do not downgrade Fix.5.7 or Fix.5.8.

Do not claim the cause is known.

### 8.4 `docs/memory/manifest.json`

Regenerate from the final intended memory tree using the existing repository
procedure. Do not hand-edit hashes or byte counts.

### 8.5 Track this contract

Place this exact contract at:

```text
docs/memory/phases/e8-4-left-lowe-mm-yield-closure-memory-reanchor-task-contract.md
```

and include it in the regenerated manifest.

---

## 9. Scientific ownership and boundaries

This task changes memory only.

Preserve:

- baseline production as authoritative;
- `no_empirical_residual`;
- all random/dummy/slow-proton subtraction;
- pion component models/fits/windows/amplitudes;
- Method-A factors and candidate identities;
- Method-B diagnostic-only role;
- canonical binning;
- yield formulas and uncertainty propagation;
- SIMC production/normalization;
- efficiencies and acceptance;
- L/T separation and cross sections.

Do not infer a correction from the visual discrepancy.

Do not alter Method-A promotion status.

Do not explain the discrepancy as cancellation unless later direct runtime
evidence proves it.

---

## 10. Positive checks

After editing, confirm all of the following:

1. `CURRENT.md` names the Left/lowe MM-versus-yield numerical closure as the
   active gate.
2. `CURRENT.md` no longer presents the canonical-five deferral decision as the
   immediate substantive next action.
3. Fix.5.7 remains narrowly `CLOSED / RUNTIME VALIDATED`.
4. Fix.5.8 remains narrowly `CLOSED / RUNTIME VALIDATED`.
5. Existing accepted F.6.3/E.8.4 Left/lowe branch-execution evidence is retained.
6. The unresolved closure is clearly distinguished from:
   - the historical E.8.3/F.6.1 lineage mismatch;
   - the absolute-SIMC unit blocker;
   - presentation/legibility acceptance.
7. Signed cancellation remains explicitly unproven.
8. Canonical-five expansion remains `DEFERRED`.
9. Final E.8/F.6.4 remain `BLOCKED`.
10. Method A remains detached/non-production.
11. Method B remains diagnostic-only and numerically excluded.
12. No farm run is authorized.

---

## 11. Negative and regression checks

Search the final candidate memory state for stale or misleading live wording.

At minimum audit for:

```text
resolve the existing user-deferred canonical-five
select the next approved roadmap item
canonical-five E.8/F.6.3 expansion decision
```

Historical occurrences may remain only when clearly historical and not
competing with CURRENT.

Also verify that no new text claims any of the following:

```text
cancellation explains the yield shift
MM/Y closure passed
SIMC normalization is wrong
Method A is production
Method B contributes numerically
canonical-five is now active
F.6.4 is unblocked
```

Preserve every accepted artifact hash, farm source identity, and runtime scope.

Run:

```text
git diff --check
```

---

## 12. Memory health

Discover the repository-selected Python interpreter using existing `TOOLS.md`
guidance.

After all versionable memory changes:

1. regenerate `docs/memory/manifest.json`;
2. verify the manifest;
3. run the ordinary health gate:

```text
<PYTHON> -B tools/check_memory_health.py --root .
```

4. because this is a memory re-anchor correcting material active-state/NEXT
   ambiguity, also run the zero-warning audit:

```text
<PYTHON> -B tools/check_memory_health.py --root . --fail-on-warning
```

Report:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
warning classification: <none or blocking/nonblocking>
manifest check: PASS | FAIL
```

Any hard failure or warning that leaves source identity, accepted evidence,
scientific ownership, blocker, or NEXT ambiguous is a hard stop.

---

## 13. Actual-diff review bundle

Do not stage merely for review.

Produce one complete temporary review bundle in the repository root:

```text
kaonlt_review.diff
```

It must contain the complete tracked diff for all modified tracked allowlisted
files and a complete `git diff --no-index /dev/null ...` addition for this new
contract.

The review bundle itself is temporary and must not be tracked.

Return the exact changed paths, memory-health report, stale-reference audit
summary, and upload the complete `kaonlt_review.diff` for independent ChatGPT
review.

Stop before commit/push.

---

## 14. Acceptance criteria

This memory re-anchor passes only if:

- live source still starts from the required HEAD;
- no scientific/executable source changed;
- the accepted narrow runtime closures remain intact;
- the unresolved Left/lowe MM-versus-yield numerical consistency issue is
  restored as the active gate;
- CURRENT owns one exact substantive NEXT centered on that closure;
- canonical-five remains `DEFERRED`;
- final E.8/F.6.4 remain `BLOCKED`;
- the absolute-SIMC issue remains separate;
- cancellation remains unproven;
- manifest and both memory-health gates pass;
- `git diff --check` passes;
- independent ChatGPT actual-diff review is still pending.

---

## 15. Hard stop

Stop and report `BLOCKED` without editing if:

- branch/HEAD/local `origin/test` do not match the required starting state;
- there is an unrelated tracked modification;
- the existing investigation no longer contains the direct MM-versus-yield
  observation or unresolved closure requirement;
- preserving accepted Fix.5.7/Fix.5.8/F.6.3/E.8.4 evidence would require
  downgrading or rewriting scientific/runtime history;
- a required change falls outside the allowlist;
- manifest or memory health fails materially.

Do not invent a replacement procedure. Do not run the farm. Do not commit or
push.
