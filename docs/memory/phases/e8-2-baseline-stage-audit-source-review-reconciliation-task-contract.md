# KaonLT E.8.2 — Post-review source-status reconciliation

## Purpose

Reconcile durable repository memory after independent ChatGPT actual-diff review
of the complete cumulative E.8.2 baseline-stage-audit implementation, including
Fix.1, Fix.2, and Fix.3.

ChatGPT has now inspected the complete cumulative review bundle and finds the
E.8.2 source implementation **PASS** for source review.

This task is **memory-only**. It must not alter the reviewed E.8.2 source,
focused tests, previous E.8.2 task contracts, accepted E.8.1/F.6.2 evidence,
production physics, or any farm/runtime machinery.

No ROOT/PyROOT, full `main.py`, procedure-PDF farm rendering, or production
runtime validation is claimed by this reconciliation.

---

## Exact committed base and current worktree

Required branch:

`test`

Required committed HEAD:

`5b2ba7f92dc51b303257dbf209df687415364602`

Commit message:

`Correct E8 no-empirical-residual chain`

The existing uncommitted E.8.2 + Fix.1 + Fix.2 + Fix.3 worktree is intentional
and is the source-reviewed candidate. Do not recreate it.

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
git diff --name-only
git diff --check
```

Required conditions:

- branch is exactly `test`;
- HEAD is exactly
  `5b2ba7f92dc51b303257dbf209df687415364602`;
- the dirty worktree contains only the already reviewed cumulative E.8.2
  candidate, this reconciliation contract, and temporary review artifact(s);
- there are no unrelated user changes.

Do not reset, stash, clean, discard, checkout over, or recreate the existing
candidate.

Do not commit, push, switch branches, update remote refs, or run the Jefferson
Lab farm.

Hard stop if branch/HEAD or worktree scope differs.

---

## Required startup reading

Read repository memory in the normal required order:

1. `docs/memory/AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read:

- `docs/memory/CODEX.md`
- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/decisions/e8-no-empirical-residual-chain-correction.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/e8-2-baseline-stage-audit-task-contract.md`
- `docs/memory/phases/e8-2-baseline-stage-audit-fix1-task-contract.md`
- `docs/memory/phases/e8-2-baseline-stage-audit-fix2-task-contract.md`
- `docs/memory/phases/e8-2-baseline-stage-audit-fix3-task-contract.md`
- `docs/memory/phases/e8-2-baseline-stage-audit.md`
- this reconciliation contract.

Then inspect the complete current cumulative diff and confirm the reviewed
source/test files have not changed since the supplied Fix.3 review bundle.

---

# Independent ChatGPT actual-diff result

ChatGPT independently reviewed the complete cumulative E.8.2 + Fix.1 + Fix.2 +
Fix.3 diff anchored to committed HEAD
`5b2ba7f92dc51b303257dbf209df687415364602`.

Result:

`PASS`

The reviewed implementation preserves the accepted baseline production branch
and adds only detached E.8.2 observation/presentation state.

## Reviewed baseline path

The accepted source-reviewed E.8.2 chain is:

```text
same authoritative event traversal
    -> true pre-proton wide-MM auxiliary capture
    -> existing random subtraction
    -> existing dummy subtraction
    -> existing slow-proton event-factor removal
    -> actual post-proton / pre-prune production snapshot
    -> unchanged production prune_hist treatment
    -> actual post-prune production snapshot
    -> exact accepted baseline pion input
    -> exact B_pi^0 built from accepted w0
    -> final baseline K_0 / MM_0(t,phi)
    -> existing Y_0(t,phi)
```

## Reviewed scientific ownership

The actual diff review found:

- no second scientific tree traversal;
- ordinary production fills remain the existing proton-factor-weighted fills;
- true pre-proton diagnostics are captured in the existing authoritative
  traversal;
- random/dummy normalization ownership is unchanged;
- slow-proton closure is against the actual post-proton/pre-prune production
  snapshot;
- the existing production `prune_hist(...)` treatment is neither moved,
  changed, rerun diagnostically, nor relabeled as proton subtraction;
- the post-prune production state must match the exact accepted pion input;
- the pion stage consumes the exact accepted
  before/template/after production objects;
- baseline pion ownership remains `b_j^0 = s_j * w0_j`;
- final `Y_0` and its existing total production error are observed, not
  recomputed;
- the E.8.2 private statistical error is taken from the existing
  `integral_with_stat_error(final_hist)` result;
- existing public `stage_window_yields`, public `Y_0`, and public yield error
  remain unchanged;
- valid-child stage inventory fails closed to historical Fit 1/Fit 2,
  empirical-residual, Method-A, and Method-B numerical-stage injection;
- active profile remains literally `no_empirical_residual`;
- Method A is numerically absent from E.8.2;
- Method B is numerically absent from E.8.2;
- accepted F.6.2 artifacts/fingerprints remain unchanged.

## Reviewed procedure finalization

The reviewed finalizer:

- retains detached Step-3 presentation state;
- waits until Step-6 data yields exist;
- rerenders the complete procedure PDF to temporary paths without rerunning
  scientific analysis;
- treats renderer failures as finalization failures;
- completes the temporary PDF and page manifest before installation;
- creates recovery copies of the preliminary PDF/manifest pair;
- installs the two final names sequentially;
- restores the preliminary pair on a handled partial-installation failure;
- reports recovery failure truthfully;
- releases transient render state only after both final artifacts are
  installed successfully.

This is a pair-safe recovery transaction, not a literal multi-file filesystem
atomic replacement.

## Reviewed page and presentation contract

The E.8.2 presentation:

- remains after existing E.8/E.8.1 content and before the terminal handoff;
- keeps the handoff last;
- exposes prompt/random;
- exposes dummy subtraction;
- exposes slow-proton removal with distinct pre-prune and post-prune production
  states;
- exposes exact baseline pion input, `B_pi^0`, and `K_0`;
- exposes final canonical `MM_0(t,phi)`, Lambda window, `Y_0`, statistical
  uncertainty, and existing total uncertainty;
- exposes producer-owned diagnostic stage-window integrals;
- never makes dormant Fit 1/Fit 2, Method A, or Method B a baseline production
  stage.

---

# Test/check boundary

Codex-reported deterministic checks are accepted as source-review evidence but
remain:

`NOT RUN by ChatGPT`

The final reviewed phase record reports:

- focused E.8.2 suite — PASS, 23 tests after Fix.3;
- `testing.test_full_background_subtraction_plots` — PASS, 102 tests with 18
  expected PyROOT/superseded-contract skips;
- `testing.test_pion_t_bin_particle_stage` — PASS, 7 tests;
- `testing.test_pion_hgcer_phase_f_runtime_contract` — PASS, 3 tests;
- `testing.test_proton_cleaning_runtime_status` — unchanged pre-existing
  34-pass/one-error result, reproduced at exact starting HEAD, with no E.8.2
  ownership;
- required `py_compile` checks — PASS;
- memory manifest regeneration/check — PASS;
- memory health/bootstrap checks — PASS;
- `git diff --check` — PASS.

The unchanged proton-cleaning fixture error remains outside E.8.2 scope. Do not
repair, suppress, or reinterpret it in this reconciliation.

These checks do not establish ROOT/PyROOT integration, actual procedure-PDF
rendering on the Jefferson Lab farm, full `main.py` execution, or production
runtime acceptance.

---

# Status reconciliation

Update durable state to:

- E.8 — `ACTIVE`
- E.8.1 narrow accepted runtime closures — unchanged
- E.8.2 baseline full-analysis stage audit, including Fix.1/Fix.2/Fix.3 —
  `SOURCE REVIEWED`
- E.8.3 detached Method-A reweighting audit — `NEXT`
- F.6.3 — `BLOCKED` pending source-reviewed E.8.3; E.8.2 prerequisite is now
  satisfied
- E.8.4 — `BLOCKED` pending F.6.3
- final E.8 — `BLOCKED` pending E.8.4 and the later runtime/visual gate
- F.6.4 — `BLOCKED` pending completed production-impact evidence
- lifecycle-hook dispatch — `BLOCKED / DEFERRED`

Do not mark E.8.2 `CLOSED / RUNTIME VALIDATED`.

Do not create an E.8.2 farm-evidence record.

The approved roadmap explicitly permits E.8.2 -> E.8.3 development through
source review without a farm run at this point.

---

# Allowed files

Modify only:

- `docs/memory/CURRENT.md`
- `docs/memory/phases/e8-2-baseline-stage-audit.md`
- `docs/memory/manifest.json`

This reconciliation contract is also an intended tracked file:

- `docs/memory/phases/e8-2-baseline-stage-audit-source-review-reconciliation-task-contract.md`

No other file may change.

---

# Frozen reviewed files

The following must remain byte-for-byte unchanged during reconciliation:

- `src/binning/calculate_yield.py`
- `src/cuts/full_background_subtraction_plots.py`
- `src/cuts/rand_sub.py`
- `src/main.py`
- `testing/test_e8_2_baseline_stage_audit.py`
- `docs/memory/phases/e8-2-baseline-stage-audit-task-contract.md`
- `docs/memory/phases/e8-2-baseline-stage-audit-fix1-task-contract.md`
- `docs/memory/phases/e8-2-baseline-stage-audit-fix2-task-contract.md`
- `docs/memory/phases/e8-2-baseline-stage-audit-fix3-task-contract.md`
- proton-cleaning source/tests;
- pion scientific source outside the already reviewed E.8.2 diff;
- `src/utility/utility.py`
- `src/utility/background_config.py`
- E.8.1 profiles/collectors/wrappers;
- accepted F.6.2 artifacts/evidence;
- all farm scripts and launchers.

If any reviewed source/test file has changed since the PASS bundle, hard stop.

---

# Required CURRENT.md reconciliation

Update `docs/memory/CURRENT.md` so it records one unambiguous current state:

1. E.8 remains `ACTIVE` and presentation-only.
2. E.8.2 is now `SOURCE REVIEWED` from independent ChatGPT inspection of the
   complete cumulative actual diff.
3. The reviewed E.8.2 chain includes the true same-traversal pre-proton
   diagnostic, actual pre-/post-prune production snapshots, exact baseline pion
   before/template/after objects, and existing `Y_0`.
4. Production physics, `w0`, proton factors, pruning, random/dummy
   normalization, canonical binning, public yields/errors, Method A, Method B,
   and F.6.2 science remain unchanged.
5. Codex-reported local checks were `NOT RUN by ChatGPT`.
6. No ROOT/PyROOT/full-analysis/farm/runtime acceptance is claimed.
7. E.8.3 is now unblocked and is `NEXT`.
8. F.6.3 remains `BLOCKED` pending source-reviewed E.8.3.
9. No farm action is required before beginning E.8.3.

Set exact NEXT to:

```text
user-controlled commit/push of the reviewed E.8.2 source/test/memory set
-> ChatGPT pushed-state review
-> E.8.3 detached Method-A reweighting audit
```

For E.8.3, preserve the roadmap boundary:

- detached/non-production only;
- consume accepted F.4/F.5/F.6.1/F.6.2 results;
- show `b_j^0 = s_j*w0_j` to `b_j^A = s_j*w0_j*C_j`;
- no child `(t,phi)` renormalization;
- no Method-B numerical input;
- no empirical residual Fit 1/Fit 2;
- do not construct the parallel full production branch; F.6.3 alone owns that.

---

# Required E.8.2 phase-record reconciliation

Update:

`docs/memory/phases/e8-2-baseline-stage-audit.md`

so that:

- status becomes `SOURCE REVIEWED`;
- it records independent ChatGPT review of the complete cumulative
  E.8.2/Fix.1/Fix.2/Fix.3 actual diff as PASS;
- it preserves the exact accepted baseline path and pre-/post-prune ownership;
- it preserves the pair-safe finalization semantics;
- it preserves the public-output non-regression result;
- Codex-reported deterministic checks are explicitly `NOT RUN by ChatGPT`;
- the known proton fixture error remains documented as pre-existing and
  out-of-scope;
- no ROOT/PyROOT/farm/runtime validation is claimed;
- no farm run is requested;
- NEXT is user commit/push -> pushed-state review -> E.8.3.

Do not rewrite historical chronology or accepted runtime evidence.

---

# Validation

Regenerate:

`docs/memory/manifest.json`

Run the established repository-memory checks:

- manifest regeneration/check;
- memory health/integrity;
- memory bootstrap check;
- `git diff --check`.

Also verify:

```bash
git status --short
git diff --stat
git diff --name-only
```

Explicitly confirm that reconciliation introduced no new diff in any reviewed
source or test file.

No scientific or regression suite needs to be modified. This is memory-only.

---

# Required complete review bundle

Refresh the root temporary:

`kaonlt_review.diff`

It must contain the **complete cumulative E.8.2 worktree**, now including this
source-review reconciliation, from committed HEAD `5b2ba7f...`.

Include:

1. complete tracked `git diff`;
2. complete `git diff --no-index -- /dev/null ...` additions for every intended
   untracked file.

At minimum include complete additions for:

- `docs/memory/phases/e8-2-baseline-stage-audit-task-contract.md`
- `docs/memory/phases/e8-2-baseline-stage-audit-fix1-task-contract.md`
- `docs/memory/phases/e8-2-baseline-stage-audit-fix2-task-contract.md`
- `docs/memory/phases/e8-2-baseline-stage-audit-fix3-task-contract.md`
- `docs/memory/phases/e8-2-baseline-stage-audit-source-review-reconciliation-task-contract.md`
- `docs/memory/phases/e8-2-baseline-stage-audit.md`
- `testing/test_e8_2_baseline_stage_audit.py`

Do not add `kaonlt_review.diff` to Git.

---

# Acceptance criteria

The reconciliation is ready for final ChatGPT review only if:

1. branch remains `test`;
2. committed HEAD remains exactly
   `5b2ba7f92dc51b303257dbf209df687415364602`;
3. the reviewed E.8.2 source/test candidate is byte-for-byte unchanged;
4. only the three allowlisted existing memory files plus this new reconciliation
   contract change during reconciliation;
5. E.8.2 status is `SOURCE REVIEWED`;
6. no E.8.2 runtime/farm acceptance is claimed;
7. E.8.3 is `NEXT`;
8. F.6.3 remains `BLOCKED` pending source-reviewed E.8.3;
9. Codex-reported checks are explicitly distinguished from ChatGPT execution;
10. the pre-existing proton regression is not repaired or suppressed;
11. Method A remains detached/non-production;
12. Method B remains diagnostic only and numerically absent from E.8.2;
13. active background profile remains `no_empirical_residual`;
14. Fit 1/Fit 2 remain dormant;
15. manifest/memory checks and `git diff --check` pass;
16. a complete cumulative `kaonlt_review.diff` is regenerated;
17. no commit, push, or farm run occurs.

---

# Hard stop

Stop and report rather than broadening scope if:

- any reviewed source/test file differs from the PASS bundle;
- reconciliation would require a source/test edit;
- any physics status or runtime claim is ambiguous;
- E.8.3 would require production Method-A application;
- Method B would become numerical;
- empirical Fit 1/Fit 2 would be reintroduced;
- the known proton regression changes;
- unrelated worktree changes are discovered.

Codex must not commit, push, update remote refs, or run the Jefferson Lab farm.
