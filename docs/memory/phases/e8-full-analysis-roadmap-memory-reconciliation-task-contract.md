# KaonLT E.8 Full-Analysis Roadmap — Memory Reconciliation Task Contract

## Objective

Reconcile tracked repository memory with the newly agreed E.8 program before any further implementation work.

This is a **documentation/memory-only** task. It must record:

1. the fresh accepted E.8.1.Fix.5/Fix.6 `Q4p4W2p74 / Left / lowe` runtime evidence already reviewed by ChatGPT;
2. the user's decision to **defer the remaining E.8.1 canonical-five farm expansion for now**;
3. the expanded role of E.8 as the complete visual audit of the kaon missing-mass/yield analysis chain;
4. the explicit Method-A reweighting presentation:
   - baseline pion background using `w0`;
   - reweighted pion background using `w0*C`;
   - baseline vs reweighted comparison, signed difference, and ratio where meaningful;
   - canonical `(t,phi)` redistribution and parent closure;
5. the new dependency order:
   - E.8.2 baseline full-analysis audit;
   - E.8.3 detached Method-A reweighting audit;
   - F.6.3 parallel full procedure with Method A;
   - E.8.4 baseline-vs-Method-A production-impact audit;
   - final E.8 closure;
   - F.6.4 explicit production-promotion decision;
6. the agreed validation cadence: local/source-reviewed development may proceed through coherent implementation steps without a Jefferson Lab farm run after every update; farm validation remains mandatory at defined runtime milestones.

No analysis source, tests, profiles, renderers, collectors, launchers, or production behavior may change in this task.

---

## Exact starting source identity

Authoritative remote `test` HEAD at contract creation:

`0ec29d4e1bb345eb37e8cca35b8b7e5cbe1b4d5`

Commit message:

`Re-pin E8.1 profile to Fix5 source`

Parent:

`53fd262b730af8f1254e411a38231aebeb6a1da3`

Codex must first verify:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
```

Required branch: `test`.

Required HEAD: exactly `0ec29d4e1bb345eb37e8cca35b8b7e5cbe1b4d5`.

The only permitted pre-existing worktree change is this newly placed task-contract file itself, if it is untracked at task start. Any other tracked or untracked change is a hard stop unless the user explicitly authorizes it.

Do not pull, reset, stash, clean, switch branches, commit, or push.

---

## Required startup reading

Before editing, read in this order:

1. `docs/memory/AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read only these task-relevant records:

- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/phase-f6-method-a-production-promotion.md`
- `docs/memory/phases/PHASE_HISTORY.md`
- `docs/memory/decisions/e8-f6-2-figure-library-implementation-contract.md`
- `docs/memory/phases/e8-1-fix5-overlay-pdf-geometry.md`
- `docs/memory/phases/e8-1-fix6-bundle-profile-repin.md`
- `docs/memory/evidence/e8-1-fix4-left-lowe-overlay-blocker.md`
- `docs/memory/evidence/f6-1-runtime-closure.md`
- `docs/memory/evidence/f6-2-scientific-runtime-closure.md`
- `docs/memory/evidence/f6-2-fix5-presentation-runtime-closure.md`
- `docs/memory/manifest.json`
- repository memory manifest/integrity helper documentation needed to regenerate/check `manifest.json`.

Do not broadly re-read or redesign unrelated phases.

---

## Evidence precedence and accepted runtime facts to record

Use current source/diff first, then the following fresh farm evidence, which ChatGPT independently inspected from:

`KaonLT_E8_1_Fix6_Q4p4W2p74_Left_lowe_20260923-090850.zip`

### Bundle/source provenance

- bundle/profile commit:
  `0ec29d4e1bb345eb37e8cca35b8b7e5cbe1b4d5b`
- required analysis/procedure source:
  `53fd262b730af8f1254e411a38231aebeb6a1da3`
- required source ancestor check: passed
- unexpected committed files after required source: none
- source worktree: clean in bundle provenance
- bundle: complete
- bundle errors: none

### Frozen F.6.2 authority preserved

- accepted F.6.2 JSON SHA-256:
  `5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1`
- accepted F.6.2 artifact fingerprint:
  `ee713b70de898bad8fa61164cbc8a54886af1ea1eb1e712ed95de17df8890cd0`
- accepted F.6.2 validation fingerprint:
  `7edc73fce20ad7dc7622b8c23a7ba7e8986e367ace605884e010945595370b3b`

### Fresh procedure PDF

- PDF SHA-256:
  `d809a80fb1a46602dbb6cb930ed1d6c8932c2583ce8c5158ac82a78e6e6301e5`
- page-manifest SHA-256:
  `d4bc89f49ce374f63055afd974d734bf782cea6e55cefcc26caf3c2f83028dc1`
- page count: 47
- `renderer_failures=[]`

### Farm tests/checks

Fresh farm evidence ran the real ROOT/PyROOT conditional regression:

- `testing.test_full_background_subtraction_plots`: 102 tests, 14 skips, PASS
- real E.8 overlay PDF regression produced `e8-overlay-pdf-regression.pdf`
- Method-A acceptance contract tests: 7 PASS
- Phase-F runtime contract tests: 3 PASS
- generic collector tests: 28 PASS
- `git diff --check`: PASS

### Independent PDF review

ChatGPT independently inspected pages 37–47.

- page 37 context: PASS
- pages 38 / 41 / 44 persisted overlays: PASS
- pages 39 / 40 / 42 / 43 / 45 / 46 maps: PASS
- page 47 handoff: PASS
- pages 38 / 41 / 44 visibly contain:
  - E.8 header;
  - all L/B/A legend entries;
  - complete first plot row;
  - no top clipping;
  - first canonical child title.
- independent text extraction confirmed:
  - `L: upstream 0 < NPE <= 2 diagnostic reference`
  - `B: physical pion control, NPE > 2, baseline w0`
  - `A: same B population, w0*C`
  - first child `Left-lowe phi0 [-180, -140)`

### Status warranted by this evidence

Record only the narrow statuses supported:

- **E.8.1.Fix.5 — CLOSED / RUNTIME VALIDATED**
  for the persisted-overlay geometry repair.
- **E.8.1.Fix.6 — CLOSED / RUNTIME VALIDATED**
  for the profile/provenance gate.
- **E.8.1 canonical-five expansion — DEFERRED**
  at the user's explicit request. Do not claim canonical-five E.8.1 closure.
- **E.8 overall — ACTIVE**.

Do not convert a single `Left / lowe` gate into canonical-five acceptance.

---

## New E.8 governing presentation requirement

The prior narrow E.8 figure-library contract remains historical evidence for the completed F.6.2 presentation work, but it no longer owns the forward E.8 roadmap.

Create a new durable architectural/roadmap decision that establishes:

> E.8 is the complete visual audit of how the kaon missing-mass spectrum and signal-region yield evolve through the analysis. Every substantive stage must expose the authoritative input spectrum, the component/treatment applied, and the authoritative output spectrum. The presentation culminates in the final canonical `(t,phi)` missing-mass spectra and extracted yields.

Presentation remains non-authoritative and must only consume authoritative upstream runtime objects, persisted snapshots, accepted detached Method-A artifacts, or later F.6.3 branch outputs. It must never recompute a fit, factor, correction, normalization, or yield merely for plotting.

The baseline production branch remains unchanged.

Method B remains diagnostic/cross-check only and is numerically excluded from every Method-A correction/application.

---

## Complete approved roadmap to record

### E.8.2 — Baseline full-analysis stage audit

Status after this memory reconciliation:

**NEXT**

Purpose: expose the complete existing baseline analysis as a coherent before/component/after spectrum and yield chain without changing physics.

#### E.8.2a — Prompt/random subtraction

Show:

- prompt kaon MM input;
- random contribution removed;
- after-random MM output;
- difference/overlay where useful.

#### E.8.2b — Dummy subtraction

Show:

- after-random input;
- normalized dummy contribution;
- after-dummy output;
- difference/overlay where useful.

#### E.8.2c — Slow-proton subtraction

Show:

- after-random/dummy kaon input;
- proton contamination estimate;
- proton-cleaned kaon output;
- before/after comparison;
- existing PID/weight diagnostics as supporting pages.

Important known implementation gap:

The current yield-level `_process_yield_data_tree()` applies accepted proton cleaning during initial `(t,phi)` filling, so existing per-cell raw/random/dummy snapshots must not be relabeled as pre-proton. E.8.2 implementation must capture or consume a true authoritative pre-proton stage snapshot rather than reconstructing or guessing one.

#### E.8.2d — Baseline pion subtraction

Show:

- proton-cleaned kaon input;
- baseline pion background actually removed;
- baseline pion-subtracted kaon output;
- before/after comparison;
- difference.

Baseline pion contribution remains:

`b_j = s_j * w0_j`.

#### E.8.2e — Residual background Fit 1

Show:

- post-pion input MM;
- actual production Fit-1 background component removed;
- post-Fit-1 MM output;
- before/after comparison;
- existing fit/window information as appropriate.

Do not refit for presentation.

#### E.8.2f — Residual background Fit 2

Show:

- Fit-1 output entering Fit 2;
- actual production Fit-2 background component removed;
- final background-subtracted MM;
- before/after comparison.

Do not refit for presentation.

#### E.8.2g — Final baseline canonical `(t,phi)` MM spectra

For every canonical `t`, show all nine canonical phi children, including explicit empty/invalid status.

Each populated final panel must expose:

- canonical `t` identity;
- phi index/range;
- final MM spectrum actually used for yield extraction;
- Lambda signal/integration window;
- relevant final fitted/background information;
- extracted yield;
- statistical uncertainty.

#### E.8.2h — Baseline stage-yield audit

Expose the existing stage-yield progression per `(t,phi)`:

- prompt;
- after random;
- after dummy;
- after proton;
- after pion;
- after Fit 1;
- after Fit 2/final.

Summarize final baseline:

`Y0(t,phi)`

versus phi for each canonical t.

Do not recalculate a yield in the renderer when the authoritative producer can persist/provide it.

---

### E.8.3 — Detached Method-A reweighting audit

Status after this memory reconciliation:

**BLOCKED** pending E.8.2 source-reviewed completion.

This remains detached/non-production and consumes accepted F.4/F.5/F.6.1/F.6.2 results only.

#### E.8.3a — Reweighting operation

Explicitly show:

`b_j^0 = s_j * w0_j`

to

`b_j^A = s_j * w0_j * C_j`.

Method A changes only the accepted pion event contribution by the accepted F.4 parent-preserving factor.

No child-by-child renormalization.

#### E.8.3b — Baseline pion background before/after reweighting

For every canonical t parent show:

- baseline pion background `B_pi^0(MM)`;
- Method-A-reweighted pion background `B_pi^A(MM)`;
- overlay;
- signed difference:
  `Delta B_pi = B_pi^A - B_pi^0`;
- ratio:
  `R_pi = B_pi^A / B_pi^0`
  only where numerically well-defined and presentation-safe.

This baseline-vs-reweighted pion-background comparison is a first-class required E.8 result, not a legend-level implication.

#### E.8.3c — Canonical `(t,phi)` redistribution and parent closure

For every canonical t show all nine children:

- `B_pi,tphi^0`;
- `B_pi,tphi^A`;
- `Delta B_pi,tphi`.

Explicitly display parent preservation:

`sum_phi B_pi,tphi^A ~= sum_phi B_pi,tphi^0`

with stored/accepted closure residual and relative residual.

Never independently normalize a child.

#### E.8.3d — Existing F.6.2 acceptance/HGCer explanation

Retain the accepted useful diagnostics:

- L/B/A normalized shapes;
- `delta x xptar`;
- `delta x yptar`;
- MM x acceptance maps;
- support/OOD;
- effective statistics;
- kaon-window refinement metrics.

Their role becomes explanatory: they explain the observed `B_pi^0 -> B_pi^A` redistribution. They do not replace the direct before/after reweighting plots.

---

### F.6.3 — Parallel full procedure plus Method A

Status after memory reconciliation:

**BLOCKED** pending E.8.2 and E.8.3 source-reviewed prerequisites, **not** pending final E.8 closure.

Preserve the current full yield calculation unchanged as the baseline branch.

The parallel Method-A branch changes exactly:

`w0_j -> w0_j * C_j`

while filling the pion-subtraction template.

Everything else remains frozen and identical between branches:

- random/dummy subtraction;
- slow-proton treatment;
- pion components/fits/windows/amplitudes except the event-level Method-A factor;
- canonical binning;
- SIMC;
- Fit-1 and Fit-2 algorithms/configuration;
- yield extraction;
- efficiencies;
- acceptance;
- L/T separation;
- cross-section formulas.

Downstream residual background fits rerun normally on the altered Method-A post-pion spectrum with algorithms/configuration frozen.

Baseline branch must remain reproducible and unchanged when Method A is disabled.

---

### E.8.4 — Baseline-vs-Method-A production-impact audit

Status after memory reconciliation:

**BLOCKED** pending F.6.3.

This presentation consumes the authoritative two-branch outputs produced by F.6.3. It must not construct those outputs itself.

#### E.8.4a — Actual pion-subtraction consequence

Using the same proton-cleaned input compare:

- `B_pi^0`;
- `B_pi^A`;
- `B_pi^A - B_pi^0`;
- baseline post-pion kaon `K_pi^0`;
- Method-A post-pion kaon `K_pi^A`;
- `K_pi^A - K_pi^0`.

#### E.8.4b — Fit-1 consequence

Compare baseline and Method-A branches at:

- Fit-1 input;
- fitted background component;
- Fit-1 output;
- branch difference.

#### E.8.4c — Fit-2 consequence

Compare baseline and Method-A branches at:

- Fit-2 input;
- fitted background component;
- final Fit-2 output;
- branch difference.

#### E.8.4d — Final canonical `(t,phi)` MM comparison

For every populated canonical cell show:

- baseline final MM;
- Method-A final MM;
- overlay;
- signed difference;
- identical signal/integration window.

#### E.8.4e — Final yield comparison

For every canonical `(t,phi)`:

- `Y0(t,phi)`;
- `YA(t,phi)`;
- `DeltaY = YA - Y0`;
- `DeltaY/Y0` where defined.

Summarize:

- versus phi for each t;
- versus t for each setting;
- later across all five canonical Q4p4W2p74 settings.

The Method-A minus baseline shift is a correction effect, not automatically a systematic uncertainty.

---

### Final E.8 closure

Status:

**BLOCKED** pending E.8.4 and later required runtime/visual validation.

Final E.8 closure requires the full visual chain:

1. prompt;
2. random subtraction;
3. dummy subtraction;
4. slow-proton subtraction;
5. baseline pion subtraction;
6. Fit 1;
7. Fit 2;
8. final `MM(t,phi)`;
9. final baseline `Y(t,phi)`;
10. stage-by-stage yield evolution;
11. detached `w0 -> w0*C` presentation;
12. baseline pion background vs reweighted pion background;
13. signed difference/ratio;
14. `(t,phi)` redistribution;
15. parent closure;
16. F.6.2 acceptance/HGCer explanation;
17. actual F.6.3 baseline-vs-Method-A post-pion spectrum;
18. baseline-vs-Method-A Fit 1;
19. baseline-vs-Method-A Fit 2;
20. final baseline-vs-Method-A `MM(t,phi)`;
21. final baseline-vs-Method-A `Y(t,phi)`;
22. `DeltaY` and fractional `DeltaY`.

Do not claim final E.8 runtime closure from the already accepted Left/lowe Fix.6 overlay gate.

---

### F.6.4 — Explicit production-promotion decision

Status remains:

**BLOCKED**

It follows the completed F.6.3/E.8.4 evidence and is the only phase that may decide whether Method A becomes production pion-background treatment.

No automatic promotion follows from detached validation or presentation.

---

## Validation cadence to record

The user explicitly does not want a Jefferson Lab farm run after every development update.

Record the workflow as:

- memory reconciliation: local/documentation only;
- E.8.2 implementation: deterministic local/source review first;
- E.8.3 implementation: deterministic local/source review first;
- F.6.3/E.8.4 implementation: deterministic local/source review through a coherent end-to-end branch;
- farm remains mandatory for ROOT/PyROOT/full-runtime claims, but is used as a **milestone gate**, not an edit-by-edit loop;
- after the coherent F.6.3 + E.8.4 branch is source-reviewed, use one narrow `Q4p4W2p74 / Left / lowe` end-to-end farm gate;
- inspect fresh artifacts, then PASS or perform one coherent repair;
- broaden only after the narrow gate passes.

Do not write that local tests substitute for farm validation.

---

## Required memory edits

### Existing files to update

1. `docs/memory/CURRENT.md`
   - reconcile live source identity with pushed Fix.6 commit;
   - record fresh Left/lowe Fix.6 runtime evidence and narrow closure;
   - mark remaining E.8.1 canonical-five expansion DEFERRED;
   - replace stale NEXT with E.8.2 baseline full-analysis audit as NEXT;
   - update blockers/dependencies exactly as above.

2. `docs/memory/MEMORY.md`
   - add the durable E.8 principle that presentation is a complete visual audit of authoritative upstream stage objects/snapshots;
   - preserve presentation-only/non-recalculation ownership;
   - preserve Method-A/Method-B boundaries and parent-preserving rule;
   - record milestone farm cadence only if it belongs here; otherwise keep it in USER/roadmap ownership and avoid duplication.

3. `docs/memory/roadmap/STATUS.md`
   - replace the stale E.8/F.6.3 dependency with the approved E.8.2 -> E.8.3 -> F.6.3 -> E.8.4 -> final E.8 -> F.6.4 sequence;
   - preserve closed F.1–F.6.2 statuses.

4. `docs/memory/phases/phase-f6-method-a-production-promotion.md`
   - reconcile the detailed Phase-F roadmap with the new dependency order;
   - retain every frozen scientific boundary;
   - expand E.8.2/E.8.3/E.8.4 content enough to preserve the complete agreed plot inventory;
   - make clear E.8.3 direct reweighting plots precede F.6.3;
   - make clear E.8.4 consumes actual F.6.3 outputs.

5. `docs/memory/phases/PHASE_HISTORY.md`
   - concise durable summary only;
   - do not duplicate the full roadmap.

6. `docs/memory/phases/e8-1-fix5-overlay-pdf-geometry.md`
   - update status to `CLOSED / RUNTIME VALIDATED`;
   - cite the new evidence record;
   - preserve source identity and repair scope.

7. `docs/memory/phases/e8-1-fix6-bundle-profile-repin.md`
   - update status to `CLOSED / RUNTIME VALIDATED`;
   - record pushed commit `0ec29d4e1bb345eb37e8cca35b8b7e5cbe1b4d5b`;
   - cite the new evidence record;
   - remove stale “no Fix.6 commit exists yet” wording.

8. `docs/memory/decisions/e8-f6-2-figure-library-implementation-contract.md`
   - do not rewrite historical details;
   - add a concise top-level supersession/status note stating that this remains the historical narrow frozen-F.6.2 figure-library contract, while future E.8 roadmap ownership has moved to the new full-analysis procedure roadmap decision.
   - do not delete it.

9. `docs/memory/manifest.json`
   - regenerate after all versioned memory edits.

### New durable records to add

10. `docs/memory/evidence/e8-1-fix6-left-lowe-runtime-closure.md`
    - record the exact fresh bundle/source/artifact/PDF/test/visual evidence above;
    - scope closure narrowly to Fix.5 geometry and Fix.6 provenance;
    - explicitly state E.8.1 canonical-five is not closed.

11. `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
    - own the complete forward architectural/presentation roadmap;
    - include E.8.2a–h, E.8.3a–d, F.6.3 handoff, E.8.4a–e, final E.8 closure, and farm milestone cadence;
    - explicitly preserve presentation-only ownership and frozen production/scientific boundaries.

12. This task contract:
    `docs/memory/phases/e8-full-analysis-roadmap-memory-reconciliation-task-contract.md`

Do not add unrelated memory records.

---

## Frozen files / forbidden edits

This task must not modify anything outside `docs/memory/`.

In particular, do not modify:

- `src/`
- `testing/`
- `run_Prod_Analysis.sh`
- `main.py`
- any configuration
- any JSON analysis artifact
- any PDF
- any bundle collector/profile
- any physics code
- any local-only root `AGENTS.md`
- `.codex/`

No production correction, normalization, cut, fit, template, prior, binning, or analysis-state behavior may change.

---

## Required wording/ownership constraints

Preserve these scientific facts exactly:

- baseline production remains authoritative until explicit F.6.4 promotion;
- Method B remains diagnostic/cross-check only and never adjusts pion weights;
- Method A uses exactly the accepted F.4 parent-preserving correction;
- no child `(t,phi)` independent renormalization;
- E.8 presentation must not recompute weights, fits, factors, yields, or corrections;
- F.6.3 is the only step here that constructs the actual parallel Method-A full-analysis branch;
- Method-A minus baseline is a correction effect, not automatically an uncertainty;
- source review does not equal farm validation;
- the accepted F.6.2 scientific JSON and fingerprints remain frozen.

Do not change the status of F.1 through F.6.2.

Do not describe the deferred four E.8.1 settings as failed.

---

## Deterministic local validation

Run the repository's established memory checks, including the actual local commands used by this repository for:

- memory manifest regeneration;
- memory manifest check;
- memory integrity/health check;
- memory bootstrap check;
- `git diff --check`.

Also run an exact path audit proving every changed/new file is under `docs/memory/`.

No ROOT/PyROOT, analysis runtime, farm, bundle, or PDF execution is required or authorized.

If a memory helper fails because the roadmap introduces a legitimate new cross-reference, repair only the memory/reference issue inside the allowed scope. Do not broaden into source changes.

---

## Diff audit requirements

Before finishing, report:

```bash
git status --short
git diff --stat
git diff -- docs/memory/
git diff --check
```

If the new files are still untracked, include their full contents in the review material.

If the complete diff is too large for terminal review, write a root-level temporary review bundle:

`kaonlt_review.diff`

containing:

1. the complete tracked `git diff`;
2. for every intended untracked file:
   `git diff --no-index -- /dev/null <file> >> kaonlt_review.diff || true`

Do not add the temporary review bundle to Git.

---

## Acceptance criteria

The memory reconciliation passes only if:

1. live starting HEAD is exactly the required commit;
2. only `docs/memory/` changes;
3. the fresh Fix.5/Fix.6 Left/lowe runtime evidence is recorded accurately;
4. E.8.1 remaining canonical-five expansion is explicitly DEFERRED, not closed or failed;
5. the complete agreed plot roadmap is preserved, including:
   - random;
   - dummy;
   - proton;
   - baseline pion;
   - Fit 1;
   - Fit 2;
   - final `(t,phi)` MM;
   - final/stage yields;
   - explicit `w0 -> w0*C`;
   - baseline pion background vs reweighted pion background;
   - signed difference/ratio;
   - `(t,phi)` redistribution and parent closure;
   - F.6.2 supporting acceptance/HGCer diagnostics;
   - F.6.3 actual parallel branch;
   - post-pion/Fit1/Fit2/final-MM/final-yield baseline-vs-Method-A comparisons;
6. F.6.3 is no longer blocked on **final** E.8 closure, but remains blocked on E.8.2/E.8.3 prerequisites;
7. E.8.4 remains blocked on F.6.3;
8. F.6.4 remains blocked pending full production-impact evidence;
9. memory manifest is regenerated and all established memory checks pass;
10. no scientific/runtime/source file changes.

---

## Hard stop

Stop immediately and report without editing if:

- branch is not `test`;
- HEAD is not exactly `0ec29d4e1bb345eb37e8cca35b8b7e5cbe1b4d5b`;
- worktree contains unexpected pre-existing changes;
- a requested memory update would require modifying analysis/test/profile/runtime source;
- the repository memory tooling reports a structural problem that cannot be repaired within this narrow memory-only scope.

Codex must not commit, push, run the Jefferson Lab farm, or claim runtime validation beyond the supplied and already reviewed evidence.
