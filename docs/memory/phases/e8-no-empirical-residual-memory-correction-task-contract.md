# KaonLT E.8 Current-Analysis-Chain Memory Correction Task Contract

## Objective

Correct the newly reconciled E.8 roadmap memory before any E.8.2 source implementation.

The pushed roadmap at:

`942444fb9586dc53e847039c97c6c3b52fe9a351`

incorrectly incorporated the historical empirical residual-background **Fit 1 / Fit 2**
machinery into the active E.8 baseline and Method-A chains.

That is wrong.

The current accepted analysis uses the active background profile:

`BG_OPT_ACTIVE_PROFILE = "no_empirical_residual"`

with both empirical residual scales forced to zero. The old Fit-1/Fit-2 code remains
in source as dormant historical machinery but is **not part of the accepted current
KaonLT scientific procedure and must not appear in the E.8 roadmap, E.8 procedure
presentation, F.6.3 Method-A branch, or E.8.4 production-impact comparison.**

The corrected current scientific/presentation chain is:

```text
selected kaon
    -> random subtraction
    -> dummy subtraction
    -> slow-proton cleaning
    -> baseline pion background treatment (w0)
    -> Method-A pion-background reweighting (w0 -> w0*C)
    -> cleanest kaon MM(t,phi)
    -> extracted Y(t,phi)
```

For branch comparison, baseline and Method A share the same proton-cleaned input:

```text
baseline:
K_proton - B_pi^0(w0) = K_0 -> MM_0(t,phi) -> Y_0(t,phi)

Method A:
K_proton - B_pi^A(w0*C) = K_A -> MM_A(t,phi) -> Y_A(t,phi)
```

There is **no empirical residual Fit 1 or Fit 2 after pion subtraction** in this
accepted chain.

This task is memory/documentation-only. Do not modify source or runtime behavior.

---

## Exact starting state

Required branch:

`test`

Required HEAD:

`942444fb9586dc53e847039c97c6c3b52fe9a351`

Commit message:

`Reconcile E8 full-analysis roadmap`

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
```

The worktree must be clean except for this newly placed untracked contract.

Hard stop on any different branch, HEAD, or unrelated worktree change.

Do not pull, reset, stash, clean, switch branches, commit, push, or run the farm.

---

## Required startup reading

Read in the repository-mandated order:

1. `docs/memory/AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read only the task-relevant records:

- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/phases/phase-f6-method-a-production-promotion.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/PHASE_HISTORY.md`
- `docs/memory/phases/e8-full-analysis-roadmap-memory-reconciliation-task-contract.md`
- `docs/memory/phases/e8-full-analysis-roadmap-memory-reconciliation-fix1-task-contract.md`
- `docs/memory/manifest.json`

Inspect the current source only to verify the specific runtime gate facts below:

- `src/utility/background_config.py`
- the Fit-1/Fit-2 conditional blocks in `src/binning/calculate_yield.py`

Do not edit source.

---

## Source facts already independently verified

At HEAD `942444fb9586dc53e847039c97c6c3b52fe9a351`,
`src/utility/background_config.py` contains:

```python
BG_STAT_SCALE1 = 0.0
BG_STAT_SCALE2 = 0.0

BG_OPT_ACTIVE_PROFILE = "no_empirical_residual"
```

and the active profile contains:

```python
"no_empirical_residual": {
    "force_bg_stat_scale1": 0.0,
    "force_bg_stat_scale2": 0.0,
}
```

The historical Fit-1/Fit-2 subtraction code in
`src/binning/calculate_yield.py` is conditional on the corresponding resolved
background scale being greater than zero.

Therefore:

- the old empirical residual-fit code still exists in source;
- its mere existence does **not** make it part of the active scientific procedure;
- the active accepted profile disables it;
- E.8 must not resurrect it as a presentation or Method-A stage.

This is a source-level statement. Do not claim a new farm validation from this memory task.

---

## User scientific clarification to preserve

The accepted current chain is:

```text
random/dummy
    -> slow proton
    -> pion background
    -> Method-A reweighting of the pion background
    -> cleanest kaon sample
    -> canonical MM(t,phi)
    -> yields
```

The historical empirical Fit 1 / Fit 2 procedure predates the current pion-model
background treatment and must not be reintroduced.

Method A is **reweighting, not rebinning**:

```text
b_j^0 = s_j * w0_j
b_j^A = s_j * w0_j * C_j
```

Canonical t/phi binning remains frozen.

Method A redistributes the pion contribution across the fixed accepted phase space while
preserving the accepted parent-t normalization contract. It does not introduce a new
empirical residual fit.

---

# Required memory corrections

## 1. `docs/memory/CURRENT.md`

Correct the `Next Action` and any other active-state text that currently includes
Fit 1 / Fit 2.

E.8.2 must now describe the **baseline current chain only**:

1. prompt/random;
2. dummy;
3. slow proton;
4. baseline pion treatment;
5. final baseline canonical `MM_0(t,phi)`;
6. baseline final/stage yields.

Remove Fit 1 / Fit 2 from E.8.2.

Retain the known true-pre-proton snapshot blocker/requirement.

Preserve:

- E.8 `ACTIVE`;
- E.8.2 `NEXT`;
- E.8.3 blocked on source-reviewed E.8.2;
- F.6.3 blocked on source-reviewed E.8.2/E.8.3;
- E.8.4 blocked on F.6.3;
- F.6.4 blocked pending production-impact evidence;
- all existing narrow Fix.5/Fix.6 runtime statuses.

Add a concise durable clarification that legacy empirical residual Fit 1/Fit 2 are
disabled by the active `no_empirical_residual` profile and are outside the current E.8
scientific chain.

---

## 2. `docs/memory/MEMORY.md`

Add one concise durable scientific-boundary statement under production ordering or
diagnostics:

- legacy empirical residual Fit 1 / Fit 2 are historical dormant machinery;
- the accepted active profile is `no_empirical_residual`;
- they are not part of the current production/presentation chain and must not be
  reintroduced by E.8, Method A, diagnostics, or presentation work.

Do not expand MEMORY into a phase log.

---

## 3. `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`

Correct the roadmap comprehensively.

### E.8.2 — baseline full-analysis audit

The corrected baseline page chain is:

#### E.8.2a — prompt/random
Show:
- prompt input;
- random component removed;
- after-random output.

#### E.8.2b — dummy
Show:
- after-random input;
- normalized dummy component removed;
- after-dummy output.

#### E.8.2c — slow proton
Show:
- after-dummy/pre-proton input;
- accepted proton contamination removed;
- proton-cleaned output;
- existing PID/weight diagnostics as supporting evidence.

Preserve the existing requirement that current yield-level raw/random/dummy snapshots
must not be mislabeled as pre-proton if they were filled after accepted proton weighting.

#### E.8.2d — baseline pion treatment
Show:
- proton-cleaned kaon input;
- actual baseline pion background `B_pi^0` using accepted `w0`;
- baseline clean kaon output:
  `K_0 = K_proton - B_pi^0`;
- before/after and signed difference where useful.

#### E.8.2e — final baseline canonical MM
For each canonical t, show every authoritative phi child:
- final baseline clean kaon MM spectrum;
- t/phi identity;
- signal/integration window;
- final baseline extracted yield and statistical uncertainty;
- explicit empty/invalid status.

There is no Fit-1/Fit-2 stage between baseline pion subtraction and this final spectrum.

#### E.8.2f — baseline stage-yield audit
Expose:
- prompt;
- after random;
- after dummy;
- after proton;
- after baseline pion / final baseline clean sample;
- authoritative final baseline `Y0(t,phi)`.

Do not include Fit-1/Fit-2 stage-yield entries in E.8.

### E.8.3 — detached Method-A reweighting audit

Retain the existing correct content:

- `w0 -> w0*C`;
- `B_pi^0` versus `B_pi^A`;
- signed difference;
- ratio where safe;
- all canonical `(t,phi)` children;
- parent closure;
- accepted F.6.2 HGCer/acceptance diagnostics as explanation.

Make explicit that Method A operates on the accepted pion background and is the next
scientific refinement after the proton-cleaned/baseline-pion picture. It does not invoke
legacy empirical residual fits.

### F.6.3 — parallel Method-A full-analysis branch

Correct the frozen downstream list.

F.6.3:

- keeps baseline branch unchanged;
- uses the same proton-cleaned input;
- changes exactly the accepted pion contribution:
  `w0_j -> w0_j*C_j`;
- produces the Method-A-clean kaon spectrum directly after pion subtraction;
- carries that spectrum through the existing yield extraction with all ordinary current
  non-empirical machinery frozen.

Explicitly state:

- legacy empirical Fit 1/Fit 2 remain disabled;
- F.6.3 must not activate, rerun, tune, or compare them;
- `BG_OPT_ACTIVE_PROFILE = "no_empirical_residual"` remains frozen unless a separate
  future scientific contract explicitly changes it.

### E.8.4 — baseline-versus-Method-A impact

Correct the branch comparison to:

1. same proton-cleaned input;
2. `B_pi^0` versus `B_pi^A`;
3. baseline clean kaon `K_0` versus Method-A clean kaon `K_A`;
4. final canonical `MM_0(t,phi)` versus `MM_A(t,phi)`;
5. `Y0`, `YA`, `DeltaY`, and `DeltaY/Y0`.

Delete the Fit-1 and Fit-2 comparison subsections entirely.

### Final E.8 closure

Correct the closure inventory to require:

- prompt/random;
- dummy;
- slow proton;
- baseline pion;
- final baseline canonical MM/yields;
- explicit `w0 -> w0*C`;
- baseline versus reweighted pion background;
- signed difference/ratio;
- fixed-bin `(t,phi)` redistribution and parent closure;
- F.6.2 explanatory HGCer/acceptance diagnostics;
- actual F.6.3 baseline-versus-Method-A clean kaon MM comparison;
- final `Y0` versus `YA`;
- `DeltaY` and fractional shift.

No Fit 1 / Fit 2.

---

## 4. `docs/memory/phases/phase-f6-method-a-production-promotion.md`

Apply the same correction to the detailed Phase-F roadmap.

Preserve every accepted F.1-F.6.2 scientific/validation detail already restored by the
previous memory repair.

Only correct the forward E.8.2/F.6.3/E.8.4 chain and add the explicit
`no_empirical_residual` boundary.

Do **not** remove historical references to empirical-fit machinery if they are clearly
describing historical source architecture. The correction is that those fits are not
part of the current accepted scientific chain.

Do not rewrite closed F.6.1/F.6.2 records.

---

## 5. `docs/memory/roadmap/STATUS.md`

Correct the concise dependency/status descriptions:

### E.8.2
`NEXT` — baseline chain:
prompt/random -> dummy -> proton -> baseline pion -> final canonical MM/yields.

### E.8.3
unchanged conceptually:
detached `w0 -> w0*C` audit.

### F.6.3
parallel Method-A pion branch only; legacy empirical Fit 1/Fit 2 remain disabled.

### E.8.4
post-pion/final-MM/final-yield baseline-versus-Method-A comparison only.

Remove Fit-1/Fit-2 from the current roadmap summaries.

Statuses/dependencies themselves remain unchanged.

---

## 6. `docs/memory/phases/PHASE_HISTORY.md`

Keep this concise.

If needed, add one sentence to the E.8/F.6 summary:

> Legacy empirical residual Fit 1/Fit 2 remain disabled under the accepted
> `no_empirical_residual` profile and are excluded from the E.8/F.6.3 forward chain.

Do not duplicate the detailed roadmap.

---

## 7. Prior roadmap-reconciliation task contracts

These historical contracts encoded the mistaken Fit-1/Fit-2 roadmap and must not remain
apparently current without qualification:

- `docs/memory/phases/e8-full-analysis-roadmap-memory-reconciliation-task-contract.md`
- `docs/memory/phases/e8-full-analysis-roadmap-memory-reconciliation-fix1-task-contract.md`

Do not rewrite their historical bodies.

Add a short top-level supersession/correction note to each stating:

- the document is retained as historical implementation context;
- any requirement that includes empirical residual Fit 1/Fit 2 in the active E.8/F.6.3
  chain is superseded;
- the current authority is
  `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
  as corrected by this task;
- the accepted active profile is `no_empirical_residual`.

---

## 8. New durable correction record

Add:

`docs/memory/decisions/e8-no-empirical-residual-chain-correction.md`

It must record:

- correction date/context;
- source observation HEAD `942444fb9586dc53e847039c97c6c3b52fe9a351`;
- current `background_config.py` facts:
  - `BG_STAT_SCALE1 = 0.0`;
  - `BG_STAT_SCALE2 = 0.0`;
  - `BG_OPT_ACTIVE_PROFILE = "no_empirical_residual"`;
  - active profile forces both scales to zero;
- legacy Fit-1/Fit-2 blocks are conditional dormant source machinery;
- user scientific clarification that those historical empirical fits predate the
  current pion-model method and are not part of the accepted current chain;
- corrected chain:
  random/dummy -> proton -> baseline pion -> Method-A reweighting -> cleanest kaon
  `MM(t,phi)` -> yields;
- Method A is reweighting on fixed canonical bins, not rebinning;
- no new farm/runtime validation claim.

Use project status discipline; this is a decision correction, not a new scientific
runtime closure.

---

## 9. `docs/memory/manifest.json`

Regenerate after all memory edits and the new files.

---

## 10. Track this task contract

Add:

`docs/memory/phases/e8-no-empirical-residual-memory-correction-task-contract.md`

to the memory manifest.

---

# Frozen files / forbidden edits

This task must modify only `docs/memory/`.

Do not modify:

- `src/`
- `testing/`
- `main.py`
- `run_Prod_Analysis.sh`
- validation profiles
- collectors
- wrappers
- configuration
- analysis JSON
- PDFs
- accepted evidence artifacts
- root `AGENTS.md`
- `.codex/`

No physics/runtime behavior changes.

---

# Scientific boundaries to preserve

- current accepted production profile remains `no_empirical_residual`;
- legacy Fit 1 / Fit 2 are not part of the accepted current chain;
- baseline random/dummy, proton, and pion treatments remain separate owners;
- baseline pion contribution remains `b_j^0 = s_j*w0_j`;
- Method-A candidate remains `b_j^A = s_j*w0_j*C_j`;
- canonical t/phi bins remain frozen;
- no child independent renormalization;
- Method B remains diagnostic/cross-check only and never numerical;
- no Method-A production promotion before F.6.4;
- F.1-F.6.2 closed runtime statuses remain unchanged;
- accepted F.6.2 JSON/fingerprints remain frozen;
- source review is not farm validation.

---

# Deterministic local validation

Run the established repository memory tools for:

- manifest regeneration;
- manifest check;
- memory integrity/health;
- memory bootstrap check;
- `git diff --check`.

Also verify:

```bash
git status --short
git diff --name-only
git diff --stat
```

Every changed path must be under `docs/memory/`.

Search the resulting current forward-roadmap records for stale active-chain language,
including at minimum:

```bash
grep -RniE 'Fit[ -]?1|Fit[ -]?2|empirical residual' \
  docs/memory/CURRENT.md \
  docs/memory/MEMORY.md \
  docs/memory/decisions/e8-full-analysis-procedure-roadmap.md \
  docs/memory/phases/phase-f6-method-a-production-promotion.md \
  docs/memory/roadmap/STATUS.md \
  docs/memory/phases/PHASE_HISTORY.md
```

Any remaining Fit-1/Fit-2 occurrence in a current-authority record must either:

- explicitly say it is legacy/dormant/disabled and excluded from the active chain; or
- be removed.

Historical task contracts may retain their old body only with the required supersession
note.

No ROOT/PyROOT or farm execution.

---

# Review bundle

Create root-level temporary:

`kaonlt_review.diff`

containing the complete tracked diff plus complete no-index diffs for every new/untracked
memory file, including this contract and the new decision record.

Do not add `kaonlt_review.diff` to Git.

---

# Acceptance criteria

PASS only if:

1. HEAD remains exactly `942444fb9586dc53e847039c97c6c3b52fe9a351`;
2. only `docs/memory/` files change;
3. E.8.2 no longer includes Fit 1 or Fit 2;
4. E.8.2 baseline chain is:
   prompt/random -> dummy -> proton -> baseline pion -> final baseline MM/yields;
5. E.8.3 retains explicit `w0 -> w0*C`, baseline-vs-reweighted pion background,
   fixed-bin redistribution, and parent closure;
6. F.6.3 does not activate or rerun legacy empirical fits;
7. E.8.4 goes directly from baseline/Method-A pion treatment to final MM/yield
   comparison without Fit 1/Fit 2;
8. the final E.8 closure inventory contains no active empirical-residual-fit stage;
9. MEMORY records the durable `no_empirical_residual` boundary;
10. prior mistaken contracts are clearly superseded in part;
11. statuses/dependencies remain otherwise unchanged;
12. manifest/memory checks pass;
13. no farm/runtime claim is added.

---

# Hard stop

Stop and report rather than editing if:

- current source no longer has `BG_OPT_ACTIVE_PROFILE = "no_empirical_residual"`;
- the active profile no longer forces both empirical scales to zero;
- correcting memory would require source changes;
- any unrelated worktree change exists;
- any requested edit would alter accepted F.1-F.6.2 scientific evidence.

Codex must not commit, push, or run the farm.
