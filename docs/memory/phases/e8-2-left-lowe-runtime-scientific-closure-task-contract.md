# KaonLT E.8.2 Left/lowe runtime/scientific closure — task contract

## 1. Task class and authority

This is a **closure/documentation-only source task** for the accepted
Q4p4W2p74 / Left / lowe E.8.2 baseline full-analysis audit.

It consumes already-supplied and already-reviewed farm evidence. It must not
rerun the farm, alter analysis physics, reopen accepted phases, or reinterpret
Method A/Method B ownership.

Before editing, read repository-root `AGENTS.md`, then the startup core in this
exact order:

1. `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read only task-relevant records:

- `docs/memory/MAINTENANCE.md`
- `docs/memory/CODEX.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/phases/e8-2-left-lowe-pdf-legibility-repair.md`
- `docs/memory/evidence/e8-2-left-lowe-runtime-owner-pass-visual-blocker-2026-10-06.md`

Authority remains:

```text
current source/diff
-> newest applicable farm/runtime evidence
-> newest applicable validation/handoff artifact
-> tracked durable memory
-> older chat/history
```

This contract records the prior ChatGPT acceptance. Codex does not claim to
have independently rerun or farm-validated the artifacts.

---

## 2. Required starting identity

Repository:

```text
https://github.com/trottar/lt_analysis/tree/test
```

Required branch:

```text
test
```

Required starting HEAD and local `origin/test`:

```text
706ae708ae69be0585472982a4df3d87c63f186f
```

Commit:

```text
Repair E8.2 PDF geometry and text legibility
```

Before editing, establish:

```bash
git branch --show-current
git rev-parse HEAD
git rev-parse refs/remotes/origin/test
git status --short --untracked-files=all
git diff --check
```

The only expected new task file before implementation is:

```text
docs/memory/phases/e8-2-left-lowe-runtime-scientific-closure-task-contract.md
```

Temporary root-level `kaonlt_review*.diff` files are permitted only as review
artifacts.

If branch/HEAD/local origin differs, or unrelated workstation state is present,
stop as `BLOCKED`. Do not reset, stash, clean, overwrite, or repair unrelated
user state.

No farm command is authorized in this closure task.

---

## 3. Accepted fresh runtime evidence

Package stem:

```text
KaonLT_E8_2_Left_lowe_scientific_audit_20261006-233426
```

Farm-evaluated source:

```text
706ae708ae69be0585472982a4df3d87c63f186f
```

Scope:

```text
Q4p4W2p74 / Left / lowe only
```

External artifact identities accepted by ChatGPT:

```text
ZIP
SHA-256:
464d0db14af4e09a53a418dc077136fa17f1efd1c4c36b57758247b8f083af23
bytes:
30272198

gate-status JSON
SHA-256:
900f32f5077754df27439c2df3c4243c3d4e9fc8244a399943616153a4089c1c
bytes:
21073

run-summary JSON
SHA-256:
d9d265bcb80f82f4bb3fff7abb243913abb2b20292b71a5d7b58f8796087400e
bytes:
17394

child log
SHA-256:
4e8713d8f4c399063a6feb43395caf1a581f6154c06660a3888b81561f297778
bytes:
51888750
```

Accepted packaged artifacts:

```text
Left_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction.pdf
SHA-256:
aacac5daa935eb1f2af6dacdf0ddfdcdd4f5b202e0d54fc11b8edd254576e61e
bytes:
1551786
pages:
74

Left_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction-manifest.json
SHA-256:
e51fe03cff97ce00e05c850b031206e97a9cb9d15b0c71279618fda1a9577f02
bytes:
18203

kaon_FullAnalysis_Q4p4W2p74_lowe.json
SHA-256:
6e6f30c630de92a138520bfe3a93612357d8b176211454c97f22bbb94b9e328e
bytes:
140864222

kaon_FullAnalysis_Q4p4W2p74_lowe_correction_ledger_no_empirical_residual.json
SHA-256:
cba9737c69868b8cb6a9ecbfba08bc6a25a7e326445298fd859cc52f403b9f92
bytes:
53424

kaon_FullAnalysis_Q4p4W2p74_lowe_correction_ledger_no_empirical_residual.csv
SHA-256:
3eb3aa22419e121140edb5829084f8401e9986cde91790d8c221751380c63278
bytes:
20343
```

RUNTIME VERIFIED by the accepted supplied evidence and prior ChatGPT review:

- owner `status=success`;
- owner `stage=complete`;
- null failure reason;
- child analysis return code 0;
- analysis completed;
- artifact verification completed;
- collection completed;
- ZIP verification completed;
- worktree cleanup completed;
- repaired disposable-worktree OUTPUT link passed;
- ordinary farm checkout preservation passed before and after;
- installed `ltsep` preservation passed before and after;
- source commit and `origin/test` matched
  `706ae708ae69be0585472982a4df3d87c63f186f`;
- only the two established allowlisted farm outputs were dirty in the ordinary
  checkout;
- all 18 required E.8.2 page IDs were structurally present;
- all 18 required E.8.2 rendered pages, PDF pages 47–64, were visually reviewed
  and accepted.

The repaired visual defects are closed:

- all nine phi children are physically visible on every subtraction page;
- no top-row clipping remains;
- final-MM headings/t context do not overlap panels;
- stage-yield text is fully readable;
- no right-edge truncation remains;
- no new renderer failure is visible.

Therefore:

```text
E.8.2 Q4p4W2p74 / Left / lowe:
CLOSED / RUNTIME VALIDATED
```

This is a narrow Left/lowe closure only. It does not establish canonical-five
E.8.2 closure, other-setting runtime acceptance, Method-A promotion,
absolute-SIMC amplitude correctness, or final E.8 closure.

---

## 4. Accepted E.8.2 scientific interpretation

The accepted stage-yield pages label their stage integrals explicitly as
**diagnostic Lambda-window integrals, not final extracted yields**.

For each canonical t parent, summing the nine phi cells gives the following
signed fractional change relative to the immediately preceding stage:

| t parent | Random | Dummy | Slow proton | `prune_hist` | Baseline pion |
| --- | ---: | ---: | ---: | ---: | ---: |
| t1, 0.4000–0.5667 GeV^2 | -0.68% | 0.00% | -7.75% | 0.00% | -20.95% |
| t2, 0.5667–0.7333 GeV^2 | -0.64% | -1.48% | -6.46% | 0.00% | -18.71% |
| t3, 0.7333–0.9000 GeV^2 | -0.44% | -2.59% | -8.08% | -0.03% | -4.81% |

Accepted interpretation:

- random subtraction is small in this scope, below 1% for all three parents;
- dummy subtraction is modest, reaching about 2.6% at t3;
- slow-proton cleaning is a consistent, material correction of about 6.5–8.1%;
- baseline pion subtraction is the largest diagnostic-window reduction at t1
  and t2, about 19–21%;
- at t3 the baseline pion effect is about 4.8%, smaller than the slow-proton
  effect;
- production `prune_hist` is essentially inert in the Lambda-window diagnostic.

The only visible nonzero prune change occurs at t3 / phi5:

```text
PRE  = 0.00018539
POST = 0
```

This is about 0.03% of the t3 pre-prune signed parent sum and about 0.014%
across all three pre-prune parent sums.

Preserve the interpretation boundary:

- `PI` is the diagnostic Lambda-window integral of the after-pion histogram;
- authoritative `Y0` is the stored final extracted yield;
- differences between `PI` and `Y0` are not a new subtraction stage;
- no empirical-residual stage is present under the active
  `no_empirical_residual` profile;
- small near-zero or negative final cells are observations of the baseline
  result, not by themselves evidence for a production failure or a redesign.

Do not broaden these conclusions beyond Q4p4W2p74 / Left / lowe.

---

## 5. Scientific and production boundaries

This closure task must not alter any scientific or runtime source.

Freeze all:

```text
src/
testing/
tools/
farm_env/
background_samples/
run_Prod_Analysis.sh
set_SymLinks.sh
```

except no file in those trees is allowed to change at all in this task.

Specifically preserve:

- random/dummy subtraction;
- slow-proton treatment;
- `prune_hist`;
- baseline pion treatment and `w0`;
- final baseline spectra/yields/errors;
- active `no_empirical_residual`;
- SIMC;
- cuts/binning/efficiencies/acceptance;
- L/T separation and cross sections;
- Method A detached/non-production role;
- Method B diagnostic/cross-check-only role;
- canonical-five provenance repair as `DEFERRED`.

No source code, profile, owner, collector, renderer, or test changes are
authorized.

---

## 6. Exact allowed versioned paths

Only these seven paths may change:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/roadmap/STATUS.md
docs/memory/evidence/e8-2-left-lowe-runtime-scientific-closure-2026-10-06.md
docs/memory/phases/e8-2-left-lowe-pdf-legibility-repair.md
docs/memory/phases/e8-2-left-lowe-runtime-scientific-closure-task-contract.md
docs/memory/phases/e8-2-left-lowe-runtime-scientific-closure.md
```

No other versioned path may change.

---

## 7. Required evidence record

Create:

```text
docs/memory/evidence/e8-2-left-lowe-runtime-scientific-closure-2026-10-06.md
```

It must record:

- farm source `706ae708ae69be0585472982a4df3d87c63f186f`;
- exact scope Q4p4W2p74 / Left / lowe;
- exact external artifact identities in section 3;
- packaged artifact identities in section 3;
- owner/provenance/artifact/page/ZIP acceptance;
- 18-page rendered visual acceptance;
- E.8.2 status `CLOSED / RUNTIME VALIDATED`;
- the stage table and interpretation boundaries in section 4;
- explicit narrow-scope exclusions:
  - no canonical-five closure;
  - no Method-A promotion;
  - no absolute-SIMC claim;
  - no final E.8 closure.

Attribute runtime/scientific acceptance to the supplied artifacts and prior
ChatGPT review, not to Codex.

---

## 8. Required phase closure record

Create:

```text
docs/memory/phases/e8-2-left-lowe-runtime-scientific-closure.md
```

Record:

- accepted evidence record;
- source identity;
- narrow closure status;
- closure of the PDF-legibility repair;
- concise accepted scientific conclusions;
- no source/physics changes;
- exact next gate.

Update:

```text
docs/memory/phases/e8-2-left-lowe-pdf-legibility-repair.md
```

only enough to change its presentation repair status from:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

to:

```text
CLOSED / RUNTIME VALIDATED
```

for Q4p4W2p74 / Left / lowe visual legibility, and link the new closure
evidence. Do not rewrite its historical implementation chronology.

---

## 9. CURRENT update

Update `docs/memory/CURRENT.md` compactly.

Required active meaning:

```text
E.8 remains ACTIVE.

E.8.2 Q4p4W2p74 / Left / lowe is CLOSED / RUNTIME VALIDATED.

Fresh source:
706ae708ae69be0585472982a4df3d87c63f186f

The baseline audit confirms:
- random small (<1%);
- dummy modest (<= about 2.6%);
- slow proton material (about 6.5–8.1%);
- baseline pion largest at t1/t2 (about 19–21%), about 4.8% at t3;
- prune_hist negligible in the Lambda-window diagnostic.

These are narrow Left/lowe conclusions.

Canonical-five provenance/identity repair remains DEFERRED.
Final E.8 and F.6.4 remain BLOCKED.
Absolute-SIMC amplitude remains separately BLOCKED.
Method A remains detached/non-production.
Method B remains diagnostic/cross-check only and numerically excluded.
```

Remove consumed wording that says E.8.2 visual validation is pending/blocked.

The sole ordinary NEXT must now be:

```text
NEXT — audit the existing E.8.3 detached Method-A reweighting audit against
current source and accepted authorities, then determine its narrow
Q4p4W2p74 / Left / lowe runtime-validation readiness. Do not redesign Method A
or reopen closed E.8.2/F-stage work without new evidence.
```

This NEXT is an **audit/readiness gate**, not authorization to run the farm.
No E.8.3 farm command is issued by this closure task.

Keep CURRENT concise. It was near its soft limit before this task; consolidate
consumed E.8.2 blocker/repair prose rather than appending a long new section.

---

## 10. Roadmap status update

Update only the E.8.2 entry in:

```text
docs/memory/roadmap/STATUS.md
```

from `SOURCE REVIEWED` to a narrow statement equivalent to:

```text
CLOSED / RUNTIME VALIDATED — only for Q4p4W2p74 / Left / lowe at farm source
706ae708ae69be0585472982a4df3d87c63f186f. Fresh owner/provenance/artifact
and all 18 rendered pages passed. Baseline stage interpretation is accepted
only for this scope. No canonical-five closure or production promotion follows.
```

Link the new evidence record.

Do not change E.8.3, F.6.3, E.8.4, canonical-five, final E.8, or F.6.4 status in
this task.

Do not edit `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`; its
role is the approved dependency/scientific contract, not the live status owner.

---

## 11. Manifest and deterministic checks

No ROOT/PyROOT or analysis execution is needed.

Discover the local Python interpreter according to tracked
`docs/memory/TOOLS.md`.

Run at minimum:

```text
python compile/check only for memory tools if required by tracked procedure;
tools/check_memory_health.py --root .
manifest regeneration/check according to tracked procedure
bootstrap JSON according to tracked procedure
git diff --check
```

Verify:

- exactly the seven allowed versioned paths changed;
- all `src/`, `testing/`, owner/profile/collector/runtime files are
  byte-identical to starting HEAD;
- CURRENT has one authoritative NEXT;
- no status is upgraded beyond the accepted Left/lowe evidence.

Report:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
warning classification: <none or blocking/nonblocking>
manifest check: PASS | FAIL
bootstrap: PASS | FAIL
git diff --check: PASS | FAIL
```

---

## 12. Review bundle

Create one complete root-level:

```text
kaonlt_review.diff
```

for ChatGPT actual-diff review.

It must include:

- tracked diffs for modified files;
- complete no-index additions for all new task/evidence/phase files.

Do not stage merely to create the bundle.

---

## 13. Hard stops

Stop as `BLOCKED` instead of broadening scope if any of these appears necessary:

- any scientific/source/test/owner/profile/collector change;
- any E.8.3 implementation change;
- any farm command;
- any ROOT/PyROOT run;
- reopening E.8.2 science;
- changing canonical-five provenance status from `DEFERRED`;
- changing Method A/B ownership;
- claiming canonical-five closure;
- claiming Method-A production promotion;
- claiming absolute-SIMC amplitude validation;
- claiming final E.8 closure.

---

## 14. Acceptance state and Codex stop point

After successful local closure reconciliation:

```text
E.8.2 Q4p4W2p74 / Left / lowe:
CLOSED / RUNTIME VALIDATED

E.8:
ACTIVE

NEXT:
E.8.3 current-source / accepted-authority audit and narrow Left/lowe
runtime-readiness determination
```

Codex stops after:

- closure memory/evidence updates;
- deterministic memory/integrity checks;
- complete `kaonlt_review.diff`.

Codex must not:

- stage merely for review;
- commit;
- push;
- update refs;
- run ROOT/PyROOT;
- run production analysis;
- run the farm;
- begin E.8.3 implementation or runtime execution.

Return:

- starting/ending branch;
- starting/ending HEAD and local `origin/test`;
- exact changed paths;
- concise closure summary;
- exact checks/results;
- byte-preservation result for frozen source/runtime files;
- memory-health report;
- `kaonlt_review.diff` path/size;
- explicit confirmation that no source physics, ROOT/PyROOT, production
  analysis, or farm command ran.
