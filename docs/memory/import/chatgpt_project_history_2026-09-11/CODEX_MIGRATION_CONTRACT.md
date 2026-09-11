# KaonLT Codex contract — migrate full historical project memory

## Task type

Documentation/durable-memory migration only.

Do **not** modify analysis source, tests, physics configuration, validation profiles, collectors, or production behavior.

Do not run the farm.

## 1. Repository and import

Repository:

```text
trottar/lt_analysis
branch: test
```

Historical import:

```text
docs/memory/import/chatgpt_project_history_2026-09-11/
```

Existing active memory:

```text
AGENTS.md
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/handoffs/CURRENT_HANDOFF.md
docs/memory/roadmap/CURRENT.md
```

This task reconciles the imported history into canonical durable memory.

## 2. Establish live source first

Before editing memory:

```bash
git status --short --branch
git rev-parse HEAD
git rev-parse origin/test
```

Read the active-memory files above and every file in the import directory.

The import was prepared when `test` was observed at:

```text
7cdb7847d72501cd3dc504565fafc604e22a6132
```

Do not assume it is still current.

Inspect current source/Git history wherever needed to resolve implementation-state claims.

## 3. Preserve evidence classes

Keep these distinct:

```text
[CURRENT_SOURCE]
[GIT_HISTORY]
[CHAT_HISTORY]
[HANDOFF]
[FARM_EVIDENCE]
[SCIENTIFIC_REFERENCE]
[INFERENCE]
[SUPERSEDED]
```

For implementation:

```text
current source
> newer exact source/diff
> newest authoritative handoff
> older history
> inference
```

For runtime:

```text
fresh farm evidence
> reviewed authoritative farm handoff
> source evidence
> historical summary
> inference
```

Never assign `CLOSED / RUNTIME VALIDATED` from a commit subject, source existence, unit tests, or collector alone.

## 4. Preserve scientific ownership

Keep strict separation among:

```text
random subtraction
slow-proton subtraction
pion-background treatment
HGCer Method A
HGCer Method B
SIMC comparisons
yield extraction
cross-section analysis
diagnostics/checkers
presentation-only plotting
```

Do not modify cuts, normalizations, windows, templates, priors, component definitions, binning, production formulas, efficiencies, acceptance, L/T separation, or uncertainties.

## 5. Build deep memory

Create/update only documentation under `docs/memory/`.

Use, where useful:

```text
docs/memory/
├── CURRENT.md
├── MEMORY.md
├── handoffs/CURRENT_HANDOFF.md
├── history/
│   ├── PROJECT_HISTORY.md
│   └── CHAT_INDEX.md
├── phases/
├── decisions/
├── evidence/
├── investigations/
├── sources/
│   ├── SOURCE_INDEX.md
│   └── ARTIFACT_INDEX.md
└── roadmap/CURRENT.md
```

Do not create empty filler files.

Migrate/reconcile:

```text
PROJECT_HISTORY.md
CHAT_INDEX.md
PHASE_HISTORY.md
VALIDATION_HISTORY.md
DECISION_HISTORY.md
SOURCE_INDEX.md
ARTIFACT_INDEX_SEED.md
KNOWN_GAPS.md
```

Keep the import directory as source provenance if useful.

Retain rejected approaches, failed farm gates, review mistakes, superseded plans, and corrective decisions. Do not rewrite history to make it appear linear.

## 6. Required decision records

Ensure durable records exist for at least:

1. slow-proton versus pion-background separation;
2. proton-before-pion production ordering;
3. K-Lambda preservation gate;
4. proposed versus applied proton semantics;
5. fixed-envelope pion-alignment comparison;
6. Method-A/Method-B independence;
7. same-`t` Method-B normalization;
8. adaptive Method-B `DO NOT PROMOTE`;
9. persisted-diagnostic round-trip review;
10. Phase-E presentation-only ownership;
11. generic collector + JSON profile architecture;
12. Method B diagnostic-only future role;
13. Method-A zero-NPE/absolute-probability limitation;
14. Method-A training/application separation;
15. parent-level normalization / no `(t,phi)` child renormalization;
16. no production promotion before F.6.

## 7. Required farm evidence records

### Phase C five-setting closure

Record:

```text
C.Fix.2.3:
9a66bc62d20a99172e326e915866877b65ae1e5d

later accepted pre-E.3:
e3853655db0809923cbf2326e2f779219128eda9

settings:
Left lowe
Left highe
Center lowe
Center highe
Right highe

result:
Phase C PASS
adaptive Method B DO NOT PROMOTE
30/30 pages rendered
```

Preserve the support statistics and scientific reason for non-promotion.

### E.3.Fix.2 Left-low

Record:

```text
implementation:
eb1710f4739ba6ef14f51419806e9fc5bd53c175

runtime/bundle HEAD:
bf53dac84e1396cfd7e3f4e0234749426bfbdcf4

gate:
Q4p4W2p74 Left-low

result:
CLOSED / RUNTIME VALIDATED
```

Include recorded source checks and rendered-page review.

Do **not** upgrade the remaining four E.3 settings without evidence.

### F.1 Fix.3/Fix.4 mechanics

Record:

- NumPy edge-array truth-value farm failure and repair;
- NumPy bool-scalar parity farm failure and repair;
- later five-setting mechanical run and focused test counts.

Do not interpret that gate as scientific acceptance of F.1 v1.

### F.1 population-ownership failure

Record the Method-A positive/low/control counts and F.1 `low=0` observation. Explain why this required Fix.5.

## 8. Reconcile F.1.Fix.5 against live source

Read current relevant files, including:

```text
src/cuts/pion_hgcer_method_a_acceptance_contract.py
src/cuts/rand_sub.py
src/cuts/full_background_subtraction_plots.py

testing/test_pion_hgcer_method_a_acceptance_contract.py
testing/test_pion_hgcer_phase_f_runtime_contract.py
testing/test_full_background_subtraction_plots.py
```

Check whether source actually implements:

```text
separate Method-A training records
separate application records
training from prompt/noRF/nommcuts NPE>0 records
application restricted to authoritative physical NPE>2 cache
separate summaries/fingerprints
pages 1-4 from training
page 5 from application
no Method-B numerical dependency
no correction/estimator/production mutation
```

If yes, assign only source-level status unless fresh farm evidence exists.

If not, document the discrepancy.

## 9. Respect known gaps

Read `KNOWN_GAPS.md`.

Do not invent:

- exact final Phase-D farm provenance;
- final remaining-four-setting E.3 closure;
- F.1.Fix.5 farm PASS;
- scientific meaning from generic June commit subjects.

If repository-owned evidence resolves a gap, cite that exact record and close only that gap.

## 10. Rebuild active memory after deep migration

Only after deep records are reconciled, update:

```text
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/handoffs/CURRENT_HANDOFF.md
docs/memory/roadmap/CURRENT.md
```

### CURRENT.md

Keep concise:

- live branch/HEAD;
- current phase/fix;
- supported source status;
- supported runtime status;
- frozen upstream architecture;
- unresolved/deferred gaps;
- one exact `NEXT`.

Do not retain an obsolete bootstrap `NEXT`.

### MEMORY.md

Curated durable rules/facts only. Chronology belongs in deep memory.

### CURRENT_HANDOFF.md

Enough for a completely new Codex chat to continue immediately. Link deep records instead of repeating them.

### roadmap/CURRENT.md

Preserve the approved Phase-F roadmap unless newer authoritative evidence supersedes it:

```text
F.1 event contract
F.2 representation freeze
F.3 detached Method-A map
F.4 detached parent-preserving A-only correction
F.5 detached event propagation
F.6 production promotion
```

## 11. Allowed status labels

Use only:

```text
CLOSED / RUNTIME VALIDATED
SOURCE REVIEWED
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
ACTIVE
DEFERRED
BLOCKED
NEXT
```

Historical evidence records may quote scientific/checker `PASS`, but canonical current state uses the labels above.

## 12. Audit

Before finishing:

```bash
git diff --check
git status --short
git diff --stat
git diff --name-only
```

Verify:

1. changed paths are under `docs/memory/` unless a strictly necessary `AGENTS.md` documentation update is justified;
2. no analysis source changed;
3. no test code changed;
4. no collector/profile changed;
5. no physics configuration changed;
6. live source identity is correct;
7. runtime statuses have farm provenance;
8. superseded decisions are marked;
9. known gaps remain explicit;
10. active memory is concise compared with deep history.

Do not commit unless the user explicitly asks.

## 13. Required final report

Report:

1. live starting HEAD;
2. files created/updated;
3. deep-memory records added;
4. canonical statuses and exact evidence for every runtime-validated status;
5. remaining known gaps;
6. whether live Fix.5 source matches the imported contract;
7. whether Fix.5 has farm evidence;
8. exact resulting `NEXT`;
9. `git diff --stat`;
10. confirmation no production/source/test/config file changed.

Then stop.
