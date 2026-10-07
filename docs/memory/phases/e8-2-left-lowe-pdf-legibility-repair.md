# E.8.2 Left/lowe PDF-legibility repair

Status: `CLOSED / RUNTIME VALIDATED` for Q4p4W2p74 / Left / lowe visual
legibility only. The [fresh closure evidence](../evidence/e8-2-left-lowe-runtime-scientific-closure-2026-10-06.md)
records prior ChatGPT acceptance of all 18 rendered E.8.2 pages at farm source
`706ae708ae69be0585472982a4df3d87c63f186f`. The implementation chronology and
then-pending acceptance boundary below remain historical.

## Authority and cause

Implemented the [contract](e8-2-left-lowe-pdf-legibility-repair-task-contract.md)
on branch `test`, starting HEAD/local `origin/test`
`8ae17e18e480527cf470a93a7a8fb6f0747308fd`. The
[accepted evidence record](../evidence/e8-2-left-lowe-runtime-owner-pass-visual-blocker-2026-10-06.md)
owns exact external artifact identities and prior ChatGPT acceptance: the
Left/lowe owner/isolation/provenance gate passed, but all 18 E.8.2 pages failed
visual acceptance through subtraction clipping, header overlap or text truncation.
Codex records those supplied facts without locally revalidating the artifacts.

SOURCE VERIFIED: subtraction pages used a 1800-by-3240 nine-row canvas;
the final-MM grid shared full-canvas coordinates with its header; stage-yield
text concatenated six long persisted names/values in one line. These are
presentation causes; no scientific algorithm or payload repair is justified.

## Renderer repair

Only E.8.2 presentation helpers in `src/cuts/full_background_subtraction_plots.py`
change. Each page uses a bounded 2400-by-2400 square canvas with a header pad
at y=0.90–1.00 and a distinct content pad at y=0.02–0.89. Header text is retained
in its own region, following the existing E.8 canvas/header/grid pattern.

The subtraction content retains 3x9 random/dummy/pion and 4x9 proton grids,
canonical child order, exact stage columns, shared per-child signed display
ranges, Lambda markers, phi labels and unavailable behavior. The proton page
still identifies the production prune boundary between columns 3 and 4.
Final-MM content retains the 3x3 grid; phi context and stored Y0/stat/total
are displayed in two short panel-note rows.

The pure stage-yield formatter returns three compact rows per valid child:
PP/AR/AD, PRE/POST/PI, then authoritative Y0/stat/total. A visible three-line
legend maps every abbreviation to its exact persisted stage name. The
diagnostic-window-integral warning remains visible. The formatter consumes
stored values with existing numeric precision and performs no integration,
fit, subtraction or yield extraction.

All 18 page IDs, semantic identities and both E.8.2 schemas are unchanged.
No continuation page is added. All other scientific functions, owner/profile/
collector paths, Method A/B and production interfaces remain frozen.

## Deterministic checks

Windows Python 3.12.10:

- `python -B -m py_compile src/cuts/full_background_subtraction_plots.py
  testing/test_e8_2_baseline_stage_audit.py`: PASS; cache directed to a temporary
  directory outside the repository.
- `python -B -m unittest testing.test_e8_2_baseline_stage_audit
  testing.test_full_background_subtraction_plots -v`: PASS, 128 tests,
  18 existing presentation-policy skips. This comprises 26 E.8.2 tests (no
  skips) and 102 existing presentation tests (18 skips).

Ubuntu WSL Python 3.12.3, used for existing symlink support:

- `python3 -B -m unittest testing.test_run_e8_2_left_lowe_scientific_audit_gate -v`:
  PASS, 27 tests, no skips; owner/test unchanged.

New fake-ROOT tests capture separate bounded pads, all nine children and every
column, Lambda markers, exact manifest IDs/semantics, retained prune wording,
all diagnostic/final-yield fields and abbreviation definitions. Extreme-value
formatting remains bounded. Source histogram contents/errors/directory/display
attributes are unchanged; Integral raises if called by the renderer. Existing
producer, failure/finalization and recovery tests continue to pass.

## Memory and acceptance boundary

Changed paths are renderer, E.8.2 test, CURRENT, manifest, supplied contract,
accepted-evidence record and this phase record only. CURRENT distinguishes
accepted narrow owner/isolation runtime from the still `BLOCKED` E.8.2 visual
audit and names the required fresh visual gate as the substantive NEXT.

RUNTIME VERIFIED: only the previously accepted supplied owner/isolation facts;
no new runtime acceptance follows from this source repair.
NOT VERIFIED: repaired physical PDF legibility and scientific stage/prune
interpretation. Local fake-ROOT checks establish source behavior only.
No ROOT/PyROOT, production analysis or farm command ran; no staging, commit,
push or ref update occurred. Actual-diff review, user synchronization and
farm-readiness review precede one fresh isolated owner run. All 18 pages must
then pass rendered visual inspection before scientific interpretation.
