# E.8.1.Debug.1 — Left/lowe launcher debug mode

## Status and starting state

`ACTIVE` — local implementation and deterministic source-contract checks are
complete at starting `test` / `origin/test` HEAD
`16d595edabd2a348f51b728bafe9731e17b430c0`. This is a user-requested
operational prerequisite before the existing E.8.1 canonical-five-setting farm
gate: it makes one full `Left / lowe` diagnostic analysis practical without
starting the actual full high-epsilon analysis.

No `SOURCE REVIEWED`, ROOT/PyROOT, full-analysis, farm, or runtime-validation
claim is made here. Independent ChatGPT source review, the user-controlled Git
handoff, and any separately requested Jefferson Lab execution remain pending.

## Exact implementation scope

The only production-source changes are:

- `run_Prod_Analysis.sh` — adds short option `-d`, documents its narrow kaon
  `Left / lowe` meaning, preserves the existing flagged Q2/W convention,
  rejects non-kaon debug use after particle selection, scopes
  `LT_ANALYSIS_DEBUG_LEFT_LOW=1` to the actual full low-epsilon `main.py`
  invocation, propagates a nonzero low result, and exits successfully before
  the full high-epsilon invocation only after a successful debug low run.
- `src/main.py` — retains ordinary `phisetlist = ["Center", "Left", "Right"]`;
  it accepts the dedicated enabled control only for kaon low-epsilon full
  analysis, rejects unsupported values/contexts (including canonical-capture
  mode), and then narrows only that full-analysis population to `["Left"]`.
- `testing/test_run_prod_analysis_debug_left_low.py` — new static/
  source-contract coverage which does not run the production launcher or its
  cleanup path.

No general setting selector, scientific calculation, production correction,
subtraction, Method-A/Method-B ownership, SIMC, yields, acceptance/binning
algorithm, canonical helper, profile, wrapper, or presentation source changed.

## Paired-canonical and ordinary-path preservation

The launcher retains its ordinary kaon order:

`low raw capture -> high raw capture -> shared canonical preflight -> full low -> full high`.

Neither raw capture receives `LT_ANALYSIS_DEBUG_LEFT_LOW`; both retain their
ordinary full setting population. The selector reaches only the subsequent full
low call. Source tracing of `resolve_canonical_analysis_bins_pre_subtraction()`
shows the full Left/lowe run validates and reuses the paired interval text and
sidecars published by the shared preflight (`validated_authoritative_interval_file`);
that branch does not write interval files. Therefore the singleton debug
population does not regenerate or replace the paired canonical t/phi result.

With `-d` absent, the new environment control is absent, `phisetlist` remains
ordinary, and the pre-existing full low and full high calls remain on their
ordinary paths.

## Deterministic local checks

All following local checks returned zero after implementation:

- `python -m py_compile src/main.py testing/test_run_prod_analysis_debug_left_low.py`
- `C:\Program Files\Git\bin\bash.exe -n run_Prod_Analysis.sh`
- `python -m unittest testing.test_run_prod_analysis_debug_left_low -v` — 8
  tests passed, covering option/help/state separation, flagged Q2/W parsing,
  paired-capture ordering and absence of debug control, shared preflight before
  full low, Left-only full low control, failure propagation, intentional stop
  before full high, ordinary low/high paths, and a selector with no canonical
  configuration mutation.

After the required memory updates, the integrity manifest write/check,
`tools/check_memory_health.py --root .`, and
`tools/memory_bootstrap.py --root . --json` also passed. The applicable
`python -m unittest testing.test_memory_health -v` suite passed all 35 tests,
and `git diff --check` returned zero (the local Git configuration emits
non-failing CRLF conversion warnings while inspecting the Windows worktree).

The default Windows `bash.exe` is a denied WSL shim in this local host; the
focused test intentionally uses the available Git Bash executable for its
non-executing syntax check. No production launcher, `git clean -fdx`, ROOT,
PyROOT, or farm command was run.

## Changed-file accounting and next steps

Codex modified these files for this local implementation:

- `run_Prod_Analysis.sh`
- `src/main.py`
- `testing/test_run_prod_analysis_debug_left_low.py`
- `docs/memory/CURRENT.md`
- `docs/memory/phases/e8-1-debug-1-left-low-debug-mode.md`
- `docs/memory/manifest.json` (regenerated integrity-only index)

`docs/memory/phases/e8-1-debug-1-left-low-task-contract.md` was the
pre-existing user-created uncommitted contract. Codex did not edit it; the
integrity manifest inventories it because it is a nonignored memory file.

Next: obtain independent ChatGPT source review of the actual diff, then the
user may commit and push the scoped implementation. Farm execution remains the
user's responsibility and requires fresh artifacts before any runtime claim.
The separate E.8.1 scientific gate remains `DEVELOPMENT COMPLETE, FARM
VALIDATION PENDING`; its required first detailed review is still
`Q4p4W2p74 / Left / highe` after this operational prerequisite.
