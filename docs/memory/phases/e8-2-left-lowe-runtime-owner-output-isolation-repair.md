# E.8.2 Left/lowe owner OUTPUT-isolation repair

Status: `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING` for repaired owner
infrastructure only; actual-diff review is pending.

## Authority and failed gate

Implemented the [repair contract](e8-2-left-lowe-runtime-owner-output-isolation-repair-task-contract.md)
on `test` at starting HEAD/local `origin/test`
`f860f832f3e4b05f468cb6bb636f82f4a26d726b`. The
[original owner record](e8-2-left-lowe-runtime-scientific-audit-owner.md) retains
its historical local implementation scope; this record owns the topology repair.

RUNTIME VERIFIED (supplied terminal output recorded in the contract): the first
isolated attempt reached Center DiamondPlot, reported `No valid file found!`,
and exited nonzero. The owner stopped at `analysis_failed`. This failed attempt
is `BLOCKED` at owner isolation and is diagnostic evidence, not E.8.2 acceptance.
No farm artifact was independently inspected or farm command executed here.

SOURCE VERIFIED: the copied-ltsep overlay makes child OUTPATH lexical to the
detached analysis worktree. The launcher creates OUTPATH before link setup;
debug mode skips `git clean -fdx`. A fresh worktree therefore acquired a real
local OUTPUT directory, which established `set_SymLinks.sh` does not replace
with the external OUTPUT symlink. Center inputs were unavailable through that
local topology. The repair belongs to the new owner, not Diamond physics or
the mature production launcher.

## Narrow repair

SOURCE VERIFIED: `prepare_debug_output_link(worktree, paths, outdir)` consumes
the existing probe dictionary and validated artifact root. It checks absolute
worktree/volatile authority, matching LTANAPATH, exact lexical OUTPATH,
existing external OUTPUT directory, selected artifact-root agreement and
absence of any worktree OUTPUT path, including dangling symlinks. It creates
one symlink in the disposable worktree and verifies its literal/resolved target
and the resolved child OUTPATH. Unexpected paths are never deleted/replaced.

Owner ordering is overlay -> path probes -> artifact-root authority -> OUTPUT
link -> external symlink preflight -> one debug child. The verified link record
is persisted before analysis starts and included in the run summary. Existing
cleanup/preservation and collection/ZIP gates remain unchanged.

Changed versioned paths are the E.8.2 owner/test, this record, the supplied
repair contract, CURRENT and generated manifest. The launcher, link script,
all scientific/production source, collector, E.8.2 profile, canonical-five
owner and other records remain byte-identical to the startup baseline.

## Deterministic local checks

- Windows Python 3.12.10: `python -B -m py_compile
  testing/run_e8_2_left_lowe_scientific_audit_gate.py
  testing/test_run_e8_2_left_lowe_scientific_audit_gate.py`: PASS; cache directed
  to a temporary directory outside the repository.
- Ubuntu WSL Python 3.12.3: `python3 -B -m unittest
  testing.test_run_e8_2_left_lowe_scientific_audit_gate -v`: PASS, 27 tests,
  no skips. Real temporary symlinks prove launcher-prefix mkdir preserves
  indirection and a harmless Center/kaon-like fixture remains visible. Negative
  cases preserve pre-existing directory/file/correct/wrong/dangling links and
  reject mismatched/missing authority. Flow tests assert exact ordering,
  pre-child status provenance, summary provenance and fail-closed preservation.
- Windows: `python -B -m unittest testing.test_e8_2_baseline_stage_audit -v`:
  PASS, 23 tests, unchanged suite.
- Windows: `python -B -m unittest
  testing.test_run_e8_4_fix5_canonical_five_plot_gate -v`: PASS, 26 tests.
- Ubuntu WSL: `python3 -B -m unittest
  testing.test_collect_pion_hgcer_validation_bundle -v`: PASS, 28 tests.

Windows cannot create symlinks here (`WinError 1314`); the full owner suite
therefore runs in the existing local Ubuntu environment without weakening tests
or changing system privileges. WSL lacks NumPy for the isolation regression;
that suite passes with the existing Windows dependencies. Explicit fixture
timestamps avoid filesystem-clock rounding; runtime freshness checks are frozen.
No local check runs ROOT/PyROOT, the production analysis or the farm.

## Remaining evidence boundary

NOT VERIFIED: repaired farm integration, actual E.8.2 stage magnitudes,
production-pruning impact and PDF legibility. Local topology tests do not close
the failed farm gate or establish scientific acceptance. Method A remains
detached/non-production; Method B diagnostic/cross-check only and numerically
excluded. Canonical-five provenance repair remains `DEFERRED`; final E.8/F.6.4
and absolute-SIMC blockers retain their scopes.

CURRENT owns the sole NEXT: actual-diff review, user-controlled commit/push,
pushed-state synchronization/farm-readiness review, then exactly one isolated
E.8.2 Left/lowe owner rerun and fresh evidence review. Codex stops at the complete
review bundle and issues no farm command.
