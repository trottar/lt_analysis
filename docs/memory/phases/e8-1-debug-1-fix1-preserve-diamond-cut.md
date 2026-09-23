# E.8.1.Debug.1.Fix.1 — preserve the ordinary Diamond cut

## Status and starting state

`SOURCE REVIEWED` — Fix.1 starts from `test` / `origin/test` committed HEAD
`4c4a18a098d184fe0b13dd1fc57010ac7f5b30c0`. Independent review of the
original E.8.1.Debug.1 source found that its selector narrowed `phisetlist`
before Step 2, creating the Diamond-cut ownership regression.

ChatGPT independently inspected the actual local Fix.1 diff and found the
repair PASS: only the downstream debug restriction moved, while the ordinary
Diamond stage remains intact. Codex-reported deterministic checks were `NOT RUN
by ChatGPT`. This remains no claim of ROOT/PyROOT, full-analysis, farm, or
runtime validation. User-controlled commit/push and targeted Jefferson Lab
debug evidence remain pending.

## Source-review finding and preserved ownership

The ordinary Step-2 `DiamondPlot()` path runs Center first over
`["Center", "Left", "Right"]`. Its Center call produces the common polygon as
`cut_mode="poly"` and `poly_points` in `inpDict`; subsequent Left processing
consumes that established cut. The original debug selector instead removed
Center before Step 2, leaving Left to derive its direct-call low-epsilon
fallback polygon. That changed a frozen scientific cut.

Fix.1 changes no Diamond algorithm, cut definition, acceptance, canonical
binning, subtraction, HGCer method, SIMC, background optimization, yield,
efficiency, L/T separation, uncertainty, or production correction. It leaves
`run_Prod_Analysis.sh`, `src/cuts/diamond.py`, canonical helpers, profiles,
collectors, and historical Debug.1 records unchanged. Analysis source outside
`src/main.py` is unchanged by Fix.1.

## Repair implementation and source trace

`src/main.py` retains early validation of `LT_ANALYSIS_DEBUG_LEFT_LOW`: only
explicit enabled values, kaon low-epsilon context, and no canonical-capture
misuse are accepted. That early block no longer changes `phisetlist`.

The ordinary Step-2 loop therefore runs Center, Left, and Right in its existing
order, registers all three ordinary Diamond PDFs, and lets Center publish the
common polygon into `inpDict`. Only after that registration, before the
`shift_prep` setting loop, Fix.1 applies `phisetlist = ["Left"]` when the debug
control is enabled. The remaining shift/SIMC/subtraction/binning/yield path is
therefore Left-only while consuming the ordinary Center-produced Diamond cut.

The frozen launcher still executes:

`low raw capture -> high raw capture -> shared canonical preflight -> full low -> full high`.

Neither capture receives the debug environment control. The full debug low call
receives it only after the paired preflight; canonical resolution continues to
validate and reuse the published authoritative interval pair without writing a
replacement. A successful Left/lowe return still causes the launcher to stop
before full high epsilon, and a failure still exits nonzero.

## Deterministic local checks

The following checks were reported as passed by Codex after implementation;
they were `NOT RUN by ChatGPT` during the independent actual-diff review:

- `python -m py_compile src/main.py testing/test_run_prod_analysis_debug_left_low.py`
- `C:\Program Files\Git\bin\bash.exe -n run_Prod_Analysis.sh`
- `python -m unittest testing.test_run_prod_analysis_debug_left_low -v` — 9
  tests passed, including the ordinary-Diamond-before-restriction ordering,
  downstream-before-`shift_prep` placement, unchanged `diamond.py`, paired
  unrestricted capture, scoped debug call, success stop, failure propagation,
  and ordinary low/high path checks.
- required manifest write/check, memory health, bootstrap, applicable
  `testing.test_memory_health`, and `git diff --check` are run after the memory
  updates below.

The default Windows `bash.exe` is a denied WSL shim on this host, so the
available Git Bash executable is used for the non-executing shell check. No
production launcher, `git clean -fdx`, ROOT/PyROOT, farm command, commit, or
push was run.

## Changed-file accounting and next steps

Codex changed only:

- `src/main.py`
- `testing/test_run_prod_analysis_debug_left_low.py`
- `docs/memory/CURRENT.md`
- `docs/memory/USER.md` (explicitly user-authorized durable collaboration
  preference)
- `docs/memory/CODEX.md` (corresponding actual-diff review workflow)
- `docs/memory/phases/e8-1-debug-1-fix1-preserve-diamond-cut.md`
- `docs/memory/manifest.json` (regenerated integrity-only index)

`docs/memory/phases/e8-1-debug-1-fix1-preserve-diamond-cut-task-contract.md`
is the pre-existing user-created contract and remains unmodified; it is part of
the eventual scoped commit and integrity-manifest inventory.
`docs/memory/phases/e8-1-debug-1-fix1-source-review-reconciliation-task-contract.md`
is also an intended new task record. The USER/CODEX changes are durable
workflow reconciliation only: a complete large review bundle is a clearly
temporary root-level file (for example, `kaonlt_review.diff`), includes tracked
diffs plus complete no-index additions for all intended new/untracked files,
requires no staging, and is removed before commit/push unless explicitly
tracked. They do not alter this repair's scientific ownership or validation
state. The original Debug.1 implementation and task records remain historical
and unchanged.

Next: user-controlled commit/push of the `SOURCE REVIEWED` repair, then the
targeted farm debug execution and return of fresh artifacts. E.8.1 remains
`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; its scientific first detailed
gate remains `Q4p4W2p74 / Left / highe` after the debug prerequisite is
resolved.
