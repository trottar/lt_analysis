# E.8.1.Debug.1 — `run_Prod_Analysis.sh` Left/lowe debug-mode task contract

## Status and historical context

This is a user-requested operational prerequisite before the existing E.8.1
farm-validation gate.

Starting branch:

`test`

Exact starting HEAD:

`16d595edabd2a348f51b728bafe9731e17b430c0`

At this HEAD:

- M0–M9 repository-memory migration is SOURCE REVIEWED.
- F.6.2 and F.6.2.Fix.5 are CLOSED / RUNTIME VALIDATED.
- E.8.1 is DEVELOPMENT COMPLETE, FARM VALIDATION PENDING.
- the reviewed bundle wrapper/state reconciliation is SOURCE REVIEWED.
- the scientific E.8.1 NEXT remains the canonical-five-setting
  `Q4p4W2p74` farm gate with `Left / highe` inspected first.
- before that gate, the user requested a narrow launcher debug mode so one
  full `Left / lowe` analysis can be run quickly while developing/debugging.

This task does not reopen or alter E.8.1, F.6.2, Method A, Method B, pion
subtraction, proton subtraction, canonical binning science, yields, or final
cross sections.

The task-contract file itself may be the sole pre-existing uncommitted file
when Codex starts. Any other unexpected worktree change is a hard stop.

## Objective

Add a narrow `-d` debug mode to tracked:

`run_Prod_Analysis.sh`

so the production launcher can perform the following diagnostic execution:

1. preserve the existing paired low/high raw canonical-binning preflight;
2. run the actual full low-epsilon analysis with only the `Left` phi setting;
3. allow that `Left / lowe` `main.py` invocation to complete normally;
4. stop successfully immediately afterward;
5. do not start the actual full high-epsilon analysis.

Without `-d`, behavior must remain unchanged.

The intended user-facing invocation is the existing launcher syntax with an
additional short flag, for example:

`./run_Prod_Analysis.sh -d 4p4 2p74`

Preserve the launcher's existing short-option/positional-argument convention;
do not redesign its CLI as part of this task.

## Scientific and runtime ownership

This is runtime/orchestration debug control only.

The new mode must not modify:

- scientific cuts;
- normalization;
- random subtraction;
- slow-proton subtraction;
- pion-background treatment;
- HGCer Method A or Method B;
- SIMC templates or normalization;
- fit windows, priors, or component definitions;
- canonical t/phi definitions or support policy;
- yield calculation;
- acceptance/binning science;
- efficiency corrections;
- L/T separation;
- uncertainties;
- production correction logic.

Debug selection determines only which detector setting proceeds through the
full low-epsilon runtime.

## Audited current runtime path

The current launcher builds separate low- and high-epsilon `main.py` argument
arrays.

For kaon analysis it then performs, before either full analysis:

- low raw canonical-prepass capture using
  `LT_ANALYSIS_CANONICAL_PREPASS_CAPTURE`;
- high raw canonical-prepass capture using
  `LT_ANALYSIS_CANONICAL_PREPASS_CAPTURE`;
- `utility/shared_canonical_binning_preflight.py`.

Those capture-mode `main.py` calls intentionally exit before full
SIMC/background processing.

The ordinary launcher then performs the actual low and high full analyses.

Inside `src/main.py`, the ordinary full-analysis setting population begins as:

`["Center", "Left", "Right"]`

and the analysis proceeds stage-major over that population.

Therefore a launcher-only early exit after the low invocation would still run
all available low-epsilon phi settings. The debug mode requires a narrowly
scoped full-analysis setting selector inside `src/main.py`.

## Required implementation

### 1. `run_Prod_Analysis.sh`

Add short option `-d`.

The help text must define it unambiguously as the fast `Left / lowe` debug
mode.

Use a dedicated shell state variable such as `d_flag`; do not overload the
existing `DEBUG` plotting/splash variable.

Ensure the existing Q2/W positional-resolution logic recognizes `-d` as a
launcher flag in the same manner as the current launcher flags. Do not
otherwise refactor argument parsing.

Debug mode is intended for kaon full-analysis debugging. After particle type
has been resolved, fail explicitly if debug mode would execute a non-kaon
analysis rather than silently changing semantics.

Do not apply the new full-analysis setting selector to either canonical
prepass capture.

The paired preflight must remain:

low raw capture
→ high raw capture
→ shared canonical preflight

with the same ordinary setting/input population that would be used without
`-d`.

Only the actual full low-epsilon `main.py` invocation receives the dedicated
debug-setting control.

Prefer an invocation-scoped environment variable, for example:

`LT_ANALYSIS_DEBUG_LEFT_LOW=1 python3 main.py ...`

Do not `export` that variable globally where it could contaminate canonical
capture or later commands.

The existing nonzero failure behavior of `main.py` must be preserved. A
failed Left/lowe run must return failure; debug mode must not convert failure
to success.

After the full low-epsilon invocation returns successfully in debug mode:

- print a clear message that the requested `Left / lowe` debug analysis
  completed;
- state that full high-epsilon processing is intentionally skipped;
- exit successfully before the full high-epsilon invocation.

The high-epsilon raw canonical-prepass capture must already have occurred.

Do not change the ordinary non-debug low/high loop behavior.

Do not change cleanup/symlink behavior except where absolutely required for
the flag to follow the existing launcher path.

### 2. `src/main.py`

Add only the minimal runtime selector needed by the launcher.

The ordinary default must remain:

`phisetlist = ["Center", "Left", "Right"]`

When the dedicated Left/lowe debug environment control is enabled for a kaon
low-epsilon full analysis, restrict the full-analysis population to exactly:

`["Left"]`

Do this before detector-setting-dependent full-analysis stage loops.

Do not introduce a general arbitrary-phi production selector.

The debug selector must not redefine canonical t/phi intervals.

The full low-epsilon run must continue to consume the canonical interval
state produced by the existing paired low/high preflight. If source tracing
shows that singleton Left execution would instead silently regenerate or
replace the paired canonical result, stop and report that as a blocker rather
than inventing a fallback.

If the debug environment control is encountered in an unsupported context,
fail explicitly rather than silently broadening or changing the requested
behavior.

### 3. Ordinary path preservation

With `-d` absent:

- preflight commands remain unchanged;
- `phisetlist` remains ordinary;
- full low analysis remains ordinary;
- full high analysis remains ordinary;
- no debug environment variable affects runtime;
- no new early exit occurs.

## Allowed source/test files

Source changes are restricted to:

- `run_Prod_Analysis.sh`
- `src/main.py`
- NEW `testing/test_run_prod_analysis_debug_left_low.py`

If an existing directly equivalent dedicated test file is found, Codex may
extend that file instead of creating the named new test, but must explain the
choice and must not broaden testing cleanup/refactors.

Everything else in scientific/runtime source is frozen.

In particular, do not modify:

- canonical-binning algorithms/helpers;
- `src/utility/shared_canonical_binning_preflight.py`;
- background configuration;
- subtraction modules;
- SIMC modules;
- yield/binning algorithms;
- presentation/rendering source;
- E.8.1 collector/profile/wrapper source.

## Required deterministic tests

Tests must not execute a destructive real production analysis or trigger
`git clean -fdx`.

Provide deterministic coverage establishing at minimum:

- shell syntax remains valid (`bash -n run_Prod_Analysis.sh`);
- `-d` is accepted and documented;
- Q2/W parsing follows the existing flagged invocation convention;
- debug mode does not affect either low or high canonical-prepass capture;
- paired low/high canonical preflight remains before full low execution;
- only the full low invocation receives the debug-setting control;
- successful debug low execution stops before full high execution;
- failed debug low execution remains nonzero;
- ordinary non-debug execution retains both full low and full high paths;
- `src/main.py` defaults to the ordinary phi-setting population;
- the dedicated debug control selects exactly `Left` for supported low-kaon
  execution;
- the debug selector does not alter canonical-binning configuration.

Tests may use source/AST/shell-contract inspection where importing or running
the real ROOT-dependent analysis would violate the local/farm boundary.

Run `python -m py_compile` on any changed/new Python file for which compilation
is locally safe.

Do not represent static/source-contract tests as proof of farm execution.

## Forbidden shortcuts

Do not:

- comment out Center/Right code;
- permanently change `phisetlist` to `["Left"]`;
- remove the high canonical capture;
- derive debug-only Left canonical bins;
- skip shared canonical preflight;
- fake high-epsilon state;
- alter physics to make singleton execution succeed;
- suppress an exception/failure to reach the debug stop;
- restructure `main.py` into a setting-major pipeline;
- introduce a broad arbitrary setting-selection feature;
- change unrelated launcher behavior;
- perform cleanup/refactors outside the narrow implementation;
- run the Jefferson Lab farm;
- commit or push.

If singleton Left/lowe cannot safely traverse the existing low-epsilon full
runtime without changing scientific ownership or canonical binning, stop and
report the exact blocker.

## Memory/history updates

This standalone contract is not a replacement for durable project memory.

As part of the scoped implementation, update only warranted memory.

Required:

- `docs/memory/CURRENT.md`
- NEW
  `docs/memory/phases/e8-1-debug-1-left-low-debug-mode.md`
- `docs/memory/manifest.json`

`CURRENT.md` must record this user-requested debug task as the immediate
operational prerequisite while preserving the scientific E.8.1 state and
scientific NEXT behind it.

The permanent phase record must include:

- starting HEAD;
- historical reason for the task;
- exact implementation scope;
- preserved paired canonical preflight;
- ordinary-path preservation;
- tests actually run and results;
- explicit local-vs-farm boundary;
- actual changed-file list;
- remaining review/push/farm steps.

After local implementation, do not call the change SOURCE REVIEWED: ChatGPT
must inspect the actual diff first. Do not claim runtime validation.

Update `MEMORY.md` or `LEARNINGS.md` only if implementation establishes new
durable knowledge that belongs there. Do not edit them merely to duplicate
the phase record.

Do not modify `USER.md`, `CODEX.md`, `TOOLS.md`, or
`handoffs/CURRENT_HANDOFF.md` unless a concrete blocker proves their owned
durable information changed. None is expected from this task.

Regenerate `docs/memory/manifest.json` whenever versioned memory changes.

## Local validation and diff audit

Before finishing, Codex must:

- inspect the actual final diff;
- confirm all changed paths are allowlisted or explicitly warranted memory;
- run `git diff --check`;
- run the focused deterministic tests;
- run applicable memory-health/manifest/bootstrap checks after memory edits;
- report each command and result;
- report any tests not runnable locally;
- trace the launcher path from `-d` through paired preflight, singleton full
  low execution, and intentional stop;
- separately trace the ordinary non-debug path and confirm it remains intact.

The real ROOT/PyROOT/full-analysis runtime remains farm-only.

## Acceptance criteria

Local implementation is acceptable for ChatGPT review only if:

- `-d` has the requested narrow meaning;
- paired low/high canonical preflight is preserved;
- canonical prepass populations are not debug-restricted;
- the actual low full analysis is restricted to Left only;
- successful debug execution stops before actual high full analysis;
- failure remains failure;
- ordinary behavior is unchanged without `-d`;
- no scientific algorithm or ownership changes;
- deterministic local tests pass;
- durable memory/history truthfully records local/proposed state;
- diff is narrow and auditable.

Farm behavior remains unvalidated until the user runs it at Jefferson Lab.

## Hard stop

After implementation, tests, memory updates, and actual-diff audit:

STOP.

Do not commit.
Do not push.
Do not run the farm.

Provide:

- concise implementation summary;
- exact changed-file list;
- tests/checks with results;
- runtime-path trace;
- any remaining limitation;
- exact non-executed user Git handoff commands using scoped `git add`,
  `git commit`, and `git push origin test`.
