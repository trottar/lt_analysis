# E.8.1.Debug.1.Fix.1 — preserve ordinary Diamond cut in Left/lowe debug mode

## Starting state

Branch: `test`

Exact starting HEAD:

`4c4a18a098d184fe0b13dd1fc57010ac7f5b30c0`

That commit introduced the requested `run_Prod_Analysis.sh -d` Left/lowe
debug mode.

Independent ChatGPT review found the launcher/preflight/canonical-binning
orchestration correct, but found one source-level blocker in `src/main.py`.

This Fix.1 repairs only that blocker.

The task-contract file itself may be the sole user-created uncommitted file
when Codex starts. Any other unexpected pre-existing worktree change is a
hard stop.

## Source-review finding

At the starting HEAD, `LT_ANALYSIS_DEBUG_LEFT_LOW` changes:

`phisetlist = ["Center", "Left", "Right"]`

to:

`phisetlist = ["Left"]`

before Step 2 calls `DiamondPlot()`.

That changes scientific cut construction.

The ordinary main.py path calls `DiamondPlot()` for Center first.
`DiamondPlot()`'s Center path builds the common Diamond polygon from the
available Center/Left/Right and low/high sources and returns that common
polygon through `cut_mode="poly"` / `poly_points` into `inpDict`.

The subsequent ordinary Left call therefore consumes the already-established
common polygon.

When debug mode suppresses Center before Step 2, the Left call has no input
polygon. `DiamondPlot()` then takes its direct-call fallback path and derives a
Left-low-only polygon (`cut_poly_source="low_fallback"`).

Therefore the pushed debug mode changes the Diamond cut and violates the
original contract's frozen-cut requirement.

## Objective

Preserve the ordinary Diamond-cut construction exactly, while retaining the
requested fast debug behavior downstream:

1. paired low/high canonical preflight remains unchanged;
2. the full low invocation enters debug mode;
3. ordinary shared Diamond-cut preparation occurs before phi restriction;
4. after Diamond preparation, downstream full analysis proceeds with only
   `Left`;
5. successful Left/lowe completion stops the launcher before full high epsilon;
6. the non-debug path remains unchanged.

## Scientific ownership

This is runtime-selection repair only.

Do not change:

- Diamond algorithms or fitted/common polygon definitions;
- cuts or acceptance;
- canonical t/phi binning;
- random subtraction;
- slow-proton subtraction;
- pion subtraction;
- HGCer Method A or Method B;
- SIMC normalization/templates;
- background optimization physics;
- yields;
- efficiencies;
- L/T separation;
- uncertainty propagation;
- production correction logic.

The debug path must consume the same ordinary common Diamond polygon that an
ordinary Left/lowe run would consume.

## Required implementation

### `src/main.py`

Keep early validation of `LT_ANALYSIS_DEBUG_LEFT_LOW`:

- explicit enabled values only;
- kaon only;
- low epsilon only;
- reject canonical-capture misuse.

But do **not** narrow `phisetlist` before Step 2.

The existing ordinary Step-2 Diamond path must remain intact and execute with:

`["Center", "Left", "Right"]`

in its existing order.

Do not modify `src/cuts/diamond.py`.

Do not optimize or shorten the Diamond stage in this fix. Scientific/path
identity is more important than saving the small additional Step-2 cost.

Only after the existing Diamond calls and their ordinary artifact registration
have completed, and before `shift_prep` begins, apply the debug restriction:

`phisetlist = ["Left"]`

for the downstream analysis.

Thus the runtime path must be:

ordinary paired low/high raw canonical preflight
→ full low debug main.py
→ ordinary Center/Left/Right Diamond preparation
→ common Center-produced polygon established in inpDict
→ debug restriction to Left
→ Left-only shift/SIMC/subtraction/binning/yield downstream processing
→ main.py completes successfully
→ launcher exits before full high epsilon.

With the debug environment absent, there must be no new setting restriction
and no behavior change.

### `run_Prod_Analysis.sh`

The pushed launcher implementation is correct for this finding and is frozen.

Do not modify it unless source tracing reveals a concrete blocker. None is
currently known.

### Canonical binning

Do not alter canonical preflight or resolver source.

The existing paired preflight and validated-authoritative interval reuse must
remain unchanged.

## Allowed source/test files

Allowed implementation files:

- `src/main.py`
- `testing/test_run_prod_analysis_debug_left_low.py`

Required memory/history files:

- `docs/memory/CURRENT.md`
- NEW `docs/memory/phases/e8-1-debug-1-fix1-preserve-diamond-cut.md`
- regenerated `docs/memory/manifest.json`

This contract file is also part of the eventual scoped commit.

Everything else is frozen, including:

- `run_Prod_Analysis.sh`
- `src/cuts/diamond.py`
- `src/binning/find_bins.py`
- `src/utility/shared_canonical_binning_preflight.py`
- subtraction/background/SIMC/yield source
- E.8.1 profiles/collectors/wrapper
- existing historical Debug.1 phase/task records.

Do not rewrite the original implementation record merely to erase the
historical failed source-review state. Fix.1 owns the correction chronology.

## Required regression tests

Extend the focused debug test so it proves source ordering, not merely the
presence of `phisetlist = ["Left"]`.

At minimum verify:

1. the ordinary default is still:
   `phisetlist = ["Center", "Left", "Right"]`;
2. early debug validation does not narrow `phisetlist`;
3. the Step-2 `DiamondPlot()` loop occurs while the ordinary setting list is
   still active;
4. the debug `phisetlist = ["Left"]` restriction occurs only after Step-2
   Diamond preparation/artifact registration;
5. the restriction occurs before the `shift_prep` setting loop;
6. no change was made to `diamond.py`;
7. canonical preflight remains paired and unrestricted;
8. only full low receives the debug environment control;
9. successful debug mode still stops before full high;
10. failed full low still propagates failure;
11. ordinary non-debug full low/high behavior remains present.

Keep tests non-destructive. Do not execute the production launcher or
`git clean -fdx`.

Run:

- `python -m py_compile src/main.py testing/test_run_prod_analysis_debug_left_low.py`
- `bash -n run_Prod_Analysis.sh` using the available local Bash
- focused debug tests
- applicable memory-health/manifest/bootstrap tests
- `git diff --check`

## Memory updates

`CURRENT.md` must truthfully record:

- pushed Debug.1 source `4c4a18a098d184fe0b13dd1fc57010ac7f5b30c0`
  received independent source review;
- that review found the pre-Step-2 selector/Diamond-cut blocker;
- E.8.1.Debug.1.Fix.1 is the immediate ACTIVE repair;
- E.8.1 itself remains DEVELOPMENT COMPLETE, FARM VALIDATION PENDING;
- the scientific E.8.1 first detailed gate remains `Q4p4W2p74 / Left / highe`
  after the debug prerequisite is resolved.

The new Fix.1 phase record must preserve:

- starting HEAD;
- exact source-review finding;
- ordinary Center-produced common-cut dependency;
- repair implementation;
- tests actually run;
- changed-file accounting;
- local-vs-farm boundary;
- remaining ChatGPT review/user push/farm-debug steps.

Do not claim SOURCE REVIEWED for Fix.1 until ChatGPT audits the actual diff.
Do not claim farm/runtime validation.

Regenerate `docs/memory/manifest.json`.

Update other durable memory only if its owned knowledge genuinely changes.

## Acceptance criteria

The repair is locally acceptable only if:

- debug and ordinary Left/lowe use the same ordinary common Diamond cut path;
- Center Diamond preparation is not skipped before debug restriction;
- downstream expensive processing is Left-only in debug mode;
- paired canonical preflight remains untouched;
- launcher debug stop behavior remains untouched;
- no scientific calculation or production ownership changes;
- ordinary non-debug behavior remains unchanged;
- focused deterministic tests pass;
- memory/history truthfully records the failed review and Fix.1 state.

## Farm boundary

No local test proves ROOT/PyROOT or full `main.py` execution.

After ChatGPT source review and user push, the user will run the targeted
Jefferson Lab debug execution and return fresh evidence.

## Hard stop

After implementation, local checks, memory updates, and actual-diff audit:

STOP.

Do not commit.
Do not push.
Do not run the Jefferson Lab farm.

Provide:

- concise implementation summary;
- exact changed-file list;
- tests/checks and results;
- explicit source trace showing Diamond preparation before Left restriction;
- exact non-executed scoped user Git handoff commands.
