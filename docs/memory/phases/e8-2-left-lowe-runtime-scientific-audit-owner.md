# E.8.2 Left/lowe runtime scientific-audit owner

Status: `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`.

## Authority and scope

Implemented the [revised task contract](e8-2-left-lowe-runtime-scientific-audit-owner-task-contract.md)
at branch `test`, starting HEAD/local `origin/test`
`7b8eb3cb231d289d18de9c168960ddddcaa39254`. This is validation infrastructure
for the active E.8.2 science audit. Production/scientific source, the collector,
existing owners/profiles, accepted evidence and scientific ownership are frozen.
The supplied contract supersedes its earlier unimplemented one-setting version.

## Source audit

SOURCE VERIFIED: the existing baseline follows the authoritative traversal,
random/dummy subtraction, event-by-event slow-proton cleaning, unchanged
production `prune_hist`, baseline pion subtraction and producer-owned final
MM0/Y0. Same-traversal detached pre-/post-proton bookkeeping observes the
conceptual proton stage without changing production ordering. The E.8.2
producer/reader checks stage algebra and the post-prune pion-input handoff;
the renderer consumes producer-owned yields and uncertainties. The active
profile remains `no_empirical_residual`, with zero empirical residual scales.
Relevant source: `src/cuts/full_background_subtraction_plots.py`,
`src/utility/utility.py`, `src/utility/correction_ledger.py` and
`src/utility/background_config.py`; the unchanged E.8.2 regression checks
same-traversal capture, pruning observations and finalization boundaries.

No new E.8.2 scientific implementation is justified by this source audit.
Actual stage magnitudes and pruning impact require a narrow runtime gate.

## Owner mechanics

SOURCE VERIFIED: `testing/run_e8_2_left_lowe_scientific_audit_gate.py` reuses
the established ordinary-checkout snapshot, disposable detached worktree,
copied-ltsep overlay, path probes, external-path preflight, bounded cleanup,
preservation, detached collection and companion-delivery helpers. It owns one
debug analysis invocation in the isolated analysis worktree, checks the debug
completion markers, and rejects high-epsilon completion. It does not invoke
candidate materialization/staging, F.6.3 lineage or E.8.4 page acceptance.

The new v4 profile declares canonical five in the frozen collector order.
This authorization inventory is distinct from the owner's explicit
`phi="Left", epsilon="lowe"` request. The unchanged collector packages only
Left/lowe plus the global summary. A one-setting v4 declaration remains invalid.
The template pins the starting HEAD; the owner resolves it to the exact
reviewed/pushed full source SHA with no allowed committed source exceptions.

The owner verifies all 18 E.8.2 pages, parent/stage identities, ordered nine
phi children, presentation flags, empty invalid-child/renderer-failure lists
and terminal handoff. Extra producer setting metadata is tolerated.
E.8.3/E.8.4 unavailable pages alone do not fail E.8.2. Freshness, nonzero size,
strict JSON, PDF magic, baseline full-analysis/ledger identity and exact ZIP
source/inventory/hash/size checks fail closed. Atomic per-attempt status and
the summary retain isolation, cleanup, preservation and acceptance boundaries;
the ZIP is accompanied by the analysis log, summary and final gate receipt.

## Deterministic local validation

Python 3.12.10; no farm, ROOT/PyROOT or full-analysis execution:

- `python -B -m py_compile` for the new owner and test: PASS; cache directed
  to a temporary directory outside the repository.
- `python -B -m unittest testing.test_run_e8_2_left_lowe_scientific_audit_gate -v`:
  PASS, 23 tests with rejection subtests and mocked operational helpers.
- `python -B -m unittest testing.test_e8_2_baseline_stage_audit -v`:
  PASS, 23 tests; existing test/source unchanged.
- Collector and canonical-five profile regression modules: PASS, 31 tests.

The real generic collector is exercised on fixture artifacts with source
commands mocked: five declared settings yield one requested/packaged Left/lowe
setting, exactly five setting artifacts and one global summary. Owner-flow
fixtures cover one invocation, copied-overlay execution, child/source failure,
ordinary/installed-ltsep preservation failure, cleanup failure, filtering,
ZIP verification and companion delivery. These are local infrastructure checks.

## Evidence and next gate

RUNTIME VERIFIED: none from this implementation task.

NOT VERIFIED: actual Left/lowe stage magnitudes, pruning impact, PDF legibility
and farm integration. E.8.2 runtime/scientific acceptance awaits supplied farm
evidence and independent review. Method A remains detached/non-production and
is not required by this gate; Method B remains diagnostic/cross-check only and
numerically excluded. No production promotion, E.8.3/E.8.4/canonical-five
acceptance or absolute-SIMC amplitude claim follows from this implementation.
The canonical-five full runtime blocker and deferred provenance repair remain;
final E.8/F.6.4 and the absolute-SIMC units blocker remain unchanged.

CURRENT owns the sole ordinary next action: actual-diff review, user-controlled
commit/push and pushed-state synchronization/farm-readiness review precede one
isolated Left/lowe E.8.2 audit and scientific review of its ZIP/companions.
Codex stops at the complete review bundle; it issues no farm command.
