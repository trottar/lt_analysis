# E.8.4 Fix.5 / Fix.5.1 / Fix.5.2 — shareable Method-A impact pages

**Status:** `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`.
Implemented/continued locally at committed `test` HEAD
`da38444e7aa60efd62d6638780776344daf40276`. Independent ChatGPT actual-diff
review, user commit/push and pushed-state synchronization remain required.

## Contracts and resolved blockers

[Fix.5](e8-4-fix5-shareable-method-a-impact-pages-task-contract.md) correctly
stopped before implementation because authoritative per-child SIMC MM was
produced only after the E.8 finalizer. [Fix.5.1](e8-4-fix5-1-simc-runtime-handoff-unblock-task-contract.md)
authorized moving the existing single SIMC-yield block before E.8 finalization.
The resulting local page candidate was preserved. Its owner/profile completion
correctly stopped at the frozen collector's canonical-five schema requirement.
[Fix.5.2](e8-4-fix5-2-collector-scope-owner-completion-task-contract.md) resolves
that operational blocker through canonical-five declaration plus the existing
explicit Left/lowe selection API; no collector change is needed.

## Source and presentation ownership

`src/main.py` changes only the order of the existing Step-6 SIMC-yield block:
existing data yields -> existing SIMC yields exactly once -> E.8 finalizer per
kaon setting -> existing correction ledger. The SIMC producer in
`src/binning/calculate_yield.py` and its normalization are unchanged.

`src/cuts/full_background_subtraction_plots.py` validates and clones the
already-normalized `hist["_xsect_support_simc"]["mm"]` matrix. It rejects missing,
nonfinite or geometrically incompatible children with explicit unavailable
status. There is no setting-wide fallback, ROOT-file reload, refill or scale.
The existing E.8.2/E.8.3/E.8.4 page IDs remain.

New `full_background.e8_4.*` pages:

- `method_a_vs_simc.t1/t2/t3`: final Method-A MM versus authoritative SIMC,
  all nine children per t, Lambda markers and explicit empty status.
- `baseline_method_a_simc.t1/t2/t3`: baseline/Method-A/SIMC overlays.
- `yield_summary.t1/t2/t3`: stored Y0/YA and statistical errors, signed DeltaY,
  and defined DeltaY/Y0; undefined fractions remain explicit.
- `parent_closure`: stored signed parent sums, residual and tolerance from the
  exact hash/fingerprint-pinned current candidate F.4, with a displayed relative
  residual derived only from those stored numbers. Equality is labeled as the
  parent-normalization sanity constraint.

No production objects, event loops, yield extraction, Method-A mathematics,
Method B, cuts, binning, efficiencies or accepted authority change.

## Complete tracked farm owner

`testing/run_e8_4_fix5_left_lowe_plot_gate.py` binds only the effective profile's
required analysis SHA to the supplied exact reviewed/pushed source. It requires
branch test, HEAD/local origin/test identity, ancestry and gate-relevant source
provenance before running the unchanged debug launcher. Known model/output
artifacts are recorded, excluded from the ordinary gate diff check, never cleaned.

The owner verifies candidate F.3/F.4 hashes, analysis completion, freshness,
existing critical pages, all ten new IDs, all nine children on each new t page,
and zero renderer failures; writes provenance; then creates a temporary clean
detached worktree at the exact pushed SHA. It loads the unchanged collector
from that tree with its source checks rooted there, artifacts still from the
ordinary canonical OUTDIR, and explicit `phi="Left", epsilon="lowe"`.
The tracked profile declares the canonical five, with unchanged required
artifact inventory and base placeholder. The effective copy is temporary only.

The ZIP must contain only one Left/lowe setting and all eight required artifacts,
with exact source identity and hashes. The temporary tree never runs analysis;
bounded removal also prunes only its own Git worktree registration, including
ignored bytecode from source checks. Failure cleans up that temporary tree and
prints no ZIP path. Success prints exactly one POSIX ZIP path. There is no
separate manual packaging step.

## Deterministic checks and runtime boundary

Python 3.12.10: nine-file `py_compile` passed using temporary bytecode paths.
Main-order 2, E.8.4 19, plotting 102 (18 skipped), F.6.3 28, debug launcher 9,
profile 2, owner 9, generic collector 28, E.8.2 23 and E.8.3 21 tests passed.
Total: 243 run, 225 passed, 18 skipped. Four skips require unavailable PyROOT;
fourteen cover explicitly retired cumulative presentation contracts. Exact
names/reasons and outputs are included in the complete review bundle.

Tests include real generic collection on synthetic artifacts, one-setting
inventory, effective-profile immutability, clean detached source identity,
bounded cleanup on success/failure, source preflight, stale artifacts,
page failures, fresh-ZIP refusal and deterministic stdout. No farm, local
`main.py`, commit or push was performed. New rendered-page/owner runtime
acceptance awaits fresh farm evidence.

## Warranted milestone alignment

The supplied independently reviewed [Left/lowe runtime evidence](../evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md)
closes the previous narrow F.6.3/current-baseline and existing E.8.4 gate.
It does not validate Fix.5, close canonical-five E.8/F.6.4, or promote Method A.
CURRENT owns the sole authoritative action:

NEXT — after user commit/push and pushed-state synchronization, run the tracked Q4p4W2p74 / Left / lowe E.8.4.Fix.5 farm owner to regenerate the full-analysis procedure PDF with the authoritative per-(t,phi) SIMC comparisons and shareable Method-A yield-impact pages, verify them, and return its single fresh validation ZIP.
