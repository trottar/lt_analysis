---
memory_schema: 3
---
# Current KaonLT development state

## Active Objective

E.8 is `ACTIVE`: present the complete kaon missing-mass and signal-region
yield chain from authoritative objects and snapshots. Preserve baseline
production, frozen F.6.2 science, and detached Method-A/Method-B boundaries.

## Current Work Item

[E.8.4.Fix.5 / Fix.5.1 / Fix.5.2](phases/e8-4-fix5-shareable-method-a-impact-pages.md)
is `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`. The local candidate adds
shareable per-child authoritative SIMC overlays, stored yield-impact summaries,
and pinned parent closure. The existing single SIMC-yield call now precedes
E.8 finalization. The tracked Left/lowe owner performs run -> verify -> package;
its canonical-five profile uses explicit Left/lowe collection from a temporary
clean detached worktree at the exact pushed source. Independent actual-diff
review remains required before user commit/push and pushed-state synchronization.

## Verified State

F.4.Refresh.2 detached candidate materialization is `CLOSED / RUNTIME VALIDATED`
from accepted materialization evidence. F.6.3 current-baseline candidate-lineage
adoption and the existing E.8.4 Method-A impact audit are `CLOSED / RUNTIME VALIDATED` only for `Q4p4W2p74 / Left / lowe`, from the supplied independently
reviewed [runtime package](evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md).
The 87-page procedure manifest had no renderer failures; real child yields
changed while parent-t preservation closed. This does not validate the new
Fix.5 pages or their farm owner.

F.1-F.6.2 historical closures remain in the [roadmap](roadmap/STATUS.md).
E.8.1.Fix.5/Fix.6 retain narrow Left/lowe overlay/profile closure; canonical-five
expansion remains `DEFERRED`. E.8.2/E.8.3 retain their recorded `SOURCE REVIEWED`
status. Final canonical-five E.8 and F.6.4 remain `BLOCKED`. Method A stays
detached/non-production; Method B stays diagnostic-only and numerically excluded.

## Source / Evidence Identity

- Observed committed `test` HEAD for this local continuation:
  `da38444e7aa60efd62d6638780776344daf40276`, also the supplied Left/lowe
  farm-evaluated source; the new Fix.5 candidate is not farm-evaluated.
- Accepted Left/lowe ZIP SHA-256:
  `200fda66fe410274df1c9a8252b9e87114d8fd10a520fb7d71691ec6b3772874`.
  Candidate F.3/F.4 hashes and post-run model-output caveat are in the
  [evidence record](evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md).
- Accepted Refresh.2 ZIP SHA-256:
  `4cbdbdd2c0403b961614a8ccbb91104f13fef59194cd48e5e3570aded0954898`;
  farm HEAD `b349967c0d4210a78b144ce6134d3c1f15970245`, materializer source
  `141a3d04f9e5d07be21dba14e0e63212c3990bf1`. See the
  [adoption phase](phases/f6-3-current-baseline-candidate-lineage-adoption.md).
- Frozen F.6.2 JSON SHA-256:
  `5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1`;
  artifact `ee713b70de898bad8fa61164cbc8a54886af1ea1eb1e712ed95de17df8890cd0`;
  validation `7edc73fce20ad7dc7622b8c23a7ba7e8986e367ace605884e010945595370b3b`.
- Active background profile remains `no_empirical_residual`; legacy empirical
  residual Fit 1/Fit 2 remain dormant with both scales zero.

## Blockers

No deterministic local implementation blocker remains. Independent actual-diff
review, user commit/push, and pushed-state synchronization are required before
the tracked owner can run. New shareable-page rendering and owner integration
still require fresh farm evidence. Historical source/cache/lineage blockers are
preserved in supporting records; supplied current-baseline Left/lowe evidence
resolves the previous narrow runtime gate without promoting Method A.

## Next Action

NEXT — after user commit/push and pushed-state synchronization, run the tracked Q4p4W2p74 / Left / lowe E.8.4.Fix.5 farm owner to regenerate the full-analysis procedure PDF with the authoritative per-(t,phi) SIMC comparisons and shareable Method-A yield-impact pages, verify them, and return its single fresh validation ZIP.

## Success Criteria

The fresh Left/lowe PDF must retain existing E.8 pages and show all 27 canonical
children in Method-A/SIMC and baseline/Method-A/SIMC overlays, stored Y0/YA,
DeltaY and defined DeltaY/Y0 summaries, and parent-normalization sanity closure.
Require no renderer failures, verified artifact freshness/hashes and one fresh
ZIP from the tracked owner. Local tests do not prove ROOT/PyROOT or farm behavior.

## Do Not Reopen Without New Evidence

Do not reopen F.1-F.6.2 accepted science, change frozen authorities or baseline
production, make Method B numerical, independently normalize children, or
promote Method A. Narrow Left/lowe closure does not close canonical-five E.8
or F.6.4. E.8 consumes authoritative outputs; F.6.3 owns the private branch.

## Relevant References

- [E.8 procedure roadmap](decisions/e8-full-analysis-procedure-roadmap.md)
- [Roadmap status](roadmap/STATUS.md)
- [Fix.5 phase](phases/e8-4-fix5-shareable-method-a-impact-pages.md)
- [Left/lowe evidence](evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md)
- [F.6.3 adoption](phases/f6-3-current-baseline-candidate-lineage-adoption.md)
