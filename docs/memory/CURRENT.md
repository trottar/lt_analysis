---
memory_schema: 3
---
# Current KaonLT development state

## Active Objective

E.8 is `ACTIVE`: audit the complete kaon missing-mass and signal-region yield
chain from authoritative objects. Preserve baseline production, accepted
upstream authorities, and detached Method-A/Method-B boundaries.

## Current Work Item

[E.8.4 Fix.5.6](phases/e8-4-fix5-6-owner-farm-readiness-and-failure-provenance.md)
is `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`: the local operational owner
candidate runs existing collector/source checks in an exact clean detached
worktree before analysis and atomically persists per-attempt stage/failure
status. Existing post-analysis artifact/page/collection/ZIP checks and bundle
inventory are preserved. Owner runtime and final delivery await a farm attempt.

Fix.5.5 source is `SOURCE REVIEWED` at
`df957a6414fc9c515d1f82228517cb801dc90350`; independent ChatGPT actual-diff
and pushed-state review passed per the supplied Fix.5.6 contract. Fix.5.4 and
Fix.5.5 remain `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`: source acceptance
establishes no new numerical farm closure or actual PDF visual quality.

## Verified State

The supplied Fix.5 observation at `Q4p4W2p74 / Left / lowe` reports a fresh
97-page PDF, zero renderer failures, all ten new Fix.5 pages and prior E.8.4
pages retained. The owner returned no final ZIP: this is render evidence,
not accepted bundle closure. See the [post-farm investigation](investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md).

The E.8.3 historical accepted F.6.1 aggregate and the current F.6.3/E.8.4
candidate branch have a confirmed lineage mismatch. Historical E.8.3 remains
`SOURCE REVIEWED`; it cannot explain the current branch without exact lineage
identity. Baseline blue curves obscured by Method-A magenta are a confirmed
presentation deficiency; Fix.5.5 addresses it from the reviewed numerical source,
with real visual acceptance still pending.

F.4.Refresh.2 remains `CLOSED / RUNTIME VALIDATED`. F.6.3 and the existing
E.8.4 gate remain `CLOSED / RUNTIME VALIDATED` only for Q4p4W2p74 / Left / lowe:
branch execution,
live-cache parity, real child changes and signed parent preservation, from the
[earlier evidence](evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md).
F.1-F.6.2 historical closures remain in the [roadmap](roadmap/STATUS.md).
E.8.1 canonical-five expansion stays `DEFERRED`; E.8.2/E.8.3 stay
`SOURCE REVIEWED`. Method A is detached/non-production; Method B is
diagnostic-only and numerically absent. Final canonical-five E.8 and F.6.4
remain `BLOCKED`.

## Source / Evidence Identity

- Supplied prior Fix.5 committed/farm source:
  `f9d70732290ea461096374ca1270b47452644991`.
- Fix.5.6 starting/unchanged committed `test` HEAD on 2026-10-02:
  `df957a6414fc9c515d1f82228517cb801dc90350`, parent
  `761fbb6c03d2d7a10bb911cf84e9ba898496fab6`. The supplied contract records
  independent Fix.5.5 actual-diff/pushed-state review PASS with exactly the
  reviewed 12 paths. Fix.5.6 owner/test/memory edits are local review candidates.
- Fresh PDF/manifest: farm-visible timestamp `2026-10-01 23:34`, PDF about
  2.2 MB; inspection copies `KaonLT_E8_4_Fix5_Left_lowe_20261001-234439.pdf`
  and `KaonLT_E8_4_Fix5_Left_lowe_20261001-234439-manifest.json`.
  These supplied observations were not independently reopened by Codex here.
- Historical E.8.3 F.4 raw SHA begins `adcc0190`; current F.6.3 candidate F.4
  raw SHA is `1d545924eba89c7f9ffa28028e307aca9b434a89beec06863cf2893887b6b902`.
  Full historical hashes and audit requirements are in the investigation.
- Earlier accepted Left/lowe farm source:
  `da38444e7aa60efd62d6638780776344daf40276`; ZIP SHA-256
  `200fda66fe410274df1c9a8252b9e87114d8fd10a520fb7d71691ec6b3772874`.
  Its candidate identities and post-run model-output caveat remain unchanged.
- Accepted Refresh.2 materialization and frozen F.6.2 identities remain
  unchanged in their linked canonical records. The active background profile
  remains `no_empirical_residual`; legacy empirical residual scales stay zero.

## Blockers

Scientific interpretation of fresh Fix.5 pages remains `BLOCKED` pending farm
evidence for the new identities. Local synthetic tests establish binwise
algebra, normal-bin yield closure, signed support and current-lineage aggregate
validation; they do not explain the observed farm curves or establish signed
cancellation there. SIMC uses `iter_weight * normfac / Ncontribute`; the traced
source does not define its luminosity/effective-charge units. E.8.4 therefore
retains an available current-F.6.3 payload and structured audits but marks only
absolute-SIMC comparison pages/claims unavailable with literal provenance.
No wrong normalization or permitted conversion is established.

Fix.5.5 real visual validation remains pending. The prior missing owner ZIP
has no established exact failure cause. The old analysis log captures only the
launcher child stream, so absence of an owner failure there cannot exclude a
post-analysis failure. Fix.5.6 adds durable owner-stage provenance and early
source checks; it does not diagnose which prior stage failed or authorize a run.

## Next Action

NEXT — Fix.5.6 actual-diff review -> user commit/push -> pushed-state review ->
one narrow Q4p4W2p74 / Left / lowe tracked-owner farm gate -> fresh
ZIP/status/log/PDF/manifest evidence review. Preserve the SIMC absolute-unit
blocker; no farm command or execution is authorized in this implementation task.

## Success Criteria

Local Fix.5.6 checks: 66 tests passed, no skips; syntax and diff checks pass.
The scoped 199-file candidate manifest and ordinary memory health pass with no
warnings or hard failures.
Source and deterministic tests do not establish farm filesystem/permissions,
ROOT/PyROOT, full analysis, ZIP delivery, PDF legibility, numerical closure or
observed signed cancellation. All scientific/presentation source is frozen.

## Do Not Reopen Without New Evidence

Do not reopen historical closures beyond the affected interpretation, replace
accepted authorities, change physics or normalization for visibility, make
Method B numerical, independently normalize children, or promote Method A.
E.8 remains a consumer; F.6.3 owns the private branch. Narrow Left/lowe evidence
does not close canonical-five E.8 or F.6.4.

## Relevant References

- [Fix.5.6 phase](phases/e8-4-fix5-6-owner-farm-readiness-and-failure-provenance.md)
- [Fix.5.5 phase](phases/e8-4-fix5-5-current-lineage-visualization-clarity.md)
- [Fix.5.4 phase](phases/e8-4-fix5-4-current-lineage-identity-audit.md)
- [Post-farm investigation](investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md)
- [Fix.5 chronology](phases/e8-4-fix5-shareable-method-a-impact-pages.md)
- [Earlier Left/lowe evidence](evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md)
- [Roadmap status](roadmap/STATUS.md)
- [E.8 procedure roadmap](decisions/e8-full-analysis-procedure-roadmap.md)
