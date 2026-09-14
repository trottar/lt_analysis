# Current KaonLT handoff

Prepared: 2026-09-14 after accepted F.1.Fix.5 closure. This record
distinguishes the analysis runtime, collector reconciliation, and F.2 baseline.

## Resume identity

- Repository: trottar/lt_analysis; branch: test.
- Recheck git status --short --branch, git rev-parse HEAD, and git rev-parse
  origin/test before work.
- Reviewed F.1.Fix.5 analysis-source commit:
  dc4fc6283001739a487ec80068f951b0e388cae6.
- F.1 runtime-evaluated analysis HEAD:
  126fa22c19bd29b9952f55b33ab43d59f9727ef6.
- F.1.Validation.Fix.1 collector reconciliation commit:
  31dd034d8404e317863bf0933c51253e1d3deeb8.
- This docs/memory-only closure commit is the exact F.2 baseline. It does not
  alter or replace the F.1 runtime identity.
- F.2 source implementation started at and is pinned to
  `02f21886cff8a721df46ebf673594856583fc3c5`.

## Established evidence

CLOSED / RUNTIME VALIDATED:

- Phase-C five-setting Method-B diagnostic closure: C.Fix.2.3
  9a66bc62d20a99172e326e915866877b65ae1e5d, later accepted pre-E.3
  e3853655db0809923cbf2326e2f779219128eda9; adaptive B is DO NOT PROMOTE.
  Read evidence/phase-c-five-setting-closure.md.
- E.3.Fix.2 Q4p4W2p74 Left-low presentation: implementation
  eb1710f4739ba6ef14f51419806e9fc5bd53c175, runtime/bundle
  bf53dac84e1396cfd7e3f4e0234749426bfbdcf4. Read
  evidence/e3-fix2-left-low-runtime.md. Do not extrapolate it to four other
  settings.

CLOSED / RUNTIME VALIDATED:

- F.1.Fix.5 for all five canonical Q4p4W2p74 settings. The runtime remains
  anchored to `126fa22c19bd29b9952f55b33ab43d59f9727ef6`; C1 reconciled its
  clean collector provenance without rerunning KaonLT. The owner accepted the
  complete gate and artifact-hash continuity with the supplied v4 bundle.

## Current gate and boundaries

CLOSED / RUNTIME VALIDATED — F.1 is accepted. C1 limits only its committed
range whitespace check to the three profile-owned validation files; the global
worktree check and committed-file identity audit remain unchanged. The F.1
artifact, PDF, and page-manifest evidence remains detached and source-owned.

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING — F.2 adds only a detached
representation-audit module, CLI, and focused tests. It uses fixed candidate
bases; strict F.1 v2 validation; canonical-t-local deterministic response
probes; and cKDTree non-prompt application support. It creates no map,
correction, event score/probability, weight, normalization, or production
object. Read phases/phase-f2-method-a-acceptance-representation.md.

Keep random, slow-proton, pion, Method A, Method B, SIMC, yield, and
cross-section ownership separate. Keep proton before pion, Method-A/B
independence, same-t Method B, Phase-E presentation-only ownership, and
parent-only normalization. Method B is diagnostic only. F.1-F.5 are detached;
no production promotion before F.6.

## Gaps and exact next action

DEFERRED — final Phase-D provenance, four E.3 settings, older proton detail,
and raw historical bundle paths/hashes. Read investigations/KNOWN_GAPS.md.

NEXT — Run only the detached F.2 analyzer on the farm against the accepted five
F.1 artifacts. Inspect its JSON, four-page PDF, source hashes, per-group
response/support/OOD evidence, and recommendation. Accept a basis or make one
coherent F.2 repair; F.3 remains blocked until acceptance.
