# Current KaonLT handoff

Prepared: 2026-09-14 after F.1.Fix.5 v4 bundle inspection. This record
distinguishes serialized-artifact closure from accepted clean-source farm
validation.

## Resume identity

- Repository: trottar/lt_analysis; branch: test.
- Recheck git status --short --branch, git rev-parse HEAD, and git rev-parse
  origin/test before work.
- Repository HEAD observed during this evidence inspection (2026-09-14):
  126fa22c19bd29b9952f55b33ab43d59f9727ef6. This is a timestamped observation,
  not a permanent live-HEAD claim.
- Reviewed F.1.Fix.5 analysis-source commit:
  dc4fc6283001739a487ec80068f951b0e388cae6.
- F.1.Fix.5 candidate farm HEAD recorded in the supplied v4 archive:
  126fa22c19bd29b9952f55b33ab43d59f9727ef6. It is not accepted closure; read
  evidence/f1-fix5-v4-bundle-inspection-2026-09-14.md.
- At this reconciliation, later commits through the observed HEAD contain
  documentation/memory/repository guidance or the three exact detached F.1
  collector/profile/test validation files; they contain no F.1 analysis-source
  changes. Those commits do not by themselves invalidate the source review, but
  the actual live test HEAD and working tree must still be established before
  work.

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

SOURCE REVIEWED:

- F.1.Fix.5 has separate Part-1 prompt/noRF/nommcuts NPE>0 Method-A training
  and physical NPE>2 application populations. It maintains independent
  summaries/fingerprints; pages 1-4 are training and page 5 application.
  It is non-authoritative, Method-B independent, and production-side-effect
  free. Read phases/phase-f1-method-a-acceptance-farm-gate.md.

## Current gate and boundaries

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING — The supplied five-setting v4
archive directly validates each F.1 v2 artifact's hash, schema, detached flags,
dual-population selection, prompt identity closure, fingerprints, and 30
canonical `(t, delta)` closures. The collector repair excludes only allowed
`docs/memory/**` from its historical range whitespace check while retaining
strict analysis and path identity gates; local collector and frozen F.1 tests
pass. Per user direction, re-run only the collector from a clean farm checkout
against the existing outputs. Do not call F.1 closed or start F.2 until its
complete manifest is inspected. Read
evidence/f1-fix5-v4-bundle-inspection-2026-09-14.md.

Keep random, slow-proton, pion, Method A, Method B, SIMC, yield, and
cross-section ownership separate. Keep proton before pion, Method-A/B
independence, same-t Method B, Phase-E presentation-only ownership, and
parent-only normalization. Method B is diagnostic only. F.1-F.5 are detached;
no production promotion before F.6.

## Gaps and exact next action

DEFERRED — final Phase-D provenance, four E.3 settings, older proton detail,
and raw historical bundle paths/hashes. F.1 clean-source runtime provenance is
unresolved. Read investigations/KNOWN_GAPS.md.

NEXT — Run the repaired F.1 collector only from a clean farm checkout against
the existing five `Q4p4W2p74` outputs. Inspect its complete manifest and
provenance, then close F.1 or make one coherent repair before F.2.
