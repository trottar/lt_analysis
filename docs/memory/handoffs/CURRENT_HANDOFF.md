# Current KaonLT handoff

Prepared: 2026-09-11 after F.1.Fix.5 validation-gate reconciliation. This
record distinguishes recovered farm evidence from live source state.

## Resume identity

- Repository: trottar/lt_analysis; branch: test.
- Recheck git status --short --branch, git rev-parse HEAD, and git rev-parse
  origin/test before work.
- Repository HEAD observed at the last reconciliation (2026-09-11):
  47723e02af03b3844dfb0323560ce39f7dfebe1a. This is a timestamped observation,
  not a permanent live-HEAD claim.
- Reviewed F.1.Fix.5 analysis-source commit:
  dc4fc6283001739a487ec80068f951b0e388cae6.
- F.1.Fix.5 farm-evaluated commit: none recorded; no F.1.Fix.5 farm artifact
  is present.
- At the last reconciliation, later commits through the observed HEAD contained
  documentation/memory, repository-guidance, or archived-memory material only;
  they did not change the reviewed F.1 analysis source or tests. Those commits
  do not by themselves invalidate the source review, but the actual live test
  HEAD and working tree must still be established before work.

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

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING — The F.1 collector/profile v4
identity policy requires reviewed Fix.5 source commit
dc4fc6283001739a487ec80068f951b0e388cae6. After it, only the three detached
collector/profile/test files and docs/memory/ are allowed; tracked root
AGENTS.md plus all other analysis and test changes are identity errors. Its artifact validator
requires the F.1 v2 dual-population contract, preserving separate positive-NPE
Method-A training and physical NPE>2 application records, summaries, and
fingerprints. Local detached checks passed; this is not a farm run or physics
claim. Read investigations/2026-09-11-phase-f1-source-identity-reconciliation.md.

Keep random, slow-proton, pion, Method A, Method B, SIMC, yield, and
cross-section ownership separate. Keep proton before pion, Method-A/B
independence, same-t Method B, Phase-E presentation-only ownership, and
parent-only normalization. Method B is diagnostic only. F.1-F.5 are detached;
no production promotion before F.6.

## Gaps and exact next action

DEFERRED — final Phase-D provenance, four E.3 settings, older proton detail,
raw historical bundle paths/hashes, and all F.1.Fix.5 farm evidence. Read
investigations/KNOWN_GAPS.md.

NEXT — The user performs one fresh five-setting F.1.Fix.5 v2 farm gate: Left
lowe, Left highe, Center lowe, Center highe, and Right highe. Collect fresh
artifacts with the v4 profile, inspect provenance/checker output, and review
the five rendered F.1 pages before recording any farm result.
