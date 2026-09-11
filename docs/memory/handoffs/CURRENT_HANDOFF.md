# Current KaonLT handoff

Prepared: 2026-09-11 after historical-memory migration and F.1.Fix.5 source
review. This record distinguishes recovered farm evidence from live source
state.

## Resume identity

- Repository: trottar/lt_analysis; branch: test.
- Recheck git status --short --branch, git rev-parse HEAD, and git rev-parse
  origin/test before work.
- Migration live HEAD: b02316f83bf8ef18641fa7217b90c04bcdca10e3.
- F.1.Fix.5 analysis source is introduced by
  dc4fc6283001739a487ec80068f951b0e388cae6; later commits through live HEAD
  carry memory/import material.

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

## Current blocker and boundaries

ACTIVE — Existing F.1 collector/profile identity policy requires
d656e15761970d7d612bb028d2746d077795e9ad and rejects later Fix.5 analysis
files. This is a static source-identity mismatch, not a failed farm run or
physics defect. Read investigations/2026-09-11-phase-f1-source-identity-reconciliation.md.

Keep random, slow-proton, pion, Method A, Method B, SIMC, yield, and
cross-section ownership separate. Keep proton before pion, Method-A/B
independence, same-t Method B, Phase-E presentation-only ownership, and
parent-only normalization. Method B is diagnostic only. F.1-F.5 are detached;
no production promotion before F.6.

## Gaps and exact next action

DEFERRED — final Phase-D provenance, four E.3 settings, older proton detail,
raw historical bundle paths/hashes, and all F.1.Fix.5 farm evidence. Read
investigations/KNOWN_GAPS.md.

NEXT — Write/review a narrow F.1.Fix.5 source-identity reconciliation contract
that selects the reviewed Fix.5 analysis commit and aligns the profile/collector
rule. Only after that, the user performs one fresh five-setting F.1 v2 farm
gate with declared artifact, checker/provenance, and rendered-five-page review.
