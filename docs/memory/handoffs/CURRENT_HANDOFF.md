# Current KaonLT handoff

Prepared: 2026-09-15 after accepted F.4 closure and F.5 source implementation. This record
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
- F.2 baseline is C2 `02f21886cff8a721df46ebf673594856583fc3c5`; its source
  implementation commit is `5549098d2552b9c092b65e31aeb77edc0807ddda`.
  F.2.Fix.1 begins exactly at the latter clean checkout and is not a runtime
  or farm-validation identity.
- The profile-driven F.2 collector repair begins exactly at
  `170e6fae3d2fed1949fc6932b8eac9ad83e3e01c`; it is source-only and does not
  replace either F.2 source identity or any farm identity.
- F.5 implementation and F.5.Fix.1 starting HEAD:
  `f9ce5aadbebc31f5163a093ab32244644ad2a088`.
- F.5.Fix.1 implementation/source-reviewed HEAD:
  `46598878c102a67275b1600e73cba3b2f166dc26`. This is the farm-Python
  compatibility repair only; it is not a F.5 farm-evaluated identity.

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

CLOSED / RUNTIME VALIDATED — F.2.Fix.1 is accepted from the supplied complete
F.2 farm bundle at `8e919fc618cea900227db5090d65da728c3aa555`. The owner
accepted `hgcer3 = (SHMS_delta, P_hgcer_xAtCer, P_hgcer_yAtCer)` after all 15
groups passed frozen `hgcer3` information and support gates. The F.2 artifact
remains historically `basis_frozen = false`; this is an external human
acceptance, not a rewrite of F.2. Read evidence/f2-fix1-runtime-closure.md.

CLOSED / RUNTIME VALIDATED — F.3 is accepted at farm/source HEAD
`5382cfc1994b078c620b32c043938134c33ffa39`; its 15 parent maps have exact
F.2 support continuity and no correction or production behavior. Read
evidence/f3-runtime-closure.md.

CLOSED / RUNTIME VALIDATED — F.4 is accepted from the supplied complete F.4
farm bundle at `67e0298c51759c7a5ba693464d2c2655bf39250d`; the exact bundle,
JSON/PDF hashes, and correction/artifact fingerprints are in
evidence/f4-runtime-closure.md. It preserves full signed canonical-t sums with
one common parent normalization and has no template or production ownership.

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING — F.5 consumes accepted
F.4/F.3/F.1 artifacts only and builds detached signed canonical `(t,phi)`
aggregate templates. It must use the F.4 shared calculator, retain explicit
zero cells, and never normalize children or create event/production output.
SOURCE REVIEWED — F.5.Fix.1 starts from the F.5 implementation
`f9ce5aadbebc31f5163a093ab32244644ad2a088` and is implemented/reviewed at
`46598878c102a67275b1600e73cba3b2f166dc26`. It removes only redundant
`zip(strict=True)` after the exact existing row/factor length guard so the
detached analyzer is compatible with farm Python.

SOURCE REVIEWED — `testing/collect_pion_hgcer_validation_bundle.py` remains
the only collector. Its original F.1-specialized profile behavior is preserved;
the new generic profile mode consumes declared global and setting-scoped files
without F.1 checkpoint/PDF extraction. The F.2 profile requires only the F.2
JSON/PDF once and five F.1 v2 JSON inputs. It is read-only with respect to
analysis artifacts and does not run F.2.

Keep random, slow-proton, pion, Method A, Method B, SIMC, yield, and
cross-section ownership separate. Keep proton before pion, Method-A/B
independence, same-t Method B, Phase-E presentation-only ownership, and
parent-only normalization. Method B is diagnostic only. F.1-F.5 are detached;
no production promotion before F.6.

## Gaps and exact next action

DEFERRED — final Phase-D provenance, four E.3 settings, older proton detail,
and raw historical bundle paths/hashes. Read investigations/KNOWN_GAPS.md.

NEXT — Run F.5's single detached analyzer/unchanged generic collector gate
against accepted F.1/F.3/F.4 artifacts. F.6 remains BLOCKED pending explicit
F.5 review.
