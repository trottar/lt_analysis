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
- F.2 baseline is C2 `02f21886cff8a721df46ebf673594856583fc3c5`; its source
  implementation commit is `5549098d2552b9c092b65e31aeb77edc0807ddda`.
  F.2.Fix.1 begins exactly at the latter clean checkout and is not a runtime
  or farm-validation identity.
- The profile-driven F.2 collector repair begins exactly at
  `170e6fae3d2fed1949fc6932b8eac9ad83e3e01c`; it is source-only and does not
  replace either F.2 source identity or any farm identity.

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

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING — F.3 constructs only the
detached, relative `hgcer3` response map. It revalidates exact F.1/F.2
provenance and reproduces every F.2
same-parent cKDTree support result before making a map available. It may not
create an absolute probability, correction, weight, normalization, event
output, Method-B numerical input, or production object.

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

NEXT — Run only the detached F.3 analyzer against the accepted five F.1 JSON
artifacts and accepted F.2 representation. Run the unchanged collector with
the F.3 profile and review all 15 models, provenance, exact support continuity,
OOD masks, and detached flags. F.4 is BLOCKED until F.3 evidence is accepted.
