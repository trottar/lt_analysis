# Durable KaonLT project knowledge

## Authority and status

docs/memory is the authoritative project development record; native Codex
memory is supplemental only. Start substantial work by reading CURRENT.md,
this file, and handoffs/CURRENT_HANDOFF.md; inspect relevant deep records; then
compare their source identity with live test.

For implementation: current source > exact newer source/diff > newest
authoritative handoff > older history > inference. For runtime: fresh farm
evidence > reviewed authoritative farm handoff > source evidence > history >
inference. Use only CLOSED / RUNTIME VALIDATED, SOURCE REVIEWED, DEVELOPMENT
COMPLETE, FARM VALIDATION PENDING, ACTIVE, DEFERRED, BLOCKED, and NEXT for
work-state claims.

## Source-identity semantics

- Record the reviewed analysis-source commit, the repository HEAD observed at
  the last reconciliation, and the farm-evaluated commit as separate values.
  A source review may be anchored to an analysis-source commit older than the
  observed repository HEAD.
- An observed repository HEAD is a timestamped checkout observation, never a
  permanent assertion that the same commit remains live. Every substantial task
  must establish the actual live `test` HEAD and working-tree state before work.
- Record a farm-evaluated commit only when the applicable farm artifact/evidence
  identifies it; otherwise state that no farm-evaluated commit is recorded.
- Documentation/memory-only commits after a reviewed analysis-source commit do
  not by themselves invalidate that source review. Reconcile the intervening
  diff and re-review only if the relevant analysis source or test scope changed.
- F.1.Fix.5 validation profile/bundle schema v4 pins
  dc4fc6283001739a487ec80068f951b0e388cae6. Its source-identity exception is
  deliberately narrow: three exact collector/profile/test files and the
  docs/memory/ prefix only. It must continue to reject tracked AGENTS.md and
  all other src/ and testing/ paths. The v4 collector validates the producer's v2 training and
  application populations separately; local collector checks are not farm
  evidence.

## Farm validation-bundle procedure

Use the authoritative [farm validation-bundle procedure](decisions/farm-validation-bundle-procedure.md)
for Jefferson Lab `tcsh` packaging. It preserves canonical farm paths,
provenance-safe temporary detached worktrees when a normal checkout has
advanced, bundle-only non-rerun behavior, collision handling, and the rule
that ordinary unrelated farm-local files must not be cleaned merely for a
detached bundle.

## Production and scientific boundaries

- Preserve the separate scientific owners for random subtraction, slow-proton
  PID contamination, pion-production background, HGCer diagnostics, SIMC,
  yields, and cross sections.
- Preserve random/dummy -> frozen binning -> slow proton -> pion subtraction.
  A setting-wide K Lambda gate controls whether proposed proton weights become
  applied; do not partially commit per-t results.
- Pion alignment comparisons use fixed evaluation envelopes. HGCer diagnostics
  consume the frozen baseline rather than redefine it.
- Preserve cuts, templates, priors, component definitions, normalizations,
  binning, efficiencies, acceptance, L/T separation, and uncertainties unless
  a narrow approved contract owns the change.

## HGCer boundaries

- Method A and Method B remain independent. Method B uses frozen Phase-A
  records and same-canonical-t relative closure: no Method-A numbers, cross-t
  pooling/interpolation, or absolute neutron normalization.
- Adaptive Method B is DO NOT PROMOTE: legacy Method B remains Phase-D context;
  neither changes production pion subtraction.
- Persisted diagnostics require producer -> serializer/checkpoint ->
  checkpoint-first payload -> consumer -> renderer review.
- Phase D compares frozen A/B states. Phase E is presentation-only and cannot
  recompute a diagnostic or construct a correction.
- Method B is diagnostic/cross-check/historical comparison only. Method A,
  which excludes NPE=0, is merely the future candidate numerical HGCer input
  and remains a positive-response relative diagnostic.
- Method-A response training and downstream application are separate
  populations. Train from prompt/noRF/nommcuts NPE>0 records; application uses
  the authoritative physical NPE>2 cache. Keep separate provenance and
  fingerprints. Future response application is event-level and parent-t
  normalized; never renormalize (t,phi) children separately.
- F.1-F.5 are detached. No production promotion before a separately validated
  and explicitly approved F.6 decision.

## Current evidence boundary

Phase C five-setting closure, E.3.Fix.2 Left-low, and F.1.Fix.5 are the
recovered phase-level CLOSED / RUNTIME VALIDATED results.

CLOSED / RUNTIME VALIDATED — F.1.Fix.5 validates the detached v2
dual-population artifacts for all five canonical Q4p4W2p74 settings. Preserve
four identities: reviewed analysis source
`dc4fc6283001739a487ec80068f951b0e388cae6`; runtime-evaluated analysis HEAD
`126fa22c19bd29b9952f55b33ab43d59f9727ef6`; collector reconciliation commit
`31dd034d8404e317863bf0933c51253e1d3deeb8`; and this docs/memory-only
closure commit, which is the F.2 baseline. The owner accepted the C1 clean
collector gate and continuity of the existing JSON/PDF/page-manifest evidence
with the supplied v4 runtime bundle. The C1 range whitespace check covers only
the profile's three explicit validation files; the global worktree check and
independent committed-file identity audit remain strict. F.1 remains detached;
F.2-F.5 remain detached and F.6 is the only possible production promotion.

CLOSED / RUNTIME VALIDATED — F.2 was implemented from C2
`02f21886cff8a721df46ebf673594856583fc3c5` as a pure Python/numpy/scipy
five-artifact representation audit. Its fixed candidates are `delta_only`,
`track3`, `hgcer3`, and non-promotable `full5_reference`; response probes are
canonical-t local, deterministic five-fold class-balanced L-BFGS-B diagnostics.
The accepted farm bundle ran at `8e919fc618cea900227db5090d65da728c3aa555`;
all 15 `hgcer3` parents passed fixed information/support gates and the owner
accepted `hgcer3 = (SHMS_delta, P_hgcer_xAtCer, P_hgcer_yAtCer)` for F.3.
The historical F.2 artifact remains `basis_frozen = false`; F.3 freezes its
own basis only from the explicit acceptance. See evidence/f2-fix1-runtime-closure.md.

CLOSED / RUNTIME VALIDATED — F.3 is a detached,
aggregate-only support-aware relative response map over the accepted ordered
`hgcer3` basis. It locally repeats strict F.1 v2 authority validation, verifies
the F.2 wrapper/representation fingerprints and accepted recommendation, then
fits one robust-scaled class-balanced L-BFGS-B model per setting × canonical-t
parent. Its support/OOD calculation must reproduce every F.2 `hgcer3` result
within `1e-12`; a single invalid parent makes the global map unavailable. It
persists coefficients/scalers/support aggregates only, while review grids are
transient cKDTree-masked PDF rendering. No absolute probability, correction,
event output, normalization, Method-B input, or production behavior exists.
The accepted F.3 bundle ran at `5382cfc1994b078c620b32c043938134c33ffa39`;
all 15 parents had exact F.2 continuity and the accepted map is frozen for
F.4. See evidence/f3-runtime-closure.md.

CLOSED / RUNTIME VALIDATED — F.4 is the sole owner of the detached
parent-preserving Method-A correction. Its accepted bundle ran at
`67e0298c51759c7a5ba693464d2c2655bf39250d`, pins accepted F.3/F.1 evidence,
uses raw `A=1` outside frozen support, and applies exactly one full signed
canonical-t parent normalization. All 15 parents passed support and closure;
no event-level correction, template, yield, normal-analysis, or production
behavior exists. See evidence/f4-runtime-closure.md.

CLOSED / RUNTIME VALIDATED — F.5 is a detached signed `(t,phi)` propagation
review of accepted F.4. Its accepted bundle ran at
`3c6a66b7df9bf17e5a428458a2281a80831f001a`: 135 canonical cells (115
occupied, 20 explicit empty), exact F.1 baseline/F.4 child propagation, and
all parent/setting closures passed. Its scientific fingerprint is
`d11b728d1089301a12c29e7f8b1798c6e0b6021ac47bd5d2afc2b62b47a1effa`; see
evidence/f5-runtime-closure.md. SOURCE REVIEWED — F.5.Fix.1 began from F.5
implementation `f9ce5aadbebc31f5163a093ab32244644ad2a088` and was implemented/
reviewed at `46598878c102a67275b1600e73cba3b2f166dc26`; it repaired only
farm-Python `zip(strict=True)` compatibility after the unchanged exact
row/factor-length gate.

CLOSED / RUNTIME VALIDATED — F.5.2 is the accepted presentation-only farm
rerender at `6634e9cb470cf35f21f5d475ec6ce33b524cd233`, packaged as
`KaonLT_PhaseF5_2_Fix1_validation_Q4p4W2p74.zip` with SHA-256
`323d5d50604d846692c2b0bb0d1455c3751b2033601046dfcc3e07ccc3ca3ebd`.
The accepted F.5 scientific fingerprint remains exactly
`d11b728d1089301a12c29e7f8b1798c6e0b6021ac47bd5d2afc2b62b47a1effa`; the
propagation payload was exactly unchanged. The 12-page review passed visual
acceptance: physics-readable terms and page-12 physical phi intervals replaced
internal indices. No scientific calculation, payload, authority, geometry,
normalization, production behavior, or F.5 result changed. F.5.2.Fix.1 at
`6634e9cb470cf35f21f5d475ec6ce33b524cd233` remains distinct from the earlier
F.5.Fix.1 source-reviewed compatibility repair. See
evidence/f5-2-runtime-closure.md.

F.6 promotion is deliberately staged. CLOSED / RUNTIME VALIDATED — F.6.1 and
F.6.1.Validation.1 are accepted from
`KaonLT_PhaseF6_1_validation_Q4p4W2p74.zip` (SHA-256
`22fb594f79a61e5378e76570df46cab9cf5195b6dc0faa04b0238919f66738a3`) at
farm/bundle HEAD `ffb7d4dc251610f6034eac11038a9841e3ef9f58`, with reviewed
scientific source pin `bfc4fe421f9fc9139a1992a0ec92e31aa101b86c`. The complete
bundle has exact accepted F.4/F.5 reproduction, five settings, 15 parents,
135 cells, all frozen F.1 hashes, and the explicit training-radians to
application-degrees phi closure. It is aggregate-only and non-production; its
27-page review was manually inspected. Method A is technically closed but has
parent-dependent shape effects: the 68/135 toward and 67/135 away comparison
counts are descriptive only. See evidence/f6-1-runtime-closure.md.

CLOSED / RUNTIME VALIDATED — F.6.2 is the detached acceptance-correlated
Method-A refinement scientific validation accepted from the previously reviewed
complete farm bundle. It persists raw discrepancy, refinement, alignment,
independent acceptance/MM, kaon-window, support, and bootstrap evidence without
a composite score or frozen numerical thresholds. It remains detached: Method B
is excluded, no production mutation is authorized, and no final-yield
uncertainty reduction is claimed. See
evidence/f6-2-scientific-runtime-closure.md.

CLOSED / RUNTIME VALIDATED — F.6.2.Fix.5 accepted the presentation-only
rerender and closed F.6.2 overall. The frozen scientific JSON SHA-256
`5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1` and
artifact/validation fingerprints remain unchanged. Keep the reviewed scientific
source `0b37af2a2927b08bdeaf897c545f290b55329cea` distinct from Fix.5 renderer
source `c88ed65cb18ba6a37358897292b77696016312d1` and bundle/profile commit
`b789d203e11f0927deb59ebcae9dc59fe8add4ae`, which are presentation-artifact
provenance only. See evidence/f6-2-fix5-presentation-runtime-closure.md.

NEXT — E.8 may receive a separate sidecar figure-library contract only after
the user commits/pushes this reconciliation and a new session re-establishes
live Git state. It must read the frozen F.6.2 JSON only, verify its SHA before
parsing, render persisted quantities only, and never recompute Method-A,
bootstrap, support/OOD, cuts, weights, binning, yields, or scientific
acceptance. F.6.3 is BLOCKED pending E.8; F.6.4 is BLOCKED pending F.6.3
evidence. See decisions/f6-2-acceptance-refinement-measurement-contract.md and
phases/phase-f6-method-a-production-promotion.md.

F.2.Fix.1 starts from the clean F.2 implementation commit
`5549098d2552b9c092b65e31aeb77edc0807ddda`. It is an authority/provenance and
deterministic-failure repair only: independently reconstruct the frozen F.1 v2
training, application, feature-metadata, child-projection, and full-contract
fingerprints; require exact F.1 authority metadata, `post_proton_noRF`, and
finite strictly increasing geometry; bind deterministic filenames to declared
payload identities before hashing; preserve the four approved recommendation
states; mark zero non-prompt rows sparse; and reject JSON/PDF path collisions.
It must not change candidates, probe numerics, information/OOD gates, or any
production/normal-analysis path.

SOURCE REVIEWED — The single validation collector is profile-driven for simple
future evidence packages. The original F.1 profile remains in its specialized
validation/PDF-extraction mode; generic profiles declare only `global` and
`settings` artifact lists with a deterministic basename, `json`/`file` kind,
and required flag. The F.2 profile pins source `170e6fae3d2fed1949fc6932b8eac9ad83e3e01c`
and packages one global representation JSON, one global four-page PDF, and the
five setting-scoped F.1 v2 JSON inputs. It never runs the analyzer or imports
analysis physics. Its source-identity exception is limited to the collector,
collector test, two collector profiles, and `docs/memory/`.
