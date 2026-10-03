---
memory_schema: 3
---
# Current KaonLT development state

## Active Objective

E.8 is `ACTIVE`: audit the kaon missing-mass and signal-region yield chain
from authoritative objects. Preserve baseline production, accepted upstream
authorities and detached Method-A/Method-B boundaries.

## Current Work Item

The external PID/background-methods research memory checkpoint is `ACTIVE` for
Q4p4W2p74 / Left / lowe Method-A scientific validity. Gates 1–4 are consumed.
The external review was explicitly user-authorized and is completed enough to
alter planning; no implementation or new farm evidence occurred. The
[research investigation](investigations/e8-4-left-lowe-method-a-external-pid-background-methods-review-2026-10-03.md)
preserves stable bibliography, source hierarchy and access limitations. Prior
[Gate-3 findings](investigations/e8-4-left-lowe-method-a-detector-response-source-science-audit-2026-10-03.md)
remain separate from EXTERNAL LITERATURE.

Stale scalar/wrong histogram is not the primary scientific issue. Fail-closed
checks, prior arithmetic and accepted narrow evidence support bookkeeping
consistency for this purpose, not production correctness. Large low-t/t1
redistribution, global-weight differences and pion-background/yield impact
remain scientifically unresolved.

EXTERNAL LITERATURE: the strongest planning lesson is to separate sample
purification, detector response, accepted `w0` physics transfer and normalization/
preservation. No external method validates a specific KaonLT correction.
F.4 signed normalization is a priority diagnostic question, not a proven defect.
Test F.3 weak-positive response as a relative pion-to-kaon mis-ID topology proxy.
HGC-free pion tagging is preferred if available; truncated/censored models are
optional model-dependent cross-checks.

## Verified State

SOURCE VERIFIED (current PID masks; Gate-3 ownership retained): NPE=0 is the
kaon PID category, NPE>0 the pion tree, NPE>2 physical pion control. Training
and application differ; P(NPE=0 | true pion,x) is pion-to-kaon HGC mis-ID into
the kaon-selected sample, not an observed pion class. F.3
`hgcer3` represents relative response, not absolute leakage probability. F.4
preserves signed canonical-t parents, not each child or MM subregion. Baseline
`w0` already owns the established pion-control-to-background transfer.
The configured HGCer hole is excluded from relevant populations. Existing
zero-photoelectron pion transfer machinery is non-authoritative, diagnostic-only
and production-side-effect-free; it is not an approved Method-A replacement.
Slow-proton proposed/applied architecture is analogue-only, not a pion prescription.

RUNTIME VERIFIED (prior supplied accepted reviews, not new artifact inspection):
the current-baseline comparator reproduced F.2/F.3 scientific payloads exactly;
F.4 was the first scientifically changed stage. F.4.Refresh.2 remains
`CLOSED / RUNTIME VALIDATED`. F.6.3/E.8.4 remain `CLOSED / RUNTIME VALIDATED`
only for Left/lowe branch execution, live-cache parity, real child changes and
signed parent preservation. See [branch evidence](evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md).

Fix.5.7 remains `CLOSED / RUNTIME VALIDATED` only for Left/lowe owner/checker
setting provenance through final ZIP; [Fix.5.8](evidence/e8-4-fix5-8-left-lowe-runtime-closure-2026-10-03.md)
remains `CLOSED / RUNTIME VALIDATED` only for Left/lowe presentation legibility.
Its 97-page bundle, structural checks, parent closure and visual review remain
accepted at their recorded scope. External literature adds no runtime validation.

Historical F.1-F.6.2 closures retain their accepted scopes in the
[roadmap](roadmap/STATUS.md). E.8.2/E.8.3 and workflow hardening remain
`SOURCE REVIEWED`; historical E.8.3 F.6.1 and current F.6.3/E.8.4 are distinct
lineages. No accepted authority is replaced.

## Source / Evidence Identity

- Observed research-checkpoint startup branch `test`, HEAD/local `origin/test`:
  `f7934f459b42b190235ffc63e486e7ba0e068e79`; timestamped observation,
  not permanent source or farm authority.
- Accepted Fix.5.8 farm source: `2ddeab47d55edb57d2f313022a948c4376730c19`.
  ZIP `KaonLT_E8_4_Fix5_Left_lowe_Q4p4W2p74_20261003-032934-987132.zip`,
  30851374 bytes; SHA-256
  `966e667b36b2626b5c16f00fc684a099a5d24e5d2015bb547612a1e63eb2a25e`.
- Current F.6.3 F.4 candidate SHA-256:
  `1d545924eba89c7f9ffa28028e307aca9b434a89beec06863cf2893887b6b902`.
  Historical authority hashes and prior source/failed-owner identities remain
  in the canonical evidence and [post-farm investigation](investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md).
- Reviewed Fix.5.4 numerical source: `761fbb6c03d2d7a10bb911cf84e9ba898496fab6`;
  Fix.5.5 presentation source: `df957a6414fc9c515d1f82228517cb801dc90350`.
  Workflow hardening synchronized through `8a6cad9c4c85bbe8fe9b32c15cc31402d4758af9`.

## Blockers

Physical validity and detector-response origin of the large low-t/t1 Method-A
redistribution remain `BLOCKED`. Accepted comparator evidence concentrates the
question at relative HGCer response, current physical application, baseline
`w0` transfer and signed F.4 normalization. INFERENCE: cancellation could amplify
normalization sensitivity; this is not an accepted explanation. NOT VERIFIED:
current-lineage t1 signed/absolute support, source decomposition, raw-response
and correction tails, child/MM-region redistribution, coordinate and support/OOD
dependence, weak-positive proxy validity and HGC-free pion-tag availability.
No particular PMT, mirror, optical, track-geometry or hardware cause is established.

Absolute-SIMC interpretation remains separately `BLOCKED`:
`SIMC_normfac_luminosity_and_charge_units_not_source_proven`. Existing
`iter_weight * normfac / Ncontribute` lacks source-proven luminosity/effective-
charge units. Only the two absolute-SIMC page families/claims are unavailable;
valid current-F.6.3 data/identity/yield/parent-closure payloads remain available.
No incorrect normalization, conversion or amplitude conclusion is established.

## Next Action

NEXT — after independent review, user commit/push and pushed-state
synchronization: define a detached current-lineage diagnostic-measurement
contract with separate tests of (1) F.3 weak-positive response as a relative
pion-to-kaon mis-ID topology proxy across hgcer3 and (2) Left/lowe t1 F.4
signed-normalization sensitivity, before correction redesign or production.

This checkpoint records research; it authorizes no farm run, scientific-source
modification, correction choice or diagnostic implementation contract here.
Canonical-five expansion remains `DEFERRED` by user decision; final E.8/F.6.4
remain `BLOCKED`.

## Success Criteria

Compact allowlisted memory preserves consumed gates, stable external sources,
evidence distinctions and accepted narrow scopes. Manifest, ordinary memory
health and diff checks must pass; independent actual-diff review and
user-controlled synchronization remain required.

## Do Not Reopen Without New Evidence

Keep `no_empirical_residual` and zero legacy residual scales. Method A remains
detached/non-production; Method B diagnostic/cross-check only and numerically
excluded. Freeze random/dummy/slow-proton/baseline-pion subtraction, SIMC,
weights, yields, uncertainties, cuts, templates, priors, binning, efficiencies,
acceptance, L/T and cross sections. Never normalize children independently,
replace accepted authorities or infer promotion from diagnostic acceptance.
E.8 remains a consumer; F.6.3 owns the private branch. Failed artifacts are not
accepted evidence; owner success alone is not scientific/visual acceptance.

## Relevant References

- [Gate-3 source/science investigation](investigations/e8-4-left-lowe-method-a-detector-response-source-science-audit-2026-10-03.md)
- [Current-baseline comparator](evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md)
- [Prior Fix.5.7 evidence](evidence/e8-4-fix5-7-left-lowe-runtime-and-fix5-visual-gate-2026-10-03.md)
- [E.8 procedure roadmap](decisions/e8-full-analysis-procedure-roadmap.md)
