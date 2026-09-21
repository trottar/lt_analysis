---
memory_schema: 2
active_objective: Farm-validate the final E.8.1 streamlined full-background-subtraction procedure PDF
current_work_item: Run the canonical-five-setting E.8.1 farm gate and inspect Left/highe first
active_status: DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
next_action: Run the five-setting Q4p4W2p74 E.8.1 farm gate with the new profile, collect the bundle, and inspect Left/highe first
baseline_commit: c25f9d8248f2acfb8b0443a482fc76781eb405ca
source_commit: 0bf1281ebf44868b68ff808090436653eea6ed60
bundle_profile_commit: c25f9d8248f2acfb8b0443a482fc76781eb405ca
---
# Current KaonLT development state

## Active Objective

Farm-validate the final E.8.1 streamlined ordinary full-background-subtraction
procedure PDF through the reviewed canonical-five-setting evidence
bundle/profile.

## Current Work Item

F.6.2 and F.6.2.Fix.5 are CLOSED / RUNTIME VALIDATED; see the
[scientific closure record](evidence/f6-2-scientific-runtime-closure.md) and
[Fix.5 closure record](evidence/f6-2-fix5-presentation-runtime-closure.md).
The standalone E.8 frozen-F.6.2 figure library was completed as an intermediate
presentation artifact. The final E.8.1 ordinary procedure-PDF source
`0bf1281ebf44868b68ff808090436653eea6ed60` and canonical-five-setting
bundle/profile infrastructure `c25f9d8248f2acfb8b0443a482fc76781eb405ca` are
SOURCE REVIEWED.

## Verified State

- F.6.2 scientific validation is CLOSED / RUNTIME VALIDATED from its accepted
  scientific bundle; see the [scientific closure record](evidence/f6-2-scientific-runtime-closure.md).
- F.6.2.Fix.5 presentation repair is CLOSED / RUNTIME VALIDATED from the
  accepted fresh presentation bundle; see the [Fix.5 closure record](evidence/f6-2-fix5-presentation-runtime-closure.md).
- The accepted yield remains the baseline. Method B remains diagnostic only;
  Method A remains detached and non-production.

## Source / Evidence Identity

The reviewed F.6.2 scientific source is
`0b37af2a2927b08bdeaf897c545f290b55329cea`. The accepted scientific JSON
SHA-256 is
`5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1`, with
artifact/validation fingerprints `ee713b70de898bad8fa61164cbc8a54886af1ea1eb1e712ed95de17df8890cd0` /
`7edc73fce20ad7dc7622b8c23a7ba7e8986e367ace605884e010945595370b3b`.

The accepted Fix.5 presentation source
`c88ed65cb18ba6a37358897292b77696016312d1` and bundle/profile commit
`b789d203e11f0927deb59ebcae9dc59fe8add4ae` are upstream presentation-artifact
provenance only; they do not replace the scientific source identity.

The E.8.1 initial integration source is
`33c114e003421690fa9b15480a4200634accf8c3`; the semantics repair is
`b02d939d9187ec656aa629f03e528fde0ca0cc0d`; and the sparse-legend repair is
the reviewed procedure-PDF source `0bf1281ebf44868b68ff808090436653eea6ed60`.
The canonical-five-setting E.8.1 bundle/profile commit is
`c25f9d8248f2acfb8b0443a482fc76781eb405ca`; it is SOURCE REVIEWED
evidence-packaging infrastructure and does not replace the reviewed analysis
source. No E.8.1 farm-evaluated commit is recorded.

## Blockers

E.8 overall is DEVELOPMENT COMPLETE, FARM VALIDATION PENDING. Source review is
not farm readiness: before requesting a farm run, verify that the current
validation bundle/profile pins the reviewed source and declares every artifact
required for the gate. That requirement is satisfied by E.8.1 bundle/profile
commit `c25f9d8248f2acfb8b0443a482fc76781eb405ca`. F.6.3 remains BLOCKED
pending the targeted E.8.1 farm gate; F.6.4 remains BLOCKED pending F.6.3
evidence.

## Next Action

NEXT — run the established five-setting `Q4p4W2p74` E.8.1 farm production with
`testing/pion_hgcer_validation_bundle_profile_e8_1.json`, collect the bundle,
and inspect `Left / highe` first. Its ordinary procedure PDF and page manifest
must show D.6-D.9, the compact E.8 context plus three atlas pages per persisted
canonical-t parent and handoff, a visible L/B/A key, sparse-state pages without
parent suppression, and the explicit F.6.3 handoff before review broadens to
the remaining settings.

## Success Criteria

The bundle must contain the frozen accepted F.6.2 JSON and ordinary procedure
PDF/page-manifest pairs for `Left / lowe`, `Left / highe`, `Center / lowe`,
`Center / highe`, and `Right / highe`. The first detailed acceptance inspection
remains `Q4p4W2p74 / Left / highe`; it does not redefine the collector/profile
as a one-setting package.

## Do Not Reopen Without New Evidence

Do not reopen accepted F.1 through F.6.2 scientific validation, alter the
baseline yield, or promote Method A. E.8 must consume the frozen F.6.2 JSON
only; it cannot recalculate Method-A science, support/OOD, weights, or yields.
No Method-A production promotion has occurred.

## Relevant References

- [F.6.2 scientific runtime closure](evidence/f6-2-scientific-runtime-closure.md)
- [F.6.2.Fix.5 presentation runtime closure](evidence/f6-2-fix5-presentation-runtime-closure.md)
- [F.6.2 measurement contract](decisions/f6-2-acceptance-refinement-measurement-contract.md)
- [F.6 production-promotion phase record](phases/phase-f6-method-a-production-promotion.md)
- [Farm validation-bundle procedure](decisions/farm-validation-bundle-procedure.md)
