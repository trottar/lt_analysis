# F.6.2 scientific runtime closure

## Scope and status

CLOSED / RUNTIME VALIDATED — the detached F.6.2 acceptance-correlated
Method-A refinement scientific gate was accepted from the previously reviewed
complete farm bundle. This record documents that accepted scientific gate only.
It does not accept the later F.6.2.Fix.4 presentation-only rerender or promote
Method A into production.

## Accepted scientific identity and continuity

- Reviewed F.6.2 scientific source:
  `0b37af2a2927b08bdeaf897c545f290b55329cea`.
- Accepted F.6.2 scientific bundle/profile HEAD:
  `eb8723bc047af7369865ea1c16b97286173b0a59`.
- Accepted F.6.2 scientific JSON SHA-256:
  `5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1`.
- Accepted F.6.2 artifact fingerprint:
  `ee713b70de898bad8fa61164cbc8a54886af1ea1eb1e712ed95de17df8890cd0`.
- Accepted F.6.2 validation fingerprint:
  `7edc73fce20ad7dc7622b8c23a7ba7e8986e367ace605884e010945595370b3b`.
- Reviewed uploaded scientific bundle SHA-256:
  `a116d502a8dc637e0d2e3b58dd69ffa34a6b25fbd4c65cb27895d01e2f7b09b3`.
- Pre-Fix.4 accepted review-PDF SHA-256 (presentation identity, not scientific
  payload identity):
  `4d07c2989c217d18d088b181390e9478fae51e4168dba2d64faea690b40f0f85`.

No analyzer farm-execution HEAD was explicitly recorded, so none is inferred
here. The accepted identities above are also retained in the
[farm validation-bundle procedure](../decisions/farm-validation-bundle-procedure.md).

## Later presentation-only provenance

F.6.2.Fix.4 source `b929815517643edaa75949f7e61ce46c6e8f0d63` was superseded
by the accepted Fix.5 presentation rerender. Fix.5 used renderer source
`c88ed65cb18ba6a37358897292b77696016312d1` and bundle/profile commit
`b789d203e11f0927deb59ebcae9dc59fe8add4ae`; its accepted PDF SHA-256 is
`f33406927b78a7b27d863c3dce6d967956fcfe1e9326b5bd2f09e3ec5bbda408`.
These are presentation-artifact provenance only, not replacements for the
scientific source or scientific bundle identity. See the separate
[Fix.5 presentation closure](f6-2-fix5-presentation-runtime-closure.md).

## Preserved boundaries and downstream state

The accepted scientific validation remains detached and aggregate-only: it
does not alter production pion weights, yields, ROOT objects, Method-B input,
or final-yield uncertainty. F.6.2 overall is CLOSED / RUNTIME VALIDATED after
the accepted Fix.5 presentation gate. E.8 is NEXT and must consume the frozen
scientific JSON/fingerprints without recalculation; F.6.3 remains BLOCKED
pending E.8, and F.6.4 remains BLOCKED pending F.6.3 evidence.
