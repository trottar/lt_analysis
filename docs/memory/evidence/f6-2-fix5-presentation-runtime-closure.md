# F.6.2.Fix.5 presentation runtime closure

## Scope and status

CLOSED / RUNTIME VALIDATED — the fresh F.6.2.Fix.5 presentation-only rerender
was accepted from the reviewed complete farm package. It closes the F.6.2
presentation gate and therefore F.6.2 overall. It is not a scientific rerun,
does not promote Method A, and does not replace the accepted scientific source
or payload provenance.

## Accepted package and presentation provenance

- Accepted package: `KaonLT_PhaseF6_2_Fix5_validation_Q4p4W2p74_20260918-021700.zip`.
- Delivered ZIP SHA-256:
  `168ff1b682b9e9dea85570323e0ca1071fd3efd1ef5767eee21569782a429aca`.
- Fix.5 renderer source:
  `c88ed65cb18ba6a37358897292b77696016312d1`.
- Fix.5 bundle/profile commit:
  `b789d203e11f0927deb59ebcae9dc59fe8add4ae`.
- Superseded Fix.4 presentation source:
  `b929815517643edaa75949f7e61ce46c6e8f0d63`.
- Prior Fix.4 PDF SHA-256:
  `a014ba89d36a38d1387448b5176b17e932f1257ab0c8fa6b605271009758358a`.
- Accepted Fix.5 PDF SHA-256:
  `f33406927b78a7b27d863c3dce6d967956fcfe1e9326b5bd2f09e3ec5bbda408`.

The original farm archive path was not recovered and is not inferred. The ZIP
identity above is the reviewed delivered bundle identity.

## Frozen scientific continuity

The presentation rerender preserved the accepted F.6.2 scientific JSON
SHA-256 exactly:

`5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1`.

It also preserved the accepted artifact fingerprint
`ee713b70de898bad8fa61164cbc8a54886af1ea1eb1e712ed95de17df8890cd0` and
validation fingerprint
`7edc73fce20ad7dc7622b8c23a7ba7e8986e367ace605884e010945595370b3b`.
The reviewed F.6.2 scientific source remains
`0b37af2a2927b08bdeaf897c545f290b55329cea`; Fix.5 source/profile identities
are upstream presentation-artifact provenance only.

## Bundle and visual review

- The supplied manifest was complete with `errors=[]`; source-ancestor and
  integrity checks passed, the detached worktree was clean, and no unexpected
  committed files were reported.
- All 15 parent overview pages passed review: Left low/high and Center were
  populated, Right-high sparse states rendered correctly, labels and footers
  remained readable, and no clipping, overlap, or raw unavailable reason
  appeared in overview tables.
- The accepted PDF has 132 pages. Regression pages 17, 36, 54, 81, 108, 131,
  and 132 were pixel-identical at 150 dpi between Fix.4 and Fix.5.
- Farm checks accepted the full-background-subtraction plot suite (90), the
  Method-A acceptance contract (7), the Phase-F runtime contract (3), and the
  collector (25), alongside compile and diff checks.

## Conclusion and boundary

F.6.2 scientific validation is CLOSED / RUNTIME VALIDATED, F.6.2.Fix.5
presentation repair is CLOSED / RUNTIME VALIDATED, and F.6.2 overall is
CLOSED / RUNTIME VALIDATED. E.8 is NEXT. Any E.8 renderer must read the frozen
scientific JSON only, verify its SHA before parsing, and render persisted
quantities without recomputing Method-A factors, bootstrap results, support or
OOD, cuts, weights, binning, yields, or scientific acceptance.
