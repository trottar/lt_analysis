# F.6.2 scientific runtime closure

## Scope and status

CLOSED / RUNTIME VALIDATED — the detached F.6.2 acceptance-correlated
Method-A refinement scientific gate was accepted from the previously reviewed
complete farm bundle. This record documents that accepted scientific gate only.
It does not accept the later F.6.2.Fix.4 presentation-only rerender or promote
Method A into production.

## Reviewed identity and continuity

- Reviewed F.6.2.Fix.4 presentation source:
  `b929815517643edaa75949f7e61ce46c6e8f0d63`.
- Reviewed clean bundle/profile commit:
  `c3252f1ec43dcca0357aebe8c6ca75d7d198f99a`.
- Accepted F.6.2 scientific JSON SHA-256:
  `5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1`.
- Accepted F.6.2 artifact fingerprint:
  `ee713b70de898bad8fa61164cbc8a54886af1ea1eb1e712ed95de17df8890cd0`.
- Accepted F.6.2 validation fingerprint:
  `7edc73fce20ad7dc7622b8c23a7ba7e8986e367ace605884e010945595370b3b`.

The recovered repository record does not contain the accepted bundle archive
path, archive SHA-256, or farm execution HEAD. They are **not recovered** here
rather than inferred. The accepted identities above are also retained in the
[farm validation-bundle procedure](../decisions/farm-validation-bundle-procedure.md).

## Preserved boundaries and remaining gate

The accepted scientific validation remains detached and aggregate-only: it
does not alter production pion weights, yields, ROOT objects, Method-B input,
or final-yield uncertainty. F.6.2 overall remains ACTIVE because
F.6.2.Fix.4 is DEVELOPMENT COMPLETE, FARM VALIDATION PENDING. The next gate is
a fresh presentation-only rerender bundle that preserves the accepted F.6.2
scientific JSON and fingerprints unchanged.
