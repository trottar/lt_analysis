# Phase F.3 Method-A support-aware acceptance map

CLOSED / RUNTIME VALIDATED — F.3 started from clean source baseline
`8e919fc618cea900227db5090d65da728c3aa555`, after the accepted F.2 bundle
recorded in `evidence/f2-fix1-runtime-closure.md`.

F.3 is a detached, aggregate-only relative-response map over the accepted
`hgcer3` basis. It consumes exactly five F.1 v2 contracts and the accepted F.2
representation, revalidates their authority/fingerprints, builds exactly 15
independent setting × canonical-t parents, and must reproduce F.2 support/OOD
summaries at absolute tolerance `1e-12`.

The model is a full-parent robust-scaled, class-balanced logistic fit with
analytic-gradient L-BFGS-B (`lambda=1e-3`, `maxiter=1000`, `gtol=1e-8`). Its
map is `exp(beta dot z)` without the fitted intercept, so the training median
defines relative response one. cKDTree p99 support masks review-only grids;
no grid, event probability, correction, weight, normalization, Method-B input,
or production application is persisted or constructed.

Local focused F.3 module/CLI tests, frozen F.2/F.1 tests, collector tests,
`py_compile`, profile JSON parsing, and `git diff --check` pass. Farm validation
ran only the detached F.3 analyzer and unchanged generic collector profile.
The accepted runtime bundle is recorded in evidence/f3-runtime-closure.md; F.4
may consume its frozen `hgcer3` map without retuning it.

F.3 can become available only if all 15 parents satisfy response support,
scaling, optimizer, finite-parameter, non-sparse support, and F.2 continuity
requirements. No pooling, fallback basis, threshold relaxation, or partial map
is allowed. The output is one deterministic JSON plus one seven-page matplotlib
review PDF. F.4 is DEVELOPMENT COMPLETE, FARM VALIDATION PENDING; it may
consume the accepted F.3 map without retuning it.
