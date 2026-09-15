# Phase F.5 Method-A signed `(t, phi)` propagation

CLOSED / RUNTIME VALIDATED — F.5 consumed the accepted F.4 correction, F.3
map, and five F.1 v2 contracts, then produced detached signed baseline and
F.4-adjusted canonical `3 x 9` `(t, phi)` aggregate templates. Its accepted
farm bundle ran at `3c6a66b7df9bf17e5a428458a2281a80831f001a`; 135 cells
(115 occupied, 20 explicit empty), F.1/F.4 continuity, and all parent/setting
closures passed. See evidence/f5-runtime-closure.md.

F.5 reuses F.4's public shared calculator and requires full F.4 reproduction
before it pairs transient correction factors to stable F.1 application-row
identities. Every occupied phi child must agree with F.4; an F.4-omitted child
is permitted only when F.1 has zero events, and is materialized as an explicit
zero cell. Parent-only normalization is frozen: no child renormalization,
clipping, smoothing, interpolation, rebinning, pooling, or new uncertainty
model exists.

The result is detached and non-authoritative. It persists aggregate matrices
and closure summaries only; it creates no event record, correction output,
ROOT histogram, Method-B input, yield, cross section, or production object.
F.5.2 is DEVELOPMENT COMPLETE, FARM VALIDATION PENDING as a presentation-only
PDF terminology rerender. It may not alter F.5 JSON, scientific fingerprint,
template data, authority, formulas, geometry, normalization, or any detached/
production boundary. Its farm gate is one rerender of accepted inputs and
visual review; the unchanged collector is optional for a new review bundle.
F.6 remains BLOCKED pending F.5.2 acceptance and a separate promotion decision.

SOURCE REVIEWED — F.5.Fix.1 begins from
`f9ce5aadbebc31f5163a093ab32244644ad2a088` after the first F.5 farm run found
that its Python version rejects `zip(strict=True)`. The repair removes only
that redundant keyword and is implemented/source-reviewed at
`46598878c102a67275b1600e73cba3b2f166dc26`; the preceding exact
row/factor-length, finite, and positive-factor gate remains unchanged and is
explicitly regression-tested. It changes no authority, ordering, geometry,
arithmetic, F.4 behavior, or detached/production boundary. The required next
gate is the same detached F.5 analyzer command followed by the unchanged F.5
generic collector.
