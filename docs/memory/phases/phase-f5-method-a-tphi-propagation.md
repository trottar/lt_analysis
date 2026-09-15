# Phase F.5 Method-A signed `(t, phi)` propagation

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING — F.5 begins from accepted F.4
farm/source HEAD `67e0298c51759c7a5ba693464d2c2655bf39250d`. It consumes only
the accepted F.4 correction artifact, accepted F.3 map, and five F.1 v2
contracts, then creates detached signed baseline and F.4-adjusted canonical
`3 x 9` `(t, phi)` aggregate templates.

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
Farm validation is one detached F.5 analyzer run plus the unchanged generic
collector using the F.5 profile. F.6 remains BLOCKED pending explicit F.5
runtime acceptance and a separate promotion decision.

SOURCE REVIEWED — F.5.Fix.1 begins from
`f9ce5aadbebc31f5163a093ab32244644ad2a088` after the first F.5 farm run found
that its Python version rejects `zip(strict=True)`. The repair removes only
that redundant keyword; the preceding exact row/factor-length, finite, and
positive-factor gate remains unchanged and is explicitly regression-tested.
It changes no authority, ordering, geometry, arithmetic, F.4 behavior, or
detached/production boundary. The required next gate is the same detached F.5
analyzer command followed by the unchanged F.5 generic collector.
