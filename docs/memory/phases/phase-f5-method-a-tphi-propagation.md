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
CLOSED / RUNTIME VALIDATED — F.5.2 is the accepted presentation-only farm
rerender at `6634e9cb470cf35f21f5d475ec6ce33b524cd233`, packaged as
`KaonLT_PhaseF5_2_Fix1_validation_Q4p4W2p74.zip` with SHA-256
`323d5d50604d846692c2b0bb0d1455c3751b2033601046dfcc3e07ccc3ca3ebd`.
Its JSON and PDF SHA-256 values are respectively
`143e3af6b1c69560e5bf351155b570d05f19e7e0ceb14665f5448c351ad501be` and
`e71bf819014b7410c8e403c54e5114460626484f21dab763323b13dac0111aa4`.
F.5.2 was presentation-only: the F.5 scientific propagation payload and its
scientific fingerprint
`d11b728d1089301a12c29e7f8b1798c6e0b6021ac47bd5d2afc2b62b47a1effa` were
exactly unchanged. The accepted 12-page review uses physics-readable terms and
physical phi intervals rather than internal phi indices. No scientific or
production behavior changed. F.6 is NEXT, pending a separate promotion
contract.

SOURCE REVIEWED — F.5.2.Fix.1 is the distinct presentation-label repair at
`6634e9cb470cf35f21f5d475ec6ce33b524cd233`; it is separate from F.5.Fix.1
below and does not replace the accepted F.5 scientific runtime identity.

SOURCE REVIEWED — F.5.Fix.1 begins from
`f9ce5aadbebc31f5163a093ab32244644ad2a088` after the first F.5 farm run found
that its Python version rejects `zip(strict=True)`. The repair removes only
that redundant keyword and is implemented/source-reviewed at
`46598878c102a67275b1600e73cba3b2f166dc26`; the preceding exact
row/factor-length, finite, and positive-factor gate remains unchanged and is
explicitly regression-tested. It changes no authority, ordering, geometry,
arithmetic, F.4 behavior, or detached/production boundary. The historical
repair was validated by the accepted F.5 farm bundle above. F.5.2 subsequently
closed as its distinct presentation-only rerender; neither repair altered F.5
science or production ownership.
