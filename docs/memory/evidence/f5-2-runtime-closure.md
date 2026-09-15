# F.5.2 accepted runtime closure

F.5.2 is `CLOSED / RUNTIME VALIDATED` from the accepted presentation-only
farm rerender. It did not change F.5 scientific or production behavior.

## Accepted bundle

- ZIP: `KaonLT_PhaseF5_2_Fix1_validation_Q4p4W2p74.zip`.
- ZIP SHA-256:
  `323d5d50604d846692c2b0bb0d1455c3751b2033601046dfcc3e07ccc3ca3ebd`.
- Farm/source HEAD: `6634e9cb470cf35f21f5d475ec6ce33b524cd233`.
- Manifest: `pion_hgcer_validation_bundle/v4`, profile
  `phase_f5_tphi_propagation_farm_review/v1`, `complete = true`,
  `errors = []`, and `unexpected committed files = []`.
- F.5.2 JSON SHA-256:
  `143e3af6b1c69560e5bf351155b570d05f19e7e0ceb14665f5448c351ad501be`.
- F.5.2 PDF SHA-256:
  `e71bf819014b7410c8e403c54e5114460626484f21dab763323b13dac0111aa4`.
- F.5 scientific fingerprint:
  `d11b728d1089301a12c29e7f8b1798c6e0b6021ac47bd5d2afc2b62b47a1effa`.
- F.5 artifact fingerprint:
  `261968ee7d9590d7d0afe0cd15155ef745a4169f95f52ffffd392aadd320de63`.

## Accepted result

The fresh propagation payload is exactly unchanged from accepted F.5. Only
`generated_at_utc` and `git_head` changed in JSON provenance. Event counts,
baseline and Method-A pion-background yields, signed changes, fractional
redistribution, t/phi geometry, parent and setting closure, authority, F.4
correction, scientific fingerprint, and artifact fingerprint remain unchanged.

The collector source checks passed and the bundle is complete. The review PDF
has 12 pages: one authority/physical-interpretation page, five
pion-background-yield-vs-phi pages, five `(t,phi)` physical-map pages, and one
fifteen-parent closure page. Visual review accepted physics-readable
presentation terminology and physical phi intervals such as `[-100, -60] deg`
in place of internal indices. No scientific change or production change
occurred.
