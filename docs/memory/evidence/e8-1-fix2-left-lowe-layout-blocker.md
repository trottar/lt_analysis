# E.8.1.Fix.2 Left/lowe procedure-PDF layout blocker

## Status

This is a direct farm-evidence blocker record. E.8.1 remains `ACTIVE`; this is
not E.8.1 acceptance or a claim about Fix.3 farm/runtime execution.

## Supplied bundle facts

- Bundle: `KaonLT_E8_1_Fix2_Q4p4W2p74_Left_lowe_20260923-040427.zip`.
- ZIP SHA-256:
  `52e0d283acdb57b56ed4c34c59f98dfe077f2a89e47fd0ced60fdcc91d5a5fd4`.
- Bundle `git_head`: `9f7094b19f3cdb0887c95f1f972d98471f2ecfdc`.
- Required analysis commit: `8985d9a212799c021c4ad1a759a689ea3826e0ea`;
  it is an ancestor and no unexpected committed files followed it.
- Manifest: `complete=true`, `errors=[]`.
- Frozen F.6.2 JSON SHA-256:
  `5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1`.
- Procedure-PDF SHA-256:
  `0745172508e2131e9256425b05bc44b19ffa6bcaa649224328679a75b7e18007`.
- Page-manifest SHA-256:
  `d4bc89f49ce374f63055afd974d734bf782cea6e55cefcc26caf3c2f83028dc1`.
- Renderer failures: `[]`. The page manifest has 47 pages: 36 retained D.6-D.9
  pages, E.8 context, three E.8 pages each for t1/t2/t3 (persisted overlays,
  delta-acceptance maps, missing-mass-acceptance maps), and E.8 handoff. No
  `full_background.e8.unavailable` page is present.
- Farm source module `testing.test_full_background_subtraction_plots`: 100
  tests, 14 skips, return code 0.

## Narrow successful gates

The available E.8 section establishes that the Fix.1 frozen persisted-setting
reader accepts the actual accepted F.6.2 artifact. The complete bundle,
ancestry/range checks, and admitted-file result establish the Fix.2 profile
provenance re-pin. Those narrow successes do not close E.8.1.

## Independent PDF review blocker

PDFium and Poppler review found context-page right-edge clipping, including
provenance and L/B/A text; ROOT em-dash mojibake on that context page; right-edge
clipping of the handoff D.10-omission sentence; and clipping of the overlay
top/header region. The sparse canonical overlay rows remain required and were
not implicated. E.8.1.Fix.3 owns only that presentation-layout repair.
