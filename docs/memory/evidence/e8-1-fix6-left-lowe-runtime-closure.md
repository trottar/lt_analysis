# E.8.1.Fix.5/Fix.6 Left/lowe runtime closure

## Scope and status

`CLOSED / RUNTIME VALIDATED` — fresh Jefferson Lab evidence closes E.8.1.Fix.5 only for the persisted-overlay ROOT-PDF geometry repair and E.8.1.Fix.6 only for its narrow validation-profile/provenance gate. It does not close E.8.1 as a canonical-five-setting program, alter production physics, promote Method A, or replace frozen F.6.2 science.

The remaining E.8.1 canonical-five expansion is `DEFERRED` by explicit user decision. It is neither failed nor accepted by this single `Q4p4W2p74 / Left / lowe` gate.

## Accepted package and provenance

- Package: `KaonLT_E8_1_Fix6_Q4p4W2p74_Left_lowe_20260923-090850.zip`.
- Bundle/profile commit: `0ec29d4e1bb345eb37e8cca35b8b7e5cbe1b4d5b`.
- Required analysis/procedure source: `53fd262b730af8f1254e411a38231aebeb6a1da3`.
- Required-source ancestor check: passed; unexpected committed files after that source: none.
- Captured source worktree: clean. Bundle manifest: `complete=true`, `errors=[]`.

## Frozen F.6.2 authority

- Accepted JSON SHA-256: `5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1`.
- Accepted artifact fingerprint: `ee713b70de898bad8fa61164cbc8a54886af1ea1eb1e712ed95de17df8890cd0`.
- Accepted validation fingerprint: `7edc73fce20ad7dc7622b8c23a7ba7e8986e367ace605884e010945595370b3b`.

All three identities were preserved by the fresh package.

## Procedure PDF, tests, and visual review

- Procedure PDF SHA-256: `d809a80fb1a46602dbb6cb930ed1d6c8932c2583ce8c5158ac82a78e6e6301e5`.
- Page-manifest SHA-256: `d4bc89f49ce374f63055afd974d734bf782cea6e55cefcc26caf3c2f83028dc1`.
- Page count: 47; `renderer_failures=[]`.
- Farm checks: `testing.test_full_background_subtraction_plots` PASS, 102 tests with 14 skips; the real E.8 overlay PDF regression produced `e8-overlay-pdf-regression.pdf`; Method-A acceptance contract PASS, 7 tests; Phase-F runtime contract PASS, 3 tests; generic collector PASS, 28 tests; `git diff --check` PASS.

ChatGPT independently inspected pages 37--47. Context page 37, persisted-overlay pages 38/41/44, all map pages 39/40/42/43/45/46, and handoff page 47 passed. The overlay pages visibly contain the E.8 header, all L/B/A legend entries, complete first plot row, no top clipping, and the first canonical child title.

Independent text extraction confirmed:

- `L: upstream 0 < NPE <= 2 diagnostic reference`
- `B: physical pion control, NPE > 2, baseline w0`
- `A: same B population, w0*C`
- `Left-lowe phi0 [-180, -140)`

## Boundary and successor

The accepted Fix.5 source remains distinct from the Fix.6 bundle/profile source. This evidence supports neither a broader E.8.1 setting claim nor final E.8 closure. The next approved work is the presentation-only E.8.2 baseline full-analysis stage audit; the full dependency chain is owned by [the E.8 procedure roadmap](../decisions/e8-full-analysis-procedure-roadmap.md).
