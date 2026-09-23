# E.8.1.Fix.4 Left/lowe persisted-overlay PDF blocker

## Status

This is direct fresh Jefferson Lab farm evidence. E.8.1 remains ACTIVE.
It closes E.8.1.Fix.4 only for its narrow bundle-profile provenance re-pin; it
does not accept E.8.1 or the persisted-overlay presentation.

## Bundle provenance and frozen input

- Bundle: KaonLT_E8_1_Fix4_Q4p4W2p74_Left_lowe_20260923-062807.zip.
- ZIP SHA-256:
  794c719d994444b49a58b0b9b326d5de019305a2e876615bf989f8c11bc01d59.
- Bundle schema: pion_hgcer_validation_bundle/v4; profile:
  phase_e8_1_full_background_procedure_pdf_farm_review/v1.
- Bundle git_head: fdd368f2084c9b9508e7ce679b8f7391f5b556f1.
- Required analysis/procedure source:
  350c34c55b2de33ad01011559dc6d8ed84d9c8a7; it is an ancestor, and
  there are no unexpected committed files.
- Manifest: complete=true, errors=[]; captured source worktree is clean.
- Frozen F.6.2 JSON SHA-256:
  5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1.
- Frozen artifact fingerprint:
  ee713b70de898bad8fa61164cbc8a54886af1ea1eb1e712ed95de17df8890cd0.
- Frozen validation fingerprint:
  7edc73fce20ad7dc7622b8c23a7ba7e8986e367ace605884e010945595370b3b.

## Passed gates

- Fresh procedure PDF SHA-256:
  9a00fc6179891e936523613293dc40ebe5790f7ac49acf37aa942aea0333b98c.
- Fresh page-manifest SHA-256:
  d4bc89f49ce374f63055afd974d734bf782cea6e55cefcc26caf3c2f83028dc1.
  It records 47 pages with renderer failures [].
- Provenance, frozen-artifact, checker, page-inventory, context page 37,
  handoff page 47, and map pages 39/40/42/43/45/46 pass review. The context
  and handoff text are readable; the map renderers require no redesign.
- Farm source checks passed: py_compile;
  testing.test_full_background_subtraction_plots (101 tests, 14 skips);
  Method-A acceptance contract (7); Phase-F runtime contract (3); collector
  (28); git diff --check; and required-range diff.

## Persisted-overlay blocker

Independent PDFium and Poppler review found the same defect on persisted-overlay
pages 38, 41, and 44: the E.8 persisted L/B/A overlays - tN header and
three-column L/B/A legend are absent, the first plot row is clipped at the
physical PDF top, and pdftotext does not contain those header/legend lines.
On page 38, Left-lowe phi0 [-180, -140) and Left-lowe phi1 [-140, -100) are
absent; later titles, including Left-lowe phi2 [-100, -60), are extractable.
This is a real farm PDF-output failure, not a viewer artifact.

E.8.1.Fix.5 owns only the persisted-overlay ROOT-PDF geometry repair. The
accepted F.6.2 artifact, frozen reader authority, persisted values, canonical
page IDs/order, context/handoff pages, map pages, profile, collector, wrapper,
and production physics remain outside that repair.
