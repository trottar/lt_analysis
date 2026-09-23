# E.8.1.Fix.5 - persisted-overlay ROOT-PDF geometry repair

## Status

SOURCE REVIEWED — this narrow source repair starts from committed test HEAD
fdd368f2084c9b9508e7ce679b8f7391f5b556f1. ChatGPT independently inspected
the complete refreshed actual diff and found the renderer/test repair PASS.
Codex-reported deterministic checks were NOT RUN by ChatGPT. No ROOT/PyROOT,
Jefferson Lab farm, or runtime acceptance is claimed for Fix.5, and the user
has not claimed a Fix.5 farm/runtime result.

## Pushed source identity

The user pushed this independently source-reviewed Fix.5 repair as
53fd262b730af8f1254e411a38231aebeb6a1da3. Independent ChatGPT pushed-state
source review passed for that exact source and its parent
fdd368f2084c9b9508e7ce679b8f7391f5b556f1: the reviewed square 3600 x 3600
overlay canvas, canonical three-by-nine grid, grid-before-header draw order,
legend semantics/colors, and pre-Print Modified()/Update() behavior remain
intact. That source review is distinct from ROOT/PyROOT, Jefferson Lab farm,
or runtime validation, none of which is claimed for Fix.5.

## Farm-owned blocker

The fresh Fix.4 Q4p4W2p74 / Left / lowe bundle passed provenance, frozen
F.6.2 identity, checker, page inventory, context page, handoff page, and all
six map-page gates. Its only presentation failure is the persisted-overlay
ROOT-PDF output on pages 38/41/44: the header and L/B/A legend are absent and
the first grid row is clipped. The direct bundle identities and review facts
are preserved in [the Fix.4 overlay blocker evidence](../evidence/e8-1-fix4-left-lowe-overlay-blocker.md).

## Repair boundary

Fix.5 changes only the E.8 persisted-overlay renderer and its focused test:

- use the square 3600 x 3600 overlay canvas while retaining the separate,
  non-overlapping header and canonical three-by-nine grid pads;
- render the grid first, then populate the header/legend as the top-level
  sibling immediately before PDF emission, and update the canvas before print;
- add a conditional real-ROOT plus pdftotext ordinary-PDF regression for
  the header, all L/B/A legend labels, and the first child title.

The persisted numerical values, L/B/A meanings/colors/proxies, first-child
title, all 27 tiles and their order, page IDs/order, reader authority,
context/handoff/map renderers, profile, collector, wrapper, accepted F.6.2
artifact, production physics, Method-A/Method-B boundaries, cuts,
normalization, templates, and binning are unchanged.

## Deterministic local validation

- python -B -m py_compile src/cuts/full_background_subtraction_plots.py testing/test_full_background_subtraction_plots.py — PASS.
- python -B -m unittest testing.test_full_background_subtraction_plots.FullBackgroundSubtractionE8Tests -v — PASS, 12 tests with 1 expected skip (PyROOT unavailable).
- python -B -m unittest testing.test_full_background_subtraction_plots -v — PASS, 102 tests with 18 expected skips.
- Memory manifest write/check, memory health, memory bootstrap, and git diff
  --check — PASS.

The conditional real-PDF regression did not run locally because PyROOT is
unavailable. These are Codex-reported source checks only, were NOT RUN by
ChatGPT, and do not replace fresh farm PDF inspection.

## Next action

NEXT — E.8.1.Fix.6 is the separate narrow profile provenance re-pin that must
require the pushed Fix.5 source 53fd262b730af8f1254e411a38231aebeb6a1da3.
It requires independent actual-diff review, user-controlled commit/push, and
ChatGPT pushed-state review before a fresh Q4p4W2p74 / Left / lowe farm/PDF
gate returns pages 37--47 for inspection. Only after persisted-overlay pages
38/41/44 visibly contain the header, L/B/A legend, first canonical child
title, and an unclipped first row may validation return to detailed Left /
highe; broader coverage remains gated on that result.
