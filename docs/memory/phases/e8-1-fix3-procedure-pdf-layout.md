# E.8.1.Fix.3 - procedure-PDF layout repair

## Status

SOURCE REVIEWED — this narrow presentation-only repair started from committed
test HEAD 9f7094b19f3cdb0887c95f1f972d98471f2ecfdc. ChatGPT independently
inspected the complete actual diff and found the source/test repair PASS.
Codex-reported deterministic checks were NOT RUN by ChatGPT. The user pushed
that reviewed source as 350c34c55b2de33ad01011559dc6d8ed84d9c8a7.

The fresh Fix.4 Left/lowe farm bundle directly confirms that Fix.3 repaired the
context and handoff pages; the context/handoff text and all six map pages pass
independent review. It does not close the persisted-overlay visual subrepair:
pages 38/41/44 remain blocked. This status stays SOURCE REVIEWED and makes no
E.8.1 acceptance claim.

## Farm-owned blocker and preserved boundary

The original Fix.2 farm evidence found context/handoff clipping, context
em-dash mojibake, and overlay clipping. The later Fix.4 evidence confirms the
context/handoff repairs while isolating the remaining defect to the
persisted-overlay ROOT-PDF geometry. The direct current facts are in [the
Fix.4 overlay blocker evidence](../evidence/e8-1-fix4-left-lowe-overlay-blocker.md).

Fix.3 preserves all full SHA-256 and fingerprint values, L/B/A
meaning/colors/proxies, persisted values, sparse children, page IDs/order, and
map pages. The accepted F.6.2 JSON and its numerical/scientific content remain
frozen. Method A remains detached/non-production, Method B remains diagnostic
only, and no production weight, yield, cut, normalization, template, binning,
efficiency, cross section, or uncertainty changes.

## Historical deterministic local validation

- python -B -m py_compile src/cuts/full_background_subtraction_plots.py testing/test_full_background_subtraction_plots.py — PASS.
- python -B -m unittest testing.test_full_background_subtraction_plots.FullBackgroundSubtractionE8Tests -v — PASS, 11 tests.
- python -B -m unittest testing.test_full_background_subtraction_plots -v —
  PASS, 101 tests with 17 expected skips; PyROOT was unavailable locally.
- Manifest write/check, memory health, and memory bootstrap — PASS.

These historical source checks remain distinct from the later direct farm
evidence and from E.8.1 closure.

## Current successor

E.8.1.Fix.5 is ACTIVE for the sole remaining persisted-overlay geometry defect.
It must receive independent actual-diff review, then user-controlled
commit/push and pushed-source review. A separate profile re-pin to its future
pushed source precedes a fresh Left/lowe gate; detailed Left / highe and broader
coverage remain gated on visual passage of the repaired overlay pages.
