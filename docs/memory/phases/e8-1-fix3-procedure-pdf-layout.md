# E.8.1.Fix.3 - procedure-PDF layout repair

## Status

`SOURCE REVIEWED` — this narrow presentation-only repair starts from committed
`test` HEAD `9f7094b19f3cdb0887c95f1f972d98471f2ecfdc`. ChatGPT independently
inspected the complete actual diff and found the source/test repair PASS.
Codex-reported deterministic checks were `NOT RUN by ChatGPT`. No ROOT/PyROOT,
farm, or runtime validation is claimed, no Fix.3 commit exists yet, and this
record makes no E.8.1 acceptance claim.

## Farm-owned blocker

The complete targeted `Q4p4W2p74 / Left / lowe` Fix.2 bundle proved the
persisted-setting reader and profile-provenance gates, but independent PDF
review found right-edge clipping on the E.8 context and handoff pages, ROOT
em-dash mojibake on the context page, and overlay top/header clipping. The
direct facts are in [the blocker evidence](../evidence/e8-1-fix2-left-lowe-layout-blocker.md).

## Repair boundary

Fix.3 changes only E.8 final presentation layout: explicitly bounded ASCII-safe
context/handoff text, and a reserved overlay header/legend pad above the
complete three-by-nine canonical grid. It retains all full SHA-256 and
fingerprint values, L/B/A meaning/colors/proxies, persisted values, sparse
children, page IDs, page order, and map pages.

The accepted F.6.2 JSON and its numerical/scientific content remain frozen.
Method A remains detached/non-production, Method B remains diagnostic only, and
no production weight, yield, cut, normalization, template, binning, efficiency,
cross section, or uncertainty changes.

## Deterministic local validation

- `python -B -m py_compile src/cuts/full_background_subtraction_plots.py testing/test_full_background_subtraction_plots.py` — PASS.
- `python -B -m unittest testing.test_full_background_subtraction_plots.FullBackgroundSubtractionE8Tests -v` — PASS, 11 tests.
- `python -B -m unittest testing.test_full_background_subtraction_plots -v` —
  PASS, 101 tests with 17 expected skips; PyROOT is unavailable locally.
- Manifest write/check, memory health, and memory bootstrap — PASS.

## Scope and next action

The intended change set is only `src/cuts/full_background_subtraction_plots.py`,
`testing/test_full_background_subtraction_plots.py`, this phase record, the
Fix.1/Fix.2 reconciliation records, `CURRENT.md`, the supplied blocker evidence,
the manifest, and the user-created Fix.3 task contract.

Next: user-controlled commit/push of this `SOURCE REVIEWED` repair, followed by
ChatGPT pushed-state review. A separate narrow profile re-pin must then target
that pushed Fix.3 source and receive review before a fresh `Left / lowe` PDF
gate. Detailed `Left / highe` and broader coverage remain gated on independent
visual passage of the repaired Left/lowe PDF.
