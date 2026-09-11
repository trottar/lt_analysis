# E.3.Fix.2 Q4p4W2p74 Left-low farm gate

## Status

CLOSED / RUNTIME VALIDATED for this exact setting and presentation gate only.

## Evidence and identities

[FARM_EVIDENCE]

- Required implementation:
  eb1710f4739ba6ef14f51419806e9fc5bd53c175.
- Runtime/bundle head:
  bf53dac84e1396cfd7e3f4e0234749426bfbdcf4.
- Gate: Q4p4W2p74 Left-low.

The recovered bundle contained manifest.json, source_state.txt,
source_checks.txt, a Left_lowe checkpoint JSON, and E.3 validation PDF. The
source procedure PDF had 61 pages; the collector extracted pages 59-61 as E.3
t1, t2, and t3.

## Recorded inspection

Source checks were required-analysis ancestry PASS, py_compile PASS,
full-background tests 61 OK, HGCer-refinement tests 63 OK, collector tests 16
OK, and git diff --check PASS. Checkpoint provenance spot checks matched.

Visual review passed complete 5x2 grids with all ten delta cells per page, no
clipping/header overlap, visible NPE=2 and low/control labels, retained negative
signed content, explicit unavailable cells, and no Method-B or correction
language.

## Boundary

This record validates independent Method-A presentation for this one setting.
It does not establish the remaining Left-high, Center-low, Center-high, or
Right-high E.3 closures. Their runtime status remains DEFERRED.

