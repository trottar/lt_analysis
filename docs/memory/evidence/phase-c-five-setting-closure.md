# Phase-C five-setting Method-B closure

## Status

CLOSED / RUNTIME VALIDATED for this specific Phase-C diagnostic closure.

## Evidence and evaluated source identity

[FARM_EVIDENCE] [HANDOFF]

- C.Fix.2.3 implementation:
  9a66bc62d20a99172e326e915866877b65ae1e5d.
- Later accepted pre-E.3 source:
  e3853655db0809923cbf2326e2f779219128eda9.
- Settings: Left lowe, Left highe, Center lowe, Center highe, Right highe.
  Right-low is not part of this gate.

## Recovered review

All 30 C.Fix.2 pages rendered correctly. The recorded scientific conclusion was
Phase C PASS. Same-canonical-t normalization converted raw Q with median about
3.687 and maximum about 277.83 to Qtilde with median about 1.013 and maximum
about 4.012.

Adaptive Method B was stable where it was measurable, with maximum chi2/ndf
1.5573, maximum absolute log pull 2.486, largest approximate low/high
difference about 1.29 sigma, and shared-cell cross-setting differences below
about 1.5 sigma. Its support was insufficient for a robust setting-wide role:

- 36 of 150 candidate-bearing cells (24%);
- 21 of 150 multi-slice cells (14%);
- 114 of 150 unavailable cells (76%);
- Left-high had 2 of 30 candidates, both single-slice.

## Conclusion and scope boundary

Adaptive Method B is DO NOT PROMOTE. Legacy Method B remains the Phase-D
reference. This evidence does not alter production pion subtraction, does not
authorize C.Fix.2.4 or D.Fix.1, and does not supply a final Phase-D farm
manifest.

Raw bundle filenames, hashes, and paths were not recovered by the import; the
reviewed handoff and page/result record are the imported farm provenance.

