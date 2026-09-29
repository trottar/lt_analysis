# E.8.4 pion-alignment determinism blocker — Left / lowe

## Status

Direct Jefferson Lab runtime evidence of a determinism blocker; this record is
not runtime closure or production acceptance.

## Artifact and setting

- Artifact: `KaonLT_E8_4_Left_lowe_alignment_20260929-102321.json`.
- Setting: `Q4p4W2p74 / Left / lowe`.
- The compared September 24 and September 29 setting-wide pion-control records
  have identical content checksum
  `0dc9c642ade833a0b8fd727be1ba55dd79a09392082c07e12afafa82cbdb0ab1` and
  identical axes, but different generated ROOT histogram names.
- The September 29 store reported `persistence_status="rejected_stale_then_created"`
  with `pion control histogram identifier mismatch`, rather than reusing the
  semantically identical alignment.

## Coupled candidate instability

The active configuration records `minimum_template_integral=1.0`,
`renormalize_shifted_templates=true`, and `interpolation_mode=linear`.

- September 24 accepted the `pi_n` residual-shift candidate `-0.003 GeV` and
  rejected its `-0.004 GeV` neighbor for insufficient template integral.
- September 29 made the reciprocal choice: it rejected `-0.003 GeV` for
  insufficient template integral and accepted `-0.004 GeV`.
- Because pion components resolve sequentially against the residual, the
  downstream parent `pi_delta` solution changed from `+0.008 GeV` to
  `+0.016 GeV`.
- The parent baseline scores were effectively identical
  (`73.77888464083372` and `73.77888471147615`), and the common setting shift
  differed by only about `4.7e-12 GeV`.

## Conclusion

The observation establishes a cache-provenance and machine-roundoff
determinism blocker in existing baseline pion alignment, not a changed
pion-control population or a scientific-model result. It requires the narrow
E.8.4.Fix.3 repair and independent source review before the blocked later
F.6.3/E.8.4 runtime gate can resume. No artifact checksum is asserted here
because none was measured for this record.
