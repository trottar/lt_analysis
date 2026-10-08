# Detached F.1 forensic reproducibility diagnostic review — 2026-10-08

## Authority and diagnostic-only scope

This receipt records user-supplied files and prior independent ChatGPT evidence
review under the [yield-first contract](../phases/e8-yield-first-memory-consistency-task-contract.md).
Codex did not locally access, hash or rerun these external raw artifacts.
The already-pushed CLI source is
`6cdd2fcaf60d0ee8fe8974dc9af8d59a625534b5`; its
[implementation record](../phases/e8-4-f1-baseline-reproducibility-diagnostic.md)
owns source/local-check history. SOURCE VERIFIED: the forensic parser is
hard-pinned to reviewed F.4 and never grants canonical equivalence acceptance.

## Supplied identities and observations

All identities/numbers here are attributed to the supplied files/prior review,
not local measurements or a fabricated farm receipt.

| Diagnostic/provenance input | SHA-256 |
| --- | --- |
| `E8_4_F1_baseline_reproducibility_Left_lowe_6cdd2fca.json` (63,467 bytes) | `a11b4ddd96a020a0cecff0961cc67a2c52470e2b8515518e88080c303a78e71c` |
| Current Left-lowe F.1 | `b5d0b10cb888b87915b3b4090e14d2c03436f93616fd5c42e79c729ecc0748e9` |
| Reviewed F.4 (supplied raw hash independently matched in prior review) | `79e7ceda7221cbeeead4ed5bc306b0e0e670741a27beaa980e22349c555e96d7` |
| Historical materialization manifest | `e7c37f55b24739e8be0fb778a6c3491e17df746263e6d9d71a4f31ec65cd7908` |
| Comparison input | `4c3fdab5d05965a1b3f1dc835c9d8c529896e6eaf2da8eb5c8ac8fc26e95257a` |

Diagnostic flags: `non_authoritative=true`,
`failed_owner_diagnostic_only=true`, `runtime_acceptance_granted=false`.
Reference F.1 is absent; original reviewed raw F.1 was not recovered.

| Canonical t index | Current/reviewed F.4-consumed count | Signed parent sum exact equal | Absolute parent sum exact equal | Current minus reviewed absolute sum |
| --- | ---: | --- | --- | ---: |
| 0 | 11019 / 11019 | false | false | +2.6627020033309634e-11 |
| 1 | 19809 / 19809 | false | false | -2.4133234299839046e-10 |
| 2 | 24552 / 24552 | false | false | -8.55021609069695e-12 |

The counts total 55,380 in each input; all six signed/absolute comparisons
are not exactly equal. Relative differences are at the scale of `10^-10`;
cause is NOT VERIFIED. Per-source/phi signed aggregates show broadly
distributed, primarily prompt-sourced differences, not an isolated phi-bin
count change. This does not identify coefficient, fitter, ROOT, template or
missing-mass lookup drift.

## Failed invariant, limits and user priority

The newer canonical-five owner at
`ebf3a0be2f50c788c22defe3601fd5b7fee48644` remains `BLOCKED` before
analysis at
`f6_3_scientific_equivalence_mismatch:f4:$.parents[0].absolute_baseline_parent_sum`.
The [failed-comparison contract](../phases/e8-4-f1-baseline-reproducibility-diagnostic-task-contract.md)
owns the supplied F.2/F.3 exact PASS and F.4 exact FAIL, with t0 absolute
reviewed/current sums `0.12072578658633845 / 0.12072578661296547`.
Small magnitude grants no negligible-roundoff exemption or runtime acceptance.
Historical second-full-run equality has a separate scope and does not resolve
this failure. No complete farm run or demonstrated physics defect follows.

Original reviewed F.1 and upstream full fit/template/weight-array/bin-edge
evidence are unavailable; event-level cause and full archival availability
remain NOT VERIFIED. Further F.1 archive/numerical archaeology is `DEFERRED`
by user priority unless concrete evidence makes it material or unavoidable.
Full-procedure proper extracted-yield readiness comes first; CURRENT alone
owns the substantive audit gate. This receipt completes no yield-readiness
milestone and authorizes no full-owner retry, bypass, tolerance change or farm
command. Final E.8/F.6.4 remain `BLOCKED`; accepted narrow Left/lowe closures
and pre-analysis-only preflight remain intact. Absolute-SIMC interpretation
stays separately `BLOCKED` at
`SIMC_normfac_luminosity_and_charge_units_not_source_proven`.
Method A remains detached/non-production and parent preserving; Method B
remains diagnostic only/numerically excluded. Physics and production are frozen.
