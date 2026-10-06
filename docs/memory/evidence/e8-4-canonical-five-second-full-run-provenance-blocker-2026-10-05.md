# Canonical-five second full-run provenance blocker

## Scope and recording authority

The second Q4p4W2p74 canonical-five full runtime gate is `BLOCKED`.
This record carries the supplied/inspected farm diagnosis and read-only
post-run comparator output specified by the
[closure/reconciliation contract](../phases/e8-science-first-reanchor-after-canonical-five-provenance-blocker-task-contract.md).
Codex records those facts; it did not rerun the farm or independently inspect
the raw farm artifacts in this checkpoint. Failed-run artifacts are diagnostic
evidence only, inadmissible for canonical-five runtime or visual closure.

## RUNTIME VERIFIED — supplied second owner attempt

Farm-evaluated source: `ace8688a71431d13b40ed19713a27746f3da6a8e`.

Gate-status basename:
`KaonLT_E8_4_Fix5_CanonicalFive_Q4p4W2p74_20261005-202554_1227600-gate-status.json`

Gate-status SHA-256:
`ddfadaf9946e262c35795689ab8d2d52a306ede9e775d98fb05df4e15924e91f`.

```text
status = failed
stage = verify_artifacts
failure_reason = new_page_missing_or_duplicate
analysis_started = true
analysis_returncode = 0
analysis_completed = true
artifact_verification_completed = false
collection_completed = false
zip_verification_completed = false
```

The normal owner lineage preflight passed before analysis, including shared
F.4 reproduction for all five settings. Analysis-child success is distinct
from owner success: the later artifact/page gate failed, and neither collection
nor ZIP verification completed. This does not establish canonical-five closure.

## RUNTIME VERIFIED — fresh page-manifest diagnosis

| Canonical setting | Page count | renderer_failures | E.8.4 inventory |
| --- | ---: | --- | --- |
| Left-lowe | 74 | `[]` | exactly one `full_background.e8_4.unavailable` |
| Left-highe | 74 | `[]` | exactly one `full_background.e8_4.unavailable` |
| Center-lowe | 74 | `[]` | exactly one `full_background.e8_4.unavailable` |
| Center-highe | 74 | `[]` | exactly one `full_background.e8_4.unavailable` |
| Right-highe | 74 | `[]` | exactly one `full_background.e8_4.unavailable` |

All five unavailable reasons were identical:

`f6_3_f4_shared_reproduction_failed:f3_fingerprint_input_content_mismatch`

INFERENCE: the generic owner failure `new_page_missing_or_duplicate` is the
downstream page-inventory consequence of F.6.3 becoming unavailable during
real analysis. It is not evidence of a duplicate-page renderer defect.

## RUNTIME VERIFIED — post-run F.1 identities

Full analysis rewrote all five F.1 artifacts after successful pre-analysis
lineage reproduction. These raw SHA-256 values and stable content fingerprints
differ from the lineage that passed the
[accepted separate preflight](e8-4-canonical-five-lineage-preflight-runtime-closure-2026-10-05.md).
Identity change alone does not establish a Method-A scientific payload change.

| Setting | Post-run raw SHA-256 | Post-run stable F.1 content fingerprint |
| --- | --- | --- |
| Left-lowe | `34e98696fab7807e9c37580d145b762fdfee3e53025c32dba586dabbef6d4f36` | `69c1b044abb143a187dd4471efba4ca430b61122204c41d6ad4225ba70d152ba` |
| Left-highe | `3648346330ab01cd52e7d58fb3002d12f758e4f9f14634017ec1e5a02e5b4a9c` | `e492e9ffb50b868855fd548ab99a97ff366ddb1a06fbac81a8be1535c3f8c3fd` |
| Center-lowe | `4fd5cc5c691d2f181ffa5134fcf803328e824b7bb94ae04fb8e91a7c8144c724` | `3e3fb1be62074a4e44c408b23674493a2c90f26ff93f6e1cb45bb7f5862ff7cc` |
| Center-highe | `fc5f340db724a0bd4dd9b2e527dd67003b088a8751bc53cce71956f28652c251` | `4cd7924d2a136e89ca508d3b4d093263fdace20b9cb62a7f540471152a4312e7` |
| Right-highe | `81a09c6dca80cbe5c5dc477ebcece19f297aadced33c7cb30cbf8492e3d7775a` | `2b41d1c72a2b0f3b62f4d4a64de0cc7731e465e3c405d085952d067a2ba6fdb8` |

## RUNTIME VERIFIED — tracked post-run comparator

The existing tracked read-only comparator
`testing/compare_method_a_current_baseline_authority.py` was run on the five
post-run F.1 artifacts against the reviewed current-baseline materialization.
Supplied farm output established:

```text
f2_scientific_payload_match = true
f3_scientific_payload_match = true
f4_scientific_payload_match = true
first_changed_stage = none
```

The full run regenerated a different F.1 raw/stable identity set, while the
tracked detached comparator reconstructed matching F.2/F.3/F.4 scientific
payloads. The [first-run comparison](e8-4-canonical-five-fresh-f1-lineage-materialization-2026-10-05.md)
had F.4 as its first changed stage against the preceding lineage; that distinct
historical result must not be attributed to this second-run comparison.

## INFERENCE — interpretation and current decision

The immediate failure is a provenance/identity coupling in the validation/
runtime path rather than evidence, from this gate, of a newly changed Method-A
F.2/F.3/F.4 scientific payload. No new Method-A scientific defect is established.

Canonical-five provenance/identity repair is `DEFERRED` by current user decision
to stop another minor-repair loop and return to the E.8 science objective.
The blocker is neither repaired nor irrelevant to later final closure.

## NOT VERIFIED — explicit non-closures

- Exact low-level cause of the run-to-run F.1 identity change.
- Whether future repair should use scientific-equivalence matching, producer
  canonicalization or another reviewed design.
- Final canonical-five runtime closure and canonical-five PDF/visual closure.
- Production promotion or production correctness.

The separate lineage preflight remains `CLOSED / RUNTIME VALIDATED` only for
its pre-analysis scope. Accepted Left/lowe branch, Fix.5.7/Fix.5.8 and earlier
F-stage closures retain their scopes without downgrade or reopening. E.8 is
`ACTIVE`; final E.8 and F.6.4 remain `BLOCKED`. Absolute-SIMC interpretation
remains separately `BLOCKED` by
`SIMC_normfac_luminosity_and_charge_units_not_source_proven`.

Method A remains detached/non-production. Method B remains independent,
diagnostic/cross-check only and numerically excluded. The baseline production
sequence and `no_empirical_residual` profile remain frozen. No pin refresh,
executable fallback, scientific redesign or farm command is authorized by this
checkpoint. CURRENT alone owns the live objective and sole ordinary NEXT.
