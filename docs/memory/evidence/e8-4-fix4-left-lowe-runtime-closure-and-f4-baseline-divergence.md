# E.8.4.Fix.4 Left/lowe runtime closure and baseline divergence

## Direct farm evidence supplied for the fresh gate

The current `Q4p4W2p74 / Left / lowe` setting-wide persisted alignment record has `alignment_semantics_version = pion_component_dynamic_alignment_semantics/v2` and `persistence_status = rejected_stale_then_created`. Its stale-record rejection reasons include `alignment_semantics_version mismatch`. This directly validates the narrow E.8.4.Fix.4 semantic-version cache rejection and recomputation at runtime. It does not validate F.6.3 or E.8.4.

The current Left/lowe F.1 raw file SHA-256 is `10086d7d4c42389c9fdd16471c49d59cd189a30980767bcdb87b85914169ef95`. A streaming frozen/current comparison found:

| Projection | Records | Identity mismatches | Projection mismatches |
| --- | ---: | ---: | ---: |
| F.3 training inputs | 52,397 | 0 | 0 |
| F.3 application inputs | 55,380 | 0 | 0 |
| F.4/F.6.3 baseline numerics | 55,380 | 0 | 55,380 |

The first baseline mismatch is identity `["prompt", 0]`:

| Field | Frozen | Current |
| --- | ---: | ---: |
| `analysis_MM` | 0.9150867495356926 | 0.9150867495404144 |
| `baseline_pion_weight_w0` | 0.0434198336219789 | 0.04235726703122994 |
| `signed_baseline_event_contribution` | 8.049822627351833e-06 | 7.852828031293563e-06 |

The fresh procedure E.8.4 gate remained unavailable with `f6_3_f4_shared_reproduction_failed:f3_fingerprint_input_content_mismatch`.

## Source interpretation and next evidence gate

F.3 numerically consumes sanitized F.1 acceptance-coordinate training and application projections; the direct comparison shows those Left/lowe projections unchanged. F.4 parent normalization consumes `signed_baseline_event_contribution` as `b` in `B = sum(b)`, `U = sum(b * raw_factor)`, and `U/B`. Every Left/lowe application baseline contribution changed. Existing broad provenance checks must remain fail closed.

The detached F.4.Refresh.1 comparator must rebuild current candidate F.2/F.3/F.4 and identify the first scientific payload difference. These observations alone do not decide which accepted stage needs refresh. Historical F.1 through F.6.2 closures remain intact; F.6.3 and E.8.4 have no runtime acceptance.
