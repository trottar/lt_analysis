# Detached Method-A aggregate JSON inspection — 2026-10-11

## Scope and attribution

Read-only aggregate checkpoint for `Q4p4W2p74 / Left / lowe`, parents
`t0/t1/t2`, primary `t1`. This is not another runtime closure. RUNTIME VERIFIED
refers only to the [October 4 accepted farm diagnostic](e8-4-left-lowe-method-a-current-lineage-diagnostic-runtime-closure-2026-10-04.md).
The user-supplied [recording contract](../phases/e8-method-a-weight-findings-memory-and-reweighting-next-task-contract.md)
reports ChatGPT's independent recalculation of aggregate comparisons from that
accepted uploaded JSON. Codex checked those identities against the October 4
record and current field definitions; it did not rehash the raw artifact,
reproduce F.4 from original inputs, run the farm, inspect new ROOT output or
newly validate PDF rendering. Historical F.4 reproduction belongs to the
accepted artifact and prior review.

## Artifact identity and checks

JSON: `Q4p4W2p74_Left_lowe_kaon_pion-background_hgcer-method-a-current-lineage-diagnostic.json`.

| Identity | Value |
| --- | --- |
| Bytes | 690916 |
| File SHA-256 | `ecd547de15aa19a50b73bfd8d9dd1c159027c6438ca38b79bc45ff69a8d0e8c7` |
| Diagnostic fingerprint | `2c68804d60aaea423cf3dffcd55c3e324c812b97298149df8599ba193ee23616` |
| Artifact fingerprint | `7b6c1b244eb15d5d262391ca4fb288d896b9776989c9f658f2877b05dd02205e` |
| Artifact-reported git_head | `aad27a4d1639eef188dc61563fd3615682835835` |
| Candidate F.4 source SHA-256 | `1d545924eba89c7f9ffa28028e307aca9b434a89beec06863cf2893887b6b902` |
| Accepted F.4 correction fingerprint | `bce0cc12ef283b0c361818b436c4e190fc4ea5165436f6c36912cfbb6de53368` |

Prior ChatGPT inspection recomputed both canonical JSON fingerprints with
`sha256(json.dumps(core, sort_keys=True, separators=(',', ':'), allow_nan=False))`;
both matched. File size/hash matched accepted evidence. Recorded flags:
`accepted_authority_match=true`, `f4_exact_reproduction.payload_identical=true`,
`non_authoritative=true`, `production_objects_mutated=false`.

## Scientific definitions checked against current source

SOURCE VERIFIED: [F.3](../../../src/cuts/pion_hgcer_method_a_acceptance_map.py)
fits relative response `r` on accepted `hgcer3` coordinates from positive-response
training (weak-positive `0<NPE<=2` versus control `NPE>2`). It is not absolute
pion-to-kaon leakage probability. [F.4](../../../src/cuts/pion_hgcer_method_a_parent_preserving_correction.py)
uses each canonical-t signed application measure:

```text
b_j = s_j*w0_j; B = sum b_j; U = sum b_j*r_j
N_signed = U/B; C_j = r_j/N_signed
b_rw,j = s_j*w0_j*C_j; sum b_rw,j = B
```

`w0` owns baseline pion-control-to-background transfer. Fitted `r` is independent;
normalized `C` depends on `w0` and the signed source measure. There is no child
renormalization. The [F.6.3 producer](../../../src/binning/calculate_yield.py)
fills `B_pi_A`, forms `MM_A=pion_input-B_pi_A`, and extracts stored `YA` using
the existing error helper. `YA/MM_A` is already the reweighted-baseline final
branch; background `b/b_rw` is distinct from final kaon `Y0/YA`. Standalone
Method-A-only final yield/MM remains NOT VERIFIED. All formulas/selections stay
unchanged.

## Recalculated aggregate findings from the accepted JSON

The following are prior independent ChatGPT read-only aggregate comparisons,
recorded here from the contract, not a new farm validation.

| Quantity | t0 | t1 | t2 |
| --- | ---: | ---: | ---: |
| Application contributions | 11019 | 19809 | 24552 |
| Signed baseline parent B | 0.1081517597541286 | 0.3374726818838243 | 0.10645779207740595 |
| N_signed = U/B | 4.114018603391793 | 0.8866714623885792 | 2.0555330576308606 |
| Median raw r | 0.990491 | 0.729158 | 0.971664 |
| Median normalized C | 0.240760 | 0.822355 | 0.472706 |
| p99 normalized C | 6.326228 | 3.370358 | 4.326288 |
| max normalized C | 24.828515 | 60.147723 | 35.827893 |
| OOD count | 13 | 21 | 33 |
| OOD fraction | 0.0011797803793447414 | 0.001060124185976119 | 0.0013440860215053752 |
| Signed-vs-abs response normalization relative difference | +0.2700% | -1.0245% | +3.1298% |
| sum_phi(abs(DeltaB_phi))/(2*B) | 14.1207% | 8.5076% | 17.9891% |

`r/C` percentiles are rounded, unweighted marginals; JSON owns exact decimals.
The phi metric describes signed background redistribution, not fit uncertainty
or an extracted-yield change. The [diagnostic source](../../../src/cuts/pion_hgcer_method_a_current_lineage_diagnostic.py)
persists `baseline_signed_sum`, `adjusted_signed_sum`, `signed_delta`,
`absolute_baseline_support` and counts; fixed coordinate bins include flows.

| Phi interval (degrees) | t0 DeltaB/B | t1 DeltaB/B | t2 DeltaB/B |
| --- | ---: | ---: | ---: |
| [-180,-140) | -12.007% | +6.513% | +12.040% |
| [-140,-100) | -2.109% | +1.950% | +5.937% |
| [100,140) | +6.316% | -3.939% | -8.079% |
| [140,180) | +7.786% | -4.544% | -9.772% |

This is a representative subset of nine canonical children. The other five
are `[-100,-60)`, `[-60,-20)`, `[-20,20)`, `[20,60)`, `[60,100)`;
their numerical changes are not supplied by this recording contract and are
NOT VERIFIED by Codex. They are not absent or zero. Full-parent signed changes
sum to approximately zero (residuals around 10^-17). Prior aggregate recalculation
reported that phi, analysis_MM, SHMS_delta, P_hgcer_xAtCer and P_hgcer_yAtCer
baseline/adjusted sums, deltas and counts reconstruct parent values within
floating-point precision. Separate marginals do not establish joint `(phi,MM)`.

| Fixed diagnostic MM interval (GeV) | t0 DeltaB | t1 DeltaB | t2 DeltaB |
| --- | ---: | ---: | ---: |
| [0.905,0.935) | -0.03367030 | -0.00434208 | +0.00423250 |
| [1.100,1.145) | +0.02400379 | +0.00568501 | +0.00101008 |

For t0, the supplied ratios are approximately -71.21% and +165.29% of the
respective local signed baseline sum. Numerators above retain the supplied
precision. The denominators are sums of `baseline_signed_sum` over these
fixed bins; their exact decimals were not supplied in the contract and remain
NOT VERIFIED here. Do not reverse-engineer them from rounded ratios or use
these ratios as uncertainties. Small/cancelled denominators can amplify local
percent changes (INFERENCE). These are signed background contributions, not
final kaon yields or changed Lambda physics windows. Parent preservation does
not require MM-region or phi-cell stability.

Representative fixed-bin coordinate changes, as supplied:

| Parent / coordinate | Interval | DeltaB |
| --- | --- | ---: |
| t0 HGC x | [-8,-4) | +0.02095294 |
| t0 HGC x | [12,16) | -0.02546252 |
| t1 HGC y | [-8,-4) | +0.01236655 |
| t1 HGC y | [8,12) | -0.01392961 |
| t2 HGC x | [4,8) | +0.01312707 |
| t2 HGC x | [-16,-12) | -0.00926777 |

These describe acceptance-correlated redistribution, not a specific hardware
cause. Positive-response proxy evidence does not directly calibrate
`P(NPE=0 | true pion,x)`.

## Interpretation and limitations

Parent preservation coexists with nontrivial redistribution and finite long
`r/C` tails, which remain unclipped. Current t0/t2 cannot silently inherit
historical F.6.2 t1 continuity. The JSON persists marginal `w0/r/C/b/abs_b`
percentiles, but no joint event-level `w0*C` distribution, tail-contribution
decomposition based on `abs(s*w0*C)`, per-cell effective sample size, or joint
factor correlations with `(phi,MM)`.

NOT VERIFIED: top-tail C dominance of signed/absolute background, statistically
justified reweighted-contribution uncertainty, absolute true-pion NPE=0 mis-ID,
direct t0/t2 cause, new production readiness or canonical-five current-lineage
validation. Correction effect is not automatically systematic uncertainty.
Any later effective-statistics study must specify a nonnegative reference
measure (e.g. abs(b) or abs(b_rw)), moments and cancellation context; signed
weights are not probabilities or Poisson counts. No metric is approved here.

[CURRENT](../CURRENT.md) owns the sole scientific NEXT: read-only comparison of
actual w0 versus w0*C and b versus b_rw, emphasizing t0/t2 tails, source classes,
signed/absolute support, squared-weight/effective statistics and joint locality
where supported. Audit existing authority-checked sources first; only if
indispensable propose a separately contracted detached aggregate-only diagnostic.
No synthetic event reconstruction, new normalized yield branch, F.4 change,
owner retry or automatic farm work. Microscopic F.4/F.1 work stays DEFERRED;
the exact private-owner mismatch and absolute-SIMC blocker in CURRENT remain
binding. Method A remains detached/non-production; Method B remains independent,
diagnostic only and numerically excluded. This checkpoint creates no phase,
promotion or additional scientific/runtime acceptance.
