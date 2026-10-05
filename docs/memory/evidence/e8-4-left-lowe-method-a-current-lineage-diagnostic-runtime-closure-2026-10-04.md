# E.8.4 Left/lowe current-lineage Method-A diagnostic closure — 2026-10-04

## Scope and evidence attribution

`CLOSED / RUNTIME VALIDATED` only for the detached aggregate diagnostic at
`Q4p4W2p74 / Left / lowe`, parents `t0`, `t1`, `t2`, primary target `t1`.
The user supplied the farm JSON/PDF to ChatGPT, which independently reviewed
and accepted their provenance, payloads, pages and scientific interpretation
before this Codex closure task. Codex records that accepted evidence under the
[corrected closure contract](../phases/e8-4-left-lowe-method-a-current-lineage-diagnostic-runtime-closure-task-contract.md);
it did not rerun the farm or independently locate, hash or revalidate the raw
artifacts. The contract is recording authority, not runtime evidence. The
runtime evidence is the accepted artifact pair identified below. No separate
external acceptance-review file is required.

## Accepted artifacts and provenance

- JSON: `Q4p4W2p74_Left_lowe_kaon_pion-background_hgcer-method-a-current-lineage-diagnostic.json`,
  690916 bytes; SHA-256
  `ecd547de15aa19a50b73bfd8d9dd1c159027c6438ca38b79bc45ff69a8d0e8c7`.
- PDF: `Q4p4W2p74_Left_lowe_kaon_pion-background_hgcer-method-a-current-lineage-diagnostic.pdf`,
  71818 bytes; SHA-256
  `0f3e0d016ce895c6ff1aa390f28c3f3e85c99e2cc8733c7fa0a176aec48f0a73`.
- JSON `git_head`: `aad27a4d1639eef188dc61563fd3615682835835`.
- JSON `runtime_claim`: `detached diagnostic only; no production promotion`.

The JSON recorded a dirty farm checkout, not a globally clean one:

```text
 M src/models/xmodel_kaon_pl.f
?? src/kaon/functions/Q4p4W2p74.model
```

Those unrelated paths were not consumed by the detached JSON diagnostic.
Acceptance rests on its exact F.1/F.3/F.4 input authority and hashes. No new ZIP
or bundle identity is inferred.

| Accepted identity | Exact value |
| --- | --- |
| Candidate validation source head | `b349967c0d4210a78b144ce6134d3c1f15970245` |
| Left/lowe F.1 SHA-256 | `10086d7d4c42389c9fdd16471c49d59cd189a30980767bcdb87b85914169ef95` |
| Candidate F.3 source SHA-256 | `eeb480d9b7ddffab1e97c7f9f27c3099d0780c2bbb4d42a0ee68c077a5752b3d` |
| Candidate F.3 map fingerprint | `3a9787fc58d26cc0816012bd1b637ad0c8f201b448d54cb1625a841131154728` |
| Candidate F.4 source SHA-256 | `1d545924eba89c7f9ffa28028e307aca9b434a89beec06863cf2893887b6b902` |
| Candidate F.4 correction fingerprint | `bce0cc12ef283b0c361818b436c4e190fc4ea5165436f6c36912cfbb6de53368` |
| Candidate F.4 artifact fingerprint | `4c935271a0b2723b58b02cc34d82a1e7d5757f7cd36b7707f400c2b28893668a` |
| Diagnostic fingerprint | `2c68804d60aaea423cf3dffcd55c3e324c812b97298149df8599ba193ee23616` |
| Diagnostic artifact fingerprint | `7b6c1b244eb15d5d262391ca4fb288d896b9776989c9f658f2877b05dd02205e` |

`accepted_authority_match = true`; shared F.4 reproduction reported
`payload_identical = true` with the exact correction fingerprint above.

## Accepted JSON/PDF gates

RUNTIME VERIFIED (prior independent ChatGPT review): JSON bytes matched SHA-256;
schema/version, `non_authoritative = true`, recomputed artifact and diagnostic
fingerprints, exact scope, parent inventory and primary target all passed.
Exact current-lineage F.1/F.3/F.4 authority, accepted-authority matching, shared
F.4 reproduction and persisted payload identity passed. The finite aggregate
payload passed source-decomposition closure, child/MM/coordinate aggregate
reconstruction and signed canonical-parent closure. There was no Method-B
numerical dependency, production-object mutation, alternative correction
construction or forbidden persisted per-event factor/identity table.

RUNTIME VERIFIED: PDF bytes matched SHA-256; exactly eight pages were independently
rendered/reviewed and legible, with no missing page, clipping, overlap or broken
plot/glyph. Their roles were: authority/scope/limitations; weak-positive/control
response; `hgcer3` topology; signed/absolute parent normalization; t1 source
decomposition/tails; t1 child/MM sensitivity; t1 coordinate sensitivity; and
interpretation boundary. Smaller weak-positive curves reflect lower statistics,
not a rendering failure.

## Measurement A — relative-response populations

RUNTIME VERIFIED:

| Parent | Control-response count | Weak-positive count | Matched control count | Training-only | Physical-only | Signed-vs-absolute normalization difference |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| t0 | 10716 | 71 | 10716 | 0 | 0 | +0.0027001080906590147 (~+0.27%) |
| t1 | 18749 | 283 | 18749 | 0 | 0 | -0.010245292038082221 (~-1.02%) |
| t2 | 22356 | 222 | 22356 | 0 | 0 | +0.03129771404630625 (~+3.13%) |

Persisted matched-population coordinate/response quantile differences are zero.
Approximate t1 raw F.3 response median/p90/p99: weak-positive
`0.872 / 1.629 / 3.295`; control-response `0.729 / 1.373 / 3.022`.
These are relative response/topology measurements, not absolute
`P(NPE=0 | true pion, x)`.

Exact accepted limitation flags:

```text
operational_kaon_pid_npe_zero = true
direct_hgc_free_true_pion_tag_present_in_consumed_artifacts = false
direct_pion_to_kaon_misid_calibration_performed = false
absolute_misid_probability_constructed = false
hgc_free_pion_tag_available_in_consumed_artifacts = false
proxy_validity_established = false
proxy_validity_evidence_label = NOT VERIFIED
```

These limit absolute mis-identification interpretation; they do not invalidate
the accepted relative acceptance-correlated refinement.

## Measurement B — t1 signed-normalization sensitivity

RUNTIME VERIFIED:

| Quantity | Accepted value |
| --- | ---: |
| B_signed | 0.3374726818838243 |
| A_abs | 0.3845152767785784 |
| cancellation | 0.8776574098983241 |
| U_signed | 0.29922739636212625 |
| V_abs | 0.34446789697425667 |
| N_signed | 0.8866714623885792 |
| N_abs | 0.8958497042306518 |
| N_signed / N_abs | 0.9897547079619178 |
| Relative normalization difference | -0.010245292038082221 (~-1.02%) |
| Signed-parent closure residual | 0.0 |

`N_abs` is a descriptive positive-support comparator only, never an applied
alternative correction. The same-setting table above shows t1 does not have
the largest normalization contrast.

| t1 source | Count | Absolute-support fraction |
| --- | ---: | ---: |
| prompt | 18749 | 0.9383929586334123 (~93.84%) |
| dummy | 93 | 0.04915800899474474 (~4.92%) |
| random | 963 | 0.012013286056093246 (~1.20%) |
| dummy-random | 4 | 0.0004357463157497171 |

Each source class is internally sign-pure; cancellation is between source
classes. Application support: `in_support_count = 19788`, `ood_count = 21`,
`ood_fraction = 0.001060124185976119` (~0.106%).

| t1 tail | p50 | p90 | p99 | max |
| --- | ---: | ---: | ---: | ---: |
| Raw r | 0.7291584557619877 | 1.3661220249087134 | 2.988399963863005 | 53.331269127265315 |
| Accepted C_signed | 0.8223547127565471 | 1.5407307924725082 | 3.3703576698100095 | 60.14772256636942 |

The finite long tail remains real; it is not clipped, capped, winsorized,
hidden or reinterpreted as a different correction.

## Accepted scientific continuity and decision

SOURCE VERIFIED: F.3 construction uses `SHMS_delta`, `P_hgcer_xAtCer`,
`P_hgcer_yAtCer`. Baseline `w0` owns pion-control -> kaon-background transfer;
F.4 preserves signed canonical-t parents, not individual children or MM regions.
Method A is an acceptance-correlated PID refinement, not merely an HGCer curve.

RUNTIME VERIFIED (historical accepted scopes):
[F.6.2](f6-2-scientific-runtime-closure.md) remains
`CLOSED / RUNTIME VALIDATED`. Its accepted detached scientific gate includes
normalized low-response/baseline/Method-A shapes, F.3 construction variables,
independent acceptance variables `SHMS_xptar`/`SHMS_yptar`, phi and missing mass,
MM x acceptance localization, support/OOD, diagnostic variance/effective
statistics, fixed kaon-window pion-change metrics and paired bootstrap intervals.
It is not an absolute HGC leakage calibration or production promotion.

The accepted [current-baseline comparator](f4-refresh1-current-baseline-authority-comparator-runtime-closure.md)
found `F.2 scientific payload match = true`, `F.3 scientific payload match = true`,
and first changed stage F.4. Current Left/lowe t1 F.4 differs from historical
accepted F.4 only at floating-point scale; t0/t2 carry substantive changes.
Historical F.6.2 acceptance/MM evidence is therefore directly relevant to t1;
that exact continuity claim must not be inherited by current-lineage t0/t2.
The [F.6.3/E.8.4 Left/lowe branch closure](f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md)
retains its execution, live-cache parity, child-change and parent-preservation scope.

INFERENCE / scientific decision: the current t1 diagnostic does not support
strong signed-normalization pathology, broad F.3 OOD application or
training-control -> physical-control population mismatch as explanations of
the large redistribution. Together with accepted F.6.2 acceptance/MM validation
and t1 current-lineage continuity, no concrete evidence from this diagnostic
warrants Method-A redesign or a replacement correction for t1. The redistribution
is compatible with the accepted detached acceptance-correlated PID/MM refinement
picture; this is not a measured detector hardware cause.

## Remaining limits and downstream dependency

NOT VERIFIED: direct absolute `P(NPE=0 | true pion, x)` calibration; absolute
pion-to-kaon HGC mis-ID probability from F.1/F.3/F.4; proof that weak-positive
response uniquely determines latent zero response; specific PMT/mirror/optical/
hardware cause; RF corroboration in this diagnostic; Method-A production
correctness outside accepted detached scopes; current-lineage t0/t2 inheritance
of historical t1 F.4 continuity; canonical-five current-lineage closure;
production promotion; and absolute-SIMC amplitude interpretation.

RF may provide additional corroboration where experimentally available at low
epsilon only; high-epsilon data lack RF. It is not a universal Method-A requirement,
and this returned diagnostic performed no RF corroboration. Missing direct
absolute calibration is an absolute-interpretation limitation, not a foundational
failure of relative acceptance-correlated refinement.

No correction-design implementation is warranted by this diagnostic. Method A
remains detached/non-production; Method B diagnostic/cross-check only and
numerically excluded. Production physics, accepted historical authorities and
`no_empirical_residual` remain unchanged. Canonical-five current-lineage expansion
remains `DEFERRED`; final E.8/F.6.4 remain `BLOCKED`. Absolute-SIMC interpretation
remains separately `BLOCKED` by
`SIMC_normfac_luminosity_and_charge_units_not_source_proven`.

After independent closure diff review, user commit/push and pushed-state
synchronization, the dependency returns to user-deferred canonical-five
reconsideration. No new source/farm work is authorized until the user explicitly
lifts `DEFERRED`; if authorized, first audit exact five-setting current-lineage
authority/validation scope before any implementation or farm contract. Final
E.8 then precedes the explicit F.6.4 production-promotion decision. CURRENT owns
the sole ordinary next action.
