# Parent-t reconstruction Part 1 — November 4, 2026

> November 4 scientific-objective correction: These are historical setting-level integrated candidate diagnostics, not exact selected-Local-t producer yields or a test of final combined data/SIMC ratio flatness. Original numbers/source identities remain; CURRENT owns the ratio-input gate.

Status: `SOURCE REVIEWED` for the user-authorized offline source/numerical audit.
Producer-equivalent Local-t parent closure and full uncertainty: `BLOCKED` on
missing inputs. Not `CLOSED / RUNTIME VALIDATED`.

## Authorization and preserved scope

Parallel extension within existing six-step roadmap: Part 1 baseline parent-t
reconstruction/bin quality; Parts 2 acceptance-aware, 3 Method-A-assisted and
4 validation/cross-section consequences remain `DEFERRED`, not executed.
Local-t selected, not production-adopted; Global/Local-(t,phi) historical.
Outcome B independent Center/high/t3 component constraint remains unresolved.
E.8.2 135/135 data and 135/135 SIMC closure stays accepted and was not rerun.
Method A detached/non-production, B diagnostic/numerically excluded;
no_empirical_residual and all physics/source conventions frozen.

## Source and numerical identity

Local/origin/live public test at startup:
`646ef3e1061f438a7367ee879eec9fb91f493e40`. Tracked worktree/index clean;
two unrelated untracked review files preserved. Only memory edited.
Low-epsilon numerical source `38307c486326faf159cf300fb7c60d6c3758a9f7`;
high-epsilon `9172aba853fd5b5b744c19add992fa4f27466180`.

SOURCE VERIFIED retained November 3 model_rows.csv: same-run Global
producer-log anchor plus exact-allcuts Global-minus-Local-t applied-background
difference. Retains printed precision; not direct serialized Local-t Y0.
Continuous Lambda selection 1.10 <= MM < 1.16 GeV. Fine-grid MM_0 unavailable
in these comparison inputs; coarse FindBin observers are not substituted.
Same-cohort SIMC printed values are code diagnostics, not absolute validation.

External Work handoff: `KaonLT_Parent_t_Reconstruction_Part1_20261104.zip`.
Contains original-data figures, 135-cell inventory, 15-parent tables, selection
sensitivities, original inputs/provenance, source snapshots, reproduction
scripts, actual memory diff and health receipt. ZIP identity in external
receipt; no workstation path belongs in this durable record.

| Numerical file | SHA-256 |
| --- | --- |
| retained model_rows.csv | `aa155dad64ea99922ff2a2f11d95ca701843fb544f11b9b64e66d75dacb80a98` |
| retained cell_inventory.csv | `2eaac7b5dad4e0d4d92a91df08f91cc895a8fd8115229c33ec71f0d3d28379ba` |
| retained high_log_anchors.json | `ae329eb17c4a82b75ef6fc45930d8510d8c5dfce42b0548e0b6b3bfca46c242e` |
| retained source_support.csv | `b292daf83663f7405b4098d79039379c43be4446b4a6cda1ed50f537ee6ffaab` |
| new parents_15.csv | `21f95bbba7d7cbab006cd694b1127bf121afd8f33ebcef9491a877c0c7846d7c` |
| new bin_quality_135.csv | `2adda7317140550d7c7c5e6bccbfd458e0b16804269221756e37fa2b653d8509` |
| new selection_sensitivity.csv | `69c31adc9f184e2bb823a9574c8396a351dd843a3e75f892ffeda20b02e47214` |

## Actual source behavior

SOURCE VERIFIED: remove_bad_bins is a no-op; prune_hist zeros fractional-error
outlier MM bins and resets low-entry/nonpositive-integral histograms. Finite
negative data yields survive producer measurement; explicit parent skip_bin
omits invalid children, distinct from applied-zero or single-scale fallback.
SIMC producer zeroes negative/invalid scalars. Python ratio zero-denominator
or nonfinite results become zero. Fortran average_ratios sums theta-setting
data and MC before division, clips tiny/negative ratios and uses error
sentinels. calc_xsect suppresses nonpositive results; lt_2D_fit selects x!=0,
and displayed graph mean/RMS is not parent yield/propagated uncertainty.
Canonical raw support uses pre-particle records, not signed significance.
No production mechanism changed.

## New numerical findings

All 135 Local-t diagnostics finite; 117 full-MM/control-supported cells,
99 with kaon Lambda events, 18 completely empty plus 18 supported without
kaon Lambda events. Six negative signed values; 16 non-production |Y|<2sigma
reference flags; 33 same-cohort zero printed SIMC denominators. Flags overlap.
Physical detector acceptance and current detached-row production skip status
are NOT VERIFIED. Protected-subfit rejection is not a child-exclusion mask.

| Setting | t | Kaon-populated phi slots | Signed nine-cell sum | Diagonal sigma reference |
| --- | --- | --- | --- | --- |
| Left/lowe | 1 | 1,2,8,9 | 0.078514804 | 0.003897619 |
| Left/lowe | 2 | 1,2,7,8,9 | 0.606673605 | 0.011315577 |
| Left/lowe | 3 | 1,2,3,5,7,8,9 | 0.552642071 | 0.010822316 |
| Left/highe | 1 | 1,2,8,9 | 0.126848081 | 0.013366964 |
| Left/highe | 2 | 1,2,3,7,8,9 | 0.881634216 | 0.035225732 |
| Left/highe | 3 | 1,2,3,7,8,9 | 0.747205916 | 0.033550694 |
| Center/lowe | 1 | 1,2,3,4,5,6,7,8,9 | 0.798197665 | 0.017350110 |
| Center/lowe | 2 | 1,2,3,4,5,6,7,8,9 | 0.715691852 | 0.017306113 |
| Center/lowe | 3 | 1,2,3,4,5,6,7,8 | 0.140045800 | 0.007361857 |
| Center/highe | 1 | 1,2,3,4,5,6,7,8,9 | 1.015515471 | 0.035098158 |
| Center/highe | 2 | 1,2,3,4,5,6,7,8,9 | 0.920311002 | 0.034161896 |
| Center/highe | 3 | 2,3,4,5,6,7,8 | 0.271520947 | 0.019652647 |
| Right/highe | 1 | 3,4,5,6,7,8 | 0.690282761 | 0.026986946 |
| Right/highe | 2 | 3,4,5,6,7 | 0.546483008 | 0.026881665 |
| Right/highe | 3 | 3,4,5,6,7 | 0.403270912 | 0.023769209 |

Existing counts/effective-mC convention. All parent sums positive. Deleting
negative children changes affected parents by at most 0.128%; dropping
no-kaon-Lambda cells by at most 0.0115%. Deleting |Y|<2sigma reference cells
lowers Center/high/t3 by 0.0184511 (6.80%). These are sensitivities, not approved
selections. Largest child supplies about 30-50% of absolute child-yield sum
in Left, 34-43% in Right. Populated-slot span is not physical acceptance.

Diagonal parent relative-error references 1.87%-10.54%, versus parent-specific
median child references 5.18%-27.31%: INFERENCE of conditional counting-noise
gain, not calibrated coverage or reduced bias. Known common data-charge
covariance f_Q^2 Y_i Y_j explicit; fit/template/shared-sample/dummy cross-child
covariance and full total sigma unavailable. Full-sigma fields null; no
independence assumption promoted.

All 15 summed child applied-background scalars agree with retained same-cut/
weight source-partition parent sums within 1e-8. This is arithmetic background
partition, not producer Y0 or physical closure. All 37 scoped hash/grid/finite/
partition checks passed; five authentic-data figures visually inspected.
No fitting, synthetic experiment, ROOT/farm execution or E.8.2 gate rerun.

Ratio of sums differs from mean child ratios on identical positive-printed-MC
support: Center/low/t3 2.57012 vs 5.34433; Left/high/t1 2.92904 vs 4.61901.
Code diagnostics only. Original October 9 producer tables separately summed,
not substituted for current Local-t or used as cross-cohort closure.

## Evidence boundary and recommendation

Useful observed-support diagnostic, not accepted full-2pi extraction. No
equivalent independent direct-parent Local-t Y0/MM_0 in inspected retained
comparison/log inputs. Historical NPZ keys describe child arrays. Direct
parent closure NOT VERIFIED; coarse observers cannot invert fine-grid pruning.
Physical acceptance, full covariance and Outcome B remain unresolved.

First obtain/review one retained same-cohort Local-t producer/covariance
record: exact child MM_0/Y0 and statistical/error components, valid/skip mask,
source/fit/normalization identity, and shared covariance. Any direct parent
requires equivalent cuts/window/pruning/producer semantics. No farm rerun
is authorized. Part 2 merits consideration after independent review and input
sufficiency; it additionally needs angular acceptance/response uncertainties,
matching epsilon information and covariance to distinguish measured coverage
from model completion. Parts 2-4 not executed. CURRENT owns the single action;
no production gate accepted, bypassed or reopened.
