# E.8 baseline-five extracted-yield scientific QA checkpoint — 2026-10-11

## Authority and scope

This evidence-only checkpoint records ChatGPT's completed detached, read-only
inspection of the previously accepted E.8.2 ZIP and procedure pages under the
[checkpoint contract](../phases/e8-baseline-five-scientific-qa-checkpoint-task-contract.md).
It is the existing roadmap's science-memory bracket, not a new scientific phase,
farm run, producer, model acceptance or independent measurement by Codex.
Codex did not ingest or rerun external raw artifacts. All inspection quantities
below are attributed to that ChatGPT examination transferred by the contract.

Accepted farm source: `8e68d2acae91fd97622793df563aba0296403a48`.
Attempt: `Q4p4W2p74_E8_2_baseline_five_8e68d2ac_20261009T223015Z`.
Validation ZIP: 111,668,823 bytes, SHA-256
`df5a64563c7322416f0895c52557bea69f49d73ac0cb47836b0635a27e7660c4`.
The [accepted runtime record](e8-2-canonical-five-baseline-yield-runtime-closure-2026-10-09.md)
owns exact companion status/summary/log identities and the original acceptance.

RUNTIME VERIFIED (prior acceptance, attributed): authoritative baseline
`Y0/MM_0`, exactly `Left/lowe`, `Left/highe`, `Center/lowe`, `Center/highe`,
`Right/highe` (no Right/lowe), three canonical t parents × nine phi children
per setting = 135 cells, and 135/135 data plus 135/135 SIMC producer-table/MM
checks. The detached QA reports all 15 parent policies `valid`.
Accepted baseline numerical yield/support readiness remains
`CLOSED / RUNTIME VALIDATED` only in its existing narrow E.8.2 scope.

## Attributed descriptive inspection

Four signed-negative printed data yields and 36 printed SIMC zeros remain.
Four-decimal printed zeros do not necessarily have empty full-precision MM support:

| Printed-zero sample | Printed zero cells | Empty MM support | Small nonzero rounded entries |
| --- | ---: | ---: | ---: |
| SIMC | 36 | 33 | 3 |
| Data | 36 | 34 | 2 |

These representation distinctions create no new numerical acceptance gate.
Coordinates below use one-based producer bins, not zero-based internal indices.

For Center/lowe t bin 1, `0.4000 <= |t| <= 0.5667 GeV^2`, summing nine phi
children in the cut-selected 1.10–1.16 GeV MM window gives:

| MM quantity | Value |
| --- | ---: |
| Net authoritative Y0 | 0.199397587332 |
| Positive MM-bin content | 0.461284788130 |
| Negative MM-bin magnitude | 0.261887200797 |
| Cancellation fraction | 0.567734309772 |
| Sum of nine printed Y0 table values | 0.1995 |

Full-precision and printed sums differ in representation precision, not producer
consistency. No data value, error or MM object is adjusted by this checkpoint.

For the same parent, existing diagnostic Lambda-window stage integrals are
`POST = 0.85726` (after proton, post-prune) and `PI = 0.14890` (after baseline
pion), with difference approximately 0.70836 or 82.63% of POST. Neighboring
t bins have roughly 8.68% and 0.33% diagnostic removals. Phi bins 4 and 7
have roughly 95.3% and 102.7%, respectively. Their authoritative positive
final yields and reported uncertainties are:

| Producer phi bin | Y0 | Reported uncertainty |
| --- | ---: | ---: |
| 4 | 0.003995 | 0.004859 |
| 7 | 0.005313 | 0.006617 |

These stage integrals are diagnostic objects; `PI` must never replace final
extracted `Y0` or be interpreted as another subtraction stage.

The stored Center/lowe t1 parent fit was applicable under the frozen production
policy. Its reported kaon-region `chi2/ndf` is approximately 246.49, versus
18.33 and 5.86 for t2/t3. INFERENCE: application acceptance alone does not prove
physical fit adequacy. Physical over-subtraction, cause and fit error are
NOT VERIFIED. Existing `w0(MM)` ratio/model application remains frozen;
these observations authorize no new fit-adequacy phase or w0 study.

## Error ownership and interpretation limits

SOURCE VERIFIED against current
[yield producer](../../../src/binning/calculate_yield.py):
`integral_with_stat_error` sums normal MM-bin contents and bin-error variance;
`_measure_yield_with_current_uncertainty` combines that quadrature with specified
data/dummy charge-normalization and existing background-fit fractional terms.
The public yield path consumes that producer-owned arithmetic.
SOURCE VERIFIED against [background configuration](../../../src/utility/background_config.py):
`no_empirical_residual` is active and forces both empirical residual scales
to zero. The supplied inspection records empirical residual fits as inactive.
No uncertainty is recalculated or replaced here.

Complete pion-component parameter-covariance propagation into `Y0` is
NOT VERIFIED. Absolute SIMC normalization/units remain independently `BLOCKED`
at `SIMC_normfac_luminosity_and_charge_units_not_source_proven`.
INFERENCE: descriptive shape-only phi fractions establish no absolute data/SIMC
agreement, correlated systematic uncertainty or chi-squared fit acceptance.
The QA adds caveats only, no production physics acceptance or private `YA` result.

The supporting detached inspection comprises a 15-row parent CSV, 135-row cell
CSV, inspection PDF and targeted w0 audit PDF in the prior ChatGPT session.
They are inspection aids, not authoritative producers or new farm acceptance;
no binary/CSV/raw-event artifact is added to repository memory.

## Existing roadmap dependency and preserved boundaries

The approved order remains E.8.2 baseline -> E.8.3 detached Method A -> F.6.3
parallel procedure -> E.8.4 impact audit -> final E.8 -> F.6.4 promotion decision.
Prior Left/lowe E.8.2/E.8.3/F.6.3/E.8.4 and other closed scopes remain intact.
Method A stays private/non-production; Method B stays diagnostic-only and
numerically excluded. Subtraction, yields/errors, SIMC, cuts, binning, models,
acceptance and all frozen scientific ownership remain unchanged.

Private canonical-five F.4 equivalence remains `BLOCKED` at
`f6_3_scientific_equivalence_mismatch:f4:$.parents[0].absolute_baseline_parent_sum`.
The [prior forensic evidence](e8-4-f1-forensic-reproducibility-diagnostic-review-2026-10-08.md)
reports three Left/lowe parent application-event counts identical, all six
signed/absolute sums not exactly equal and reviewed F.1 absent. It grants no
canonical acceptance, presumed roundoff exemption or tolerance relaxation.
Final E.8/F.6.4 remain `BLOCKED`; absolute-SIMC remains separately blocked.

[CURRENT](../CURRENT.md) owns the sole substantive successor, conditional on
actual-diff acceptance, user publication and pushed-state synchronization:
audit the first blocking exact private F.4 invariant and available admissible
evidence against the existing roadmap contract; determine the earliest valid
resolution gate without changing acceptance criteria or retrying the owner.
Further F.1 archival/numerical archaeology remains `DEFERRED` unless that
prerequisite audit proves it materially unavoidable. This checkpoint executes
no successor investigation, source change, Method-A implementation or farm work.
