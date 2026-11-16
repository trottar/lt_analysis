# Lambda bin exclusion and reverse reconstruction — task contract

## Authority, starting point and ownership

Starting public/local test: `77c44eb60192f38c6979e71a2ae0cd8a920ada76`.
User handoff: November 12 colleague-inspired four-criterion detached scientific test.
Preserve all pre-existing dirty changes and the previous physical-assessment packet.
The October 28/30 Local-t baseline, unchanged iteration-zero KLambda SIMC,
completed three-component fit, 1024 bootstrap replicas and 116 model responses
are scientific inputs, not editable implementation outputs. Outcome B is unresolved
and is not a prerequisite to this observer comparison. No physical acceptance is implied.

## Exact allowed paths

- testing/analyze_lambda_bin_exclusion.py (new ROOT-independent consumer/CLI)
- testing/test_analyze_lambda_bin_exclusion.py (new deterministic tests)
- this contract
- docs/memory/CURRENT.md (active task/NEXT, preserving prior pending checkpoint)
- docs/memory/investigations/full-t-acceptance-edge-check1-2026-11-05.md (new result)
- docs/memory/manifest.json (regenerate)

All other tracked/untracked bytes, raw index and refs are frozen. No other
production, fit, bootstrap, utility, producer or roadmap source may change.
Artifact-only scientific plotting and packaging scripts remain outside the repository
and accompany the review bundle. No new runtime integration is owned.

## Source/data audit and preserved runtime

calculate_yield_data measures final child histograms through the existing helper,
retains errors and skips invalid frozen-parent children; it does not own this mask.
_xsect_support_data and _xsect_support_simc retain normalized MM children.
average_ratios.f sums setting data and SIMC separately and their reported variances
before division. Reuse that aggregation order and existing setting normalization;
do not import ROOT-dependent producers merely to implement a detached linear sum.
Do not copy legacy sentinel zero/error conventions into unavailable observer values.
The ordinary runtime path is untouched; there are no new production callers.

Inputs: the completed physical-assessment ZIP with exact original 135 MM arrays,
source-bound selected-event acquisition ZIP, conditional covariance factors and all
116 existing response vectors. Verify member manifests and original event ZIP hash,
unique axes, counts and source-level histogram closure before interpretation.
No inverse reconstruction of cells from parents is permitted.

## Trial definition frozen before SIMC evaluation

Signal [1.114,1.126), left [1.102,1.114), right [1.126,1.138) GeV.
Each width is 12 MeV; endpoints are actual retained experimental MM edges. The
signal includes the already-observed peak-containing [1.114,1.120) and
[1.120,1.126) bins. This is an explicit Lambda trial, not an omega window or
fitted mass/resolution result. Use exact retained edge values and half-open
intervals; do not fractionally fill a boundary or round events into a window.

Total entries: all retained kaon nommcuts rows in that physical t/phi cell,
including source streams and zero weights. Signal entries: raw kaon allcuts rows
in the signal interval. Retain per-source counts, nonzero-weight counts and
absolute-weight effective counts separately; raw entries are not net signal or
significance. Signed normalized cleaned signal/left/right yields use original
kaon minus applied Local-t pion arrays, preserving random/dummy and low-epsilon
slow-proton/RF weights. No per-cell charge factor exists or is introduced.

Four criteria, all evaluated and all failure reasons retained:
N_total>=T; N_signal>=S; Y_signal>=Y_left; Y_signal>=Y_right.
Primary provisional pair (T,S)=(30,20). Predeclared sensitivity pairs:
(0,0), (15,10), (45,30), (60,40); no scan or selection based on SIMC agreement.
(0,0) retains the two shape checks and is not the all-retained control.
Explicit unavailable criteria produce indeterminate/excluded-for-comparison status,
never invented counts. Existing empty/weak/nonpositive/finite-MC-absence states
are reported separately. Negative yield or poor SIMC agreement alone cannot be
an exclusion criterion. No SIMC value may enter mask construction.

## Required consumer and outputs

One reusable CLI accepts assessment ZIP, source-event ZIP and fresh output directory.
No embedded workstation paths, farm calls, fits, resampling, model/normalization
changes or source-data deletion. Output every cell and complete criterion/reason
records for each trial. Fixed masks apply identically to data and MC.

Reverse aggregate original/retained/rejected children into 15 parents, five settings,
54 epsilon/t/phi cells and six angular-ratio distributions. Preserve common
normalizations; no scaling remaining cells to recover lost yield or acceptance.
Output counts, signed yields, original background, finite-MC states, tail fractions,
peak-bin/shape metrics where supported, source/support fractions, and all coherent
model variations under the same mask. Missing ratios or shapes are null.

All-retained control must reproduce original 135 cells, 15 parents and six combined
observer-window distributions to retained precision (rtol 1e-9, atol 2e-12).
A control failure stops interpretation. These are observer window integrals, not
new producer-owned Y0 or final acceptance.

For each fixed aggregation A, background replicas/covariance use the same selection:
A*C*A.T, checked against aggregate retained replicas. Retain cross-setting/background
correlations where supplied. Native-bin background factors use only complete bins;
partial-edge covariance remains unavailable separately. Original fixed-weight data
variance remains separate from conditional background bootstrap variance. Do not add
them or invent kaon/background/PID covariance. Model scenarios vary one parent as
published; do not invent simultaneous cross-parent scenarios or probabilistic weights.
No refit/bootstrap rerun, post-selection significance or complete physical one-sigma
claim. Selection bias and acceptance loss must be discussed explicitly.

## Deterministic local verification

Positive: all four logic checks, exact equal-width/window binding, source raw counts
versus signed yields, frozen mask independence from MC, parent/setting/combined sums,
retained+rejected=original, signed/negative observation preservation, and covariance
with nonzero off-diagonals. Model masking must use the same A.
Negative: malformed axes/duplicates, unavailable criterion, invalid thresholds,
nonfinite input, incompatible covariance, source/hash mismatch, no retained cells,
zero MC and nonpositive signed shape support. All failed reasons must survive JSON.
Boundary regressions: half-open lower/upper inclusion; no arbitrary fractional fill.
Actual-data checks: all-retained closure; all 135 exact source histogram/count
bindings; all 116 model vectors and 1024 replicas retained; coherent masking and
covariance replay; no generated criterion/mask depends on MC agreement.
Run existing adjacent bootstrap/propagation suites if interface imports change;
otherwise focused consumer tests and actual-data closure are the relevant gate.
Compile and Python 3.9 grammar checks; manifest and ordinary memory health;
complete tracked/untracked actual diff and byte-exact patch reconstruction.

## Scientific report and acceptance criteria

Embedded original-data figures and tables cover all five settings and 15 parents,
original/retained/excluded MM, shape versus yield magnitude, tail residuals, pion
sensitivity, six angular ratios, support/statistical cost and all exclusion reasons.
Prioritize Center/low t3, Center/high t1/t2 and Center/high t3. Report whether apparent
improvement survives comparable support; no optimization toward MC, invented purity,
full acceptance, fitted peak widths or calibrated significance. Histogram-bin peak
and width descriptors must be identified as such; unavailable fitted quantities stay
unavailable. No 0+/-10000 convention. An improvement with severe support loss is
insufficient for adoption; absent an accepted physical-acceptance loss threshold,
quantify cost and leave acceptability for independent scientific review.

## Hard stop and next gate

Missing exact required input, failed all-retained closure or preservation invariant
stops the affected interpretation; no farm acquisition/rebuild workaround.
Stop after local implementation, science comparison, complete actual diff and review
ZIP. NEXT is independent ChatGPT actual-source and scientific review. No commit,
staging, push, farm execution, production adoption or Method A/B change. The six
roadmap statuses and all prior scoped closures remain unchanged.
