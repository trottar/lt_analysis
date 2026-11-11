# Detached joint pion MM-delta fit -- November 11, 2026

## Authority, objective and starting state

Starting authoritative test HEAD: f346579713168bb9c36b48128df51117c83059ae.
The November 11 user handoff accepts the regenerated 15-role/45-parent
projection and fixed-model diagnostic. It authorizes revised relative neutron,
Delta and SIDIS amplitudes as the primary scientific objective. Template or
acceptance mismatch is a model-uncertainty concern, not a prerequisite that
postpones this bounded detached fit. Local-t selection is settled.

The existing fixed-comparison memory edits are pending and must be retained.
All other pending work, index bytes and refs are frozen. The six roadmap
milestones/statuses, Outcome B and Method A/B ownership remain unchanged.

## Allowlist and frozen interfaces

- testing/fit_pion_mm_delta.py (new, ROOT-independent detached diagnostic CLI)
- testing/test_fit_pion_mm_delta.py (new deterministic tests)
- docs/memory/phases/joint-pion-mm-delta-fit-2026-11-11-task-contract.md (this contract)
- docs/memory/CURRENT.md
- docs/memory/MEMORY.md (only durable priority/uncertainty clarification)
- docs/memory/investigations/full-t-acceptance-edge-check1-2026-11-05.md
- docs/memory/manifest.json

All production src/, existing testing code, extraction/reference builders,
roadmap, no_empirical_residual, cuts/events/axes/alignment/normalization,
production coefficients/weights, Method A/B and final-yield owners are frozen.

## Audited path and before/after behavior

The November 11 fixed-comparison ZIP retains source-bound numerical_inputs.npz,
the original fit diagnostics/alignment extracts, signed source moments, hashes,
all 15 parents and complete experimental support classifications. The supplied
regenerated projection supplies fresh raw per-delta channel histograms,
alignment operators and same-cohort closure references. No ROOT execution is
needed for this offline analysis.

src/cuts/pion_component_fits.py::_solve_joint_template_amplitudes and the
current full-range regularized joint implementation use SciPy bounded
least squares with nonnegative amplitudes and staged targets/recorded widths.
testing/pion_local_t_transfer_variance_comparison.py documents exact original
staged targets, prior widths, fit-bin lists and conditional covariance limits.
Their ROOT-oriented owners remain unmodified. Reuse their bounded least-squares
mathematics and recorded priors through the same SciPy numerical primitive;
do not import fake ROOT, AST-extract production functions, redesign staging,
or claim a new production integration.

Previously there were fixed-B residuals only. The new diagnostic fits one
nonnegative amplitude per channel per canonical-t parent, shared across ten
delta bins; settings/parents remain independent. Default covers all 15, with
Center/high t3 evaluated first. Explicit setting/parent selection is allowed
for deterministic development; it must be recorded, never called full scope.

## Numerical definition and scientific ownership

Preserve the original parent normalization and native MM bins. Fit only
occupied joint cells with positive finite signed counting Sumw2 in the
recorded original MM fit-bin domain. Keep empty data and absent finite-MC
support as distinct exported states, without invented variance or sparse
threshold. Zero/negative signed observations remain observations.

Use the recorded staged amplitude targets and regularization widths, with
nonnegative bounds and no cross-parent equality. Verify replay of the original
marginal joint objective/coefficients before fitting new dimensions.

The primary candidate is regularized generalized least squares with counting
plus conditional finite-MC covariance frozen at the original B amplitudes.
This defines an explicit fixed quadratic objective and permits fixed/candidate
comparison under the same weights. MC covariance must differentiate common
stochastic parent normalization and fixed interpolation, including off-diagonal
terms. Counting-only and covariance-re-evaluated-at-candidate fits are separate
statistical sensitivities, not silently changed objectives. No same-data prior
is relabeled independent evidence.

Report conditional coefficient covariance using the linear estimator's
sampling derivative, with separate counting and MC contributions and the
retained deterministic regularization. Also report penalized curvature,
data-only rank, correlations, constrained boundary flags and parameter profile
contours. Near-boundary symmetric errors are not calibrated confidence limits.
Use delta-objective=1 profile contours only as conditional penalized contours.

Quantify model sensitivity by leave-one-occupied-delta-bin-out refits and by
retained, non-rejected alignment candidates/baseline shifts, normalized once
at the same parent. The central alignment remains fixed. Report actual envelopes
and surviving residuals; these variations have no invented probability law,
are not added to statistical errors in quadrature, and do not exhaust physical
acceptance/shape uncertainty. Retain a no-regularization and Lambda-excluded
fit as explicit robustness studies, with no adoption by default.

Compare fixed versus candidate in broad native MM and complete native Lambda
bins; preserve exact allcuts observations and boundary-bin integral bounds.
Report model mass in empty/unsupported cells, residual quadratics and original
signed source support. Distinguish in-sample improvement, conditional precision,
model sensitivity, true component identifiability and physical validation.
Evaluate the pion-control denominator and beneath-Lambda mixture only; keep
kaon-side amplitudes and applied production pion weights frozen.

## Verification and forbidden shortcuts

Positive checks: real 15-role/45-parent closure, original marginal fit replay,
actual fixed arrays unchanged, signed-source closure, correlated-normalization
Jacobian/PSD, known deterministic coefficient recovery, covariance decomposition,
nonnegative boundary/profile, all 15 numerical fits and regional arithmetic.
Negative checks: malformed/nonfinite arrays, nonpositive variances on fitted
rows, missing provenance/prior fields, member/hash/reference mismatch,
template closure mismatch, no supported fit cells, invalid selections,
nonconverged solver, invalid output and duplicate publication.
Regression: default remains detached; no ROOT/production imports or changes;
scale equivariance and recorded priors, signed negative/empty-cell retention,
original historical hashes not imposed on regenerated events.

Never tune normalization, vary data cuts, relabel source events as species,
floor empty-cell variance, invent acceptance/alignment covariance, add a
production fallback, change source physics or promote Outcome B/Method A/B.
No web research or extra farm evidence is necessary for this bounded objective.

## Local validation, deliverables and hard stop

Create this single contract directly in the repository before implementation.
Short Codex launch instruction: read root AGENTS/startup memory and this
contract; implement only the detached diagnostic, tests and allowed checkpoint;
run deterministic checks and all-15 real-array evaluation; preserve source and
pending work; return full actual diffs and stop for independent review.

Run deterministic local numerical tests, original fitting/extractor regression
checks, compilation and Python 3.9 grammar. Regenerate/check manifest and run
ordinary memory health; report sizes and classify existing warnings nonblocking,
batched at the next checkpoint. No farm/runtime claim follows from local tests.

Provide complete cumulative and incremental actual diffs, full new files and
before/after snapshots, real input identities, result arrays/tables/correlations,
uncertainty/sensitivity breakdown and original project-data plots embedded in
a report. Inspect complete diffs and verify index/refs/unrelated bytes. Include
exact conditional user-only Git commands without executing them. Source review
is independent after this implementation; stop at that gate. No commit, push,
farm operation, automatic production adoption or additional phase follows.
