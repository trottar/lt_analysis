# Source-aware detached pion bootstrap -- November 11, 2026

## Authority and objective

Starting local/public test HEAD: f346579713168bb9c36b48128df51117c83059ae.
The accepted joint-fit source is local/unpublished and matches the reviewed
kaonlt_joint_mm_delta_20261111_review_bundle.zip byte-for-byte. Reuse it unchanged.
The user authorizes source-aware data/SIMC bootstrap and correlated model
uncertainty treatment, all five settings/15 parents, Center/high t3 first.
Do not wait for externally owned PID systematics or promote production weights.

## Allowlist and frozen files

- testing/bootstrap_pion_mm_delta.py (new detached numerical CLI)
- testing/test_bootstrap_pion_mm_delta.py (new deterministic tests)
- docs/memory/phases/pion-bootstrap-transfer-2026-11-11-task-contract.md
- docs/memory/CURRENT.md
- docs/memory/MEMORY.md (short durable covariance/PID-interface clarification only)
- docs/memory/investigations/full-t-acceptance-edge-check1-2026-11-05.md
- docs/memory/manifest.json

All existing source/tests, accepted fit/extractor/contracts, prior pending edits,
roadmap milestones, production physics, normalization/cuts/alignments/PID,
Method A/B, Outcome B, refs and raw index bytes are frozen. Cumulative review
diffs must include the accepted but unpublished joint-fit additions; incremental
diffs isolate this bootstrap extension without changing accepted source.

## Audited path and limitations

Reuse testing/fit_pion_mm_delta.py's loader, normalized MC covariance and bounded
solver. Reuse actual experimental rows from the source-bound November 9 Check 1
ZIP and detector-coordinate/source-entry bindings from the validated ten-bin
packet. The regenerated projection retains source entry, t/phi/MM/delta bins and
event weight for all 15 channel files. Retain hashes and all common row identities.
Source-aware Poisson event multiplicities are an established local diagnostic
convention; do not independently fluctuate overlapping histograms.

The original transfer owner is build_simc_shape_pion_control_weights in
src/cuts/pion_component_subtraction.py, called by the Local-t observer adapter.
It uses the protected applied kaon amplitudes (excluding signal templates), the
pion-control fitted denominator, recorded floor/lower/upper clipping and zero
flow weights. snapshot(None) explains the stored clip_max unavailable sentinel
when the producer explicitly emits that key as None; verify this exact source
meaning rather than guessing a clipping policy.

Post-cleaning kaon rows exist and reconstruct the frozen kaon MM input, but raw
learning populations, original physical event IDs common across PID trees and
the complete staged/protected/learned-PID refit path are not reconstructed by
this detached estimator. Full numerator-fit and upstream-learned covariance
remain unavailable. Do not call all kaon data missing. Complete control covariance
and conditional transfer with protected A fixed; numerator and denominator use
the SAME bootstrapped channel templates, retaining their known SIMC correlation.
The SAME pion event draws feed fitted B and transferred background, retaining
their known covariance. Do not independently draw a numerator covariance.

## Bootstrap definition

One Poisson(1) multiplicity per unique retained source/tree/entry identity in
each actual experimental cohort; reuse it across selections, parent labels,
MM/delta/phi bins, fit inputs and application. A row is never keyed by entry alone.
SIMC files/channels are independent from data and each other, keyed by recorded
file identity and entry. Shared selected MC rows use one draw for parents and phi
children. Unknown cross-PID/shared-learned/cross-setting systematic relationships
remain explicitly unavailable. Do not invent shared raw-event equivalence.

Reconstruct raw Sumw=sum k*w and Sumw2=sum k*w^2 (not k^2*w^2), native arrays,
fixed alignment and stochastic common parent normalization for each replica.
Fit three nonnegative amplitudes with the recorded staged targets, prior widths
and original MM domain. Recompute counting and finite-SIMC precision for each
replica at the original B convention; positive finite replica Sumw2 controls
fit membership, without an added sparse-bin threshold. Stored priors, alignment
and upstream PID coefficients remain fixed/conditional, not silently included.

Run paired data-only, MC-only and joint ensembles with reproducible seeds and
replica IDs. Default 1024 replicas, with exact first-half/full covariance and
batch stability diagnostics, failure counts and boundary frequencies. Inputs
from disjoint statistical sources are drawn independently; report actual joint
covariance and nonlinear interaction/cross terms, not an unjustified quadrature
sum. Failed/unsupported replicas keep their identifiers/reasons. Covariance
requires at least two common supported replicas; never silently zero-fill missing
results, drop one tail or substitute analytic covariance after failure.

## Transfer and correlated model response

For each supported replica use the established MM-only numerator/denominator
ratio and actual signed pion events to compute MM/phi background, exact Lambda
selection and integrated yield. The numerator A is fixed but its channel
templates are shared with the denominator in the MC replica. Verify baseline
transfer and signed application against original cached rows before bootstrapping.
Record production floor/clipping/flow behavior and unsupported denominators.
Label this conditional transfer covariance, not full physical subtraction error.

Keep deterministic alignment/leave-delta influence cases and residual-derived
shape directions distinct from bootstrap statistics. Export full correlated
response vectors and their outer-product kernels, plus amplitude/background
envelopes and discrepancy coverage. Kernels have no assigned probability law;
model covariance itself remains null until justified nuisance covariance is
provided. Do not equate scan ranges or data-trained discrepancy directions to
one sigma, alter the central templates to force agreement, or assign PID zero.
Keep PID covariance null/external, with axes/units/provenance and cross-covariance
requirements explicit for later combination. No default independence assumption.

## Validation, deliverables and hard stop

Before implementation place this single contract in the repository. Inline
Codex launch: read root AGENTS, startup memory and this contract; implement only
the detached bootstrap/tests/checkpoint; validate real input closure and all-15
replica results; deliver actual cumulative/incremental diffs and stop for review.

Test known shared-event covariance, repeated-identity coherence, signed weights,
Sumw2 multiplicity, parent normalization, independently sampled MC channels,
exact transfer, source/array/hash failures, boundaries, reproducible replicas,
unsupported denominators and covariance/interaction algebra. Check frozen source
identity, original regression suites, compilation/Python 3.9 grammar, manifest
and ordinary memory health. Existing size warnings are nonblocking and batch
at next milestone. Never fabricate ROOT/farm or physical coverage validation.

Produce actual all-15 amplitude/bootstrap covariance, correlations/boundaries,
conditional Lambda/MM/phi background covariance (low-rank factors permitted
with exact axes and supported-replica masks), model response kernels, unsupported
terms and a null external-PID interface. Include readable original plots embedded
in a scientific report, convergence/closure evidence, complete actual diffs,
full snapshots and exact conditional user-only Git commands in one review ZIP.
Preserve all unrelated bytes/index/refs. Stop for independent source/scientific
review; no commit, push, farm execution, production adoption or extra phase.
