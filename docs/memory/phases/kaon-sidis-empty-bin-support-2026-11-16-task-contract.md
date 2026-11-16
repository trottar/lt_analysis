# D1.a — original-support exclusion in kaon staged SIDIS WLS

## Authority, base and stop

User: November 16 D1.a Empty-Bin Statistical-Weight Correction request.
Exact live/local test base: `c567382027cca3d6aa8c1019c298459d046a33ee`.
This is one separately identifiable, default-off local candidate. Implement only
after the source audit and fixed-weight eight/seven-bin replay. Stop after local
deterministic checks and complete actual-diff review bundle. No commit, push,
remote update, ROOT/main.py/farm run, production-yield regeneration or D2 advance.
Independent actual-diff review precedes any user-only publication handoff.

## Ownership and observed before/after

Before: kaon SIDIS residual WLS can accept a bin with no entries in any original
signed stream and zero experimental Sumw2 because a fitted neutron template
creates a small nonzero residual error. Candidate B Center/high t3 bin 62
([1.066,1.072) GeV) dominates the score. Eight-bin G=-38.7723381448,
H=250.768868521, unconstrained=-0.154613841716 and accepted boundary zero.
Direct exclusion of only 62 with retained weights yields G=10.6094147226,
H=48.4742198501 and amplitude=0.218867157746; this is conditional diagnostic
evidence, not physical background acceptance or a target to force.

After, only when `inpDict['d1a_kaon_sidis_support_enabled'] is True`: original
event support reaches the kaon-side staged least-squares SIDIS step. A bin is
excluded only if all four original source entry counts and all four original
source Sumw2 are exactly zero AND its residual variance is positive from prior
template subtraction. Record native indices and exact reason. Supported signed
cancellation, zero/negative residuals and small errors remain eligible.
Preserve all original histogram bins, contents, errors, plots and yield ownership.
No error floor, positive amplitude floor or alternative estimator is introduced.

## Source audit and interface

The direct proton event loop already owns signed coefficients, prepared MM/t,
nommcuts/hole acceptance and frozen final RF selection. Collect additive original
counts/Sumw2 there, independently of cleaned/proton event weight, with the exact
final no-MM-cut/RF population. Include flow bins and exact MM edges. Original
entries whose cleaned factor is zero still count as observations. Require all
four source streams with valid nonzero coefficients and trees; missing support
is an explicit candidate failure, never an inferred empty population.

Carry metadata setting/global and canonical-t -> rand_sub parent handoff ->
parent builder/detached comparison -> component builder -> kaon entry -> staged
pass. Never reconstruct authoritative kaon parents by rereading trees. Bind
metadata to the source histogram axes; validate lengths, nonnegative integer
counts, finite nonnegative Sumw2 and source completeness. No inference from
residual content/error or cleaned histogram error is permitted.

The generic single-template builder/solver may gain optional objective exclusions
and minimum-bin requirements, with unchanged defaults. Only the enabled A/pi_sidis
least-squares call supplies them. Require at least two retained fit bins for a
one-parameter fit with residual degrees of freedom and positive finite template
denominator. Missing/zero template support or insufficient bins must fail
explicitly and cannot be converted to an accepted zero by downstream refinement.
Use an explicit candidate RuntimeError at the staged-pass boundary on failure.
Do not apply the rule to window-integral SIDIS fits or other stages.

## Allowed paths

- `src/cuts/pion_component_fits.py`: optional objective mask/support plumbing,
  validation, diagnostics and fail-closed candidate behavior only.
- `src/cuts/proton_contamination_weights.py`: default-off additive original
  support metadata in existing event loop; no physics/fill/selection changes.
- `src/cuts/rand_sub.py`: opt-in producer flag and metadata handoff only.
- `src/cuts/pion_t_bin_parents.py`: additive support handoff/parent metadata only.
- `src/cuts/pion_t_template_comparison.py`: forward existing parent support to
  detached clone fit only; no baseline mutation or new alignment scan.
- `testing/test_kaon_sidis_original_support.py` and one compact source-bound
  fixture under `testing/fixtures/d1a_kaon_sidis_original_support.json`.
- This standalone contract at
  `docs/memory/phases/kaon-sidis-empty-bin-support-2026-11-16-task-contract.md`.
- `docs/memory/CURRENT.md`, one concise phase record at
  `docs/memory/phases/kaon-sidis-empty-bin-support-2026-11-16.md`, and generated
  `docs/memory/manifest.json` for warranted active-state/history changes only.

All other paths and pending local work are frozen. Preserve original index,
unrelated tracked/untracked bytes, refs and the seven pre-existing status flags.
Never stage merely to make additions reviewable. If another scientific producer
must change or original support is ambiguous, stop with the exact blocker.

## Scientific freeze

No changes to neutron/Delta/pion-control/joint/protected-Delta algorithms,
templates/alignment/cuts/normalizations, random/dummy/proton subtraction, PID
transfer formula, efficiencies, Local-t ownership, yields, SIMC or L/T. Later
stages can respond numerically to a changed SIDIS coefficient through their
existing path; do not change their estimator or acceptance policy. D1.a does
not resolve protected-Delta failures, Outcome B, full covariance or physical
pion transfer. Method A remains detached/non-production and B diagnostic-only.
The four-criterion 30/20 Lambda (t,phi) mask and reviewed 67/135 cohort are a
different owner from this fit-objective exclusion. Unity/phi-flatness at iteration
zero remain nonrequirements. D2-D5 and six roadmap statuses remain unchanged.

## Validation and deliverables

Tests A-G: original eight-bin boundary; sole unsupported bin exclusion and exact
seven weights/amplitude; signed cancellation; observed zero/negative residual;
missing/zero template; insufficient information; default-off and other stage/
prefix/mode regression; unchanged histogram bytes; producer zero-cleaned-weight
entries retained; missing/malformed support failure. Exercise actual numerical
functions with detached histogram adapters and relevant existing regression
suites in Linux/WSL where available. No local test proves ROOT/farm behavior.

Inspect existing evidence for all 15 candidate-B canonical parents under the
same original-support rule. Retain original weights; report support/exclusion,
scores/amplitudes/failures. Plot native-bin scores, K/neutron/residual, conditional
SIDIS/full pion-model change (other coefficients explicitly frozen for direct
comparison) and 15-parent summary. Give per-native-bin Lambda-window changes;
no fractional-bin fill for nonaligned window edges. Preserve rejected protected
Delta statuses; no claim of complete downstream refit or regenerated yields.

Re-anchor CURRENT narrowly: user reports independently reviewed D1 completed;
numerical origin is source-reviewed, full boundary covariance not validated.
Immediate priority is this D1.a local candidate, not D2. Final status may be
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING only after deterministic checks;
independent source review still pending, no runtime/scientific promotion.
Regenerate/check memory manifest and ordinary health, report nonblocking warnings.

Return one ZIP with request, contract, complete tracked diff AND complete new
files, before/after files, deterministic code/results/logs, 15-parent table,
native arrays/figures, provenance, integrity manifest/receipt and full health.
Acceptance: targeted behavior proven locally, defaults/other stages preserved,
no frozen path changed, all additions reviewable. Hard stop at independent review.

Candidate path safeguards from the source audit: bypass existing canonical-parent
reuse only when explicitly enabled, so an older frozen fit cannot silently bypass
the candidate. Require staged_only kaon mode and fixed alignment (no live residual
shift scan); alignment scoring calls remain unchanged. Unsupported candidate
paths fail explicitly rather than silently broadening alignment/joint behavior.
The verified producer source key dummy_prompt is normalized to dummy in metadata.
