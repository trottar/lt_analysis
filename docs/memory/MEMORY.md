# Durable KaonLT project knowledge

## Authority and evidence

`docs/memory/` is the repository-owned durable continuity layer. Native
assistant/Codex memory is supplemental only and cannot override current source,
repository memory, or supplied farm evidence. `CURRENT.md` is the sole
ordinary active-state authority. The handoff records exceptional transfer state
only and cannot override CURRENT.

Work can and should directly update authorized repository memory, including
scientific priorities, investigation/accountability records and workflow
preferences. Memory-only updates do not require a separate Codex task,
contract, launch prompt or farm run. Preserve pending changes, inspect the
actual diff, regenerate the manifest and run memory-health checks. Reserve
the Codex implementation workflow for source-changing work. User-only
publication and farm authorization remain separate. This is the user's
November 4 clarification and supersedes the former Codex-only memory rule.

At substantial-work start, establish the actual branch, HEAD, and worktree.
For implementation, current source and relevant diff outrank stale summaries;
for runtime acceptance, fresh farm/runtime evidence outranks source review.
Source review never establishes farm validation. Historical and chat material
is lower authority than canonical current, evidence, phase, and decision
records.

## Source and provenance semantics

Interpreter code aliases need both lexical provenance under an interpreter-owned
root and canonical containment; canonical containment alone admits arbitrary
external aliases, while unconditional lexical/canonical equality rejects genuine
lib64-to-lib environments. Loader-only source/extension reads remain distinct
from package data, bytecode and scientific reads. Review adjacent read, finder,
namespace and identity guards together. The
[October 22 incident](investigations/kaonlt-matched-cohort-farm-readiness-failure-2026-10-22.md)
records inadequate assistant readiness judgment and lost scientific throughput.
Replaying a Python 3.9 path shape on another version does not validate Python 3.9
or the installed farm dependency closure; local checks never establish runtime
acceptance. CURRENT owns the unresolved acquisition gate.

Reviewed scientific source, repository observation, farm-evaluated source,
renderer source, bundle/profile source, and artifact provenance are distinct
roles. A source review may be anchored to an analysis source older than the
repository observation. A stored observed HEAD is a timestamped observation,
not a permanent statement of live identity.

Record a farm-evaluated identity only when direct applicable farm evidence
supports it. Documentation or memory-only work does not automatically
invalidate an older reviewed analysis source; reconcile the intervening diff
and re-review only when the relevant analysis source or test scope changed.
Never collapse provenance roles merely for convenience.

Exact serialized lineage identity and scientific-equivalence identity are
distinct. A permitted bridge reconstructs the current lineage through the
same public scientific builders, compares an explicit closed scientific
projection exactly, preserves both provenance roles, and consumes only
current-lineage transient factors. It never silently interchanges those identities.

## Scientific ownership and production ordering

Preserve separate scientific owners for random subtraction, slow-proton PID
contamination, pion-production background, HGCer diagnostics, SIMC, yields,
and cross sections. The production ordering remains random/dummy -> frozen
binning -> slow proton -> pion subtraction. A setting-wide K-Lambda gate
controls whether proposed proton weights become applied; never partially
commit per-t proton results.

Legacy empirical residual Fit 1/Fit 2 are dormant historical machinery. The
accepted active profile is `no_empirical_residual`, which forces both empirical
residual-background scales to zero; those fits are not part of the current
production or presentation chain and must not be reintroduced by E.8, Method A,
diagnostics, or presentation work.

Pion-alignment comparisons use fixed evaluation envelopes. HGCer diagnostics
consume the frozen baseline rather than redefine it. Cuts, templates, priors,
component definitions, normalizations, binning, efficiencies, acceptance, L/T
separation, and uncertainty propagation remain frozen unless a narrow approved
contract explicitly owns a change.

### Pion statistical interpretation

SIMC pi_n/pi_delta/pi_sidis shapes estimate pion contamination. Subtract once
at iteration zero; later K+Lambda SIMC iterations reuse it. Local-t is the
user-selected scientific baseline, not production-adopted. Global/Local-(t,phi)
remain historical comparisons; initial K+Lambda SIMC is provisional.
The fitted pion-control denominator is essential: pion-only numerator/denominator
weights act on signed experimental pion events. Small subtraction, positive
yield, agreement or a near-unity data/SIMC ratio cannot validate that transfer.
High epsilon has no RF slow-proton subtraction. Scalar empty-error fallback
fails equivalent-unit invariance. Source-aware refits establish conditional
behavior, not universal coverage or physical bias/MSE. Staged priors and the
residual Lambda gauge reuse data; shapes do not fix kaon-PID composition.
Protected Delta acceptance cannot validate earlier pion removal. Judge applied
background beneath Lambda, not full-MM chi2. Missing Sigma0 input is a separate
improvement. Diagnostics are not corrections; see the November 3 investigation.

### Clean-yield, angular-ratio and SIMC-iteration objective

The established scientific objective is a physically defensible clean
background-subtracted K+Lambda experimental yield and its corresponding
K+Lambda SIMC comparison, with data/SIMC ratios flat versus phi in each t bin
at both epsilons. At fixed epsilon combine complementary spectrometer settings
by summing data and SIMC yields separately before division. For branch b and
SIMC iteration n:

\[
R_{\epsilon,t,\phi}^{(n),b}
=\frac{\sum_{s\in S_\epsilon}Y_{\mathrm{data},s,t,\phi}^{(0),b}}
{\sum_{s\in S_\epsilon}Y_{\mathrm{SIMC},s,t,\phi}^{(n)}}.
\]

At Q4p4W2p74, S_low is Left+Center and S_high is Left+Center+Right: six
epsilon-by-t distributions, 54 nominal cells. Do not average setting ratios.
Apply experimental corrections and selected pion subtraction once at iteration
zero. Immediate angular QA is agreement with a constant versus phi using
valid covariance where available; unity is not required at iteration zero.
Keep those experimental yields/subtraction fixed in later SIMC cross-section
model iterations; established iteration machinery should bring ratios toward
unity while reproducing data angular dependence. Never retune pion background,
area-match or add an artificial normalization to force unity.

Method-A evaluation compares selected Local-t versus selected Local-t + A on
the same cohort, fitted Local-t baseline, cuts, bins, normalization and common
iteration-zero SIMC denominator. Its primary observable is the change in ratio
versus phi, not yield increase, parent-sum appearance or smoother error bars.
Distinguish angular flatness, normalization, physical pion-background support
beneath Lambda and statistical/correlated/background uncertainty. Neither a
flat ratio nor unity proves purity; incorrect subtraction can compensate for a
SIMC-model error. Invent no thresholds, covariance or fitting algorithm.
Accepted setting-local E.8 debugging is not final combined-yield acceptance.
Missing matched producer evidence stops that comparison, not a reason to
commission another historical or conditional substitute.

### Observable-chain and combined-setting review

Before commissioning models or investigations, trace producer -> serialized
Y0/MM_0/error/validity records -> physics_lists yield files/list.settings ->
average_ratios setting combination -> the actual observable. At Q4p4W2p74,
combine Left+Center low epsilon and Left+Center+Right high epsilon into 54
(t,phi,epsilon) cells before judging final angular gaps. Sum data and SIMC
separately; do not average setting ratios. Prove branch/cohort ownership:
t_bin configuration alone does not identify the selected Local-t template
branch. Missing offline records are not missing functionality. Method A is a
separate subsequent comparison; an off-target conditional fit clears no gate.

## HGCer and Method-A/Method-B boundaries

Method A and Method B are independent. Method B uses frozen upstream records
and same-canonical-t relative closure; it uses no Method-A numbers, no cross-t
pooling or interpolation, and no absolute neutron normalization. Adaptive
Method B remains DO NOT PROMOTE. Method B is diagnostic/cross-check/historical
comparison only and never changes production pion weights.

Presentation cannot recompute a diagnostic or construct a correction. Method A
remains a candidate numerical HGCer path only after an explicit, validated F.6
production-promotion decision. Its positive-response training and physical
application populations remain distinct: training uses the established
prompt/noRF/nommcuts NPE>0 population, while physical application uses the
authoritative NPE>2 cache. Keep their provenance and fingerprints separate.
Future response application is event-level and parent-t normalized; never
renormalize `(t,phi)` children independently. Detached F-stage diagnostics do
not become production merely because they are accepted.

## Diagnostics and runtime validation

Scientific evaluation starts with full-MM component reconstruction under verified
common selection, coordinates, binning and charge/efficiency/acceptance
normalization, then signal-region peak/width/tail/residual closure, then the
unbiased background-subtracted K+Lambda data/SIMC yield ratio for matching
`(t,phi,epsilon)` bins and established L/T extraction. Show spectra and
signal-region zooms; background-fit chi2 or a large control-domain mismatch alone
does not establish physical contamination or signal bias. Document any floated
background amplitude or physically justified MM offset as diagnostic shape
nuisances, distinct from authoritative absolute normalization. Never force unity
by area-matching final signal SIMC. Unproven absolute SIMC units still block
absolute interpretation; diagnostic shifts/scales authorize no production change.

Historical October initiating observation, not November model preference:
the user's visual comparison places deterioration during canonical-t model construction/
refinement, before final pion subtraction; the setting-wide model appears more
sensible. This observation and the report of earlier working behavior establish
neither accepted model quality, an exact first-bad operation, identical canonical-t
inputs nor a working SHA. The setting-wide full procedure remains in the October
`pion-fit-debug` PDF; fewer main-PDF pages reflect presentation routing, not
removal of global fitting. Audit all three SIMC components through model,
transfer and signed application; one amplitude (`B_n=0`), algebraic closure,
large SIMC totals or statistical support cannot establish scientific model
adequacy or the initiating regression/root cause. The
[investigation](investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md)
retains attributed evidence and uncertainties; [CURRENT](CURRENT.md) owns active
scope and priority. No code change or physics acceptance follows from this lesson.

User-requested future scientific intent, not an accepted correction: preserve
the complete setting-wide staged `pi_n`/`pi_delta`/`pi_sidis` procedure at t and
later `(t,phi)`. Anchor t fits to global without suppressing genuine potentially
large t dependence; consider full child refinement anchored to its own t parent
only with meaningful local support and covariance. Parent and child share events,
so the parent fit cannot be treated as an independent measurement. Current
parent alignment, local-stage regularization and inherited parent event weights
are distinct from such amplitude/covariance conditioning. No fixed percentage
clamp, numerical prior width, covariance model or automatic fallback algorithm
is approved; hierarchy efficacy and production acceptance remain NOT VERIFIED.

When applicable, trace persisted diagnostics through producer ->
serializer/checkpoint -> checkpoint-first payload -> consumer -> renderer.
Local deterministic tests do not establish farm integration. Farm, ROOT/PyROOT,
and full `main.py` behavior require direct farm evidence.

E.8 presentation is the complete visual audit of the kaon missing-mass and
signal-region yield chain: every substantive stage shows its authoritative
input spectrum, applied component/treatment, and authoritative output spectrum,
culminating in canonical `(t,phi)` missing-mass spectra and extracted yields.
It only consumes authoritative runtime objects, persisted snapshots, accepted
detached Method-A artifacts, or later branch outputs. It never recomputes a
fit, factor, correction, normalization, or yield for plotting. Method-A
reweighting remains parent-preserving (`w0 -> w0*C`) with no independent child
normalization; Method B remains numerically excluded.

The [four-category presentation intent](decisions/e8-full-analysis-procedure-roadmap.md#user-clarified-final-presentation-requirement)
requests baseline, Method A, reweighted baseline and K+Lambda SIMC. Their
scientific roles differ: F.3 independently derives relative detector response
`r` from positive-response training on accepted `hgcer3`, not an absolute
pion-mis-ID probability or extracted kaon yield. F.4 uses signed baseline
`b_j=s_j*w0_j`, `B=sum b_j`, `U=sum b_j*r_j`, and valid-parent
`C_j=r_j/(U/B)` to preserve `sum b_j*C_j=B`. Parent-normalized `C` depends
explicitly on `w0` and signed source coefficients; no child normalization.
F.6.3 fills `B_pi_A` with `s_j*w0_j*C_j`, forms `MM_A=pion_input-B_pi_A`,
and extracts `YA` with the existing uncertainty helper. Existing `YA/MM_A`
already is the Method-A-reweighted baseline final branch; do not invent a
missing independent reweighted-baseline producer or treat it as Method-A-only.
Retain six model/refinement presentations: Local-t versus Local-t + Method A
is primary; Global and Local-(t,phi), each before/after A, are the other four
historical combinations. Independent Method-A response diagnostics are separate;
no standalone A-only final yield or seventh yield series. `Y_rw/MM_rw`
remain request notation, not new keys. `B_pi_A` is background, not final kaon
signal. Current E.8.4 MM overlays remain three-way with `Y0/YA` summaries;
Method-A response/weight evidence may be shown as diagnostics, without false
yield units or a promise of four independent final series. CURRENT retains
the sole NEXT; this semantic correction creates no phase or runtime acceptance.

Correction-factor percentile tails alone do not measure tail dominance of
`s*w0*C`. Compare signed and absolute baseline/adjusted contributions; assess
effective statistics with an explicit nonnegative support measure and moment
definitions, retaining signed cancellation context. Marginal factor percentiles
cannot reconstruct joint event-level `(phi,MM)` support or an uncertainty.

The current-baseline F.6.3 Method-A full-analysis path has accepted direct
runtime evidence for Q4p4W2p74 / Left / lowe only: real canonical child yields
change while signed parent-t normalization is preserved. Parent equality is a
constraint, not the magnitude of the redistribution effect. See the
[Left/lowe evidence](evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md).
No canonical-five closure or Method-A promotion follows.

SOURCE VERIFIED at the yield-first checkpoint: E.8.2 final canonical MM
consumes after-baseline-pion `after_pion_final` and stored baseline `Y0`;
Lambda-window stage integrals (including `PI`) are diagnostics, not extracted
yields or an additional subtraction stage. Historical E.8.3 F.6.1 is separate
from private current F.6.3/E.8.4 `MM_0/MM_A` and stored `Y0/YA`. INFERENCE:
cross-branch or diagnostic-integral comparisons may create apparent discrepancies;
that possibility resolves none without exact producer-owned closure.

E.8 shareable pages consume the existing normalized per-child SIMC MM support
from `hist["_xsect_support_simc"]["mm"]` after the existing SIMC-yield producer.
Clone only for display; never reload, refill, renormalize, or substitute a
setting-wide SIMC shape. Yield summaries consume stored branch values, and
parent closure consumes hash/fingerprint-pinned candidate F.4 records.

An explanatory aggregate page and an applied downstream branch cannot be
interpreted as the same correction lineage without exact source, artifact and
fingerprint identity. Historical accepted E.8.3 F.6.1 aggregates and the
current-baseline F.6.3/E.8.4 candidate are distinct; individual page correctness
does not establish cross-page linkage. Renderer success does not establish
scientific linkage consistency.

E.8.3 has accepted direct farm evidence only for Q4p4W2p74 / Left / lowe at
`63eb2bf6`: provenance/preservation passed,
authority-page full-value legibility is repaired, and pages 66–72 are unchanged
(pixel-identical). The [accepted closure](evidence/e8-3-left-lowe-runtime-closure-2026-10-07.md)
attributes acceptance to supplied artifacts/prior ChatGPT review. This narrow
runtime/presentation closure does not link historical F.6.1 to current F.6.3,
promote Method A, validate absolute SIMC or close canonical-five/final E.8/F.6.4.

Current-branch yield-impact presentations require explicit histogram-scalar
closure against exact producer-owned integration semantics. When fractional
signed-yield changes appear disproportionate to spectrum overlays, report
signed integral, positive-bin support, negative-bin support and absolute
support separately; cancellation is a hypothesis until those numbers exist.
For the present Left/lowe work item, fail-closed identity checks, prior arithmetic
checks and accepted narrow evidence support treating bookkeeping as internally
consistent. Stale scalar/wrong histogram is no longer the primary unresolved
question. Arithmetic consistency does not validate the scientific origin or
production correctness of Method-A redistribution. The accepted current-lineage
[Left/lowe diagnostic](evidence/e8-4-left-lowe-method-a-current-lineage-diagnostic-runtime-closure-2026-10-04.md)
now measures t1 signed/absolute support, source decomposition, finite response/
correction tails, child/MM redistribution, coordinate dependence and support/OOD.
Its ~-1.02% signed-versus-positive-support normalization difference, ~0.106% OOD
fraction, exact matched control populations and zero recorded parent residual
do not support strong normalization pathology, broad OOD or population mismatch
as explanations of t1. INFERENCE: with accepted F.6.2 acceptance/MM validation
and comparator t1 continuity, no correction redesign is warranted by this
diagnostic. Direct absolute calibration and hardware cause remain NOT VERIFIED;
the general histogram/scalar rule remains.
SIMC/data amplitude interpretation requires proven same-cell, same-object,
same-normalization and same-unit provenance, including the relevant
integration-window integral. Cloning existing SIMC without renderer scaling
alone does not prove those identities. See the [post-farm investigation](investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md).

Fix.5.4 source audits use normal-bin, unweighted-by-width integration with
flows excluded; exact content/error ownership fingerprints include flow bins.
Current pion-template aggregates use current F.6.3 children only. F.4 parent
application sums cover broader MM support and cannot be equated to cut-window
template integrals. SIMC `iter_weight * normfac / Ncontribute` has no proven
luminosity/effective-charge unit authority in the traced source. Preserve that
absolute-comparison blocker without inventing a scale; data/SIMC agreement or
disagreement cannot be inferred from the current overlay amplitude.
Unresolved absolute SIMC units block the two absolute-SIMC page families and
their interpretation only. They do not make a valid current-F.6.3 E.8.4
data-only/identity/yield/parent-closure payload unavailable.

Fix.5.4 numerical source `761fbb6` (full identity in CURRENT and the phase record) passed
independent ChatGPT actual-diff and pushed-state review per the supplied
Fix.5.5 contract; this is source acceptance, not new farm closure.
Fix.5.5 presentation consumes stored child signed/absolute support and stored
current-F.6.3 aggregates without histogram reintegration or child reaggregation.
Baseline-last solid/open-marker blue and wider dashed magenta encode identity
beyond color. Historical E.8.3 must visibly identify its accepted F.6.1 lineage
as separate from current F.6.3/E.8.4. Fake ROOT checks establish draw commands and
immutability. The later Fix.5.8 farm gate validates the narrow Left/lowe
owner/presentation result and actual PDF legibility; canonical-five closure
and absolute-SIMC interpretation remain unresolved. See the
[accepted closure](evidence/e8-4-fix5-8-left-lowe-runtime-closure-2026-10-03.md).

Fix.5.5 presentation source `df957a6` passed independent actual-diff and
pushed-state review per the Fix.5.6 contract; full identity is in CURRENT and its
phase record. The later Fix.5.8 Left/lowe structural/visual farm gate is accepted;
it does not supply broader numerical interpretation, other-setting acceptance
or production promotion.
An analysis subprocess log contains only its child stream, not later owner
artifact/collection/ZIP failures. Absence of a gate-failure message there is
not evidence of owner success. Fix.5.6 persists atomic per-attempt owner status
outside the unchanged bundle and runs the unchanged detached collector's source
checks before analysis, retaining final collection/rechecks and ZIP verification.

A completed bundle is insufficient by itself: inspect provenance, checker
gates, structured payloads or logs, and rendered pages as applicable. When a
checker disagrees with raw evidence or implementation behavior, inspect the raw
evidence and implementation rather than forcing data to fit the checker. Use
one narrow gate -> one targeted farm run -> fresh evidence -> inspect -> one
coherent repair. See [LEARNINGS.md](LEARNINGS.md) for generalized failure and
review lessons.

Coherent implementation steps may proceed through deterministic local checks
and source review without a farm run after every edit. Farm validation remains
mandatory for ROOT/PyROOT or full-runtime claims and is scheduled at explicit
milestones; inspect fresh artifacts before broadening a narrow runtime gate.

Before a farm handoff, audit input authority, producer (if any), checker (if
any), collector, invocation owner, and returned artifact. Every required
executable step must be tracked, deterministic locally where possible,
independently source reviewed, pushed, and pushed-state reviewed. Unowned
multi-step orchestration blocks farm readiness; one reviewed CLI may own a
complete direct operation. After substantial implementation, review, closure,
or pushed-state handoff, Codex and ChatGPT report memory health. Unresolved
warnings block progress only for material active-state/provenance ambiguity
or a hard integrity violation. Record nonblocking warnings and batch their
repair at the next checkpoint or milestone audit; hard failures remain blocking.

For source-changing implementation, the default loop is memory checkpoint ->
substantive science/implementation ->
deterministic local test -> ChatGPT actual-diff review -> user commit/push ->
ChatGPT pushed-state synchronization review -> narrow farm/runtime validation
when required -> fresh evidence review -> memory checkpoint. Commit/push and
pushed-state review are required synchronization stages, not independent
scientific phases. A matching pushed candidate with materially accurate CURRENT
proceeds directly to the substantive gate, without a new reconciliation contract
merely to restate the push.

A memory inconsistency blocks only when it creates concrete ambiguity about
current source identity, accepted evidence, frozen scientific interfaces, active
scientific ownership, or the exact next operation. Otherwise record it, continue
the active milestone, and batch correction at a meaningful checkpoint. Every
substantive scientific loop must produce a plot, yield table, validated numerical
comparison, accepted runtime artifact, or one directly evidenced scientific/
runtime blocker with one coherent repair. Documentation alone is not scientific
completion. See the [workflow decision](decisions/memory-bracketed-scientific-throughput.md).

### Review discipline across projects

Prefer mature working source and the shortest science-producing path. A reviewed
diff and passing synthetic suites do not establish farm readiness. Trace CLI ->
preflight -> child -> validation/collection through actual farm path aliases,
symlink/cwd/import geometry, cheap-to-expensive check ordering and expected work
before approving a farm handoff. Expose unverified hard invariants and stop at
them; do not substitute repeated speculative infrastructure fixes for evidence.
The [seventh-failure record](investigations/kaonlt-matched-cohort-farm-readiness-failure-2026-10-22.md#seventh-failure-and-repeated-review-failures--october-23)
records repeated assistant review failures and the withdrawn direct-patch advice.
For source-changing implementation, preserve audit -> one standalone Markdown
contract -> short Codex launch ->
Codex implementation/local checks -> independent actual-diff review -> user
commit/push -> pushed-state synchronization -> farm only if proven ready.
Publication and farm execution remain user-only. Unreviewed source-changing
patch/application instructions bypass that workflow and are forbidden; this
does not prohibit authorized direct Work memory updates or contract-authorized
Codex source implementation.

The [ninth-attempt accountability checkpoint](investigations/kaonlt-matched-cohort-farm-readiness-failure-2026-10-22.md#ninth-failure-io-cost-and-assistant-accountability--october-24)
adds a durable review constraint: match verification cost to actual read/write
semantics, distinguish existence-only inputs from scientific byte consumers, and
expose whole-directory hash/copy scope and expected scale before readiness.
Diagnose a failed preservation/performance gate without blindly weakening it or
ordering another run to discover the next invariant; first justify the minimum
missing scientific evidence. Repeated premature readiness approvals and replacing
the tracked memory workflow with a writing block are assistant process failures.
Use the established three-check receipt; this lesson creates no support phase.

### Science priority and evidence admissibility

A real safety/provenance hard gate stops the specific unsafe operation; it
neither authorizes unrelated repair cycles nor automatically blocks independent
safe evidence collection. Anchor every support gate to the exact requested
fresh artifact. After two consecutive non-evidence supporting gates, stop new
support work and professionally re-anchor: name that artifact, the first
non-bypassable blocker and one minimal safe operation. Preserve failed gates;
no bypass or speculative farm command follows. Model algebra, historical plots,
synthetic tests, installed-package attestations and fresh farm signed-event
science are distinct evidence classes. The [October 15–17 incident](investigations/kaonlt-science-gate-displacement-2026-10-17.md)
records the assistant's failure to apply these existing priorities; it changes
no scientific ownership, acceptance or Method-A/Method-B role.

### Prevent redundant Work before authorization

The [November 3 accountability correction](investigations/kaonlt-assistant-roadmap-displacement-2026-11-03.md)
records ChatGPT's failure to apply already-read safeguards: settled choices and
Outcome B were reopened/repeated, mature SIMC machinery was underweighted, and
assistant-created studies displaced approved priorities. Before commissioning
Work, check implemented source, applicable accepted runtime gates, completed
investigation conclusions, the exact new evidence-producing outcome and its
dependency on the user-approved roadmap. Stop duplicates before spending resources.
An independently investigable topic is not automatically an approved next task;
priority alterations require concrete new evidence and user authorization.
Missing offline evidence is not missing production functionality. Preserve
scoped acceptance and the exact absolute-unit limitation; neither a conditional
formula nor archive omission demonstrates a needed normalization correction.
Reports, receipts and this accountability correction are not scientific milestones.
Enforce the existing displacement alarm prospectively; do not add another support
cycle to document its repeated violation. Responsibility is ChatGPT's planning
and review judgment, not the user's actions or Work's lack of farm access.

### Method-A detector-response scientific lessons

The [Gate-3 source/science audit](investigations/e8-4-left-lowe-method-a-detector-response-source-science-audit-2026-10-03.md)
preserves the evidence and limits behind these durable constraints:

- A relative detector-response representation is a different scientific object
  from an absolute leakage probability; accepted F.3 `hgcer3` is the former.
- Baseline `w0` already owns the established pion-control-to-background transfer.
  Multiplying another response/transfer object into it requires distinct
  scientific ownership and evidence that normalization/transfer physics is not
  double-counted; this is not an approved replacement algorithm.
- Signed parent-preserving normalization is not a positive probability-measure
  normalization. Where cancellation can matter, explicitly decompose signed
  and absolute support, source classes and response/correction tails. Possible
  amplification requires measurement, not an assumed explanation. The accepted
  Left/lowe t1 diagnostic measured only ~1.02% signed-versus-positive-support
  normalization contrast and supports no strong signed-normalization pathology.
- NPE=0 is the operational kaon PID category, NPE>0 the pion tree and NPE>2
  physical pion control. P(NPE=0 | true pion,x) is pion-to-kaon HGC mis-ID into
  the kaon-selected sample, not an observed pion class. Existing zero-photoelectron
  pion transfer machinery remains diagnostic-only and non-authoritative, with
  no production side effects or automatic replacement role.
- Slow-proton event probability -> proposed effect -> support/preservation
  diagnostic -> applied-effect decision is a structural analogue, not a
  numerical pion prescription; HGC-positive control alone cannot establish
  true-pion mis-ID into the kaon-selected sample.
- Accepted current-baseline comparator evidence reproduced F.2/F.3 scientific
  payloads exactly, with F.4 the first scientifically changed stage. This
  records floating-point-scale t1 F.4 differences versus substantive t0/t2
  changes. Accepted F.6.2 is acceptance-correlated refinement validation using
  normalized shapes, independent `SHMS_xptar`/`SHMS_yptar`, phi/MM and MM x
  acceptance localization, support/OOD, effective statistics, fixed kaon-window
  pion-change metrics and paired bootstrap intervals. That accepted evidence
  remains directly relevant to current t1 through this continuity; t0/t2 must
  not automatically inherit the exact t1 claim. No specific hardware cause or
  production promotion follows.
- Lack of an HGC-free true-pion tag in consumed F.1/F.3/F.4 limits direct absolute
  `P(NPE=0 | true pion,x)` calibration and unique weak-positive-to-zero-response
  proxy claims. It does not invalidate accepted relative acceptance-correlated
  refinement. RF can corroborate where experimentally available at low epsilon;
  high-epsilon data lack RF, so it is not a universal Method-A requirement. The
  accepted returned diagnostic did not perform RF corroboration.

### External PID/background-methods lessons

The [external review](investigations/e8-4-left-lowe-method-a-external-pid-background-methods-review-2026-10-03.md)
is EXTERNAL LITERATURE and planning inference, not KaonLT runtime evidence:

- Sample purification, conditional detector response, physics transfer and
  normalization/closure are distinct roles. Signed component-separation weights
  can be valid statistical weights without being physical event probabilities.
- Detector-response transport should use a physically meaningful positive
  population measure, with calibration/application support, composition,
  detector coordinates and run-condition matching checked. This design lesson
  does not approve replacing KaonLT's signed parent-preservation rule.
- The external review's absolute-proxy question is whether F.3 weak-positive pion response
  (0<NPE<=2 versus >2) is a valid relative proxy for pion-to-kaon HGC mis-ID
  topology across hgcer3. Independent HGC-free pion tagging is the preferred
  direct calibration if available. Zero-truncated/censored response models are
  optional model-dependent validation/cross-check machinery, not the foundation
  or a mandatory solution; Poisson remains a hypothesis requiring validation.
- Coupled cross-species PID may require a response/confusion matrix; its
  calibration, unidentified state and uncertainty propagation must be explicit.
  Literature does not mandate a KaonLT matrix or correction.
- External methods motivated separate F.4 signed-normalization measurements,
  not a proven defect or explanation of t1. The accepted current-lineage diagnostic
  supplies those measurements without establishing direct absolute proxy validity
  or warranting a replacement correction. Literature remains distinct from
  accepted KaonLT runtime evidence.

### Full-analysis runtime isolation

The ordinary full launcher intentionally mutates repository-local working
state through `git clean`, symlink setup, `set_sig_fortran` model rewrites and
generated model/cache files. Canonical validation must execute it only inside
an owner-created detached disposable analysis worktree under a sibling temporary
root, never in the ordinary farm checkout. Copy the actually imported ltsep
package to a separate PYTHONPATH overlay; uniquely match its PATH_TO_DIR config
against the ordinary baseline and change only copied LTANAPATH. Caller changes
alone do not override fixed config values. Prove overlay import identity and
equal launcher/set-sig/link-setup path dictionaries before execution. Preserve
external environment except the three debug/bypass variables; source background
config with resolved ltsep USER and require external links need no repair.
Record and recheck ordinary identity/status/farm-local bytes and original ltsep
path/config hashes even after child failure; never restore mismatches or copy
runtime mutations back. Cleanup is bounded to owner-created siblings. The
canonical-five owner/profile is a
local orchestration candidate, not runtime acceptance or production promotion.

## Canonical record ownership

- [CURRENT.md](CURRENT.md) owns the sole active objective, blockers, and exact
  next action; the handoff is exceptional transfer state only.
- `evidence/` owns accepted farm/runtime gate facts and artifacts; `phases/`
  owns implementation and fix chronology; `decisions/` owns scientific and
  architectural contracts.
- `roadmap/STATUS.md` owns approved dependency/status structure, not the
  exact active action. `LEARNINGS.md` owns generalized reusable lessons.
- `TOOLS.md`, `COMMUNICATION.md`, and `CODEX.md` own operational references,
  farm-delivery communication, and source-changing workflow respectively.

Use the authoritative [farm validation-bundle procedure](decisions/farm-validation-bundle-procedure.md)
for Jefferson Lab `tcsh` packaging rather than duplicating its procedure here.
