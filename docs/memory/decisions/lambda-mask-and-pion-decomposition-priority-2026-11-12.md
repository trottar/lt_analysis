# Lambda exclusion adoption and pion-decomposition priority — November 12, 2026

## Authority and status boundaries

Authority: the user's explicit November 12 memory-synchronization/scientific-priority
instruction. This records scientific decisions only. Public/local test at checkpoint:
`77c44eb60192f38c6979e71a2ae0cd8a920ada76`; pending source/memory changes remain
local. User adoption does not certify source review, publication, production
integration, final yield acceptance or farm/runtime validation. Those remain
separate; no phase is promoted here. CURRENT owns the single active NEXT.

## Adopted exclusion procedure and reviewed cohort mask

The user adopts the four Lambda (t,phi) criteria:

1. At least 30 total experimental entries.
2. At least 20 signal-region entries.
3. Signed cleaned signal-region yield at least the equal-width left-region yield.
4. Signed cleaned signal-region yield at least the equal-width right-region yield.

The reviewed diagnostic cohort retains 67/135 physical cells. Counts and signed
normalized yields are distinct. Its signal [1.114,1.126), left [1.102,1.114) and
right [1.126,1.138) GeV windows use the retained native experimental edges and
half-open inclusion. Total counts are retained kaon nommcuts source rows; signal
counts are allcuts rows in the signal interval, including signed-source streams
and zero-weight records. These are the reviewed definitions, not net significance.
Every failed criterion and empty/weak/nonpositive/finite-MC support state remains.
No extra negative-yield or poor-data/SIMC-agreement exclusion is adopted.

Preserve authoritative signed data, unchanged setting-level normalization,
no bin-level charge correction, and identical retained experimental/SIMC cells.
No renormalizing the surviving cells or replacing excluded measurements by
0 +/-10000. Preserve source/cohort/mask identity; a later correction must not
silently retune the mask against agreement or replace this reviewed support.
The 67-cell count describes this cohort, not guaranteed acceptance for other data.

Mask evidence: `kaonlt_lambda_bin_exclusion_20261112_review_bundle.zip`, SHA-256
`cf18308e93c9bc1444b2b7d73f139a97a2aff9e66c8bfa500e2a304d9182d9e5`,
member `report/primary_135_cell_decisions.csv`. The SHA-256 of the sorted retained
ID list encoded as compact JSON is
`504c22d76d71473ea7635c2b4e2e42b65c49edb2df707d179693b8078bfce09d`.
The review receipt includes the exact member hash and all 67 IDs. Existing support
losses and remaining MM discrepancies are retained, not erased by user adoption.
Production integration and farm validation of this mask are pending.

## Iteration-zero interpretation

Neither unity nor phi-flatness of data/SIMC ratios is an acceptance requirement
at iteration zero. Inspect these ratios diagnostically; do not reject an otherwise
supported initial result solely for nonunity or angular variation. The SIMC
cross-section parameterization is iterated subsequently with the initial pion
subtraction fixed. Do not tune background, area-match or alter normalization to
force either property. Preserve physical MM/yield requirements and uncertainties;
this is not a declaration of purity or final cross-section acceptance.
This supersedes earlier wording treating immediate constant-ratio agreement as
a prerequisite, without rewriting historical numerical results.

## All-15-parent MM presentation requirement

Every experimental/SIMC parent MM comparison uses identical bin edges. Never
mix coarse and fine histograms in the same comparison. Produce separate
unit-normalized shape and original yield-normalized comparisons; display shape
normalization never changes authoritative yields. Only exact bin aggregation or
source-event-backed binning can provide a common grid: no assumed fractional-bin
contents, interpolation or invented within-bin distribution. Report unsupported
mapping honestly; no new acquisition or replot is authorized by this memory task.
Earlier mixed-grid plots remain historical evidence, not fulfillment of this
requirement. No existing plot, histogram or analysis source is changed here.

## Immediate scientific priority and planned sequence

The fitted zero kaon-side SIDIS component is physically concerning and requires
an acceptance-informed assessment, particularly Center/high t3. Use the established
pion-control sample and existing neutron/Delta/SIDIS channel SIMC decomposition.
Examine accepted phase-space support, channel shapes/overlap, component degeneracies,
physically supported constraints and available uncertainties. Distinguish pion-control
amplitudes from kaon-side amplitudes and fitted zero from proven physical absence.
Do not impose a positive SIDIS lower bound without evidence or repeat completed
fit/bootstrap/model-response investigations as a new prerequisite.

Preserve accepted source-aware work, the 1024-replica conditional covariance and
116 model responses within their existing scopes. Complete physical covariance,
external PID/cross terms and the formal independently constrained Outcome B gate
remain unresolved, not zero or resolved by this decision. Neither Outcome B nor
unrelated normalization/Method-A gates displace this user's immediate assessment;
no replacement evidence-acquisition exercise is commissioned.

Following decomposition assessment and any separately authorized correction,
plan the established Local-t pion-subtraction chain -> (t,phi) binning -> adopted
exclusion mask -> Lambda spectra and yields. Preserve existing SIMC iteration,
L/T machinery and random/dummy/slow-proton ownership. Method A remains detached,
Method B diagnostic/excluded. This records the intended sequence, not present
permission to implement a correction or execute farm analysis.

## Review boundary

This memory-only update stops at independent ChatGPT actual-diff review. User-only
commit/push and subsequent published-state synchronization remain required. The
substantive priority remains the decomposition assessment, not another memory
reconciliation phase. Existing six roadmap statuses/dependencies and scoped
source/runtime closures are preserved. No physics source, Codex implementation,
farm run, refit, new bootstrap, commit or push occurs in this task.

## November 17 controlled population-correction comparison

The user's corrective instruction now requires independently executed OLD and
CORRECTED fits for all 15 parents through the same unchanged fitter. Only the
population constructor may differ; experimental/source events, all scientific
settings and procedures, alignment and primary D1.a default-off state are identical.
The complete binding requirements and parity/output obligations are in the
[D2.a governing record](../phases/event-population-authority-2026-11-16.md#november-17-binding-controlled-comparison-requirement).
Historical Global A/raw-pion B outputs and population/overlap plots do not answer
this controlled fitting question. Retain full signal-plus-background, pion-only
models and actual event-applied subtraction as distinct objects. This requirement
changes no mask, scientific model, roadmap status or source/runtime acceptance.

### November 17 correction: runtime location is not a source defect

Needing Jefferson Lab's ROOT environment does not establish a required code or
runner update. Existing acquired inputs, retained rows, published population
correction and shared fitter must be considered before implementation is scoped.
Demonstrate the exact missing existing-code capability separately; do not infer
it from absent Work PyROOT. The prepared matched-arm contract is DEFERRED and
unpublished. Its necessity and any future authorization must be independently
established. No new fit, farm update or acceptance is implied by this correction.

### November 17 correction: numerical fitting is distinct from ROOT execution

The user's NumPy clarification supersedes making an existing ROOT session the
scientific NEXT merely because the public fitter imports ROOT. Retained genuine
NumPy records and histogram arrays are numerical fitting inputs. Work must assess
whether the existing numerical routines/staged procedure can operate through an
array interface, preserving all mathematics, settings, signal protection, gates,
uncertainty treatment and failure policies between OLD and CORRECTED.

This does not establish that a ready equivalent array fitter exists. Verify the
actual reused code and every relevant histogram/coordinate/error operation; do
not substitute another statistical model, minimizer or historical amplitudes.
Record any specific irreducible ROOT operation separately. A literal import
failure is not evidence that offline fitting cannot be done. ROOT/main.py/farm
integration acceptance remains separate and still needs its own runtime evidence.
No production change, new farm contract or ROOT installation follows by inference.

### November 17 final D2.a assessment and hybrid investigation boundary

The user authorizes final scientific diagnostics A->B->C->D before Method A:
MM-boundary/protected-failure attribution, source-aware uncertainty/prediction,
parent-level STRICT/FORCED hybrid feasibility, then all15/native135 Lambda
peak/tail and fixed67-mask comparisons. This authorizes investigation only,
not adoption. Keep STRICT unchanged; never optimize choices by individual phi
or use the intended KLambda MC shape validator to define the hybrid selector.
Iteration-zero unit-normalized Lambda shape is a primary diagnostic, separately
from unchanged original yields. Later SIMC iterations retain the selected pion
subtraction; initial unity/phi-flatness and MC rescaling are not requirements.

The final diagnostic has now been executed. No independent/stable hybrid was
demonstrated; source rejections, conditional covariance limits, mode-specific
unavailable proposals and native variation remain. The [final assessment](../phases/event-population-authority-2026-11-16.md#november-17-final-scientific-assessment-before-method-a)
owns measured findings and exact input/evidence identities. Independent ChatGPT
review must decide whether a scientifically supported pion baseline can be
frozen before any detached Method A comparison. No Method A application, new
contract,acquisition,farm operation or production baseline promotion follows
automatically. CURRENT alone owns the actual next gate. All six milestone
statuses,Local-t,67/135,OutcomeB BLOCKED,D3-D5 DEFERRED and MethodA/B boundaries
remain unchanged.

### November 17 acceptance/constraint study and decision boundary

The user authorized Work-owned detached sensitivity calculations,not adoption.
Iteration-zero goals remain:determine pion contamination;freeze a subsequently
selected subtraction before SIMC cross-section iteration;prioritize Lambda peak/
relevant tail with separate original-yield and unit-shape views;require neither
initial data/MC unity nor phi-flatness. Assess variance against model dependence,
native yields and matched-epsilon response;no per-cell or MC-shape-only tuning.
Physical cross sections require actual matched accepted calibration/producer inputs.

The [completed measured study](../phases/event-population-authority-2026-11-16.md#november-17-completed-acceptance-and-model-dependence-sensitivity)
records the approved grid and actual150 central/19200 paired protected refits.
A continuous numerical plateau is established but no better physically demonstrated
full-cohort tradeoff is accepted. STRICT unchanged;FORCED diagnostic;no baseline,
threshold or hybrid adoption. Original67/135 mask,Local-t,six milestones,OutcomeB
BLOCKED,D3-D5 DEFERRED,MethodA detached and MethodB excluded remain. CURRENT owns
the sole NEXT:independent ChatGPT actual numerical/diff review then user baseline
decision. No new acquisition,farm,Codex or Method A task follows automatically.

### Corrected three-stage research sequence

The prior handoff incorrectly merged acceptance,model constraints and windows
and omitted acceptance-gate interactions. Previous findings retain their actual
configuration scope. ACTIVE Stage1 is the complete acceptance-threshold factorial;
independent numerical/diff review and user Stage2 authorization precede any
Lambda-constraint combinations;window combinations are a later separate stage.
No threshold,baseline,hybrid or Method A treatment is adopted. CURRENT owns NEXT;
see the existing D2.a phase for the exact source audit and declared grid. All
six roadmap milestones,Local-t,fixed67 mask,OutcomeB BLOCKED,D3-D5 DEFERRED,
MethodA detached and MethodB excluded retain their established scope/status.

Stage1 is now numerically complete:480 raw combinations with explicit mathematical/
empirical redundancy;source-owned gauge softening and protected gates verified.
Native Center/lowt3 responses remain heterogeneous despite parent shape improvement;
no policy or baseline adopted. Previous findings retain actual tested-policy scope.
The existing D2.a phase owns measurements. Independent Stage1 evidence/diff review
is NEXT;user authorization is required before Stage2. Constraint/window research,
Method A and production adoption do not follow automatically.

### November 18 reviewed Stage 1 and authorized Stage 2

Independent ChatGPT Stage1 scientific/published-memory review is complete within
its numerical scope. The user explicitly authorizes Work-owned Stage2 Lambda-
constraint combinations next;no acceptance policy or production baseline adopted.
The existing D2.a phase owns exact54-case scope and results. F120 remains primary;
a single frozen F280 secondary diagnostic is not cross-family optimization.
Stage3 windows and Stage4 cross-family work remain separately reviewed/unauthorized.
CURRENT owns NEXT;stop after Stage2 and return scientific/memory bundles for
independent review. Preserve fixed67/135,Local-t,six milestones and all scientific,
Method A/B and source/runtime boundaries. User alone publishes accepted changes.

Stage2 is now numerically complete on all54 combinations/all128 paired draws.
A Gaussian/retention.90 stability plane exists;retention/mode error reductions
carry conditional-mean and source-availability dependence. No improved total
physical tradeoff,baseline,hybrid,normalization or Method A adoption is established.
The existing D2.a phase owns measurements. NEXT is independent ChatGPT Stage2
numerical/diff review before the next scientific decision. Stage3 windows and
Stage4 cross-family optimization remain separately authorized future scopes.

### November 18 Stage2 review complete;Stage3 authorized

Independent ChatGPT Stage2 scientific/published-memory review is complete within
its offline scope;no physical baseline adopted. User authorizes Stage3 protected
MM-window/boundary factorial only. Existing D2.a phase owns evidence/scope. L26,
F120 primary/F280 separate fixed diagnostic,D1.a OFF,fixed67 mask and scientific
inputs/source stay unchanged. CURRENT owns NEXT:Stage3 then independent review.
STOP after Stage3;Stage4,Method A,source changes,Codex,farm and publication are
not authorized. Preserve all existing scientific decisions and six milestones.

Stage3 is now numerically complete:12 protected-MM windows,128 paired source
draws and32 heldout partitions perparent. Near-full-window numerical stability
is scoped;larger truncations change gate/fallback and native yield without
demonstrated total physical improvement. No baseline adopted. Existing D2.a
phase owns exact results;CURRENT NEXT is independent ChatGPT Stage3 science/diff
review before user decision. STOP after Stage3;no automatic Stage4 or Method A.
All scientific/mask/roadmap/source/runtime and publication boundaries preserved.
