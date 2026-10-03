# KaonLT — external PID/background-methods research memory checkpoint

## 1. Task class and objective

This is a **memory-only scientific-research reconciliation task** under
`docs/memory/CODEX.md`.

It follows the completed Gate-4 memory reconciliation and the user-authorized
scientific-direction decision to perform a deep external review before any
implementation contract.

Objectives:

1. preserve the external literature/methods review in one detailed repository
   investigation record with stable bibliographic provenance;
2. preserve only durable cross-phase lessons in `MEMORY.md`;
3. advance `CURRENT.md` and the roadmap to one precise diagnostics-first NEXT.

This task does **not** implement a pion correction, probability model,
normalization, transfer factor, response matrix, or production path. It does
not choose one external method as the KaonLT solution. It does not run the farm
or alter scientific source.

## 2. Exact starting identity

Required branch:

```text
test
```

Required starting local `HEAD` and local `origin/test`:

```text
f7934f459b42b190235ffc63e486e7ba0e068e79
```

Before editing, report:

```bash
git branch --show-current
git rev-parse HEAD
git rev-parse origin/test
git status --short --untracked-files=all
```

Hard stop on any identity mismatch or unexpected unrelated worktree state that
cannot be preserved safely. Do not reset, clean, stash, overwrite, commit,
push, update refs, or run the farm.

## 3. Mandatory startup reads

Read in exact repository order:

1. `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read:

- `docs/memory/MAINTENANCE.md`
- `docs/memory/CODEX.md`

Then task-relevant records:

- `docs/memory/investigations/e8-4-left-lowe-method-a-detector-response-source-science-audit-2026-10-03.md`
- `docs/memory/phases/e8-4-left-lowe-method-a-gate4-final-memory-reconciliation-task-contract.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/sources/SOURCE_INDEX.md`
- `docs/memory/evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md`
- `docs/memory/evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md`

External literature informs design but cannot override actual KaonLT source or
accepted farm evidence.

## 4. KaonLT baseline to preserve

### 4.1 Populations

```text
operational kaon PID category:
  P_hgcer_npeSum == 0 (not a true-species label)

pion tree:
  P_hgcer_npeSum > 0

training:
  prompt / noRF / nommcuts / P_hgcer_npeSum > 0

training weak-positive response:
  0 < P_hgcer_npeSum <= 2

training control-response:
  P_hgcer_npeSum > 2

physical application:
  authoritative physical-pion-control population with P_hgcer_npeSum > 2
```

NPE=0 is the operational kaon PID category, not a standalone observed pion
population/class. `P(NPE=0 | true pion,x)` is specifically the pion-to-kaon HGC
misidentification probability into the kaon-selected sample. Positive-response
training alone cannot directly calibrate it.

### 4.2 F.3

Accepted basis:

```text
hgcer3 = (SHMS_delta, P_hgcer_xAtCer, P_hgcer_yAtCer)
```

F.3 is a support-aware **relative response representation**, not an absolute
leakage probability, zero-photoelectron probability, detector efficiency, or
production correction.

### 4.3 F.4

Within each canonical-`t` parent:

```text
b_j = signed_source_coefficient_j * w0_j
r_j = relative Method-A response
B   = sum_j b_j
U   = sum_j b_j * r_j
C_j = r_j / (U / B)
```

thus:

```text
sum_j b_j * C_j = B
```

This preserves the signed parent, not each child/MM region/Lambda yield/event.

### 4.4 `w0`

`w0` already owns the accepted pion-control -> kaon-background transfer.
No detector-response object may silently duplicate this transfer.

### 4.5 Unresolved issue

The Left/lowe low-`t` / `t1` redistribution is real in the accepted private
branch, but its physical detector-response origin is not established.
Signed-normalization amplification remains `INFERENCE`, not an accepted cause.
Absolute-SIMC provenance remains a separate blocker.

## 5. External source hierarchy

### Tier A — peer-reviewed primary sources

#### A1. LHCb PID calibration

R. Aaij et al., “Selection and processing of calibration samples to measure the
particle identification performance of the LHCb experiment in Run 2,”
*EPJ Techniques and Instrumentation* **6**, 1 (2019).

- DOI: `10.1140/epjti/s40485-019-0050-z`
- arXiv: `1803.00824`

Supported points:

- probe identity established without the PID observable being calibrated;
- residual background may be removed using signed sWeights;
- PID response depends on kinematics, occupancy and detector conditions;
- calibration samples are reweighted separately to match reference populations;
- sample-purification weights and PID-efficiency transport are distinct;
- calibration/reference mismatch is an analysis-specific systematic.

#### A2. COMPASS RICH purity-efficiency matrix

G. D. Alexeev et al. (COMPASS), “Multiplicities of positive and negative
pions, kaons, and unidentified hadrons from deep-inelastic scattering of muons
off a liquid hydrogen target,” *Phys. Rev. D* **112**, 012002 (2025).

- DOI: `10.1103/q4rb-bhcg`
- arXiv: `2410.12005`
- CERN: `https://cds.cern.ch/record/2914029/files/document.pdf`

Supported points:

- `P_ij` diagonal = ID efficiencies, off-diagonal = mis-ID probabilities;
- data control samples: `K0_S -> pi pi`, `phi -> K K`, `Lambda -> p pi`;
- response depends primarily on momentum and RICH entrance angle;
- about 30 matrices form a 2D response grid;
- true yields recovered by matrix inversion;
- PID/mis-ID uncertainty propagated through matrix-element variations.

#### A3. HERMES RICH/unfolding

B. Hommez, “Hadron identification with the HERMES RICH,”
*Nucl. Instrum. Meth. A* **502**, 294–299 (2003).

- DOI: `10.1016/S0168-9002(03)00291-2`

A. Airapetian et al. (HERMES), “Multiplicities of charged pions and kaons from
semi-inclusive deep-inelastic scattering by the proton and the deuteron,”
*Phys. Rev. D* **87**, 074029 (2013).

- DOI: `10.1103/PhysRevD.87.074029`
- arXiv: `1212.5407`

Supported points:

- cross-species PID represented through response/probability matrices;
- true `(pi,K,p)` populations unfolded using inverse matrices;
- response can depend on momentum and event topology;
- unidentified state considered before truncation in multiplicity analysis;
- RICH simulation included/tuned mirror reflectivity, PMT response and mirror
  roughness; decay samples constrained the model.

#### A4. sPlot

M. Pivk and F. R. Le Diberder, “sPlot: a statistical tool to unfold data
distributions,” *Nucl. Instrum. Meth. A* **555**, 356–369 (2005).

- DOI: `10.1016/j.nima.2005.08.106`
- arXiv: `physics/0402083`

Signed event weights can be legitimate component-separation weights under the
method assumptions; they are not physical event probabilities.

#### A5. COWs

H. Dembinski, M. Kenzie, C. Langenbruch, M. Schmelling,
“Custom Orthogonal Weight functions (COWs) for event classification,”
*Nucl. Instrum. Meth. A* **1040**, 167270 (2022).

- DOI: `10.1016/j.nima.2022.167270`

Supported points:

- generalizes sWeights to broader discriminating/control relationships;
- component-separation weights and nonuniform efficiencies remain distinct;
- covariance treatment matters for weighted inference.

#### A6. Fake Factor / Matrix Method

K. Lehmann and B. Stelzer, “The Fake Factor Method and its relation to the
Matrix Method,” *Nucl. Instrum. Meth. A* **1054**, 168376 (2023).

- DOI: `10.1016/j.nima.2023.168376`

Supported point: data-driven control-to-signal transfer for misidentified objects
can be built from explicit efficiencies/matrix relations rather than an arbitrary
shape normalization.

#### A7. ATLAS Universal Fake Factor

ATLAS Collaboration, “Estimation of backgrounds from jets misidentified as
tau-leptons using the Universal Fake Factor method with the ATLAS detector,”
*Eur. Phys. J. C* **85**, 1441 (2025).

- DOI: `10.1140/epjc/s10052-025-14916-1`
- arXiv: `2502.04156`

Supported points:

- fake factors measured in dedicated data regions;
- distinct fake sources have distinct mis-ID probabilities;
- target-region factor combines source-specific factors with source composition;
- composition mismatch is an explicit systematic;
- validation regions are required.

#### A8. PMT-response caution

“On estimating the photoelectron yield and the resultant inefficiency of a
photomultiplier-based detector,” *Nucl. Instrum. Meth.* **225**, 153–163 (1984).

- DOI: `10.1016/0167-5087(84)91352-8`

Measured pulse-height response need not be a simple Poisson distribution in
photoelectron number; response spreading can affect inferred inefficiency.

E. H. Bellamy et al., “Absolute calibration and monitoring of a spectrometric
channel using a photomultiplier,” *Nucl. Instrum. Meth. A* **339**, 468–476
(1994).

- DOI: `10.1016/0168-9002(94)90183-X`

Realistic PMT response functions can be required for quantitative inference.

### Tier B — direct Jefferson Lab / Hall C detector analogues

#### B1. SHMS HGC commissioning

G. R. Ambrose, “Blinded by the Light: Commissioning of the Hall C SHMS Heavy
Gas Cherenkov Detector,” M.Sc. thesis, University of Regina (2018).

- JLab: `JLAB-PHY-19-2849`
- `https://misportal.jlab.org/ul/publications/view_pub.cfm?pub_id=15789`

Supported points:

- HGC photoelectron response explicitly characterized;
- response/efficiency spatially nonuniform;
- central/mirror-boundary geometry is important;
- efficiency/contamination can use independently identified species;
- sub-threshold response can arise from detector processes.

#### B2. SHMS spectrometer performance

S. F. Ali et al., “The SHMS 11 GeV/c spectrometer in Hall C at Jefferson Lab,”
*Nucl. Instrum. Meth. A* **1083**, 171070 (2026).

- DOI: `10.1016/j.nima.2025.171070`
- arXiv: `2503.08706`

Supported points:

- HGC performance is reported through efficiency and contamination for
  independently selected species;
- detector geometry/threshold response is distinct from physics-background
  normalization.

#### B3. CLAS12 HTCC spatial efficiency

V. Klimenko, Ph.D. thesis, section “Cherenkov Counter Efficiency Correction.”

- `https://misportal.jlab.org/sti/publications/22810/attachments/12991/KlimenkoPhD.pdf`

Supported points:

- HTCC NPE response mapped versus mirror-impact `(x,y)`;
- NPE spectra fit bin-by-bin;
- threshold efficiency inferred from response model;
- data and MC efficiencies treated separately;
- correction = data/MC detector-efficiency ratio;
- spatial nonuniformity associated with mirror/PMT geometry effects.

Close detector-level analogue to HGC `xAtCer,yAtCer`, but not a direct solution
to KaonLT's pion-to-kaon HGC mis-ID calibration or F.3 proxy-validity question.

#### B4. CLAS12 complex kaon sWeights

CLAS12 Note 2026-003, “Kaon ID with Two-dimensional Sequential S-weights.”

- `https://misportal.jlab.org/mis/physics/clas12/viewFile.cfm/2026-003.pdf?documentId=186`

Supports use of statistical component weights for complex kaon backgrounds;
does not turn such weights into detector-efficiency probabilities.

### Tier C — optional truncated-response statistical cross-check references

UCLA Statistical Consulting Group:

- `https://stats.oarc.ucla.edu/sas/dae/zero-truncated-poisson-regression/`
- `https://stats.oarc.ucla.edu/sas/dae/zero-truncated-negative-binomial/`

Supported points:

- if zeros cannot be observed, likelihood must condition on `Y>0`;
- ordinary Poisson/NB is inappropriate simply because only positive counts are
  observed;
- overdispersion must be assessed before Poisson is accepted.

These are methodological references for optional model-dependent cross-checks,
not detector validation or the foundation of the KaonLT response question.

### Tier D — project chronology/context only

Public KaonLT/Hall-C Redmine notes may be cited only as project chronology:
HGC hole studies, detector-efficiency studies, complementary-PID discussions,
and HGC-cut corrections. They are lower authority than peer-reviewed sources
and current repository evidence.

## 6. External-review findings to preserve

Use an explicit `EXTERNAL LITERATURE` label in the new investigation, separate
from KaonLT `SOURCE VERIFIED` and `RUNTIME VERIFIED`.

### 6.1 Distinct scientific roles

Across LHCb, COMPASS, HERMES, CLAS12 and Hall C, detector efficiency or
misidentification response is its own calibrated object.

Keep distinct:

```text
sample/component purification
detector response / PID probability
physics control-to-background transfer
normalization / closure constraint
```

For KaonLT, `w0` already owns physics pion-control -> kaon-background transfer.

### 6.2 Independent tagging is preferred

Strong experimental methods establish species identity independently of the PID
observable being calibrated. Independent HGC-free pion tagging is the preferred
direct KaonLT calibration if available. HGC-selected weak-positive versus >2
response is a candidate relative mis-ID topology proxy whose validity must be
tested, not a reason to mandate truncation/selection modeling.

### 6.3 Detector coordinates matter

Response is parameterized in detector-driving variables: momentum, entrance
angle/position, occupancy/topology, and run conditions. This supports the
physical relevance of `hgcer3`, especially `xAtCer,yAtCer`, but does not validate
current F.3/F.4 numerics.

### 6.4 Signed weights are not probabilities

LHCb provides a direct separation:

```text
signed sWeights -> calibration-sample purification
reference/calibration population weights -> response transport
PID efficiency -> conditional detector-response object
```

This is highly relevant because F.4 normalizes over
`b_j = signed_source_coefficient_j * w0_j`.

### 6.5 Response/confusion matrices are mature alternatives

COMPASS/HERMES model diagonal efficiencies, off-diagonal mis-ID, and sometimes
an unidentified state. This is a mature alternative representation when
cross-species leakage is coupled, but does not mandate a KaonLT matrix.

### 6.6 Weak-positive response as a candidate relative mis-ID topology proxy

The research question is whether F.3 weak-positive pion response (`0<NPE<=2`
versus `NPE>2`) is a valid relative proxy for pion-to-kaon HGC mis-ID topology
across hgcer3. NPE=0 is the kaon-selected category; conditional on true pions,
`P(NPE=0 | true pion,x)` is their HGC mis-ID probability into that sample.
The positive-response classes alone do not identify an absolute mis-ID
probability or establish proxy validity. Independent HGC-free pion tagging is
the preferred direct calibration if available. Zero-truncated/censored models
are optional model-dependent validation/cross-check machinery, not the
conceptual foundation or a mandatory solution.

### 6.7 Poisson is a hypothesis

CLAS12/Hall-C studies make Poisson-like Cherenkov response plausible, but PMT
literature shows measured response can depart from naïve Poisson. Future
optional truncated/censored modeling must compare response families and perform
predictive checks; this does not require such modeling.

### 6.8 Source composition matters

ATLAS demonstrates transfer-factor dependence on source composition. KaonLT
prompt/random/dummy source classes and signs cannot be hidden inside one global
factor without explicit composition/support checks.

### 6.9 F.4 signed normalization is a priority question, not a proven defect

External methods generally separate component subtraction weights from positive
detector-response probabilities. This makes F.4 signed normalization a priority
diagnostic target.

It does **not** prove that F.4 is wrong, that signed cancellation causes `t1`,
or that positive normalization is automatically correct. Label this `INFERENCE`.

## 7. KaonLT conceptual decomposition from the research

Preserve as a research result, **not an approved implementation**:

```text
A. component/sample identity or purification
B. weak-positive HGC response and its candidate relative pion-to-kaon mis-ID proxy
C. accepted pion-control -> kaon-background transfer w0
D. normalization / closure / physics-preservation rule
```

The present Method-A path couples B and D through F.4 after consuming C.
The next task must separately test B proxy validity and D signed-normalization
sensitivity before choosing a replacement.

## 8. Diagnostics-first NEXT

After this memory checkpoint passes independent actual-diff review, user
commit/push, and pushed-state synchronization, the sole substantive NEXT must be:

```text
NEXT — define a detached current-lineage diagnostic-measurement contract that
separately tests (1) F.3 weak-positive response as a relative pion-to-kaon HGC
mis-ID topology proxy across hgcer3 and (2) Left/lowe t1 F.4 signed-normalization
sensitivity, before any correction redesign or production implementation.
```

The future diagnostic should separately test:

1. weak-positive response versus relative pion-to-kaon HGC mis-ID topology
   across hgcer3, preferably calibrated with independent HGC-free pion tagging
   if available, with population matching and support/OOD checked;
2. F.4 signed-normalization sensitivity through at least these measurements:

```text
B = sum b_j
sum |b_j|
|B| / sum |b_j|
U = sum b_j r_j
U / B

source-separated authoritative classes

distributions:
  raw r_j
  final C_j
  w0_j
  signed b_j
  absolute |b_j|

dependence:
  canonical t
  phi child
  MM region
  SHMS_delta
  P_hgcer_xAtCer
  P_hgcer_yAtCer
  F.3 support / OOD

comparisons:
  positive-measure versus signed-normalization summaries
  with authoritative source classes and w0 ownership preserved
```

It should also ask whether an independently tagged pion sample exists in current
KaonLT/Hall-C data that can calibrate HGC response without using HGC itself.

Existing detached zero-photoelectron transfer diagnostics remain non-authoritative
cross-checks. Truncated/censored response models are optional model-dependent
validation only, not required for either primary test. Unavailable HGC-free tags
leave direct calibration NOT VERIFIED rather than mandating a model.

No correction choice is approved here.

## 9. Durable MEMORY lessons

Add only compact reusable lessons:

1. sample purification, detector response, physics transfer, and
   normalization/closure are distinct roles;
2. signed component-separation weights can be valid statistical weights without
   being probabilities;
3. detector response should be transported on a physically meaningful positive
   population measure with support/composition checked;
4. NPE=0 is the kaon category; true-pion zero response is pion-to-kaon HGC
   mis-ID, with HGC-free tagging preferred for direct calibration if available;
   truncated/censored modeling is an optional model-dependent cross-check;
5. cross-species PID may require a response/confusion matrix;
6. detector coordinates and calibration/application population matching matter;
7. pure Poisson is a detector-response hypothesis requiring validation;
8. separately test F.3 weak-positive proxy validity and F.4 signed-normalization
   sensitivity; neither proxy validity nor a normalization cause is established.

Link to the detailed investigation; do not duplicate the bibliography in MEMORY.

## 10. New investigation record

Create:

```text
docs/memory/investigations/e8-4-left-lowe-method-a-external-pid-background-methods-review-2026-10-03.md
```

It must contain:

- KaonLT question and frozen boundaries;
- source hierarchy;
- bibliographic table with citation, DOI/arXiv/stable URL, methodological point,
  transferability, and limitation;
- detailed Hall C/SHMS, CLAS12, COMPASS, HERMES, LHCb, sPlot/COWs,
  fake-factor/matrix, zero-truncated, and PMT-response sections;
- comparison table against current KaonLT architecture;
- labels:
  `SOURCE VERIFIED`, `RUNTIME VERIFIED`, `EXTERNAL LITERATURE`, `INFERENCE`,
  `NOT VERIFIED`;
- diagnostics-first NEXT.

Never label an external paper `RUNTIME VERIFIED`.

## 11. SOURCE_INDEX update

Update:

```text
docs/memory/sources/SOURCE_INDEX.md
```

Add a compact `External PID/background-methods research — 2026-10-03` section
that points to the investigation and lists stable primary bibliographic anchors
grouped by:

- Hall C / JLab detector sources;
- RICH/PID response-matrix sources;
- calibration/transport sources;
- statistical weighting/transfer sources;
- truncated/PMT-response references.

The source index remains provenance/navigation only and owns no status/NEXT.

## 12. CURRENT update

Compact rather than append.

It must state:

- Gate 4 consumed;
- external review explicitly user-authorized and completed enough to alter
  planning;
- no implementation occurred;
- strongest external lesson is separation of purification / detector response /
  `w0` transfer / normalization-preservation;
- no external method validates a specific KaonLT correction;
- F.4 signed normalization is a priority diagnostic question, not proven defect;
- NPE=0 is the operational kaon category, NPE>0 the pion tree, NPE>2 physical
  pion control; conditional true-pion zero response is pion-to-kaon HGC mis-ID;
- F.3 weak-positive proxy validity is the primary response question; HGC-free
  tags are preferred if available, truncated/censored modeling optional;
- no detector hardware cause accepted;
- Method A detached/non-production;
- Method B diagnostic only/numerically excluded;
- canonical-five `DEFERRED`;
- final E.8/F.6.4 `BLOCKED`;
- absolute-SIMC separately `BLOCKED`.

Sole NEXT:

```text
NEXT — after independent review, user commit/push and pushed-state
synchronization: define a detached current-lineage diagnostic-measurement
contract with separate tests of (1) F.3 weak-positive response as a relative
pion-to-kaon mis-ID topology proxy across hgcer3 and (2) Left/lowe t1 F.4
signed-normalization sensitivity, before correction redesign or production.
```

## 13. Roadmap update

Minimally update dependency to:

```text
accepted Left/lowe evidence
-> Gates 1–4 consumed
-> external analogous-experiment / PID-background methods review consumed
-> detached current-lineage diagnostic-measurement contract
-> evidence-based method-selection decision
-> implementation contract only if warranted
-> later canonical-five reconsideration
-> final E.8
-> F.6.4 explicit production-promotion decision
```

Preserve all frozen statuses.

## 14. Allowed versioned paths

Only:

```text
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/roadmap/STATUS.md
docs/memory/sources/SOURCE_INDEX.md
docs/memory/manifest.json
docs/memory/investigations/e8-4-left-lowe-method-a-external-pid-background-methods-review-2026-10-03.md
docs/memory/phases/e8-4-left-lowe-method-a-external-pid-background-methods-research-memory-task-contract.md
```

The final two are intended new tracked files.

`VALIDATION_HISTORY.md` is not allowlisted: literature is not farm evidence.

## 15. Frozen paths

Everything else is frozen, especially:

```text
src/
testing/
tools/
farm_env/
run_Prod_Analysis.sh
AGENTS.md
docs/memory/USER.md
docs/memory/CODEX.md
docs/memory/MAINTENANCE.md
docs/memory/COMMUNICATION.md
docs/memory/handoffs/CURRENT_HANDOFF.md
docs/memory/evidence/
docs/memory/decisions/
all prior phase contracts
all prior investigations
```

## 16. Positive checks

Confirm:

```text
exact starting HEAD/origin/test
only seven allowlisted paths changed
Gate 4 consumed
external review recorded with stable sources
SOURCE VERIFIED distinct from EXTERNAL LITERATURE
no external paper labeled RUNTIME VERIFIED
signed weights not called probabilities
w0 ownership preserved
relative response distinct from absolute leakage probability
operational PID categories and conditional true-pion mis-ID distinguished
weak-positive proxy validity separated from F.4 normalization sensitivity
HGC-free tagging preferred if available; response modeling optional
Poisson treated as model hypothesis
F.4 signed normalization diagnostic priority, not proven defect
Method A detached
Method B numerically excluded
canonical-five DEFERRED
final E.8/F.6.4 BLOCKED
absolute-SIMC separate
one diagnostics-first NEXT
```

## 17. Negative checks

Confirm:

```text
no src/testing/tool changes
no farm commands
no ROOT/PyROOT claims
no probability/response-matrix implementation
no w0/F.4 modification
no production correction
no child renormalization
no Method-B numerical use
no claim signed cancellation causes t1
no claim F.4 is wrong
no claim Poisson is correct HGC model
no NPE=0 standalone observed pion population/class
no mandatory truncated/censored solution
no PMT/mirror hardware cause claimed
no external paper promoted to KaonLT runtime evidence
no implementation contract
```

## 18. Memory health and manifest

Discover `<PYTHON>` through `docs/memory/TOOLS.md`.

After edits:

1. regenerate/check manifest;
2. run:
   `<PYTHON> -B tools/check_memory_health.py --root .`
3. run:
   `git diff --check`

Report:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
warning classification: <none or classifications>
manifest check: PASS | FAIL
git diff --check: PASS | FAIL
```

No runtime claim follows.

## 19. Actual-diff review bundle

Create temporary repository-root:

```text
kaonlt_review.diff
```

Include complete tracked diffs plus complete `git diff --no-index /dev/null ...`
additions for both intended new Markdown files, using the POSIX/Git-Bash procedure
in `CODEX.md`.

Do not stage for review. Stop before staging, commit, push, farm execution, or
implementation.

## 20. Acceptance criteria

Pass only if:

- exact identity and allowlist pass;
- literature is recorded accurately/conservatively;
- primary sources distinguished from theses/notes/method references;
- no external method is promoted into KaonLT;
- diagnostics-first conclusion is preserved;
- all frozen ownership remains unchanged;
- manifest/health/diff checks pass;
- complete review bundle produced;
- no farm/source implementation/commit/push occurred.

## 21. Hard stop

Return `BLOCKED` if:

- repository identity differs;
- unexpected worktree state cannot be preserved;
- current source contradicts Gate-3 ownership;
- a source claim cannot be supported by the supplied citation;
- external literature would have to be treated as KaonLT runtime evidence;
- Method-A/Method-B/baseline ownership would have to change;
- memory health has an unresolved hard failure.

Do not invent a workaround.

This contract records research and advances to a diagnostics-first scientific
gate. It does not authorize implementation.
