# External PID/background-methods review — 2026-10-03

## Question, authority and scope

The user explicitly authorized external research after Gate 4. This memory-only
checkpoint records that review under the [research contract](../phases/e8-4-left-lowe-method-a-external-pid-background-methods-research-memory-task-contract.md),
at branch `test`, observed HEAD/local `origin/test`
`f7934f459b42b190235ffc63e486e7ba0e068e79`.
Gates 1–4 are consumed. The question is how to diagnose the physical meaning of
the large Q4p4W2p74 / Left / lowe low-t/t1 Method-A redistribution before choosing
a correction. The primary response question is whether F.3 weak-positive pion
response (0<NPE<=2 versus >2) is a valid relative proxy for pion-to-kaon HGC
mis-ID topology across hgcer3, separately from F.4 signed-normalization sensitivity.
External methods inform design; they cannot override KaonLT source
or accepted farm evidence, establish a detector cause, or validate a correction.

No scientific source, probability model, response matrix, transfer, normalization,
correction or production path was implemented. No farm run occurred. This record
does not write the later diagnostic-measurement implementation contract.

## SOURCE VERIFIED — KaonLT baseline

The [Gate-3 investigation](e8-4-left-lowe-method-a-detector-response-source-science-audit-2026-10-03.md)
owns the detailed source audit. Targeted population/F.4 checks at this checkpoint
and the memory-only intervening commit range preserve its ownership. The repair
rechecked `PION_HGCER_PID_MASK_CONFIG` in
[`background_config.py`](../../../src/utility/background_config.py):

| Object | Current KaonLT meaning |
| --- | --- |
| Operational kaon PID category | `P_hgcer_npeSum == 0`; selected category, not a true-species label |
| Pion tree | `P_hgcer_npeSum > 0` |
| Training | prompt / noRF / nommcuts / `P_hgcer_npeSum > 0` |
| Weak-positive response class | `0 < P_hgcer_npeSum <= 2` |
| Control-response class | `P_hgcer_npeSum > 2` |
| Physical application | authoritative physical-pion-control population, `P_hgcer_npeSum > 2` |
| F.3 | support-aware relative response on `hgcer3 = (SHMS_delta, P_hgcer_xAtCer, P_hgcer_yAtCer)` |
| Baseline `w0` | accepted pion-control -> kaon-background physics transfer |
| F.4 | signed canonical-t parent preservation, without independent child normalization |

NPE=0 is the operational kaon PID category, not a standalone observed pion
population/class. Conditional on true species, `P(NPE=0 | true pion,x)` is
specifically the pion-to-kaon HGC misidentification probability into the
kaon-selected sample. Positive-response training cannot directly calibrate that
probability. F.3 weak-positive versus >2 response may provide a relative topology
proxy, whose validity remains NOT VERIFIED; it is not an absolute mis-ID
probability, detector efficiency or production correction. F.4 uses:

```text
b_j = signed_source_coefficient_j * w0_j
r_j = relative Method-A response
B = sum_j b_j
U = sum_j b_j r_j
C_j = r_j / (U / B)
sum_j b_j C_j = B
```

This preserves the signed parent, not each phi child, MM region, Lambda yield
or event contribution. `w0` transfer ownership cannot be silently duplicated
by a second detector-response object. Existing zero-photoelectron transfer
machinery remains detached, non-authoritative and diagnostic-only. The configured
HGCer hole is already excluded from relevant populations; no simple known-hole
correction explanation follows. Slow-proton architecture remains analogue-only.

## RUNTIME VERIFIED — existing KaonLT evidence only

The [current-baseline comparator](../evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md)
reproduced F.2/F.3 scientific payloads exactly, with F.4 the first changed stage.
The [F.6.3/E.8.4 evidence](../evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md)
supports only the accepted Left/lowe private branch, current-lineage application,
live-cache parity, real child changes and signed parent preservation. Observed
large low-t effects are real at that scope; their detector origin is not known.
F.4.Refresh.2 and earlier detached F-stage closures retain their recorded scopes;
Fix.5.7/Fix.5.8 retain narrow owner/checker/presentation closures.

These are accepted supplied farm-review records, not raw archives newly reopened
or fresh runtime validation in this task. No external paper is KaonLT runtime
evidence. The post-run model-file provenance caveat remains unchanged.

## EXTERNAL LITERATURE — hierarchy and verification scope

Tier A: peer-reviewed primary detector/analysis/method papers. Tier B: direct
JLab/Hall-C analogues, with peer-reviewed SHMS performance distinguished from
theses and collaboration notes. Tier C: institutional statistical-method
references; these are not detector validation. Tier D: public project chronology,
including Redmine, has context-only authority and is not used to validate a method.

Stable anchors below were checked through primary publisher, author manuscript,
institutional or laboratory records. Full-text LHCb, COMPASS, HERMES multiplicity,
SHMS and ATLAS material was available; publisher abstracts support the weighting,
fake-factor and PMT cautions. Some JLab PDF endpoints failed to reopen through
the web tool. Their detailed results are preserved as the contract-supplied
external review, explicitly identified below, rather than a new independent
full-text check. Metadata access is not a detector-validation claim.

### Stable bibliography

| ID/tier | Citation and stable anchor | Methodological point | Transferability to KaonLT | Limitation |
| --- | --- | --- | --- | --- |
| A1 | R. Aaij et al., “Selection and processing of calibration samples to measure the particle identification performance of the LHCb experiment in Run 2,” EPJ Techniques and Instrumentation 6, 1 (2019). [DOI 10.1140/epjti/s40485-019-0050-z](https://doi.org/10.1140/epjti/s40485-019-0050-z); [arXiv:1803.00824](https://arxiv.org/abs/1803.00824), Secs. 4–5 | Independent tagging, purification and response transport separated | Population matching and support checks | LHCb detector/calibration samples differ |
| A2 | G. D. Alexeev et al. (COMPASS), “Multiplicities of positive and negative pions, kaons, and unidentified hadrons from deep-inelastic scattering of muons off a liquid hydrogen target,” Phys. Rev. D 112, 012002 (2025). [DOI 10.1103/q4rb-bhcg](https://doi.org/10.1103/q4rb-bhcg); [arXiv:2410.12005](https://arxiv.org/abs/2410.12005); [CERN published PDF](https://cds.cern.ch/record/2914029/files/document.pdf), PID and systematic-uncertainty sections | Calibrated purity-efficiency matrices | Explicit efficiency/mis-ID roles | RICH multi-species calibration is not available merely from F.3 |
| A3a | B. Hommez, “Hadron identification with the HERMES RICH,” NIM A 502, 294–299 (2003). [DOI 10.1016/S0168-9002(03)00291-2](https://doi.org/10.1016/S0168-9002(03)00291-2); [author institution record](https://biblio.ugent.be/publication/316360) | RICH response/model characterization | Detector-model calibration context | Detailed optical-model findings below are supplied-review findings |
| A3b | A. Airapetian et al. (HERMES), “Multiplicities of charged pions and kaons from semi-inclusive deep-inelastic scattering by the proton and the deuteron,” Phys. Rev. D 87, 074029 (2013). [DOI 10.1103/PhysRevD.87.074029](https://doi.org/10.1103/PhysRevD.87.074029); [arXiv:1212.5407](https://arxiv.org/abs/1212.5407), Sec. III.3.2 | Species unfolding with unidentified state | Response-matrix alternative | MC calibration/model dependence; no mandatory KaonLT matrix |
| A4 | M. Pivk and F. R. Le Diberder, “sPlot: a statistical tool to unfold data distributions,” NIM A 555, 356–369 (2005). [DOI 10.1016/j.nima.2005.08.106](https://doi.org/10.1016/j.nima.2005.08.106); [arXiv:physics/0402083](https://arxiv.org/abs/physics/0402083) | Statistical component weights | Signed weights can be legitimate without being probabilities | Model and component-wise discriminant/control assumptions |
| A5 | H. Dembinski, M. Kenzie, C. Langenbruch, M. Schmelling, “Custom Orthogonal Weight functions (COWs) for event classification,” NIM A 1040, 167270 (2022). [DOI 10.1016/j.nima.2022.167270](https://doi.org/10.1016/j.nima.2022.167270); [arXiv:2112.04574](https://arxiv.org/abs/2112.04574) | Generalized component weighting and covariance | Purification/efficiency separation | Still requires an identified statistical model |
| A6 | K. Lehmann and B. Stelzer, “The Fake Factor Method and its relation to the Matrix Method,” NIM A 1054, 168376 (2023). [DOI 10.1016/j.nima.2023.168376](https://doi.org/10.1016/j.nima.2023.168376) | Explicit efficiency relations for fake transfer | Control-to-target transfer reasoning | Lepton fake definitions are not HGC pion leakage |
| A7 | ATLAS Collaboration, “Estimation of backgrounds from jets misidentified as tau-leptons using the Universal Fake Factor method with the ATLAS detector,” EPJ C 85, 1441 (2025). [DOI 10.1140/epjc/s10052-025-14916-1](https://doi.org/10.1140/epjc/s10052-025-14916-1); [arXiv:2502.04156](https://arxiv.org/abs/2502.04156), Sec. 5 | Source composition and dedicated-region validation | Composition/support diagnostics | Jet/tau sources cannot be mapped numerically onto prompt/random/dummy |
| A8a | “On estimating the photoelectron yield and the resultant inefficiency of a photomultiplier-based detector,” NIM 225, 153–163 (1984). [DOI 10.1016/0167-5087(84)91352-8](https://doi.org/10.1016/0167-5087(84)91352-8) | Pulse-height spreading affects inefficiency inference | Caution against naïve Poisson inference | No KaonLT response-family selection |
| A8b | E. H. Bellamy et al., “Absolute calibration and monitoring of a spectrometric channel using a photomultiplier,” NIM A 339, 468–476 (1994). [DOI 10.1016/0168-9002(94)90183-X](https://doi.org/10.1016/0168-9002(94)90183-X) | Realistic PMT response function | Quantitative response calibration | Instrument calibration is not a pion-background normalization |
| B1 | G. R. Ambrose, “Blinded by the Light: Commissioning of the Hall C SHMS Heavy Gas Cherenkov Detector,” M.Sc., University of Regina (2018), JLAB-PHY-19-2849. [JLab publication 15789](https://misportal.jlab.org/sti/publications/15789) | HGC NPE calibration and efficiency | Direct SHMS detector analogue | Thesis; detailed spatial findings supplied by review, PDF not reopened |
| B2 | S. F. Ali et al., “The SHMS 11 GeV/c spectrometer in Hall C at Jefferson Lab,” NIM A 1083, 171070 (2026). [DOI 10.1016/j.nima.2025.171070](https://doi.org/10.1016/j.nima.2025.171070); [arXiv:2503.08706](https://arxiv.org/abs/2503.08706), HGC performance/Table 5 | Independently selected species efficiency/contamination | Closest peer-reviewed detector context | Different settings/thresholds, no validation of current t1 |
| B3 | V. Klimenko, “Differential Cross Sections from CLAS12 RG-A Inclusive Electron Scattering,” Ph.D. thesis (2024), section “Cherenkov Counter Efficiency Correction,” JLAB-PHY-24-4212. [JLab PDF](https://misportal.jlab.org/sti/publications/22810/attachments/12991/KlimenkoPhD.pdf); [OSTI 2459311](https://www.osti.gov/biblio/2459311) | Spatial response/threshold efficiency modeling | Detector-coordinate diagnostic analogue | Electron HTCC, not pion-to-kaon HGC mis-ID calibration; detailed findings supplied-review only |
| B4 | A. Acar, M. Bashkanov, D. Watts, N. Zachariou, “Kaon ID with Two-dimensional Sequential S-weights,” CLAS12 Note 2026-003 (2026). [JLab note PDF](https://misportal.jlab.org/mis/physics/clas12/viewFile.cfm/2026-003.pdf?documentId=186) | Statistical kaon component separation | Complex-background weighting context | Collaboration note, not peer-reviewed HGC response evidence |
| C1 | UCLA Statistical Consulting Group, “Zero-Truncated Poisson Regression.” [Institutional method reference](https://stats.oarc.ucla.edu/sas/dae/zero-truncated-poisson-regression/) | Likelihood conditioned on positive counts | Optional model-dependent selection cross-check | Count-model example, not analog pulse-height detector validation |
| C2 | UCLA Statistical Consulting Group, “Zero-Truncated Negative Binomial.” [Institutional method reference](https://stats.oarc.ucla.edu/sas/dae/zero-truncated-negative-binomial/) | Truncation with overdispersion | Optional model-dependent response-family comparison | Does not establish HGC negative-binomial behavior |

## EXTERNAL LITERATURE — detailed methods and limits

### Hall C / SHMS

Ambrose B1 characterizes NPE response and detector efficiency; the supplied
review additionally identifies spatial nonuniformity, central/mirror boundaries,
independently selected species and sub-threshold detector response. B2 defines
HGC efficiency and contamination using non-HGC identification of clean species,
at specified thresholds and settings. These motivate detector-level observables
separate from physics-background normalization. They establish no mirror, PMT,
alignment or hardware cause for KaonLT t1. [B1](https://misportal.jlab.org/sti/publications/15789),
[B2](https://arxiv.org/html/2503.08706v2).

### CLAS12

The supplied B3 review describes HTCC NPE spectra fitted by mirror-impact `(x,y)`
bins, model-derived threshold efficiencies, separate data/MC efficiencies and
their ratio. Spatial response tied to optics/PMTs is a close coordinate analogue,
but the electron HTCC study cannot calibrate true-pion mis-ID into KaonLT
kaon-selected NPE=0 events or validate the F.3 relative topology proxy.
B4 addresses kaon backgrounds with sequential statistical weights, not detector
probabilities. B3's PDF was not independently reopened; B4's laboratory index
and indexed primary abstract confirm identity/context, not a fresh full-note audit.
[B3](https://misportal.jlab.org/sti/publications/22810/attachments/12991/KlimenkoPhD.pdf),
[B4](https://misportal.jlab.org/mis/physics/clas12/viewFile.cfm/2026-003.pdf?documentId=186).

### COMPASS

A2 uses decay samples `K0_S -> pi pi`, `phi -> K K`, `Lambda -> p pi` to
determine diagonal ID efficiencies and off-diagonal mis-ID probabilities.
About 30 matrices cover momentum and RICH entrance angle; inverse matrices
recover species yields, and variations of matrix elements propagate PID
uncertainty. This is an established representation when leakage is coupled,
not permission to construct a KaonLT response matrix from selected F.3 classes.
[A2, Sec. 3.2 and PID systematics](https://arxiv.org/pdf/2410.12005).

### HERMES

A3b relates true pi/K/p populations to identified species plus an unidentified
state; momentum, charge and event topology index its response. Removing the
unidentified row then inverting yields event-counting weights. The supplied
A3a review describes simulation tuning of mirror reflectivity, PMT response
and mirror roughness constrained by decay samples. Those detailed optical
claims were not independently reopened here. Neither matrix inversion nor
optical tuning is a KaonLT correction choice.
[A3b, Sec. III.3.2](https://arxiv.org/html/1212.5407),
[A3a](https://doi.org/10.1016/S0168-9002(03)00291-2).

### LHCb calibration and transport

A1 establishes probe species without using the calibrated PID requirement;
residual background is statistically subtracted. Detector response depends on
kinematics, occupancy and operating conditions. Calibration/reference transport
uses population matching separately from sample purification; mismatch remains
analysis-specific. The conceptual separation is signed purification weights,
reference/calibration population weights and conditional PID efficiency.
Transferability depends on adequate support and detector-variable coverage,
not simply a global parent normalization.
[A1, calibration samples and Sec. 5](https://arxiv.org/html/1803.00824).

### sPlot and COWs

A4 provides component-separation weights that may be signed; their interpretation
is statistical, not a physical event probability. Its model and component-wise
discriminant/control independence assumptions matter. A5 generalizes weighting
to broader modeled relationships and treats covariance of weighted inference.
These support keeping purification and efficiency distinct. They do not prove
that KaonLT subtraction coefficients are sWeights, justify ignoring correlations,
or select a replacement for F.4.
[A4](https://arxiv.org/abs/physics/0402083),
[A5](https://doi.org/10.1016/j.nima.2022.167270).

### Fake-factor / matrix transfer

A6 derives fake-factor relations from matrix-method efficiency relations for
misidentified leptons. A7 measures dedicated-region factors for different fake
sources, combines them according to target composition and evaluates composition
mismatch and validation tests. The useful lesson is explicit source composition
and support, not imported tau/lepton factors. Prompt/random/dummy signs and
classes must remain explicit when diagnosing KaonLT; statistical subtraction
classes are not automatically the physical fake categories used by ATLAS.
[A6](https://doi.org/10.1016/j.nima.2023.168376),
[A7, Sec. 5](https://arxiv.org/html/2502.04156v2).

### Optional truncated/censored and PMT-response cross-checks

C1/C2 describe likelihoods conditional on `Y>0`; a positive-only sample should
not be fitted as an untruncated Poisson/NB sample. Overdispersion needs assessment.
For KaonLT, `P(NPE=0 | true pion,x)` means pion-to-kaon HGC mis-ID into the
kaon-selected sample. Independent HGC-free pion tagging, if available, is the
preferred direct calibration: it can test mis-ID topology across hgcer3 without
selecting the tag on HGC response. The weak-positive versus >2 classes alone
do not identify an absolute mis-ID probability or establish proxy validity.
Truncated/censored models are optional model-dependent validation/cross-check
machinery for that topology question, not its conceptual foundation or a
mandatory solution.
[C1](https://stats.oarc.ucla.edu/sas/dae/zero-truncated-poisson-regression/),
[C2](https://stats.oarc.ucla.edu/sas/dae/zero-truncated-negative-binomial/).

A8a warns that measured pulse-height distributions can depart from naïve
photoelectron Poisson assumptions through response spreading; A8b motivates
realistic PMT response functions. Thus Poisson is a plausible starting hypothesis,
not an accepted HGC model. If optional modeling is pursued, candidate-family
comparisons and predictive checks must distinguish latent counts, measured NPE
response, truncation and threshold/censoring. No family or modeling requirement
is selected here.
[A8a](https://doi.org/10.1016/0167-5087(84)91352-8),
[A8b](https://doi.org/10.1016/0168-9002(94)90183-X).

## INFERENCE — comparison with KaonLT architecture

The research motivates this decomposition, not an approved implementation:

```text
A. component/sample identity or purification
B. weak-positive HGC response and its candidate relative pion-to-kaon mis-ID proxy
C. accepted pion-control -> kaon-background transfer w0
D. normalization / closure / physics-preservation rule
```

The current path couples B's relative representation and D through F.4 after
consuming C. The next gate must separately test B proxy validity and D signed-
normalization sensitivity before selecting a method.

| External role | KaonLT comparison | Diagnostic implication | Unestablished conclusion |
| --- | --- | --- | --- |
| Independent species tag | Training is already HGC-positive selected | Prefer HGC-free pion tagging, if available, for direct mis-ID calibration | No such sample is established |
| Conditional detector response | F.3 represents weak-positive versus >2 relative `hgcer3` response | Test its relative proxy validity against pion-to-kaon HGC mis-ID topology | Neither proxy validity nor absolute mis-ID is established |
| Positive population transport | F.4 sums signed coefficient times `w0` | Compare positive-measure and signed summaries with source composition explicit | No authorization to replace signed normalization |
| Purification weights | KaonLT prompt/random/dummy source coefficients | Preserve authoritative classes/signs and separate their contributions | They are not automatically sWeights or probabilities |
| Physics control-to-background transfer | Baseline `w0` already owns transfer | Inspect `w0` jointly with response/correction tails | No extra transfer may silently double-count it |
| Coupled species response matrix | No approved KaonLT matrix in this task | Consider representational requirements only after measurements | Literature does not mandate matrix implementation |
| Optional response likelihood | NPE=0 is kaon-selected; true-pion mis-ID is conditional on species | Model-dependent truncated/censored cross-check only if useful | No mandatory model, selected family or validated absolute mis-ID |
| Validation/closure | F.4 signed parent equality is enforced | Separate parent integrity from child/MM/yield sensitivity | Parent equality is not scientific correction validity |

Signed-normalization amplification is a priority hypothesis: cancellation in
`B` and `U` could increase sensitivity of `U/B`. External separation of statistical
weights and positive response probabilities motivates measuring it. It proves
neither that F.4 is wrong nor that signed cancellation causes t1; replacing it
with positive normalization would itself need a scientific preservation decision.

## NOT VERIFIED — open questions

No detector hardware cause is accepted. The exact current-lineage t1 decomposition,
relative pion-to-kaon mis-ID proxy validity, HGC-free tag availability, optional-model
predictive performance and production suitability remain unverified. Literature
cannot resolve the separate absolute-SIMC unit/provenance blocker. Independent
full-text reopening of B1/B3 and the detailed A3a/B4 supplied review remains an
access limitation, not new experimental evidence or a change to KaonLT ownership.

## Diagnostics-first conclusion and NEXT

NEXT — after independent actual-diff review, user commit/push and pushed-state
synchronization: define a detached current-lineage diagnostic-measurement
contract with separate tests of (1) F.3 weak-positive response as a relative
pion-to-kaon mis-ID topology proxy across hgcer3 and (2) Left/lowe t1 F.4
signed-normalization sensitivity, before correction redesign or production.
This record does not create that contract.

The later measurement should retain exact source/artifact/population identity
and authoritative source classes, with two distinct tests:

1. **Response/proxy validity:** compare F.3 weak-positive response (0<NPE<=2
   versus >2) with relative pion-to-kaon HGC mis-ID topology across hgcer3,
   including coordinates, population matching and support/OOD. Prefer direct
   calibration with independent HGC-free pion tagging if available; establish
   its availability and species purity rather than assume them. An unavailable
   tag leaves direct calibration NOT VERIFIED, not a mandate to fit a model.
2. **F.4 signed-normalization sensitivity:** measure current-lineage Left/lowe
   t1 signed/absolute support, source decomposition, response/correction tails
   and child/MM-region redistribution, separately from proxy validity:

```text
B = sum b_j
sum |b_j|
|B| / sum |b_j|
U = sum b_j r_j
U / B
source-separated authoritative classes
distributions: raw r_j, final C_j, w0_j, signed b_j, absolute |b_j|
dependence: canonical t, phi child, MM region, SHMS_delta,
            P_hgcer_xAtCer, P_hgcer_yAtCer, F.3 support/OOD
comparisons: positive-measure versus signed-normalization summaries
            with authoritative source classes and w0 ownership preserved
```

Existing detached zero-photoelectron transfer diagnostics may be inspected only
as non-authoritative cross-checks. Truncated/censored response models remain
optional model-dependent validation, with assumptions and predictive checks
explicit; they are not required for either primary test.

These are planning targets, not implemented diagnostics or a correction choice.
Keep baseline production authoritative, `no_empirical_residual` and zero legacy
residual scales. Random/dummy/slow-proton/baseline-pion subtraction, SIMC,
yield/uncertainty formulas, cuts, templates, priors, binning, efficiencies,
acceptance, L/T and cross sections remain frozen. Method A is detached/non-production;
Method B diagnostic/cross-check only and numerically excluded. Canonical-five
expansion is `DEFERRED`; final E.8/F.6.4 and absolute-SIMC interpretation remain
`BLOCKED` at their separate scopes. No correction redesign, response-matrix
implementation, production promotion or farm command is authorized here.
