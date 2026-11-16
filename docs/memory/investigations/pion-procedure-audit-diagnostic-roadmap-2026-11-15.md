# November 15 pion-procedure audit and bounded diagnostic sequence

## Authority and evidence status

This record began with the November 15 Scientific Memory Update and D1
Preparation handoff. The November 16 D1/D1.a Memory Closure and D2 Preparation
user handoff now records D1 and D1.a as SOURCE REVIEWED within their bounded
scopes and authorizes read-only D2 after memory review/publication/synchronization.
CURRENT remains the sole active-priority/exact-NEXT owner. This operation updates
memory only: no D2, source implementation, production reanalysis or farm run.

November 15 historical public/local test observation:
`d811c7e6f9c2c4e113238fc578f136dacb9aae32`. Its one contract-manifest CRLF/LF
discrepancy recorded 8662 bytes/SHA-256
`a95917beb1e79e513ff5efee877fe52dfe5d13dd65da0e38169ee2837e06b7ad`
instead of the actual 8661 bytes/SHA-256
`a6bcf9b79ac802a0d883ff5e199cd7978092adeb8b702ddce4c7f38886ebc19b`.
That historical integrity discrepancy was repaired in the published D1-preparation memory at
`c567382027cca3d6aa8c1019c298459d046a33ee`, the current live/local base observed
at this checkpoint. The manifest now passes. D1.a source/test changes remain
local/uncommitted; source review is not publication or runtime acceptance.

Audit artifact: `KaonLT_Full_Pion_Background_Procedure_Audit_20261115_review_bundle.zip`,
6030458 bytes, SHA-256
`f1e27b578d7d9cdfb28719430d9b5e112702767cbd45f46c8cca93fc7047951c`.
All 59 members/CRC and sealed member identities were checked for this checkpoint.
`REPORT.md`/`REPORT.html` contain six project-data plots; `source/` records exact
source identities/lines, and `numerical/` records the actual saved operands.
The audit executed 1120 deterministic replay assertions, not independent
physics experiments. No new fit, bootstrap/model campaign, ROOT/main.py run,
farm execution, production correction or runtime acceptance occurred. This
handoff records the user's specified source-bound findings; it does not infer
independent source-review acceptance from bundle existence or promote a status.

## Established scientific findings and limits

### Center/high t3 SIDIS boundary

The saved raw canonical-t candidate has A_n=0.182202244064,
A_Delta=0.569000638411, A_SIDIS=0 and B_SIDIS=2.984506898807.
The eight selected SIDIS-stage MM centers span 1.063 to 1.105 GeV in native
6-MeV bins; the retained anchor is
`[1.0601751288600272, 1.1101751288600272]` GeV. After the neutron stage,
the SIDIS weighted residual numerator is -38.7723381448 and denominator
250.768868521, giving unconstrained amplitude -0.154613841716. The solver's
nonnegative boundary correctly returns zero. The unweighted residual sum is
positive; its sign is not the sign of this inverse-variance estimator.

The native bin `[1.066,1.072)` GeV, centered at 1.069, has zero original kaon
rows in every prompt/random/dummy/dummy-random source, not signed cancellation.
Its negative neutron-subtracted residual contributes -49.381752867 to the
weighted numerator. Its retained residual sigma, 0.000231675886, is the
fixed-amplitude neutron-template MC error. The early fitted neutron amplitude
has conditional stage sigma 0.162155076037, but the complete learned-neutron/
data/template covariance is not included in that residual treatment.
This establishes the numerical origin, not stability under a justified complete
uncertainty treatment and not physical absence of SIDIS. Do not delete this
physical histogram bin, invent an independent covariance or impose a positive SIDIS floor.
Exact operands: `numerical/Center_highe_t3_sidis_stage_bins.csv` and
`numerical/empty_bin_support.json` in the audit.

November 16 D1 closure: the remaining seven bins contribute +10.6094147226.
Zero experimental entries and zero experimental Sumw2 in all four signed
streams establish genuine absence of observed support, not cancellation. The
small template-only residual error assigns excessive inverse-variance weight
in this particular staged WLS procedure. The earlier fixed-weight neutron
amplitude sign-crossing argument is not a valid physical stability test when
the implemented variance dependencies are restored. D1 is SOURCE REVIEWED,
not proof of physical SIDIS absence or complete covariance. The user accepts
this defect as established; do not restart a general D1 investigation.

The [D1.a implementation record](../phases/kaon-sidis-empty-bin-support-2026-11-16.md)
owns the exact reviewed opt-in support rule, contract, source identities,
eight/seven-bin comparison and test limits. Its seven-bin G=10.6094147226,
H=48.4742198501 and nonnegative A_SIDIS=0.218867157746 follow from unchanged
retained weights, without a positive floor. Physical histogram bins remain.
Only Center/high t3 changes in the retained fifteen-parent fixed-weight replay;
the 13 focused tests pass and the broader inherited error remains documented.
D1.a is SOURCE REVIEWED, unpublished/default-off and not farm validated.

### Accepted templates and production ownership

Ordinary production parents inherit setting-wide pion shapes: main's loader
does not pass t/phi edges. In the detached scoped path, adding independently
unit-normalized nonempty phi children produces equal-phi normal-bin mixtures,
not event-population-weighted canonical-t parents. Both differ from existing
raw event-weighted canonical-t templates. For Center/high t3 the normal-bin
L1(raw versus equal-phi) distances are 0.37110741 neutron, 0.17867978 Delta,
0.16469883 SIDIS; these are shape distances, not significances or efficiencies.
All 45 raw parent/child closure and retained alignment identities were checked.
Baseline and raw-local candidates both retain A_SIDIS=0 for Center/high t3;
template replacement alone is not the demonstrated cure. The audit's proposed
template contract is unapproved and is not the current next operation.

### Protected Delta and transfer interpretation

The retained audit reports 11/15 quality-rejected protected-Delta candidate
fits and 13/15 baseline fits. The configured `zero_pi_delta` failure policy can
set the applied Delta amplitude to zero while the outer component result stays
accepted. This is policy execution, not proof of absent Delta. Center/high t3
itself has a successful declared Lambda-only protected fallback, distinct from
the generic outer fallback flag; early neutron/SIDIS coefficients stay fixed.
Outer acceptance does not validate each component physically.

The pion-control decomposition defines the background-channel model; its transfer
to kaon-selected data requires physical PID/acceptance support and properly
propagated uncertainty. The implemented sum(A_i T_i)/sum(B_i T_i) is an
arithmetically reproducible fitted-mixture ratio; independent A estimates are
not measured channel-specific detector transfer probabilities. Do not simply
copy B into A or force positivity. Pion-control NPE>2 and kaon-selected NPE==0
are different detector populations. Existing 135 signed transfer/Sumw2 replay
closures are not independent physical-purity or predictive-closure tests.

Preserve the accepted/reviewed source-aware 1024-replica bootstrap, 116 coherent
model responses and their scopes. Fixed protected-kaon numerator, learned
targets/alignment, unavailable full training/PID/model/cross-setting covariance
and empty-cell classifications remain explicit. Formal Outcome B is BLOCKED:
no qualifying independent joint non-SIDIS anchor/pion-gauge constraint has been
established. This is separate from the authorized read-only diagnostic sequence;
it does not commission a new evidence-acquisition search.

## Bounded diagnostic sequence and D2 scope

These are diagnostics, not new production milestones. CURRENT alone owns the
exact conditional NEXT; the six production-step labels/dependencies are frozen.

| Diagnostic | Recorded state and boundary |
| --- | --- |
| D1 | **SOURCE REVIEWED**: zero-SIDIS numerical origin and source-aware conditional covariance limitations established; no physical-absence claim or repeat general study. |
| D1.a | **SOURCE REVIEWED**: opt-in original-support weighted-fit correction; user publication and narrow runtime validation remain separate. |
| D2 | **NEXT after memory publication and independent synchronization**: accepted neutron/Delta/SIDIS component shapes and canonical-t phase-space support; read-only reuse of retained evidence. |
| D3 | **DEFERRED**: physical HGCer pion-control-to-kaon PID transfer. |
| D4 | **DEFERRED**: protected-Delta rejection and fallback interpretation. |
| D5 | **DEFERRED**: independent or suitably held-out pion-background predictive closure; same-data fit agreement is not independent closure. |

**D2 question:** Do the accepted neutron, Delta and SIDIS SIMC templates correctly
represent the canonical-t populations used to model pion backgrounds? Here
accepted denotes the selected/accepted simulation population, not demonstrated
physical adequacy of the background model.

Reuse already retained raw canonical-t projections and earlier comparisons;
do not repeat a completed study or launch another unrestricted fit campaign.
Inspect across all five settings/fifteen canonical-t parents:

1. Accepted channel populations versus t, phi, SHMS delta, momentum and missing mass.
2. Finite-SIMC support and effective statistics.
3. Production setting-wide shapes versus raw event-weighted canonical-t shapes.
4. Equal-phi aggregation effects, keeping detached diagnostic and production ownership distinct.
5. Neutron/Delta/SIDIS shape overlap and degeneracy.
6. Retained component-model and Lambda-background sensitivity, without new unrestricted fits.

The earlier Center/high t3 raw-versus-setting-wide normal-bin L1 differences
are neutron 0.639572 (about 0.640), Delta 0.301043 and SIDIS 0.399584 (about
0.400). These are shape differences, not significances or automatic corrections.
Distinguish actual construction defects from finite-MC limits and physical
acceptance differences. The earlier proposed canonical-t template contract is
not automatically implemented; D2 is diagnostic, not source-changing production.

D1.a publication and farm validation are required before claiming its accepted
production behavior, but do not block separately authorized read-only D2.
Formal Outcome B remains BLOCKED as the distinct independent-component-constraint
requirement. No D2, D1.a farm operation or source implementation occurs in this
memory checkpoint. Stop at its review bundle; publication/synchronization follow.

## Literature references for diagnostics, not production prescriptions

Bibliographic identities were checked against the primary records/papers below.
These are methodological references, not KaonLT calibration, numerical priors,
detector-response measurements or implementation approval. Detector-specific
numbers cannot be imported from another experiment.

- **Hall C pion analysis:** H. P. Blok et al., “Charged pion form factor between
  Q2=0.60 and 2.45 GeV2. I. Measurements of the cross section for the
  1H(e,e'pi+)n reaction,” Physical Review C 78, 045202 (2008),
  [arXiv:0809.3161](https://arxiv.org/abs/0809.3161).
  Relevance: experimental calibration, physical selection and systematic analysis.
- **Hall C HGCer:** George Ryan Ambrose, “Blinded by the Light: Commissioning
  of the Hall C SHMS Heavy Gas Cherenkov Detector,”
  [JLab publication record 15789](https://misportal.jlab.org/sti/publications/15789).
  Official indexed title/author verified; year/degree metadata are not asserted
  here because direct record/PDF retrieval was unavailable. The indexed primary
  thesis describes momentum-dependent pion efficiency and response geometry;
  it supplies no measured KaonLT transfer coefficient for this task.
- **LHCb PID calibration:** R. Aaij et al., “Selection and processing of
  calibration samples to measure the particle identification performance of
  the LHCb experiment in Run 2,” EPJ Techniques and Instrumentation 6, 1 (2019),
  [arXiv:1803.00824](https://arxiv.org/abs/1803.00824),
  [publisher record](https://link.springer.com/article/10.1140/epjti/s40485-019-0050-z).
  Relevance: true-particle calibration samples and kinematically dependent PID
  efficiency; not direct amplitude copying or a transferable detector number.
- **CMS control regions:** CMS Collaboration, “A search for new phenomena ...
  using the alpha_T variable” (title abbreviated), EPJ C 77, 294 (2017),
  [publisher article](https://link.springer.com/article/10.1140/epjc/s10052-017-4787-8).
  Relevance: control-to-signal transfer factors, uncertainty and closure tests
  with appropriate validation regions; not same-data fit closure by relabeling.
- **Finite MC:** R. Barlow and C. Beeston, “Fitting using finite Monte Carlo
  samples,” Computer Physics Communications 77, 219-228 (1993),
  [DOI:10.1016/0010-4655(93)90005-W](https://doi.org/10.1016/0010-4655(93)90005-W).
  Title/authors/journal/pages and DOI verified in indexed primary-paper copies;
  direct PDF retrieval failed. Finite-template sampling must be accounted for;
  an unweighted-Poisson template construction is not automatically a likelihood
  for KaonLT signed random/dummy-subtracted data.
- **Weighted statistics:** G. Bohm and G. Zech, “Statistics of weighted Poisson
  events and its applications,” NIM A 748, 1-6 (2014),
  [arXiv:1309.1287](https://arxiv.org/abs/1309.1287).
  Relevance: weighted-event statistical structure. Applicability to signed,
  learned-weight populations and their cross covariance must be demonstrated.
- **Boundaries/nuisances:** G. Cowan, K. Cranmer, E. Gross and O. Vitells,
  “Asymptotic formulae for likelihood-based tests of new physics,” EPJ C 71,
  1554 (2011), [arXiv:1007.1727](https://arxiv.org/abs/1007.1727).
  Relevance: physical boundaries and nuisance-parameter treatment. Asymptotic
  results do not guarantee coverage for eight sparse signed bins, and zero
  fitted coefficients do not prove absent physical components.

## Preserved progression and ownership

The six existing [roadmap milestones](../roadmap/STATUS.md) and their exact
work-state labels/dependencies remain unchanged. Intended progression:
physically defensible pion subtraction -> selected Local-t -> native (t,phi)
reconstruction -> adopted mask -> Lambda spectra/yields -> existing SIMC
cross-section iterations -> existing L/T separation. Diagnostics D1-D5 do not
promote production Step 1 or reopen accepted phases.

Preserve the [November 12 decision](../decisions/lambda-mask-and-pion-decomposition-priority-2026-11-12.md):
four criteria, >=30 total and >=20 signal entries plus both neighboring signed-
yield comparisons; reviewed 67/135 retained-cell cohort; matched experimental/
SIMC support; unchanged setting normalization. Mask adoption is separate from
integration/farm acceptance. Every data/SIMC parent MM comparison uses identical
native or exactly aggregated edges, with separate unit-normalized shape and
original yield-normalized views; no fabricated fractional-bin content.
Iteration-zero data/SIMC unity or phi-flatness is not an acceptance requirement.
Initial pion subtraction stays fixed during subsequent SIMC iteration.
Random/dummy/proton ownership stays frozen. Method A remains detached and
non-production; Method B remains diagnostic-only and numerically excluded.

## Review boundary and maintenance

This November 16 memory-only candidate records the user's D1/D1.a source-review
closure and conditional D2 decision. It preserves all pending implementation,
tests, coefficients, histograms and scientific/runtime boundaries. No numerical
diagnosis, D2 execution or farm/source operation was repeated. The next sequence
is independent memory actual-diff review -> user-only memory publication ->
independent pushed-state synchronization -> read-only D2. CURRENT alone owns NEXT.
Existing CURRENT/MEMORY soft-size warnings are nonblocking; batch consolidation
at the next checkpoint/milestone audit without growing MEMORY.md here.
