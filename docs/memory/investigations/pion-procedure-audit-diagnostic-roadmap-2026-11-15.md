# November 15 pion-procedure audit and bounded diagnostic sequence

## Authority and evidence status

User authority: the Scientific Memory Update and D1 Preparation handoff following
the November 15 full pion-background procedure audit. This is a memory-only
checkpoint. CURRENT is the sole active-priority/exact-NEXT owner. D1 has not
begun. Independent memory-diff review, user publication and independent
pushed-state verification precede the next scientific operation.

Public/local test observed at this checkpoint:
`d811c7e6f9c2c4e113238fc578f136dacb9aae32`. The November 12 memory decision is
published, including the adopted exclusion and required earlier task contract.
The committed manifest has one verified integrity defect: it records that
contract as 8662 bytes/SHA-256
`a95917beb1e79e513ff5efee877fe52dfe5d13dd65da0e38169ee2837e06b7ad`, while
the committed and current local contract is 8661 bytes/SHA-256
`a6bcf9b79ac802a0d883ff5e199cd7978092adeb8b702ddce4c7f38886ebc19b`.
The mismatch is the previously identified final CRLF/LF byte difference.
This memory operation regenerates the manifest from actual bytes; it makes
no change to the contract. A passing local manifest is not a claim that this
correction is already published.

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
empty bin, invent an independent covariance or impose a positive SIDIS floor.
Exact operands: `numerical/Center_highe_t3_sidis_stage_bins.csv` and
`numerical/empty_bin_support.json` in the audit.

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

## Bounded D1-D5 sequence

This sequence describes diagnostic scope and order, not competing active-state
authority. CURRENT alone owns exact NEXT and its review/publication conditions.

1. **D1 — SIDIS-boundary stability and source-aware covariance.** Initially
   read-only, Center/high t3, using existing retained eight-bin and source-event
   evidence. Reconstruct signed prompt/random/dummy/dummy-random contributions,
   original observed/model populations, neutron-subtracted residuals, per-bin
   weighted scores and conditional variances. Trace fitted neutron-amplitude,
   data/template and training-covariance dependencies. Where the existing
   evidence supports a statistically justified uncertainty treatment, assess
   boundary sensitivity; otherwise identify precisely what is unidentifiable.
   Preserve empty observations. Do not remove bins, impose a positive SIDIS
   floor, invent missing covariance, replace the production estimator, rerun
   completed fit/bootstrap/model studies or start farm acquisition. Determine
   stability or inability to identify the coefficient; do not presume either.
2. **D2 — Accepted channel support.** Reuse completed raw canonical-t comparisons
   for neutron/Delta/SIDIS accepted populations, effective statistics, momentum,
   SHMS delta/phi coverage, shapes and channel overlap; no repeated study.
3. **D3 — Physical pion PID-transfer support.** Examine existing measurements/
   diagnostics of HGCer response, phase-space overlap, momentum/delta, detector
   position and available run dependence. State supported transfer constraints
   and unavailable terms; invent no calibration or covariance.
4. **D4 — Protected-Delta interpretation.** Classify recorded rejection reasons
   across 15 parents, distinguishing fit quality, component support, degeneracy
   and configured fallback. The existing rejection count is not a new fit.
5. **D5 — Predictive background closure.** Use independent or suitably held-out
   evidence only where available and supported. Same-data fit agreement is not
   independent closure; preserve formal Outcome B limitations.

No D1 execution, new estimator, production correction or source-changing task
contract is created by this memory update. Any justified correction requires
a separate standalone repository contract, independent review and explicit
user authorization. Stop this operation at the memory-review bundle.

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

Memory-only candidate; independent actual-diff review and user publication are
pending. No source/test/contract/fit/cut/normalization/weight/yield changes and
no D1/farm execution. CURRENT owns the conditional D1 NEXT. No MEMORY.md growth
is needed: lasting lessons and diagnostic references are owned by this record.
Existing CURRENT/MEMORY soft-size warnings are nonblocking; consolidate at the
next checkpoint/milestone audit, not through a new maintenance gate now.
