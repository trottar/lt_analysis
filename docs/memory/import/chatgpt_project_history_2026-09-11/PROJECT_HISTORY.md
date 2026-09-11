# KaonLT project development history

This is the chronological deep-memory record. Current-state rules belong in canonical active memory only after migration/reconciliation.

## 1. Long-running KaonLT analysis foundation

`[SCIENTIFIC_REFERENCE] [CHAT_HISTORY]`

KaonLT is a mature Jefferson Lab Hall C analysis of E12-09-011. The broader analysis includes HMS/SHMS coincidence reconstruction, PID/cointime cuts, efficiencies, kinematic offsets, random subtraction, fixed acceptance/binning, SIMC validation, yield/average-kinematics extraction, unseparated cross sections, L/T separation, and systematic uncertainties.

The user’s dissertation is the primary historical scientific reference for KaonLT apparatus, analysis architecture, kinematics, cross-section extraction, and the SHMS HGCer inefficiency context. Ali Usman’s dissertation is a complementary reference for the same experimental data set, Hall C/SIMC behavior, pion PID, SHMS response, SIMC resolution, and pion-production backgrounds.

These references supply physics context. Current repository source remains authoritative for 2026 implementation behavior.

## 2. Early 2026 pion-background development

`[GIT_HISTORY] [CHAT_HISTORY]`

By June 2026, Git history contains repeated empirical pion-background work on `pi-n` and `pi-Delta` subtraction for `Q2=3.0` at `W=2.32` and `W=3.14`. Representative chronology anchors include:

```text
395451f351d43b2610ab212ad20f3186428a6320
6cfdd47c39c2efe47fb39056c5cada73d74b1c4e
2a19bde89ace00ead2eafa5d899d79c2dd59fc0b
9db45f2ad72999d223e0790b2ce9a6f18ac1c9a2
64ea04e7f004e12cea1a68679b4396fe72f3e714
```

These commit subjects are chronology only; they do not independently prove final production behavior.

The model evolved toward staged pion components:

```text
pi_n -> pi_delta -> pi_sidis
```

with later joint refinement while preserving the scientific identity of each component.

## 3. SIMC background implementation

`[GIT_HISTORY] [CHAT_HISTORY]`

June 2026 source history then added/iterated SIMC background machinery. Useful anchors include:

```text
7a7584910b3a3932e8f2e8e5cfabef8ed55223c5
  Added Ali's background simc input files

e1640390bf67afad6b9bb5360cabe607d3c04fdc
ab43cac0ace94f63405fe50a5ec5ab5340ae1a4b
8284f087e0420fcb79d6b46d972d76e19c936b5d
277f269f64692a5a1972775c8526403ed7d659db
```

The mature production ownership became:

```text
src/cuts/rand_sub.py
    setting-wide ordering/runtime integration

src/cuts/proton_contamination_weights.py
    slow-proton treatment

src/cuts/pion_component_fits.py
src/cuts/pion_component_subtraction.py
    pion component model and baseline pion weights

src/binning/calculate_yield.py
src/binning/ave_per_bin.py
    downstream bin propagation

src/utility/background_config.py
    defaults and setting/phi overrides
```

## 4. Dynamic pion/SIMC alignment

`[CHAT_HISTORY] [HANDOFF]`

The pion model gained dynamic component-alignment scans. The core design separated the candidate fit region from a fixed evaluation envelope so candidates could not improve their score merely by moving away from difficult bins.

Durable rules included:

- rank all candidates on common evaluation bins;
- always evaluate an explicit baseline candidate;
- retain offset/window-expansion/support/lost-integral/localization diagnostics;
- reject unsupported or boundary/pathological scans rather than forcing a result;
- preserve staged component order instead of a full uncontrolled Cartesian scan;
- reuse the existing shifted-template implementation.

Preferred staged policy:

```text
pi_n:
    global + fine scans

pi_delta:
    global + fine scans when support permits

pi_sidis:
    inherit parent/current alignment by default
    scan disabled initially
```

This alignment machinery belongs to the pion-background model and is separate from the later HGCer Method-A/Method-B diagnostic program.

## 5. Slow-proton contamination becomes a separate correction

`[CHAT_HISTORY] [SCIENTIFIC_REFERENCE]`

A broad low-missing-mass contamination was identified in kaon-selected data, especially in some low-epsilon/left settings.

The central architecture decision was:

```text
slow proton:
PID contamination within kaon-selected data

pion background:
separate pion-production events mapped from pion-control data
```

They cannot be collapsed into one empirical MM background.

Production ordering became:

```text
random/dummy correction
-> freeze bin-counting spectra
-> slow-proton cleaning
-> pion component subtraction
-> optional residual diagnostics
```

## 6. Rejected global slow-proton proof-of-concept

`[CHAT_HISTORY] [SUPERSEDED]`

An early broad/global timing-aerogel proof-of-concept was useful diagnostically but not trusted as a production model. It left structured/non-Gaussian residual behavior and risked excessive subtraction in the protected `K Lambda` region.

Durable lesson:

> Agreement in the broad low-MM region is insufficient if genuine kaon signal is removed.

The method therefore moved toward local momentum/timing event-level treatment.

## 7. Local timing/PID slow-proton architecture

`[CHAT_HISTORY] [HANDOFF]`

The slow-proton procedure evolved to event-level probabilities constrained by timing and PID observables whose separation varies across SHMS momentum/acceptance.

A low-aerogel timing-center fallback hierarchy was established:

```text
nominal offset fit
-> expanded offset fit
-> stable-center fallback
```

A failed optional offset fit must not automatically reject a statistically supported cell.

No hard aerogel `5 NPE` production cut was introduced merely because diagnostics use a guide line.

## 8. K-Lambda preservation gate

`[HANDOFF] [FARM_EVIDENCE]`

A setting-wide safety gate was designed for timing-based proton cleaning:

```text
build proposed proton probabilities
-> measure proposed removal in protected K-Lambda window
-> compare to threshold
-> PASS: commit proposed probabilities
-> FAIL: bypass cleaning for entire setting
```

Default maximum Lambda removal:

```text
10%
```

Per-`t` values are diagnostic only; no partial application.

On failure:

```text
applied proton probability = 0
cleaned factor = 1
```

while proposed-model quantities remain available for diagnostics.

This also forced a terminology rule: `applied` refers only to committed post-gate production quantities.

## 9. Background procedure consolidation

`[SCIENTIFIC_REFERENCE] [HANDOFF]`

By late July/early August 2026, the production note documented two distinct corrections:

1. event-level slow-proton cleaning;
2. MM-dependent pion-background subtraction.

The reproducibility record was expected to retain configuration, sample provenance, dynamic pion alignment, serialized proton result, accepted/bypassed state, and diagnostic pages.

## 10. HGCer local pion-refinement question

`[CHAT_HISTORY]`

Once the baseline proton/pion treatment was working, the new question became whether the established pion event weight:

```text
w_pi,e^0
```

had an additional unmodeled local dependence across SHMS/HGCer acceptance.

The work was explicitly diagnostic before correction.

Two independent observables were defined:

```text
Method A:
detector-response driven
low-HGCer behavior in pion-control data

Method B:
outcome/closure driven
local MM closure of frozen pion prediction
against actual kaon host
```

Their independence became a core invariant.

## 11. Phase A — exact event-level baseline contract

`[HANDOFF]`

Phase A exposed the existing pion prediction and actual host state without changing production.

Later Phase C was explicitly forbidden to:

```text
reopen ROOT trees
rebuild event selection
rebuild pion weights
rerun proton cleaning
```

It had to consume frozen Phase-A records.

## 12. Phase A.1/A.2/A.3 — provenance and host hardening

`[HANDOFF]`

Recorded predecessor sequence:

```text
Phase A.1 — noRF / frozen-weight / provenance repairs
Phase A.2 — identity-host runtime wiring
Phase A.3 — independent upstream noRF identity-host closure
```

These repairs made host provenance explicit before Method B.

## 13. Phase B/B.1 — Method A

`[HANDOFF]`

Phase B established Method A; B.1 hardened statistics/diagnostics.

Approved source anchor before Phase C:

```text
18b06ec1aa8ba42859dd0980705ce8374d4720ab
Harden Method A statistical diagnostics
```

Later frozen Method-A response definition:

```text
positive: NPE > 0
low:      0 < NPE <= 2
control:  NPE > 2
```

with:

```text
f_low = N_low / N_positive
R_low_control = N_low / N_control
```

plus support and interval information.

Because `NPE=0` is excluded, current Method A is a positive-response relative leakage diagnostic, not yet an absolute leakage probability.

## 14. Phase C — independent Method-B local closure

`[HANDOFF]`

Method B asks whether the frozen `t`-binned pion prediction reproduces pion-sensitive MM regions equally well across `delta` within one canonical `t` parent.

Regional diagnostic:

```text
Q_ijr = H_ijr / P_ijr
```

same-`t` parent normalization:

```text
Qtilde_ijr = Q_ijr / <Q_ir>
```

It is not an absolute pion normalization.

Protected `K Lambda` and `K Sigma0` signals are excluded from scale determination.

The neutron/`pi_n` region may be a closure region but is not an absolute normalization anchor.

Method B consumes the total frozen pion prediction, not separate Method-B factors for `pi_n`, `pi_delta`, and `pi_sidis`.

No Method-A numerical values are allowed in Method B.

## 15. C.Fix.1 persistence failure

`[HANDOFF] [FARM_EVIDENCE]`

A key review failure occurred:

- live Method B produced `mm_regions` / `protected_regions`;
- the serializer did not persist them;
- the renderer was checkpoint-first;
- farm pages therefore lost the expected information.

Permanent source-review rule:

```text
producer
-> serializer/checkpoint
-> checkpoint-first display payload
-> consumer
-> renderer
```

must be traced for every persisted diagnostic.

## 16. C.Fix.2 adaptive partition diagnostic

`[HANDOFF]`

Adaptive Method B tested whether baseline-support-driven MM partitioning could improve local support without using Method A.

Partition identifier:

```text
baseline_delta_support_first_passing_protected_outward/v2
```

Candidate policy:

```text
0 slices -> unavailable
1 slice  -> Qtilde, available_single_slice
>=2      -> inverse-variance log-space, available_multi_slice
```

Multi-slice consistency was descriptive, not a veto.

C.Fix.2.2 accepted source:

```text
c7af5ada2ac6cb342c36cf351b99166ac6b1fbf1
```

and the next gate was deliberately one `Q4p4W2p74 Left-low` farm validation.

## 17. Phase-C five-setting closure

`[HANDOFF] [FARM_EVIDENCE]`

Accepted C.Fix.2.3:

```text
9a66bc62d20a99172e326e915866877b65ae1e5d
```

later accepted pre-E.3 source:

```text
e3853655db0809923cbf2326e2f779219128eda9
```

Five settings:

```text
Left lowe
Left highe
Center lowe
Center highe
Right highe
```

No Right-low.

Recorded result:

```text
all 30 C.Fix.2 pages rendered correctly
Phase C: PASS
adaptive Method B: DO NOT PROMOTE
```

Adaptive support was sparse:

```text
36/150 candidate-bearing = 24%
21/150 multi-slice       = 14%
114/150 unavailable      = 76%
```

Left-high had only 2/30 candidates, both single-slice.

Where measurable, consistency was not the problem:

```text
max chi2/ndf = 1.5573
max |log pull| = 2.486
largest approximate low/high discrepancy ~1.29 sigma
shared-cell cross-setting differences <~1.5 sigma
```

Same-`t` normalization was useful:

```text
raw Q median ~3.687, max ~277.83
Qtilde median ~1.013, max ~4.012
```

Conclusion:

```text
do not alter legacy Method B
do not alter production pion subtraction
no C.Fix.2.4
no D.Fix.1
adaptive remains diagnostic context
Phase D uses legacy Method B
Phase C scientifically closed
```

## 18. Phase D — independent A/B comparison

`[HANDOFF]`

Phase D compares the already-independent diagnostics.

Frozen states include:

```text
Both methods available
Both present; ratio undefined
Method A only
Method B only
Neither method available
```

For positive/comparable values:

```text
r_BA = B/A
Delta_AB = ln(B)-ln(A)
```

A real Method-A zero is retained and not floored to manufacture a ratio.

Phase D uses legacy Method B, not adaptive B.

Later work treats D as frozen upstream comparison/closure infrastructure. Exact final farm provenance is a known gap unless a dedicated bundle/handoff is found.

## 19. Phase E — presentation over frozen results

`[HANDOFF]`

Phase E is presentation-only. It consumes frozen upstream records/results and cannot silently recompute or promote physics.

## 20. E.3 — independent local Method-A presentation

`[HANDOFF]`

E.3 displays canonical `(t,delta)` HGCer context and stored Method-A support/`f_low`.

The signed noRF/no-MM-cut spectra shown for context are not the same population that owns the stored prompt-pion `f_low`. E.3 must not recompute `f_low` from displayed spectra.

No Method-B, adaptive-B, Phase-D numerical input, or correction is allowed.

## 21. E.3.Fix.2 farm-driven presentation repair

`[HANDOFF] [FARM_EVIDENCE]`

Initial detached tests missed:

1. header/title overlap with the first HGCer row;
2. E.3 pages interleaved rather than final.

Fix.2 accepted implementation:

```text
eb1710f4739ba6ef14f51419806e9fc5bd53c175
```

changed only:

```text
src/cuts/full_background_subtraction_plots.py
testing/test_full_background_subtraction_plots.py
```

and restored final page ordering:

```text
E.3 t1
E.3 t2
E.3 t3
```

## 22. JSON-driven validation-bundle architecture

`[HANDOFF] [FARM_EVIDENCE]`

The validation process itself was improved:

```text
generic stable collector
+
gate-specific JSON profile
```

Rules:

- collector packages the exact review artifacts automatically;
- profiles change for meaningful new gates;
- `complete=true` is not scientific PASS;
- stale checkout/source provenance must be machine-detectable;
- collector can allow validation-only commits above a required analysis commit while rejecting unexpected analysis-source changes.

Later validation-infrastructure HEAD:

```text
bf53dac84e1396cfd7e3f4e0234749426bfbdcf4
```

with no new analysis-source changes relative to E.3.Fix.2.

## 23. E.3.Fix.2 Left-low farm PASS

`[FARM_EVIDENCE]`

Runtime/bundle HEAD:

```text
bf53dac84e1396cfd7e3f4e0234749426bfbdcf4
```

required implementation:

```text
eb1710f4739ba6ef14f51419806e9fc5bd53c175
```

in ancestry.

Recorded checks:

```text
required analysis ancestry: PASS
py_compile: PASS
full-background tests: 61 OK
HGCer refinement tests: 63 OK
collector tests: 16 OK
git diff --check: PASS
```

The source procedure PDF had 61 pages; collector extracted 59–61 as E.3 `t1/t2/t3`.

Visual review passed complete 5x2 grids, no clipping/header overlap, visible `NPE=2`, low/control labeling, negative signed content, explicit unavailable cells, and no Method-B/correction language.

This specific Left-low gate is runtime validated. The recovered evidence does not by itself prove the remaining four E.3 settings passed.

## 24. Method-B role is later narrowed

`[HANDOFF] [SUPERSEDED]`

An earlier procedure note envisioned both A and B contributing to a future correction.

Later explicit scientific decision superseded that:

```text
Method B:
keep fully in analysis and keep plots
diagnostic/cross-check/historical comparison only
must not numerically adjust pion event weights

Method A:
sole candidate numerical HGCer leakage input
```

History preserves both states; the later decision is authoritative.

## 25. Phase-F roadmap

`[HANDOFF]`

Approved sequence:

```text
F.1 — detached Method-A acceptance event contract
F.2 — freeze map representation after five-setting support review
F.3 — detached Method-A map with support/OOD handling
F.4 — detached parent-preserving A-only correction
F.5 — detached event-level (t,phi) propagation
F.6 — production promotion only after F.5 validation
```

Training/response coordinates and yield coordinates are intentionally different.

Response may be learned in:

```text
canonical t parent + acceptance variables
```

but application happens event-by-event before filling downstream `(t,phi)` children.

Do not independently renormalize children; normalization belongs to the canonical `t` parent.

## 26. F.1 lineage and early repairs

`[HANDOFF]`

```text
b0f6868c17b0aae1d97293d4a103a3736fd50fd0  pre-F.1 base
b49dd1cdaf579407c5375d116e048dbf7681a4f7  F.1
4fa198e4e43b776f18ef2db7fc918df18857da27  F.1.fix.1
d656e15761970d7d612bb028d2746d077795e9ad  F.1.fix.2
81b25890ff5e600617db4040cd99e3c875cb9f7c  validation bundle infra
8b1ad5b735e1b5cc95a5d876526d321296cec16b  F.1.fix.3
eb253046ff6e23ec94c8638315c4bc712fa4f992  F.1.fix.4
```

Fix.1 repaired fingerprint/provenance and true `(t,phi)` page-5 mapping.

Fix.2 replaced unreadable page-5 text with a descriptive matrix (occupancy, delta, xptar, yptar) without adding a fit/correction.

## 27. F.1.Fix.3 farm failure

`[HANDOFF] [FARM_EVIDENCE]`

Farm exposed:

```text
ValueError:
truth value of an array with more than one element is ambiguous
```

from array-like edge geometry used in truth-value fallback expressions such as `value or ()`.

Fix.3 added safe edge materialization. Later farm output proved that mechanical regression closed.

## 28. F.1.Fix.4 farm failure

`[HANDOFF] [FARM_EVIDENCE]`

Next farm run exposed:

```text
parent_child_parity_mismatch:allcuts
```

because parent used Python `bool` while child cache could use a NumPy boolean scalar with identity-style comparison.

Fix.4 detached NumPy-like child scalars to native Python scalars at the F.1 read boundary.

Later five-setting run closed this mechanical regression.

## 29. Post-Fix.4 five-setting mechanical gate

`[FARM_EVIDENCE]`

All five settings produced mechanically valid F.1 artifacts:

```text
status=available
diagnostic_stage=complete
records nonempty
fingerprints populated
five F.1 pages rendered
```

Focused checks:

```text
full-background tests: 90 PASS
F.1 contract tests: 10 PASS
Phase-F runtime tests: 3 PASS
collector tests: 14 PASS
py_compile: PASS
git diff --check: PASS
```

This validated mechanics but exposed a scientific population-ownership problem.

## 30. F.1 v1 population-ownership bug

`[FARM_EVIDENCE] [HANDOFF]`

Method-A counts:

```text
Setting       positive  low  control
Left low      52397     576  51821
Left high     14862     285  14577
Center low    37035     568  36467
Center high   22429     479  21950
Right high    23963     736  23227
```

F.1 v1 showed:

```text
low = 0
control = exactly Method-A control
```

for every setting.

Root cause:

```text
F.1 v1 consumed Phase-A pion_records
from the authoritative downstream physical pion-control cache
which already applies NPE > 2
```

The downstream application cache was correct. The scientific mistake was using it as the Method-A response-training population.

## 31. F.1.Fix.5 dual-population architecture

`[HANDOFF] [CURRENT_SOURCE]`

Training source:

```text
pion_hgcer_tdelta_diagnostic["records"]["pion"]
+
pion_hgcer_tdelta_diagnostic["phase_e_acceptance_records"]["pion"]
```

exact identity join:

```text
(source_label, entry_index)
```

training selection:

```text
source_label == "prompt"
nommcuts == True
NPE > 0
```

classes:

```text
low:     0 < NPE <= 2
control: NPE > 2
```

Training owns F.1 pages 1–4 and future Method-A response modeling.

Application source remains Phase-A pion records + authoritative parent/child caches, physically restricted to `NPE>2`, and owns page 5 plus future event-level application.

Current pushed `test` HEAD observed during this import:

```text
7cdb7847d72501cd3dc504565fafc604e22a6132
Pion subtraction hgcer checker| Phase F.1.fix.5 || Added working memory
```

This proves Fix.5-era source exists on `test`; it does not prove farm validation.

## 32. Repository-owned active memory

`[CURRENT_SOURCE]`

The repo now contains `AGENTS.md` and `docs/memory/` active state.

Codex is required to read current memory, inspect relevant deep records, compare to live `test` HEAD, and update memory as part of completing meaningful work.

The bootstrap intentionally avoided inventing historical farm status. This import is the next step.

## 33. Current checkpoint

`[CURRENT_SOURCE]`

The next task is memory migration/audit.

No production source should be changed during migration.

After migration, a fresh Codex chat should be able to recover current state and deep history without a giant continuation prompt.
