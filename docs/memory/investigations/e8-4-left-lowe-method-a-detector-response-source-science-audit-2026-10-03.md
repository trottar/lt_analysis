# Left/lowe Method-A detector-response source/science audit — 2026-10-03

## Authority, identity and scope

This record preserves the completed read-only Gate-3 audit supplied in the
[Gate-4 contract](../phases/e8-4-left-lowe-method-a-gate4-final-memory-reconciliation-task-contract.md).
Gate 4 made targeted source/evidence confirmations at branch `test`, local
HEAD/origin/test `3a3a4217c78dde07123e71cd36e43ee54666d1e4`; it did not repeat
the full audit, change scientific source or introduce farm evidence. The
scientific question is Q4p4W2p74 / Left / lowe, especially low-t/t1 redistribution.
CURRENT alone owns active state and exact NEXT. This is an investigation record,
not an implementation prescription or production-promotion decision.

## SOURCE VERIFIED — population and response ownership

`src/cuts/pion_hgcer_method_a_acceptance_contract.py` preserves:

| Role | Selection |
| --- | --- |
| Training | prompt / noRF / nommcuts / `P_hgcer_npeSum > 0` |
| Training low-response class | `0 < P_hgcer_npeSum <= 2` |
| Training control-response class | `P_hgcer_npeSum > 2` |
| Physical application | authoritative physical-pion-control population, `P_hgcer_npeSum > 2` |

Training and application have separate records and fingerprints; application
is not the training population. The target zero-photoelectron pion population
is absent from the positive-response control sample.

`src/cuts/pion_hgcer_method_a_acceptance_map.py` constructs a support-aware
relative detector-response representation on the accepted F.3 basis:

```text
hgcer3 = (SHMS_delta, P_hgcer_xAtCer, P_hgcer_yAtCer)
relative response = exp(coefficients dot robust-scaled hgcer3), without intercept
```

The intercept is provenance only. F.3 constructs neither an absolute pion
leakage probability nor an absolute zero-photoelectron probability, production
correction or production-side weight. Do not rename it a probability map.
The accepted basis/support boundary is not reopened by this audit.

## SOURCE VERIFIED — F.4 signed parent-preserving mathematics

`src/cuts/pion_hgcer_method_a_parent_preserving_correction.py` uses physical
application contributions within each canonical-t parent:

```text
b_j = signed_source_coefficient_j * w0_j
r_j = relative Method-A response
B   = sum_j b_j
U   = sum_j b_j * r_j
C_j = r_j / (U / B)
sum_j b_j * C_j = B
```

Source requires finite positive B, U and normalization, checks parent closure,
and applies neutral raw response outside the frozen F.3 support before parent
normalization. The equality is a signed-parent constraint, not a bound on
individual corrections. It does not preserve each `(t,phi)` child, MM subregion,
Lambda-window yield or individual event contribution. Independent child
renormalization is forbidden. Full parent application support and cut-window
template support are distinct integrals.

`src/cuts/pion_hgcer_event_contract.py` builds `w0` through
`simc_shape_pion_weight_from_value` in `src/cuts/pion_component_subtraction.py`
and records the signed coefficient separately. The baseline weight already
represents the established pion-control-to-kaon-background transfer from the
accepted pion component model. A future response/transfer quantity cannot
automatically supply another independent contamination normalization or be
multiplied into `w0`: distinct scientific ownership and absence of double
counting must first be established. This constraint approves no algorithm.

## SOURCE VERIFIED — geometry and existing diagnostic machinery

The relevant diagnostic/application selections already exclude the configured
HGCer geometric hole. `src/cuts/particle_subtraction.py` applies `not
hole_rejected` to allcuts and nommcuts; the acceptance/event contracts consume
the selected populations. Separate no-hole display populations do not change
this ownership. The present redistribution must not simply be called a
correction of events inside the known hole. Response dependence elsewhere in
the HGCer plane remains possible; no hardware cause follows from exclusion.

`src/cuts/pion_hgcer_transfer.py` already contains zero-photoelectron transfer
diagnostics using zero-truncated detector-response families to infer a
control-to-zero-photoelectron response transfer. The machinery predates the
later accepted `hgcer3` program. It is non-authoritative, diagnostic-only and
production-side-effect-free, with no production subtraction entry point. Its
existence does not validate its model for the present lineage or make it an
automatic Method-A replacement. Nothing is copied, promoted or applied here.

The slow-proton architecture in `src/cuts/proton_contamination_weights.py`
separates event probability, proposed effect, support/preservation diagnostics
and applied-effect decisions. It is a useful architectural analogue only.
The pion problem has a statistically absent target population under the
positive-response selection; no slow-proton numerical model can be transplanted
into pion subtraction on this analogy.

## RUNTIME VERIFIED — prior accepted reviews only

These are existing accepted runtime scopes recovered from canonical records,
not raw artifacts newly inspected or runs performed by Gates 2–4:

- [F.1](../evidence/f1-fix5-v4-bundle-inspection-2026-09-14.md),
  [F.2](../evidence/f2-fix1-runtime-closure.md) and
  [F.3](../evidence/f3-runtime-closure.md) retain accepted detached scopes,
  population distinctions and the frozen `hgcer3` basis.
- The [current-baseline comparator](../evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md)
  reported exact F.2/F.3 scientific payload matches and F.4 as the first changed
  scientific stage. Raw serialized provenance can differ without a scientific
  payload change. This does not replace historical accepted authorities.
- F.4 detached closure and F.4.Refresh.2 materialization remain accepted at
  their recorded scopes; [F.6.1](../evidence/f6-1-runtime-closure.md) and
  [F.6.2](../evidence/f6-2-scientific-runtime-closure.md) remain detached
  validation, with scientific and presentation provenance separate.
- [F.6.3/E.8.4](../evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md)
  retain only Q4p4W2p74 / Left / lowe branch execution, current-lineage
  application/live-cache parity, real child changes and signed parent preservation.
  Candidate F.4 SHA-256 is
  `1d545924eba89c7f9ffa28028e307aca9b434a89beec06863cf2893887b6b902`.
  The packaging-time model-file caveat does not prove pre-run cleanliness or a
  scientific failure. Historical E.8.3 F.6.1 is a separate lineage.
- [Fix.5.7](../evidence/e8-4-fix5-7-left-lowe-runtime-and-fix5-visual-gate-2026-10-03.md)
  retains its narrow owner/checker setting-provenance closure;
  [Fix.5.8](../evidence/e8-4-fix5-8-left-lowe-runtime-closure-2026-10-03.md)
  retains its narrow Left/lowe structural/visual legibility closure. Neither
  establishes a detector cause, absolute-SIMC interpretation or production validity.

The comparator therefore concentrates the scientific question at the interaction
of accepted relative HGCer response, current physical application population,
baseline `w0` transfer and signed F.4 parent normalization. It does not show
F.2/F.3 themselves changing scientifically.

## MEMORY ONLY — workflow and historical drift

The supplied Gate-4 contract records Gate 2's completed read-only memory audit:
no blocking ownership/evidence/boundary contradiction, with consumed live wording
and nonblocking historical drift. Its reported outside-allowlist drift includes
KNOWN_GAPS, the dated F.1 identity investigation, older F.6 promotion/phase
chronology, the E.8 procedure decision and older F.6.2 downstream wording.
Those historical records are preserved; they do not acquire CURRENT authority.
This section records the supplied audit rather than claiming a fresh complete
repository-memory audit. Canonical evidence above was read for scope continuity.

## INFERENCE — signed-normalization sensitivity hypothesis

Because B and U are sums of signed contributions, positive/negative cancellation
could make `U/B` more sensitive than a positive-measure response average.
Parent preservation alone cannot diagnose that sensitivity. This is a plausible
hypothesis, not an accepted cause or proof that signed cancellation explains the
observed t1 effect. Relative response is not an absolute leakage probability;
neither the hypothesis nor arithmetic integrity establishes production correctness.

## NOT VERIFIED — current-lineage t1 decomposition and detector cause

The unresolved question is why the current Left/lowe candidate produces large
low-t/t1 redistribution and downstream pion-background/yield impact, and whether
a physically appropriate detector-response model supports it rather than an
inappropriate normalization/transfer construction. No cause is accepted.

Exact current-lineage diagnostics still required, without implementation here:

```text
sum b_j
sum |b_j|
|sum b_j| / sum |b_j|
sum b_j r_j
source-separated signed support and source-sign/source-class decomposition
raw r_j distribution and tails
final C_j distribution and tails
(t,phi) redistribution
missing-mass-region redistribution
hgcer3 response-coordinate dependence
support/OOD dependence
comparison to physically interpretable zero-photoelectron response transfer
```

Existing source diagnostic fields do not mean this exact current-lineage
decomposition has been evaluated. No particular mirror, PMT, optical alignment,
track geometry or hardware origin is established. Bookkeeping consistency is
not the primary blocker; no stale-scalar/wrong-histogram defect is established
as the main issue. The general histogram/scalar identity rule remains intact.

Absolute-SIMC provenance/units remain a separate blocker:
`SIMC_normfac_luminosity_and_charge_units_not_source_proven`. Valid current-F.6.3
data/identity/yield/closure payloads remain available; absolute overlay amplitude
does not prove agreement, disagreement or incorrect normalization.

## Preserved boundaries and closing checkpoint

No new runtime evidence exists from Gates 2–4. Keep baseline production
authoritative, `no_empirical_residual`, zero legacy residual scales and all
accepted ownership. Random/dummy/slow-proton/baseline-pion subtraction, SIMC,
yield/uncertainty formulas, cuts, templates, priors, binning, efficiencies,
acceptance, L/T and cross sections are frozen. Method A remains detached/
non-production; Method B diagnostic/cross-check only and numerically excluded.
Canonical-five expansion remains `DEFERRED`; final E.8/F.6.4 remain `BLOCKED`.

After Gate-4 actual-diff review, user commit/push and pushed-state synchronization,
the project holds a scientific-direction checkpoint to decide whether an external
analogous-experiment / PID-background methods review should precede any
implementation contract. This record authorizes neither research nor implementation.
