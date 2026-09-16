# F.6.2 acceptance-refinement measurement contract

F.6.2 is `NEXT`, not active. This is a pre-implementation scientific contract
for a detached acceptance-correlated Method-A refinement validation. It does
not authorize production mutation, Method-B numerical input, a new pion model,
or a claim about final kaon-yield uncertainty.

## Physical interpretation

The accepted baseline pion model remains the guide for the broad, physically
realistic missing-mass distribution, including pion leakage beneath the kaon
signal. Method A is an acceptance-dependent refinement of that baseline, not
an independent replacement model. A large Method-A/baseline departure is a
question for physical explanation, not an automatic failure.

Pion leakage can be present in kaon-bearing regions. A valid refinement need
not avoid the kaon signal region; it must be assessed there by whether its
change is consistent with observed pion-leakage topology and the corresponding
acceptance information. The low-HGCer population is an observed leakage
reference shape, not an absolute leakage probability or total-normalization
ground truth.

## Unit, populations, and normalization

The primary diagnostic unit is each canonical `(setting, t, phi)` child. For
each shape variable, define `L` as the observed low-HGCer response reference,
`B` as the baseline physical pion-control prediction, and `A` as the Method-A
refined physical pion-control prediction. Normalize `L`, `B`, and `A`
independently to unit area only for shape comparisons. Retain actual signed
event contributions separately and never mix their meaning with normalized
shape quantities.

## Per-variable interpretable measurements

For each variable `v`, persist raw metrics rather than a composite quality
score:

```text
H_B(v) = H(L_v, B_v)
H_A(v) = H(L_v, A_v)
DeltaH(v) = H_A(v) - H_B(v)

TV(L,B) = 1/2 sum_i |L_i - B_i|
TV(L,A) = 1/2 sum_i |L_i - A_i|
R(v) = TV(A,B) = 1/2 sum_i |A_i - B_i|
```

`DeltaH < 0` moves toward the low-response reference; `DeltaH > 0` moves
away. `R` measures refinement size, not failure. Define baseline residual
`r_i = L_i - B_i` and correction `c_i = A_i - B_i`, then, when both norms are
nonzero:

```text
kappa(v) = sum_i r_i c_i / sqrt[(sum_i r_i^2)(sum_i c_i^2)]
rho(v) = sqrt(sum_i c_i^2) / sqrt(sum_i r_i^2)
```

`kappa` near `+1` is aligned with the observed residual, near zero is largely
unrelated, and negative is opposite. `rho` compares correction and residual
amplitudes; neither `rho ~ 1` nor a large value is an automatic requirement or
failure. Zero-norm and unavailable cases must be persisted explicitly, never
substituted silently.

Keep construction-space variables (`SHMS_delta`, `P_hgcer_xAtCer`,
`P_hgcer_yAtCer`) separate from primary independent acceptance corroboration
(`SHMS_xptar`, `SHMS_yptar`). Other frozen matched quantities may be used only
when consistently available and scientifically justified; none becomes a
Method-A input.

## Joint MM, kaon-window, and support diagnostics

For each child, construct transient matched `(MM, SHMS_xptar)` and
`(MM, SHMS_yptar)` distributions. For each selected acceptance quantity `a`,
retain raw 2D distributions and, where defined:

```text
r_ij = L_ij - B_ij
c_ij = A_ij - B_ij
kappa_MMxa = sum_ij r_ij c_ij / sqrt[(sum_ij r_ij^2)(sum_ij c_ij^2)]
```

This tests whether an MM leakage feature and an acceptance feature identify
the same localized population. Two-dimensional total variation may supplement,
but not replace, these distributions and alignment.

Use the existing frozen kaon-signal MM window without retuning it:

```text
P_B^K = baseline pion contribution in the window
P_A^K = Method-A pion contribution in the window
DeltaP^K = P_A^K - P_B^K
f_refine^K = (P_A^K - P_B^K) / P_B^K, when P_B^K is valid
```

Persist `N_low`, `N_control`, and support/OOD counts or fractions. For positive
physical-control weights only, retain
`N_eff = (sum_j w_j)^2 / sum_j w_j^2` for `w0` and `w0*C`. Do not apply
`N_eff` to signed background contributions. For those signed contributions,
retain diagnostic proxies `V_B = sum_j (s_j*w0_j)^2`,
`V_A = sum_j (s_j*w0_j*C_j)^2`, and `R_V = V_A/V_B` when defined. These are
not final yield uncertainties.

## Bootstrap and descriptive cases

Require event-level bootstrap intervals for at least `DeltaH`, `kappa`, `rho`,
and `DeltaP^K` or `f_refine^K`. Resample the shared B/A control population as
paired data; resample the low-response reference consistently but independently
unless source ownership later proves otherwise. Freeze replica count, seed, and
confidence-interval convention in the implementation contract before final
classification results are examined. Do not tune bootstrap or binning post hoc.

Use no composite score and no numerical case thresholds. Cases are descriptive:

- **A — baseline adequate / refinement small:** small discrepancy and `R`, with
  no statistically meaningful movement.
- **B — acceptance-supported refinement:** appreciable discrepancy and
  refinement, positive physical alignment, preferably movement toward `L`,
  independent acceptance corroboration, and MM×acceptance localization. It
  does not require Method A to stay close to baseline.
- **C — statistically inconclusive:** inadequate low/control support, broad
  intervals, unresolved movement, or insufficient support for a conclusion.
- **D — physically concerning / unexplained refinement:** substantive change
  with weak/opposite alignment, no independent corroboration, or no physical
  MM/acceptance explanation. Variance inflation can strengthen concern but is
  not alone sufficient.

The objective is an overall physically sound, acceptance-supported, explainable
pion-background refinement, especially in contamination-dominated regions,
while retaining the baseline as the broad-MM guide. Not every child or metric
must improve. A later F.6.3 full-procedure phase alone may test final extraction
and propagated uncertainty effects.
