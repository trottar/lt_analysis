# F.6.1 accepted runtime closure

F.6.1 and F.6.1.Validation.1 are `CLOSED / RUNTIME VALIDATED` from the
accepted detached farm bundle. This evidence closes validation only; it does
not promote Method A, alter production pion subtraction, construct yields, or
replace the accepted baseline.

## Accepted bundle and source identity

- ZIP: `KaonLT_PhaseF6_1_validation_Q4p4W2p74.zip`.
- ZIP SHA-256:
  `22fb594f79a61e5378e76570df46cab9cf5195b6dc0faa04b0238919f66738a3`.
- Farm/bundle HEAD:
  `ffb7d4dc251610f6034eac11038a9841e3ef9f58`.
- Reviewed F.6.1 scientific source pin:
  `bfc4fe421f9fc9139a1992a0ec92e31aa101b86c`.
- Manifest: `complete = true`, `errors = []`,
  `required_analysis_commit_is_ancestor = true`, and
  `unexpected_committed_files_after_required_analysis_commit = []`.

## Accepted artifacts and exact continuity

- F.6.1 JSON SHA-256:
  `62bad2dae0f65f1fff87eeffd071c29dbbd4cdc15a43f37d6d9853cb58b25bd6`.
- F.6.1 PDF SHA-256:
  `2de4e98c0b22962cd9dbb68e385be1735ce5977b41378d13846c257cd9b874af`.
- F.6.1 artifact fingerprint:
  `377872a218a780481347402e4410c81bd568a4fb7682cd49f0f3352b499bad41`.
- F.6.1 validation fingerprint:
  `d992d789b3434897d76691df27c190a0f51479f67c0d5b496e84835e517c77f8`.
- Reproduced upstream JSON SHA-256 values: F.3
  `04a9a576767ba3f81c8c03daab2461998894c655f0d875bd532ad30cb46c0d95`,
  F.4 `adcc01900b8adc211e296aaedaa314fdb86207d51483ebfd4b3edf69f558f188`,
  and F.5.2
  `143e3af6b1c69560e5bf351155b570d05f19e7e0ceb14665f5448c351ad501be`.
- All five accepted F.1 source hashes match their frozen accepted values.
- Accepted authority matches; F.4 and F.5 reproduce exactly. The review has
  five settings, 15 parents, 135 canonical cells, and maximum absolute parent
  closure residual approximately `2.22e-16`.
- The accepted F.5 scientific fingerprint remains
  `d11b728d1089301a12c29e7f8b1798c6e0b6021ac47bd5d2afc2b62b47a1effa`.

## Detached validation boundary

The repaired source contract closes training/Part-1 `evt.ph_q` in radians to
physical application `phi_degrees = evt.ph_q * 180/pi`. The artifact is
aggregate-only, non-authoritative, validation-only, and requires manual
review. Production application/object mutation, yield and cross-section
construction, ROOT objects, Method-B numerical dependency, child
renormalization, smoothing/interpolation, persisted event correction, and a
shape-improvement gate are all false. The 27-page review PDF rendered
successfully and was manually inspected.

The captured farm source-state record included unrelated modified/untracked
model files. This is not a closure blocker: the committed-source gate found no
unexpected files after the reviewed source pin, and the F.6.1 inputs exactly
match accepted frozen artifacts.

## Scientific outcome

F.6.1 is a runtime-validation PASS, not evidence that Method A uniformly
improves every shape. Across nine one-dimensional comparisons in each of 15
parents, 68 of 135 moved toward the low-response reference and 67 moved away.

- Construction-space variables (`SHMS_delta`, `P_hgcer_xAtCer`,
  `P_hgcer_yAtCer`): 19/45 improve, 26/45 worsen, mean
  `H_MethodA - H_baseline` approximately `+0.0340`.
- Independent/secondary variables (`SHMS_xptar`, `SHMS_yptar`, `phi`): 23/45
  improve, 22/45 worsen, mean approximately `-0.00311`.
- Missing-mass/kinematic variables (`analysis_MM`, `Q2`, `W`): 26/45 improve,
  19/45 worsen, mean approximately `+0.0105`.

These counts are descriptive only, not an acceptance score or threshold. The
accepted conclusion is that Method A is technically closed, source-traceable,
authority-consistent, and exactly preserves accepted F.4/F.5 closure, but its
shape effect is parent-dependent. Agreement with the baseline pion model is
not itself the acceptance criterion, and disagreement with it is not an
automatic failure. F.6.2 is the next detached acceptance-correlated refinement
validation; see `../decisions/f6-2-acceptance-refinement-measurement-contract.md`.
