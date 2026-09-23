# E.8 no-empirical-residual current-chain correction

## Context

On 2026-09-23, the E.8/F.6 forward-roadmap memory was corrected before E.8.2
implementation. Earlier roadmap text incorrectly treated legacy empirical
residual Fit 1/Fit 2 machinery as part of the accepted current analysis chain.
This record corrects that memory/decision boundary only.

## Source observation at correction

The observed repository HEAD was
`942444fb9586dc53e847039c97c6c3b52fe9a351`. At that source state,
`src/utility/background_config.py` specifies:

- `BG_STAT_SCALE1 = 0.0`;
- `BG_STAT_SCALE2 = 0.0`;
- `BG_OPT_ACTIVE_PROFILE = "no_empirical_residual"`;
- the active `no_empirical_residual` profile forces both scales to `0.0`.

The retained Fit-1/Fit-2 blocks in `src/binning/calculate_yield.py` are
conditional on their corresponding resolved background scale being greater
than zero. They are therefore dormant historical source machinery under the
accepted active profile, not an active scientific procedure stage.

## Corrected current chain

The user-confirmed accepted chain is:

    random/dummy
        -> slow proton
        -> baseline pion background
        -> Method-A pion-background reweighting
        -> cleanest kaon MM(t,phi)
        -> yields

For the baseline/Method-A comparison, both branches use the same proton-cleaned
input:

    K_proton - B_pi^0(w0)   = K_0 -> MM_0(t,phi) -> Y_0(t,phi)
    K_proton - B_pi^A(w0*C) = K_A -> MM_A(t,phi) -> Y_A(t,phi)

Legacy empirical residual Fit 1/Fit 2 predate the current pion-model method
and are excluded from E.8 presentation, the F.6.3 Method-A branch, and E.8.4
comparison. They must not be activated, rerun, tuned, or presented without a
separate future scientific contract.

## Preserved Method-A and evidence boundaries

Method A is reweighting, not rebinning: `b_j^0 = s_j * w0_j` becomes
`b_j^A = s_j * w0_j * C_j` on fixed canonical `(t,phi)` bins with the accepted
parent-t normalization contract. It does not introduce an empirical residual
fit or independently normalize a child. Method B remains diagnostic/cross-check
only and numerically excluded.

This correction changes no source, profile, artifact, or runtime behavior. It
does not claim new farm, ROOT/PyROOT, or runtime validation; existing accepted
F.1-F.6.2 and narrow E.8.1 Fix.5/Fix.6 runtime evidence retain their recorded
scopes.
