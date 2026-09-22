# Farm validation-bundle wrapper Fix.1 state reconciliation

## Scope

This is a post-push active-state reconciliation only. At live `test` and
`origin/test` HEAD `31abff9b8ab02f715cf7d6285a8dd2f53855a5e4`, the user
reported that ChatGPT's source review of the pushed parameterized farm
validation-bundle wrapper passed.

## Reconciled state

The wrapper source at `31abff9b8ab02f715cf7d6285a8dd2f53855a5e4` is
`SOURCE REVIEWED`. That source review is not farm execution, ROOT/PyROOT
validation, bundle acceptance, or E.8.1 runtime acceptance. The existing
[wrapper phase record](farm-validation-bundle-wrapper.md) remains an immutable
historical account of the pre-push local/proposed state and is not revised.

`CURRENT.md` now removes the obsolete review/commit/push prerequisite and
restores the sole `NEXT` to the targeted E.8.1 farm gate: use the reviewed
wrapper to package and inspect `Q4p4W2p74 / Left / highe` first, then broaden
only after that detailed gate passes.

## Preserved identities and boundaries

This reconciliation preserves the E.8.1 ordinary procedure-PDF source
`0bf1281ebf44868b68ff808090436653eea6ed60`, canonical bundle/profile source
`c25f9d8248f2acfb8b0443a482fc76781eb405ca`, and frozen F.6.2 scientific
source/artifacts. It does not modify analysis source, the wrapper, its tests,
the collector, profiles, production physics, Method-A/Method-B ownership,
accepted yield, cuts, normalization, templates, or binning.

## Validation boundary

Only `CURRENT.md`, this reconciliation record, and the integrity-only manifest
change. No farm run is requested or claimed. Fresh user-supplied farm evidence
remains required before any E.8.1 runtime conclusion.
