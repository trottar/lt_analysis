# Canonical-five current-lineage preflight runtime closure

## Status and acceptance boundary

The separate Q4p4W2p74 canonical-five `--lineage-preflight-only` farm gate is
`CLOSED / RUNTIME VALIDATED` for its lineage/pre-analysis scope only.

This checkpoint records the user-supplied receipt and its prior accepted
evidence review, as specified by the
[closure contract](../phases/e8-4-canonical-five-lineage-preflight-runtime-closure-task-contract.md).
Codex records that acceptance; it did not rerun the farm or independently
revalidate the raw receipt. Full canonical-five validation remains
`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; full runtime is NOT VERIFIED.

## RUNTIME VERIFIED — supplied accepted farm evidence

Returned receipt basename:

`KaonLT_E8_4_CanonicalFive_LineagePreflight_Q4p4W2p74_20261005-182635_1227600-gate-status.json`

Receipt SHA-256: `3bad694fa000e34d502b73d9124309c24884a5d6a5746ef9a863cb65e6238ba4`.

Farm-evaluated source commit: `b8bafd3f2853523ca9deb7aa7d6b572357584d6f`.

The receipt records:

```json
{
  "status": "success",
  "stage": "success",
  "mode": "lineage-preflight-only",
  "canonical_five_runtime_validation": false,
  "analysis_started": false,
  "analysis_completed": false,
  "artifact_verification_completed": false,
  "collection_completed": false,
  "zip_verification_completed": false,
  "failure_reason": null,
  "worktree_cleanup_completed": true
}
```

Reviewed fresh-F1 materialization verification and candidate F.3/F.4 staging
passed. The observed staged raw identities were:

| Candidate | Raw SHA-256 |
| --- | --- |
| F.3 | `c5b86452b790ecbaf5b8f0df05da67efa2fa92aab157b12b153ed7e491a38228` |
| F.4 | `79e7ceda7221cbeeead4ed5bc306b0e0e670741a27beaa980e22349c555e96d7` |

The five observed F.1 raw identities and transient factor populations were:

| Exact setting ID | F.1 raw SHA-256 | Factor count |
| --- | --- | --- |
| Left-lowe | `eb6f659da0511f0f6ea867420fda0ee6feb77e6698509ccd31f80a92cc541c07` | 55380 |
| Left-highe | `544ea08f71e74b6b59bc33d05458d4acc051e96fb0327f54ec092a331245c01e` | 15134 |
| Center-lowe | `2b193dec46b10aebc74dfb51634d19a7fef29362fd0a252d940f2894d18a0a16` | 52518 |
| Center-highe | `c857911396bed03f9e418ca45a609509ebf1eadf8084d4c0ff2f9cc3a9980941` | 24967 |
| Right-highe | `e03245494d54a6f812c426d189b220417fdc0888ba1f58a961fbe250ff097652` | 40254 |

These are runtime population counts for this gate, with no scientific
interpretation of their magnitudes.

The real five-setting F.6.3 `reconstruct_transient_factor_map` preflight passed.
For every setting, F.5 authority validation and shared F.4 persisted/recomputed
reproduction passed; `f4_shared_reproduction_passed=true`. Accepted F.3
source/map/algorithm/artifact identities and accepted F.4
source/correction/artifact identities matched the private current-lineage
candidate authority records. Historical general authorities were not replaced.
Each `accepted_f4_runtime_authority.accepted_authority_match=true` and
`current_baseline_candidate_lineage=true`. Both global
`f6_3_lineage_preflight.passed` and
`f6_3_lineage_preflight.f4_shared_reproduction_passed` were true.

No full analysis ran: `analysis_started=false`. Every setting recorded
`production_promotion_performed=false` and `event_correction_persisted=false`.

## Runtime isolation and preservation

The detached analysis worktree used the exact pushed source above.
All three LTANAPATH isolation probes passed. The external SIMC-link preflight
passed without required mutation. Installed ltsep preservation and final
ordinary-checkout preservation passed, including unchanged pre-existing
farm-local checkout state. Worktree cleanup completed. No machine-specific
absolute filesystem paths are recorded here.

## Explicit exclusions and retained scientific boundaries

This preflight does not validate:

- `run_Prod_Analysis.sh` or full `main.py` under the refreshed lineage;
- low/high completion markers or E.8.4 page production;
- PDF rendering, full-analysis JSON or correction ledgers;
- artifact freshness, collector execution or ZIP verification;
- full canonical-five runtime closure or production promotion;
- final E.8, F.6.4 or absolute-SIMC interpretation.

Method A remains detached/non-production. Method B remains independent,
diagnostic/cross-check only and numerically excluded. Scientific mathematics,
production subtraction, historical authorities and the active
`no_empirical_residual` profile remain unchanged.

The
[earlier failed gate/materialization record](e8-4-canonical-five-fresh-f1-lineage-materialization-2026-10-05.md)
retains the first full-run failure as historical diagnostic evidence. Acceptance
of this later narrow preflight does not retroactively accept those artifacts.

## Workflow consequence

The separate preflight gate and its evidence review are consumed. After
actual-diff review, user commit/push and pushed-state synchronization of this
checkpoint, the next operation is the final farm-readiness audit of the normal
canonical-five owner at the pushed source. Only a PASS may authorize one full
Q4p4W2p74 canonical-five farm run; this checkpoint does not authorize that run.
