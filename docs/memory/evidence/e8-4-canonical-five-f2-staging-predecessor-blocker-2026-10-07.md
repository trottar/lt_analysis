# Canonical-five F.2 staging-predecessor blocker

## Recording authority and scope

This record carries only supplied/accepted diagnostic facts from the
[repair contract](../phases/e8-4-canonical-five-f2-staging-predecessor-repair-task-contract.md).
Codex did not run the farm or independently inspect the raw farm artifacts.
The failed attempt supplies no canonical-five runtime or visual closure;
the full runtime gate remains `BLOCKED`.

## RUNTIME VERIFIED — supplied failed owner attempt

Attempt stem:
`KaonLT_E8_4_Fix5_CanonicalFive_Q4p4W2p74_20261007-231233_3914085`.

Farm-evaluated source: `946b5512bb00d22bb446bd96d6ab9ee35f5f989b`.

```text
status = failed
stage = candidate_staging
failure_reason = unknown_candidate_target:Q4p4W2p74_kaon_pion-background_hgcer-method-a-acceptance-representation-current-baseline-candidate.json
analysis_started = false
analysis_completed = false
artifact_verification_completed = false
collection_completed = false
zip_verification_completed = false
```

Analysis never started. Neither artifact acceptance, collection nor ZIP
acceptance occurred. This first failed invariant is distinct from the earlier
[second full-run provenance blocker](e8-4-canonical-five-second-full-run-provenance-blocker-2026-10-05.md).

## RUNTIME VERIFIED — supplied existing and reviewed F.2 identities

Observed target basename:
`Q4p4W2p74_kaon_pion-background_hgcer-method-a-acceptance-representation-current-baseline-candidate.json`.

Observed raw SHA-256:
`182433ccd2d13d3f7af80b963e69190ad8f2eb4e295c8216b36585e8ee7409b2`.

Observed schema:
`pion_hgcer_method_a_acceptance_representation_artifact/v1`.
Observed flags: `non_authoritative = true`, `probe_only = true`.

| Identity | Observed target | Reviewed candidate |
| --- | --- | --- |
| Raw SHA-256 | `182433ccd2d13d3f7af80b963e69190ad8f2eb4e295c8216b36585e8ee7409b2` | `2fa715b2b5c2d5077e38416fb9eceb814e013e72f308c43526fddff36ddd962e` |
| Representation fingerprint | `89976fd7a94a608be81a9bf0e177c5c59b5b5d6e6f39b5c704c0e719ce877585` | `e9df6e09c31dda4d7d724792ce3c40125e27530c4c76e0dbe7f5b9861d0e6216` |
| Algorithm fingerprint | `2ca2cc03ae9fae90c50c620f267d467a213a819c1e724a1b731bd88e2a0184bd` | `2ca2cc03ae9fae90c50c620f267d467a213a819c1e724a1b731bd88e2a0184bd` |
| Artifact fingerprint | `c63e7461240622fc64e9a2d566c687deb5fb620b54b0d02f99dc5485a8d2d400` | `87d9a8780cb2de9f2830afd9e388064151f796654ddb5c4a8e1938443806a0c6` |

The observed provenance was non-identifying:

```json
{
  "generated_at_utc": null,
  "git_head": null,
  "git_status_short": null,
  "input_paths": {}
}
```

The farm read-only comparison rechecked the reviewed materialization hash:
`reviewed_sha256_expected = true`.

## RUNTIME VERIFIED — supplied exact scientific comparison

The farm compared the observed and reviewed F.2 bodies using the pushed
source-owned `scientific_projection`, `first_mismatch` and `F2_PROVENANCE`
in `src/cuts/pion_hgcer_method_a_parallel_full_procedure.py`:

```text
f2_scientific_payload_match = true
first_mismatch = null
```

## Decision and retained boundaries

The supplied exact comparison authorizes the byte-pinned `182433...` state
only as a replaceable staging predecessor. It is not accepted scientific
authority, a known historical authority, reviewed runtime candidate,
materialization input authority or fallback F.2. Staging must replace it with
exact reviewed `2fa715...` bytes; every other unknown target still fails closed.
No scientific-comparison-at-staging fallback is authorized.

Canonical-five runtime remains `BLOCKED`. The preceding equivalence repair
remains `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; no farm closure follows
from this diagnosis or local owner repair. Existing closed scopes remain
unchanged. Final E.8/F.6.4 and absolute-SIMC interpretation remain `BLOCKED`.
Method A remains detached/non-production; Method B remains diagnostic/cross-check
only and numerically excluded. CURRENT alone owns the conditional next gate.
