# Canonical-five fresh-F1 lineage: failed gate and detached materialization

## Evidence and acceptance boundary

This checkpoint records the farm diagnosis and detached comparison/materialization
accepted in the user-supplied task. Codex did not rerun or independently validate
the farm artifacts. The first isolated Q4p4W2p74 canonical-five run completed
full analysis but failed the downstream E.8.4 page gate. All five page manifests
rendered only `e8_4.unavailable`, each with:

`f6_3_f4_shared_reproduction_failed:f3_fingerprint_input_content_mismatch`

Fresh F.1 raw/stable/training/application identities differed from the previously
pinned candidate lineage. Detached comparison established F.2 and F.3 scientific
equality and F.4 as the first changed stage. Detached materialization completed
without errors, accepted-authority mutation, production mutation or promotion.
Failed-run artifacts remain diagnostic-only, inadmissible for runtime closure.
Canonical-five runtime is NOT VERIFIED. Method A remains detached/non-production;
Method B remains diagnostic/cross-check only and numerically excluded.

## Exact reviewed identities

Materialization construction source head: `463d2657f696ecee33113edc3393ac51083a8944`.

| File role | Raw SHA-256 |
| --- | --- |
| comparison | `4c3fdab5d05965a1b3f1dc835c9d8c529896e6eaf2da8eb5c8ac8fc26e95257a` |
| manifest | `e7c37f55b24739e8be0fb778a6c3491e17df746263e6d9d71a4f31ec65cd7908` |
| f2 | `2fa715b2b5c2d5077e38416fb9eceb814e013e72f308c43526fddff36ddd962e` |
| f3 | `c5b86452b790ecbaf5b8f0df05da67efa2fa92aab157b12b153ed7e491a38228` |
| f4 | `79e7ceda7221cbeeead4ed5bc306b0e0e670741a27beaa980e22349c555e96d7` |

Candidate basenames:

- f2: `Q4p4W2p74_kaon_pion-background_hgcer-method-a-acceptance-representation-current-baseline-candidate.json`
- f3: `Q4p4W2p74_kaon_pion-background_hgcer-method-a-acceptance-map-current-baseline-candidate.json`
- f4: `Q4p4W2p74_kaon_pion-background_hgcer-method-a-parent-preserving-correction-current-baseline-candidate.json`

| Fresh F.1 setting | Raw SHA-256 |
| --- | --- |
| Left-lowe | `eb6f659da0511f0f6ea867420fda0ee6feb77e6698509ccd31f80a92cc541c07` |
| Left-highe | `544ea08f71e74b6b59bc33d05458d4acc051e96fb0327f54ec092a331245c01e` |
| Center-lowe | `2b193dec46b10aebc74dfb51634d19a7fef29362fd0a252d940f2894d18a0a16` |
| Center-highe | `c857911396bed03f9e418ca45a609509ebf1eadf8084d4c0ff2f9cc3a9980941` |
| Right-highe | `e03245494d54a6f812c426d189b220417fdc0888ba1f58a961fbe250ff097652` |

| Candidate | Fingerprint field | Value |
| --- | --- | --- |
| f2 | representation_fingerprint | `e9df6e09c31dda4d7d724792ce3c40125e27530c4c76e0dbe7f5b9861d0e6216` |
| f2 | artifact_fingerprint | `87d9a8780cb2de9f2830afd9e388064151f796654ddb5c4a8e1938443806a0c6` |
| f3 | map_fingerprint | `6042e9485827d9884d2a841d5c6cb5e23401e1026c8e2d4bd2022ede15843548` |
| f3 | algorithm_fingerprint | `ba29630b2f40a87cbadbe751504ce48f23e2b17a08378a2c8a219f93131cb912` |
| f3 | artifact_fingerprint | `e1cf0dd14644034f5d21df27d29f48355cde7c8e85c6123b709a6975c152a302` |
| f4 | correction_fingerprint | `71527eecd5828766ebcb3b9ef5e4e1931eb242fabfc1059f04e4eec104f6cc98` |
| f4 | artifact_fingerprint | `0b3899e5a3f24a34e04aa550ef162414c1d95af9a786e748612cbcb018dbcad4` |

F.3 candidate reconstruction farm-source sentinel: `0000000000000000000000000000000000000000`.
F.4 inherits the exact fresh F.3 raw hash and map/algorithm/artifact fingerprints.

## Materialization gates

```json
{
  "complete": true,
  "errors": [],
  "non_authoritative": true,
  "accepted_authority_mutated": false,
  "production_objects_mutated": false,
  "production_application_performed": false,
  "method_a_promoted": false,
  "f2_scientific_payload_match": true,
  "f3_scientific_payload_match": true,
  "f4_scientific_payload_match": false,
  "first_changed_stage": "F4"
}
```

## Local refresh boundary

The narrow implementation refreshes only private F.6.3 pins and adds owner
verification/staging of reviewed F.3/F.4. A five-setting preflight loads current
F.1 and staged candidate JSON and runs the real shared F.4/F.5 validators before
full analysis. Local synthetic checks are SOURCE VERIFIED only; actual-diff
review, user synchronization and a new applicable farm gate remain required.
Historical accepted authorities and all scientific mathematics remain frozen.
