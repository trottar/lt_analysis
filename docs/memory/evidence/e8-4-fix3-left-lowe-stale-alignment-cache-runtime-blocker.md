# E.8.4 Fix.3 Left/lowe stale-alignment-cache blocker

## Status

`BLOCKED` — supplied Jefferson Lab farm comparison exposes a new E.8.4/F.6.3
provenance blocker. It is not E.8.4 runtime closure or an accepted replacement
for frozen F.3/F.4 authority.

## Supplied farm and artifact observations

- The fresh `Q4p4W2p74 / Left / lowe` ordinary low-epsilon analysis completed.
  E.8.4 was unavailable with
  `f6_3_f4_shared_reproduction_failed:f3_fingerprint_input_content_mismatch`.
- Accepted F.3 file SHA-256 remains
  `04a9a576767ba3f81c8c03daab2461998894c655f0d875bd532ad30cb46c0d95`;
  accepted F.4 file SHA-256 remains
  `adcc01900b8adc211e296aaedaa314fdb86207d51483ebfd4b3edf69f558f188`.
- Four canonical settings reproduce accepted F.1 raw and stable identity. Only
  `Left / lowe` differs. Its frozen/fresh F.1 raw SHA-256 values are
  `f16c5e89f62e848ec221335bd4916265757040af9970b005d8f67cc800e0077d`
  and `10086d7d4c42389c9fdd16471c49d59cd189a30980767bcdb87b85914169ef95`.
  Its frozen/fresh stable F.1 fingerprints are
  `76bcf74ac1de6c19d402e078bc8b8ac7f1c0f821e6d7fda1ef97338541d7c988`
  and `250d95f3f210e3bca7315439777a71adbdb1c21c5c74bfeab36490e432f7dfee`.
- Unchanged F.1 partitions: `part1_config_fingerprint`,
  `method_a_event_population_fingerprint`,
  `acceptance_feature_metadata_fingerprint`, and
  `application_child_assignment_projection_fingerprint`.
- Changed F.1 partitions: `phase_a_contract_fingerprint`,
  `phase_a_pion_event_population_fingerprint`, `coordinate_fingerprint`,
  `method_a_fingerprint`, `method_a_training_population_fingerprint`, and
  `application_population_fingerprint`.
- Frozen/fresh `Left / lowe` coordinate fingerprints are
  `207456612e2f4b88daa8de7d89b55e0a56d52e01096608fb3ef7d9316d19f2b4`
  and `e41d4d004e047076c213cab7bedf5a51fe56f6f7b83c3c9e0fd5e2a2e52a4066`.
  The fresh value matches the September 29 pre-Fix.3
  [alignment blocker](e8-4-left-lowe-pion-alignment-determinism-blocker.md).

## Source audit and limit

At committed starting HEAD `8d62dcfdf2fd8d08298c075940b4ab1b28d7f079`,
`load_or_resolve_pion_component_alignment` returns a compatible persisted
record before running the candidate scan. Its cache compatibility checks
schema/config/scope/identity/axis/parent/templates/pion-control identity but
did not include the Fix.3 scan-acceptance semantic revision. Thus a pre-Fix.3
schema-v2 record can bypass the repaired scan predicate. This is source proof
combined with the supplied farm comparison; the fresh farm output did **not**
directly print `persistence_status=reused`. The accepted F.3/F.4 checks must
remain fail closed. A new narrow farm gate is required after Fix.4 review and
user-controlled source/provenance reconciliation.
