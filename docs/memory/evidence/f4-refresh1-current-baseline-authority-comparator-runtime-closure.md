# F.4.Refresh.1 current-baseline authority comparator runtime closure

## Direct supplied farm evidence

The user supplied the result of the detached farm comparator as `KaonLT_F4_Refresh1_authority_comparison_Q4p4W2p74_20260930-134654.json`, SHA-256 `c16287e9192288ed5f116b5eebf95f49261d7fd09755ee923070154e047cf7b5`. The JSON itself is not in the local worktree for this task; the following are the supplied result fields, not a new local or farm execution.

Accepted raw JSON SHA-256: F.2 `87e0a16870ea34a6611702e1bd4fc02c4dc6a8daa132e75b725a1b782a9a31da`, F.3 `04a9a576767ba3f81c8c03daab2461998894c655f0d875bd532ad30cb46c0d95`, and F.4 `adcc01900b8adc211e296aaedaa314fdb86207d51483ebfd4b3edf69f558f188`. Current Left/lowe F.1 raw SHA-256 was `10086d7d4c42389c9fdd16471c49d59cd189a30980767bcdb87b85914169ef95`; the other four F.1 raw identities remained the accepted F.2 inputs.

| Comparator result | Value |
| --- | --- |
| F.2 scientific payload match | `true` |
| F.3 scientific payload match | `true` |
| F.4 scientific payload match | `false` |
| First changed stage | `F4` |
| Candidate F.2 serialized SHA-256 | `182433ccd2d13d3f7af80b963e69190ad8f2eb4e295c8216b36585e8ee7409b2` |
| Candidate F.3 serialized SHA-256 | `eeb480d9b7ddffab1e97c7f9f27c3099d0780c2bbb4d42a0ee68c077a5752b3d` |

The first F.4 scientific mismatch was `$.parents[0].absolute_baseline_parent_sum`: accepted `0.12294577672157099`, candidate `0.12072578661296547`.

| Global maximum | Absolute difference | Relative difference |
| --- | ---: | ---: |
| Baseline parent sum | 0.0049304592644921486 | 0.044263728042183405 |
| Parent normalization | 0.06654204106342765 | 0.033455174261301984 |
| Correction factor summary | 1.1986283985416364 | 0.03237215807179576 |

Only the three Left-lowe canonical-t parents differed numerically in the F.4 detail. The `t=1` difference was at floating-point scale; `t=0` and `t=2` carried the substantive current-baseline changes.

## Interpretation and closure boundary

F.4.Refresh.1 and F.4.Refresh.1.Fix.1 are `CLOSED / RUNTIME VALIDATED` **only for the detached current-baseline comparison gate**: accepted F.2 and F.3 scientific payloads reproduced exactly, while F.4 was the first scientifically changed stage. The earlier Fix.1 list-recursive reporting is included in that gate.

Accepted F.2 and F.3 science do not need a scientific refresh on this evidence. Their current-baseline provenance lineage still needs materialization because fail-closed F.3/F.4 input identities bind to current F.1. Accepted historical F.4 remains valid for its accepted historical baseline. A distinct, non-authoritative current-baseline F.4 candidate must be materialized and reviewed before any source-owned authority decision.

This comparison does not replace accepted F.2/F.3/F.4 artifacts, validate a refreshed F.4 or F.5 authority, close F.6.3/E.8.4 runtime, alter production physics, or promote Method A.
