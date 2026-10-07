# E.8.2 Left/lowe runtime/scientific closure — accepted 2026-10-06 evidence

## Authority and scope

Status: `CLOSED / RUNTIME VALIDATED` only for Q4p4W2p74 / Left / lowe.
This record transcribes accepted supplied farm artifacts and prior ChatGPT
runtime, visual and scientific review specified in the
[closure contract](../phases/e8-2-left-lowe-runtime-scientific-closure-task-contract.md).
Codex records that acceptance; it did not rerun the farm or independently
inspect/hash/render the external artifacts. Identities below are the exact
accepted supplied identities, not newly measured local values.

## Accepted fresh runtime evidence


Package stem:

```text
KaonLT_E8_2_Left_lowe_scientific_audit_20261006-233426
```

Farm-evaluated source:

```text
706ae708ae69be0585472982a4df3d87c63f186f
```

Scope:

```text
Q4p4W2p74 / Left / lowe only
```

External artifact identities accepted by ChatGPT:

```text
ZIP
SHA-256:
464d0db14af4e09a53a418dc077136fa17f1efd1c4c36b57758247b8f083af23
bytes:
30272198

gate-status JSON
SHA-256:
900f32f5077754df27439c2df3c4243c3d4e9fc8244a399943616153a4089c1c
bytes:
21073

run-summary JSON
SHA-256:
d9d265bcb80f82f4bb3fff7abb243913abb2b20292b71a5d7b58f8796087400e
bytes:
17394

child log
SHA-256:
4e8713d8f4c399063a6feb43395caf1a581f6154c06660a3888b81561f297778
bytes:
51888750
```

Accepted packaged artifacts:

```text
Left_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction.pdf
SHA-256:
aacac5daa935eb1f2af6dacdf0ddfdcdd4f5b202e0d54fc11b8edd254576e61e
bytes:
1551786
pages:
74

Left_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction-manifest.json
SHA-256:
e51fe03cff97ce00e05c850b031206e97a9cb9d15b0c71279618fda1a9577f02
bytes:
18203

kaon_FullAnalysis_Q4p4W2p74_lowe.json
SHA-256:
6e6f30c630de92a138520bfe3a93612357d8b176211454c97f22bbb94b9e328e
bytes:
140864222

kaon_FullAnalysis_Q4p4W2p74_lowe_correction_ledger_no_empirical_residual.json
SHA-256:
cba9737c69868b8cb6a9ecbfba08bc6a25a7e326445298fd859cc52f403b9f92
bytes:
53424

kaon_FullAnalysis_Q4p4W2p74_lowe_correction_ledger_no_empirical_residual.csv
SHA-256:
3eb3aa22419e121140edb5829084f8401e9986cde91790d8c221751380c63278
bytes:
20343
```

RUNTIME VERIFIED by the accepted supplied evidence and prior ChatGPT review:

- owner `status=success`;
- owner `stage=complete`;
- null failure reason;
- child analysis return code 0;
- analysis completed;
- artifact verification completed;
- collection completed;
- ZIP verification completed;
- worktree cleanup completed;
- repaired disposable-worktree OUTPUT link passed;
- ordinary farm checkout preservation passed before and after;
- installed `ltsep` preservation passed before and after;
- source commit and `origin/test` matched
  `706ae708ae69be0585472982a4df3d87c63f186f`;
- only the two established allowlisted farm outputs were dirty in the ordinary
  checkout;
- all 18 required E.8.2 page IDs were structurally present;
- all 18 required E.8.2 rendered pages, PDF pages 47–64, were visually reviewed
  and accepted.

The repaired visual defects are closed:

- all nine phi children are physically visible on every subtraction page;
- no top-row clipping remains;
- final-MM headings/t context do not overlap panels;
- stage-yield text is fully readable;
- no right-edge truncation remains;
- no new renderer failure is visible.

Therefore:

```text
E.8.2 Q4p4W2p74 / Left / lowe:
CLOSED / RUNTIME VALIDATED
```

This is a narrow Left/lowe closure only. It does not establish canonical-five
E.8.2 closure, other-setting runtime acceptance, Method-A promotion,
absolute-SIMC amplitude correctness, or final E.8 closure.

---

## Accepted E.8.2 scientific interpretation

The accepted stage-yield pages label their stage integrals explicitly as
**diagnostic Lambda-window integrals, not final extracted yields**.

For each canonical t parent, summing the nine phi cells gives the following
signed fractional change relative to the immediately preceding stage:

| t parent | Random | Dummy | Slow proton | `prune_hist` | Baseline pion |
| --- | ---: | ---: | ---: | ---: | ---: |
| t1, 0.4000–0.5667 GeV^2 | -0.68% | 0.00% | -7.75% | 0.00% | -20.95% |
| t2, 0.5667–0.7333 GeV^2 | -0.64% | -1.48% | -6.46% | 0.00% | -18.71% |
| t3, 0.7333–0.9000 GeV^2 | -0.44% | -2.59% | -8.08% | -0.03% | -4.81% |

Accepted interpretation:

- random subtraction is small in this scope, below 1% for all three parents;
- dummy subtraction is modest, reaching about 2.6% at t3;
- slow-proton cleaning is a consistent, material correction of about 6.5–8.1%;
- baseline pion subtraction is the largest diagnostic-window reduction at t1
  and t2, about 19–21%;
- at t3 the baseline pion effect is about 4.8%, smaller than the slow-proton
  effect;
- production `prune_hist` is essentially inert in the Lambda-window diagnostic.

The only visible nonzero prune change occurs at t3 / phi5:

```text
PRE  = 0.00018539
POST = 0
```

This is about 0.03% of the t3 pre-prune signed parent sum and about 0.014%
across all three pre-prune parent sums.

Preserve the interpretation boundary:

- `PI` is the diagnostic Lambda-window integral of the after-pion histogram;
- authoritative `Y0` is the stored final extracted yield;
- differences between `PI` and `Y0` are not a new subtraction stage;
- no empirical-residual stage is present under the active
  `no_empirical_residual` profile;
- small near-zero or negative final cells are observations of the baseline
  result, not by themselves evidence for a production failure or a redesign.

Do not broaden these conclusions beyond Q4p4W2p74 / Left / lowe.

## Retained boundaries

Scientific/runtime source and production physics are unchanged by this closure.
Method A remains detached/non-production; Method B remains diagnostic/cross-check
only and numerically excluded. Canonical-five provenance/identity repair remains
`DEFERRED`. Absolute-SIMC amplitude interpretation remains separately `BLOCKED`;
final E.8 and F.6.4 remain `BLOCKED`. No canonical-five or other-setting closure,
Method-A promotion, absolute-SIMC claim or final E.8 closure follows.
