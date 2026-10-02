# E.8.4 Fix.5 Left/lowe post-farm identity, lineage and visualization — 2026-10-02

## Scope and evidence authority

`Q4p4W2p74 / Left / lowe`, pushed/farm source
`f9d70732290ea461096374ca1270b47452644991`.
The [checkpoint contract](../phases/e8-4-fix5-4-post-farm-identity-audit-memory-checkpoint-task-contract.md)
supplies the direct farm observations and independent review findings below.
Codex records them; it did not independently open the farm PDF/manifest or
perform a numerical audit in this memory-only task. Local source excerpts
confirm the different historical/current input pins. No new bundle hash or
runtime acceptance is inferred from source or the supplied render observations.

## Direct observations

The expensive full-analysis/render stage completed in the ordinary canonical
artifact directory:
`/group/c-kaonlt/USERS/trottar/lt_analysis/OUTPUT/Analysis/KaonLT`.

- `Left_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction.pdf`:
  farm-visible timestamp `2026-10-01 23:34`, approximately 2.2 MB.
- `Left_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction-manifest.json`:
  timestamp `2026-10-01 23:34`; `page_count = 97`; `renderer_failures = []`.
- Required `full_background.e8_4.*` pages were present:
  `method_a_vs_simc.t1/t2/t3`, `baseline_method_a_simc.t1/t2/t3`,
  `yield_summary.t1/t2/t3`, and `parent_closure`. Prior E.8.4 pages remained.
- User-created inspection-only Globus copies:
  `KaonLT_E8_4_Fix5_Left_lowe_20261001-234439.pdf` and
  `KaonLT_E8_4_Fix5_Left_lowe_20261001-234439-manifest.json`.

Fresh t1 pages display near-overlapping baseline/Method-A MM spectra with
substantially different stored scalar yields in several populated children:

| t parent | phi interval [degrees] | Displayed Y0 | Displayed YA |
| --- | --- | ---: | ---: |
| t1 | [-180, -140) | 0.033906 | 0.025165 |
| t1 | [140, 180) | 0.029873 | 0.016411 |

These are observation values, not universal constants or proof of which object
is wrong. Some populated cells also show a large data/SIMC amplitude difference.
Its physical meaning remains unresolved until object, normalization and unit
identity close.

## Confirmed lineage finding

**CONFIRMED:** sequential E.8.3 and E.8.4 pages describe different Method-A
lineages. E.8.3 presents the historical accepted persisted F.6.1 aggregate;
its accepted raw input SHA-256 identities are:

| Historical E.8.3 input | SHA-256 |
| --- | --- |
| F.4 | `adcc01900b8adc211e296aaedaa314fdb86207d51483ebfd4b3edf69f558f188` |
| F.5 | `143e3af6b1c69560e5bf351155b570d05f19e7e0ceb14665f5448c351ad501be` |
| F.6.1 | `62bad2dae0f65f1fff87eeffd071c29dbbd4cdc15a43f37d6d9853cb58b25bd6` |

The active F.6.3 branch consumes the separately validated current-baseline
candidate F.4 raw input:
`1d545924eba89c7f9ffa28028e307aca9b434a89beec06863cf2893887b6b902`.
Candidate F.3 remains
`eeb480d9b7ddffab1e97c7f9f27c3099d0780c2bbb4d42a0ee68c077a5752b3d`.

Source confirmation at the checkpoint HEAD: the `E8_3_F4_INPUT_SHA256`,
`E8_3_F5_INPUT_SHA256`, and `E8_3_F6_1_INPUT_SHA256` constants in
`src/cuts/full_background_subtraction_plots.py` pin the historical inputs;
`F6_3_CANDIDATE_F4_VALIDATION_AUTHORITY_BY_KINEMATIC` in
`src/cuts/pion_hgcer_method_a_parallel_full_procedure.py` pins the candidate.
The historical F.4 correction fingerprint is
`362241005c02f2149e260c391b5c3d35793287573128b42cf5ed693419d9d2f3`;
the candidate correction fingerprint is
`bce0cc12ef283b0c361818b436c4e190fc4ea5165436f6c36912cfbb6de53368`.

E.8.3 must not be interpreted as the aggregate explanation of current
F.6.3/E.8.4 without exact lineage identity. This cross-page provenance and
presentation mismatch does not invalidate historical E.8.3 `SOURCE REVIEWED`
work or replace accepted F.6.1 authority.

## Unvalidated numerical checks and hypotheses

Scientific interpretation is `BLOCKED`. Missing checks are not demonstrated
numerical failures. The following four audit requirements remain separate.

### 1. Pion-template to final-MM algebra

**NOT YET VALIDATED:** for every populated canonical child, prove bin by bin:

```text
MM_A - MM_0 = -(B_pi_A - B_pi_0)
```

Use the exact current-branch objects with compatible axes and a documented
floating-point tolerance. Current plots expose both sides without a fail-closed
numerical equality gate. No algebraic defect is asserted here.

### 2. Histogram to scalar-yield closure

**NOT YET VALIDATED:** for every populated canonical child, prove:

```text
Integral(MM_0) = Y0
Integral(MM_A) = YA
Integral(MM_A - MM_0) = YA - Y0
```

Use exact branch histograms and producer-owned integration semantics: identify
the Lambda window, selected bins/edge convention, underflow/overflow treatment,
any producer-owned bin-width or normalization factors and units. Do not choose
a plotting integral merely because it matches the displayed scalar. Trace the
stored scalar and histogram provenance through each handoff. The observed t1
visual/scalar inconsistency requires this audit; it does not prove whether a
producer, sidecar, payload or renderer object is incorrect.

### 3. Signed-cancellation diagnostics

**NOT YET VALIDATED:** report each child's signed integral, positive-bin
support, negative-bin support and absolute support separately, with explicit
window and sign conventions. The spectra contain signed random/dummy/background-
subtracted content. Those diagnostics distinguish an object mismatch from a
large fractional signed-yield change caused by positive/negative cancellation.
Cancellation remains a hypothesis; do not assert it without those numbers.

### 4. SIMC source, normalization and unit identity

**NOT YET VALIDATED:** the renderer clones existing per-child support from
`hist["_xsect_support_simc"]["mm"]` without renderer-side renormalization.
That source fact alone does not prove that the plotted child is the exact
same-cell authoritative yield-chain SIMC object in the same normalization and
units as `MM_0/MM_A`. Trace its production, existing normalization, handoff,
clone and display; record provenance and relevant integration-window integral.
Verify canonical t/phi identity and MM geometry. The visual amplitude difference
is not yet an established physics discrepancy or proof of incorrect SIMC scaling.

## Exact planned audit

Under a separate reviewed implementation contract, trace current-lineage
producer -> sidecar -> payload -> renderer for each populated canonical child.
Inventory all 27 canonical cells, including explicit empty/unavailable cases.
Carry source/artifact/fingerprint, setting, t/phi and axis identity through that
trace; compare actual inputs and stored values rather than page adjacency.
Report the four checks above with signed residuals, declared tolerances and
per-cell outcomes. Implement only numerical invariants or narrow repairs
warranted by the evidence, preserving producer-owned science. This record does
not implement or authorize those source edits.

## Confirmed visualization deficiency and dependency

**CONFIRMED presentation deficiency:** in many baseline/Method-A/SIMC panels,
the blue baseline is almost completely obscured by the magenta Method-A curve,
making near-equality and small differences hard to evaluate independently.

Visibility repair is dependency-blocked until numerical identity closes.
Future presentation changes may adjust draw order, line style/width, markers
and their size/frequency, or legend clarity. A difference/ratio support panel
is permissible only with already-validated stored objects and explicitly
presentation-owned arithmetic. Never alter histogram contents, normalization,
cuts, binning, fits, Method-A factors, SIMC scaling, yield extraction or
production objects to improve appearance.

## Farm/package caveat

The tracked owner did not return the final ZIP. Render completion and the
inspection copies are not accepted bundle closure or proof of complete
run -> verify -> package integration. The exact post-render operational cause
remains undiagnosed and separate from the numerical/visual findings. Diagnose
that packaging failure from existing artifacts/logs under its own warranted
scope; do not require an expensive analysis rerun merely to explain it.
No further farm run is authorized by this checkpoint.

The earlier accepted da38444 Left/lowe package and its post-run model-output
provenance caveat remain intact in the [earlier evidence record](../evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md).
This new observation supplies no replacement accepted ZIP identity.

## Scientific boundaries and retained statuses

F.4.Refresh.2 and historical F.6.1/F.6.2 closures remain unchanged. F.6.3's
current-baseline Left/lowe mechanics remain `CLOSED / RUNTIME VALIDATED` for
branch execution, live-cache parity, real child changes and signed parent
preservation. The existing narrow E.8.4 accepted scope is not enlarged by this
render observation. Historical E.8.3 remains `SOURCE REVIEWED`.
Method A remains detached/non-production; Method B remains diagnostic-only
and numerically absent. Final canonical-five E.8 and F.6.4 remain `BLOCKED`.
No authority replacement, scientific calculation or runtime behavior changes.

## Exact next sequence

[Fix.5.4](../phases/e8-4-fix5-4-current-lineage-identity-audit.md) is `ACTIVE`.
The sole substantive action is the current-lineage audit/fix after this
checkpoint's independent review, user commit/push and pushed-state review, and
a separate implementation contract. Dependency order:

1. this memory checkpoint -> ChatGPT PASS -> user commit/push -> pushed-state review;
2. separate Fix.5.4 numerical audit/fix contract;
3. Codex numerical implementation;
4. ChatGPT actual-diff review;
5. user commit/push;
6. pushed-state review;
7. only then write/run the visualization-only contract from that reviewed pushed numerical source, after numerical closure;
8. Codex visualization implementation;
9. ChatGPT actual-diff review;
10. user commit/push;
11. pushed-state review;
12. one narrow Q4p4W2p74 / Left / lowe farm run;
13. fresh scientific and visual evidence review.

Stop this checkpoint before audit implementation, visualization implementation,
commit, push or farm execution. CURRENT owns the exact ordinary NEXT.
