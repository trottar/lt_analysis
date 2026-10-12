# E.8 current-snapshot contribution runtime closure - 2026-10-12

Status: `CLOSED / RUNTIME VALIDATED` only for Q4p4W2p74 internal detached
candidate coherence and the descriptive Left/lowe t0/t1/t2 aggregate/PDF
diagnostic (primary t0/t2, t1 reference).

## Acceptance, attribution and identities

RUNTIME VERIFIED (attributed): user-supplied Jefferson Lab materialization ZIP,
candidate JSONs, farm seven-input preflight receipt and contribution ZIP, with
prior independent ChatGPT review, own this narrow acceptance. The
[checkpoint contract](../phases/e8-current-snapshot-contribution-runtime-evidence-checkpoint-task-contract.md)
transfers the completed acceptance. Codex records it without farm replay,
independent raw-artifact validation or possession of the five F.1 inputs.
The seven-input locator is local-only and contains no public path authority.

- Farm-evaluated/reviewed Git HEAD: `75756d7b32767f18949386c4c5864c756d5e85f6` on `test`.
- Materialization ZIP basename: `kaonlt-current-f1-snapshot-518d9d7533904bd686d6d02e112dc851.zip`; ZIP SHA-256 `5c55452ced080514c8f910ebd491effc74f2f36ff1affa28afb26a900c0b300d` (independently recalculated), four-member archive CRC PASS.
- Materialization manifest schema: `method_a_current_f1_contribution_snapshot_materialization/v1`; fingerprint `a242bea92966795c0776ae370ce079e3ca57eea97051c04016ba350156bed6b3`. It reports `internal_snapshot_consistency=PASS`, cross-stage bindings, exact F.4 reproduction, persisted readback, and all 15 parent signed closures. Its F.3 construction authority uses an explicitly non-authoritative all-zero 40-character source-HEAD sentinel. `canonical_owner_equivalence=BLOCKED`; no production application or Method-B numerical dependence.
- New F.2/F.3/F.4 candidate raw SHA-256, in order: `11e8b4a46a3a6d78c67d51761db99bd7a22b137120f1f8a8dab03e09e895b0f8`; `5b098fc92b9107d3312944b3605be6cab85becf66f5ff36fc31aafaa1756c12d`; `946680dc447bf7904db5ed30d2881037350cb8b584429701af5628af61aa18cf`.
- Five current F.1 raw SHA-256 values (these **identify inputs; they are not accepted historical authority replacements**): `Left-lowe=99dcbf00a6ed93410ec0e38b317b25db94e0651d220a837662309d4be9940e6b`; `Left-highe=ced1b5aad28bf7bd071caf3af3ef045c6931824569af1aefc64ce8e6cc4d5585`; `Center-lowe=440f92711df824d7c8b9bf4846d9e2c0c59b8c85fd9cb3124ba512a502c5c936`; `Center-highe=40d9175bc1afaaa5cacc965c5ece88b3de08aacd4d20a07ea9bf248e3c21ef96`; `Right-highe=e2d6b229c8c92bc8698f9b95e392bcc77aa4f56cdbeb999702da6930948efda3`.
- Subsequent **farm** read-only seven-input preflight: `preflight=PASS`, `internal_snapshot_consistency=PASS`, all seven fingerprint validations true, `first_failure=null`, `outputs_written=false`, canonical equivalence still `BLOCKED`, snapshot fingerprint `3644eb1b91932519a6a92f08107e77350d558ef7c1bb376cbcd15781df7085e7`.
- Full contribution ZIP basename: `kaonlt-current-snapshot-contribution-c6392716c8244bdcb56807bb8d1342bc.zip`; ZIP SHA-256 `b053ce131d1c92b202dc0e945a4d539757064b8633035f31861035a1c7622fd1` (independently recalculated from supplied bytes); contains exactly `diagnostic.json`, `diagnostic.pdf`, `manifest.json`; ZIP CRC and relevant persisted numerical/PDF checks passed earlier ChatGPT review. Manifest `pion_hgcer_method_a_current_snapshot_contribution_debug_bundle/v1`, fingerprint `b0fe81bc528ea1df381241f4631b477a613793e448f20f096a7f125b3e0fdf1c`, diagnostic fingerprint `27dee190435495fbc9c59bb4476ec17f8b5db20e4a6b5fe7b3fffb39ddb88322`. Member hashes: `diagnostic.json=57e638f578e958b62e0d9a212cb4f1d14bb152eeff4c2da9aec8c3e30aad0b5b`, `diagnostic.pdf=26a58272710b1be73a52eee5f2b061cb17e63df22e9b941d23f84b13028a7a3d`. PDF is **eight detached pages**, renderer input is persisted aggregates only, previously visually inspected by ChatGPT. Farm diagnostic provenance reported a dirty checkout; identities of modified/untracked farm paths remain NOT VERIFIED, not silently described as clean.

The materialization and passed read-only preflight resolve predecessor
`snapshot_failed:neff_moment_mismatch` and
`snapshot_failed:f3_fingerprint_input_content_mismatch` only in this detached
path. Preserve the [materialization history](../phases/e8-current-f1-coherent-detached-materialization-task-contract.md),
[checker repair](../phases/e8-neff-moment-validator-parity-repair-task-contract.md)
and [debug history](../phases/e8-current-snapshot-method-a-contribution-debug-task-contract.md).
No earlier failed artifact is retroactively accepted.

## Descriptive quantities and scientific limits

SOURCE VERIFIED: tracked materializer
`testing/materialize_method_a_current_f1_contribution_snapshot.py` retains
non-authoritative F.3 construction and exact internal F.4 reproduction;
`testing/analyze_pion_hgcer_method_a_current_snapshot_contribution_debug.py`
reads five F.1 plus F.3/F.4, writes no preflight output, and renders only
persisted diagnostic aggregates. The detached validators preserve canonical
blocking and exclude production mutation and Method-B numerical dependence.
These source checks establish no new runtime acceptance.

Define `b=s*w0`, `b_rw=s*w0*C`, signed source coefficient `s`;
`sum b_rw=sum b`. For nonnegative support measures `abs(b)` and `abs(b_rw)`
separately, `N_eff=(sum abs(b))^2/sum b^2` (and the reweighted analogue).
Top shares below rank by actual `abs(b_rw)`, not by `C` alone. Neff and
shares describe support, not uncertainties or physical event probabilities.

RUNTIME VERIFIED (attributed): persisted uploaded diagnostic quantities from
prior ChatGPT numerical review, transcribed exactly from the contract:

| `Left-lowe` parent | t0 | t1 | t2 |
| --- | ---: | ---: | ---: |
| Application records | 11,019 | 19,809 | 24,552 |
| Signed baseline and adjusted parent total (`B≈B_rw`) | 0.108151759731 | 0.337472682094 | 0.106457792086 |
| `N_eff(abs(b))` | 4,912.9217 | 9,417.0855 | 13,193.3543 |
| `N_eff(abs(b_rw))` | 865.6855 | 5,265.9988 | 4,477.3820 |
| Retained fraction | 17.62% | 55.92% | 33.94% |
| Top 1% of `abs(b_rw)` share | 26.87% | 11.46% | 13.50% |
| Top 5% of `abs(b_rw)` share | 62.20% | 24.69% | 39.38% |
| Total absolute-support change `(A_rw/A−1)` | −0.269% | +1.035% | −3.035% |
| F.3 out-of-support fraction | 0.118% | 0.106% | 0.134% |

Prior review independently checked source/phi/MM/joint aggregation and
numerical moments. Significant parent-internal signed background relocation
across `(phi, analysis_MM)` occurs despite preserved signed parents. t0 is
especially concentrated; t2 also has substantial concentration. Diagnostic MM
uses 40 fixed bins over 0.8-1.4 GeV plus flows, separate from authoritative
Lambda yield integration and final `Y0/YA`. Parent preservation is not per-cell
conservation; tail concentration alone establishes no production failure or
correction error. F.3 `hgcer3` remains relative detector response, not absolute
pion-to-kaon mis-ID.

MEMORY ONLY: historical accepted closures are retained under their own records;
this checkpoint neither reruns nor reinterprets them. INFERENCE: mathematical
redistribution alone cannot attribute a downstream stage/yield change without
matching source-owned objects, populations, lineage and units.

NOT VERIFIED: modified/untracked farm-path identities, canonical-five F.4
scientific equivalence, new current-F.6.3 production closure, hardware cause,
absolute HGC mis-ID/proxy calibration, cross-setting effect, final kaon yield,
production correctness, and SIMC-unit proof.

## Preserved blockers and successor boundary

Private canonical-five Method A and final E.8/F.6.4 remain `BLOCKED` at
`f6_3_scientific_equivalence_mismatch:f4:$.parents[0].absolute_baseline_parent_sum`.
No canonical-owner retry or historical F.1/F.4 pin replacement. F.1 archaeology
remains `DEFERRED`. Absolute SIMC remains independently `BLOCKED` by
`SIMC_normfac_luminosity_and_charge_units_not_source_proven`. Method A stays
private/non-production; Method B stays independent diagnostic/cross-check only
and numerically excluded. `no_empirical_residual` and all science remain frozen.

[Baseline-five closure](e8-2-canonical-five-baseline-yield-runtime-closure-2026-10-09.md)
uses farm source `8e68d2acae91fd97622793df563aba0296403a48`, distinct from
this detached diagnostic's `75756d7b32767f18949386c4c5864c756d5e85f6`.
No same event population, numerical equivalence or direct subtractability is
assumed. [CURRENT](../CURRENT.md) owns the sole push-stable successor: read-only
baseline-to-detached-background scientific crosswalk after actual-diff acceptance,
user publication and pushed-state review, starting Left/lowe t0/t1/t2. Audit
matching keys, cuts/windows/bins, run/lineage provenance and weighting units,
and distinguish exact same-lineage comparisons from qualitative juxtaposition.
Report the first missing comparable source-owned object/accepted artifact;
invent no producer, reconstruct no Y0/YA and retry no blocked owner.
No farm command or successor scientific implementation is authorized here.
