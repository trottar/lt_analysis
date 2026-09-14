# F.1.Fix.5 v4 bundle inspection — 2026-09-14

## Supplied evidence

- Archive: `KaonLT_F1Fix5_validation_Q4p4W2p74_full_v4.zip`.
- Archive SHA-256:
  `7bf7a4890eeca3ae00b43bebd245ee3f6ed658cf8ef223c5b3447a70eb694292`.
- Embedded `manifest.json` SHA-256:
  `54df07577fb44f2322158ca0adade03e81bde6f74c84d445e30bde1b6d2bfe2e`.
- Manifest schema/profile: `pion_hgcer_validation_bundle/v4` /
  `phase_f1_batched_farm_acceptance_review/v2`.
- Recorded farm HEAD: `126fa22c19bd29b9952f55b33ab43d59f9727ef6`.
- Required reviewed analysis commit:
  `dc4fc6283001739a487ec80068f951b0e388cae6`, recorded as an ancestor.

## Independent artifact audit

The archive contains exactly Left-lowe, Left-highe, Center-lowe,
Center-highe, and Right-highe for `Q4p4W2p74`. For every archived F.1 JSON,
direct inspection confirmed its recorded SHA-256; v2 schemas; detached flags;
frozen feature list and setting identity; exact positive-NPE training and
NPE>2 application selections; prompt application identity closure; all stored
population/feature/projection/full fingerprints; and all 30 stored canonical
`(t, delta)` Method-A and partition closures.

Each F.1 review PDF hash-matches its manifest entry and contains five PDF page
objects. Its source page manifest contains each required F.1 page ID once in
the required order, with `scope=setting` and `authoritative=false`. A local
PDF renderer was unavailable, so this is structural rather than visual PDF
inspection.

## Initial gate result

**BLOCKED** — the supplied archive alone is not F.1 runtime closure.

The bundle manifest records `complete=false` with one error:
`git_diff_check_required_analysis_commit_range` returned 2. Its output lists
only trailing-blank-line warnings in allowed `docs/memory/` files; the
committed range otherwise contains only `docs/memory/` plus the three exact
F.1 validation collector/profile/test files.

More importantly, captured farm `git status --short` was not clean: modified
`src/models/xmodel_kaon_pl.f` and `src/models/xmodel_pion_pl.f`, plus untracked
`src/kaon/functions/Q4p4W2p74.model` and
`src/kaon/xsects/unsep_Q44W274.csv`. The archive does not prove whether those
changes predated the F.1 outputs. Do not infer clean scientific-source
provenance from the matching committed HEAD.

The serialized F.1 v2 artifacts close, but the evidence did not satisfy the
complete clean-source farm gate needed to mark F.1.Fix.5
`CLOSED / RUNTIME VALIDATED` or to begin F.2.

## Source-gate repair

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING — The collector now preserves
the source-identity and path allowlists while running its historical range
whitespace check with only the allowed `docs/memory/**` pathspec excluded.
Analysis files, collector/profile/test files, and unapproved paths remain in
that check and in the independent committed-file identity audit.

Local validation passed:

- `python -m py_compile testing/collect_pion_hgcer_validation_bundle.py
  testing/test_collect_pion_hgcer_validation_bundle.py`;
- collector suite: 18 tests; and
- frozen F.1 contract/runtime suite: 10 tests.

The user explicitly chose a collector-only clean revalidation of the existing
farm outputs. Do not rerun the physics analysis for the allowed-documentation
whitespace issue. A fresh complete collector manifest from a clean checkout is
still required before closure.

## NEXT

Run the repaired collector only, from a clean farm checkout, against the
existing five `Q4p4W2p74` outputs. Inspect its complete manifest and retained
artifact/PDF provenance, then close F.1 or make one coherent repair. Do not
begin F.2 until that result is reviewed.

```tcsh
cd /group/c-kaonlt/USERS/trottar/lt_analysis
python3 testing/collect_pion_hgcer_validation_bundle.py \
  --outdir /lustre24/expphy/volatile/hallc/c-kaonlt/trottar/OUTPUT/Analysis/KaonLT \
  --kinematic Q4p4W2p74 \
  --output /lustre24/expphy/volatile/hallc/c-kaonlt/trottar/OUTPUT/Analysis/KaonLT/KaonLT_F1Fix5_validation_Q4p4W2p74_full_v4_clean-collector.zip
```

Use a new output name; this must not overwrite the supplied archive.
