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

## Historical initial gate result

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

C1 `31dd034d8404e317863bf0933c51253e1d3deeb8` repaired the collector without
changing analysis code, profile semantics, F.1 artifacts, ROOT paths, or
physics. C1 is a validation-infrastructure commit, not a runtime commit. It
retains the global worktree `git diff --check`, but scopes only the historical
committed-range whitespace check to these exact profile-owned files:

- `testing/collect_pion_hgcer_validation_bundle.py`;
- `testing/test_collect_pion_hgcer_validation_bundle.py`; and
- `testing/pion_hgcer_validation_bundle_profile.json`.

The independent committed-file identity audit is unchanged, so unexpected
`src/`, unrelated `testing/`, and other non-allowlisted paths still fail.

Local validation passed:

- `python -m py_compile testing/collect_pion_hgcer_validation_bundle.py
  testing/test_collect_pion_hgcer_validation_bundle.py`;
- collector suite: 21 tests; and
- frozen F.1 contract/runtime suite: 10 tests.

The user chose a collector-only clean revalidation of the existing farm outputs;
no physics analysis rerun occurred or was required.

## Accepted closure

The owner accepted F.1.Fix.5 as `CLOSED / RUNTIME VALIDATED` after the clean C1
collector gate. The runtime-evaluated KaonLT analysis remains
`126fa22c19bd29b9952f55b33ab43d59f9727ef6`; C1 is only the later validation
reconciliation commit. The accepted review required the fresh collector to
confirm the same source-artifact inventory and SHA-256 continuity with this
supplied runtime bundle, including the five F.1 JSON artifacts, five procedure
PDFs, and five page-manifest JSONs. The C1 archive itself was not copied into
this checkout; this record preserves the owner-accepted closure decision rather
than inventing its path or digest.

## Next

The docs/memory-only closure commit is the F.2 starting baseline. F.2 remains
detached and must not construct a map, correction, event probability, weight,
or production application.

```tcsh
cd /group/c-kaonlt/USERS/trottar/lt_analysis
python3 testing/collect_pion_hgcer_validation_bundle.py \
  --outdir /lustre24/expphy/volatile/hallc/c-kaonlt/trottar/OUTPUT/Analysis/KaonLT \
  --kinematic Q4p4W2p74 \
  --output /lustre24/expphy/volatile/hallc/c-kaonlt/trottar/OUTPUT/Analysis/KaonLT/KaonLT_F1Fix5_validation_Q4p4W2p74_full_v4_clean-collector.zip
```

Use a new output name; this must not overwrite the supplied archive.
