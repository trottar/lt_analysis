# E.8 F.6.2 figure-library implementation contract

## Status and supersession

This remains the historical narrow frozen-F.6.2 figure-library implementation
contract and preserves the completed F.6.2 presentation context. It does not
own the forward E.8 program. Future E.8 architectural, presentation, and
dependency ownership is superseded by the
[full-analysis procedure roadmap](e8-full-analysis-procedure-roadmap.md).

## Objective and exact starting state

Implement E.8 as a separate, deterministic, presentation-ready PDF figure
library over the one accepted frozen F.6.2 JSON. This contract was prepared
after the F.6.2 closure reconciliation at observed clean `test` HEAD
`73bb8c2f9c300890ada23598aeacc8e8df859b26`, descended from the reconciliation
baseline `b789d203e11f0927deb59ebcae9dc59fe8add4ae`.

E.8 is `ACTIVE`. It is an implementation task, not a scientific or farm
closure. F.6.3 remains `BLOCKED` pending accepted E.8 evidence.

This narrow contract repair starts from observed clean `test` HEAD
`0a93aaa9fd1d4ba58ec02e03d24e7c59f061044d`. After the user reviews, commits,
and pushes this repair, S1 must start from that newly pushed repair commit only
after a fresh live-HEAD check. Neither the historical `73bb8c2...` preparation
commit nor this repair-start observation is a permanent S1 source identity; the
later user-created S1 commit is the E.8 renderer-source identity.

## Frozen authority and scientific boundary

The only physics input is the accepted F.6.2 artifact JSON with SHA-256
`5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1`, artifact
fingerprint `ee713b70de898bad8fa61164cbc8a54886af1ea1eb1e712ed95de17df8890cd0`,
and validation fingerprint
`7edc73fce20ad7dc7622b8c23a7ba7e8986e367ace605884e010945595370b3b`.
The scientific source remains
`0b37af2a2927b08bdeaf897c545f290b55329cea`.

Fix.5 renderer source `c88ed65cb18ba6a37358897292b77696016312d1` and
bundle/profile commit `b789d203e11f0927deb59ebcae9dc59fe8add4ae` are upstream
presentation-artifact provenance only. They must not be called E.8 source or
profile identities. E.8 records its own renderer-source and bundle/profile
commits after those commits exist.

The implementation must never rebuild or alter Method-A factors, bootstrap
results, support/OOD, cuts, weights, binning, pion fits, random/dummy or
slow-proton treatment, Method B, yields, cross sections, or scientific
acceptance. It must not import or modify `rand_sub.py`,
`full_background_subtraction_plots.py`, the F.6.2 artifact builder, or the
accepted F.6.2 review renderer.

## Two-commit implementation and profile sequence

### S1 — renderer source

At this contract baseline, add only:

- `testing/render_pion_hgcer_method_a_acceptance_refinement_figure_library.py`
- `testing/test_render_pion_hgcer_method_a_acceptance_refinement_figure_library.py`

The renderer is a standalone Python CLI using `matplotlib` Agg and existing
standard-library facilities only. It has no ROOT dependency and produces no
analysis artifact. After local source review, the user commits and pushes S1.
The resulting S1 commit is the E.8 renderer-source commit.

### S2 — evidence profile

Only after the user supplies the pushed S1 commit, add:

- `testing/pion_hgcer_validation_bundle_profile_e8.json`
- `testing/test_pion_hgcer_validation_bundle_profile_e8.py`

Use the existing generic collector unchanged. The new profile must pin S1 as
`required_analysis_commit`, allow only the collector/profile/rendering files
and `docs/memory/` after that source, and package exactly the frozen F.6.2 JSON,
the E.8 PDF, and the E.8 page-manifest JSON. After S2 is committed/pushed, its
commit is the E.8 bundle/profile commit.

This present contract-writing change is documentation only. S1 and S2 are
separate follow-up implementation steps authorized only within the exact file,
validation, and authority boundaries stated here.

## Renderer CLI and fail-closed input contract

`main()` accepts exactly these required command-line paths/identities:

```text
--input-json PATH
--output-pdf PATH
--output-manifest PATH
--renderer-source-commit SHA1
--bundle-profile-commit SHA1
```

All paths are explicit, resolved, distinct, and the output paths must not
exist. The two commit arguments are lowercase 40-hex values. The renderer must:

1. Read the input bytes once, hash those same bytes, and require the accepted
   JSON SHA-256 before parsing them.
2. Require the F.6.2 artifact schema, `validation.available=true`, all frozen
   fingerprints above, the aggregate-only/non-production flags, 15 parents,
   135 canonical children, and exactly 115 populated plus 20 explicitly empty
   children. Any structural or authority mismatch is a nonzero failure.
3. Render from persisted JSON values only. An explicitly persisted unavailable
   L/B/A shape or metric remains a visible literal unavailable reason; it is
   not synthesized, imputed, or converted to a new decision.
4. Write PDF and manifest through private temporary files, atomically promote
   both only on success, and remove every temporary or partially promoted output
   on failure. The input file is never opened for writing.

The deterministic output basenames are
`Q4p4W2p74_kaon_pion-background_hgcer_method-a-acceptance-refinement-figure-library.pdf`
and
`Q4p4W2p74_kaon_pion-background_hgcer_method-a-acceptance-refinement-figure-library-manifest.json`.
The resolved output paths may use caller-provided directories, but their
basenames must be exactly these values; a wrong basename is a nonzero failure
before rendering. The PDF, manifest, and input paths remain distinct, and the
two output paths must not pre-exist.

## Figure-library and page-manifest contract

Render one provenance cover and one polished physics page for every populated
canonical child: exactly 116 pages from the accepted artifact. Traverse the
stored parent order, then each stored child order; do not select, rank, or omit
a populated child. Empty children are listed as omitted in the manifest only.

Each child page has:

- Row 1: L/B/A normalized overlays for `analysis_MM`, `SHMS_xptar`, and
  `SHMS_yptar`; every overlay y-axis is labeled **Normalized bin fraction**,
  and the missing-mass panel visibly marks the frozen kaon window.
- Row 2: L/B/A `SHMS_delta × SHMS_xptar` localization maps read directly from
  the persisted `SHMS_delta__SHMS_xptar` joint-distribution payloads.
- Row 3: L/B/A `SHMS_delta × SHMS_yptar` localization maps read directly from
  the persisted `SHMS_delta__SHMS_yptar` joint-distribution payloads.
- A compact persisted-metric footer limited to population counts, effective
  sample sizes, prompt/full OOD fractions, one-dimensional `analysis_MM`
  `DeltaH` and `kappa`, persisted `analysis_MM__SHMS_xptar` and
  `analysis_MM__SHMS_yptar` `kappa`, kaon-window `P_B^K`, `P_A^K`,
  `DeltaP^K`, `f_refine^K`, and `R_V`, including only their persisted
  values/intervals or literal unavailable states.

The displayed localization maps are `SHMS_delta × acceptance`, whereas the
persisted acceptance-refinement diagnostics in the footer remain
`analysis_MM × acceptance`; neither replaces or is collapsed into the other.
For each 2D row, the physical x/y axes remain their named variables and the
color scale/colorbar is labeled **Normalized bin fraction**. Each row shares
one color normalization across its available L/B/A matrices, determined from
the finite persisted values of those matrices; panels must never autoscale
independently. Unavailable populations remain visibly unavailable with their
literal persisted reason. This does not renormalize the persisted `unit_area`
values or introduce a scientific transformation.

L is the upstream low-HGCer observed reference from accepted upstream
experimental events before the downstream gate (`0 < NPE <= 2`); B is the
physical pion-control population (`NPE > 2`) with `w0`; A is that same physical
control population with `w0*C`. L is neither production kaon nor absolute
leakage truth. The page caption must state that B-to-A comparison is descriptive
acceptance refinement, not a uniform-improvement, yield, or promotion claim.

The manifest schema is
`pion_hgcer_f6_2_e8_figure_library_manifest/v1`. It records the accepted input
SHA/fingerprints and scientific source, Fix.5 upstream presentation provenance,
the CLI-supplied E.8 renderer-source and bundle/profile commits, output
basenames, page count, 15/135/115/20 inventory, and ordered page records. Each
page record identifies cover or child page ID, setting, canonical-t index, phi
index/range, and whether it was rendered or explicitly empty. It is descriptive
provenance, never a scientific acceptance decision.

## Tests, validation, and hard stops

S1 tests use temporary fixtures only and cover successful deterministic output,
116-page inventory/ordering, input immutability, manifest provenance, output
collision rollback, wrong SHA, malformed JSON/schema/fingerprint/authority,
bad geometry or matrix dimensions, invalid 15/135/115/20 inventory, and
literal persisted unavailable rendering. Static tests must also prove no ROOT,
`rand_sub`, F.6.2 builder, factor, bootstrap, support/OOD, or yield-calculation
path is imported or called.

Given byte-identical accepted F.6.2 input JSON and identical CLI provenance
identities, repeated E.8 renders must produce byte-identical PDF and
page-manifest outputs. S1 must omit any wall-clock timestamp, hostname,
absolute temporary path, process ID, random UUID, or other run-dependent
metadata; it must use fixed or omitted Matplotlib/PDF `CreationDate` and
`ModDate`, stable PDF metadata, deterministic page ordering, and deterministic
manifest serialization with stable key ordering and formatting. Focused tests
must render the same fixture twice and assert identical SHA-256 values for both
PDFs and both manifests.

S2 profile tests require the three declared global artifacts, pin the supplied
S1 commit, reject missing/invalid E.8 PDF or manifest, and reject unexpected
committed source changes. Local validation for each implementation step runs
the focused unittest module, `py_compile` for each changed Python file, the
relevant profile test, `git diff --check`, and an exact allowlist audit.

Hard-stop rather than fallback for any frozen-input hash/fingerprint mismatch,
malformed persisted payload, unknown commit identity, existing output path,
or unexpected changed file. Do not run the analyzer, regenerate F.6.2 JSON,
modify the ordinary procedure PDF, create farm artifacts, commit, push, or
claim farm/runtime validation.

## Acceptance and authority

Codex performs local edits and deterministic checks only. ChatGPT audits the
actual diff. The user alone commits/pushes S1 and S2 and runs the farm. ChatGPT
reviews the fresh E.8 ZIP/PDF/manifest before E.8 can become
`CLOSED / RUNTIME VALIDATED`; until then it remains `ACTIVE` and F.6.3 remains
`BLOCKED`.
