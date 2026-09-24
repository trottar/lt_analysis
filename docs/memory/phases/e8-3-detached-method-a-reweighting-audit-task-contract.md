# KaonLT E.8.3 — Detached Method-A Reweighting Audit

## Starting state

- Repository: `trottar/lt_analysis`
- Branch: `test`
- Exact authoritative starting remote HEAD: `91bb7809d27d84d7709a6600dbb0dc9ab514a458`
- E.8.2: `SOURCE REVIEWED`
- E.8.3: `NEXT`
- F.6.3: `BLOCKED` pending source-reviewed E.8.3
- F.6.2/Fix.5: `CLOSED / RUNTIME VALIDATED`

This contract authorizes **presentation-only E.8.3 development**. It does not authorize Method-A production promotion, a parallel production branch, Method-B numerical use, empirical residual fits, or any change to the baseline production analysis.

The current audit establishes that the accepted persisted F.4/F.5/F.6.1/F.6.2 artifacts are sufficient for E.8.3. **No new scientific persisted-source exposure phase is required.**

## 1. Mandatory startup and worktree preservation

Before editing:

1. Establish the actual local branch, HEAD, and dirty-worktree state.
2. Confirm branch `test` and that the local base contains starting HEAD `91bb7809d27d84d7709a6600dbb0dc9ab514a458`.
3. Read, in order: root `AGENTS.md`, `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, `docs/memory/USER.md`.
4. Read the directly relevant records: `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`, `docs/memory/phases/phase-f4-method-a-parent-preserving-correction.md`, `docs/memory/phases/phase-f5-method-a-tphi-propagation.md`, `docs/memory/phases/phase-f6-method-a-production-promotion.md`, `docs/memory/evidence/f6-1-runtime-closure.md`, `docs/memory/evidence/f6-2-scientific-runtime-closure.md`, `docs/memory/evidence/f6-2-fix5-presentation-runtime-closure.md`.
5. Inspect the live runtime/presentation path in `src/cuts/full_background_subtraction_plots.py`, `src/cuts/rand_sub.py`, `src/binning/calculate_yield.py`, and `src/main.py`.
6. Inspect the current F.4/F.5/F.6.1/F.6.2 producers/readers before writing new validation code.

Do not reset, stash, clean, or otherwise destroy accepted local work. If local HEAD has advanced beyond the stated starting HEAD, hard stop and report the exact new HEAD plus intervening diff before implementation.

## 2. Scientific ownership

E.8.3 is a **detached, presentation-only Method-A reweighting audit**.

The accepted baseline pion contribution remains:

```text
b_j^0 = s_j * w0_j
```

The accepted detached Method-A variation remains:

```text
b_j^A = s_j * w0_j * C_j
```

with `C_j` exactly the accepted F.4 parent-preserving correction. Parent-only normalization is frozen; no canonical `(t,phi)` child is independently renormalized.

Method B remains diagnostic/cross-check only and contributes no numerical input. Legacy empirical residual Fit 1/Fit 2 remain dormant and absent from E.8.3.

## 3. Frozen accepted input authority

E.8.3 must fail closed unless persisted detached inputs match the accepted authority.

### F.4

Accepted F.4 JSON SHA-256:

```text
adcc01900b8adc211e296aaedaa314fdb86207d51483ebfd4b3edf69f558f188
```

Validate the current accepted schema, detached/non-production flags, Q4p4W2p74 identity, 15 canonical parents, and persisted correction fingerprint/continuity. Do not persist or reconstruct event-level `C_j` factors.

### F.5

Accepted F.5/F.5.2 JSON SHA-256:

```text
143e3af6b1c69560e5bf351155b570d05f19e7e0ceb14665f5448c351ad501be
```

Frozen F.5 scientific fingerprint:

```text
d11b728d1089301a12c29e7f8b1798c6e0b6021ac47bd5d2afc2b62b47a1effa
```

F.5 owns the accepted canonical `3 x 9` baseline and Method-A-adjusted signed `(t,phi)` aggregate matrices and parent-closure summaries.

### F.6.1

Accepted F.6.1 JSON SHA-256:

```text
62bad2dae0f65f1fff87eeffd071c29dbbd4cdc15a43f37d6d9853cb58b25bd6
```

Accepted artifact fingerprint:

```text
377872a218a780481347402e4410c81bd568a4fb7682cd49f0f3352b499bad41
```

Accepted validation fingerprint:

```text
d992d789b3434897d76691df27c190a0f51479f67c0d5b496e84835e517c77f8
```

F.6.1 owns the accepted detached parent-level signed `analysis_MM` baseline, Method-A, and difference shapes.

### F.6.2

Accepted F.6.2 JSON SHA-256:

```text
5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1
```

Accepted artifact fingerprint:

```text
ee713b70de898bad8fa61164cbc8a54886af1ea1eb1e712ed95de17df8890cd0
```

Accepted validation fingerprint:

```text
7edc73fce20ad7dc7622b8c23a7ba7e8986e367ace605884e010945595370b3b
```

Reuse the existing fail-closed E.8/F.6.2 reader. Do not create a second F.6.2 calculation.

## 4. Why direct E.8.3 implementation is authorized

The current persisted artifacts already provide the required presentation quantities:

- F.4: accepted parent-preserving correction authority.
- F.5: persisted baseline and Method-A canonical `3 x 9` signed `(t,phi)` aggregates plus closure.
- F.6.1: persisted parent-level signed missing-mass baseline/Method-A/difference shapes.
- F.6.2: persisted acceptance/HGCer explanatory diagnostics already consumed by the procedure-PDF path.

Therefore E.8.3 is readers + validation + presentation only. Do not introduce a scientific producer, checkpoint, event-factor table, or production branch.

## 5. Allowed files

Source changes are restricted to:

```text
src/cuts/full_background_subtraction_plots.py
src/cuts/rand_sub.py
```

Tests may change only:

```text
testing/test_full_background_subtraction_plots.py
testing/test_pion_hgcer_phase_e_runtime_contract.py
```

A new narrowly scoped test file is allowed if useful:

```text
testing/test_e8_3_detached_method_a_reweighting_audit.py
```

Durable memory/history may change only as warranted:

```text
docs/memory/CURRENT.md
docs/memory/phases/e8-3-detached-method-a-reweighting-audit.md
docs/memory/phases/e8-3-detached-method-a-reweighting-audit-task-contract.md
docs/memory/roadmap/STATUS.md
docs/memory/manifest.json
```

## 6. Frozen files and interfaces

Frozen by this contract:

```text
src/binning/calculate_yield.py
src/main.py
src/cuts/pion_component_subtraction.py
src/cuts/pion_component_fits.py
src/cuts/particle_subtraction.py
src/cuts/pion_hgcer_method_a_parent_preserving_correction.py
src/cuts/pion_hgcer_method_a_tphi_propagation.py
src/cuts/pion_hgcer_method_a_reweighting_validation.py
src/cuts/pion_hgcer_method_a_acceptance_refinement_validation.py
testing/analyze_pion_hgcer_method_a_parent_preserving_correction.py
testing/analyze_pion_hgcer_method_a_tphi_propagation.py
testing/analyze_pion_hgcer_method_a_reweighting_validation.py
testing/analyze_pion_hgcer_method_a_acceptance_refinement_validation.py
```

All production physics, public yields, pion weights, proton factors, random/dummy normalization, pruning, cuts, templates, fit windows, priors, binning, SIMC, efficiencies, acceptance, L/T separation, and cross-section behavior are frozen.

If implementation appears to require changing any frozen file, **HARD STOP**.

## 7. Required E.8.3 architecture

Add an explicit presentation schema, following current module conventions, e.g. `full_background_subtraction_e8_3/v1`.

The payload must explicitly carry availability/reason, current setting, input paths, observed SHA-256 values, accepted fingerprints, canonical parent/t/phi identity, no-production flags, and accepted persisted display quantities only.

Add narrow fail-closed readers in `full_background_subtraction_plots.py` for accepted F.4, F.5, and F.6.1 JSON, following the existing E.8/F.6.2 reader pattern. Each reader must verify raw-byte SHA, JSON/schema, fingerprints, detached/non-production flags, Q4p4W2p74 identity, canonical settings/parents/geometry, finite persisted values, and required continuity fields. No same-schema fallback and no recomputation.

Select the current setting only through exact canonical identity:

```text
Left-lowe
Left-highe
Center-lowe
Center-highe
Right-highe
```

No aliases or cross-setting pooling.

## 8. Required presentation

E.8.3 appends a compact, physics-readable detached Method-A section to the existing full-background-subtraction procedure PDF while preserving all prior E.8/E.8.1/E.8.2 pages.

### E.8.3a — authority and reweighting operation

Render one authority page stating:

```text
b_j^0 = s_j * w0_j
b_j^A = s_j * w0_j * C_j
```

and that `C_j` is accepted F.4, parent normalization preserves the signed canonical-t parent sum, no child is independently normalized, E.8.3 is detached/non-production, Method B is absent numerically, empirical residual Fit 1/Fit 2 are inactive, and F.6.3 alone may later construct a parallel production branch. Show accepted SHA/fingerprint provenance.

### E.8.3b — baseline versus reweighted pion missing mass

For each of the three canonical t parents in the current setting, consume persisted F.6.1 signed `analysis_MM` arrays and render:

- baseline `B_pi^0(MM)`;
- Method-A `B_pi^A(MM)`;
- direct overlay;
- signed difference `Delta B_pi = B_pi^A - B_pi^0`;
- ratio `R_pi = B_pi^A / B_pi^0` only where the persisted baseline denominator is finite and nonzero.

Do not clip, cap, smooth, interpolate, renormalize, or fabricate undefined ratio bins. Use persisted F.6.1 edges exactly.

The renderer may validate that persisted delta equals persisted Method-A minus persisted baseline, but must not silently replace a persisted delta with a newly calculated scientific result.

### E.8.3c — canonical `(t,phi)` redistribution and closure

For each canonical t parent, consume persisted F.5 `3 x 9` aggregate matrices and show all nine physical phi children:

- `B_pi,tphi^0`
- `B_pi,tphi^A`
- `Delta B_pi,tphi`

Explicit zero/empty children remain visible. Display the persisted accepted parent closure `sum_phi B_pi,tphi^A ~= sum_phi B_pi,tphi^0`, including the persisted absolute residual and persisted relative residual if the accepted F.5 schema contains one. Do not invent a new scientific closure metric if it is not persisted.

Do not reconstruct child aggregates from events and do not independently normalize a child.

### E.8.3d — accepted F.6.2 explanation

Do not duplicate or recompute the accepted E.8.1/F.6.2 explanatory pages. Preserve those existing pages and add a concise textual cross-reference explaining that they provide the accepted L/B/A shapes, acceptance maps/correlations, support/OOD, effective statistics, and kaon-window refinement context for the observed E.8.3 redistribution. No new F.6.2 science.

## 9. Runtime integration

`rand_sub.py` may change only to locate and load the accepted persisted F.4/F.5/F.6.1 artifacts from established `OUTPATH`, build the E.8.3 presentation payload, and pass it through the existing final procedure-PDF rendering path.

Use deterministic filenames already owned by accepted F-stage modules. No configurable fallback path.

Preserve the current E.8.2 producer/runtime state and the existing PDF/page-manifest pair-safe finalization. E.8.3 must not create an independent PDF lifecycle.

If a required accepted artifact is missing, stale, malformed, has the wrong SHA/fingerprint, or fails authority/geometry validation, render an explicit E.8.3 unavailable section. Do not substitute a recomputed result and do not invalidate the already completed baseline production calculation.

## 10. Required tests

At minimum cover:

1. exact accepted F.4 SHA accepted; wrong SHA rejected;
2. exact accepted F.5 SHA/fingerprint accepted; wrong authority rejected;
3. exact accepted F.6.1 SHA/artifact/validation fingerprints accepted; wrong authority rejected;
4. nonfinite JSON rejected;
5. wrong kinematic rejected;
6. missing/duplicate canonical setting or parent rejected;
7. malformed F.5 `3 x 9` geometry rejected;
8. malformed F.6.1 MM edges/content lengths rejected;
9. explicit empty F.5 children preserved;
10. current-setting selector returns only the requested setting;
11. Method-B numerical fields are neither required nor consumed;
12. empirical residual Fit 1/Fit 2 are neither required nor consumed;
13. persisted F.6.1 baseline/Method-A/delta arrays pass through without normalization/reconstruction;
14. persisted F.5 child aggregates pass through without child renormalization;
15. ratio mask is defined only where baseline is finite and nonzero;
16. no clipping/capping/smoothing/interpolation occurs;
17. persisted F.5 parent closure is displayed;
18. unavailable inputs render explicit unavailable E.8.3 pages rather than raising through production;
19. prior E.8/E.8.1/E.8.2 page ordering is unchanged except intentional appended E.8.3 pages;
20. baseline production objects/yields remain untouched;
21. no Method-B builder, empirical-residual fit builder, Method-A production application, or yield recalculation is callable from the E.8.3 path.

## 11. Local validation

Run at minimum:

```bash
python -m py_compile \
  src/cuts/full_background_subtraction_plots.py \
  src/cuts/rand_sub.py

python testing/test_full_background_subtraction_plots.py
python testing/test_pion_hgcer_phase_e_runtime_contract.py
```

If added:

```bash
python testing/test_e8_3_detached_method_a_reweighting_audit.py
```

Also run relevant consumed-artifact regression suites:

```bash
python testing/test_pion_hgcer_method_a_parent_preserving_correction.py
python testing/test_pion_hgcer_method_a_tphi_propagation.py
python testing/test_pion_hgcer_method_a_reweighting_validation.py
python testing/test_pion_hgcer_method_a_acceptance_refinement_validation.py
```

Run all repository memory health/manifest checks required by repository memory. Local checks do not establish ROOT/PyROOT, full `main.py`, farm, or runtime validation.

## 12. Forbidden shortcuts

Do not construct a production Method-A template; mutate `w0`; apply `C_j` inside production pion subtraction; construct `K_A`, `MM_A`, or `Y_A`; change `K_0`, `MM_0`, or `Y_0`; derive or persist event-level correction factors; renormalize a child; recompute F.4/F.5/F.6.1/F.6.2; use Method B numerically; activate empirical residual Fit 1/Fit 2; alter cuts/windows/templates/priors/normalization/binning/SIMC/efficiency/acceptance/LT/cross sections; invent a Method-A uncertainty; treat Method-A minus baseline as an uncertainty; soften authority failures into fallbacks; or commit/push/run the farm.

## 13. Durable memory

Keep E.8.3 `ACTIVE` until independent ChatGPT actual-diff review.

Create `docs/memory/phases/e8-3-detached-method-a-reweighting-audit.md` recording exact start HEAD, presentation-only ownership, accepted input SHAs/fingerprints, changed files, runtime path, Codex-run tests, farm boundary, no promotion, no Method-B numerical dependency, no empirical-residual dependency, no F.6.3 work, and exact NEXT = independent ChatGPT actual-diff review.

Update `CURRENT.md` only enough to make E.8.3 the active local work item while preserving:

```text
F.6.3 BLOCKED pending source-reviewed E.8.3
E.8.4 BLOCKED pending F.6.3
final E.8 BLOCKED pending E.8.4/runtime visual gate
F.6.4 BLOCKED pending full production-impact evidence
```

Do not claim E.8.3 `SOURCE REVIEWED` yourself. Regenerate `docs/memory/manifest.json` whenever versioned memory changes.

## 14. Diff audit and review bundle

Before stopping:

1. print exact branch and HEAD;
2. print `git status --short`;
3. print exact changed-path list;
4. verify every changed path is allowlisted;
5. run `git diff --check`;
6. inspect all changed source/tests;
7. confirm frozen files are byte-unchanged relative to starting HEAD;
8. confirm no production-physics file outside the two allowed presentation-integration files changed;
9. confirm no Method-B or empirical-residual numerical dependency was added.

Generate root-level `kaonlt_review.diff` containing the complete tracked diff from starting HEAD plus complete `git diff --no-index /dev/null ...` sections for every intended untracked/new file, including this contract and the E.8.3 phase record. `kaonlt_review.diff` is temporary and must not be committed.

## 15. Acceptance criteria

Implementation is complete only if:

1. E.8.3 consumes accepted persisted F.4/F.5/F.6.1/F.6.2 artifacts only.
2. SHA/fingerprint authority is fail-closed.
3. No scientific correction/event-factor calculation is created.
4. `b_j^0 -> b_j^A = b_j^0*C_j` is shown explicitly.
5. Parent-level baseline/Method-A MM uses persisted F.6.1 arrays.
6. Signed MM difference is direct and physics-readable.
7. Ratio is shown only where mathematically defined.
8. Canonical `3 x 9` redistribution uses persisted F.5 aggregates.
9. All nine phi children remain visible for each t parent.
10. No child normalization occurs.
11. Persisted parent closure is displayed.
12. Existing F.6.2 explanation pages remain authoritative.
13. Method B is numerically absent.
14. Empirical residual Fit 1/Fit 2 are absent.
15. Baseline production behavior and E.8.2 outputs are unchanged.
16. F.6.3 production branch is not constructed.
17. Required deterministic tests pass or an exact blocker is reported.
18. Memory remains `ACTIVE` pending independent ChatGPT review.
19. Complete `kaonlt_review.diff` is generated.
20. No commit, push, or farm run is performed.

## 16. Hard stop

Stop and report a blocker if accepted F.4/F.5/F.6.1 persisted schemas do not actually contain the quantities described above; a required E.8.3 plot would require reconstructing event-level factors; production pion subtraction would need modification; any frozen file would need editing; accepted authority cannot be uniquely pinned; local HEAD materially differs from starting HEAD; or existing dirty work cannot be preserved safely. Do not expand scope to solve a hard-stop condition.
