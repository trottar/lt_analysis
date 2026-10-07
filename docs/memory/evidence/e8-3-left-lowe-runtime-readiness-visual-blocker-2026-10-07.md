# E.8.3 Left/lowe runtime-readiness visual blocker — 2026-10-07

## Authority and continuity

E.8.3 remains `SOURCE REVIEWED`. Its narrow runtime-readiness gate is
`BLOCKED` only at authority-page provenance legibility. This record transcribes
supplied artifacts and prior ChatGPT inspection specified in the
[repair contract](../phases/e8-3-authority-page-provenance-legibility-repair-task-contract.md);
Codex did not rerun, hash or independently render the external artifacts.

Current implementation-start source: `0ac6134c6f76efa731f94390b986b9533b53a7f6`.
SOURCE VERIFIED: the compare from farm source
`706ae708ae69be0585472982a4df3d87c63f186f` to current source contains exactly
seven `docs/memory/` paths. Farm E.8.3 source continuity is preserved.
The [accepted E.8.2 closure](e8-2-left-lowe-runtime-scientific-closure-2026-10-06.md)
owns package/owner/companion identities.

## Accepted fresh E.8.3 evidence

The fresh procedure PDF belongs to package stem:

```text
KaonLT_E8_2_Left_lowe_scientific_audit_20261006-233426
```

Farm source:

```text
706ae708ae69be0585472982a4df3d87c63f186f
```

Scope:

```text
Q4p4W2p74 / Left / lowe
```

Procedure PDF:

```text
Left_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction.pdf
SHA-256:
aacac5daa935eb1f2af6dacdf0ddfdcdd4f5b202e0d54fc11b8edd254576e61e
```

Page manifest:

```text
Left_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction-manifest.json
SHA-256:
e51fe03cff97ce00e05c850b031206e97a9cb9d15b0c71279618fda1a9577f02
```

The existing accepted E.8.2 closure evidence owns the complete package,
owner/provenance/artifact, and companion identities.

RUNTIME VERIFIED from the fresh supplied artifacts and ChatGPT inspection:

- the procedure PDF has 74 pages;
- page-manifest `renderer_failures` is empty;
- setting identity is exactly Left / lowe / Q4p4W2p74 / kaon;
- all eight E.8.3 page records are present:
  - page 65 `full_background.e8_3.authority`;
  - page 66 `full_background.e8_3.mm`, t1;
  - page 67 `full_background.e8_3.tphi`, t1;
  - page 68 `full_background.e8_3.mm`, t2;
  - page 69 `full_background.e8_3.tphi`, t2;
  - page 70 `full_background.e8_3.mm`, t3;
  - page 71 `full_background.e8_3.tphi`, t3;
  - page 72 `full_background.e8_3.f6_2_cross_reference`;
- pages 66–72 are visually complete and legible;
- all three MM pages show persisted baseline, Method-A, delta, and gap-safe
  defined-bin-only ratio points;
- all three `(t,phi)` pages show all nine children, persisted EMPTY/POPULATED
  identity, baseline/Method-A/delta values, and persisted parent closure;
- the F.6.2 cross-reference page is legible;
- no E.8.3 renderer failure is present.

This evidence establishes that the existing E.8.3 runtime route is active in
the real farm procedure path. It does **not** close E.8.3 because page 65 fails
visual provenance legibility.

---

## 4. First failed invariant: authority-page SHA-256 clipping

The sole observed E.8.3 visual blocker is page 65:

```text
full_background.e8_3.authority
```

The page is otherwise legible, including:

- historical accepted F.6.1 lineage label;
- detached/non-production statement;
- Method-B exclusion;
- dormant empirical residual statement;
- F.6.3/E.8.4/F.6.4 ownership text;
- setting identity;
- the three fingerprint lines.

But these two current lines place two 64-character SHA-256 values on one line:

```text
Accepted input SHA-256: F.4=<64 hex> F.5=<64 hex>
F.6.1=<64 hex> F.6.2=<64 hex>
```

The physical A4 farm PDF clips the right-hand hashes:

- the F.5 SHA-256 is visibly truncated;
- the F.6.2 SHA-256 is visibly truncated.

PDF text extraction confirms the same right-edge loss.

SOURCE VERIFIED cause in:

```text
src/cuts/full_background_subtraction_plots.py
```

`_e8_3_render_authority_page()` passes those long strings directly into the
non-wrapping `_e8_text_page()` / `_e8_add_text()` path.

This is a presentation-only defect. The persisted authority values themselves
are correct and already pass source/runtime authority checks.

Therefore the current E.8.3 readiness state is:

```text
E.8.3:
SOURCE REVIEWED

E.8.3 Left/lowe runtime-readiness gate:
BLOCKED at rendered authority-page provenance legibility
```

Do not reinterpret or recompute the E.8.3 scientific arrays.

## Retained boundary

The [local repair](../phases/e8-3-authority-page-provenance-legibility-repair.md)
changes displayed provenance layout only. Fresh physical PDF acceptance remains
NOT VERIFIED. No E.8.3 runtime closure, new scientific magnitude conclusion,
Method-A promotion or production correctness follows. Historical authorities,
scientific arrays and production remain frozen. Method A stays
detached/non-production; Method B diagnostic/cross-check only and numerically
excluded. Canonical-five provenance repair remains `DEFERRED`; final E.8
and F.6.4 remain `BLOCKED`.
