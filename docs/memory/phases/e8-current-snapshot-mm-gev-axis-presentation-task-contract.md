# KaonLT E.8 — Current-snapshot missing-mass axis in GeV (presentation-only)

## Authority, starting state, and objective

**Task class:** one narrowly scoped presentation-source change. Current published `test` HEAD observed for this audit: `75756d7b32767f18949386c4c5864c756d5e85f6`. Re-establish actual remote/local branch, HEAD, index and complete worktree; do not rely on this stored observation as live state.

The user requests that the **accepted detached current-snapshot contribution diagnostic PDF** display **actual analysis missing mass in GeV**, not missing-mass bin indices. This applies to page 5 (three missing-mass marginal plots) and pages 6–8 (all joint phi-versus-MM heatmaps). It is a plotting-coordinate/axis-only improvement, **not** a new scientific measurement, rebinning, normalization, correction, authority change, yield calculation, or production promotion.

Read root `AGENTS.md` and the four mandatory startup files in exact order, then `docs/memory/MAINTENANCE.md`, `CODEX.md`, `TOOLS.md`, `COMMUNICATION.md`, the current-snapshot diagnostic task contract and this contract. Audit live renderer/owner, frozen aggregate/geometry validator, and tests before editing. The previous accepted October 12 detached ZIP and PDF remain `CLOSED / RUNTIME VALIDATED` **only as already rendered**. New renderer pages are pending validation until the user supplies a fresh farm-rendered ZIP/PDF. Canonical-five F.4 equivalence and absolute-SIMC remain `BLOCKED`.

## Explicit pending-unpublished checkpoint constraint

**Do not push the previously accepted four-file evidence checkpoint yet.** The user's workstation holds an independently accepted but **uncommitted** checkpoint diff: `kaonlt_review_e8_current_snapshot_checkpoint_2026_10_12.diff` (previous uploaded revision 38,409 bytes, SHA-256 `c2d1cb53914aaa600471eea3db686c7bdd976956bc501387591b7b4c68f535f5`). The pending changes are `docs/memory/CURRENT.md`, `docs/memory/manifest.json`, new `docs/memory/evidence/e8-current-snapshot-contribution-runtime-closure-2026-10-12.md` and new `docs/memory/phases/e8-current-snapshot-contribution-runtime-evidence-checkpoint-task-contract.md`. The original review-file name/path is temporary, and is **not** scientific authority.

Preserve both previously accepted new Markdown records **byte-for-byte**, the corrected t1/t0/t2 statement, scientific NEXT intent, historical closures, unrelated dirty files and index state. Confirm the pending four-file candidate is present and materially matches accepted content; **stop** rather than overwriting/reverting/recreating it if missing, changed or conflicting. Subsequent actual-diff review must include the complete **combined** uncommitted review bundle (not only the incremental renderer diff), including all new files. No staging or Git publication by Codex.

## Exact additional path allowlist

May edit:

1. `testing/analyze_pion_hgcer_method_a_current_snapshot_contribution_debug.py` — `write_review_pdf` and minimal pure presentation helpers in this same tracked owner only.
2. `testing/test_pion_hgcer_method_a_current_snapshot_contribution_debug.py` — deterministic tests of the persisted-JSON renderer and MM-axis geometry.
3. `docs/memory/CURRENT.md` — only a brief status note to distinguish the already accepted original diagnostic from the GeV-axis re-render pending farm visual inspection, and make the **one NEXT** the narrow visual rerender/verification gate, with subsequent baseline/background science retained as successor.
4. `docs/memory/manifest.json` — generated integrity metadata following memory edits.

Create **only**:

5. `docs/memory/phases/e8-current-snapshot-mm-gev-axis-presentation-task-contract.md` — this task contract.

The two previously accepted, pending checkpoint Markdown additions are allowed as **unchanged carry-forward** files in the complete review bundle; they are not editable by this task. No change to numerical aggregate builders, diagnostic JSON schemas/fingerprints, `src/`, frozen F.1–F.4 builders, Method B, private F.6.3, production, historical renderer, yield/MM extraction, existing evidence, or user environment records. Do not create an unrelated standalone PDF/script as a substitute for the tracked renderer.

## Required display semantics

Existing persisted `diagnostic.json` contains `display_geometry.mm_edges`: **41** physical bin edges defining **40 normal bins** in `[0.800, 1.400] GeV` with width `0.015 GeV`. Each displayed parent has `mm`: **42 entries**, and `joint`: **9 phi rows × 42 MM entries**. The first and last entries are categorical **underflow (UF)** and **overflow (OF)**, **not** physical MM measurements. Use only the validated stored edges; do not reconstruct event-level data or alter their contents.

- On **page 5**, plot the 40 normal marginal cells at their **true physical bin centers in GeV** (first center `0.8075`, last `1.3925`) with readable x-axis labels and suitable ticks near `0.8, 0.9, ..., 1.4`. Preserve baseline, reweighted, signed-delta curves exactly. Render UF/OF values distinctly (isolated markers or separated display strips), explicitly labeled `UF` and `OF`; never connect them as if they were measurements at real masses.
- On **pages 6–8**, retain all 42 joint columns. Map the 40 normal columns to the **physical MM edges** and use only clearly separated **categorical display strips** for UF/OF outside the normal range, labeled `UF`/`OF`. Matplotlib `imshow(..., extent=...)` (or equivalent pixel-preserving method) may put side strips of one normal-bin display width on either side; these are a **visual accommodation**, not claims of physical MM coordinates. Clearly mark the true `0.8`/`1.4` physical boundaries (e.g., separator lines); tick labels for normal columns should be in GeV. Preserve original matrix ordering, raw values, 9 phi rows, color limits/semantics and colorbars.
- All MM x-axis labels should indicate `Analysis missing mass [GeV]` and explain that `UF`/`OF` are outside-range categories. Do **not** use `MM bin`, `0..41`, or index numbers as the primary x-axis. Keep axes legible, avoid crowded ticks, maintain all eight existing pages, exact page-title inventory and descriptive/non-authoritative disclaimers.
- The change is **render-only**. No numerical recomputation, bin-width multiplication, edge drift, interpolation, clipping, flow elimination, altered tails/Neff, histogram mutation, diagnostic JSON rewrite, changed fingerprint/ZIP manifest schemas, or production effect.

## Deterministic local gates

Add a focused test for physical centers/edges/ticks using the source-defined/persisted `mm_edges` (including the first and last centers), explicit separated `UF`/`OF` semantics, and unchanged 42-value/9×42 visualization coverage. Validate the complete persisted-only eight-page renderer under `Agg`, existing title/page count, and no transient-F.4/data reload. The original diagnostic and snapshot authority tests must still pass. Inspect rendered page 5 and pages 6–8 with actual PDF extraction/rendering where feasible; **local** rendering is a regression check, not fresh Jefferson Lab farm acceptance.

Run at minimum:

```text
python -B -m unittest testing.test_pion_hgcer_method_a_current_snapshot_contribution_debug testing.test_pion_hgcer_method_a_reweighted_contribution_debug testing.test_materialize_method_a_current_f1_contribution_snapshot
```

Use the locally discovered available interpreter, not a hardcoded version; report if dependencies cannot run. Check Python syntax, `git diff --check`, allowed paths and preservation, memory manifest **write then check**, bootstrap/ordinary memory health. A CURRENT 8 KiB **soft** warning is nonblocking if no hard failures or material drift; do not broaden the task into unrelated memory refactoring.

## Delivery, review and next gate

Return one complete actual-diff review bundle covering the four already accepted, pending checkpoint paths plus this change, including **complete `--no-index` additions** for each of the **three** untracked Markdown documents (two previously accepted and this new contract). Include a concise changeset report, deterministic test receipt, renderer screenshots/page-text evidence if available, preservation status, exact remaining warning(s) and full KaonLT health receipt. **Do not stage, commit, push or run the farm.** ChatGPT independently reviews the actual diff; only the user can publish after acceptance.

After push and pushed-state synchronization, **one** user-controlled, source-owned full detached current-snapshot diagnostic rerun is permitted only after verifying the same seven-input snapshot and exact preflight fingerprint, and existing external output roots. Review ZIP manifest/JSON identity and new PDF **page 5 and pages 6–8** before granting *visual* acceptance. Do not reinterpret numerical science or reopen the canonical-owner blocker. Then return to the read-only baseline missing-mass/yield versus detached-background scientific crosswalk (with lineage/normalization comparability checked before quantitative comparison).

## Hard stops

Stop on an unexpected source/branch/worktree conflict, changed prior accepted checkpoint, unknown MM geometry or flow layout, inability to map normal cells consistently from persisted edges, new scientific content or PDF role, changed numerical arrays/schemas/fingerprints, unreviewed executable steps, or failed manifest/health/renderer checks. Diagnose the first invariant. No automatic farm command, canonical-owner retry, Method-A promotion or new acceptance.
