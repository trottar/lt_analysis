# KaonLT E.8 — paired low-epsilon scientific scope and GeV-axis acceptance memory checkpoint

## Authority and objective

Starting published `test` HEAD for ChatGPT audit: `c8e830f38cf59a22b7b6d9753b1a847fe1124e70`. Reverify actual `test` HEAD, `origin/test`, branch, index, full worktree, and relevant tracked files. Read root `AGENTS.md`, then in exact order `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, `docs/memory/USER.md`, then `docs/memory/MAINTENANCE.md`, `CODEX.md`, and the directly applicable evidence. Stop if live source identity materially differs or a scoped target has unexpected local modifications.

**Task class:** one memory-only checkpoint requested explicitly by the user, **before** the next scientific analysis. No numerical, renderer, production, or analysis-source changes.

Record the user's standing **E.8 quick-debugging preference**:

> Analyze **Left/lowe and Center/lowe together by default**, not only Left/lowe. The low-epsilon runs are adjacent and quick to analyze. Left/lowe has restricted phi support; Center/lowe covers a wider phi range. Comparisons must report actual populated phi coverage for each setting. Preserve Left-only accepted historical conclusions as setting-local, not representative of all low epsilon.

Do not reinterpret this as a claim of complete physical acceptance, equality of cross-setting normalization, a command to include other settings, or a new high-epsilon analysis requirement. Do not infer or sum event populations across settings without valid independent authority. Do not change existing cuts or sources.

## Previously accepted evidence to checkpoint

The user-supplied October 13 GeV-axis rerun was independently accepted in the preceding ChatGPT review, for **presentation-only** closure, with exact identities rechecked from the supplied ZIP:

- New ZIP: `kaonlt-current-snapshot-contribution-b7719f070b9048debd28d3c1ba666c13.zip`
- ZIP SHA-256: `610c6414a1ea9a59219516b4e2473348a702d794ecfbdbfb71ea45347aeb3e6d`
- Farm source HEAD: `c8e830f38cf59a22b7b6d9753b1a847fe1124e70`
- Seven-input snapshot fingerprint: `3644eb1b91932519a6a92f08107e77350d558ef7c1bb376cbcd15781df7085e7`
- Manifest fingerprint: `d94df1fc25f7acd255e219546897a528436bb7bce4a1d439a75a6b32103b42ab`
- Diagnostic fingerprint: `2faa10d6532297082ac312c93e0e243fc8ebe238eebbaec15027ad6b3f9b4584`
- PDF SHA-256: `753df68284f7dfc85a1f3949ca4ff281d15bfe4406d5a806ceadf9eb1bfed9b1`
- Exactly three ZIP members, CRC/member digest/manifest gates PASS; all eight PDF pages inspected; pages 5–8 show physical 0.8–1.4 GeV coordinates and distinct UF/OF categories. The diagnostic numerical objects are identical to the previously accepted `c6392716...` original except new Git provenance and consequent digest.

Attribute **RUNTIME VERIFIED** to the supplied ZIP and prior independent ChatGPT review, not Codex execution. Preserve October 12 original diagnostic acceptance and all earlier scoped closures unchanged.

## Exact allowlist

Edit only:

1. `docs/memory/CURRENT.md` — compactly correct obsolete "farm visual rerun pending" text; mark **GeV-axis presentation** `CLOSED / RUNTIME VALIDATED`; state E.8 scientific work remains `ACTIVE` and the default **paired Left/lowe + Center/lowe** quick-debugging scope. Set one push-stable scientific NEXT: **read-only paired low-epsilon phi-coverage and baseline/reweighted pion-contribution inspection using already accepted five-setting F.4 materialization, starting with each setting's t0/t1/t2 source-owned phi cells**. Separate within-setting results; inspect cross-setting overlap cautiously. If physical MM-resolution Center aggregates are missing, report the precise missing source-owned data/owner gate instead of inventing results. State microscopic F.1/F.4 forensic work stays `DEFERRED` and canonical-five owner/absolute SIMC stays `BLOCKED`. Do not send NEXT back to those deferred tasks.
2. `docs/memory/USER.md` — one durable, narrowly phrased collaboration-preference note under the existing preferences section for paired Left/Center low-epsilon diagnostics and acceptance-sensitive comparisons.
3. `docs/memory/evidence/e8-current-snapshot-contribution-runtime-closure-2026-10-12.md` — short appended GeV-axis farm visual-acceptance addendum, scoped to **Left/lowe presentation only**, with above exact identity and numerical-preservation boundary. Do not mislabel Center/lowe as newly farm-validated for MM/PDF.
4. `docs/memory/manifest.json` — generated, not hand-edited.

Add exactly one new tracked contract:
5. `docs/memory/phases/e8-paired-lowe-scientific-scope-memory-checkpoint-task-contract.md` — this file, unchanged in substance.

`MEMORY.md`, `CURRENT_HANDOFF.md`, `LEARNINGS.md`, roadmaps, other evidence, all `src/`, `testing/`, production/runtime input paths and user environment are frozen. Preserve unrelated dirty files, index and untracked review bundles. Do not stage, commit, push, or run the farm.

## Acceptance and validation

- Verify the four modified paths plus single contract are the **only** task changes and the new contract is included **in full** in the review bundle. Do not stage to capture untracked additions.
- Preserve accepted October 12 Left-only F.4 numbers and scientific interpretation; do not promote the detached correction or claim a current same-lineage final `YA/MM_A`.
- Wording must explicitly say `Left/lowe + Center/lowe` is the **default for forthcoming paired quick-debug science**, not historical coverage of earlier Left-only artifacts.
- The exact single `NEXT` must stay accurate before/after publication and should yield science rather than more memory-only chores.
- `python -B tools/update_memory_manifest.py --root . --write` then `--check`, with a discovered working interpreter; ordinary `python -B tools/check_memory_health.py --root .`, bootstrap as required, `git diff --check`, path scope, SHA/manifest integrity, and numeric-preservation inspection. Report `CURRENT`/`MEMORY`/`CURRENT_HANDOFF` byte counts and all warnings. The CURRENT 8-KiB soft warning is nonblocking unless hard/materiel ambiguity arises; make compact substitutions rather than expanding substantially, but do not perform unrelated compaction.
- Produce **one complete combined actual-diff bundle** in the repository root, including all four tracked modifications plus the entire new contract as a `git diff --no-index /dev/null` addition; report exact size/SHA.
- ChatGPT reviews actual diff. User alone performs commit/push. ChatGPT then verifies pushed source and begins the paired analysis. No JLab farm operation is authorized by this memory task.

**Hard stop:** If scoped local changes conflict with already accepted evidence, if target files are not present, if identity/manifest gates fail, or if a purported scientific comparison requires non-comparable lineages or absent source-owned data, stop and report the first failed invariant without broadening scope.
