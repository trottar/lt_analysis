# E.8 four-way final-yield and missing-mass presentation — memory consistency checkpoint task contract

## 1. Authority, starting identity, and task class

**Task class:** one narrowly scoped **memory-only clarification/checkpoint** of an explicitly stated user presentation requirement. This is not a physics implementation, new E.8 phase, proof of an independent reweighting branch, accepted runtime evidence, or farm authorization.

Authoritative repository: `trottar/lt_analysis`, branch `test`.
**ChatGPT-verified starting remote HEAD:** `1924bc1a8d8e218968a27a8327c0528d5f96eaa6`.

Before making changes, establish local branch, `HEAD`, `origin/test`, `git status --short --branch`, full tracked/untracked/index state and preservation baseline. Require local `test`, local HEAD and `origin/test` to equal the starting SHA. A dirty worktree is not automatically grounds to discard or reset anything: preserve every unrelated byte, index entry, existing review bundle and untracked path, and stop `BLOCKED` if the proposed allowlist overlaps prior unrelated changes or their provenance cannot be separated. Never pull, reset, stash, clean, stage, commit, push, or change refs as part of the Codex task.

Read exactly in tracked root-`AGENTS.md` startup order: (1) `AGENTS.md`, (2) `docs/memory/CURRENT.md`, (3) `docs/memory/MEMORY.md`, (4) `docs/memory/handoffs/CURRENT_HANDOFF.md`, (5) `docs/memory/USER.md`. Then read `docs/memory/MAINTENANCE.md`, `docs/memory/CODEX.md`, `docs/memory/TOOLS.md`, `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`, the E.8 part of `docs/memory/roadmap/STATUS.md`, and relevant E.8.4 source (`src/cuts/full_background_subtraction_plots.py`, producer `src/binning/calculate_yield.py`) **read-only**. Follow CURRENT direct references only as needed. The tracked records and live source outrank this handoff if material drift is established; stop and report it instead of rewriting science.

## 2. User-requested clarification — exact scope

On October 11, 2026, the user corrected an earlier three-way presentation summary. The **desired final comparison contains four named categories for each canonical `(t,phi)` cell, both in extracted yields and in missing-mass spectra**:

| Requested presentation category | Yield label (notation for request only) | MM label (notation for request only) |
| --- | --- | --- |
| Baseline | `Y0(t,phi)` | `MM_0(t,phi)` |
| Method A | `YA(t,phi)` | `MM_A(t,phi)` |
| **Reweighted baseline** | `Y_rw(t,phi)` | `MM_rw(t,phi)` |
| `K+Lambda` SIMC | `Y_SIMC(t,phi)` | `MM_SIMC(t,phi)` |

These are **four requested presentation slots**, not evidence that four independent producer results exist. In particular, `Y_rw` and `MM_rw` are **notation for an unresolved requested branch**, NOT newly established output keys, production formulas, files, authority objects, or scientific values. Preserve the user's separate reweighted-baseline label; do not silently merge it into Method A, but do not manufacture a distinction where the current source provides none. An eventual independent source/definition audit must determine whether the requested reweighted baseline is separately defined and supported, whether it overlaps existing `w0 -> w0*C` Method-A semantics, or requires a further **user-reviewed scientific definition**. Record that determination as `NOT VERIFIED`, not as a new NEXT or new implementation phase.

The requested eventual comparison includes (when authoritative outputs are available and scientifically comparable): four MM curves with common canonical bin geometry and Lambda window, four extracted-yield series versus phi per canonical t and summaries across t/settings; existing uncertainties without fabricated `DeltaY` errors; impact differences and defined fractions with exact reference names fixed only after identifying the reweighted-baseline source. Show unavailable/undefined states rather than substituting lookalike or independently scaling SIMC. Existing `B_pi^A` is a **reweighted pion-background component** and is not automatically a separate final `MM_rw` or `Y_rw` branch.

**SOURCE VERIFIED from precontract audit:** the current E.8.4 renderer implements baseline/Method-A/SIMC three-way MM overlays, and its yield summary uses `Y0`/`YA` and `DeltaY`; the tracked E.8 roadmap specifies the same implemented scope. There is no separately defined/consumed `Y_rw`/`MM_rw` final output in those reviewed paths. The user-requested expansion is **not** claimed to be implemented or farm validated.

## 3. Approved minimal edits / exact allowlist

1. `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md` — **primary specification home**. Add one short, clearly labeled *user-clarified final presentation requirement / unresolved source-identity boundary* near E.8.4 or final-E.8 closure. Preserve the existing **historical** E.8.4c/d two-/three-series deliverable descriptions as what the implemented code actually does. Record the four labels for BOTH yields and MM, limits, source-definition uncertainty, and eventual renderer prerequisites, with explicit `NOT VERIFIED` for a separate reweighted-baseline authority.
2. `docs/memory/CURRENT.md` — add at most a **small concise pointer** near the active E.8/current-work context to the above requirement and unresolved fourth-branch source definition. **Keep the sole exact scientific NEXT unchanged**: audit the first blocking private F.4 scientific-equivalence invariant to determine earliest valid resolution, conditional on ChatGPT diff acceptance, user commit/push, and pushed-state synchronization. Do not make this clarification a NEXT, an E.8 phase, or a Method-A blocker resolution. Preserve accepted F.1 deferral.
3. `docs/memory/MEMORY.md` — at most one compact durable clarification/pointer if necessary to ensure the distinction survives a new chat: four *requested* slots, separate `rw` definition unverified, existing E.8.4 renderer three-way, no forged fourth result. No redundant roadmap rewrite.
4. `docs/memory/roadmap/STATUS.md` — at most one compact cross-reference, keeping the approved E.8.2 -> E.8.3 -> F.6.3 -> E.8.4 -> final E.8 -> F.6.4 order and every accepted/blocked status unchanged. Do not put a competing active NEXT here.
5. `docs/memory/manifest.json` — update exclusively through tracked manifest tool once all edits/new contract are placed.
6. `docs/memory/phases/e8-four-way-presentation-memory-consistency-task-contract.md` — this exact supplied contract, preserved byte-for-byte. **This contract is the only new Markdown file.**

No other path is allowed. In particular **do not edit** `src/`, `testing/`, `AGENTS.md`, `docs/memory/USER.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, closed evidence records, F.4/F.6 artifacts, producer/renderer code, farm ownership, or unrelated memory. `CURRENT_HANDOFF` is exceptional-only; normal cross-chat resumption follows CURRENT and the canonical roadmap.

## 4. Before/after and scientific ownership

Before: E.8 roadmap/current plot producer require/use `Y0`/`YA`, `MM_0`/`MM_A`, `K+Lambda` SIMC comparisons with the existing absolute-SIMC restrictions; the explicit *separate reweighted-baseline final yield/MM presentation requirement* is absent. The current E.8 status is `ACTIVE`; baseline-five `Y0/MM_0` is `CLOSED / RUNTIME VALIDATED`; private canonical-five F.4 is `BLOCKED`; final E.8/F.6.4 `BLOCKED`.

After: the **user's four-slot request is durably explicit and discoverable**, and the lack of an independently validated fourth final output is clearly identified. Existing two-/three-series implemented pages, accepted narrow results, and all scientific/runtime status remain exactly as before. This checkpoint **does not define** a correction, a new yield, a fourth MM histogram, a new normalization, a model-fit gate, or a new experimental acceptance criterion. In the new chat the current F.4 prerequisite audit remains the sole NEXT; the four-way requirement is a future presentation constraint to resolve when its scientific definition and upstream outputs can be source-verified.

Preserve current `no_empirical_residual` and all frozen production/subtraction/yield/error/acceptance/binning/SIMC inputs. Method A is private/non-production, parent-t preserving `w0 -> w0*C`; Method B remains diagnostic only and numerically excluded; empirical residual fits remain disabled. Original F.4 exact equality must not be bypassed, renormalized, or relaxed. F.1 archaeological search remains conditionally `DEFERRED` unless the exact prerequisite audit makes it materially necessary. Absolute-SIMC luminosity/effective-charge units remain `BLOCKED`: do not claim quantitative four-way amplitude agreement or normalize SIMC for display.

## 5. Positive, negative and regression acceptance checks

**Positive:**

- Both *yields* and *MM distributions* explicitly list all four requested categories under an E.8-specific canonical roadmap section, with the fourth result marked unverified as to authority, producer, and independence.
- Existing E.8.4 baseline/Method-A/SIMC three-way and `Y0`/`YA` code/documentation are accurately described as the **current implementation**, not retroactively four-way. Distinguish `B_pi^A` reweighted *pion background* from an independent final `MM_rw`/`Y_rw` branch.
- The roadmap pointer, MEMORY pointer (if added), CURRENT pointer and STATUS pointer are consistent: exactly one authoritative NEXT, owned by CURRENT; no implicit new stage or farm gate.
- Baseline-five accepted 135-cell `Y0/MM_0` closure, Left/lowe E.8.2/E.8.3/F.6.3/E.8.4 accepted narrow scopes, F.4 exact mismatch, absolute-SIMC units blocker and F.6.4 nonpromotion all remain intact.
- Source-identity language uses `SOURCE VERIFIED`, `RUNTIME VERIFIED`, `INFERENCE`, `MEMORY ONLY`, `NOT VERIFIED` correctly; user-defined presentation intent is not mislabeled as farm verification.

**Negative:** reject any claim that `Y_rw`/`MM_rw` currently exist in authoritative producer output, that they are independently extracted, that they are mathematically distinct from Method A, that `B_pi^A` is already such a final output, that current pages are four-way, or that absolute SIMC normalization was verified. Reject any fabricated fourth histogram/yield, child normalization, new `DeltaY` uncertainty, new science formula, automatic production promotion, changed closure, or altered CURRENT NEXT. Do not treat the user's notation as existing JSON or ROOT field names. Do not relabel prior three-way validated pages as failing due to the clarified future request.

**Regression:** audit all affected docs for inconsistent uses of *three-way*, *four-way*, *reweighted* and *baseline*, distinguishing historical delivered from future requested, and verify all relative links and memory manifest entries. Check the `no_empirical_residual`, Method A/B boundaries and F.4 gate text remain unchanged. Do not refactor adjacent records merely for wording cleanup. If a real scientific-ownership conflict is discovered, report it and stop rather than silently choosing one interpretation.

## 6. Deterministic local and memory health checks

Use a discovered working local Python executable `<PYTHON>`; do not assume a path or run ROOT/farm. After editing:

```sh
<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/check_memory_health.py --root .
<PYTHON> -B tools/memory_bootstrap.py --root . --json

git diff --check
```

Report exact pass/fail and any warning, `CURRENT.md`, `MEMORY.md` and `CURRENT_HANDOFF.md` byte counts. `CURRENT.md` had previously approached the 8-KiB soft boundary: maintain brevity, consolidate only narrowly if required to avoid an actual warning, and do not sacrifice active-state accuracy for cosmetic limits. Hard errors or material ambiguity STOP `BLOCKED`.

Verify only allowlisted paths changed, the supplied contract is unchanged, unrelated state is byte-for-byte preserved, Git HEAD/refs/index did not move, and no compiler/generated artifacts remain. Capture initial/final status and manifest/bootstrap outcomes. Do **not** claim farm/ROOT/PyROOT validation.

## 7. Actual-diff handoff and hard stop

Stop after local edits and deterministic checks for **ChatGPT actual-diff review**. Produce a single complete temporary review bundle in the repository root, including tracked changes in all allowlisted existing paths and the **entire new-file addition** of the supplied contract via `git diff --no-index /dev/null <contract>` (not merely ordinary `git diff`), with no staging merely for review. Preserve any existing review bundle by choosing a new name. Report exact file inventory, byte count and SHA-256, plus full health and preservation results. Do not commit, push, run the Jefferson Lab farm, or initiate the current F.4 audit in this checkpoint.

After independent ChatGPT diff PASS, **only the user** may commit and push accepted scoped changes. After push, ChatGPT must verify live SHA, gate-relevant blobs and `CURRENT`/NEXT continuity. Only then does the next analysis chat proceed directly to CURRENT's F.4 prerequisite audit without a redundant memory-reconciliation cycle.

**STOP `BLOCKED`** on wrong branch/identity, material CURRENT/source drift, any unauthorized file edit, failed manifest/health/integrity check, ambiguity requiring a new scientific definition, acceptance-status mutation, competing NEXT, fabricated fourth producer, or required change beyond this memory-only allowlist. No workaround or automatic science redesign.
