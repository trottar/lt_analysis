# E.8.2 canonical-five baseline yield/MM runtime closure — repository task contract

## 1. Classification, authority and hard prerequisites

**Task class:** evidence-only repository closure of an already accepted Jefferson Lab farm gate. It is not a new scientific calculation, source fix, Method-A lineage repair, or instruction to rerun the farm.

Authoritative repository: `trottar/lt_analysis`, branch `test`. At the time of this ChatGPT audit, live `test` HEAD was:

`8e68d2acae91fd97622793df563aba0296403a48`

Before any edit, establish actual local branch/HEAD, `origin/test`, full worktree/index state, and whether tracked/untracked local changes exist. Require the exact above baseline and its matching `origin/test`; otherwise **STOP BLOCKED** and report the difference. Do not reset, fetch-and-overwrite, clean, stash, stage, change file modes, or discard unrelated local work.

Read the tracked startup core in exactly this order: repository-root `AGENTS.md`, `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, `docs/memory/USER.md`. Then follow `docs/memory/MAINTENANCE.md`, `docs/memory/CODEX.md`, the E.8 roadmap and the [baseline-five owner task contract](e8-2-canonical-five-baseline-yield-gate-task-contract.md). Use tracked evidence records only as needed for historical-scope distinctions. Root `AGENTS.md` and live source outrank this explanatory handoff.

**Accepted gate:** the independent E.8.2 `Q4p4W2p74` baseline-five yield/MM owner, from one full user-run farm attempt at the source above. Its scope is strictly baseline `Y0/MM_0` and corresponding authentic K-Lambda SIMC yields/support, not private `YA/MM_A` or production-promotion evidence. ChatGPT already accepted the uploaded runtime evidence in the immediately preceding review and rechecked the supplied archive, status, receipt and hashes. This contract authorizes recording that acceptance, not reassessing it or demanding farm artifact access in Codex.

## 2. Exact accepted farm evidence and identities

Attempt stem: `Q4p4W2p74_E8_2_baseline_five_8e68d2ac_20261009T223015Z`

| Supplied artifact | Size (bytes) | SHA-256 |
| --- | ---: | --- |
| Validation ZIP | 111668823 | `df5a64563c7322416f0895c52557bea69f49d73ac0cb47836b0635a27e7660c4` |
| Companion gate-status JSON | 20273 | `a8cdfa6ee35a461f3313f63c1b8faa7c4db5416c9dfffafff5f87b4f986fad29` |
| Companion run-summary JSON | 965570 | `856128f6c41f27e5a09a52aff01c37825c8d9de7279cda36d11009f3e0d2f583` |
| Companion child log | 292966650 | `f630c33b5d759c3c7541fe552f83524a63edb923ce7c9e51c05badeaeba6c282` |

The accepted owner receipt and independent ZIP review establish:

- `status=success`, `stage=complete`, `failure_reason=null`, analysis subprocess return code 0; analysis, artifact checks, collection and ZIP verification complete; disposable worktree cleanup complete.
- Primary checkout identity/ref/index/content and original installed-`ltsep` preservation passed, including final checks. The original full-analysis launcher executed only in the disposable isolated worktree. No claim that the ordinary checkout is clean now.
- Exact canonical settings: `Left/lowe`, `Left/highe`, `Center/lowe`, `Center/highe`, `Right/highe`; **no Right/lowe**. Two epsilon analyses with active profile `no_empirical_residual` and the existing 0th-iteration cut-selected MM support.
- ZIP CRC and manifest/source/provenance verification passed. Exactly **29 required evidence artifacts**: five procedure PDFs, five corresponding page manifests, two epsilon full-analysis JSON files, two JSON and two CSV correction ledgers, two source NPZ MM support files, ten byte-for-byte producer yield `.dat` tables, and one canonical owner receipt. The ZIP has 38 entries including directories and three collector control entries. All 29 packaged artifact SHA-256 identities matched its manifest; the packaged owner receipt matches the uploaded companion summary byte-for-byte.
- **Five PDFs, 97 pages each (485 pages)**. The five page manifests match the PDF counts, show zero renderer failures and all **18 E.8.2 stage/final pages per setting** (90 such records across five). Representative final-MM/stage-yield pages for each setting were rendered and inspected; do not claim all 485 pages underwent independent pixel-by-pixel visual acceptance.
- Exactly **135 `(t_index,phi_index)` cells** (three canonical t parents × nine phi children × five settings). All producer data and SIMC table-to-cut-selected-NPZ numerical gates passed: **135/135 data**, **135/135 SIMC**, with one-based producer `(phi_bin,t_bin)` mapped to zero-based NPZ child indices, normal-bin integration and exact existing four-decimal text formatting. A separate independent reintegration of the supplied receipt found **zero discrepancies**; source table uncertainties were not changed or rederived as a scientific replacement.
- Four signed-negative producer data yields were retained: `Left/lowe t2 phi7 = -0.0001`, `Left/highe t2 phi7 = -0.0017`, `Right/highe t1 phi2 = -0.0003`, `Right/highe t1 phi8 = -0.0001`. The existing nonnegative SIMC-yield convention gives **36 zero SIMC cells**. These are observed producer results, not grounds for invented normalization or a demonstrated physics defect. All 135 parent policies are recorded `valid`; no accepted skip/zero parent policy was needed in this attempt.
- Cut-selected MM window: `1.10–1.16 GeV`, **100 normal bins** per child. These are integrated numerical-support objects, not the wider E.8.2 stage-display histograms. Earlier diagnostic Lambda-window stage integrals are not authoritative extracted `Y0`.

The accepted source `src/binning/calculate_yield.py`, `src/physics_lists.py`, `src/binning/xsect_support.py`, existing page producers, frozen profile, and `testing/run_e8_2_canonical_five_baseline_yield_gate.py` remain unchanged.

Optional identification within the evidence record, if useful for future archive checking:

- Low-epsilon NPZ SHA-256: `c38968915a1e58a822986cea317b1246201210f669baca0e121b7f8cdee6cfd4`
- High-epsilon NPZ SHA-256: `d342f6b6aab40486f956e77128cd452392066f831661d556c04921c9d7b2f705`
- Procedure PDF SHA-256 by setting: Left/lowe `5e6bcdb838bb052a5285f965347ae1874809fd1c83941cbe0680dc47f1f05a86`; Left/highe `746b69faf30df08927468f18cb5ff0302813fa6619a0d29df0f36737b6e3e3e7`; Center/lowe `b9180adc9fb7890ba1526c11f61db0f43d0dae6f0319007a1848b099cb2b7dce`; Center/highe `c3b9941973fb3e52ddce6829743a86229d82b60a658fb005da62044f28cbfe81`; Right/highe `a4bdeadbee5d3d6b8b999fccf2badee34ad10d6754bee50c032bf05cf5601470`.

## 3. Accepted scientific scope and explicitly preserved blockers

Record `CLOSED / RUNTIME VALIDATED` **only** for the independent **E.8.2 canonical-five baseline extracted-yield and matching MM numerical/procedure gate** at the exact source/attempt above. This expands baseline-five runtime coverage beyond the earlier accepted Left/lowe-only E.8.2 closure but does not modify that historical record or retroactively close any other component of E.8.

Preserve every existing unrelated accepted phase and scientific ownership, especially:

1. Original **private Method-A canonical-five F.4 scientific-equivalence** mismatch at `$.parents[0].absolute_baseline_parent_sum`: `BLOCKED`. The original five-setting Method-A owner is **not retried**. No `YA/MM_A` comparison, current/reviewed F.4 equivalence, numerical Method-A improvement or production promotion is accepted by baseline-only success.
2. **Absolute SIMC amplitude/units** remains independently `BLOCKED`: `SIMC_normfac_luminosity_and_charge_units_not_source_proven`. Same-model SIMC-to-data numerical support exists, but no absolute-scale agreement, rate or luminosity conclusion is authorized.
3. Final E.8 and F.6.4 remain `BLOCKED`, with Method A private/nonproduction and Method B diagnostic-only, numerically excluded. The active `no_empirical_residual` profile and frozen production/normalization/yield/uncertainty/cuts remain intact.
4. The accepted prior E.8.2 Left/lowe, E.8.3 Left/lowe, F.6.3/E.8.4 Left/lowe, and other closed scopes are not reopened or downgraded. E.8.4 unavailable-page records in this bundle do not establish private branch availability or acceptance.
5. Spot-checked PDF readability, structural page inventories, table/NPZ closure and runtime owner success are distinct evidence claims. Do not overextend the result to cross-section accuracy, total systematics, detector-response causality or physical interpretation of the negative/zero yields.

## 4. Exact allowed repository changes

**New files** (this contract is the sole task contract):

- `docs/memory/phases/e8-2-canonical-five-baseline-yield-runtime-closure-task-contract.md` — this text, unchanged after placement.
- `docs/memory/evidence/e8-2-canonical-five-baseline-yield-runtime-closure-2026-10-09.md` — concise durable accepted-evidence record using exact verified hashes, numerical quantities, and scope limits above, explicitly crediting user farm artifacts and independent ChatGPT review (Codex did not execute the farm or accept evidence independently).

**Modified files**, only as necessary:

- `docs/memory/CURRENT.md` — record this new narrow closed baseline-five result and link evidence; replace the consumed farm-gate NEXT with a push-stable substantive science NEXT.
- `docs/memory/roadmap/STATUS.md` — reflect baseline-five E.8.2 gate closure while distinguishing the blocked broader Method-A canonical-five/full E.8 gate.
- `docs/memory/manifest.json` — generated manifest update for changed repository memory.

No other files are allowlisted. In particular: no `src/`, `testing/`, model, procedure page, collector/profile/owner, prior evidence, `AGENTS.md`, `MEMORY.md`, handoff, or unrelated historical document changes. Do not create a separate acceptance-review artifact or new design contract. Preserve all pre-existing worktree/index/untracked state.

## 5. CURRENT NEXT: one substantive scientific gate after closure

Use one sole, push-stable NEXT, **conditional on ChatGPT actual-diff acceptance, user commit/push, and pushed-state synchronization**:

**NEXT — E.8 baseline-five extracted-yield scientific QA:** Consume the now accepted baseline `Y0/MM_0` and K-Lambda SIMC cut-selected MM evidence, audit per-setting/per-parent `(t,phi)` yield distributions and uncertainty/reporting conventions, and characterize four signed-negative data cells and 36 zero SIMC cells without adjusting the data, normalizing SIMC amplitudes or reinterpreting stage integrals as authoritative yields. Establish which baseline-only shape/coverage or scientific comparisons are meaningful and what extra independently authorized evidence would be needed for private `YA/MM_A` comparisons. Return one physics-oriented plot/yield comparison/readiness result or one directly evidenced blocker. Do not perform or authorize a private Method-A F.4/YA comparison until its independent lineage/equivalence gate is resolved.

The NEXT must not be another paperwork-only reconciliation task, full Method-A owner retry, physics reimplementation or unsupported claim of complete E.8 closure. `CURRENT.md` should remain compact; its pre-edit size is approximately **8063 UTF-8 bytes**, near its **8192-byte soft warning**. Consolidate only replaced/consumed local wording when necessary, not unrelated chronology. A nonblocking soft warning is recorded per MAINTENANCE; do not start a broad cleanup.

## 6. Deterministic local checks and handoff

After the scoped memory edits, run the tracked manifest write then check commands, ordinary memory health (`tools/check_memory_health.py`), bootstrap, Markdown/link/sole-NEXT checks, and `git diff --check`. Discover an actual working Python interpreter (no assumed workstation path). Independently audit the exact allowlist and changed-file content, confirm accepted evidence identity and source/no-physics-diff boundary, check expected Git modes `100644`, record before/after HEAD/ref/index/status and unrelated file-byte preservation, and classify any memory warnings under MAINTENANCE.

Report exact `CURRENT.md`, `MEMORY.md`, `CURRENT_HANDOFF.md` UTF-8 byte counts, manifest/health/bootstrapping results and any unresolved warning. Never turn the previous farm result into a claim of a new farm run or revalidation by Codex. No staging, commit, push, remote-ref updates, Jefferson Lab farm execution, or scientific-source changes.

Deliver **one complete actual-diff review bundle** containing tracked changes and full no-index contents of the new closure evidence and this contract; preserve other review bundles and unrelated dirty state. **Stop for ChatGPT actual-diff review.** The user alone performs accepted Git publication after that review; ChatGPT then checks pushed-state continuity before the next scientific milestone.
