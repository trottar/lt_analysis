# E.8 baseline-only canonical-five extracted-yield validation owner — task contract

## 1. Objective and starting authority

**Task class:** one narrowly scoped operational/validation implementation, not a physics change or a Method-A repair.

Repository: `trottar/lt_analysis`, branch `test`.
Audited remote HEAD: `af45126899c0d89707b29b95f7138b3388d432cd` (2026-10-08). Before work, Codex must establish its actual **local** branch, HEAD, `origin/test`, and complete worktree/index status. Require this exact HEAD and local `origin/test`, or **STOP BLOCKED**, return the mismatch and request a fresh audit. Do not reset, clean, stash, stage, overwrite or discard the user's unrelated dirty changes (notably the three previously reported tracked docs under `docs/memory/`).

First read the exact tracked root-`AGENTS.md` startup core in order: `AGENTS.md`, `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, `docs/memory/USER.md`. Follow `docs/memory/MAINTENANCE.md`, `docs/memory/CODEX.md`, `docs/memory/TOOLS.md`, `docs/memory/COMMUNICATION.md`, and this contract. Check relevant accepted E.8.2 Left/lowe runtime evidence and E.8 full-procedure roadmap, without reopening their closed scopes.

**Goal:** A single tracked, deterministic owner/profile/checker/collector path that, once independently reviewed, pushed, and farm-executed by the user, can validate the existing ordinary full **Q4p4W2p74 baseline** analysis for the **exact canonical five** and preserve the **same producer-owned data and K-Lambda SIMC yields and MM support** passed to the existing L/T-separation chain. It must not require or claim an accepted Method-A F.4 lineage or Method-A candidate comparison. This closes no runtime gate until applicable farm evidence passes.

## 2. SOURCE VERIFIED audit basis (must recheck before editing)

- `testing/run_e8_2_left_lowe_scientific_audit_gate.py` owns the accepted isolated **Left/lowe only** E.8.2 farm test and `-d` debug launcher; its profile declares five potential settings but its requested set and run/checks cover only Left/lowe.
- `testing/run_e8_4_fix5_canonical_five_plot_gate.py` owns a **full five-setting** isolated run (`./run_Prod_Analysis.sh 4p4 2p74`), with mandatory F.2/F.3/F.4 Method-A candidate/equivalence preflight before analysis. Its F.4 failure remains independently `BLOCKED`; do not skip or modify that gate in its owner.
- `src/main.py` Step 6 calls `find_yield_data`, `find_yield_simc`, then E.8.2 post-yield presentation finalization; Step 6 ratio feeds ordinary cross-section workflow; the ordinary high-epsilon Step 7 invokes existing L/T procedures.
- `src/binning/calculate_yield.py` produces baseline Y0 from its cut-selected final MM; its E.8.2 *wide/no-MM-cut* display is not that integrated histogram. The private F.6.3 builder fails closed to an unavailable parallel branch when Method-A candidate authority cannot be reconstructed; it does not replace public Y0.
- `src/physics_lists.py` writes downstream data/SIMC per-setting `src/kaon/yields/yield_data.*.dat` and `yield_simc.*.dat`, with **four-decimal formatting** of yield/error and `(phi_bin,t_bin)` coordinates. Do not treat rounded text as unrounded scientific precision.
- `src/binning/xsect_support.py` emits `kaon_xsect_support_Q4p4W2p74_lowe.npz` and `_highe.npz` (0th-iteration case), with data/SIMC `mm` arrays, error arrays and bin edges from existing yield producers. Preserve existing iteration names/semantics if actual runtime differs; do not silently take another iteration's cached data. The SIMC source comes from the existing `iter_weight` model for K-Lambda in `src/main.py`.
- `src/normalize/get_eff_charge.py` takes the existing SIMC `normfac/Ncontribute`; **absolute** SIMC-to-data yield units remain **BLOCKED** as `SIMC_normfac_luminosity_and_charge_units_not_source_proven`. No synthetic scale, implicit floating normalization or claim of absolute agreement is allowed.

## 3. Exact implementation scope

Permitted NEW tracked files:

- `testing/run_e8_2_canonical_five_baseline_yield_gate.py`
- `testing/pion_hgcer_validation_bundle_profile_e8_2_canonical_five_yields.json`
- `testing/test_run_e8_2_canonical_five_baseline_yield_gate.py`
- `testing/test_pion_hgcer_validation_bundle_profile_e8_2_canonical_five_yields.py`
- `docs/memory/phases/e8-2-canonical-five-baseline-yield-gate-task-contract.md` (this file, unchanged)

Permitted MODIFIED tracked files **only if needed to keep one accurate, push-stable NEXT**:

- `docs/memory/CURRENT.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/manifest.json`

Do not modify `src/`, `run_Prod_Analysis.sh`, `set_SymLinks.sh`, `farm_env/`, existing owners or profiles, collector, particle-subtraction/simulation/model files, production yield formats, tests outside this allowlist, `MEMORY.md`, or historical evidence. If a required feature cannot be implemented within the allowlist, STOP `BLOCKED` with the earliest concrete reason; do not broaden the contract implicitly.

## 4. Intended after-behavior: independent baseline owner

Implement one new, source-pinned owner using existing proven operational isolation/collector primitives, with no scientific source change. Source reuse should be narrow and explicit; do not copy or invoke the candidate F.4 preflight, candidate staging, candidate restoration, or candidate scientific-equivalence validation. Do not modify the existing canonical-five Method-A owner.

Fixed setting inventory, no Right/lowe: `Left/lowe`, `Center/lowe`, `Left/highe`, `Center/highe`, `Right/highe`. Validate both epsilon analyses and the exact two ordinary full completion markers. Launch only **`./run_Prod_Analysis.sh 4p4 2p74`** inside a newly created **detached disposable worktree**, not the user's primary farm checkout, and without `-d` or new launcher flags.

The owner must preserve all already established safety invariants: exact pushed source, primary checkout status/ref/content fingerprints, isolated `ltsep` copied configuration and `LTANAPATH` caller probes, non-mutating external SIMC-link preflight, bounded owner-created worktrees/cleanup even on failure, full subprocess log and per-attempt status, success/failure receipt, freshness and SHA-256, completion markers, distinct collector worktree/source identity, complete ZIP verification, primary checkout and installed-ltsep preservation. Never repair external links or mutate the primary checkout. Do not copy any analysis-produced model/parameter/working files back to primary checkout.

**Important:** the ordinary analysis may generate an `E.8.4 unavailable` diagnostic because Method-A F.4 equivalence is blocked. This is *allowed only in the new baseline owner's scope*, after checking that baseline Y0/E.8.2 survived independently. E.8.4 availability/unavailability is not a baseline-production success criterion and its reason must be recorded honestly. The original Method-A five-setting owner remains unchanged and `BLOCKED`.

## 5. Producer-to-evidence yield and MM gate

Require and validate, at minimum:

1. Each of five settings' fresh baseline procedure PDF and strict page manifest; for each canonical t parent, all nine phi indices explicitly represented in six E.8.2 page families, with exact setting/t/phi semantics, full inventory, `no_empirical_residual`, no renderer failures, and no silently missing/invalid populated children. Explicit frozen skip/zero child states must be recorded and judged under their **existing** accepted policy; do not relabel an invalid child as an accepted yield.
2. Fresh, strict full-analysis JSON for each epsilon with exact setting inventory and runtime configuration, plus matching `no_empirical_residual` correction ledgers (JSON and CSV) and provenance.
3. Both existing 0th-iteration `kaon_xsect_support_Q4p4W2p74_{lowe,highe}.npz` outputs with exactly source-provided canonical t/phi edges and setting matrices `data_mm_<setting>_values/errors`, `simc_mm_<setting>_values/errors`, `mm_edges`. Reject wrong/missing/ambiguous iteration, epsilon, setting, shapes, nonfinite bins, or mismatched axes. Do not modify the existing NPZ writer.
4. Five exact **data** and five exact **SIMC** downstream `src/kaon/yields/*.dat` producer files written by `src/physics_lists.py` in the detached worktree during the same analysis. Resolve their real names using the source's existing polarization, epsilon, Left/Center/Right conventions; require unique matches, fresh content, exact canonical `(phi_bin,t_bin)` coordinates, valid/explicit skips and nonfinite rejection. These source files exist inside the disposable worktree and must be **copied byte-for-byte into uniquely named owner-owned evidence artifacts** under the already resolved analysis OUTPATH **before worktree cleanup**. Log each source-relative basename, copied basename, SHA-256 and byte count. Do not use the user's ordinary checkout as a yield-file source.
5. A **read-only numerical consistency audit**: for each valid `(t,phi)` compare the producer `.dat` baseline Y0 against the normal-bin integral of the matching *cut-selected* `data_mm` support; compare its SIMC table yield against the matching SIMC MM support with the existing nonnegative SIMC-yield convention, respecting the known four-decimal `.dat` formatting as a bounded textual quantization error. Detect coordinate swaps, duplicates, omissions, negative/undefined cases, bin/window mismatch and unexpected differences. This computation is **checker-only**, never a new physics producer or replacement yield. Do not alter producer yields, uncertainties or models. If a source-established numerical rule is unknown (especially excluded/zero bins), stop and identify the invariant instead of picking an arbitrary tolerance.
6. Retain the per-cell data and SIMC **MM content/error and bin edges** from the existing NPZ, and the source-persisted table Y0, YSIMC, errors, identities, numerical gate outcomes and SHA references in one compact, versioned baseline-five owner receipt (JSON). Record the actual model-input/iteration identities available from the existing full analysis, `.hist` or model provenance; do not claim a freeze or parameter identity that was not verified. Do not generate interpretive agreement plots or scale SIMC to data in this first gate.

The independent validation receipt should indicate `SOURCE VERIFIED` only for checks actually performed in the owner, not assert full 3D/ROOT validation beyond the sourced farm artifacts. Ensure the reviewer can identify which per-setting data and SIMC yields will later be paired with private YA and future reweighted yields **without changing the existing cross-section inputs**.

## 6. Collector and deliverable

Define a distinct v4 compatible `generic_artifacts` profile for exactly five settings, scoped to this baseline owner. Include required baseline PDF/manifest, existing low/high full-analysis JSON and correction ledgers, the two NPZ support files, ten byte-for-byte producer yield tables under deterministic owner-owned names, one neutral owner summary/receipt and any required status/log companion evidence. Ensure collector handles one epsilon-level shared artifact appearing for more than one setting without identity ambiguity or duplicate-path failure. Do not invent a new collector mode. Verify source SHA, manifest, requested setting inventory, all required artifacts, summaries, hashes and final ZIP. If required NPZ/yield files cannot be delivered correctly by the unchanged collector, STOP `BLOCKED` and return the exact unsupported contract instead of weakening artifact completeness.

The expected returned object after a *later, user-run farm gate* is a self-contained SHA-pinned baseline-five validation ZIP, plus companion owner status/log/summary as required by the existing operational pattern. A failed gate produces only failed diagnostic evidence, not an accepted bundle or scientific yield comparison.

## 7. Deterministic local validation

Add tests with synthetic owner/collector/filesystem fixtures only. Cover: exact full launcher/cwd and completion markers; exact five settings/no Right-lowe; wrong branch/HEAD/origin and unexpected dirty state; isolation/path/ltsep/symlink preservation; source-check/ZIP failure; stale/wrong-epsilon/wrong-setting PDFs, manifests, JSON/CSV/NPZ; missing/duplicate yield tables; four-decimal quantization; error propagation not modified; cell `(phi,t)` ordering; support axes/matrix identity; negative/undefined SIMC bins and accepted skip/zero handling; Method-A candidate missing/mismatched cases **not** incorrectly blocking baseline-only validation; no calls to F.2/F.3/F.4 Method-A preflight; success and failure cleanup; owner summary/provenance accuracy. Use current regression suites for left-lowe E.8.2, canonical-five isolation, collector/profile and producer impact. Mock farm subprocesses; no real farm, ROOT/PyROOT or external environment repair.

Run source syntax/AST checks, all focused tests, required existing regression tests, `git diff --check`, complete diff/allowlist/other-state byte preservation. Record tests exactly rather than claiming complete runtime validation.

## 8. Gate sequence, memory, review and hard stops

Operational chain: approved immutable source+inputs -> source/checkout/collector preflight -> detached worktree+ltsep path isolation+external-link safety -> unchanged full Q4p4W2p74 launcher -> full low+high completion -> captured authentic yield `.dat` and existing MM support -> numerical owner checker -> baseline E.8.2 PDF/manifest/full-analysis/ledger provenance -> owner status/summary -> detached generic collector -> final source/checkout preservation -> verified ZIP and bounded cleanup.

A missing owner stage, checker, supported collector path or input contract means `BLOCKED`; no farm instruction before actual-diff review, user commit/push, and pushed-state synchronization. The exact canonical-five Method-A F.4 equivalence failure stays `BLOCKED` and must not be retried or waived. Absolute-SIMC amplitude comparability stays independently `BLOCKED`. Do not claim Method A improved SIMC agreement. Method A remains private/nonproduction, Method B diagnostic-only, no reweighted data branch is invented.

Update only permitted memory records for a concise, **push-stable sole NEXT**: after ChatGPT actual-diff acceptance, user publication and pushed-state synchronization, perform an authorized *narrow baseline-five farm gate* using the new tracked owner; returned evidence will determine whether baseline yields are accepted. Historical Left/lowe accepted closures remain unchanged and no canonical-five runtime closure is claimed. Keep `CURRENT.md` compact (<8192 bytes if possible), manifest regenerated/checked, ordinary health and bootstrap PASS. Report CURRENT/MEMORY/CURRENT_HANDOFF exact UTF-8 sizes and warnings/classification per MAINTENANCE.

Stop for ChatGPT actual-diff review. Supply one complete review bundle (tracked diff plus full no-index additions of new files), do not overwrite unrelated review bundles. Do not commit, push, run farm, add new live instrumentation, change scientific code, change the original F.4 owner, or authorize scientific interpretation/production promotion.
