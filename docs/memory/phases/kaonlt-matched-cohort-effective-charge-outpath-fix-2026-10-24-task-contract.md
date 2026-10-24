# KaonLT — Eighth-failure effective-charge OUTPATH correction (one bounded fix)

**Task:** Repair the demonstrated `baseline_root_outpath_mismatch:effective_charge` and ensure the **same real no-mode effective-charge CSV writer** is confined to the already-owned private output in the matched-cohort child. One narrow correction, one combined local review/test receipt. No new science phase, helper framework, separate smoke task, background cleanup, physics adjustment, farm run, or commit/push by Codex.

**Anchor:** GitHub `test` at `92af408d0431d2955fc6e8544edebad3b45df5dc`; re-establish actual local `test`, `origin/test`, HEAD and complete worktree/index state before touching files. Source and farm provenance take priority over the prior assistant's erroneous `Farm readiness: PASS`. Farm worktree state is **NOT VERIFIED** in this contract.

## Farm evidence and first failed invariant

The user supplied the **eighth** attempt's original gate-status JSON and read-only installed-`ltsep` probe. These are **RUNTIME VERIFIED as supplied**, not accepted science:

- Source argument `92af408d0431d2955fc6e8544edebad3b45df5dc`; failure `baseline_root_outpath_mismatch:effective_charge` at `testing/run_pion_local_t_matched_cohort_gate.py::probe_private_callers`, inside `overlay_path_validation`; `analysis_started=false`, `evidence_admissible=false`.
- Installed Python `/usr/bin/python3.9`; source/collector preflight passed. Timings: source ~20.89 s, installed package ~3.11 s, overlay validation ~13.39 s, final preservation ~0.48 s; no shared-output snapshot, ROOT staging, child, or collection started. `shared_inputs_outputs_preserved` remained **NOT VERIFIED** in the failed receipt.
- Actual ordinary no-mode Root for `src/setup/findEffectiveCharge.py`: `LTANAPATH=/group/c-kaonlt/USERS/trottar/lt_analysis`, `OUTPATH=/group/c-kaonlt/USERS/trottar/hallc_replay_lt/UTIL_KAONLT/OUTPUT/Analysis/None`, `ANATYPE=Kaon`.
- Actual ordinary Plot-mode Root: `OUTPATH=/group/c-kaonlt/USERS/trottar/lt_analysis/OUTPUT/Analysis/KaonLT`. After the **exact** `findEffectiveCharge.py` `None` substitution, the no-mode OUTPATH is absolute and resolves to the Plot output, but is **not lexically equal** to `LTANAPATH/OUTPUT/Analysis/KaonLT`. The supplied Python diagnostic returned `absolute=true`, `lexical_equal=false`, `baseline_resolved_equal=true`, `derived_resolved_equal=true`. `/group` canonically maps to `/u/group`; the artifact outdir canonicalizes to the Lustre volatile output tree.

**SOURCE VERIFIED at the anchor:** `probe_private_callers()` enforces the erroneous baseline literal equality at its `base_raw` guard and has a separate private `raw` guard. `src/setup/findEffectiveCharge.py` calls `create_ltsep_root(os.path.realpath(__file__))`, substitutes `None`, then writes `OUTPATH/table_kaon_<KIN>.csv` through `getDataTable.py`/`utility.py`. `prepare_runtime_overlay()` changes only **copied** `LTANAPATH`; baseline `UTILPATH` stays stable. Consequently, **private no-mode writer destination is NOT VERIFIED and may still resolve to shared output**. The fixture `Root` incorrectly constructs *both* modes' OUTPATH from `LTANAPATH`; it does not model the real no-mode `UTIL_KAONLT` anchor, explaining why previous local passes did not establish this farm path contract.

## Workflow and scope

1. Read tracked root `AGENTS.md`, then exactly `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, `docs/memory/USER.md`; then `MAINTENANCE.md`, `CODEX.md`, the existing eighth-failure incident/owner contract, and only relevant source. Emit full KaonLT health check. Preserve October 20 accepted science and all closed statuses; update CURRENT to record the eighth failed gate and keep the scientific NEXT (fresh same-cohort Center/Left lowe canonical t1/t2/t3 A/B plus signed events, existing seven-strength C and six-parent MM/variance).
2. **Baseline resolver:** For no-mode, evaluate the same actual `OUTPATH` normalization consumed by `findEffectiveCharge.py`; require an absolute, existing/authorized canonical destination equal to ordinary Plot/`baseline['OUTPATH']` and authorized shared artifact root. Do **not** demand the no-mode lexical origin be `LTANAPATH`. Retain no-mode spelling/source role in the receipt so non-equivalent configurations fail. Keep Plot-mode strict and retain caller, cwd, installed/copy package identity, unique config, all stable input roles, and installed-source preservation.
3. **Private resolver is a distinct hard gate:** Determine from the *copied* ltsep normal import, no-mode and Plot probes, and the real writer expression where `table_kaon_<KIN>.csv` would be written by the isolated child. Require that both modes' actual writable destinations canonically equal the **owner-created private** output directory and cannot escape via symlink/alias. Do not make private accept the ordinary shared directory just because ordinary baseline canonical paths agree. Verify parent topology and source->sink; no child until all checks pass. Do not change `VOLATILEPATH`, shared external roles, installed ltsep, installed config, external links, or trust environment-only redirects.
4. **If the existing copied-`LTANAPATH` overlay does not redirect no-mode:** choose the smallest justified safe remedy **within this one task** based on real `findEffectiveCharge.py` source and existing frozen-source/provenance rules. Prefer an owner-contained private resolution consistent with all existing constraints, but do not insert monkeypatches, rewrite shared links/config, use direct manual farm patches, or weaken byte-preservation. A narrowly scoped path-only correction to `src/setup/findEffectiveCharge.py` is **conditionally permitted only if** Codex first establishes that no owner-only solution preserves the contract, demonstrates ordinary production resolves to the same destination and returns/format/physics are unchanged, and updates all genuinely required source-identity/test allowlists without defeating frozen evidence or the historical source check. If those proofs cannot be made, **stop `BLOCKED`** with the first exact writer/authority invariant and no speculative change. Never silently modify other frozen production paths.
5. **Consolidated deterministic tests:** emulate the *actual* two distinct ltsep Root modes (ordinary no-mode anchored in `UTIL_KAONLT/OUTPUT/Analysis/None` with symlink to normal shared output, ordinary Plot anchored in `LTANAPATH`); exercise baseline/pass, wrong canonical baseline/fail, private no-mode-to-shared/fail, private no-mode-to-private/pass, bad alias/package/config/fail, and the full acquisition CLI pre-child path with mock child proving no science execution on failures. Use real subprocess imports where supported. Include a source-level check of the effective-charge table writer, input preservation, cheap-to-expensive order, current timing/status behavior and the existing regression suite. Mark JLab filesystem/installed `ltsep` and physical ROOT behavior **NOT VERIFIED** locally.
6. **One consolidated result:** return changed-path inventory, complete tracked diff and full new-file additions in one temporary root `kaonlt_review.diff`, checks/test counts/skips, exact source-path chain CLI -> baseline/private no-mode -> `findEffectiveCharge` CSV sink -> child -> preservation/collector, and all remaining farm-only unknowns. No independent checker or second contract. No farm command or automatic retry. User alone commits/pushes after ChatGPT independent full-diff review; pushed-state synchronization is required before any next farm readiness decision.

**Primary source scope:** `testing/run_pion_local_t_matched_cohort_gate.py`, `testing/test_run_pion_local_t_matched_cohort_gate.py`. Conditional path-only `src/setup/findEffectiveCharge.py` (plus directly necessary focused tests and owner source-provenance/allowlist changes) only under item 4's proof. Memory: `docs/memory/CURRENT.md`, existing `docs/memory/investigations/kaonlt-matched-cohort-farm-readiness-failure-2026-10-22.md`, `docs/memory/manifest.json`, and **this single task contract**. Do not create a new framework or widen scientific source modifications. Preserve source-dependent values/formatting, effective charge calculation, CSV row contents, normalization, cuts, models, weights, uncertainties, A/B/C algorithm, and non-production Method A/independent excluded Method B roles.

**Hard stop:** If the private no-mode writer cannot be proven private, this is a genuine output-isolation blocker, **not** a lexical comparison to bypass. Report it together with the narrowest source-supported next remedy rather than authorize a ninth acquisition.


## Local implementation and consolidated receipt

Status: `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`. Starting HEAD and local
origin/test remain `92af408d0431d2955fc6e8544edebad3b45df5dc`; independent complete
diff review is NEXT. No commit, push, remote-ref update or farm command occurred.

SOURCE VERIFIED: the baseline guard accepts existing canonical no-mode/Plot
equivalence without demanding a no-mode LTANAPATH lexical anchor. Plot stays
strict. Stable UTILPATH cannot redirect the no-mode writer through a copied
LTANAPATH-only config: the minimal permitted producer correction derives its
CSV OUTPATH from LTANAPATH/OUTPUT/Analysis/<ANATYPE>LT. The owner evaluates the
actual audited producer AST, checks ordinary equivalence to the prior sink,
and requires the private sink to equal the owner-created matched-output tree.
Shared or redirected backing paths fail before the child. Read-only no-mode
helpers never use Root.OUTPATH; their inputs remain unchanged. Only this producer
was added to the current owner's committed-successor allowlist. The historical
frozen resolver/source check remains unchanged and grants no new attestation.
The existing writer-audit status also records the producer's source identity.

The three consolidated checks:

1. Pathing and imports: WSL Python 3.12.3 owner suite, 63 tests in 59.917 s:
   62 PASS, one skip (real scientific dependencies unavailable in WSL).
   Command: `python3 -B -m unittest
   testing.test_run_pion_local_t_matched_cohort_gate -q`.
   Actual distinct no-mode replay-UTIL_KAONLT and Plot-LTANAPATH patterns,
   /group-to-/u/group-shaped ancestor aliases, real subprocess imports and
   cwd/package/config negatives are covered. An extracted unchanged real CSV
   sink writes only inside the owner-created private output. Wrong baseline,
   shared writer and private backing alias escapes fail. Windows Python 3.12.10
   ordinary numpy/pandas/dateutil/ctypes/zoneinfo imports separately PASS; that
   import check is not a replacement for the skipped WSL dependency fixture.
2. Code infrastructure and preservation: WSL frozen resolver/history, collector
   and memory-health regression suites, 107/107 PASS in 70.410 s.
   Commands: `python3 -B -m unittest
   testing.test_attest_pion_chain_ltsep_default_root
   testing.test_collect_pion_hgcer_validation_bundle
   testing.test_memory_health -q`.
   Windows event-owner suite plus the focused producer AST preservation test,
   15/15 PASS in 1.807 s: `python -B -m unittest
   testing.test_run_pion_local_t_event_validation_gate
   testing.test_run_pion_local_t_matched_cohort_gate.PrivateBoundary.test_effective_charge_change_is_path_only_and_old_sink_is_rejected -q`.
   After removing only the old/new output-path statement, the producer AST is
   identical to its frozen predecessor. Charge calculations, stdout, input paths,
   CSV row/format logic and all physics remain unchanged. Syntax compilation of
   the three changed Python files, with caches outside the repository, and
   `git diff --check` PASS. No standalone historical attestation is newly granted.
3. Existing status/log debugging: included in the owner suite above, not counted
   again. The actual acquisition CLI parser -> execute_gate -> resolver preflight
   runs real baseline/private probes with a mock child. Shared writer failure
   returns 1, names the existing status file, records overlay_path_validation,
   analysis_started=false and evidence_admissible=false, performs preservation,
   and never reaches shared snapshot, staging, child, collector or delivery.
   Existing stage-timing, traceback, cleanup, full-byte comparisons and status-I/O
   failure tests remain. Fixture collector stdout is not a farm execution.

Execution chain and the owner-only-remedy proof are recorded in the existing
incident's October 24 section. The review bundle contains all actual tracked
diffs and this complete new-file addition; no staging is needed for review.
Changed paths: the matched owner and its test, src/setup/findEffectiveCharge.py,
CURRENT, the existing October 22 incident, manifest, and this supplied contract.
Only those six tracked files may differ from the 985-file startup byte snapshot.
The final preservation check and manifest/health results are reported in chat.

JLab filesystem/installed ltsep, actual Python 3.9 and dependency closure,
physical ROOT, full child/runtime preservation, fresh matched scientific outputs,
PDFs and production acceptance remain NOT VERIFIED. October 20 accepted science,
Method A's detached role and independent/excluded Method B are preserved.
No acquisition rerun or farm-readiness approval follows from this local receipt.
