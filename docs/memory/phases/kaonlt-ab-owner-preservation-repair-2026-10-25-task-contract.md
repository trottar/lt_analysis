# KaonLT A/B owner: reviewed successor and private-output repair

Status: `ACTIVE`. One bounded repair of the established A/B owner; no new scientific phase, acquisition owner, farm operation or publication.

Repository destination: `docs/memory/phases/kaonlt-ab-owner-preservation-repair-2026-10-25-task-contract.md`.

## Starting authority and demonstrated failures

Starting `test` HEAD and local `origin/test`: **`1144fe56b68b0980f7f7ce59887987b10e0f3d9a`**, independently checked against live GitHub `test` on October 25, 2026. The workstation checkout was clean before contract placement. Farm worktree state is NOT VERIFIED here. Re-establish branch, HEAD, index, full tracked/untracked state and pending contracts before editing; read root `AGENTS.md`, `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, `docs/memory/USER.md` in that order. Then read MAINTENANCE and the two existing October 19 comparison contracts plus the October 25 fresh-source contract. Stop on an incompatible HEAD or overlapping local change. Do not fetch, reset, clean, stash, stage, commit or push.

The user supplied two blockers, corroborated by current source:

1. `preflight()` compares against historical `554ad84630ad77920ce17954d99f73df8ed2e6ac` and permits only the owner/profile/test plus memory. Live history includes these 13 legitimate additional paths:

```text
src/binning/calculate_yield.py
src/setup/findEffectiveCharge.py
testing/pion_local_t_event_validation_capture.py
testing/pion_local_t_event_validation_gate_profile.json
testing/pion_local_t_event_validation_replay.py
testing/pion_local_t_transfer_variance_comparison.py
testing/run_pion_local_t_event_validation_gate.py
testing/run_pion_local_t_matched_cohort_gate.py
testing/test_pion_local_t_event_validation.py
testing/test_pion_local_t_matched_cohort_science.py
testing/test_pion_local_t_transfer_variance_comparison.py
testing/test_run_pion_local_t_event_validation_gate.py
testing/test_run_pion_local_t_matched_cohort_gate.py
```

2. `comparison_artifact_preexists`: all four accepted Center/Left lowe A/B JSON/PDFs already occupy their ordinary names. Preserve them in place. Their existence is supplied evidence; their bytes and current farm identities were not independently inspected in this audit.

## Feasibility finding and first implementation gate

Current `--outdir` is **not an isolation switch**: `execute_gate()` requires the configured ordinary artifact root, and frozen `debug.prepare_debug_output_link()` requires the worktree OUTPUT link to target `VOLATILEPATH/OUTPUT`. The producer writes directly to its supplied `outpath`. Arbitrary `--outdir`, VOLATILEPATH changes, or a bare private link would fail existing checks or lose required inputs. The original owner's boundary explicitly allowed ordinary production writes; it cannot establish the stronger preservation requested now.

A source-level private topology exists in `testing/run_pion_local_t_matched_cohort_gate.py`: `prepare_private_output`, `input_plan`/`stage_input_files`, `diamond_selections`, and `audit_default_output_writer`. This establishes a bounded implementation pattern, **not farm compatibility or permission to run that failed acquisition**. `findEffectiveCharge.py` at this HEAD directs its CSV into LTANAPATH-local OUTPUT; `shift_prep` can update analysed ROOT containers, which therefore require private copies. Diamond selection and both epsilon prepasses require more than the four lowe data/dummy files.

**Before implementing the repair**, trace the actual unchanged `-d` path and make a short source-backed read/write map in this contract's outcome section. Show how owner-private OUTPUT can retain all required inputs, default/Plot path resolution, both prepasses, and independent clean-source collection without any write-through to accepted outputs or original mutable inputs. Identify every reached external writer/link repair and its guard. Test the topology with real temporary Linux/WSL symlinks and harmless files; no production imports or execution.

If that map cannot close within the allowlist below, **STOP `BLOCKED` with the first exact writer/path/check**. Do not implement a partial source-gate relaxation and call the owner repaired. Do not bypass the frozen helper's check, patch it, or claim that merely redirecting the four comparison files preserves other scientific outputs. No second contract or speculative farm command follows this stop.

## Allowed implementation after that gate

1. **Reviewed source successor.** Retain the historical baseline identity as history/ancestry evidence. Add the exact starting HEAD above as the explicit reviewed successor anchor for this owner's future committed-range check and tracked profile. Verify historical ancestry. From that anchor to the requested exact pushed SHA, allow only the existing owner, its test, its profile and `docs/memory/`; reject every other committed source change. Thus the 13 paths are accepted at their exact reviewed checkpoint bytes, not placed on an open future-edit allowlist. Keep branch/HEAD/origin equality, dirty-source rejection, detached source preflight, exact committed SOURCE_FILES hashes, profile equality, and collector source checks. The effective collector profile still pins the exact requested SHA with `allowed_committed_files=[]`. Never rewrite historical A/B artifacts or C's archival source pin.

2. **Private runtime outputs in this owner.** Keep `--outdir` as the ordinary configured authority/receipt location and `--output` as the existing fresh ZIP destination. Create a unique private OUTPUT tree only inside the established owner-created disposable sibling area, with an exclusive worktree OUTPUT link. Add a narrowly named owner-local private-topology check in place of calling the shared-output link helper for this mode. Explicitly verify ownership, containment, absence of preexisting paths including dangling links, no overlap with shared output, literal/resolved link identity, and runtime OUTPATH resolution. The frozen shared helper remains unchanged and retains its original regression tests. This is an explicit stronger owner invariant, not an ignored exception.

3. Preserve the copied-ltsep mechanism: **only copied LTANAPATH changes**; installed package/config and stable external roles, including VOLATILEPATH, remain unchanged. Prove launcher/setup/link caller agreement and reached default/Plot writer paths, including the actual effective-charge CSV expression. Scope TMPDIR and writable runtime caches to owner-private directories. External symlink preflight must prove no repair is needed; existing external setup directories must not need mutation. Unresolved caller/dependency/path behavior is a pre-child hard stop, not a reason for a new resolver framework or smoke farm run.

4. Stage only the source-derived consumed merged ROOT/histogram inputs as independent byte-verified private copies, including both epsilon prepasses and the exact Diamond candidate namespace/tie behavior. Preserve ordinary setting availability and selection; do not silently restrict normal all-low debug processing. No hard links or mutable-input symlinks. Record original and runtime-container identities separately because frozen derived-branch writes can change private ROOT bytes. Inspect the established staging logic as reference; do not invoke/import the matched-cohort orchestration, its observer hooks, event capture, resolver instrumentation, recursive shared inventory or expensive external byte-hashing gate. Inventory existence-only skims by the actual existence predicate, never whole-file hashes. Expose the exact copy list and total byte count before copying; missing/colliding inputs fail before expensive work. No whole-directory hash/copy operation.

5. The four preexistence/freshness checks apply to the **private generation targets**. Independently record and recheck the four accepted ordinary artifacts' identities and links without moving or rewriting them. All ordinary scientific writes must be excluded by the traced private topology, not assumed safe because those four hashes match. Retain exclusive receipt/log/ZIP destination checks. Enable only the existing A/B observer for exactly one unchanged `./run_Prod_Analysis.sh -d 4p4 2p74` child in its disposable cwd; environment conflicts still fail closed. Preserve current all-low completion markers and no full high-epsilon execution.

6. Verify private artifacts with the existing payload, frozen-parent, source, PDF, availability/negative-evidence and freshness gates. Collect from the private artifact root using the unchanged collector imported from a **separate clean detached source** at the exact pushed SHA. Keep private files alive through verification and collection; do not exit their lifetime before the collector reads them. Preserve exact four-artifact ZIP inventory, CRC, hashes, manifest, log/summary/status semantics. Finish bounded cleanup and all preservation checks before success delivery. Failure retains diagnostic receipts/logs and runs independent preservation checks; no failed evidence is delivered as accepted. Do not copy runtime scientific outputs back to ordinary output.

## Strict edit allowlist and frozen interfaces

- `testing/run_pion_t_template_comparison_gate.py`
- `testing/test_run_pion_t_template_comparison_gate.py`
- `testing/pion_t_template_comparison_gate_profile.json`
- `docs/memory/CURRENT.md`: shortest material blocker/NEXT correction; keep under 16 KiB.
- `docs/memory/manifest.json`: generated only.
- This single placed contract: append a concise implementation/outcome record, including any hard stop.

All other files are frozen, including the 13 checkpoint paths, shared owners/helpers, collector, launcher, link setup, installed/copied Python code, A/B producer, C comparison and tests, science/configuration and accepted evidence. No new CLI, observer, renderer, event acquisition, external path configuration or general framework. No changes to A/B/C mathematics, seven-strength grid, fixed K_B, SUMW2, fits, alignment, templates, weights, priors, cuts, normalization, binning, subtraction, yields/errors or presentation. Method A stays detached/non-production; Method B diagnostic and numerically excluded. Historical October 20 acceptance and the nine failed acquisition attempts retain their scopes.

## Required local evidence and acceptance

Run only deterministic source/pure-Python tests; actual ROOT/PyROOT, main.py, production launcher, farm and PDF rendering are forbidden here.

- Reproduce the original 13-path rejection and ordinary four-artifact collision. Show the anchor change accepts exactly the reviewed history and rejects a subsequent edit to any of those 13 paths, an unrelated committed path, wrong branch/SHA/origin, dirty source and tampered profile. Keep committed-source hash and final collector rejection tests.
- Real Linux/WSL topology fixtures: shared accepted artifacts plus representative ordinary outputs remain byte/link-identical; child inputs remain visible; required private copies have distinct inodes; harmless writes representing CSV, PDFs and derived ROOT updates land privately. Cover default/Plot callers, lexical/canonical aliases, both prepasses and Diamond selection parity. Wrong/outside/dangling/existing links, missing inputs, path escapes and external writers stop before the mocked child.
- Mocked end-to-end existing owner: exactly one child, collection from private output while it exists, independent clean-source checks, four fresh objects/strict ZIP, cleanup before delivery. Inject source/path/copy/child/collector/ZIP/cleanup/preservation failures and assert no success delivery, no accepted-file mutation and retained failure evidence. Freshness, partial-unavailable scientific evidence and six-page PDF structure gates remain strict. Never accept stale artifacts or manipulate runtime timestamps to pass.
- Run the affected owner suite and existing debug-link, debug-launcher, collector and template-comparison regression suites where supported. Report exact commands, interpreter/OS, counts, skips and failures; a missing required test is a hard stop, not a PASS. Compile changed Python without leaving caches; run `git diff --check`.
- Regenerate then check memory manifest and run ordinary memory health. Report CURRENT/MEMORY/handoff bytes and classify warnings. Startup audit: manifest PASS; CURRENT 16221 bytes (nonblocking soft-size warning), MEMORY 29898, handoff 323. Replace short CURRENT text rather than grow it. Do not edit supplemental native memory.

Return the complete actual tracked diff plus full new-file additions without staging, in one temporary repository-root review bundle if needed. Show unchanged unrelated worktree/index/refs and frozen-file identities. Include a concise source-only feasibility verdict, unresolved installed-environment facts and measured copy scope; tests never establish farm readiness. Supply a non-executed scoped user Git handoff only after successful local completion, excluding temporary review files.

NEXT is independent ChatGPT actual-diff review, then user-controlled publication and pushed-state review. Any later farm operation requires a separate evidence-based readiness decision; **this contract supplies none and no farm command**. The eventual scientific target remains fresh Center/Left t1/t2/t3 A/B followed by unchanged offline C histogram JSON/PDF. Matching source HEAD alone proves neither same execution nor population. Signed-event application, final Lambda yields/variance, physical bias truth, C selection and production adoption remain NOT VERIFIED.

## Local implementation outcome â€” October 25, 2026

`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; independent actual-diff review
is NEXT. Starting branch/HEAD/local origin matched the exact checkpoint above.
Only this placed contract was initially untracked; staged content was empty.
The initial manifest discrepancy was exactly that new contract, with no changed
or missing existing entry. No farm, ROOT/PyROOT, actual main.py, production child
or PDF rendering was executed. Test children/backends were synthetic mocks.

### First-gate source read/write map

The source map and two real temporary Linux experiments preceded owner edits.
They used harmless files and pure path/filename expressions only; no production
modules or matched-cohort orchestration were imported or executed.

| Reached source/path | Read/write boundary and owner guard |
| --- | --- |
| run_Prod_Analysis.sh normal -d | Keeps both raw epsilon prepasses and ordinary all-low selection, then exits before full high; no -c merge/error-log writer. TMPDIR holds the paired captures privately. Command and completion markers are unchanged. |
| launcher mkdir and set_SymLinks.sh | Private OUTPUT is created exclusively before the child. Stable skim and shared background directories must already exist, so their mkdir calls are no-ops. Frozen external_symlink_preflight rejects missing/broken external SIMC OUTPUTS/input/worksim links and any background worksim repair. External links/directories are rechecked before the child and at preservation gates. Local simc/log/src-data/background links belong to the disposable tree; the reached source consumes their external referents read-only. |
| findEffectiveCharge.py -> getDataTable.py -> utility.data_to_csv | The actual AST expression LTANAPATH/OUTPUT/Analysis/<ANATYPE>LT is checked, independently of no-mode Root.OUTPATH. Its CSV lands in private OUTPUT. The other three no-mode helpers use caller-supplied CSV or return scalars, not Root.OUTPATH. |
| main.py, diamond.py and Plot_LTSep consumers | Source-only import traversal finds 24 caller/cwd/mode combinations. Baseline/copied package, canonical caller, stable external roles and private Plot output are checked before copying/child execution. Diamond's low/mid/high candidate namespace, case behavior, exclusions, fallback and first equal-length tie are compared after staging. |
| shift_prep.py -> shift_MM.add_derived_branches_to_file | Data/dummy analysed ROOT containers are independent verified copies, including both prepasses and all ordinarily available settings. MM_shift/t_shift writes affect private containers; original and runtime identities are retained separately. Original file and link identities are independently rechecked, including failed-copy/child cases. |
| set_sig_fortran.py, find_bins.py, model activation, tshift compiler, iteration stamp | Model rewrites, interval/metadata publication, param_active copy, compiler outputs and stamp stay under copied LTANAPATH in the disposable worktree. Analysis PDFs/ROOT/JSON and cache_transfer stay under private OUTPATH. Writable Python/XDG/Matplotlib caches and temporary files are owner-private. High-epsilon separation is not reached; the main.py jput block is dormant text. |
| external efficiency/runlists/PID logs, background/Sigma0 templates and skims | Resolver roles and configured sources remain unchanged. No external link repair, recursive inventory or whole-directory byte hash/copy is introduced. Skims are inventoried with the actual exists predicate, including absence; their bytes are never opened or hashed. |
| independent collector | The effective profile pins the requested exact SHA and allows no committed source differences. The unchanged collector is imported from a separate clean detached source, while private inputs/artifacts remain alive. Four-artifact ZIP checks precede both bounded cleanups and independent preservation; only then are success companions delivered. Ordinary scientific outputs are never copied back. Receipt/log destinations remain the existing exclusive ordinary/transfer roles, so the whole output root is not claimed immutable. |

The first Linux topology experiment copied 21 files / 975 fixture bytes with
distinct inodes and preserved shared bytes and symlinks after private CSV/PDF/
derived-ROOT stand-in writes. The second executed frozen pure Diamond selectors
and preserved 12 low/mid/high, fallback and tie selections. The implemented
owner fixture copies 21 files / 1017 bytes. These are synthetic measurements:
actual farm input inventory/bytes, installed dependency roots, real Root resolver
behavior, native ROOT updates, performance and farm preservation remain NOT VERIFIED.
The owner records and emits its exact source-derived copy list and byte total
before any copy; it fails on missing inputs, collision or changed plan.

### Repair and local evidence

Historical BASE_HEAD is retained as ancestry evidence; REVIEWED_HEAD and the
tracked profile pin the exact starting checkpoint. Only owner/test/profile and
memory are allowed after it; the 13 checkpoint paths are not future allowances.
Tests execute checkpoint owner code to reproduce unexpected_committed_source
and comparison_artifact_preexists with four accepted targets intact. New tests
reject each of the 13 future edits and an unrelated path, as well as branch,
SHA/origin, dirty source, profile, path/copy/child/collector/ZIP/cleanup and
independent preservation failures. Existing source-hash, freshness, partial-
unavailable, frozen-parent and six-page structure validators remain unchanged.
No artifact was retimed to pass; mock-child fixtures make genuinely fresh writes.
One initial positive CSV fixture had inconsistent LTANAPATH/OUTPATH; correcting
that fixture retained the guard and its strict negative test.

Local commands (WSL Ubuntu, Linux Python 3.12.3):

```text
python3 -B -m unittest testing.test_run_pion_t_template_comparison_gate -v
python3 -B -m unittest testing.test_run_e8_2_left_lowe_scientific_audit_gate testing.test_run_prod_analysis_debug_left_low testing.test_collect_pion_hgcer_validation_bundle -v
```

Final affected suite: 29 passed. Frozen regression suites: 27 + 15 + 28 passed.
Windows Python 3.12.10:

```text
python -B -m unittest testing.test_pion_t_template_comparison -v
```

23 passed; total 122 distinct tests, zero skips/failures in final results.
Both changed Python files passed py_compile with bytecode directed to temporary
storage. git diff --check passed. AST comparison against starting HEAD confirms
only resolved_profile, preflight and execute_gate changed among existing owner
functions; scientific, payload, source-hash, completion and ZIP validators are
unchanged. All frozen/unrelated tracked file bytes and refs are checked against
the startup snapshot. The staged index remains empty as at startup; raw index
bytes changed during Git inspection, so byte-identical index metadata is not
claimed and no restoration/staging was attempted.

CURRENT is minimally updated; generated manifest and ordinary health results
accompany the review bundle. CURRENT's existing soft-size warning is nonblocking
and deferred to the next checkpoint/milestone audit. This is source/local proof,
not farm readiness, same-run evidence, scientific acceptance or production
adoption. Historical acceptance and the nine acquisition failures retain scope.


### Narrow optional-SIMC compatibility follow-up - October 25, 2026

Independent review reported 99 passing WSL tests but reproduced an uncovered
copy_plan regression: data/dummy presence made every archived SIMC mandatory.
The new WSL regression first failed with
required_simc_input_missing:Prod_Coin_Q4p4W2p74right_lowe.root.
SOURCE VERIFIED: frozen shift_prep returns an unavailable result when archived
SIMC is absent; main.py filters those settings after Diamond processing.

Only copy_plan changed in executable owner code (five added lines). Missing
archived ROOT now follows that skip for optional Right/lowe and high-epsilon
settings. Center/Left lowe data/dummy, archived ROOT and histogram stay mandatory.
A present optional lowe ROOT still requires its histogram; an existing invalid
ROOT is rejected rather than hidden. The unchanged complete Diamond glob still
stages raw ROOT candidates, including skipped settings and both epsilon scopes.
No source checks, private confinement, accepted-artifact guards, scientific
source, model/fit mathematics, profile, launcher or producer were changed.

Three new tests reproduce the failure, execute only frozen availability/filter
AST without ROOT, check shared/private selection parity, preserve the complete
raw namespace, cover each optional high setting with other high inputs present,
and retain mandatory/malformed-input failures. Final WSL Ubuntu Python 3.12.3:

```text
python3 -B -m unittest testing.test_run_pion_t_template_comparison_gate testing.test_run_e8_2_left_lowe_scientific_audit_gate testing.test_run_prod_analysis_debug_left_low testing.test_collect_pion_hgcer_validation_bundle -q
python3 -B -m unittest testing.test_pion_t_template_comparison -q
```

102 + 23 = 125 tests passed, zero skips/final failures. Both affected Python
files passed py_compile with temporary bytecode; whitespace checks passed.
The prior review bundle and all pre-existing unrelated edits are byte-preserved.
Refs and empty staged contents are preserved; raw index metadata changed during
read-only Git inspection. CURRENT already names the same review/readiness gate
and is left unchanged. Manifest is regenerated for this appended outcome.

Local/source verdict: PASS for this narrow compatibility repair; updated diff
still needs independent review. Farm compatibility/readiness is NOT VERIFIED.
After actual-diff acceptance and user publication/pushed-state synchronization,
the smallest remaining readiness evidence is the intended farm interpreter and
actually imported ltsep/package/config identities plus real baseline/copied
caller/cwd path observations (default CSV and Plot outputs), stable external
roles/no-repair links and existing setup directories, and the bounded consumed
input/Diamond namespace, copy-byte and destination/capacity observations. Include
an available high-epsilon prepass population and mandatory Center/Left lowe
inputs; synthetic fixtures cannot supply those facts. Applicable existing
receipts can supply them; no new support contract, acquisition wrapper, science
smoke cycle or whole-directory hash gate is proposed. Native ROOT writes,
actual procedure/comparison PDFs and runtime preservation remain evidence of a
later separately justified user operation, not claims from these local tests.
No farm operation, publication, or additional contract was executed/created.
