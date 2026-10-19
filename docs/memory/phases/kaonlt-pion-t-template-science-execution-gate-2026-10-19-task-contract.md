# KaonLT — One controlled canonical-t template A/B farm evidence gate

**Task class:** narrow farm-bound owner implementation, **not** a new scientific phase.  
**Observed source:** GitHub `test` at `554ad84630ad77920ce17954d99f73df8ed2e6ac` (2026-10-19).  
**Work state:** `ACTIVE` contract; existing diagnostic `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`.  
**Sole scientific NEXT:** one actual Center/lowe and Left/lowe fixed-alignment A/B comparison, with independent returned evidence review.

## 1. Scientific target and smallest scope

Produce **fresh, farm-generated**, detached diagnostic evidence using the *already published* `src/cuts/pion_t_template_comparison.py` and its opt-in `LT_ANALYSIS_PION_T_TEMPLATE_COMPARISON=1`:

- `Q4p4W2p74 / Center / lowe`, canonical `t1,t2,t3`, primary scientific test; and
- `Q4p4W2p74 / Left / lowe`, canonical `t1,t2,t3`, separate, acceptance-limited control.

The sole controlled fit variable is the three pion-component SIMC template shapes: **A** is the frozen authoritative parent with its selected templates; **B** calls the *existing* fitter using genuine raw-weighted canonical-t templates and fixed parent alignment, without modifying A. Evaluate model construction and residuals **before** pion subtraction. The three components are `pi_n`, `pi_delta`, and `pi_sidis`; protected signals and data are unchanged. A failed, sparse, or unavailable B is an informative explicitly rejected scientific outcome, not permission to fall back, reweight, retune, fit children, or replace A.

The producer already writes exactly these four evidence files directly to the runtime `OUTPATH`:

```text
Center_kaon_Q4p4W2p74_lowe_pion-t-template-comparison.json
Center_kaon_Q4p4W2p74_lowe_pion-t-template-comparison.pdf
Left_kaon_Q4p4W2p74_lowe_pion-t-template-comparison.json
Left_kaon_Q4p4W2p74_lowe_pion-t-template-comparison.pdf
```

**No new diagnostic, physics adjustment, wide owner/collector program, or extra evidence package.** The ZIP may also contain the collector's provenance files and manifest, and the owner may return a separate run log/receipt; it must not demand unrelated E.8 or F.6 material as a condition for this A/B test.

## 2. Authority, source and local-state gate

Before editing, obtain live `origin/test`, local `test` HEAD, complete staged/unstaged/untracked worktree state, existing pending task contracts and temporary diffs. Read **in full, in order**: `AGENTS.md`, `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, `docs/memory/USER.md`; then the directly relevant `MAINTENANCE.md`, `CODEX.md`, `COMMUNICATION.md`, `TOOLS.md`, the current baseline-pion investigation and the published controlled-comparison contract. Preserve all unrelated bytes, index state, prior task contracts and review files. Stop `BLOCKED` if the expected baseline is incompatible with current source, if an intended new file already exists, or if another pending change overlaps this task. The observed SHA is an audit anchor, not authority to overwrite newer work.

The earlier proposed **collector-expansion** contract (`kaonlt-pion-t-template-four-artifact-profile-2026-10-19-task-contract.md`) is **superseded and must not be implemented** here. A harmless pre-existing untracked copy may be left untouched; never stage it accidentally or silently resolve an overlapping implementation.

## 3. Frozen science and run boundary

No changes to: `run_Prod_Analysis.sh`, `set_SymLinks.sh`, `src/main.py`, `src/cuts/rand_sub.py`, `src/cuts/pion_t_template_comparison.py`, `src/utility/pion_component_shapes.py`, `src/cuts/pion_component_fits.py`, `src/cuts/pion_t_bin_parents.py`, any production source/configuration, `testing/collect_pion_hgcer_validation_bundle.py`, or previously accepted owner/helper code. Preserve random/dummy/proton/pion ordering, all frozen data and MC weights, PID, templates, cuts, priors, cross sections, histogram coordinates, parent fits, children, yields, and `no_empirical_residual`. Method A stays detached and Method B numerically excluded. Neither A/B status nor a successful optimizer is scientific acceptance.

Codex performs local implementation and deterministic tests **only**. It must not execute the farm, launch the actual production script, commit/push, delete or relocate farm outputs, mutate the ordinary checkout, patch installed `ltsep`, or claim ROOT/PyROOT or rendered-PDF validation.

## 4. Existing verified operations: reuse rather than reinvent

Relevant published source at observed HEAD:

- `run_Prod_Analysis.sh -d 4p4 2p74` skips `git clean`, runs ordinary low-epsilon available settings, stops before full high epsilon, **but still** invokes `set_SymLinks.sh`, `set_sig_fortran.py`, and normal `main.py`; the runner is **not** harmless in the ordinary checkout.
- `testing/run_e8_4_fix5_canonical_five_plot_gate.py` supplies `owned_worktree`, `prepare_runtime_overlay`, `probe_paths`, `external_symlink_preflight`, checkout/installed-`ltsep` preservation checks, detached collector import and gate-status/delivery utilities. This owner has other E.8-specific physics/provenance gates that **must not** be imported into the new scientific test.
- `testing/run_e8_2_left_lowe_scientific_audit_gate.py::prepare_debug_output_link` handles the early `-d` output-link requirement inside a disposable analysis worktree. Reuse this narrow, proven path helper, **not** that owner's historical Left-only completion, page or full-analysis gates.
- `testing/collect_pion_hgcer_validation_bundle.py` already accepts `generic_artifacts` with **nonempty global declarations** and **empty `settings` artifact declarations**. Its profile schema requires the original five-setting list; the four explicit comparison artifacts can be declared as **global entries** with literal basenames. This avoids *any* collector code change. The five empty setting metadata folders are schema machinery, **not** scientific results or five analysis settings.
- The existing generic collector validates source ancestry and committed range, checks JSON syntax, records hashes and manifests and fails on absent required files; it **does not** establish comparison-schema quality or PDF legibility. The owner must check those separately.

Recheck the actual code in Codex rather than relying solely on this inventory. **One tracked owner CLI** must own source preflight → disposable worktree/copy-`ltsep` path isolation → one debug analysis invocation → completion/diagnostic identity/freshness checks → independent-source generic collection → ZIP/manifest verification → returnable receipt. There must be **no interactive multi-command farm workaround**. No need to add new generalized wrappers.

## 5. Allowlisted changes

**New, first verify absent:**

1. `testing/run_pion_t_template_comparison_gate.py` — the sole narrowly owned execution/verification/collection CLI.
2. `testing/test_run_pion_t_template_comparison_gate.py` — deterministic, dependency-free focused owner/profile/negative-path tests.
3. `testing/pion_t_template_comparison_gate_profile.json` — a declarative `generic_artifacts` profile of exactly four required objects; no global or per-setting extras.
4. `docs/memory/phases/kaonlt-pion-t-template-science-execution-gate-2026-10-19-task-contract.md` — **this transferred contract**.

**Existing, only if required for a material source-state update:**

5. `docs/memory/CURRENT.md` — at most a short, accurate, push-stable mention of the owner candidate without changing the sole scientific NEXT or claiming runtime readiness; no standalone reconciliation.
6. `docs/memory/manifest.json` — regenerate only if versionable memory changes (the new contract itself requires manifest update).

No alteration of `MEMORY.md`, handoff, historical incident, roadmap or other source is authorized. If an already published owner has a concrete invariant that cannot be satisfied without changing a frozen path, stop `BLOCKED` and report that exact blocker; do not silently add a new dependency gate.

## 6. Exact profile semantics

The profile must pass the **existing, unchanged** collector's v4 schema:

```json
{
  "schema_version": "pion_hgcer_validation_bundle_profile/v4",
  "validation_profile": "pion_t_template_comparison_center_left_lowe/v1",
  "collection_mode": "generic_artifacts",
  "settings": [
    {"phi": "Left", "epsilon": "lowe"},
    {"phi": "Left", "epsilon": "highe"},
    {"phi": "Center", "epsilon": "lowe"},
    {"phi": "Center", "epsilon": "highe"},
    {"phi": "Right", "epsilon": "highe"}
  ],
  "artifacts": {
    "global": [
      {"key": "center_comparison_json", "basename_template": "Center_kaon_Q4p4W2p74_lowe_pion-t-template-comparison.json", "kind": "json", "required": true},
      {"key": "center_comparison_pdf", "basename_template": "Center_kaon_Q4p4W2p74_lowe_pion-t-template-comparison.pdf", "kind": "file", "required": true},
      {"key": "left_comparison_json", "basename_template": "Left_kaon_Q4p4W2p74_lowe_pion-t-template-comparison.json", "kind": "json", "required": true},
      {"key": "left_comparison_pdf", "basename_template": "Left_kaon_Q4p4W2p74_lowe_pion-t-template-comparison.pdf", "kind": "file", "required": true}
    ],
    "settings": []
  },
  "source_identity": {
    "required_analysis_commit": "554ad84630ad77920ce17954d99f73df8ed2e6ac",
    "allowed_committed_files": [],
    "allowed_non_analysis_path_prefixes": ["docs/memory/"]
  }
}
```

At the eventual single-run gate, accept **only the exact pushed SHA** as `--source-commit`. Like the existing reviewed E.8 owners, the new owner verifies the checked-in profile's frozen baseline SHA and constructs one **temporary** effective profile with `required_analysis_commit` set to the just-verified exact pushed SHA, `allowed_committed_files=[]`; do not edit the tracked profile, broaden allowed source changes, change the collector, or substitute stale base identity. The detached collection worktree must be checked at that exact SHA before packaging. Assert exactly four evidence-file members plus the collector's known standard metadata/directory members.

## 7. One owner invocation, isolation, and fail-closed order

Implement a CLI accepting `--source-commit` plus the minimal existing-owner input pattern (`--repo`, `--outdir`, `--output`), all validated against the configured JLab owner/profile boundary. Do **not** bake private workstation/farm paths into new public memory. An eventual user-controlled `tcsh` command is supplied **only after** actual-diff acceptance, publication, pushed-state review and farm readiness review.

Before invoking the child:

1. Validate requested kinematics are fixed to `Q4p4W2p74`; profile is exactly as contracted; source commit is full 40-hex, local `test` and `origin/test` equal the supplied pushed commit; no unexpected changed or untracked source; source identity and collector source checks PASS from a separate owned, clean, detached worktree. **Do not** interpret a generic `complete=true` as scientific validity.
2. Snapshot ordinary checkout/source status and known farm-local modifications. Do not clean, reset, stash, check out over or restore the ordinary checkout. Record the preexisting identities of four diagnostic target files and log/receipt/ZIP destinations. Require unique new ZIP, log, receipt and companion destinations. Do not erase or silently overwrite previously accepted detached comparison evidence; if a target already exists without explicit authorized overwrite semantics, fail `comparison_artifact_preexists` rather than deleting it. Ordinary production outputs may be written by the established `-d` analysis; explicitly record this limitation—**not** an immutable whole-output-root guarantee.
3. Create one uniquely owner-controlled disposable analysis worktree at exactly the source commit, using the existing bounded helper. Create and verify its copied `ltsep` configuration with **only copied LTANAPATH changed**, and prove identical non-`LTANAPATH`/`OUTPATH` path roles. Prepare the existing `-d` output link **inside that worktree** after path verification; separately verify external SIMC/background links will require no repair. No `git clean` in ordinary checkout, no installed-package edit or symlink repair outside the disposable worktree.
4. Set the diagnostic environment variable to exactly `1` **only for this owner's child analysis**; remove or reject inherited conflicting `LT_ANALYSIS_DEBUG_LEFT_LOW`, unpaired-binning/prepass debug overrides, or other injected opt-ins that would change this experiment. Do not change global shell/ordinary environment. Run the unchanged `./run_Prod_Analysis.sh -d 4p4 2p74` **exactly once** from the disposable worktree with its verified environment, capturing a fresh full child log.
5. Require exit 0 and the **current** all-low completion markers (`All available low-epsilon settings completed in debug mode.` and `Full high-epsilon processing is intentionally skipped in -d debug mode.`), with no high-epsilon full-analysis completion. Do **not** reuse the older `Left / lowe debug analysis completed.` gate, because that is not this launcher mode.
6. Verify ordinary checkout and installed `ltsep` preservation **even on any child, path, renderer, artifact or cleanup failure**; never copy back mutations, never mark a child-success receipt as owner-success. The owned temporary-worktree cleanup remains bounded to owner-created paths. Preserve failure receipts and logs for diagnosis; no accepted ZIP delivery after a failed gate.

Do not add an E.8 F.4 equivalence, unrelated page-manifest, Method-A/Method-B, absolute SIMC normalization or high-epsilon gate: those are different questions.

## 8. Structural evidence acceptance **before collection**

For **each** Center and Left artifact pair, require:

- Exact expected file basename, regular nonempty newly generated file, fresh timestamp/changed identity relative to pre-run snapshot and exact computed SHA-256; strict JSON with no nonfinite values; `%PDF-` signature, readable page structure and **exactly six pages** (two pages each for pion- and kaon-side model/residual displays for `t1,t2,t3`), using a verified existing PDF-reading backend or fail closed if none can verify the count.
- `schema_version == "pion_t_template_comparison/v1"`; `diagnostic_only == true`, `production_applied == false`, `scientifically_accepted == false`; `source_HEAD == --source-commit`; `enable_environment_variable == "LT_ANALYSIS_PION_T_TEMPLATE_COMPARISON"`; `phi_setting` equals the setting; correct `Q2`, `W`, `EPSSET`, `ParticleType`, and four canonical edges/t indices. The three parent records must be unique `t_bin_index` 0,1,2 with frozen A and unchanged-authority checks, explicit fixed-alignment test mode, frozen input identities, all three component provenance records, and B result or explicit per-t rejection reasons. Guard mixed/duplicate/missing settings/t bins.
- All required `source_files` identities are available, present, SHA-256 verified against **the exact committed source bytes** (not only against another file created by the child). `source_identity_rejection_reasons == []`. Independently check that the JSON's `pdf.status` and `pdf.sha256` agree with the freshly verified PDF. Source identities and A-parent preservation must never be silently downgraded to nonblocking warnings.
- Preserve distinction between producer **availability** and scientific model quality. A `partial_unavailable` result with explicit component/fit rejection reasons and intact provenance/authoritative invariants is **valid negative diagnostic evidence**; do not force `available` or invent missing candidate shapes. Treat structural failure, missing expected parent, or unexplained unavailable state as an owner failure.

A readable six-page PDF is **not** proof of visual legibility or numerical fit adequacy. Those are reviewed by ChatGPT after a successful returned evidence package.

The owner writes an atomic, independently inspectable gate-status receipt recording precise failure stage/reason, HEAD, source/preflight, output identities, analysis completion, preservation and collection status, and a concise run summary/log identity. A failed invariant stops forward motion immediately; no automatic rerun, copy of failed artifacts as accepted evidence, or scientific interpretation.

## 9. Package and verify one narrow ZIP

Only after source + child + structural + provenance checks PASS, use the existing generic collector from the exact detached pushed source with the runtime-instantiated profile. Package **only the four required comparison outputs** as declared global artifacts, plus its standard `manifest.json`, `source_state.txt`, `source_checks.txt` and schema-required empty directory entries. Keep source checks and field-level semantic checks distinct.

Check ZIP CRC, exact allowlisted archive members, manifest `complete == true`, empty `errors`, `git_head == required_analysis_commit == --source-commit`, correct `validation_profile`, no extra scientific files, and independently match each member's bytes/length/SHA-256 to the freshly verified files. Recheck ordinary checkout and installed-package preservation after collection and before reporting success. Provide one fresh ZIP, plus uniquely named log and owner status/run-summary receipts via existing safe companion-delivery helpers. Do not use `testing/package_pion_hgcer_validation_bundle.tcsh` to launch a new analysis: it is bundle-only and cannot own this workflow.

## 10. Local deterministic positive/negative/regression tests

Implement focused standard-library-only tests with synthetic artifacts, mock subprocesses, isolated fake worktrees and a mocked collector module. No ROOT, farm data, network or production launch during tests. At minimum cover:

1. Exact four literal basenames and profile shape with empty setting-artifacts and five schema-required metadata settings; profile accepted by **unchanged** v4 collector; global-only archiving with correct four hashes; no collector edits.
2. Source SHA/branch/remote mismatch; dirty source; edited required source; stale or unexpected committed content; untrusted environment overrides; missing/unsafe external link; wrong or preexisting OUTPUT link; inherited paths; `ltsep` copy/import/path mismatch.
3. Exactly **one** real-subprocess *mock* child invocation with diagnostic env `1`, `-d` parameters and disposable `cwd`; current all-low markers, no high epsilon; error path for exit 1 or missing markers.
4. Missing, stale, empty, malformed JSON, invalid status/setting, incorrect t inventory, wrong parent/authority/normalization/provenance/committed source hash, mixed setting, invalid PDF hash/page count and six-page success. If candidate is unavailable with explicit valid rejection reasons and all invariants pass, accept it as negative evidence and do not relabel it successful.
5. Failure at each gate prevents downstream collection and success delivery; test retention of failure receipt and preservation after child failure and after source/ZIP failure. Unexpected old target files or output collision cannot be deleted or overwritten. No need to rerun on failed gates.
6. ZIP manifest, identity, hash, completeness and exact contents; no high-epsilon/Right diagnostic accidentally included; sidecar log/receipt bounded delivery.
7. Regression tests for `testing/test_run_e8_2_left_lowe_scientific_audit_gate.py`, `testing/test_run_prod_analysis_debug_left_low.py`, `testing/test_collect_pion_hgcer_validation_bundle.py`, and `testing/test_pion_t_template_comparison.py` where runnable locally. Do not rework those existing tests. If a specific external dependency is missing, report its test as **unexecuted**, not passed or a reason to invent general environment repairs; a required owner test that cannot run is a bounded hard stop.

Run compilation/syntax checks of changed Python files, focused unit tests, `git diff --check`, manifest generation **then** manifest check, and `<PYTHON> -B tools/check_memory_health.py --root .`. Discover `<PYTHON>` locally. Report warning byte counts for CURRENT, MEMORY, CURRENT_HANDOFF and classify warnings as blocking/nonblocking; don't enter a cosmetic memory-reconciliation cycle.

## 11. Actual-diff review, hard stops and delivery

Codex must not commit, push or run Jefferson Lab farm. Return:

- exact local/remote identity and worktree inventory before edits;
- exact paths changed, actual complete tracked diff **and complete new-file additions**; if lengthy create a single temporary repository-root `kaonlt_review.diff` including both `git diff -- <tracked>` and `git diff --no-index /dev/null <new>` for each new file, without staging merely for review;
- focused tests and missing tests, source/identity/cleanup gates and memory-health receipt;
- an exact **non-executed** `git add` / `git commit` / `git push origin test` handoff only for the task's accepted paths, excluding prior superseded contracts and any temporary diff;
- the single next step: independent ChatGPT **actual-diff** review, then user publication, ChatGPT pushed-state synchronization, **one** user-run scientific gate, independent validation and physics interpretation.

**Hard stops:** Any unexpected source change or output path, absent contract-required file, failure to prove the existing owner helper safe for this use, destructive external write, unresolved provenance/identity, missing required local test, owner/checker/manifest/PDF/artifact failure, incompatible pending local change, or need to expand frozen source requires a precise `BLOCKED` report. Do not convert the stop into speculative infrastructure, a second contract, an alternate manual farm command, broad reconciliation, or weakening a safety check. Preserve older closed scopes and the sole CURRENT scientific NEXT.
