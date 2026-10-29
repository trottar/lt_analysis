# KaonLT — Q4p4W2p74 High-Epsilon Three-Model Evidence Acquisition

**Task class:** narrowly scoped source-changing farm acquisition; **state:** `ACTIVE` (contract, not implementation approval).  
**Repository:** `trottar/lt_analysis`, branch `test`.  
**Source-audited starting HEAD (2026-10-29):** `38307c486326faf159cf300fb7c60d6c3758a9f7`.  
**Scientific output:** fresh Center/highe, Left/highe, Right/highe canonical-t and child-(t,phi) Global/Local-t/Local-(t,phi) control/transfer evidence at `Q4p4W2p74` from one controlled JLab execution.  
**Runtime boundary:** ONLY the user runs the Jefferson Lab farm, after independent actual-diff review, user commit/push, pushed-state synchronization and readiness approval. This contract is not a farm command.

## 1. Startup and established state

Before edits: verify actual `test` branch, HEAD, `origin/test`, and **full worktree state**. Read the tracked root `AGENTS.md`, then `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, and `docs/memory/USER.md`, in exactly that order. Read `docs/memory/MAINTENANCE.md`, the relevant baseline-pion investigation, the current comparison owner/profile/tests, and the October 22–24 farm incident before deciding on source changes. Follow any CURRENT direct references needed for this task.

Do **not** silently treat the audited SHA as the current HEAD; stop and re-anchor if it has moved. Do not overwrite untracked/dirty work. Tracked CURRENT's pre-October-28 NEXT is materially stale relative to supplied October 28–29 low-epsilon evidence; document this drift and make only the smallest accurate memory update as part of this milestone. Do not reopen closed phases or change unrelated scientific status.

Accepted low-epsilon Work evidence consists of the October 28–29 Q4p4W2p74 Center/lowe and Left/lowe Global/Local-t/Local-(t,phi) numerical comparisons and the corrected 60-page visual review. It favors Local-t **as a candidate**, not an adopted production change. Existing historical high-epsilon Phase-C/Method-B plots are **not** the requested high-epsilon three-model evidence. Method A remains detached/non-production; Method B remains independent and excluded numerically.

## 2. Exact audited source restrictions

The current restriction is source-visible; do not perform another exploratory high-epsilon inventory:

- `src/cuts/pion_t_template_comparison.py::observe_template_comparison`, around lines 428–504, only accepts Q4p4W2p74 Center/Left lowe, and writes a `..._lowe_...` JSON/PDF.
- `src/cuts/pion_t_phi_template_comparison.py::requested/observe`, around lines 27–34 and 204–225, only accepts Q4p4W2p74 Center/Left lowe and writes `..._lowe_...` JSON. Its child results, unlike inherited parent curves, are independently fitted only when accepted.
- `run_Prod_Analysis.sh`, around lines 1707–1740, performs **both** canonical epsilon prepasses but its `-d` path intentionally exits after low epsilon; the normal path processes high epsilon. Around lines 338–346 the normal path also executes a `git clean` and symlink setup, and must never be launched against the ordinary checkout for this gate.
- `testing/run_pion_t_template_comparison_gate.py` is explicitly a Center/Left lowe execution owner with hardcoded settings, artifact basenames, low-only completion markers, A/B and child identity checks, source/checkout pins, private input staging, detached worktree, copied ltsep overlay, cleanup/preservation and strict collector ZIP verification. Its `copy_plan` currently treats high as **raw prepass only** and requires high `.root`, not high `.hist` and merged high pion controls.
- `testing/pion_t_phi_template_comparison_gate_profile_v2.json` includes five setting entries for collection context but only declares **six low-epsilon global artifacts**. Listing high settings is not proof of high analysis or collection.
- `src/binning/calculate_yield.py` already invokes the child observer without a low-epsilon-specific condition when the explicit opt-in is set. Preserve its numerical physics and do not modify it unless source inspection establishes an essential, contract-scoped defect.

## 3. Bounded implementation objective

Implement a **separately selectable, explicit high-epsilon diagnostic acquisition mode** for exactly Center/highe, Left/highe and Right/highe at Q4p4W2p74. Reuse the existing parent/child fitters, transfer code, physical cuts, frozen canonical bins, signed random/dummy bookkeeping, slow-proton cleaning, `no_empirical_residual`, source/charge/SIMC conventions and established control/yield definitions. Existing comparison diagnostic fits are expected; **do not introduce new models, fitting algorithms, tuning, weights, regularization, event selection or production adoption**.

Favor extending the existing owner behind an explicit high-epsilon option or a small high-epsilon wrapper that shares its validated internals. Preserve the existing low-epsilon owner CLI/profile and historical output behavior **bit for bit in default mode**. Avoid duplicating its ~1,280 lines or creating a separate orchestration framework. One tracked owner must control the complete high-epsilon sequence, including failures.

For the launcher: retain the existing `-d` low-only behavior by default. Add a narrow environment/owner opt-in that **only for kaon Q4p4W2p74, from the approved private worktree** permits `-d` to execute the existing high-epsilon `main.py` invocation after both canonical-binning prepasses and the low-epsilon phase. Ensure unrelated `-d` invocations remain low-only. Do not weaken the normal-launcher `git clean`/path protections; do not use the normal `-c` merge or run the full launcher in the ordinary farm checkout. Do not change the canonical binning physics, accepted low-epsilon analysis or protected signal templates.

For observers: accept the three exact high-epsilon settings only under the high-epsilon diagnostic opt-in and current physical fitting context; name all files with the real epsilon token, independently record setting/epsilon/source HEAD/input hashes and failure status, and preserve default low-epsilon behavior. The Center/Left/Right high output cannot silently reuse low input, output, fit, or provenance. Capture all **nine parent t comparisons** and the canonical 27 child slots per high setting (81 slots total), with independently fitted child models only where support and existing fitter acceptance allow. Preserve `empty_child`, `fit_failed`, `fit_rejected`, template-fallback and unavailable records; **never** synthesize zero-background estimates or mark failures as successful fits.

The high-epsilon parent comparison produces 6 diagnostic pages per setting under the existing parent renderer, if applicable. The child observer produces numerical JSON, not an invented production PDF. The established Work presentation pipeline can build the detailed high-epsilon four-group figures from the authenticated JSON plus any required archived numerical source objects after successful farm delivery; do not claim those figures already exist.

## 4. Actual runtime input/output and safety boundary

Keep all scientific execution in the **existing owner-created detached disposable worktree and privately isolated output namespace**. The actual installed ltsep `PATH_TO_DIR` must be sourced from the existing resolver and probed for baseline/private caller parity; the ordinary checkout, installed ltsep, accepted artifacts, external links, originals and shared outputs must remain unchanged. Do not create/fix external shared links or silently mutate the ordinary repository. Inspect the actual source of every writer (`findEffectiveCharge`, setup helpers, launcher, `main.py`, collector) and verify both lexical and resolved paths, including symlink aliases and PATH/`PYTHONPATH` import identities, using the actual farm interpreter at the preflight stage.

Before any child run, the owner must confirm **all three requested high-epsilon settings can actually execute**, including existing merged Data/Dummy **kaon and pion** inputs, their exact consumed SIMC `.root` and `.hist` files, any required matching runlist/efficiency references, and configured external K+Sigma0 inputs where applicable. Frozen `shift_prep` and common binning prepasses must resolve to the private copies and preserve selection/coordinate identity. Check exact consumed inputs and output collisions **before** expensive copying, work or full analysis; do not silently skip a missing high-epsilon setting or mix it with low-epsilon sources. Fail with the first precise missing role/path if a truly required input is absent.

**I/O cost is an acceptance condition.** Produce a declared inventory *before* heavy work showing explicit files, semantic role (byte-consumed versus existence-only), proposed copied bytes, expected full-read/hash bytes, and output-space requirements. No recursive hashing of the shared analysis tree; no full-content hash of existence-only skim ROOTs. Check skim existence with the same stat/existence semantics as the actual scientific consumer. Use cryptographic identity for scientifically consumed and staged inputs and retain content-preservation checks for those exact files. Do not weaken input protection to hide a source mutation. The owner must print/record file counts and predicted I/O scale, rejecting unbounded scope or shortage before processing.

Carry over strict name/provenance/canonical phi and t closure, signed coefficients, independent source classes, parent/child values and Sumw2, accepted/refused child fits, required PDF page count, payload status, actual-produced timestamp, JSON schema, `source_HEAD`, source-file identity and ZIP CRC/manifest/SHA-256 checks. Require successful high-epsilon completion markers and all **three** required settings, not a return code alone. Failed owner/checker/collector/preservation/manifest/page gates stop delivery; partial products remain inadmissible and must not be repackaged as successful evidence. Preserve a status receipt that states the *first failed invariant*, stage timestamps, actual environment, planned and actual byte costs, and preservation status even on failure. Do not add an excessive independent support/test system.

## 5. Prevent already documented farm failures

The independent actual-diff review must verify that the implementation handles each relevant historical failure class using the **real source call path**, not just isolated mocks:

1. `uproot`/Python 3.9 `lib64` package alias and namespace-module resolution: do not assume workstation Python 3.12 matches farm imports; preserve actual installed-package probe and provenance.
2. `caller_path_invalid` for `src/setup/findEffectiveCharge.py` and `/group` versus `/u/group` path topology: verify lexical and canonical caller roots and actual writer destinations.
3. `baseline_root_outpath_mismatch:effective_charge`: prove ordinary/private output resolution for each writer and every requested epsilon, including fallback and implicit output paths.
4. `source_changed_during_hash`: do not hash a 336.7-GB whole tree or existence-only skims; scope expensive reads to necessary consumed bytes with predictable I/O and stable-content checks.
5. Lost outputs, stale artifacts, path alias loops, untrusted environment variables, missing high-epsilon pion/kaon/SIMC inputs, wrong source checkout, source pin failure, and PDF/ZIP provenance inconsistency: all must fail before claims of acceptance.
6. The October 27 `child_failed_fit_evidence` gate: ensure the actual serialized missing-fit record matches the owner/checker's exact unavailable schema, and rejected child fits do not carry applied weights.

No known failed invariant may be silently bypassed. Where current environment testing cannot reproduce a farm-only dependency, report `NOT VERIFIED` and do not claim a farm preflight passed.

## 6. Required deterministic tests and review gate

Implement tests using the existing test framework, including:

- All five logical setting/epsilon selections and wrong-mode rejection; low-only behavior unchanged with high opt-in absent.
- A source-level or controlled command test showing both canonical prepasses, low completion, then high completion **only** under the new explicit opt-in; proof no normal checkout `git clean`, external writer, or unreviewed production path is exercised.
- All three high input plans, including high `.hist`, pion Data/Dummy and case-sensitive output names; negatives for missing paired inputs, optional vs required inputs, source/epsilon contamination, stale or collided outputs, symlink/alias escape.
- Existing parent/child JSON and histogram provenance invariants, source/child closure, signed source weights, Sumw2 and accepted-versus-rejected/failed/empty fits. Include an **actual serialized** failed/empty child round-trip, not only hand-written fixtures.
- Correct owner/profile settings, exact required artifact counts, actual PDF page checks, manifest and ZIP inventory, source pin and publication guards, and fail-closed on any absent Center/Left/Right high setting.
- Local-source parity for previously approved low-epsilon gate, complete preservation on both success and injected failure, and regression tests for the six failure categories above.
- `py_compile` on modified Python files; `bash -n` on launcher; relevant existing regression suites; `git diff --check`; and `python -B tools/check_memory_health.py --root .`.

**Review must inspect the actual full changed-source/diff**, including ROOT-call wiring, shell path, dependency pin changes, private I/O path, source allowlist, and parent-child result semantics. A Codex summary or synthetic unit-test pass is not approval. A new high-epsilon owner is NOT farm-ready until executable farm-appropriate smoke/readiness evidence at the reviewed/pushed SHA covers installed interpreter, package, callers, input inventory, output permissions/aliases and byte scope. No farm command is included in this task contract.

## 7. Explicit change scope

Only if required by the audited minimal high-epsilon acquisition, permit:

- `run_Prod_Analysis.sh` — opt-in continuation only, no default change.
- `src/cuts/pion_t_template_comparison.py` and `src/cuts/pion_t_phi_template_comparison.py` — epsilon/setting gating, output names and metadata only.
- `testing/run_pion_t_template_comparison_gate.py` — reuse/parameterize the established owner if the old mode stays identical; **or** add `testing/run_pion_high_epsilon_comparison_gate.py` as a narrow composition wrapper without duplicating its machinery.
- New `testing/pion_t_phi_template_comparison_highe_gate_profile_v1.json` for high settings; leave the existing low-profile `v2` frozen.
- `testing/test_run_pion_t_template_comparison_gate.py`, `testing/test_pion_t_template_comparison.py`, `testing/test_pion_t_phi_template_comparison.py`, plus a narrowly named new test file only if an additive high owner requires it.
- `docs/memory/CURRENT.md`, and the exact baseline-pion investigation and this task contract for the scoped state update. Other memory changes only if unavoidable and individually justified.

**Out of scope:** `src/main.py`, `src/binning/calculate_yield.py`, fitters, physics templates/config, cuts, normalization, HGCer methods, Method A/B, production outputs, scientific acceptance thresholds, broad launcher refactor, generic packaging infrastructure, unrelated source/test owners, and any farm execution. If a real blocker demonstrably requires an out-of-scope edit, **stop**, show the first invariant and propose a narrow amended contract; do not silently expand.

Before adding new paths, verify absence in the *current* checkout. Preserve the historical reviewed blob pins for unaffected files. Do not remove/relax an allowlist or source check to make the new commit appear valid; extend only the exact reviewed expected paths and identities. Preserve the existing low-epsilon owner and its regression tests.

## 8. Deliverable, decision and user-only publication

Codex delivers:

1. One actual complete review bundle: tracked diff plus complete additions for new/untracked files, and exact changed-path list. Do not stage merely for review. If large, provide an exact command that writes one complete temporary review file in the repository root for the user to upload.
2. Deterministic test commands and exact PASS/FAIL output, with real interpreter/environment identities, observed default low-epsilon regression result, and input/exit/output/copy/hash scope summary.
3. The MAINTENANCE full health check, including live source identity, worktree, source/CURRENT consistency, and the precise gate state.
4. A compact technical handoff identifying only unverified farm-specific facts; **no farm-ready claim** and **no farm command**.

**STOP** after implementation and deterministic checks. ChatGPT reviews the actual diff. Only after approval does the user review and perform a scoped commit/push; ChatGPT then verifies pushed state and the reviewed execution path. Prefer **one** user-operated acquisition with fail-fast checks when published source and directly applicable real-environment evidence already justify it. If a material farm-only interpreter, installed-package, path or I/O invariant remains `NOT VERIFIED`, authorize at most one narrow no-child smoke before acquisition, review its actual receipt, and stop on the first failed gate. Never collapse source review and pushed-state review, claim farm readiness from mocks, or rerun a failed gate automatically.

**Success criterion for this task:** a reviewable, complete, low-epsilon-regression-safe high-epsilon acquisition implementation with a realistic, bounded JLab preflight contract. It is not yet accepted high-epsilon physics.
