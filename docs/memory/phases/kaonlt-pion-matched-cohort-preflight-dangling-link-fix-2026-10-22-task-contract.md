# KaonLT — Repair unrelated dangling-link preflight without delaying matched-cohort physics

**Task class:** One narrowly scoped failed-farm-gate repair within the existing matched-cohort acquisition. **Not a new phase.**

**Audit anchor:** Live remote `test` at `30e029c744ca2638733a1f4c070a71968641f157` on 2026-10-22. Codex must establish its actual local `test` HEAD, `origin/test`, index, staged/unstaged/untracked state and source/diff before editing. Preserve unrelated dirty state; stop on a conflicting edit or source discrepancy.

**Immediate evidence:** The user-supplied `KaonLT_matched_cohort_Q4p4W2p74_20261022-133634-gate-status.json` reports `status=failed`, `stage=source_preflight`, `analysis_started=false`, `evidence_admissible=false`, source HEAD `30e029…`, and `failure_reason="[Errno 2] No such file or directory: '/lustre24/expphy/volatile/hallc/c-kaonlt/trottar/OUTPUT/Analysis/HeeP'"`. No event evidence, accepted ZIP or science was produced. This HeeP path is *a reported failure target*, not a new input or a directory to create/repair.

**Single scientific goal retained:** One user-operated, owned `Q4p4W2p74` matched-cohort acquisition for **Center/lowe and Left/lowe**, canonical `t1/t2/t3`, producing fresh, provenance-matched A/B templates and signed-event packets, then the already-implemented offline seven-strength C refits and event-applied full-MM and Λ-region pion-background/variance comparisons. Do not select λ, infer a physical MSE minimum, promote production or claim final yields.

## 1. Authority, source-level diagnosis, and boundary

Read root `AGENTS.md`, followed **in order** by `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, and `docs/memory/USER.md`. Then read `docs/memory/MAINTENANCE.md`, `docs/memory/CODEX.md`, `docs/memory/TOOLS.md`, `docs/memory/COMMUNICATION.md`, the existing `docs/memory/phases/kaonlt-pion-matched-cohort-attempt-isolation-2026-10-21-task-contract.md`, and only its direct active references needed for this fix. Emit the full KaonLT health check.

**SOURCE VERIFIED at audit anchor:** `testing/run_pion_local_t_matched_cohort_gate.py::inventory` calls `p.resolve(strict=True)` for **every** symlink encountered by broad directory inventory; a dangling unrelated symlink therefore raises `FileNotFoundError` before `p.is_file()` can distinguish a required consumable file. `execute_gate` calls `inventory(shared)` from `shared_state()` and `external_input_inventory()` calls `inventory(root)` over `CUTPATH` and `PARAMPATH`. A source-preflight failure from an unrelated link is therefore reproducible by a targeted synthetic fixture. The supplied gate status does **not** give a stack trace proving which invocation encountered the actual farm HeeP link; distinguish reproduced code behavior from that unverified farm caller.

Trace the earliest source-preflight read that could have raised the observed missing-HeeP exception, including the owner and relevant frozen delegates. Confirm by an executable deterministic regression (not a conjecture alone). If the failed invariant is not this broad-inventory symlink dereference, **stop `BLOCKED`** with the first verified actual caller rather than guessing, creating HeeP directories or widening the repair.

## 2. Minimal repair, preserving strict science safety

Modify only the *broad preservation/metadata inventory* treatment of symlinks whose referent does not exist. Preserve their lexical path, literal link target, and explicit unresolved/missing-target state **without dereferencing a missing target**. Include the record in the before/after equality check so retargeting, disappearance, replacement with a file, or becoming resolvable during the child is caught, not ignored. Leave normal regular-file hashing, intact symlink referent identity where already hashed, and scope of the shared-output preservation inventory unchanged.

A dangling path is acceptable **only as an irrelevant metadata record**; do **not** treat it as an available or safe input. Maintain existing strict checks for every actually consumed input, configured Sigma0, original/private ROOT files, background and SIMC sources, fixed output destinations, required directory existence, installed/copied `ltsep`, external links required by the real child, source checks, writer audit, and preservation after failure. Never downgrade an `identity(path)` failure, `required_external_input_missing`, `external_link_changed`, or output collision into success. If the root problem instead involves a required KaonLT link, fail closed and identify its actual role.

Do **not** create/modify HeeP or any shared/external directory, remove/relink a farm symlink, weaken original versus private-source lineage, fabricate input files, change `VOLATILEPATH` or turn off the owner's preservation gate. Do not modify launcher, ROOT/PyROOT, production physics, generic collector, installed `ltsep`, existing A/B/C algorithms, signal definitions, cuts, weights, templates, priors, SIMC, yields, normalization, efficiencies, or cross sections. No change to the 7-strength sweep or downstream replay/plots.

## 3. Tests — minimum directly relevant evidence

Add focused deterministic synthetic tests to the **existing** `testing/test_run_pion_local_t_matched_cohort_gate.py`:

1. Broad shared-output inventory containing an unrelated symlink pointing to an absent `OUTPUT/Analysis/HeeP` must record it consistently, without creating/reading the missing target or stopping acquisition preflight.
2. A similar dangling symlink under a broad `CUTPATH` or `PARAMPATH` scan must not block otherwise sufficient configured KaonLT inputs.
3. Before/after inventory must detect a dangling link being changed, removed, replaced or becoming resolvable during the child; existing regular-file byte mutations must still fail.
4. A required KaonLT input (e.g., configured Sigma0 or one of the consumed neutron/delta/sidis files) resolving through a dangling symlink must **still fail** with the existing required-input invariant. Configured input mutation and original-file-preservation tests remain strict.
5. Existing no-child-on-hard-failure, no-collection-on-failure, exact one-child invocation, protected output boundary, live-link preservation and source-validation tests must remain green.

Run targeted test module(s), the source-owned suites affected by the edit, `py_compile`, and `git diff --check`. If an optional environment dependency prevents a test, report it as not run, not passed, without starting another infrastructure repair. These synthetic checks do not establish farm ROOT validation.

## 4. Strict allowlist and memory checkpoint

Allowed source edits:
- `testing/run_pion_local_t_matched_cohort_gate.py`
- `testing/test_run_pion_local_t_matched_cohort_gate.py`

Allowed memory edits for this **one failed gate**:
- `docs/memory/CURRENT.md` (replace stale pre-push NEXT and identify blocked preflight / conditional science goal succinctly)
- `docs/memory/investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md` (one brief failure/repair checkpoint, with the supplied failure distinguished from source inference)
- `docs/memory/manifest.json` (generated only)
- `docs/memory/phases/kaonlt-pion-matched-cohort-preflight-dangling-link-fix-2026-10-22-task-contract.md` (this single contract)

**No other files.** The existing `CURRENT.md` is close to its 16-KiB hard limit: replace obsolete text or compact in-place; do not append enough prose to break memory health. No unrelated memory consolidation, generic path subsystem, extra testing framework, new science phase or adjacent helper edit. If the narrow safe repair requires a frozen file or expanded allowlist, stop and report the blocker.

After memory edits, run `<PYTHON> -B tools/update_memory_manifest.py --root . --write` and `--check`, then `<PYTHON> -B tools/check_memory_health.py --root .` (ordinary nonstrict form). Report exact bytes and any nonblocking warnings; run bootstrap if required by the existing local workflow.

## 5. Acceptance and user handoff

Codex should produce a complete actual diff and all new-file content for independent ChatGPT review. Report the reproduced failing invariant, exact repaired behavior, tests run, unchanged protected source and memory health. **Codex does not commit/push or execute the JLab farm.**

Required progression: targeted regression and fix → independent ChatGPT actual-diff review → **user** commit/push → ChatGPT pushed-state synchronization → **one** user-run guarded farm acquisition using reviewed owner → fresh provenance/manifest/ZIP verification → offline A/B/C six-parent signed-event MM and variance report → scientific interpretation. No farm command is authorized by this contract before the independent pushed-state review. If a gate fails, stop on that invariant; do not automatically rerun, collect or interpret inadmissible artifacts.

**Definition of success for this fix:** The unrelated HeeP dangling-link condition no longer blocks metadata-only inventory, while required-input and preservation failures remain hard stops. The actual science remains to be established by a subsequently accepted farm artifact.
