# KaonLT A/B owner — actual-farm volatile alias compatibility correction

## Authority, active science, and entry gate

- Authoritative repository: `https://github.com/trottar/lt_analysis/tree/test`.
- Reviewed starting pushed `test` HEAD: **`df07ecedf01b2d299c6fc61ac293c1fb4101ffd2`**. This is the five-path volatile-runtime-placement repair on top of `942fab7ed06ba53bff18ce053093a89a2cd55896`.
- **First establish actual local `test` HEAD, `origin/test`, index and full tracked/untracked worktree state.** Follow repository-root `AGENTS.md`; read in exact order `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, `docs/memory/USER.md`; then `docs/memory/MAINTENANCE.md`, the existing `kaonlt-ab-owner-volatile-runtime-placement-2026-10-25-task-contract.md`, and only relevant source/tests. Stop on overlapping edits, unexpected identity, or a genuinely incompatible local state. Preserve all unrelated work.
- **Science remains unchanged:** the requested next deliverable is a *single* private fresh Center/Left low-epsilon, canonical t1/t2/t3 A/B diagnostic at `Q4p4W2p74`, followed separately by the unchanged seven-strength C scientific comparison. Method A remains detached/non-production; Method B remains diagnostic and numerically excluded. No production promotion, physics modification, or new acquisition phase is authorized.

## First failing invariant — established before another farm command

The October 25 pushed owner (`testing/run_pion_t_template_comparison_gate.py`) correctly relocated its disposable `analysis`, `ab-output`, independent ROOT/hist copies, private `ltsep` overlay, caches and temporary files away from the ordinary `/group` repository. **However, its volatile-base preflight is inconsistent with the existing JLab farm's infrastructure aliases:**

```python
VOLATILE_MOUNT = Path("/volatile")

# volatile_base():
require(mount.resolve(strict=True) == mount and
        not mount.is_symlink() and mount.is_mount(), "volatile_mount_invalid")
identity = directory_identity(Path(baseline["VOLATILEPATH"]))

# directory_identity() currently requires path.resolve(strict=True) == path
# and disallows every symlink in its ancestry.
```

Historical real farm receipt `Q4p4W2p74_E8_2_baseline_five_8e68d2ac_20261009T223015Z-gate-status.json` establishes the **installed resolver's** `VOLATILEPATH=/volatile/hallc/c-kaonlt/trottar`; ordinary `OUTPUT` and SIMC-worksim symlinks resolve to `/lustre24/expphy/volatile/hallc/c-kaonlt/trottar/...`, and the no-repair external-link preflight passed. The prior October 16 farm `namei` evidence records legitimate aliases including `/volatile -> /lustre/enp/volatile`, `/lustre -> /lustre24`, and `/lustre24/enp -> expphy`. The configured logical `/volatile` path therefore need not itself be a physical mount point or a non-symlink directory. The October 25 `df` on the same farm reported approximately **4.4 PiB filesystem-wide free** for `/volatile` (NOT proof of the account's quota). Earlier on-ifarm checks verified installed Python 3.9 `replay_lt_env`, `ltsep`, owner/source preflight, 32 owner tests, and source-derived **19 input files / 8.85 GiB**; these precede the alias repair and do not prove the corrected owner runs.

**Hard invariant:** preserve the *literal* resolver-provided `VOLATILEPATH` and established external `OUTPUT` links, but validate its **canonical filesystem referent** and create all new owner-controlled data paths beneath the canonical, non-group volatile subtree. Reject genuine group-device aliases, mutable/escaping user-controlled symlinks, and changes in trusted anchor/base/owner-created identities. Do not weaken safety or add arbitrary new filesystem roots.

## Exact, minimum source correction

1. **Modify only the A/B owner's volatile alias handling.** In `testing/run_pion_t_template_comparison_gate.py`, adjust `volatile_base()`, `directory_identity()` and narrowly dependent placement/check functions as needed. Do not require the *lexical* `/volatile` alias to satisfy `resolve()==path`, `not is_symlink()` or `Path.is_mount()`. `Path.is_mount()` is not evidence that the logical alias itself is a mount, and the physical Lustre mount can sit above the canonical `.../volatile` subtree.

2. **Validate the actual canonical anchor, not a guessed destination.** Obtain the raw base **only** from unchanged `isolation.path_fields(repo, repo)` (`baseline["VOLATILEPATH"]`); require an absolute nonempty configured path. Resolve the existing `/volatile` infrastructure entry point and the configured base with `resolve(strict=True)`. Require the canonical base to be a **strict descendant** of the canonical volatile anchor, to exist as an accessible/writable directory and to be owned by the account as appropriate. Confirm the canonical base/anchor device and relevant filesystem identity are distinct from the actual `/group`/repository filesystem; reject canonical overlap with the repository, group, shared output, or any other forbidden authority. Use resolved filesystem identities/containment, not simple string-prefix checks. Do **not** posit that the canonical volatile anchor is itself `is_mount()`. A synthetic test may substitute the trusted anchor, but the real code must use the existing `/volatile` alias and resolver base.

3. **Differentiate trusted infrastructure aliases from owner-controlled symlinks.** Accept only the established root/infrastructure alias chain to the canonical volatile anchor; validate, snapshot and recheck the identity of the anchor and resolver base (including alias referents) before creating anything, during placement validation and before owner-scoped cleanup. A trusted infrastructure alias is not permission to accept arbitrary symlinks inside user-writable directories or the fresh private parent. The canonical base, fresh parent and all child data destinations must be represented by canonical, non-symlink paths, protected by the existing ownership/mode/device/inode/containment checks. Reject swapped links, escape through user-controlled aliases, wrong mount/device, changed owner or access, or unexpected alias rebinding. If any identity cannot be established safely, fail **before** large copying or `-d` child execution; unsafe cleanup fails closed without traversal of external targets.

4. **Do not modify the reviewed runtime topology.** Preserve the unique owner-created 0700 parent and the current `analysis`, `ab-output/Analysis/KaonLT`, `runtime-python/ltsep`, private TMP/XDG/MPL/pycache directories, exclusive `analysis/OUTPUT` symlink and literal/canonical verification. Preserve `analysis_worktree()`'s safe, bounded removal and failure receipts. Parent and children must be created under the verified **canonical** base (the supplied raw alias remains unchanged in the baseline/child `ltsep` fields). Existing ordinary repo `OUTPUT`, external SIMC links, archived ROOT inputs and accepted PDFs/JSON must be untouched. Do not move data back onto `/group` or create a second source checkout with data under `/group`.

5. **Keep copied-`ltsep` identical.** Do not edit `prepare_runtime_overlay()`'s unique installed config selection through `SetPath`, imported-module identity, `LTANAPATH`-only copied-byte substitution, untouched stable roles (particularly `VOLATILEPATH`), `private_callers` actual default/Plot/cwd probes or installed-package preservation. The frozen shared `isolation.prepare_runtime_overlay()` and shared `owned_worktree()` remain unchanged. Keep `copy_plan`, optional SIMC availability, required Center/Left lowe inputs, high-epsilon prepasses, Diamond namespace and first-tie parity, byte-verified independent copies, `skim_inventory`, artifacts, PDF/payload validators, independent clean-source collection, CRC/ZIP/manifest, status/log/summary and all scientific code untouched.

6. **Retain transparent evidence.** Existing `runtime_placement` records must expose the resolver's literal `VOLATILEPATH`, canonical volatile anchor and base, relevant device/UID/inode/alias identity, actual created parent/output paths and 19-file source-derived intended copy size. Do not mask a preflight failure, fake farm compatibility, or modify the owner to fall back to `/group` when volatile validation fails.

## Required regression coverage and acceptance

- In `testing/test_run_pion_t_template_comparison_gate.py`, create a genuine temporary Linux **symlink-chain fixture** representing the known logical-to-canonical alias shape: e.g. a synthetic `/volatile` alias through one or more `/lustre/...` aliases ending in a real `.../expphy/volatile` subtree. Patch only the test's configured anchor/group paths where necessary; **do not mock away resolution or `is_mount()` to conceal the reproduced failure**. First reproduce rejection by the pushed `df07ece` behavior, then prove the corrected `volatile_base` accepts the authentic alias chain and chooses the canonical owner-private parent. Use `Path.resolve(strict=True)`, symlinks and real synthetic Linux files where possible; no farm/PyROOT or full analyzer.
- Verify actual owner `analysis_worktree`, copied harmless `ltsep` overlay, private `ab-output`, selected 21 synthetic ROOT/hist copies (NOT replacements for real 19 farm inputs), caches, input-original preservation, independent `source` collector worktree and bounded success cleanup. Preserve exactly one child on mocked success; no child on any preflight failure.
- Negative fixtures: alias to group or repository; configured base outside resolved volatile anchor; user-writable or swapped alias; missing/dangling base; wrong device/owner/mode; changed anchor/base identity; worktree, output, cache or parent symlink swap; and unsafe cleanup. Assert no accepted/shared output modification, no staged ROOT copy before rejection, no recursive deletion outside the owned private parent.
- Run the affected Linux 40+ A/B tests, existing debug-owner/launcher/collector suites and detached template suite. Keep previous tests; explain any unavoidable platform symlink fixture limitation without claiming full PASS where tests did not run. Check both Python files with `py_compile`, `git diff --check`, memory manifest write/check, ordinary `tools/check_memory_health.py --root .` (nonblocking CURRENT soft-size warning may remain), actual diff and all changed-file identities. Preserve Python 3.9 syntax and installed-`ltsep` compatibility. All tests are deterministic local evidence, **not corrected ROOT/PDF or farm validation**.

## Strict edit allowlist

- `testing/run_pion_t_template_comparison_gate.py` — only alias/canonical volatile placement and dependent invariants/evidence.
- `testing/test_run_pion_t_template_comparison_gate.py` — real-alias positive/negative regressions.
- `docs/memory/phases/kaonlt-ab-owner-volatile-alias-compatibility-2026-10-25-task-contract.md` — this standalone contract, plus concise implementation results once completed.
- `docs/memory/CURRENT.md` — minimum necessary current work/NEXT continuity; do not reopen closed science or broaden status.
- `docs/memory/manifest.json` — regenerate and check.

**Forbidden:** `src/`, `run_Prod_Analysis.sh`, `farm_env/`, `testing/run_e8_4_fix5_canonical_five_plot_gate.py`, `testing/pion_t_template_comparison_gate_profile.json`, `ltsep` installed/config files, production outputs, unrelated tests/memory, source physics, user Git publication, and Jefferson Lab farm execution. If a narrowly necessary deviation is identified, stop with an exact blocker rather than editing outside this contract.

## Completion and handoff

1. Return the **complete actual tracked diff** and full new-contract contents for independent ChatGPT review, plus exact test/health results and unchanged-source assertions. Preserve all unrelated local worktree/index changes. Do not stage, commit or push.
2. After independent review, the **user** commits/pushes; ChatGPT verifies the actual pushed five-path successor against `df07ecedf01b2d299c6fc61ac293c1fb4101ffd2` and CURRENT/NEXT continuity.
3. Only then decide one **bounded read-only farm readiness gate** based on historical `VOLATILEPATH`, symlink/mount identity, *current* capacity/quota (4.4 PiB `df` is not user-quota proof), interpreter and input identities. All farm commands must be valid **interactive-safe `tcsh`**: never `exit`, `logout`, `exec`, unreviewed cleanup, or an automatically chained run after failure. Farm execution remains solely user-owned.
4. If source review and farm readiness both pass, run **one** private A/B owner, inspect fresh status/log, four A/B JSON/PDF artifacts, ZIP and preservation; then return immediately to the seven-strength C physics comparison. Corrected ROOT behavior, actual PDFs and runtime preservation remain **NOT VERIFIED** until supplied farm evidence passes.

**NEXT:** One alias-compatibility correction; no new science phase, no matched-cohort wrapper, no repeated support-contract cycle.


## Local implementation outcome - October 25, 2026

`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`. Starting test HEAD and local
origin/test matched df07ecedf01b2d299c6fc61ac293c1fb4101ffd2; staged content
was empty. The startup manifest discrepancy was exactly this new contract;
no existing entry was changed or missing. Pre-existing local work is preserved.

SOURCE VERIFIED: only directory_identity and volatile_base changed among
existing owner definitions. A small alias-identity helper snapshots the trusted
entry point's link chain, literal targets, referents, UID/device/inode and
protected parent identities. Every trusted link retains the anchor entry's
owner; writable user parents and extra base aliases are rejected. The canonical
anchor need not itself be a mount. The base is a strict canonical descendant
on the same device, distinct from the repository/group filesystem. The resolver
literal remains unchanged. Existing parent/worktree/output/cache guards and
cleanup automatically recheck the expanded placement record before writes and
cleanup. No new root, resolver, fallback, topology or acquisition was added.

All other existing function/class source bytes match the starting SHA,
including prepare_runtime_overlay, execute_gate, staging, optional-SIMC and
Diamond behavior, caller probes, preservation and artifact/collector gates.
The copied-ltsep algorithm is byte-for-byte unchanged. All 986 frozen/unrelated
versionable file byte identities, user branch/remote refs and empty staged
contents match startup. Internal automatic Codex refs are excluded from the
user-ref comparison. No staging, commit, push or farm action was performed.

### Deterministic Linux evidence

WSL Ubuntu, Python 3.12.3, account UID 1000. The fixtures use real /tmp and
/dev/shm filesystems with distinct device identities, real symlinks and strict
resolution. Synthetic infrastructure represents volatile -> lustre/enp/volatile,
lustre -> lustre24, and enp -> expphy. Protected fixture aliases are owned by
the fixture account and read-only; real infrastructure trust follows the actual
configured anchor owner. No resolve or is_mount behavior is mocked. The exact
pushed df07ece guard first reproduces volatile_mount_invalid on this topology;
the corrected guard then accepts its canonical, non-mount volatile subtree.

The actual owner worktree context, harmless copied installed-package/import/
SetPath probes, 21 independent ROOT/hist fixture copies (1017 bytes), caches
and exclusive OUTPUT link remain canonical and volatile. Raw VOLATILEPATH is
identical in baseline/child roles. The mocked end-to-end owner performs exactly
one child and independent source collection while private data exist, then
bounded cleanup before delivery. Git worktree operations and analysis are
mocked; installed-package probe subprocesses are real harmless fixture code.

Negative cases cover group/repository targets and devices, outside/missing/
dangling bases, writable/untrusted aliases, same-referent link rebinding,
changed base/worktree identity, wrong UID/mode, output/cache/parent swaps and
unsafe cleanup. Mutable alias and group placement fail before ROOT staging or
child execution. Shared accepted objects and original fixture inputs are
preserved; unsafe cleanup stops without traversing an external target.

```text
python3 -B -m unittest testing.test_run_pion_t_template_comparison_gate testing.test_run_e8_2_left_lowe_scientific_audit_gate testing.test_run_prod_analysis_debug_left_low testing.test_collect_pion_hgcer_validation_bundle testing.test_pion_t_template_comparison -q
```

139 passed, zero skips/failures: 46 owner (all prior 40 retained), 27 debug-owner,
15 launcher, 28 collector and 23 detached template tests. Final combined run:
28.832 seconds. Both touched Python files compile with temporary bytecode and
parse with Python 3.9 syntax; Python 3.9 execution is NOT VERIFIED. Whitespace,
generated memory-manifest and ordinary health checks accompany the review.
CURRENT's pre-existing soft-size warning is nonblocking and deferred to the
next checkpoint/milestone audit. The supplied real 19-file / 8.85 GiB inventory
remains historical farm evidence, not replaced or remeasured by these fixtures.

No corrected ROOT/PDF/farm validation or scientific acceptance is claimed.
NEXT: independent actual-diff review, user publication, pushed-state review,
then the contract's single bounded read-only farm readiness gate before any
private A/B run. No farm command or new supporting phase is introduced.
