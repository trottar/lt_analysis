# E.8.4 Fix.5.7 — page-manifest setting provenance / owner repair

**Status:** `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`.

## Source and scope

Started on `test` at `cdc6ead47be3c46987a8418b796d07d87b80f829` under the
[contract](e8-4-fix5-7-page-manifest-setting-provenance-owner-repair-task-contract.md).
Accessible Git checks established no pre-existing tracked modifications; the
supplied contract was the sole untracked file. Initial sandbox-denied reads
produced an incorrect dirty-file observation, corrected by accessible read-only
Git checks. The startup manifest mismatch was solely that new contract:
no existing entry changed or disappeared. Authorized manifest regeneration
restored manifest and ordinary health PASS before substantive editing.

## SOURCE VERIFIED root cause and repair

`rand_sub.py` takes the complete checkpoint setting; the checkpoint producer
serializes kinematic token, Q2, W, semantic epsilon, filename epsilon, phi and
particle. Render-state capture detaches that setting, and the page-manifest
builder preserves its complete mapping. The owner previously compared exact
dictionary equality to four keys, rejecting valid extra provenance. The old
four-key test fixture masked that integration mismatch.

The verifier now checks a mapping and exact Q4p4W2p74 / lowe / Left / kaon
identity plus semantic `epsilon_setting=low`, allowing producer Q2/W metadata
without rewriting the manifest. The regression builds the real pure-Python
checkpoint setting and passes it through the real page-manifest builder for
both focused verification and the complete synthetic owner path. Wrong or
missing identity and invalid setting types fail with the unchanged literal
`page_manifest_setting_invalid`; successful and failed verification leave the
manifest unchanged. Existing page and owner safety gates remain intact.

Only verifier/tests and warranted memory change. Scientific/presentation source,
launcher, collector/tests, profile, candidate identities and eight-artifact
inventory remain frozen. Method A remains detached/non-production; Method B
remains diagnostic-only and numerically absent. Absolute-SIMC provenance and
all prior accepted runtime scopes remain unchanged.

## Verification and NEXT

The implementation was `ACTIVE` before completing deterministic checks.
Python 3.12: the required command was:

```text
python -B -m unittest testing.test_run_e8_4_fix5_left_lowe_plot_gate testing.test_full_background_subtraction_plots testing.test_pion_hgcer_refinement_checkpoint testing.test_collect_pion_hgcer_validation_bundle testing.test_pion_hgcer_validation_bundle_profile_e8_4_fix5
```

159 tests ran: 141 passed, 18 existing skips (4 unavailable-PyROOT and
14 superseded historical procedure-tail assertions); no failures/errors.
Owner-only suite: 16 passed with no skips. Syntax checks use in-memory
`compile(...)` for both changed Python files. AST comparison to starting HEAD
confirms `verify_pages` is the only changed owner function/class; only a Mapping
import changes elsewhere in that owner. `git diff --check` passes.
All 164 frozen tracked files match starting HEAD after checkout line-ending
normalization; accessible Git diff also confirms no frozen file changed.

The complete synthetic path retains detached source preflight/rechecks,
collection/ZIP verification, status success flags, eight-artifact inventory,
status exclusion, one-path success stdout and empty failure stdout. The page
tests retain missing/duplicate, scope, phi-inventory and renderer-failure gates,
and explicitly exercise schema and t-identity failures. No launcher runs and
no real worktree/ref operation is performed by these mocked owner tests.

Final manifest check and ordinary memory health pass with one nonblocking
soft-size warning: CURRENT is 8593 bytes (limit 8192). MEMORY is 12489 bytes;
CURRENT_HANDOFF is 323 bytes. No active-state/provenance ambiguity or hard
integrity failure remains. Batch CURRENT compaction at the next memory
checkpoint/milestone audit. The manifest includes only tracked memory and the
two intended Fix.5.7 additions. The complete temporary root review diff includes
all five changed tracked files and both new files; remove it before commit.

NEXT — after independent actual-diff review, user commit/push and ChatGPT
pushed-state synchronization, one narrow Q4p4W2p74 / Left / lowe tracked-owner
farm gate and fresh status/ZIP/log/PDF/manifest review. CURRENT owns this sole
next action; no farm command is authorized before synchronization.
No commit, push, remote-ref update, farm execution or farm command is authorized
or performed. ROOT/PyROOT, full analysis, PDF rendering/legibility, farm
filesystem behavior, real freshness/delivery and numerical interpretation
remain NOT VERIFIED here.
