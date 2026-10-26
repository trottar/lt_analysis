# A/B private staging: restore the consumed pion Data/Dummy inputs

## Authority and audit

Starting branch `test`, HEAD and local `origin/test` must equal
`b5769e51c01655d9756481c91b258f029e219211`; preserve all unrelated state.
Read root AGENTS, CURRENT, MEMORY, CURRENT_HANDOFF, USER in that order, then
MAINTENANCE. User alone publishes Git and executes farm validation.

The supplied October 26 ifarm2402 attempt
`KaonLT_pion_t_template_AB_Q4p4W2p74_b5769e5_20261026T145834Z-gate-status.json`
failed `analysis_failed` after about 4m51: private Center pion Analysed
Data/Dummy lowe files were absent. This is supplied failure evidence, not a
locally inspected complete receipt or scientific acceptance.

SOURCE VERIFIED audit: `run_Prod_Analysis.sh -d 4p4 2p74` defaults to kaon;
merging existing skim trees is confined to `-c`, which this owner never invokes.
The missing files are pre-existing merged pion inputs, not products of the
kaon debug child. `src/main.py` runs paired raw prepasses, returns early for
capture, then executes full available low settings. `rand_sub` sets
SubtractedParticle=pion and opens both files via
`_open_subtracted_particle_tree_bundle` before component fitting. Current
`copy_plan` includes kaon and archived SIMC plus the kaon-only Diamond namespace;
it omits every pion input. `stage_copies` copies exactly that plan.

At accepted October 19 owner source
`6d13abb4681b938c56def2a1f455244b0880096c`, OUTPUT linked to ordinary shared
outputs, making existing pion inputs available. The launcher, main and rand_sub
are byte-identical between that commit and this starting HEAD. Do not restore
shared writes. Existing shared pion file presence/contents on today's farm are
NOT VERIFIED locally; fail early if they are absent or invalid.

## One repair and exact allowlist

1. `testing/run_pion_t_template_comparison_gate.py`: add pion Data/Dummy lowe
   inputs to the existing copy inventory for each setting that passes the
   existing kaon/archived-SIMC availability guards. Center and Left remain
   mandatory. Right lowe also needs both pion files when available; Right with
   absent archived SIMC stays optional and skipped downstream. High epsilon
   remains raw prepass only and must not acquire a pion dependency.
   Require each consumed pion input to be a regular file before copying any
   inputs or invoking analysis. Report the exact missing filename. Preserve
   copy-byte accounting, hashes, distinct inodes, provenance and final checks.
2. `testing/test_run_pion_t_template_comparison_gate.py`: extend existing input
   fixtures and regression coverage; keep prior cleanup/path/mode/source guards.
3. `docs/memory/phases/kaonlt-ab-stage-pion-inputs-2026-10-26-task-contract.md`:
   this single contract and concise measured local outcome.
4. `docs/memory/CURRENT.md`: minimally record supplied failure, repair and
   substantive NEXT conditional on independent review/publication/synchronization.
5. `docs/memory/manifest.json`: regenerate/check.

All other files are frozen, especially src/, launcher, farm_env, shared owner
helpers, profiles, collectors, installed ltsep, external links and old attempts.
No physics, cuts, weights, fits, priors, binning, normalizations, subtraction,
yield/error or scientific ownership changes. Method A stays detached and
Method B stays diagnostic and numerically excluded. No redesign or fallback.

## Representative verification

- Execute the frozen pure file-opener function extracted by AST, replacing only
  ROOT opening/tree retrieval with fixture handles. Demonstrate the old plan
  misses the exact Center pair, and repaired staged Center/Left/available Right
  pairs reach the unchanged opener with private paths and equal bytes.
- Each Data/Dummy file absent separately and both absent; directory and dangling
  symlink input; mandatory Center/Left; available Right; unavailable Right with
  no pion files; high epsilon with no pion files. Missing required inputs must
  fail before staging/child, not after minutes of analysis.
- Plan/copy race, collision, source mutation, private inode/hash identities,
  original/shared preservation and Diamond namespace remain strict. Existing
  mocked full-owner flow must assert pion inputs at child first use.
- Run Linux affected owner, debug-owner, launcher, collector and detached
  comparison tests. Run applicable shared-helper tests using a working existing
  interpreter if WSL lacks NumPy; report exact platform and missing tests.
  Parse changed Python with Python 3.9 grammar, compile without science imports,
  run diff whitespace checks, manifest write/check and ordinary memory health.

## Acceptance, review and hard stop

Deliver only the five allowed paths, complete actual tracked diff plus full new
contract addition without staging, exact local results and frozen-source proof.
Preserve index, branch and remote refs; never commit, push, fetch, clean or run
farm. Existing failed directories/artifacts are untouched. If frozen source or
other infrastructure needs changing, stop with the exact blocker under this
contract; no second support contract or speculative farm rerun.

NEXT: independent ChatGPT actual-diff review, then user publication and
pushed-state verification, then decide one private Center/Left lowe t1/t2/t3
A/B gate using the unchanged owner CLI. Verify updated copy_plan/status contains
the pion pairs, fresh four A/B files, existing preservation and collection
checks; then resume the unchanged C comparison. No farm command before those
reviews. Local tests establish neither ROOT/tree schema, rendered PDFs, farm
preservation nor scientific acceptance. Independent review is not self-review.


## Local implementation outcome — October 26, 2026

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING. GitHub live test, local HEAD and
origin/test were confirmed at b5769e51c01655d9756481c91b258f029e219211 before
implementation; worktree/index were clean. No publication, farm or cleanup of
previous attempts occurred. The implementation adds seven lines to copy_plan;
stage_copies, source/history/blob checks, launcher and science are unchanged.

Before the fix, the new two-test Linux run reproduced the exact Center missing
file at the frozen production opener and nine missing-pair subcase failures.
After the fix, all 154 affected Linux tests passed (Python 3.12.3): A/B owner,
debug owner, launcher, collector and detached comparison. Real temporary copies,
symlinks, separate group/data filesystem devices, prior-attempt preservation,
02700 mode, caller selection, cleanup, freshness and provenance guards remain
covered. Six new test methods cover this repair, including failure before any
stage_copies or child call and unchanged first-use opening for Center/Left/Right.
The complete fixture inventory is 27 files/1332 bytes, previously 21/1017;
actual farm inventory sizes remain unknown and are reported by the real plan.

Shared helper: 36 Linux operational/cleanup tests passed in 5.530 s via the
established exact-test-AST route, supplying only accepted_pages fixture data
generated with the existing Windows interpreter. No test methods/guards were
rewritten. The one frozen scientific-lineage reconstruction test was NOT RUN
successfully: direct WSL import lacks NumPy. A direct Windows 37-test attempt
also failed (43 subtest failures, 16 errors), principally because its POSIX
cleanup safety/UID interfaces are unavailable. These attempts are not passes;
the unchanged scientific-lineage test remains NOT VERIFIED. No environment
installation or safety bypass was performed.

Changed Python files pass Python 3.9 grammar parsing and py_compile; this is not
Python 3.9 execution. Diff whitespace, manifest write/check and ordinary memory
health are required before delivery. Existing CURRENT soft-size warning is
nonblocking and batched; no scientific state is promoted. The contract and
CURRENT preserve the fresh A/B then C goal and one ordinary NEXT.

Author actual-source/diff audit: the only executable change is consumed lowe
pion-pair inventory after existing availability guards; high prepass, optional
Right, complete Diamond kaon namespace and original/copy verification stay
unchanged. Current launcher/main/rand_sub equal accepted October 19 source.
Independent ChatGPT actual-diff review remains NEXT; this author's audit is
not that independent acceptance. Complete tracked diff and full new contract
addition accompany delivery. Farm ROOT trees/schema, present shared pion files,
copy cost/quota, actual current preservation and rendered PDFs are NOT VERIFIED.
No farm command is issued. After independent review, user publication and
pushed-state review, decide the existing single private A/B gate; missing pion
inputs now name the exact file before staging/analysis instead of failing late.
