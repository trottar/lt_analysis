# KaonLT profile-ledger consumer repair

## Objective and authority

Implement one bounded bookkeeping repair, starting clean test/HEAD/origin/test and live GitHub test `19dbac361c08e6aec88e3c3bcc4cd77216ad83ad`. Recheck state; preserve any unrelated work. Follow AGENTS and CURRENT/MEMORY/handoffs/CURRENT_HANDOFF/USER in order, then MAINTENANCE, the October 29 high-epsilon acquisition/final-parent-publication contracts and the relevant investigation. This contract amends the earlier main.py freeze only for the exact paths/order described here. User authorization to repair the ongoing acquisition persists; publication and farm execution remain user-only.

The goal remains fresh validated Global/Local-t/Local-(t,phi) high-epsilon evidence through the existing owner. No new method, fitting, normalization, subtraction, model acceptance or infrastructure framework is authorized.

## Actual failure and source audit

The supplied excerpt reaches high Step 6 yields and ledger export, then completes Step 7 separation but raises FileNotFoundError opening generic `kaon_FullAnalysis_Q4p4W2p74_lowe_correction_ledger.json` at main.py:1700. Total 1 h 22 m 11.66 s; owner analysis_failed. Failed status basename: `KaonLT_HighEpsilon_Acquisition_20261029T221118-gate-status.json`. Actual executed SHA, full status, cleanup/preservation and surviving artifact bytes remain NOT VERIFIED; do not assign the current source SHA to that run without its receipt. Internal fit/completion messages do not establish owner or scientific acceptance.

Source: write_correction_ledger calls write_json_with_aliases with json and json_profile paths. The primary is json_profile. For no_empirical_residual the generic alias is intentionally not written. main.py nevertheless chooses ["json"] at the Step-6 debug/context assignment and both Step-7 low/high readers. The high generic reader would fail next after the low reader. This is a writer/consumer namespace mismatch, not a missing scientific producer or a reason to create an alias. main_auto/main_iter already prefer profile ledgers; they are outside this change. Low/high bg summaries are directly written at the same paths main reads; the input bundle uses its existing profile-aware loader. Step-8 missing optional archive items already warn/skip. Preserve those paths and semantics.

The earlier high-path audit should have caught this source-visible mismatch. The existing owner attempts guarded cleanup of its private data parent on failure; no successful artifact recovery can be assumed without the receipt. Do not manipulate remote attempts or introduce resume/retention infrastructure in this repair.

## Required implementation

1. Replace only the three main.py correction-ledger selections with the existing ["json_profile"] primary path. No generic fallback for a missing nonnominal profile ledger. For None/empty profile this path equals the generic path; nominal_weighted already writes both. Preserve all payload calculations and writer/alias policy.
2. Move the existing zeroth-bundle/low-high-entry/path selection and four JSON input reads to the start of the existing high Step-7 block, before Fortran/model activation and run_xsect. Retain the same inputs, errors, validation and variables. Keep separation execution, model activation order and final-summary timing/math unchanged. In particular leave the final-summary stage_start after separation; do not accidentally time separation twice. Missing or malformed finalization inputs must stop before expensive separation.
3. Audit remaining finalization and archive consumers for the same naming mistake; amend only this evidenced ledger path. Do not waive comparison/summary validators, invent a ledger, copy a stale generic file, create generic aliases, zero-fill missing metadata or skip finalization.
4. Update exact existing owner committed-source allowlists for main.py, the two scoped test paths and this contract. Preserve all ancestry/blob pins and unrelated-source rejection. Existing negative tests that used main.py as an unrelated committed change must now use a truly frozen unrelated source; dirty-main rejection remains required.

## Allowed paths

- src/main.py: three ledger path keys and relocation of the existing input-read block only. No numerical/physics change, new run mode, CLI change or skipping separation.
- testing/test_main_correction_ledger_paths.py: one new ROOT-independent producer/consumer test file.
- testing/test_e8_4_fix5_main_order.py: retain producer/finalizer ordering checks; extend the exact AST-normalization guard only for these three ledger keys and this one input-block relocation, preserving its prior SIMC-block allowance. Do not replace it with an unrestricted main.py exception.
- testing/run_pion_t_template_comparison_gate.py and its existing test file: exact source allowlists, temporary committed-history coverage and only necessary unrelated-source negative-fixture changes.
- docs/memory/CURRENT.md, the existing baseline-pion investigation, this contract under docs/memory/phases, and manifest.json.

Frozen: correction_ledger.py, frozen_manifest.py, epsilon comparison and final-summary calculations/validators, main_auto/main_iter, all scientific fits/cuts/templates/priors/normalizations/weights/uncertainties/yields, parent/child observers, launchers, profiles, installed packages, source isolation/cleanup, production adoption and accepted evidence. Stop at an essential out-of-scope blocker rather than expanding silently.

## Deterministic validation

Before the fix, reproduce the exact low generic-path FileNotFoundError using the actual unchanged ledger writer and AST-extracted actual main Step-6/Step-7 consumer statements. Stub only unused ROOT integration imports; do not replace file writing, path resolution, JSON encoding/reading or the consumer with a toy function. Control-flow fixtures are synthetic and establish no physics.

Positive cases: no_empirical_residual, nominal_weighted, None and empty profile; both low and high. Actual written primary JSON must equal main's recorded and consumed path, decoded data must match writer payloads, and generic aliases must remain absent for nonnominal profile. Seed a wrong/stale generic ledger to prove it is ignored. For missing or malformed required profile files, fail explicitly without fallback. Verify low/high bg summaries and bundle selection reach the same consumer context. Exercise the existing ledger comparison consumer using the resolved paths where appropriate, without changing thresholds.

Order test: the existing JSON reads precede run_xsect in the real main high block; final-summary timer and existing numerical producer order remain intact. The AST integrity guard must accept only the three path-key substitutions, exact input-block relocation and its pre-existing SIMC allowance, and still reject another source mutation. Test actual temporary Git histories with the newly allowed paths, and ensure unrelated frozen source changes still fail. Do not commit to or modify user refs/index for fixtures.

Run focused ledger/path tests with the available numerical Windows Python, AST/main-order tests, existing owner/POSIX tests under WSL and memory tests. Report versions, genuine skips and final outcomes. Compile changed Python sources without repo bytecode; Bash syntax only; whitespace including untracked additions; manifest, bootstrap and memory hard checks. No dependency installation or farm/ROOT/full-main execution.

## Memory, preservation and handoff

Replace only active high-acquisition paragraphs/NEXT in CURRENT, keeping it below 16 KiB; append the actual failure, precise namespace invariant and local evidence to the existing investigation. Retain previous failed attempts and evidence boundaries; no new maintenance phase. Regenerate manifest for the scoped memory/new files.

Preserve source outside the allowlist, USER, Git index/refs, installed ltsep, external symlinks, original inputs and all historical attempts. Produce a complete actual tracked diff plus full additions for this contract and the new test, a concise full-health handoff, exact local tests and limitations. Stop for independent ChatGPT actual-diff review. No staging/commit/push, no farm command, no automatic retry and no claim that prior outputs survived or are accepted. User publication and pushed-state review precede any further runtime decision.
