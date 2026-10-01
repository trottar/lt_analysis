# F.4.Refresh.2.Validation.2/Fix.1 — Final Pre-Push Source-Review Reconciliation

## 1. Purpose

Perform the required memory/status reconciliation after independent ChatGPT
actual-diff/source-runtime-path review **PASSED** the repaired
F.4.Refresh.2.Validation.2/Fix.1 tracked execution-owner candidate.

Reviewed cumulative bundle:

```text
kaonlt_review(20261001-092421).diff
```

Required committed base:

```text
c406c138285727503b12115177dfb8bc7efcb7fe
Harden post-push CURRENT continuity
```

Independent ChatGPT review established that:

- live remote `test` remained exactly
  `c406c138285727503b12115177dfb8bc7efcb7fe`;
- the cumulative candidate is based on that exact committed HEAD;
- the tracked execution owner is orchestration-only;
- scientific candidate construction remains delegated to the unchanged reviewed
  F.4.Refresh.2 materializer;
- package/collector ownership remains delegated to the unchanged reviewed
  package wrapper and generic collector;
- the owner path is:
  preflight -> reviewed materializer subprocess -> explicit materialization
  completion/provenance verification -> reviewed package wrapper -> returned-ZIP
  verification;
- the frozen materializer source pin remains
  `141a3d04f9e5d07be21dba14e0e63212c3990bf1`;
- the frozen F.4.Refresh.1 comparison SHA-256 remains
  `c16287e9192288ed5f116b5eebf95f49261d7fd09755ee923070154e047cf7b5`;
- Fix.1 correctly admits only the two already-reviewed post-materializer
  hardening paths omitted by the first Validation.2 candidate:
  `testing/test_memory_health.py` and `tools/check_memory_health.py`;
- those two files remained byte-identical to pushed `c406c138...` with Git blob
  identities:
  `8f7454dea48c00148dd379a8c141a005c2367d59` and
  `2736916c9d1ae63741c1da920ad6be7c14aa0bab`;
- both the profile and execution owner now use exactly the same six-path
  post-materializer committed allowlist:
  `testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json`,
  `testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py`,
  `testing/run_f4_refresh2_materialize_verify_package.py`,
  `testing/test_run_f4_refresh2_materialize_verify_package.py`,
  `testing/test_memory_health.py`,
  `tools/check_memory_health.py`;
- `docs/memory/` remains the only allowed non-analysis committed-path prefix;
- arbitrary additional `testing/`, `tools/`, materializer, or `src/` changes
  still fail closed;
- the repaired owner retained fresh-output/no-overwrite behavior, exact
  comparison identity, exact kinematic gate, immutable five-output packaging,
  and returned-ZIP hash verification;
- the original Validation.2 owner source differed from the repaired owner only
  by the two exact hardening allowlist entries;
- no scientific/runtime source, accepted authority, production correction,
  Method-A scientific builder, Method-B path, or production physics changed;
- no farm execution, ROOT/PyROOT, full `main.py`, production analysis,
  candidate authority acceptance, F.6.3/E.8.4 runtime closure, or Method-A
  production promotion was established.

The existing pushed delegated component identities independently observed during
review were:

```text
testing/materialize_method_a_current_baseline_authority.py
bd00b14338d6549840685258caaa54f26fb784d8

testing/collect_pion_hgcer_validation_bundle.py
e5eb836ad65226a81ad98573c6af4ab948fd1ed3

testing/package_pion_hgcer_validation_bundle.tcsh
0896db4ff7b3d029eae284af55eb2d34fea6ed16
```

Codex-reported local deterministic checks in the reviewed bundle were:

```text
py_compile owner + focused test                         PASS
owner suite                                            11 tests OK, 0 skips
F.4.Refresh.2 profile suite                             4 tests OK, 0 skips
unchanged materializer suite                           13 tests OK, 0 skips
generic collector suite                                28 tests OK, 0 skips
package-wrapper suite                                   3 tests, 2 OK, 1 existing tcsh-unavailable skip
memory-health suite                                    36 tests OK, 0 skips
manifest write/check                                   PASS
memory bootstrap                                       PASS
strict memory health --fail-on-warning                 PASS, zero warnings
git diff --check                                       PASS
farm / ROOT / PyROOT / full analysis                  NOT RUN
```

Those unit suites were **NOT RUN by ChatGPT**. ChatGPT reviewed the actual
cumulative diff, current pushed source identities, source/provenance path, and
the supplied deterministic-check output.

This task is **memory/status reconciliation only**. It must not alter the
independently reviewed substantive source candidate.

---

## 2. Exact starting state

Required branch:

```text
test
```

Required committed HEAD:

```text
c406c138285727503b12115177dfb8bc7efcb7fe
```

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short --untracked-files=all
git log -1 --oneline
```

The cumulative source-reviewed candidate before adding this reconciliation
contract must contain exactly these intended versionable paths:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/f4-refresh2-validation2-fix1-post-hardening-source-allowlist-task-contract.md
docs/memory/phases/f4-refresh2-validation2-fix1-post-hardening-source-allowlist.md
docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner-task-contract.md
docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner.md
docs/memory/roadmap/STATUS.md
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
testing/run_f4_refresh2_materialize_verify_package.py
testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py
testing/test_run_f4_refresh2_materialize_verify_package.py
```

This newly placed reconciliation contract is additionally permitted:

```text
docs/memory/phases/f4-refresh2-validation2-final-pre-push-source-review-reconciliation-task-contract.md
```

Temporary repository-root files matching:

```text
kaonlt_review*.diff
```

may remain untracked. They remain outside the candidate and must not be staged,
rewritten, or deleted by this task.

Root `AGENTS.md` and `.codex/` remain local-only/untracked as established by the
repository workflow.

If committed HEAD differs or unrelated tracked changes exist, **STOP**.

Do not reset, stash, clean, discard, commit, push, run the farm, materialize
farm artifacts, mutate accepted authority, or run the full analysis.

---

## 3. Mandatory repository-memory startup

Read in this exact order:

1. root `AGENTS.md` if present;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read only task-relevant records:

```text
docs/memory/roadmap/STATUS.md
docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization.md
docs/memory/phases/f4-refresh2-validation1-farm-materialization-bundle-profile.md
docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner-task-contract.md
docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner.md
docs/memory/phases/f4-refresh2-validation2-fix1-post-hardening-source-allowlist-task-contract.md
docs/memory/phases/f4-refresh2-validation2-fix1-post-hardening-source-allowlist.md
docs/memory/phases/memory-health-operational-completeness-hardening.md
docs/memory/phases/post-push-current-continuity.md
docs/memory/investigations/f4-refresh2-validation1-operational-readiness-failure.md
```

Do not reopen adjacent architecture or scientific phases.

---

## 4. Files that must remain byte-identical

Do not edit the independently source-reviewed substantive candidate:

```text
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py
testing/run_f4_refresh2_materialize_verify_package.py
testing/test_run_f4_refresh2_materialize_verify_package.py
```

Do not edit either implementation contract:

```text
docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner-task-contract.md
docs/memory/phases/f4-refresh2-validation2-fix1-post-hardening-source-allowlist-task-contract.md
```

Do not edit the frozen delegated or hardening files:

```text
testing/materialize_method_a_current_baseline_authority.py
testing/test_materialize_method_a_current_baseline_authority.py
testing/compare_method_a_current_baseline_authority.py
testing/test_compare_method_a_current_baseline_authority.py
testing/collect_pion_hgcer_validation_bundle.py
testing/test_collect_pion_hgcer_validation_bundle.py
testing/package_pion_hgcer_validation_bundle.tcsh
testing/test_package_pion_hgcer_validation_bundle_tcsh.py
testing/test_memory_health.py
tools/check_memory_health.py
```

Also freeze every scientific/runtime/production source path, including all
`src/` files, `src/main.py`, `run_Prod_Analysis.sh`, accepted F.2/F.3/F.4
authorities, F.5/F.6.3/E.8.4 runtime logic, Method B, and production physics.

Before memory edits, record the exact `git hash-object` identities of the four
source-reviewed substantive candidate files. After reconciliation, verify the
same four identities again. If any substantive candidate byte changes, **STOP**.

Also verify the two hardening files remain:

```text
testing/test_memory_health.py
8f7454dea48c00148dd379a8c141a005c2367d59

tools/check_memory_health.py
2736916c9d1ae63741c1da920ad6be7c14aa0bab
```

and the delegated pushed source blobs remain:

```text
testing/materialize_method_a_current_baseline_authority.py
bd00b14338d6549840685258caaa54f26fb784d8

testing/collect_pion_hgcer_validation_bundle.py
e5eb836ad65226a81ad98573c6af4ab948fd1ed3

testing/package_pion_hgcer_validation_bundle.tcsh
0896db4ff7b3d029eae284af55eb2d34fea6ed16
```

If any frozen identity differs, **STOP**.

---

## 5. Allowed reconciliation changes

Modify only:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner.md
docs/memory/phases/f4-refresh2-validation2-fix1-post-hardening-source-allowlist.md
docs/memory/roadmap/STATUS.md
```

Create exactly:

```text
docs/memory/phases/f4-refresh2-validation2-final-pre-push-source-review-reconciliation-task-contract.md
```

No other versionable file may change.

This task may only reconcile status, evidence boundary, and push-stable NEXT.
It must not edit implementation or tests.

---

## 6. Required status reconciliation

### 6.1 F.4.Refresh.2.Validation.2

Advance the local execution-owner status from:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

to:

```text
SOURCE REVIEWED
```

Record that independent ChatGPT actual-diff/source-runtime-path review of:

```text
kaonlt_review(20261001-092421).diff
```

passed.

Record the exact source-reviewed execution chain:

```text
explicit validated inputs
    ->
tracked F.4.Refresh.2 execution owner
    ->
unchanged reviewed materializer
    ->
explicit manifest/hash completion verification
    ->
unchanged reviewed package wrapper
    ->
unchanged generic collector + F.4.Refresh.2 profile
    ->
exact returned ZIP verification
```

Record that this is source review only and does not establish farm execution.

### 6.2 F.4.Refresh.2.Validation.2.Fix.1

Advance Fix.1 from:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

to:

```text
SOURCE REVIEWED
```

Record that the repaired candidate:

- preserved `required_analysis_commit =
  141a3d04f9e5d07be21dba14e0e63212c3990bf1`;
- allowed exactly the six reviewed post-materializer paths;
- preserved `docs/memory/` as the only allowed non-analysis prefix;
- retained exact pushed hardening-file blobs;
- rejected arbitrary extra `testing/`, `tools/`, materializer, and `src/`
  changes;
- changed no scientific/runtime/production source.

### 6.3 Preserve upstream/downstream status

Keep:

```text
F.4.Refresh.1 / Fix.1       CLOSED / RUNTIME VALIDATED
F.4.Refresh.2 / Fix.1       SOURCE REVIEWED
F.4.Refresh.2.Validation.1  SOURCE REVIEWED
F.6.2 / Fix.5               CLOSED / RUNTIME VALIDATED
E.8                          ACTIVE
E.8.2                        SOURCE REVIEWED
E.8.3                        SOURCE REVIEWED
F.6.3                        SOURCE REVIEWED
E.8.4                        SOURCE REVIEWED
final E.8                    BLOCKED
F.6.4                        BLOCKED
```

Preserve the narrow historical E.8.4.Fix.4 and E.8.1.Fix.5/Fix.6 closures
exactly as already recorded.

Do not reopen F.1-F.6.2 accepted science.

Do not promote Method A.

Method B remains diagnostic-only.

---

## 7. Farm-readiness boundary

After this reconciliation, farm readiness is still:

```text
BLOCKED
```

because the source-reviewed owner has not yet been:

1. user committed/pushed;
2. independently pushed-state reviewed.

Do not provide or run a farm command in this reconciliation.

Do not claim:

- farm materialization;
- farm package-wrapper execution;
- ROOT/PyROOT validation;
- full `main.py` validation;
- candidate-authority acceptance;
- F.5 authority refresh;
- F.6.3 runtime acceptance;
- E.8.4 runtime acceptance;
- production acceptance;
- Method-A production promotion.

After a future successful user push and pushed-state review, the next substantive
gate is one narrow F.4.Refresh.2 `Q4p4W2p74` materialize -> verify -> package
farm operation using explicit validated canonical-five F.1, accepted F.2/F.3/F.4,
and reviewed comparison inputs. That later farm command must be prepared only
after pushed-state review verifies the complete tracked chain.

---

## 8. CURRENT.md push-stable NEXT

CURRENT must retain exactly one ordinary NEXT.

Replace the pre-review NEXT with a push-stable substantive NEXT equivalent to:

```text
NEXT — after user-controlled commit/push and pushed-state review of the SOURCE REVIEWED F.4.Refresh.2.Validation.2/Fix.1 execution owner, prepare the single narrow Q4p4W2p74 F.4.Refresh.2 materialize -> verify -> package farm gate; do not begin F.6.3/E.8.4 until the returned F.4.Refresh.2 evidence is reviewed.
```

The wording may be line-wrapped but must preserve that meaning exactly.

Do not make commit/push itself the NEXT.

Do not describe already pushed `c406c138...` work as local/unpushed.

Do not say the new owner is pushed before the user actually pushes it.

---

## 9. Validation/evidence boundary to record

Record accurately:

```text
ChatGPT actual-diff/source-runtime-path review          PASS
reviewed bundle                                         kaonlt_review(20261001-092421).diff
live remote test observed during review                 c406c138285727503b12115177dfb8bc7efcb7fe

Codex-reported owner suite                              11 tests OK, 0 skips
Codex-reported profile suite                             4 tests OK, 0 skips
Codex-reported materializer suite                       13 tests OK, 0 skips
Codex-reported collector suite                          28 tests OK, 0 skips
Codex-reported package-wrapper suite                     3 tests, 1 tcsh-unavailable skip
Codex-reported memory suite                             36 tests OK, 0 skips
Codex-reported strict memory health                     PASS, zero warnings
Codex-reported manifest check                           PASS
Codex-reported git diff --check                         PASS

ChatGPT unit tests                                      NOT RUN
farm materialization                                    NOT RUN
farm packaging                                          NOT RUN
ROOT/PyROOT                                             NOT RUN
full main.py                                            NOT RUN
```

Do not silently convert the local synthetic collector test's printed synthetic
ZIP path into a farm-runtime claim.

---

## 10. Manifest and strict memory-health checks

After all versionable memory changes:

```bash
python -B tools/update_memory_manifest.py --root . --write
python -B tools/update_memory_manifest.py --root . --check
python -B tools/memory_bootstrap.py --root .
python -m unittest testing.test_memory_health
python -B tools/check_memory_health.py --root . --fail-on-warning
git -c core.safecrlf=false diff --check
```

Strict health must return zero warnings.

Required report:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
manifest check: PASS | FAIL
```

Do not modify `CURRENT_HANDOFF.md` absent a real exceptional transfer state.

No substantive source test rerun is required for a byte-identical memory-only
reconciliation. If deterministic suites are rerun, report them accurately and
do not change the evidence boundary.

---

## 11. Cumulative diff audit

Before stopping:

```bash
git status --short --untracked-files=all
git -c core.safecrlf=false diff --stat
git -c core.safecrlf=false diff --no-ext-diff
```

The cumulative intended versionable candidate must contain exactly:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/f4-refresh2-validation2-final-pre-push-source-review-reconciliation-task-contract.md
docs/memory/phases/f4-refresh2-validation2-fix1-post-hardening-source-allowlist-task-contract.md
docs/memory/phases/f4-refresh2-validation2-fix1-post-hardening-source-allowlist.md
docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner-task-contract.md
docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner.md
docs/memory/roadmap/STATUS.md
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
testing/run_f4_refresh2_materialize_verify_package.py
testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py
testing/test_run_f4_refresh2_materialize_verify_package.py
```

Temporary `kaonlt_review*.diff` files may remain untracked and are outside the
candidate inventory.

Do not use `git add -A`.

---

## 12. Required final pre-push review bundle

Create one fresh repository-root bundle:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

It must contain byte-faithfully:

1. branch;
2. committed HEAD;
3. full `git status --short --untracked-files=all`;
4. cumulative `git diff --stat`;
5. complete cumulative tracked `git diff --no-ext-diff`;
6. complete `git diff --no-index /dev/null ...` for every intended untracked
   candidate file, including this final reconciliation contract;
7. pre/post `git hash-object` identities for the four source-reviewed
   substantive candidate files;
8. frozen hardening and delegated-component blob checks;
9. exact reconciliation-check commands/results;
10. memory-health report;
11. final cumulative changed-path inventory.

Do not stage the review bundle.

---

## 13. Acceptance criteria

The reconciliation is acceptable only if:

1. committed HEAD remains
   `c406c138285727503b12115177dfb8bc7efcb7fe`;
2. the four source-reviewed substantive candidate files are byte-identical to
   the candidate reviewed in `kaonlt_review(20261001-092421).diff`;
3. both implementation task contracts remain byte-identical;
4. hardening files retain exact pushed blob identities;
5. materializer, collector, package wrapper, scientific/runtime source, and
   accepted authorities remain unchanged;
6. Validation.2 is recorded `SOURCE REVIEWED`;
7. Validation.2.Fix.1 is recorded `SOURCE REVIEWED`;
8. F.4.Refresh.2/Fix.1 and Validation.1 remain `SOURCE REVIEWED`;
9. historical runtime closures remain unchanged;
10. F.6.3/E.8.4 remain `SOURCE REVIEWED` and runtime blocked;
11. farm readiness remains `BLOCKED` pending user push and pushed-state review;
12. CURRENT contains one push-stable substantive NEXT, not commit/push itself;
13. no farm/runtime/production claim is added;
14. manifest check passes;
15. strict memory health passes with zero warnings;
16. `git diff --check` passes;
17. one fresh complete cumulative review bundle is created.

---

## 14. Hard stop

After memory/status reconciliation, manifest regeneration/check, strict
memory-health checks, cumulative diff audit, and creation of the fresh
pre-push review bundle:

**STOP.**

Do not edit substantive implementation source.

Do not commit.

Do not push.

Do not run the farm.

Do not provide a farm command.

Do not materialize candidate authority.

Do not change accepted authority.

Do not begin F.6.3/E.8.4.

Return the fresh cumulative review bundle for independent ChatGPT final
pre-push reconciliation review.
