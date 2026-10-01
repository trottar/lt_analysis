# F.4.Refresh.2.Validation.2.Fix.1 — post-hardening source-allowlist repair Codex contract

## Objective

Repair one source-provenance defect found during independent ChatGPT actual-diff review of:

```text
kaonlt_review(20260930-211453).diff
```

The F.4.Refresh.2.Validation.2 execution-owner architecture is otherwise retained. The defect is narrow:

- the frozen F.4.Refresh.2 source boundary remains the reviewed materializer commit
  `141a3d04f9e5d07be21dba14e0e63212c3990bf1`;
- after that commit, the independently reviewed/pushed memory-health hardening changed two tracked non-scientific operational files:
  - `testing/test_memory_health.py`;
  - `tools/check_memory_health.py`;
- the current F.4.Refresh.2 profile and the new execution owner do not allow those two exact paths in their post-materializer committed-range provenance rule;
- therefore the owner would fail its own `git diff --name-only 141a3d04...<bundle commit>` preflight on the current branch even before materialization.

Repair only that mismatch. Do not redesign the owner, change the materializer source pin, broaden to a directory prefix, or alter scientific/runtime code.

The exact committed base remains:

```text
c406c138285727503b12115177dfb8bc7efcb7fe
Harden post-push CURRENT continuity
```

Remote `test` was independently checked by ChatGPT at that exact HEAD before this repair contract.

---

## Why this repair is warranted

Independent review of the pushed range:

```text
141a3d04f9e5d07be21dba14e0e63212c3990bf1
    ...
c406c138285727503b12115177dfb8bc7efcb7fe
```

established that, outside `docs/memory/`, the already-pushed range contains:

```text
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
testing/test_memory_health.py
testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py
tools/check_memory_health.py
```

The two memory-health files were part of the independently source-reviewed hardening change subsequently pushed at:

```text
3fd4edfcd05fd4cc25fb4e2119fb5f212ee916a1
```

and inherited unchanged by `c406c138...`.

Their current pushed Git blob identities are:

```text
testing/test_memory_health.py
8f7454dea48c00148dd379a8c141a005c2367d59

tools/check_memory_health.py
2736916c9d1ae63741c1da920ad6be7c14aa0bab
```

These files are repository-memory health infrastructure. They are not F.1/F.2/F.3/F.4 scientific builders, analysis runtime, accepted authority, or production physics.

The repair is therefore to permit these **two exact already-reviewed paths** in the F.4.Refresh.2 post-materializer committed-range rule. They must remain byte-identical. This is not permission to edit them.

---

## Required starting state

Before editing, read repository memory in the required order:

1. root `AGENTS.md` if present;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read only task-relevant records/source, including:

- `docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner-task-contract.md`
- `docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner.md`
- `docs/memory/phases/memory-health-operational-completeness-hardening.md`
- `docs/memory/phases/post-push-current-continuity.md`
- `docs/memory/phases/f4-refresh2-validation1-farm-materialization-bundle-profile.md`
- `docs/memory/investigations/f4-refresh2-validation1-operational-readiness-failure.md`
- `testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json`
- `testing/run_f4_refresh2_materialize_verify_package.py`
- their focused tests.

Verify:

```bash
git branch --show-current
git rev-parse HEAD
git status --short --untracked-files=all
git log -1 --oneline
```

Required committed state:

```text
branch: test
HEAD: c406c138285727503b12115177dfb8bc7efcb7fe
```

The existing Validation.2 candidate may contain exactly the paths already present in the reviewed bundle:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner-task-contract.md
docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner.md
docs/memory/roadmap/STATUS.md
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
testing/run_f4_refresh2_materialize_verify_package.py
testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py
testing/test_run_f4_refresh2_materialize_verify_package.py
```

plus this newly placed repair contract:

```text
docs/memory/phases/f4-refresh2-validation2-fix1-post-hardening-source-allowlist-task-contract.md
```

Root-level temporary `kaonlt_review*.diff` files may remain untracked and must not be staged, rewritten, or deleted by this task. Root `AGENTS.md` and `.codex/` remain local-only/untracked.

If committed HEAD differs, or unrelated tracked changes are present, **STOP**.

---

## Mandatory pre-edit provenance audit

Before modifying the candidate, run:

```bash
git diff --name-only \
  141a3d04f9e5d07be21dba14e0e63212c3990bf1..c406c138285727503b12115177dfb8bc7efcb7fe \
  | grep -v '^docs/memory/' || true

git rev-parse c406c138285727503b12115177dfb8bc7efcb7fe:testing/test_memory_health.py
git rev-parse c406c138285727503b12115177dfb8bc7efcb7fe:tools/check_memory_health.py
```

The non-memory committed-range inventory must be exactly:

```text
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
testing/test_memory_health.py
testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py
tools/check_memory_health.py
```

The two hardening blobs must be exactly:

```text
testing/test_memory_health.py
8f7454dea48c00148dd379a8c141a005c2367d59

tools/check_memory_health.py
2736916c9d1ae63741c1da920ad6be7c14aa0bab
```

If any additional non-memory path exists in that pushed range, or either hardening blob differs, **STOP** and report it. Do not broaden the repair.

---

## Allowed substantive changes

Modify only:

```text
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
testing/run_f4_refresh2_materialize_verify_package.py
testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py
testing/test_run_f4_refresh2_materialize_verify_package.py

docs/memory/CURRENT.md
docs/memory/roadmap/STATUS.md
docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner.md
docs/memory/manifest.json
```

Create exactly:

```text
docs/memory/phases/f4-refresh2-validation2-fix1-post-hardening-source-allowlist-task-contract.md
docs/memory/phases/f4-refresh2-validation2-fix1-post-hardening-source-allowlist.md
```

The original Validation.2 task contract must remain byte-identical:

```text
docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner-task-contract.md
```

No other tracked file may change.

---

## Frozen files and scientific ownership

Do not edit:

```text
testing/test_memory_health.py
tools/check_memory_health.py

testing/materialize_method_a_current_baseline_authority.py
testing/test_materialize_method_a_current_baseline_authority.py
testing/compare_method_a_current_baseline_authority.py
testing/test_compare_method_a_current_baseline_authority.py
testing/collect_pion_hgcer_validation_bundle.py
testing/test_collect_pion_hgcer_validation_bundle.py
testing/package_pion_hgcer_validation_bundle.tcsh
testing/test_package_pion_hgcer_validation_bundle_tcsh.py
```

Also freeze:

- every `src/` file;
- `src/main.py`;
- `run_Prod_Analysis.sh`;
- accepted F.2/F.3/F.4 authority artifacts/constants;
- F.1 scientific artifacts/builders;
- F.5/F.6.3/E.8.4 scientific/runtime logic;
- Method B;
- production physics/corrections.

The two memory-health paths are being **allowed**, not edited or re-reviewed scientifically.

---

## Exact source-identity repair

### Profile

Keep:

```text
required_analysis_commit =
141a3d04f9e5d07be21dba14e0e63212c3990bf1
```

Keep:

```json
"allowed_non_analysis_path_prefixes": ["docs/memory/"]
```

Do not move the required analysis commit forward.

The profile's `allowed_committed_files` must become exactly this six-path set:

```text
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py
testing/run_f4_refresh2_materialize_verify_package.py
testing/test_run_f4_refresh2_materialize_verify_package.py
testing/test_memory_health.py
tools/check_memory_health.py
```

No directory wildcard/prefix is allowed.

Do not add:

```text
testing/
tools/
src/
run_Prod_Analysis.sh
```

or any other blanket path rule.

### Execution owner

`testing/run_f4_refresh2_materialize_verify_package.py` must use the same exact six-path post-materializer committed allowlist as the profile.

Preserve all existing owner behavior from the reviewed candidate:

```text
preflight
  -> reviewed materializer subprocess
  -> post-materialization completion/provenance verification
  -> existing package wrapper
  -> returned-ZIP verification
```

Do not alter:

- materializer command construction;
- frozen comparison SHA;
- F.1 alias inventory;
- candidate output names;
- materialization verification;
- immutable packaging;
- ZIP verification;
- stale-state behavior;
- kinematic restriction;
- canonical artifact/Globus roots.

The only owner-source change should be the exact committed-range allowlist needed to recognize the already-reviewed hardening files.

---

## Focused regression tests

Update the focused tests so they prove all of the following.

### Positive provenance

A synthetic committed-range inventory containing:

```text
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py
testing/run_f4_refresh2_materialize_verify_package.py
testing/test_run_f4_refresh2_materialize_verify_package.py
testing/test_memory_health.py
tools/check_memory_health.py
docs/memory/CURRENT.md
```

passes the owner's committed-range check.

The profile test must assert the exact same six-path allowlist.

### Negative provenance

Each of these must still fail closed as an unexpected committed path:

```text
testing/unreviewed_helper.py
tools/unreviewed_helper.py
testing/materialize_method_a_current_baseline_authority.py
src/cuts/unreviewed.py
```

Do not weaken the existing scientific/materializer negative cases.

### Frozen hardening files

Tests must not modify or monkeypatch the real hardening source files themselves. Their permission is path provenance only.

---

## Local deterministic validation

Run at minimum:

```bash
python -m py_compile \
  testing/run_f4_refresh2_materialize_verify_package.py \
  testing/test_run_f4_refresh2_materialize_verify_package.py

python -m unittest testing.test_run_f4_refresh2_materialize_verify_package
python -m unittest testing.test_pion_hgcer_validation_bundle_profile_f4_refresh2
python -m unittest testing.test_materialize_method_a_current_baseline_authority
python -m unittest testing.test_collect_pion_hgcer_validation_bundle
python -m unittest testing.test_package_pion_hgcer_validation_bundle_tcsh
python -m unittest testing.test_memory_health
```

If local `tcsh` is unavailable, retain and report the existing wrapper-test skip exactly. That is not farm validation.

Re-run the pushed-range audit:

```bash
git diff --name-only \
  141a3d04f9e5d07be21dba14e0e63212c3990bf1..c406c138285727503b12115177dfb8bc7efcb7fe \
  | grep -v '^docs/memory/' || true
```

and re-check the frozen hardening blobs:

```bash
git rev-parse HEAD:testing/test_memory_health.py
git rev-parse HEAD:tools/check_memory_health.py
```

They must still be:

```text
8f7454dea48c00148dd379a8c141a005c2367d59
2736916c9d1ae63741c1da920ad6be7c14aa0bab
```

---

## Memory/status updates

Create:

```text
docs/memory/phases/f4-refresh2-validation2-fix1-post-hardening-source-allowlist.md
```

Record:

- independent ChatGPT review of `kaonlt_review(20260930-211453).diff` found one narrow source-provenance blocker;
- the owner architecture itself was not rejected;
- the blocker was the omitted already-reviewed/pushed memory-health paths after the frozen materializer source;
- no scientific/runtime file or accepted authority changed;
- the repair admits exactly those two historical operational paths;
- `required_analysis_commit` remains the reviewed materializer source;
- no farm operation was run;
- repaired source still requires independent ChatGPT actual-diff review before final pre-push reconciliation.

Keep Validation.2 status:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

for the repaired local candidate.

Do **not** mark it `SOURCE REVIEWED` in this task.

Update `CURRENT.md` and `roadmap/STATUS.md` only as needed to state the repair accurately.

The sole ordinary CURRENT NEXT after this repair must be:

```text
NEXT — independent ChatGPT actual-diff/source-runtime-path review of the repaired F.4.Refresh.2.Validation.2/Fix.1 execution-owner candidate; if it passes, separately contract final pre-push reconciliation.
```

Do not change `CURRENT_HANDOFF.md` unless an actual exceptional transfer state appears.

---

## Manifest and strict memory health

After versionable memory changes:

```bash
python -B tools/update_memory_manifest.py --root . --write
python -B tools/update_memory_manifest.py --root . --check
python -B tools/memory_bootstrap.py --root .
python -m unittest testing.test_memory_health
python -B tools/check_memory_health.py --root . --fail-on-warning
git -c core.safecrlf=false diff --check
```

Required report:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
manifest check: PASS | FAIL
```

Any warning under the strict command blocks completion.

---

## Diff audit and review bundle

Before stopping:

```bash
git status --short --untracked-files=all
git -c core.safecrlf=false diff --stat
git -c core.safecrlf=false diff --no-ext-diff
```

The cumulative candidate must contain only the original Validation.2 candidate plus this Fix.1 contract/phase and the warranted narrow edits listed above.

Create one fresh repository-root:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

It must contain:

1. branch and committed HEAD;
2. full `git status --short --untracked-files=all`;
3. cumulative tracked diff stat;
4. complete cumulative tracked `git diff --no-ext-diff`;
5. complete `git diff --no-index /dev/null ...` for every intended new/untracked candidate file;
6. the exact pushed-range provenance audit output;
7. the exact two frozen hardening blob checks;
8. all deterministic test/check commands and results;
9. memory-health report;
10. final cumulative changed-path inventory.

The review bundle remains untracked and must not be staged.

---

## Acceptance criteria

This repair passes locally only if all are true:

1. committed base remains `c406c138285727503b12115177dfb8bc7efcb7fe`;
2. the pushed-range audit proves the pre-existing non-memory range contains only the four expected paths;
3. `testing/test_memory_health.py` and `tools/check_memory_health.py` retain their exact pushed blob identities;
4. `required_analysis_commit` remains `141a3d04f9e5d07be21dba14e0e63212c3990bf1`;
5. profile and owner allow exactly six committed files after that source: original profile/test, owner/test, and the two hardening files;
6. no blanket `testing/`, `tools/`, `src/`, or other prefix is introduced;
7. arbitrary other testing/tools/scientific/materializer changes still fail closed;
8. owner runtime path remains materializer -> verify -> wrapper -> ZIP verify;
9. all frozen scientific/runtime/package/materializer/hardening files remain unchanged;
10. focused and regression tests pass, with any existing local `tcsh` skip reported accurately;
11. strict memory health and manifest checks pass with no warnings;
12. no farm, ROOT/PyROOT, full analysis, accepted-authority update, production correction, or Method-A promotion occurs;
13. repaired candidate remains `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`;
14. one fresh cumulative review bundle is produced for independent ChatGPT review.

---

## Hard stop

After the narrow repair, deterministic checks, memory update, manifest regeneration/check, strict memory-health report, diff audit, and fresh review bundle:

**STOP.**

Do not perform final pre-push reconciliation.

Do not commit.

Do not push.

Do not run the farm.

Do not provide a farm command.

Do not mutate accepted authority.

Do not begin F.6.3/E.8.4.

Return the fresh review bundle for independent ChatGPT actual-diff/source-runtime-path review.
