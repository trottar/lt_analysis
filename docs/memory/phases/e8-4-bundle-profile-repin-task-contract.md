# KaonLT E.8.4 — Narrow Farm-Bundle Profile Re-pin to Pushed E.8.4 Source

## Purpose

Prepare the existing E.8.1 generic full-background procedure-PDF validation
bundle for the first narrow E.8.4 farm gate.

The pushed E.8.4 source is:

```text
1aa1fd4184a6f8b20043e00ebb1ed3e9505a4935
Add E8.4 Method-A production-impact audit
```

The current profile
`testing/pion_hgcer_validation_bundle_profile_e8_1.json` is still pinned to the
older Fix.5 analysis source:

```text
53fd262b730af8f1254e411a38231aebeb6a1da3
```

and therefore intentionally rejects the later E.8.2/E.8.3/F.6.3/E.8.4 source
history as unexpected committed analysis changes.

This task changes **validation source provenance only**. It does not alter the
analysis, renderer, collector, wrapper, artifact definitions, accepted frozen
F.6.2 JSON, Method A/B, production physics, or farm launcher.

---

# 1. Exact starting state

Required branch:

```text
test
```

Required committed HEAD:

```text
1aa1fd4184a6f8b20043e00ebb1ed3e9505a4935
Add E8.4 Method-A production-impact audit
```

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git log -1 --oneline
git status --short
```

The worktree must be clean except for local-only repository instructions such
as root `AGENTS.md` / `.codex/` and this task contract after the user copies it
into `docs/memory/phases/`.

Hard stop on any unrelated tracked modification. Do not reset, stash, clean,
commit, push, or run the farm.

---

# 2. Mandatory startup reading

Read in order:

1. local root `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read only:

- `docs/memory/phases/e8-1-fix6-bundle-profile-repin.md`
- `docs/memory/phases/e8-4-production-impact-audit.md`
- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/roadmap/STATUS.md`
- this task contract
- `testing/pion_hgcer_validation_bundle_profile_e8_1.json`
- `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`
- `testing/collect_pion_hgcer_validation_bundle.py`
- `testing/package_pion_hgcer_validation_bundle.tcsh`

Do not reopen settled scientific architecture.

---

# 3. Required profile change

In:

```text
testing/pion_hgcer_validation_bundle_profile_e8_1.json
```

change only:

```text
source_identity.required_analysis_commit
```

from:

```text
53fd262b730af8f1254e411a38231aebeb6a1da3
```

to:

```text
1aa1fd4184a6f8b20043e00ebb1ed3e9505a4935
```

Preserve exactly:

- schema version `pion_hgcer_validation_bundle_profile/v4`;
- validation-profile identity
  `phase_e8_1_full_background_procedure_pdf_farm_review/v1`;
- collection mode `generic_artifacts`;
- ordered canonical-five setting inventory;
- global frozen F.6.2 acceptance-refinement-validation JSON;
- per-setting full-background procedure PDF;
- per-setting full-background page-manifest JSON;
- `allowed_committed_files` limited to:
  - `testing/pion_hgcer_validation_bundle_profile_e8_1.json`
  - `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`
- `allowed_non_analysis_path_prefixes = ["docs/memory/"]`.

Do not add E.8.4 source files to the allowlist. The point of the re-pin is that
the pushed E.8.4 analysis source becomes the new required analysis base, while
the later bundle/profile commit remains a distinct provenance commit.

---

# 4. Focused test change

In:

```text
testing/test_pion_hgcer_validation_bundle_profile_e8_1.py
```

update the exact reviewed-source constant/expectation to:

```text
1aa1fd4184a6f8b20043e00ebb1ed3e9505a4935
```

Preserve the tests that prove:

- exact five-setting inventory;
- exact artifact inventory;
- complete synthetic packaging;
- missing/invalid required artifacts produce incomplete bundles;
- required source ancestry is mandatory;
- only the profile/test pair and `docs/memory/` are allowed after the required
  analysis commit;
- later analysis changes remain fail-closed;
- output collision is fail-closed.

Do not weaken the source-identity test to admit arbitrary later `src/` changes.

---

# 5. Frozen implementation

Remain byte-unchanged:

```text
testing/collect_pion_hgcer_validation_bundle.py
testing/package_pion_hgcer_validation_bundle.tcsh
src/cuts/full_background_subtraction_plots.py
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
src/binning/calculate_yield.py
src/cuts/rand_sub.py
src/cuts/pion_component_subtraction.py
src/main.py
run_Prod_Analysis.sh
```

Also do not alter the accepted frozen F.6.2 artifact or any production physics.

---

# 6. Allowed changes

Substantive validation files:

```text
testing/pion_hgcer_validation_bundle_profile_e8_1.json
testing/test_pion_hgcer_validation_bundle_profile_e8_1.py
```

Warranted durable memory/history only:

```text
docs/memory/CURRENT.md
docs/memory/USER.md
docs/memory/manifest.json
docs/memory/phases/e8-4-bundle-profile-repin-task-contract.md
docs/memory/phases/e8-4-bundle-profile-repin.md
docs/memory/roadmap/STATUS.md
```

Do not change other files unless a concrete blocker requires stopping and
reporting it.

---


Also update:

```text
docs/memory/USER.md
```

with the durable workflow preference established by the user in this chat:

```text
Do not provide or authorize farm-run commands until the exact validation
bundle/profile required for that gate is ready, independently source-reviewed,
pushed, and pushed-state reviewed. Advance one workflow gate at a time; do not
provide commands for later gates before the current gate has passed.
```

This is a collaboration/workflow preference, not a scientific status change.

# 7. Deterministic local checks

Run:

```bash
python -B -m py_compile testing/test_pion_hgcer_validation_bundle_profile_e8_1.py
python -m json.tool testing/pion_hgcer_validation_bundle_profile_e8_1.json >/dev/null
python -B -m unittest testing.test_pion_hgcer_validation_bundle_profile_e8_1 -v
python -B -m unittest testing.test_collect_pion_hgcer_validation_bundle -v
```

Then memory integrity:

```bash
python -B tools/update_memory_manifest.py --root . --write
python -B tools/update_memory_manifest.py --root . --check
python -B tools/check_memory_health.py --root .
python -B tools/memory_bootstrap.py --root . --json
python -B -m unittest testing.test_memory_health -v
git -c core.safecrlf=false diff --check
```

These checks are local source/provenance checks only. They do not establish
ROOT/PyROOT, farm, procedure-PDF, or runtime validation.

---

# 8. Memory/status

Create:

```text
docs/memory/phases/e8-4-bundle-profile-repin.md
```

Record:

- starting pushed E.8.4 source `1aa1fd...`;
- previous stale profile source `53fd...`;
- new required analysis source `1aa1fd...`;
- profile/test-only allowlist preserved;
- artifact inventory unchanged;
- collector/wrapper unchanged;
- deterministic checks actually run by Codex;
- no farm run;
- no scientific or production change.

Until independent ChatGPT review, status is:

```text
ACTIVE
```

Do not change E.8.4 itself from `SOURCE REVIEWED`.

Final E.8 remains `BLOCKED`.
F.6.4 remains `BLOCKED`.

The immediate next gate after this task is only:

```text
independent ChatGPT review of the complete bundle-profile re-pin diff
```

Do not provide, prepare, or execute farm-run or packaging commands as part of
this task. Later gates are handled only after the current gate has passed.

---

# 9. Diff audit

Before stopping:

```bash
git status --short
git diff --stat
git diff -- testing/pion_hgcer_validation_bundle_profile_e8_1.json
git diff -- testing/test_pion_hgcer_validation_bundle_profile_e8_1.py
git diff -- docs/memory
git -c core.safecrlf=false diff --check
```

Confirm every frozen implementation file is byte-unchanged.

Create a complete timestamped review bundle containing the tracked diff plus
no-index sections for any intended new/untracked memory files.

Do not stage merely to create the review bundle.

---

# 10. Acceptance criteria

The re-pin is locally implementation-complete only if:

1. starting committed HEAD is exactly `1aa1fd...`;
2. profile required analysis commit is exactly `1aa1fd...`;
3. profile/test-only committed-range allowlist is unchanged;
4. canonical-five setting inventory is unchanged;
5. artifact inventory is unchanged;
6. collector and tcsh wrapper are unchanged;
7. focused profile tests and collector regressions pass;
8. memory checks and `git diff --check` pass;
9. E.8.4 remains `SOURCE REVIEWED`;
10. bundle-profile re-pin remains `ACTIVE` pending independent review;
11. `docs/memory/USER.md` records the one-gate-at-a-time / no-premature-farm-command preference;
12. no commit, push, farm run, or production change occurs.

---

# 11. Hard stop

Stop after the narrow provenance re-pin, deterministic checks, warranted memory
update, manifest regeneration, diff audit, and review-bundle creation.

Do not commit.
Do not push.
Do not run the farm.

The user alone controls commit/push.

The exact NEXT is only:

```text
independent ChatGPT review of the refreshed complete bundle-profile re-pin diff
```
