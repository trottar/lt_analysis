# KaonLT F.4.Refresh.2.Validation.1 — Farm Materialization Bundle/Profile Gate

## 1. Purpose

Create one narrow, source-reviewed farm-validation profile for the already-pushed
F.4.Refresh.2 current-baseline Method-A candidate materializer.

This gate must **not** change the materializer, comparator, scientific builders,
accepted authorities, production analysis, or any Method-A/Method-B physics.

The pushed materializer source is:

```text
141a3d04f9e5d07be21dba14e0e63212c3990bf1
F4 Refresh.2: add current-baseline candidate materializer
```

The materializer is already independently source reviewed. It has not yet been
run on the farm.

This task adds only the declarative generic-bundle profile and its focused test,
plus warranted durable memory. The profile must package the five deterministic
F.4.Refresh.2 materialization outputs and must pin the materializer source above.
It must reuse the existing generic v4 validation-bundle collector unchanged.

The later farm gate, after this profile is independently source reviewed, pushed,
and pushed-state reviewed, will:

1. run the pushed materializer directly against the canonical-five current F.1
   artifacts, accepted F.2/F.3/F.4 artifacts, and the validated F.4.Refresh.1
   comparison JSON;
2. write the detached candidate lineage into a task-local output directory;
3. collect exactly the five materialization outputs with this profile into one
   fresh ZIP;
4. return that ZIP for independent evidence review.

This task does **not** run any of those farm steps.

---

## 2. Exact starting state

Required branch:

```text
test
```

Required committed starting HEAD:

```text
141a3d04f9e5d07be21dba14e0e63212c3990bf1
```

Required commit subject:

```text
F4 Refresh.2: add current-baseline candidate materializer
```

Before editing, run:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
git log -1 --oneline
```

Hard requirements:

- branch must be `test`;
- committed HEAD must be exactly
  `141a3d04f9e5d07be21dba14e0e63212c3990bf1`;
- inspect the current worktree before editing;
- pre-existing temporary `kaonlt_review*.diff` files may remain untracked and
  must not be deleted, rewritten, staged, or treated as candidate files;
- root `AGENTS.md` and `.codex/` remain local-only/untracked;
- if unrelated tracked changes exist, STOP and report the blocker;
- do not reset, stash, clean, discard, commit, push, or run the farm.

The previously supplied user push showed only temporary review bundles left
untracked after the commit. Re-verify the actual local worktree rather than
assuming that state remains unchanged.

---

## 3. Mandatory repository-memory startup

Read in this exact order:

1. root `AGENTS.md`;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read only task-relevant records/source:

- `docs/memory/CODEX.md`
- `docs/memory/TOOLS.md`
- `docs/memory/COMMUNICATION.md`
- `docs/memory/decisions/farm-validation-bundle-procedure.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/phase-f6-method-a-production-promotion.md`
- `docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization.md`
- `docs/memory/evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md`
- `testing/materialize_method_a_current_baseline_authority.py`
- `testing/test_materialize_method_a_current_baseline_authority.py`
- `testing/collect_pion_hgcer_validation_bundle.py`
- `testing/test_collect_pion_hgcer_validation_bundle.py`
- `testing/pion_hgcer_validation_bundle_profile_f6_1.json`
- `testing/test_pion_hgcer_validation_bundle_profile_f6_1.py`

Use the F.6.1 profile/test only as an architectural example of the existing
generic collector. Do not copy F.6.1 scientific content or expand this task.

---

## 4. Established evidence and status entering this task

Preserve these facts:

```text
F.4.Refresh.1         — CLOSED / RUNTIME VALIDATED
F.4.Refresh.1.Fix.1   — CLOSED / RUNTIME VALIDATED
                         detached comparison gate only

F.4.Refresh.2         — SOURCE REVIEWED
F.4.Refresh.2.Fix.1   — SOURCE REVIEWED

historical F.1–F.6.2  — CLOSED / RUNTIME VALIDATED

F.6.3                 — SOURCE REVIEWED
E.8.4                 — SOURCE REVIEWED
E.8                   — ACTIVE
final E.8             — BLOCKED
F.6.4                 — BLOCKED
lifecycle hook        — BLOCKED / DEFERRED
```

F.4.Refresh.1 direct farm evidence remains:

```text
comparison file:
KaonLT_F4_Refresh1_authority_comparison_Q4p4W2p74_20260930-134654.json

comparison SHA-256:
c16287e9192288ed5f116b5eebf95f49261d7fd09755ee923070154e047cf7b5

F.2 scientific payload match = true
F.3 scientific payload match = true
F.4 scientific payload match = false
first changed stage = F4

accepted F.2 raw SHA-256:
87e0a16870ea34a6611702e1bd4fc02c4dc6a8daa132e75b725a1b782a9a31da

accepted F.3 raw SHA-256:
04a9a576767ba3f81c8c03daab2461998894c655f0d875bd532ad30cb46c0d95

accepted F.4 raw SHA-256:
adcc01900b8adc211e296aaedaa314fdb86207d51483ebfd4b3edf69f558f188

candidate F.2 serialized SHA-256:
182433ccd2d13d3f7af80b963e69190ad8f2eb4e295c8216b36585e8ee7409b2

candidate F.3 serialized SHA-256:
eeb480d9b7ddffab1e97c7f9f27c3099d0780c2bbb4d42a0ee68c077a5752b3
```

The pushed F.4.Refresh.2 materializer remains detached and non-authoritative. It
requires explicit canonical-five current F.1 inputs, accepted F.2/F.3/F.4,
the reviewed comparison JSON plus expected raw SHA, an explicit source HEAD,
and an existing output directory. It rebuilds and validates the candidate
lineage before publishing its manifest last.

No accepted authority has changed.

---

## 5. Scientific ownership and frozen boundaries

This task owns only:

- one declarative generic validation-bundle profile;
- one focused pure-Python profile test;
- durable status/history needed to record this source-development gate.

This task does **not** own scientific calculation.

Do not modify:

- `testing/materialize_method_a_current_baseline_authority.py`;
- `testing/test_materialize_method_a_current_baseline_authority.py`;
- `testing/compare_method_a_current_baseline_authority.py`;
- `testing/test_compare_method_a_current_baseline_authority.py`;
- `testing/collect_pion_hgcer_validation_bundle.py`;
- `testing/test_collect_pion_hgcer_validation_bundle.py`;
- any `src/` file;
- any F.1/F.2/F.3/F.4 scientific builder;
- F.5;
- F.6.1/F.6.2 accepted science;
- F.6.3;
- E.8.4;
- random/dummy subtraction;
- slow-proton subtraction;
- pion component fitting/subtraction/alignment;
- Method B;
- SIMC;
- cuts/windows/templates/priors/binning;
- efficiencies/acceptance;
- yield extraction;
- L/T separation;
- cross sections;
- active `no_empirical_residual` behavior.

Do not update accepted F.3/F.4 runtime-authority constants or replace accepted
F.2/F.3/F.4 artifacts.

Do not promote Method A.

---

## 6. Allowed substantive files

Create exactly:

```text
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py
docs/memory/phases/f4-refresh2-validation1-farm-materialization-bundle-profile-task-contract.md
docs/memory/phases/f4-refresh2-validation1-farm-materialization-bundle-profile.md
```

Allowed existing durable-memory updates:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization.md
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/roadmap/STATUS.md
```

No evidence record is created in this task because no new farm/runtime evidence
exists yet.

If implementation appears to require any collector, materializer, comparator,
scientific, runtime, or production source change, STOP and report the blocker.

---

## 7. Exact profile contract

Create:

```text
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
```

It must use:

```text
schema_version:
pion_hgcer_validation_bundle_profile/v4

validation_profile:
phase_f4_refresh2_current_baseline_candidate_materialization_farm_review/v1

collection_mode:
generic_artifacts
```

The profile must retain the collector-required canonical-five setting inventory
in this exact order:

```text
Left / lowe
Left / highe
Center / lowe
Center / highe
Right / highe
```

This profile has **no setting-scoped artifacts**. Therefore:

```json
"settings": []
```

inside the `artifacts` object is required.

The top-level profile `settings` array still contains the canonical five above,
because the existing generic collector validates that fixed setting inventory.
Do not change the collector merely to remove this existing requirement.

---

## 8. Exact packaged artifact inventory

The profile must declare exactly five **global**, **required**, `json`
artifacts.

Use these exact keys and basename templates:

### 8.1 Reviewed comparison-input copy

```text
key:
f4_refresh1_comparison_input

basename_template:
{kinematic}_kaon_pion-background_hgcer-method-a-current-baseline-authority-comparison-input.json
```

### 8.2 Candidate F.2

```text
key:
f4_refresh2_candidate_f2

basename_template:
{kinematic}_kaon_pion-background_hgcer-method-a-acceptance-representation-current-baseline-candidate.json
```

### 8.3 Candidate F.3

```text
key:
f4_refresh2_candidate_f3

basename_template:
{kinematic}_kaon_pion-background_hgcer-method-a-acceptance-map-current-baseline-candidate.json
```

### 8.4 Candidate F.4

```text
key:
f4_refresh2_candidate_f4

basename_template:
{kinematic}_kaon_pion-background_hgcer-method-a-parent-preserving-correction-current-baseline-candidate.json
```

### 8.5 Materialization manifest

```text
key:
f4_refresh2_materialization_manifest

basename_template:
{kinematic}_kaon_pion-background_hgcer-method-a-current-baseline-authority-materialization-manifest.json
```

For all five:

```text
kind = json
required = true
```

Do not package:

- accepted F.2/F.3/F.4 files redundantly;
- the five current F.1 files redundantly;
- F.5/F.6.1/F.6.2/F.6.3/E.8.4 artifacts;
- procedure PDFs;
- Method B;
- production outputs.

The five-file package is the smallest complete F.4.Refresh.2 output set because
the copied comparison input contains the reviewed detailed comparison and input
identities, while the materialization manifest records the verified source/input
lineage and the candidate-output identities.

If later evidence review finds that an additional pre-existing accepted artifact
is genuinely required, that is a later evidence request, not a reason to broaden
this source gate preemptively.

---

## 9. Exact source-identity boundary

The profile must pin the already-reviewed and pushed materializer source:

```text
required_analysis_commit:
141a3d04f9e5d07be21dba14e0e63212c3990bf1
```

This is intentionally the **materializer source commit**, not the later
profile/bundle commit.

The materializer's future farm invocation must likewise receive:

```text
--source-head 141a3d04f9e5d07be21dba14e0e63212c3990bf1
```

Do not infer or substitute the later profile commit as materializer provenance.

The profile's exact `allowed_committed_files` must be:

```text
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py
```

The exact allowed non-analysis prefix remains:

```text
docs/memory/
```

No collector, materializer, comparator, or scientific source file may be
allowlisted after the required source commit.

This ensures any later modification of the materializer or scientific path
fails the bundle source-provenance gate rather than being silently admitted.

---

## 10. Preserved runtime / packaging path

The source/runtime path after this change must be:

```text
reviewed pushed materializer source at 141a3d...
    -> later direct farm materializer invocation
    -> five deterministic non-authoritative output JSON files
    -> unchanged generic v4 validation-bundle collector
    -> new F.4.Refresh.2 generic profile
    -> fresh ZIP with:
         global/<five exact JSON files>
         manifest.json
         source_state.txt
         source_checks.txt
         canonical setting directories with no declared artifacts
```

The collector remains read-only with respect to source artifacts.

The collector's bundle `complete=true` means only that the declared package and
source checks succeeded. It does **not** accept the materialized candidate
scientifically or as source-owned authority.

Independent evidence review must later inspect both:

1. the collector `manifest.json`; and
2. the packaged F.4.Refresh.2 materialization manifest.

No new semantic checker should duplicate the materializer's scientific logic in
this task.

---

## 11. Focused profile tests

Create:

```text
testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py
```

Use the existing F.6.1 generic-profile test as an architectural pattern only.

Tests must use synthetic temporary files and the existing collector. Do not
invoke the real materializer builders or farm artifacts.

At minimum test:

1. profile schema version is the existing v4 profile schema;
2. validation-profile identifier is exact;
3. collection mode is `generic_artifacts`;
4. top-level setting inventory is exactly canonical five in the frozen order;
5. global artifact declarations are exactly the five required JSONs in Section 8;
6. artifact-scoped `settings` declarations are exactly `[]`;
7. required analysis/source commit is exactly
   `141a3d04f9e5d07be21dba14e0e63212c3990bf1`;
8. exact allowed committed files are only the new profile and focused profile test;
9. allowed non-analysis prefix is only `docs/memory/`;
10. a complete synthetic package returns `returncode == 0`, bundle
    `complete == true`, and contains each of the five global payloads exactly
    once plus `manifest.json`, `source_state.txt`, and `source_checks.txt`;
11. each of the five synthetic JSON records is reported `status == exists` and
    `json_status == valid`;
12. the five canonical setting manifest entries contain no declared artifacts;
13. missing any required global artifact produces bundle `complete == false`;
14. malformed JSON for any required global artifact produces bundle
    `complete == false`;
15. required source commit not being an ancestor produces
    `required_analysis_commit_not_present`;
16. an unexpected committed scientific/runtime/materializer file after
    `141a3d...` produces
    `unexpected_committed_files_after_required_analysis_commit`;
17. allowed profile/test paths plus `docs/memory/` changes do not become
    unexpected committed files;
18. an already-existing output ZIP is refused rather than overwritten.

Do not test scientific equality numerics here; that belongs to the unchanged
materializer and comparator tests.

---

## 12. Negative and regression boundaries

The implementation must demonstrate that:

- the generic collector is reused unchanged;
- the materializer is unchanged;
- comparator source is unchanged;
- scientific/runtime source is unchanged;
- accepted authority constants are unchanged;
- the profile cannot silently package an undeclared extra artifact;
- required artifacts cannot be omitted without `complete=false`;
- malformed required JSON cannot yield `complete=true`;
- an unreviewed committed materializer/scientific change after `141a3d...`
  causes the source-provenance gate to fail;
- docs/memory changes remain non-analysis provenance only.

Do not loosen any existing collector source-identity rule to make this profile
pass.

---

## 13. Required local deterministic checks

Discover the local interpreter and use it consistently as `<PYTHON>`.

Run:

```bash
<PYTHON> -B -m py_compile \
  testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py

<PYTHON> -B -m unittest \
  testing.test_pion_hgcer_validation_bundle_profile_f4_refresh2 -v

<PYTHON> -B -m unittest \
  testing.test_collect_pion_hgcer_validation_bundle -v

<PYTHON> -B -m unittest \
  testing.test_materialize_method_a_current_baseline_authority -v

<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/check_memory_health.py --root .
<PYTHON> -B tools/memory_bootstrap.py --root . --json
<PYTHON> -B -m unittest testing.test_memory_health -v

git -c core.safecrlf=false diff --check
```

Report exact pass/fail counts and skips.

A pre-existing `CURRENT.md` soft-size warning remains nonfatal if it is the only
memory-health warning and all required checks otherwise pass.

These checks do not establish farm/runtime validation.

---

## 14. Durable memory/status updates

Create the phase record:

```text
docs/memory/phases/f4-refresh2-validation1-farm-materialization-bundle-profile.md
```

During Codex implementation, record:

```text
F.4.Refresh.2.Validation.1 — ACTIVE
```

pending independent ChatGPT actual-diff/source-provenance review.

Do **not** advance it to `SOURCE REVIEWED`; that requires the later independent
ChatGPT review and final-pre-push reconciliation gate.

Preserve:

```text
F.4.Refresh.2         — SOURCE REVIEWED
F.4.Refresh.2.Fix.1   — SOURCE REVIEWED
F.4.Refresh.1         — CLOSED / RUNTIME VALIDATED
F.4.Refresh.1.Fix.1   — CLOSED / RUNTIME VALIDATED
historical F.1–F.6.2  — CLOSED / RUNTIME VALIDATED
F.6.3                 — SOURCE REVIEWED
E.8.4                 — SOURCE REVIEWED
E.8                   — ACTIVE
final E.8             — BLOCKED
F.6.4                 — BLOCKED
lifecycle hook        — BLOCKED / DEFERRED
```

Update the existing F.4.Refresh.2 phase, Phase-F roadmap, and roadmap status only
as needed to record that the pushed materializer is now the pinned source for an
ACTIVE profile/bundle candidate.

Do not create farm evidence.

Regenerate `docs/memory/manifest.json` whenever versioned memory changes.

---

## 15. CURRENT.md exact NEXT after implementation

Set exactly one ordinary NEXT:

```text
NEXT — independent ChatGPT actual-diff/source-provenance review of the F.4.Refresh.2.Validation.1 farm materialization bundle/profile candidate.
```

Do not add another NEXT.

Do not say the profile is source reviewed before that independent review.

---

## 16. Farm-validation boundary

This task must not:

- run the materializer on farm artifacts;
- run the collector on the farm;
- create a farm ZIP;
- invoke ROOT/PyROOT;
- run `main.py`;
- rerun full analysis;
- update accepted F.2/F.3/F.4 authority;
- re-pin F.5/F.6.3/E.8.4 authority;
- promote Method A.

The later farm operation is a detached pure-Python materialization plus
bundle-collection gate, not a full scientific rerun and not a production run.

After this source candidate passes:

1. independent ChatGPT actual-diff review;
2. the established separate final-pre-push reconciliation task;
3. user-controlled commit/push;
4. independent pushed-state review;

only then may ChatGPT provide the exact `tcsh` farm command.

That later command must use the canonical farm paths from repository memory,
must not disturb the ordinary farm checkout, and should use a detached worktree
at the reviewed profile/bundle commit when needed for provenance isolation.
The exact farm command is deliberately **not** part of this task.

---

## 17. Required diff audit

Before stopping, run:

```bash
git status --short
git -c core.safecrlf=false diff --stat
git -c core.safecrlf=false diff --no-ext-diff
```

Audit every intended new/untracked file with:

```bash
git diff --no-index /dev/null <path> || true
```

The `|| true` exception is allowed only for `git diff --no-index` returning 1
when it displays a difference. Do not suppress test or implementation failures.

Expected candidate inventory is exactly:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization.md
docs/memory/phases/f4-refresh2-validation1-farm-materialization-bundle-profile-task-contract.md
docs/memory/phases/f4-refresh2-validation1-farm-materialization-bundle-profile.md
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/roadmap/STATUS.md
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py
```

If an allowed existing memory file needs no actual change, it need not appear in
the final diff. No path outside the allowlist may appear.

Temporary `kaonlt_review*.diff` files remain outside the candidate inventory.

---

## 18. Required cumulative review bundle

Create one fresh temporary repository-root review bundle:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

It must contain byte-faithfully:

1. branch;
2. committed HEAD;
3. `git status --short`;
4. complete cumulative `git diff --stat`;
5. complete tracked
   `git -c core.safecrlf=false diff --no-ext-diff`;
6. complete `git diff --no-index /dev/null ...` output for every intended
   new/untracked candidate file;
7. exact local validation commands/results;
8. final changed-path inventory.

The review bundle is temporary, must remain untracked, and is not part of the
candidate change set.

Return the fresh bundle for independent ChatGPT review.

---

## 19. Acceptance criteria

The local candidate is acceptable for independent source review only if:

1. committed starting HEAD remains exactly
   `141a3d04f9e5d07be21dba14e0e63212c3990bf1`;
2. only the profile, focused profile test, and allowed memory/history change;
3. materializer/comparator/collector/scientific/runtime source remain byte
   unchanged;
4. the new profile reuses `generic_artifacts`;
5. canonical-five top-level settings remain exact;
6. the profile declares exactly the five deterministic F.4.Refresh.2 output
   JSONs and no setting-scoped artifacts;
7. every declared artifact is required;
8. source identity pins the pushed materializer source `141a3d...`;
9. only profile/test are allowed committed files after that source, plus
   `docs/memory/`;
10. focused positive, negative, and provenance tests pass;
11. existing collector and materializer tests pass;
12. memory manifest/check/bootstrap/health checks pass subject only to known
    nonfatal soft-size warning;
13. `git diff --check` passes;
14. F.4.Refresh.2 remains `SOURCE REVIEWED`;
15. F.4.Refresh.2.Validation.1 remains `ACTIVE` pending independent ChatGPT
    review;
16. no farm/runtime/authority/promotion claim is made;
17. one fresh cumulative review bundle is produced.

---

## 20. Hard stop

After implementation, deterministic local checks, warranted memory/history
updates, manifest regeneration, diff audit, and creation of the fresh review
bundle:

**STOP.**

Do not commit or push.

Do not run the farm materializer.

Do not run the farm collector.

Do not create the farm ZIP.

Do not update accepted authority.

Do not modify F.5/F.6.3/E.8.4 authority.

Do not run full analysis.

Do not promote Method A.

NEXT is independent ChatGPT actual-diff/source-provenance review of the
F.4.Refresh.2.Validation.1 farm materialization bundle/profile candidate.
