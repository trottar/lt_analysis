# KaonLT E.8.4.Fix.4 — Validation-Bundle Profile Re-Pin Task Contract

## 1. Purpose

Re-pin the existing E.8.4/E.8.1 generic validation-bundle profile to the newly
pushed and independently reviewed E.8.4.Fix.4 analysis source before any new
Jefferson Lab farm validation is authorized.

Independent pushed-state review established that remote `test` is now exactly:

```text
6e2adf7a37ac9e79cad99242686804cf51701644
```

with subject:

```text
E8.4 Fix.4: version persisted alignment semantics
```

and parent:

```text
8d62dcfdf2fd8d08298c075940b4ab1b28d7f079
```

The pushed commit contains exactly the reviewed ten-path Fix.4 candidate and
does not contain either temporary `kaonlt_review*.diff` file.

The live validation profile still requires the older Fix.3 analysis source:

```text
29d7b7f9635db899939efeb3508e941e994e8928
```

and its focused profile test still pins the same value in `REVIEWED_SOURCE`.

The profile's committed-range allowlist remains intentionally narrow:

```text
testing/pion_hgcer_validation_bundle_profile_e8_1.json
testing/test_pion_hgcer_validation_bundle_profile_e8_1.py
```

with only:

```text
docs/memory/
```

accepted as a non-analysis prefix.

Because the newly pushed Fix.4 commit contains a substantive change to:

```text
src/cuts/pion_component_fits.py
```

and that path is not an allowed post-required-analysis exception, attempting to
collect a new Fix.4 validation bundle while the profile still requires
`29d7b7f...` would correctly fail source-provenance validation.

This task changes only the profile's reviewed analysis-source pin and the
matching focused test constant, plus warranted repository memory/history. It
does not change the collector, wrapper, artifact inventory, canonical-five
setting declaration, scientific/runtime source, accepted F.6.2 science, or any
analysis behavior.

---

## 2. Exact starting state

Required branch:

```text
test
```

Required committed starting HEAD:

```text
6e2adf7a37ac9e79cad99242686804cf51701644
```

Required commit subject:

```text
E8.4 Fix.4: version persisted alignment semantics
```

Required parent:

```text
8d62dcfdf2fd8d08298c075940b4ab1b28d7f079
```

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
```

Hard requirements:

- branch must be `test`;
- committed HEAD must be exactly
  `6e2adf7a37ac9e79cad99242686804cf51701644`;
- inspect the existing worktree before editing;
- temporary root `kaonlt_review*.diff` files may remain untracked and must not
  be staged, deleted, rewritten, or treated as candidate files;
- root `AGENTS.md` and `.codex/` are local-only/untracked and must not be
  changed or removed;
- if unrelated tracked edits already exist, STOP and report the blocker;
- do not reset, stash, clean, discard, commit, push, package, rerender, or run
  the Jefferson Lab farm.

---

## 3. Mandatory repository-memory startup

Read in this exact order:

1. root `AGENTS.md`;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read only task-relevant records/source, including:

- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/e8-4-fix4-alignment-cache-semantics.md`
- `docs/memory/phases/e8-4-fix4-final-pre-push-source-review-reconciliation-task-contract.md`
- `docs/memory/phases/e8-4-fix3-bundle-profile-repin.md`
- `testing/pion_hgcer_validation_bundle_profile_e8_1.json`
- `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`
- `testing/collect_pion_hgcer_validation_bundle.py`
- `testing/package_pion_hgcer_validation_bundle.tcsh`

Read additional files only if required to verify the exact profile/source
provenance contract.

---

## 4. Current reviewed state to preserve

Preserve these statuses:

- F.1 through F.6.2, including F.6.2.Fix.5 —
  **CLOSED / RUNTIME VALIDATED**
- E.8 — **ACTIVE**
- E.8.2 — **SOURCE REVIEWED**
- E.8.3 — **SOURCE REVIEWED**
- F.6.3 — **SOURCE REVIEWED**
- E.8.4 — **SOURCE REVIEWED**
- E.8.4.Fix.3 — **SOURCE REVIEWED**
- E.8.4.Fix.4 — **SOURCE REVIEWED**
- final E.8 — **BLOCKED**
- F.6.4 — **BLOCKED**
- lifecycle-hook dispatch — **BLOCKED / DEFERRED**

The fresh `Q4p4W2p74 / Left / lowe` E.8.4 gate remains **BLOCKED** by:

```text
f6_3_f4_shared_reproduction_failed:f3_fingerprint_input_content_mismatch
```

That blocker is not runtime-closed by the Fix.4 source push. A fresh narrow farm
gate is still required after source-provenance preparation is complete.

Create:

- E.8.4.Fix.4 validation-bundle profile re-pin — **ACTIVE** pending independent
  ChatGPT actual-diff/source-provenance review.

---

## 5. Pushed-state review facts

Record, where warranted, that independent GitHub pushed-state review established:

```text
remote test HEAD:
6e2adf7a37ac9e79cad99242686804cf51701644

subject:
E8.4 Fix.4: version persisted alignment semantics

parent:
8d62dcfdf2fd8d08298c075940b4ab1b28d7f079
```

The pushed commit contains exactly the reviewed Fix.4 candidate paths:

```text
docs/memory/CURRENT.md
docs/memory/evidence/e8-4-fix3-left-lowe-stale-alignment-cache-runtime-blocker.md
docs/memory/manifest.json
docs/memory/phases/e8-4-fix3-pion-alignment-determinism.md
docs/memory/phases/e8-4-fix4-alignment-cache-semantics-task-contract.md
docs/memory/phases/e8-4-fix4-alignment-cache-semantics.md
docs/memory/phases/e8-4-fix4-final-pre-push-source-review-reconciliation-task-contract.md
docs/memory/roadmap/STATUS.md
src/cuts/pion_component_fits.py
testing/test_pion_component_dynamic_alignment.py
```

No temporary review bundle is part of the pushed commit.

The live pushed source still contains the reviewed semantic-version gate:

```text
pion_component_dynamic_alignment_semantics/v2
```

and the live pushed focused alignment tests contain the reviewed stale-cache,
wrong/malformed-version, current-reuse, and stale-parent coverage.

This pushed-state review is source/provenance evidence only. It is not
ROOT/PyROOT, full-analysis, procedure-PDF, farm, runtime, production, or
Method-A promotion evidence.

---

## 6. Allowed substantive changes

Only:

```text
testing/pion_hgcer_validation_bundle_profile_e8_1.json
testing/test_pion_hgcer_validation_bundle_profile_e8_1.py
```

No other substantive production/test/source file may change.

Allowed existing memory/history files:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/e8-4-fix4-alignment-cache-semantics.md
docs/memory/roadmap/STATUS.md
```

Allowed new records:

```text
docs/memory/phases/e8-4-fix4-bundle-profile-repin-task-contract.md
docs/memory/phases/e8-4-fix4-bundle-profile-repin.md
```

If another substantive file is required, STOP and report the blocker.

---

## 7. Frozen files and interfaces

Do not modify:

```text
src/cuts/pion_component_fits.py
testing/test_pion_component_dynamic_alignment.py
testing/collect_pion_hgcer_validation_bundle.py
testing/package_pion_hgcer_validation_bundle.tcsh
src/cuts/full_background_subtraction_plots.py
src/cuts/rand_sub.py
src/cuts/pion_hgcer_method_a_acceptance_contract.py
src/cuts/pion_hgcer_method_a_acceptance_map.py
src/cuts/pion_hgcer_method_a_parent_preserving_correction.py
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
src/cuts/pion_component_subtraction.py
src/cuts/particle_subtraction.py
src/binning/calculate_yield.py
src/binning/ave_per_bin.py
src/main.py
run_Prod_Analysis.sh
src/utility/background_config.py
```

Also freeze:

- profile schema version;
- validation-profile name;
- collection mode;
- canonical five settings;
- global artifact declarations;
- setting artifact declarations;
- F.6.2 frozen JSON basename;
- `allowed_committed_files`;
- `allowed_non_analysis_path_prefixes`;
- collector source;
- wrapper source;
- all scientific artifacts;
- all production physics.

---

## 8. Required profile change

In:

```text
testing/pion_hgcer_validation_bundle_profile_e8_1.json
```

change exactly:

```text
source_identity.required_analysis_commit
```

from:

```text
29d7b7f9635db899939efeb3508e941e994e8928
```

to:

```text
6e2adf7a37ac9e79cad99242686804cf51701644
```

No other profile field may change.

In:

```text
testing/test_pion_hgcer_validation_bundle_profile_e8_1.py
```

change exactly:

```python
REVIEWED_SOURCE = "29d7b7f9635db899939efeb3508e941e994e8928"
```

to:

```python
REVIEWED_SOURCE = "6e2adf7a37ac9e79cad99242686804cf51701644"
```

No other test constant, artifact declaration, setting inventory, allowlist, or
test behavior may change.

---

## 9. Required before/after behavior

### Before

The profile requires:

```text
29d7b7f9635db899939efeb3508e941e994e8928
```

The current pushed analysis source is:

```text
6e2adf7a37ac9e79cad99242686804cf51701644
```

The committed range includes the substantive Fix.4 source
`src/cuts/pion_component_fits.py`, which is outside the profile's allowed
post-required-analysis exceptions. The collector should therefore fail closed
rather than silently package the new analysis source under the old pin.

### After

The profile requires exactly:

```text
6e2adf7a37ac9e79cad99242686804cf51701644
```

The focused test pins the same reviewed source.

A later profile/bundle commit may differ from that required analysis source only
through the already-frozen exact exceptions:

```text
testing/pion_hgcer_validation_bundle_profile_e8_1.json
testing/test_pion_hgcer_validation_bundle_profile_e8_1.py
docs/memory/
```

The newly pushed Fix.4 source is therefore the explicit reviewed analysis
authority rather than an unexpected committed source change.

---

## 10. Required deterministic local validation

Use the actual repository-local Python interpreter.

Run at minimum:

```bash
python -B -m py_compile \
  testing/test_pion_hgcer_validation_bundle_profile_e8_1.py \
  testing/collect_pion_hgcer_validation_bundle.py

python -B -m json.tool \
  testing/pion_hgcer_validation_bundle_profile_e8_1.json \
  >/dev/null

python -B -m unittest \
  testing.test_pion_hgcer_validation_bundle_profile_e8_1 -v

python -B -m unittest \
  testing.test_collect_pion_hgcer_validation_bundle -v

python -B -m unittest \
  testing.test_pion_component_dynamic_alignment -v

python -B tools/update_memory_manifest.py --root . --write
python -B tools/update_memory_manifest.py --root . --check
python -B tools/check_memory_health.py --root .
python -B tools/memory_bootstrap.py --root . --json
python -B -m unittest testing.test_memory_health -v

git -c core.safecrlf=false diff --check
```

If PyROOT-dependent tests skip, report the exact count/reason.

These are local deterministic checks only. Do not claim farm validation.

---

## 11. Required durable memory/history

### 11.1 `docs/memory/phases/e8-4-fix4-alignment-cache-semantics.md`

Append only warranted pushed-state provenance:

- pushed commit
  `6e2adf7a37ac9e79cad99242686804cf51701644`;
- independent pushed-state review PASS;
- exact source/provenance boundary;
- profile re-pin remains a distinct subsequent gate.

Do not upgrade runtime status.

### 11.2 New profile-repin phase record

Create:

```text
docs/memory/phases/e8-4-fix4-bundle-profile-repin.md
```

Record:

- starting HEAD `6e2adf7...`;
- exact old/new required-analysis source;
- only profile/test substantive changes;
- canonical five/artifact inventory unchanged;
- allowlists unchanged;
- collector/wrapper/scientific source unchanged;
- status **ACTIVE** pending independent ChatGPT actual-diff/source-provenance
  review;
- local test results accurately;
- no farm/runtime claim.

### 11.3 `CURRENT.md`

Update stale pre-push state to reflect that the user-controlled Fix.4 push has
occurred and independent pushed-state review passed.

Set the profile re-pin to **ACTIVE** pending independent actual-diff review.

The exact ordinary NEXT after implementation must be:

```text
NEXT — independent ChatGPT actual-diff/source-provenance review of the E.8.4.Fix.4 validation-bundle profile re-pin candidate.
```

Do not add a second NEXT.

### 11.4 `roadmap/STATUS.md`

Record the pushed Fix.4 source identity and profile-repin ACTIVE gate without
changing scientific phase ownership.

### 11.5 Manifest

Regenerate `docs/memory/manifest.json`.

Do not change `MEMORY.md`, `USER.md`, `LEARNINGS.md`, `CODEX.md`, `TOOLS.md`,
decisions, evidence, or unrelated phase records unless a concrete blocker is
found; if so STOP before expanding scope.

---

## 12. Required diff audit

Before stopping:

```bash
git status --short
git -c core.safecrlf=false diff --stat

git -c core.safecrlf=false diff -- \
  testing/pion_hgcer_validation_bundle_profile_e8_1.json \
  testing/test_pion_hgcer_validation_bundle_profile_e8_1.py \
  docs/memory/CURRENT.md \
  docs/memory/manifest.json \
  docs/memory/phases/e8-4-fix4-alignment-cache-semantics.md \
  docs/memory/roadmap/STATUS.md
```

Audit new files:

```bash
git diff --no-index /dev/null \
  docs/memory/phases/e8-4-fix4-bundle-profile-repin-task-contract.md || true

git diff --no-index /dev/null \
  docs/memory/phases/e8-4-fix4-bundle-profile-repin.md || true
```

Verify that the substantive profile/test diff contains only the two source-pin
replacements described in Section 8.

---

## 13. Required cumulative review bundle

Create one fresh temporary repository-root review bundle:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

It must contain byte-faithfully:

1. branch;
2. committed HEAD;
3. `git status --short`;
4. complete `git diff --stat`;
5. complete tracked `git -c core.safecrlf=false diff --no-ext-diff`;
6. complete `git diff --no-index /dev/null ...` output for both intended new
   phase/contract records;
7. exact local validation commands/results;
8. final changed-path inventory.

Do not normalize or rewrite raw diff text.

Temporary existing `kaonlt_review*.diff` files must remain untracked and must
not be part of the intended candidate.

Stop for independent ChatGPT review.

---

## 14. Acceptance criteria

The implementation candidate passes this task only if:

1. committed HEAD remains exactly
   `6e2adf7a37ac9e79cad99242686804cf51701644`;
2. only the two allowed substantive profile/test files change;
3. `required_analysis_commit` becomes exactly `6e2adf7...`;
4. `REVIEWED_SOURCE` becomes exactly `6e2adf7...`;
5. no other profile field changes;
6. canonical-five settings remain unchanged;
7. artifact inventory remains unchanged;
8. committed-file allowlist remains unchanged;
9. non-analysis prefix allowlist remains exactly `docs/memory/`;
10. collector/wrapper source remains unchanged;
11. Fix.4 analysis/test source remains unchanged;
12. all scientific/runtime interfaces remain frozen;
13. local profile/collector tests pass;
14. dynamic-alignment regression remains passing subject to explicit PyROOT
    skips;
15. memory checks pass, allowing only the known nonfatal CURRENT soft-size
    warning when exit code is zero;
16. Fix.4 remains **SOURCE REVIEWED**;
17. profile re-pin is **ACTIVE** pending independent ChatGPT review;
18. no farm command, bundle collection, commit, push, or production promotion
    occurs;
19. a fresh byte-faithful cumulative review bundle is produced.

---

## 15. Forbidden shortcuts

Do not:

- add `src/cuts/pion_component_fits.py` to `allowed_committed_files`;
- broaden any allowlist;
- make the collector ignore substantive source changes;
- leave `required_analysis_commit` at `29d7b7f...`;
- pin to a future/uncommitted source;
- modify the collector or wrapper;
- modify Fix.4 source/test behavior;
- change canonical settings or artifacts;
- alter scientific calculations or runtime behavior;
- regenerate accepted F.1--F.6.2 artifacts;
- run the farm;
- collect a farm bundle;
- commit or push.

---

## 16. Hard stop

After the exact profile/test source-pin replacements, warranted memory/history
updates, deterministic local checks, manifest regeneration, diff audit, and
creation of the fresh cumulative review bundle:

**STOP.**

NEXT is independent ChatGPT actual-diff/source-provenance review.

No final pre-push reconciliation, Git handoff, or farm command is authorized by
this task.
