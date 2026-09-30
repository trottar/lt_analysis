# KaonLT F.4.Refresh.2.Fix.1 — Materialization Manifest Integrity Repair

## 1. Purpose

Repair two narrow manifest/publication-integrity defects found by independent
ChatGPT actual-diff/source-runtime-path review of:

```text
kaonlt_review(20260930-104453).diff
```

The reviewed cumulative F.4.Refresh.2 candidate is otherwise correctly scoped
and follows the approved detached builder chain, but it must not advance to
final pre-push reconciliation until these two issues are repaired.

### Finding A — provenance-inclusive stage fingerprints are mislabeled as scientific fingerprints

The current materializer records:

```text
scientific_payload_fingerprint = body.get("fingerprint")
```

for accepted and candidate F.2/F.3/F.4 artifacts.

That label is not correct for the existing stage `fingerprint` fields:

- F.2 `fingerprint_inputs` include current F.1 content/source-file provenance;
- F.3 `fingerprint_inputs` include candidate F.2 source-file identity plus
  current F.1 content/source-file provenance;
- F.4 `fingerprint_inputs` include candidate F.3 source-file/map/artifact
  identity plus current F.1 content/source-file provenance.

Those fingerprints are useful stage/payload identities, but they are not the
provenance-excluded scientific equality hashes used by F.4.Refresh.1.

The materialization manifest must preserve the distinction between:

```text
stage identity / provenance-bound fingerprint
```

and:

```text
scientific equality result
```

Scientific equality remains owned only by the explicit `scientific_gate`
comparison.

### Finding B — `--overwrite` can leave a stale complete manifest over a partially replaced candidate set

The current publication order correctly writes the new manifest last, but when
`--overwrite` is used and an older manifest already exists, that old
`complete=true` manifest remains visible while the four candidate files are
being replaced.

If publication fails after one or more candidate replacements but before the
new manifest replacement, the old manifest can remain and falsely describe the
now-mixed output set as complete.

For overwrite publication, an existing manifest must therefore cease to be a
visible completion marker before any candidate target is replaced. If
publication fails after candidate replacement begins, no manifest may remain
at the final manifest path.

These are detached manifest/provenance-safety repairs only. They do not change
F.2/F.3/F.4 scientific calculations, accepted authorities, candidate numerical
content, or Method-A production ownership.

---

## 2. Exact starting state

Required branch:

```text
test
```

Required committed HEAD:

```text
08f8c5be84ab7a54278e8c623eae5d9d3c43d938
```

Required cumulative local candidate:

```text
kaonlt_review(20260930-104453).diff
```

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
```

Hard requirements:

- committed HEAD remains exactly
  `08f8c5be84ab7a54278e8c623eae5d9d3c43d938`;
- the existing cumulative F.4.Refresh.2 candidate must match the reviewed
  bundle apart from this repair;
- pre-existing `kaonlt_review*.diff` files remain temporary/untracked and must
  not be removed, rewritten, or staged;
- root `AGENTS.md` and `.codex/` remain local-only/untracked;
- do not reset, stash, clean, discard, commit, push, or run the farm;
- if unrelated tracked edits exist, STOP and report the blocker.

---

## 3. Mandatory repository-memory startup

Read in exact order:

1. root `AGENTS.md`;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read only task-relevant records/source:

- `docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization-task-contract.md`
- `docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization.md`
- `testing/materialize_method_a_current_baseline_authority.py`
- `testing/test_materialize_method_a_current_baseline_authority.py`
- `testing/compare_method_a_current_baseline_authority.py`
- `src/cuts/pion_hgcer_method_a_acceptance_representation.py`
- `src/cuts/pion_hgcer_method_a_acceptance_map.py`
- `src/cuts/pion_hgcer_method_a_parent_preserving_correction.py`

Do not redesign adjacent architecture.

---

## 4. Review findings already accepted

The following parts of the F.4.Refresh.2 candidate passed source inspection and
must remain unchanged:

- exact canonical-five F.1 inventory;
- comparison raw-SHA gate;
- accepted F.2/F.3/F.4 raw-SHA continuity;
- current F.1 raw-SHA and provenance continuity;
- public F.2 -> F.3 -> F.4 builder dependency order;
- deterministic candidate F.2/F.3 writer SHA reproduction;
- exact independent F.2/F.3 scientific equality;
- in-memory zero-head candidate-F.3 diagnostic authority override;
- exact F.4 detailed comparison/summary reproduction;
- distinct non-authoritative candidate filenames;
- task-local temporary directory;
- public F.2/F.3/F.4 writers;
- candidate byte verification;
- no accepted-authority mutation;
- no scientific/runtime source mutation;
- no production application or Method-A promotion.

Do not reopen those parts unless this repair exposes a concrete blocker.

---

## 5. Allowed substantive edits

Only:

```text
testing/materialize_method_a_current_baseline_authority.py
testing/test_materialize_method_a_current_baseline_authority.py
docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization.md
docs/memory/manifest.json
```

Create exactly one new repair contract:

```text
docs/memory/phases/f4-refresh2-fix1-materialization-manifest-integrity-task-contract.md
```

No other file may change.

In particular, keep byte-identical:

```text
testing/compare_method_a_current_baseline_authority.py
testing/test_compare_method_a_current_baseline_authority.py
docs/memory/CURRENT.md
docs/memory/evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md
docs/memory/phases/f4-refresh1-current-baseline-authority-comparison.md
docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization-task-contract.md
docs/memory/phases/phase-f6-method-a-production-promotion.md
docs/memory/roadmap/STATUS.md
```

and every file under `src/`.

If a scientific/runtime source file must change, STOP.

---

## 6. Fix A — manifest fingerprint semantics

### 6.1 Remove the misleading key

Do not emit:

```text
scientific_payload_fingerprint
```

for accepted or candidate F.2/F.3/F.4.

Do not create a new hash and call it a scientific fingerprint in this repair.

Scientific equality is already represented by:

```text
scientific_gate.f2_scientific_payload_match
scientific_gate.f3_scientific_payload_match
scientific_gate.f4_scientific_payload_match
scientific_gate.first_changed_stage
```

and remains an exact provenance-excluded comparison.

### 6.2 Persist the real source-owned stage fingerprint names

For `accepted_inputs` and `candidate_outputs`, persist stage fingerprints using
their actual semantics:

F.2:

```text
representation_fingerprint = representation["fingerprint"]
algorithm_fingerprint = representation["algorithm_fingerprint"]
```

F.3:

```text
map_fingerprint = acceptance_map["fingerprint"]
algorithm_fingerprint = acceptance_map["algorithm_fingerprint"]
accepted_basis = acceptance_map["accepted_basis"]
```

F.4:

```text
correction_fingerprint = correction["fingerprint"]
```

Keep existing:

```text
artifact_fingerprint
raw_sha256
```

and the already-required F.4 inherited F.3 identities.

Do not rename source-owned fingerprint fields inside the candidate artifacts.
This change is only the materialization manifest's labeling.

### 6.3 Scientific/provenance separation regression

The focused test must explicitly demonstrate that a provenance-bound stage
fingerprint is allowed to differ between accepted and candidate artifacts while
the corresponding:

```text
scientific_payload_match
```

is true.

At minimum cover F.2 and F.3.

The manifest must never imply that the provenance-bound stage fingerprint is
the scientific equality proof.

---

## 7. Fix B — overwrite completion-marker safety

### 7.1 Completion-marker invariant

The final manifest path is the completion marker for the five-file detached
materialization set.

Before replacing any of:

```text
comparison
F.2 candidate
F.3 candidate
F.4 candidate
```

under `--overwrite`, an already-existing final manifest must no longer be
visible at its final target path.

Do this only after:

- all in-memory scientific/provenance gates pass;
- all temporary candidate files are written;
- all temporary byte/hash checks pass;
- the new temporary manifest is built and verified.

### 7.2 Path-bounded handling

Any handling of the previous manifest must be restricted to the exact final
manifest target and the task-created temporary directory.

No wildcard deletion.

A valid implementation is:

1. atomically move the existing manifest target into a backup path inside the
   task-created temporary directory;
2. begin candidate publication;
3. publish the new manifest last.

Equivalent path-bounded semantics are acceptable.

### 7.3 Failure behavior

Track whether candidate publication has begun.

If failure occurs **before** any candidate target has been replaced:

- restoring the previous manifest is permitted and preferred if it had been
  moved aside.

If failure occurs **after** any candidate target has been replaced:

- do **not** restore the previous manifest;
- ensure the final manifest target is absent;
- allow the partial candidate files to remain visibly incomplete without a
  completion marker;
- do not attempt broad rollback or cleanup of the output directory.

The task-created temporary directory may be cleaned normally.

Never delete or rewrite input artifacts.

### 7.4 Non-overwrite behavior

Preserve existing non-overwrite behavior:

- any pre-existing target rejects before temporary writes;
- no final output is published.

---

## 8. Required focused tests

Preserve all existing F.4.Refresh.2 tests and add/strengthen tests for:

1. accepted F.2 manifest row uses `representation_fingerprint`, not
   `scientific_payload_fingerprint`;
2. candidate F.2 manifest row uses `representation_fingerprint`;
3. accepted/candidate F.3 rows use `map_fingerprint`;
4. F.4 rows use `correction_fingerprint`;
5. `scientific_payload_fingerprint` appears nowhere in the materialization
   manifest;
6. F.2 provenance-bound representation fingerprints may differ while
   `f2_scientific_payload_match == true`;
7. F.3 provenance-bound map fingerprints may differ while
   `f3_scientific_payload_match == true`;
8. successful `--overwrite` still publishes the new manifest last;
9. when a prior complete manifest exists and overwrite publication fails after
   at least one candidate replacement, the final manifest path is absent;
10. the failed overwrite never restores the stale complete manifest over the
    mixed candidate set;
11. a failure before candidate publication may restore the prior manifest;
12. input files remain byte-identical through success and failure cases.

Do not weaken any existing test.

---

## 9. Local deterministic validation

Run:

```bash
python -B -m py_compile \
  testing/materialize_method_a_current_baseline_authority.py \
  testing/test_materialize_method_a_current_baseline_authority.py

python -B -m unittest \
  testing.test_materialize_method_a_current_baseline_authority -v

python -B -m unittest \
  testing.test_compare_method_a_current_baseline_authority -v

python -B -m unittest \
  testing.test_pion_hgcer_method_a_acceptance_representation -v

python -B -m unittest \
  testing.test_pion_hgcer_method_a_acceptance_map -v

python -B -m unittest \
  testing.test_pion_hgcer_method_a_parent_preserving_correction -v

python -B -m unittest \
  testing.test_f6_3_parallel_full_procedure_method_a -v

python -B tools/update_memory_manifest.py --root . --write
python -B tools/update_memory_manifest.py --root . --check
python -B tools/check_memory_health.py --root .
python -B tools/memory_bootstrap.py --root . --json
python -B -m unittest testing.test_memory_health -v

git -c core.safecrlf=false diff --check
```

Report exact skips.

No ROOT/PyROOT, farm materialization, profile creation, packaging, or full
analysis.

---

## 10. Required phase-memory update

Append a Fix.1 section to:

```text
docs/memory/phases/f4-refresh2-current-baseline-candidate-materialization.md
```

Record:

```text
F.4.Refresh.2.Fix.1 — ACTIVE
```

pending independent review of the repaired cumulative candidate.

Record the independent review finding from:

```text
kaonlt_review(20260930-104453).diff
```

as two narrow manifest-integrity repairs:

1. provenance-inclusive F.2/F.3/F.4 stage fingerprints were incorrectly labeled
   `scientific_payload_fingerprint`;
2. overwrite publication could leave an old `complete=true` manifest visible
   over a partially replaced output set if publication failed.

Also correct the existing phase wording:

Instead of implying the uncommitted materializer is contained *at* committed
HEAD `08f8c5...`, state that the local candidate is **based on** committed
`test` HEAD:

```text
08f8c5be84ab7a54278e8c623eae5d9d3c43d938
```

Keep overall:

```text
F.4.Refresh.2 — ACTIVE
```

Do not change `CURRENT.md` or its existing exact NEXT.

Regenerate `docs/memory/manifest.json`.

---

## 11. Required diff audit

Before stopping:

```bash
git status --short
git -c core.safecrlf=false diff --stat
git -c core.safecrlf=false diff --no-ext-diff
```

Audit every intended new/untracked candidate file with:

```bash
git diff --no-index /dev/null <path> || true
```

The cumulative candidate may contain only the original reviewed F.4.Refresh.2
paths plus:

```text
docs/memory/phases/f4-refresh2-fix1-materialization-manifest-integrity-task-contract.md
```

and the permitted Fix.1 edits.

No scientific/runtime source may change.

---

## 12. Fresh cumulative review bundle

Create:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

at repository root.

It must contain byte-faithfully:

1. branch;
2. committed HEAD;
3. `git status --short`;
4. complete cumulative `git diff --stat`;
5. complete cumulative tracked diff;
6. complete no-index diffs for every intended untracked candidate file,
   including both F.4.Refresh.2 contracts;
7. exact validation commands/results;
8. final cumulative changed-path inventory.

Do not stage the review bundle.

---

## 13. Acceptance criteria

Fix.1 passes only if:

1. committed HEAD remains
   `08f8c5be84ab7a54278e8c623eae5d9d3c43d938`;
2. no scientific/runtime source changes;
3. original F.4.Refresh.2 builder/scientific gates remain unchanged;
4. `scientific_payload_fingerprint` is absent from the manifest;
5. F.2 manifest identity is labeled `representation_fingerprint`;
6. F.3 manifest identity is labeled `map_fingerprint`;
7. F.4 manifest identity is labeled `correction_fingerprint`;
8. exact scientific equality remains owned by `scientific_gate`;
9. provenance-bound fingerprint differences do not falsely fail scientific
   equality;
10. an old manifest cannot remain as a valid completion marker once overwrite
    candidate publication begins;
11. failure after partial overwrite publication leaves no final manifest;
12. non-overwrite preflight behavior remains fail-closed;
13. inputs remain unchanged;
14. public writers and candidate byte checks remain intact;
15. accepted authority constants remain untouched;
16. F.4.Refresh.2 and Fix.1 remain `ACTIVE`;
17. F.6.3/E.8.4 remain runtime blocked;
18. deterministic local suites pass subject to explicitly reported skips;
19. no farm run, profile creation, packaging, authority update, commit, push,
    full analysis, or Method-A promotion occurs;
20. a fresh cumulative review bundle is produced.

---

## 14. Hard stop

After the narrow Fix.1 implementation, focused/regression tests, phase-memory
update, manifest regeneration, diff audit, and fresh cumulative review bundle:

**STOP.**

Do not commit, push, run the farm materializer, create the validation profile,
package candidate artifacts, update accepted F.2/F.3/F.4/F.5/F.6.3 authority,
rerun the full analysis, or promote Method A.

NEXT remains independent ChatGPT actual-diff/source-runtime-path review of the
repaired cumulative F.4.Refresh.2 candidate.
