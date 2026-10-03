# KaonLT — E.8.4 Fix.5.7 page-manifest setting provenance / owner repair

## 1. Objective and exact starting HEAD

Repair the single source-verified integration mismatch that caused the tracked
E.8.4 Fix.5 owner to reject the freshly rendered page manifest with:

```text
stage = verify_artifacts
failure_reason = page_manifest_setting_invalid
```

Required branch:

```text
test
```

Required starting HEAD:

```text
cdc6ead47be3c46987a8418b796d07d87b80f829
```

Commit subject at that HEAD:

```text
Close KaonLT workflow hardening
```

This task is a narrow owner/checker integration repair. It must not change
analysis physics, renderer science, page IDs, page availability, profile scope,
collector semantics, candidate identities, artifact inventory, launcher
behavior, Method-A numerics, Method-B numerics, SIMC normalization, yields,
cuts, templates, priors, binning, efficiencies, acceptance, L/T separation, or
cross sections.

No farm run is authorized by this contract.

---

## 2. Mandatory startup and preflight

Read in exact order:

1. repository-root `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read only the task-relevant records/source:

- `docs/memory/MAINTENANCE.md`
- `docs/memory/CODEX.md`
- `docs/memory/TOOLS.md`
- `docs/memory/phases/e8-4-fix5-6-owner-farm-readiness-and-failure-provenance.md`
- `docs/memory/investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md`
- `testing/run_e8_4_fix5_left_lowe_plot_gate.py`
- `testing/test_run_e8_4_fix5_left_lowe_plot_gate.py`
- `src/cuts/rand_sub.py`
- `src/cuts/pion_hgcer_refinement_checkpoint.py`
- `src/cuts/full_background_subtraction_plots.py`
- directly relevant existing tests for those producer interfaces

Before editing, establish and report:

```bash
git branch --show-current
git rev-parse HEAD
git status --short --untracked-files=all
```

Hard requirements:

- branch is exactly `test`;
- HEAD is exactly `cdc6ead47be3c46987a8418b796d07d87b80f829`;
- preserve all unrelated worktree state;
- any unrelated tracked modification is a hard stop;
- do not reset, clean, stash, checkout over, commit, push, or run the farm;
- do not assume the ChatGPT-visible remote state proves the local worktree is clean.

---

## 3. Source-verified failure mechanism

Preserve this source finding.

### Runtime producer path

`src/cuts/rand_sub.py` sets the full-background page-manifest setting from:

```python
(histDict.get("pion_hgcer_refinement_checkpoint") or {}).get("setting")
```

The current checkpoint producer in
`src/cuts/pion_hgcer_refinement_checkpoint.py` serializes the setting with:

```text
kinematic_token
Q2
W
epsilon_setting
epsilon_filename_token
phi_setting
particle_type
```

`capture_full_background_subtraction_e8_2_render_state(...)` detaches that
setting without redefining its identity.

`build_full_background_subtraction_page_manifest_artifact(...)` calls
`_full_background_manifest_setting(...)`, which requires the four core identity
fields but copies the complete supplied setting mapping into the page manifest.

Therefore the real current page manifest legitimately carries additional
producer-owned provenance beyond the four core identity fields.

### Owner mismatch

`testing/run_e8_4_fix5_left_lowe_plot_gate.py::verify_pages(...)` currently
requires exact dictionary equality to only:

```python
{
    "kinematic_token": "Q4p4W2p74",
    "epsilon_filename_token": "lowe",
    "phi_setting": "Left",
    "particle_type": "kaon",
}
```

Exact equality rejects the valid current producer payload when the additional
`Q2`, `W`, and `epsilon_setting` fields are present.

### Test gap

`testing/test_run_e8_4_fix5_left_lowe_plot_gate.py::pages()` constructs only the
four-key owner subset. The owner success tests therefore do not exercise the
actual producer-to-owner setting interface.

This task repairs that integration mismatch only.

---

## 4. Scientific ownership and frozen source

### Scientific ownership

This is owner/checker orchestration and provenance validation only.

- Method A remains detached/non-production.
- Method B remains diagnostic/cross-check only and numerically excluded.
- The active `no_empirical_residual` profile remains unchanged.
- The unresolved absolute-SIMC unit blocker remains unchanged.
- Existing accepted F.6.3 / E.8.4 Left-low runtime evidence remains unchanged.
- Historical E.8.3 lineage remains distinct from current F.6.3/E.8.4 lineage.

### Frozen files

Do not modify these unless an unexpected source contradiction forces a hard
stop and a new contract:

```text
src/
run_Prod_Analysis.sh
testing/collect_pion_hgcer_validation_bundle.py
testing/test_collect_pion_hgcer_validation_bundle.py
testing/pion_hgcer_validation_bundle_profile_e8_4_fix5.json
```

The failure is not repaired by changing the renderer or page-manifest producer.

---

## 5. Allowed substantive files

Expected source/test changes:

```text
testing/run_e8_4_fix5_left_lowe_plot_gate.py
testing/test_run_e8_4_fix5_left_lowe_plot_gate.py
```

Allowed memory/contract changes when warranted:

```text
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/LEARNINGS.md
docs/memory/phases/e8-4-fix5-7-page-manifest-setting-provenance-owner-repair-task-contract.md
docs/memory/phases/e8-4-fix5-7-page-manifest-setting-provenance-owner-repair.md
docs/memory/investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md
docs/memory/manifest.json
```

Do not broaden the task merely because adjacent files are convenient to edit.

---

## 6. Required before/after behavior

### Before

A valid real page manifest produced from the current HGCer refinement checkpoint
can contain the producer-owned setting metadata:

```text
kinematic_token
Q2
W
epsilon_setting
epsilon_filename_token
phi_setting
particle_type
```

The owner rejects it because it compares the whole setting mapping for exact
equality with a four-key subset.

### After

The owner must validate the exact gate-relevant setting identity without
requiring the manifest mapping to contain only four keys.

At minimum require exact values for:

```text
kinematic_token = Q4p4W2p74
epsilon_filename_token = lowe
phi_setting = Left
particle_type = kaon
```

Also validate the current producer's semantic epsilon provenance:

```text
epsilon_setting = low
```

when exercising the real current producer interface.

Do not reject valid producer-owned `Q2` / `W` metadata merely because those
fields exist.

Do not mutate, strip, rewrite, or normalize the manifest setting inside the
owner. The owner is a verifier, not a producer.

A mismatch of any gate-relevant identity field must still fail closed with:

```text
page_manifest_setting_invalid
```

Do not weaken any other page-manifest gate.

---

## 7. Required regression that closes the missed interface

The prior four-key synthetic fixture is insufficient as the only positive case.

Add a deterministic regression that exercises the current producer-to-owner
setting shape. Prefer a real pure-Python producer path where practical:

1. construct or obtain a current-format HGCer refinement checkpoint setting
   containing the full producer-owned setting metadata;
2. feed that setting through the current
   `build_full_background_subtraction_page_manifest_artifact(...)` path, or an
   equivalently direct source-owned helper path;
3. populate the existing required page inventory;
4. prove `owner.verify_pages(...)` accepts that manifest;
5. prove the manifest retains its producer-owned extra metadata;
6. mutate each gate-relevant identity field and prove the owner rejects it with
   `page_manifest_setting_invalid`.

If importing the real producer in the owner test creates a genuine deterministic
dependency blocker, stop and report it. Do not silently replace this required
interface test with another four-key hand-built fixture.

The existing synthetic complete-owner test may remain, but its successful page
manifest must no longer represent a shape that contradicts the current runtime
producer.

---

## 8. Positive checks

Prove locally that:

- current full producer setting metadata is accepted by the owner;
- `kinematic_token`, `epsilon_filename_token`, `phi_setting`,
  `particle_type`, and semantic low epsilon remain strictly checked;
- extra current producer-owned metadata is preserved and tolerated;
- `verify_pages(...)` still returns the page count for a valid manifest;
- the complete synthetic owner success path still reaches collection and ZIP
  verification;
- status-sidecar success semantics are unchanged;
- success stdout remains exactly one ZIP path.

---

## 9. Negative checks

Prove fail-closed behavior for at least:

- wrong `kinematic_token`;
- wrong `epsilon_filename_token`;
- wrong `phi_setting`;
- wrong `particle_type`;
- wrong current semantic `epsilon_setting`;
- missing required core identity field;
- non-mapping/invalid `setting`;
- existing page-schema failure;
- nonempty renderer failures;
- missing/duplicate required page;
- wrong page scope / t identity / phi inventory.

Preserve the existing literal failure reasons unless a source contradiction
requires otherwise.

---

## 10. Regression checks

Preserve all Fix.5.6 owner safety behavior:

- source branch/HEAD/origin checks;
- known-dirt handling;
- no ordinary-checkout cleanup;
- clean detached collector source preflight;
- exact `./run_Prod_Analysis.sh -d 4p4 2p74` invocation;
- completion markers;
- freshness/refreshed-artifact checks;
- candidate hashes;
- strict JSON;
- PDF magic;
- page IDs/scopes/inventory;
- run summary;
- final clean detached collection;
- final collector source rechecks;
- ZIP source/scope/profile/hash/inventory verification;
- durable gate-status sidecar;
- empty stdout on failure;
- one ZIP path on success.

The gate-status sidecar remains outside the validation bundle.

---

## 11. Forbidden shortcuts

Do not:

- change `src/cuts/full_background_subtraction_plots.py` to strip metadata;
- change `src/cuts/rand_sub.py` to manufacture a four-key owner-specific setting;
- change the checkpoint schema to satisfy the owner;
- remove provenance fields from the manifest;
- weaken schema, renderer-failure, page-ID, page-scope, child-inventory, artifact,
  source, candidate, collector, or ZIP checks;
- catch and ignore `page_manifest_setting_invalid`;
- skip page-manifest verification;
- alter the validation profile or eight-artifact inventory;
- alter scientific/presentation source;
- rerun the farm;
- fabricate a replacement manifest;
- package the failed-gate artifacts as accepted evidence;
- interpret the failed-gate PDF scientifically.

---

## 12. Preserved execution chain and operational readiness

This task remains farm-bound after source review and synchronization, but does
not itself authorize farm execution.

Preserve the complete tracked chain:

```text
input authority / current candidate identities
-> run_e8_4_fix5_left_lowe_plot_gate.py owner preflight
-> clean detached collector source preflight
-> run_Prod_Analysis.sh -d 4p4 2p74
-> main.py / rand_sub.py analysis + renderer
-> page-manifest producer
-> owner completion/artifact/page verification
-> run-summary creation
-> clean detached generic collector
-> collector source rechecks
-> validation ZIP
-> owner ZIP verification
-> returned ZIP/status/log/PDF/manifest evidence
```

Source owners:

```text
owner/checker:
  testing/run_e8_4_fix5_left_lowe_plot_gate.py

owner regression:
  testing/test_run_e8_4_fix5_left_lowe_plot_gate.py

renderer/page-manifest producer:
  src/cuts/full_background_subtraction_plots.py   [frozen]

runtime handoff into producer:
  src/cuts/rand_sub.py                            [frozen]

setting authority:
  src/cuts/pion_hgcer_refinement_checkpoint.py   [frozen]

collector:
  testing/collect_pion_hgcer_validation_bundle.py [frozen]

launcher:
  run_Prod_Analysis.sh                            [frozen]
```

No new driver is required: the existing reviewed tracked owner already owns the
complete multi-step operation. The defect is inside its manifest-setting
verification contract.

Farm readiness remains blocked until:

1. Codex implements and passes deterministic checks;
2. ChatGPT reviews the actual diff;
3. the user commits/pushes;
4. ChatGPT performs pushed-state synchronization review.

Only then can a new narrow Left/lowe owner farm gate be considered.

---

## 13. Local validation

Discover the working repository Python as required by `docs/memory/TOOLS.md`;
represent it as `<PYTHON>`.

Run syntax checks for changed Python files.

Run at minimum:

```text
testing.test_run_e8_4_fix5_left_lowe_plot_gate
testing.test_full_background_subtraction_plots
testing.test_pion_hgcer_refinement_checkpoint
testing.test_collect_pion_hgcer_validation_bundle
testing.test_pion_hgcer_validation_bundle_profile_e8_4_fix5
```

Run any directly affected focused suites needed to prove the real producer-to-owner
regression.

Run:

```bash
git diff --check
```

No local test establishes ROOT/PyROOT, full `main.py`, procedure-PDF rendering,
farm filesystem behavior, real artifact freshness, ZIP delivery, PDF legibility,
numerical closure, signed cancellation, or scientific interpretation.

---

## 14. Memory health

If memory changes:

1. regenerate the versioned manifest;
2. check the manifest;
3. run the ordinary task-final health check.

Use:

```bash
<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/check_memory_health.py --root .
```

Report:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
warning classification: <none or blocking/nonblocking>
manifest check: PASS | FAIL
```

Do not use `--fail-on-warning` unless a warning is materially blocking or this
task explicitly becomes a zero-warning maintenance task.

---

## 15. Phase / CURRENT requirements

Create/update the Fix.5.7 phase record only with source-supported facts.

During implementation:

```text
ACTIVE
```

After deterministic local checks pass:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

Do not claim runtime acceptance.

Record the source-verified root cause succinctly:

- current runtime producer carries the full checkpoint setting mapping;
- owner exact-four-key equality rejected valid extra provenance;
- the old owner test fixture masked the interface mismatch;
- repair is verifier/test only.

Leave the scientific blockers unchanged, including absolute-SIMC provenance.

### Push-stable NEXT

At local completion, CURRENT's sole ordinary NEXT must remain the substantive
next gate, conditional on synchronization, e.g.:

```text
NEXT — after independent Fix.5.7 actual-diff review, user commit/push, and
ChatGPT pushed-state synchronization, run one narrow Q4p4W2p74 / Left / lowe
tracked-owner farm gate and review the fresh status/ZIP/log/PDF/manifest
evidence. No farm command is authorized before those synchronization gates.
```

Do not make commit/push itself the sole NEXT.

---

## 16. Diff audit and acceptance criteria

Before stopping, report:

- exact starting and ending local HEAD;
- branch;
- full worktree state;
- exact changed paths;
- exact new/untracked intended files;
- deterministic test commands/results;
- `git diff --check`;
- memory health/manifest results when applicable;
- whether unrelated state was preserved;
- whether any scientific or frozen file changed.

Create one complete temporary review bundle in the repository root if the diff
is not comfortably reviewable inline:

```text
kaonlt_review.diff
```

It must include:

- complete tracked diff for every intended changed tracked file;
- complete `git diff --no-index /dev/null ...` addition for every intended new
  file.

Do not stage files merely to make them reviewable.

Acceptance requires:

1. only the owner/checker integration mismatch is repaired;
2. the real producer setting shape passes;
3. wrong gate identity still fails closed;
4. no frozen/scientific source changes;
5. all deterministic checks pass;
6. memory remains materially consistent;
7. the actual diff is ready for independent ChatGPT review.

---

## 17. Execution authority and hard stop

Codex performs local edits and deterministic checks only.

Codex must not:

- commit;
- push;
- update remote refs;
- run the Jefferson Lab farm;
- claim ROOT/PyROOT/full-analysis/PDF/farm validation;
- provide a farm command as though the gate were already synchronized.

Hard stop after local implementation, deterministic checks, memory/manifest
checks, and complete review material.

Return `BLOCKED` instead of inventing a workaround if:

- local branch/HEAD do not match the required starting state;
- unexpected tracked work exists;
- the real source path contradicts the diagnosed mismatch;
- fixing the issue requires touching a frozen scientific/runtime producer;
- the required producer-to-owner regression cannot be implemented
  deterministically without a broader dependency change.
