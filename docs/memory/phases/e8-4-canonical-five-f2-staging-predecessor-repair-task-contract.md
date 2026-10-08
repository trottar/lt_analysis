# E.8.4 canonical-five F.2 staging-predecessor migration repair — task contract

## 1. Task class and starting authority

This is a **narrow source-changing owner-state migration repair** under
`docs/memory/CODEX.md`.

Authoritative repository:

`https://github.com/trottar/lt_analysis/tree/test`

Exact pushed starting identity:

```text
branch: test
HEAD: 946b5512bb00d22bb446bd96d6ab9ee35f5f989b
origin/test: 946b5512bb00d22bb446bd96d6ab9ee35f5f989b
```

The previous F.2/F.3/F.4 scientific-equivalence repair remains source-reviewed
and pushed. This task owns only the **first failed invariant** from the fresh
farm attempt: candidate staging rejected one pre-existing F.2 target whose exact
bytes were not in the owner's recognized predecessor set.

Do not run the Jefferson Lab farm.
Do not stage, commit, push, update refs, reset, clean, stash, or overwrite
unrelated local work.

---

## 2. Required startup and hard start gate

Before editing, read in exact order:

1. `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read:

- `docs/memory/MAINTENANCE.md`
- `docs/memory/CODEX.md`
- `docs/memory/COMMUNICATION.md`
- `docs/memory/TOOLS.md`
- `docs/memory/phases/e8-4-canonical-five-f2-handoff-equivalence-repair.md`
- `docs/memory/phases/e8-4-canonical-five-f2-handoff-equivalence-repair-task-contract.md`
- `docs/memory/evidence/e8-4-canonical-five-second-full-run-provenance-blocker-2026-10-05.md`
- `docs/memory/roadmap/STATUS.md`

Then inspect:

```text
testing/run_e8_4_fix5_canonical_five_plot_gate.py
testing/test_run_e8_4_fix5_canonical_five_plot_gate.py
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
```

The runtime source is audit-only and frozen for this task.

Report:

```bash
git branch --show-current
git rev-parse HEAD
git rev-parse origin/test
git status --short --untracked-files=all
```

Required branch/HEAD/origin are exactly the identity above.

Preserve unrelated local state byte-for-byte. Hard stop if its ownership is
ambiguous or this repair cannot be completed without resetting, cleaning,
stashing, deleting, or staging unrelated paths.

---

## 3. Fresh failed-gate evidence

Fresh owner attempt stem:

```text
KaonLT_E8_4_Fix5_CanonicalFive_Q4p4W2p74_20261007-231233_3914085
```

Farm-evaluated source:

```text
946b5512bb00d22bb446bd96d6ab9ee35f5f989b
```

Supplied gate-status fields:

```text
status = failed
stage = candidate_staging
failure_reason =
  unknown_candidate_target:
  Q4p4W2p74_kaon_pion-background_hgcer-method-a-acceptance-representation-current-baseline-candidate.json

analysis_started = false
analysis_completed = false
artifact_verification_completed = false
collection_completed = false
zip_verification_completed = false
```

Therefore no child analysis, artifact acceptance, collection, ZIP acceptance or
canonical-five runtime claim follows from this attempt.

### Existing target observed on the farm

Basename:

```text
Q4p4W2p74_kaon_pion-background_hgcer-method-a-acceptance-representation-current-baseline-candidate.json
```

Exact observed raw SHA-256:

```text
182433ccd2d13d3f7af80b963e69190ad8f2eb4e295c8216b36585e8ee7409b2
```

Observed wrapper/body identities:

```text
schema_version:
  pion_hgcer_method_a_acceptance_representation_artifact/v1

non_authoritative = true
probe_only = true

representation_fingerprint:
  89976fd7a94a608be81a9bf0e177c5c59b5b5d6e6f39b5c704c0e719ce877585

algorithm_fingerprint:
  2ca2cc03ae9fae90c50c620f267d467a213a819c1e724a1b731bd88e2a0184bd

artifact_fingerprint:
  c63e7461240622fc64e9a2d566c687deb5fb620b54b0d02f99dc5485a8d2d400
```

Its provenance block was non-identifying:

```json
{
  "generated_at_utc": null,
  "git_head": null,
  "git_status_short": null,
  "input_paths": {}
}
```

This object must **not** be promoted to accepted authority or described as a
known historical authority.

### Reviewed F.2 candidate

Exact reviewed candidate raw SHA-256:

```text
2fa715b2b5c2d5077e38416fb9eceb814e013e72f308c43526fddff36ddd962e
```

Reviewed representation fingerprint:

```text
e9df6e09c31dda4d7d724792ce3c40125e27530c4c76e0dbe7f5b9861d0e6216
```

Reviewed algorithm fingerprint:

```text
2ca2cc03ae9fae90c50c620f267d467a213a819c1e724a1b731bd88e2a0184bd
```

Reviewed artifact fingerprint:

```text
87d9a8780cb2de9f2830afd9e388064151f796654ddb5c4a8e1938443806a0c6
```

The reviewed materialization hash itself was rechecked on farm during the
read-only comparison:

```text
reviewed_sha256_expected = true
```

### Exact source-owned scientific comparison

The observed existing F.2 and the reviewed F.2 were compared on the farm using
the pushed source-owned functions and exclusion set from:

```text
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py

scientific_projection(...)
first_mismatch(...)
F2_PROVENANCE
```

Supplied result:

```text
f2_scientific_payload_match = true
first_mismatch = null
```

This means the two F.2 bodies are scientifically identical under the same
closed exact projection already reviewed for the private F.6.3 equivalence
bridge. Their serialized/provenance identities differ.

This evidence authorizes only the exact predecessor-state migration below.

---

## 4. Objective

Permit the canonical-five owner to replace **exactly one newly observed,
scientifically equivalent, byte-pinned pre-existing F.2 state** with the already
reviewed F.2 candidate.

The observed predecessor hash:

```text
182433ccd2d13d3f7af80b963e69190ad8f2eb4e295c8216b36585e8ee7409b2
```

may be treated only as:

```text
recognized existing target state -> safe to replace with reviewed F.2 bytes
```

It must never be treated as:

```text
accepted scientific authority
reviewed runtime candidate
materialization input authority
fallback F.2
current F.2 after staging
```

After staging, the target must still have exact reviewed SHA-256:

```text
2fa715b2b5c2d5077e38416fb9eceb814e013e72f308c43526fddff36ddd962e
```

All other unknown pre-existing F.2 hashes must continue to fail closed.

---

## 5. Required source behavior

### 5.1 Explicit predecessor constant

In:

```text
testing/run_e8_4_fix5_canonical_five_plot_gate.py
```

define one explicit source-owned predecessor-state mapping/set for candidate
staging.

It must include only the exact F.2 candidate basename and exact observed hash:

```text
Q4p4W2p74_kaon_pion-background_hgcer-method-a-acceptance-representation-current-baseline-candidate.json
  ->
182433ccd2d13d3f7af80b963e69190ad8f2eb4e295c8216b36585e8ee7409b2
```

A clear name such as:

```text
RECOGNIZED_STAGING_PREDECESSORS
```

is preferred.

Do not put this hash into:

- `MATERIALIZATION_SHA256`;
- `CANDIDATES`;
- `HISTORICAL_INPUT_SHA256`;
- F.6.3 authority constants;
- scientific-equivalence authority records;
- accepted Left/lowe `CANDIDATES`;
- any production authority.

It is a staging-state exception only.

### 5.2 Existing F.3/F.4 behavior unchanged

Current recognized F.3/F.4 predecessor behavior through
`accepted.CANDIDATES.get(name)` must remain unchanged.

Do not rewrite historical Left/lowe owner constants just to generalize this
repair.

### 5.3 Fail-closed staging

For each candidate target, allowable `before` states remain:

```text
absent
already exact reviewed candidate
already-recognized historical predecessor from existing owner behavior
plus, for F.2 only, exact 182433... observed predecessor
```

Any other pre-existing hash must raise:

```text
unknown_candidate_target:<basename>
```

before any replacement.

Symlink/non-file/source-target-alias checks remain unchanged.

The owner must continue to:

1. hash-check the reviewed source materialization;
2. copy to a temporary file;
3. hash-check the temporary copy;
4. recheck that the destination has not changed since planning;
5. atomically replace;
6. require the final destination hash equals the reviewed candidate hash.

No direct overwrite, delete-first, wildcard, directory cleanup, or unconditional
replacement is allowed.

### 5.4 Observability

The existing installation record already stores:

```text
path
before_sha256
after_sha256
action
```

That is sufficient if preserved exactly.

For the observed predecessor case the deterministic record must show:

```text
before_sha256 =
  182433ccd2d13d3f7af80b963e69190ad8f2eb4e295c8216b36585e8ee7409b2

after_sha256 =
  2fa715b2b5c2d5077e38416fb9eceb814e013e72f308c43526fddff36ddd962e

action = replaced
```

Do not add language implying the predecessor was accepted authority.

---

## 6. Scientific source frozen

Do not modify any scientific or runtime-analysis source, including:

```text
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
src/binning/calculate_yield.py
src/cuts/pion_hgcer_method_a_acceptance_contract.py
src/cuts/pion_hgcer_method_a_acceptance_representation.py
src/cuts/pion_hgcer_method_a_acceptance_map.py
src/cuts/pion_hgcer_method_a_parent_preserving_correction.py
src/cuts/pion_hgcer_method_a_tphi_propagation.py
src/cuts/full_background_subtraction_plots.py
```

Also do not modify:

```text
run_Prod_Analysis.sh
set_SymLinks.sh
testing/pion_hgcer_validation_bundle_profile_e8_4_canonical_five.json
testing/collect_pion_hgcer_validation_bundle.py
farm_env/
background_samples/
```

This task changes no:

- F.1 population;
- F.2/F.3/F.4 science;
- scientific-equivalence projection;
- thresholds/tolerances;
- parent normalization;
- Method-A factors;
- Method-B behavior;
- yields;
- cuts;
- SIMC;
- production subtraction;
- procedure-page science.

---

## 7. Allowed source/test paths

Only:

```text
testing/run_e8_4_fix5_canonical_five_plot_gate.py
testing/test_run_e8_4_fix5_canonical_five_plot_gate.py
```

may change for executable behavior.

If another executable path is required, return `BLOCKED`.

---

## 8. Required tests

Add focused deterministic tests proving:

### Positive predecessor migration

With the F.2 destination pre-populated with bytes whose SHA-256 is exactly:

```text
182433ccd2d13d3f7af80b963e69190ad8f2eb4e295c8216b36585e8ee7409b2
```

and the source materialization supplying the reviewed F.2:

- staging succeeds;
- action is `replaced`;
- recorded `before_sha256` is exactly `182433...`;
- recorded `after_sha256` is exactly `2fa715...`;
- final destination bytes/hash are the reviewed F.2;
- F.3/F.4 staging behavior remains normal.

Use an exact fixture for the observed predecessor bytes if those bytes are
already available in an existing deterministic fixture. If they are not
available locally, do **not** invent bytes that hash to `182433...`.
Instead test the policy boundary by narrowly patching the hash reader for the
pre-existing F.2 planning/recheck only, while leaving real source-copy/final-hash
behavior exercised. The test must make that mocking explicit.

### Existing valid states

Retain/prove:

- absent F.2 -> installed;
- already-reviewed F.2 -> no-op;
- existing recognized F.3/F.4 predecessor -> still replaceable.

### Negative fail-closed cases

Require failure for:

- any F.2 hash other than absent/current/exact `182433...`;
- wrong candidate source hash;
- target changes between planning and atomic replacement;
- target symlink;
- target non-file;
- source/target alias.

Do not broaden to a scientific-comparison-at-staging fallback.

The staging decision must be exact hash membership, not “scientific payload
looks equal” at runtime. The exact scientific comparison has already been
performed as diagnostic evidence and is the reason this one hash may be
recognized.

---

## 9. Regression tests

Run at minimum:

```text
testing/test_run_e8_4_fix5_canonical_five_plot_gate.py
testing/test_f6_3_parallel_full_procedure_method_a.py
testing/test_pion_hgcer_validation_bundle_profile_e8_4_canonical_five.py
testing/test_compare_method_a_current_baseline_authority.py
testing/test_collect_pion_hgcer_validation_bundle.py
```

Also run the complete owner-focused test group used in the prior repair if
readily available.

No ROOT/PyROOT or farm execution.

---

## 10. Runtime evidence record

Add:

```text
docs/memory/evidence/e8-4-canonical-five-f2-staging-predecessor-blocker-2026-10-07.md
```

Record only supplied/accepted diagnostic facts:

- failed attempt stem;
- source `946b5512...`;
- gate stopped at `candidate_staging`;
- analysis never started;
- exact failure reason;
- observed F.2 raw/body/artifact identities;
- empty/null provenance block;
- reviewed F.2 raw/fingerprint identities;
- farm read-only source-owned scientific comparison:
  `f2_scientific_payload_match=true`, `first_mismatch=null`;
- interpretation: exact `182433...` bytes are accepted only as a replaceable
  staging predecessor state, not scientific authority;
- canonical-five runtime remains `BLOCKED`.

Do not claim farm closure from this record.

---

## 11. Memory/status paths allowed

After successful local implementation/checks, only these memory paths may
change:

```text
docs/memory/evidence/e8-4-canonical-five-f2-staging-predecessor-blocker-2026-10-07.md
docs/memory/phases/e8-4-canonical-five-f2-staging-predecessor-repair-task-contract.md
docs/memory/phases/e8-4-canonical-five-f2-staging-predecessor-repair.md
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/roadmap/STATUS.md
docs/memory/manifest.json
```

Prefer not to edit MEMORY unless a genuinely reusable rule is needed. The
existing durable exact-lineage/scientific-equivalence rule already covers the
main lesson.

### CURRENT

Record the failed staging attempt as historical diagnostic evidence.

After deterministic local repair checks, set the narrow owner repair to:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

Keep the earlier F.2/F.3/F.4 equivalence repair at:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

until a fresh full canonical-five run validates the combined source.

Keep:

```text
canonical-five full runtime gate: BLOCKED
final E.8: BLOCKED
F.6.4: BLOCKED
absolute-SIMC interpretation: BLOCKED
Method A: detached/non-production
Method B: diagnostic/cross-check only, numerically excluded
```

Push-stable sole NEXT:

```text
NEXT — after ChatGPT actual-diff review, user commit/push, and pushed-state
synchronization, perform a fresh farm-readiness audit for one isolated
Q4p4W2p74 canonical-five owner retry. Only Farm readiness: PASS may authorize
the retry.
```

Do not authorize the farm from the local implementation.

---

## 12. Local health and integrity checks

Discover `<PYTHON>` according to repository convention.

After all changes:

```bash
<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/check_memory_health.py --root .
<PYTHON> -B tools/memory_bootstrap.py --root . --json
git diff --check
```

Report:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
warning classification: <none or each warning: blocking | nonblocking>
manifest check: PASS | FAIL
bootstrap: PASS | FAIL
git diff --check: PASS | FAIL
```

Local checks establish only `SOURCE VERIFIED` behavior.

---

## 13. Review bundle

Create repository-root temporary:

```text
kaonlt_review.diff
```

It must contain:

- complete diffs for every modified tracked file;
- complete additions for every new file using
  `git diff --no-index /dev/null ...`;
- no unrelated paths.

Do not stage merely for review.

Stop before commit/push.

---

## 14. Acceptance criteria

Return:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

only if:

- exact start identity passed;
- unrelated state was preserved;
- executable changes are limited to owner + owner tests;
- exact `182433...` F.2 may be replaced by exact reviewed `2fa715...` F.2;
- arbitrary unknown F.2 hashes still fail closed;
- existing F.3/F.4 predecessor policy is unchanged;
- source materialization hash checking remains unchanged;
- atomic staging/recheck/final-hash requirements remain unchanged;
- no scientific/runtime-analysis source changed;
- no equivalence projection changed;
- focused and regression tests pass;
- memory/manifest/bootstrap/diff checks pass;
- complete review bundle exists;
- no farm/staging/commit/push/ref update occurred.

---

## 15. Hard stops

Return `BLOCKED` if:

- branch/HEAD/origin differ;
- unrelated local state cannot be preserved;
- implementation requires modifying scientific/runtime-analysis source;
- implementation requires accepting a wildcard/range/class of F.2 predecessor
  hashes instead of exact `182433...`;
- implementation would classify `182433...` as accepted scientific authority;
- implementation would allow scientific-comparison-at-staging to replace
  arbitrary unknown targets;
- materialization source hash verification must be weakened;
- atomic replacement/recheck/final-hash verification must be weakened;
- another executable path is required;
- any test shows an arbitrary unknown F.2 can be overwritten;
- memory/manifest hard checks fail.

Do not run the farm from this task.
