# KaonLT — canonical-five fresh-F1 candidate-lineage refresh after failed E.8.4 gate

## 1. Objective

Repair the exact blocker exposed by the first isolated canonical-five farm run.

The full Q4p4W2p74 analysis completed successfully in the isolated runtime.
The canonical-five owner then correctly failed because every setting rendered
E.8.4 unavailable with:

f6_3_f4_shared_reproduction_failed:f3_fingerprint_input_content_mismatch

Read-only farm diagnosis established that all five fresh F.1 artifacts differ
from the previously pinned F.1 lineage. A detached comparison using those fresh
five F.1 inputs then established:

- F.2 scientific payload match = true
- F.3 scientific payload match = true
- F.4 scientific payload match = false
- first changed stage = F4

A detached materialization of that fresh-F1 lineage subsequently passed.

This task must:
1. refresh only the private/non-production F.6.3 current-baseline candidate
   authority pins;
2. give the existing canonical-five owner explicit ownership of verifying and
   staging that reviewed candidate materialization before full analysis;
3. preserve all established physics and historical accepted authorities.

No Method-A redesign, F.2/F.3 scientific redesign, production promotion, or
baseline-analysis change is allowed.

## 2. Required starting identity

Branch:

test

Required starting HEAD and origin/test:

463d2657f696ecee33113edc3393ac51083a8944

Before editing run:

git branch --show-current
git rev-parse HEAD
git rev-parse origin/test
git status --short --untracked-files=all

Hard stop on branch/HEAD/origin mismatch or unrelated tracked changes.

Do not reset, clean, stash, commit, push, update refs, or run the farm.

## 3. Startup reads

Read in exact repository-required order:

1. AGENTS.md
2. docs/memory/CURRENT.md
3. docs/memory/MEMORY.md
4. docs/memory/handoffs/CURRENT_HANDOFF.md
5. docs/memory/USER.md

Then read:

- docs/memory/MAINTENANCE.md
- docs/memory/CODEX.md
- docs/memory/COMMUNICATION.md
- docs/memory/LEARNINGS.md
- docs/memory/TOOLS.md

Task-relevant source/tests:

- src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
- src/cuts/pion_hgcer_method_a_parent_preserving_correction.py
- src/cuts/pion_hgcer_method_a_tphi_propagation.py
- testing/test_f6_3_parallel_full_procedure_method_a.py
- testing/run_e8_4_fix5_canonical_five_plot_gate.py
- testing/test_run_e8_4_fix5_canonical_five_plot_gate.py
- testing/run_e8_4_fix5_left_lowe_plot_gate.py
- testing/pion_hgcer_validation_bundle_profile_e8_4_canonical_five.json
- testing/collect_pion_hgcer_validation_bundle.py

Do not reopen closed F stages or redesign the analysis.

## 4. Accepted fresh-F1 evidence for this task

Comparison input SHA-256:

4c3fdab5d05965a1b3f1dc835c9d8c529896e6eaf2da8eb5c8ac8fc26e95257a

Materialization manifest SHA-256:

e7c37f55b24739e8be0fb778a6c3491e17df746263e6d9d71a4f31ec65cd7908

Materialization source head:

463d2657f696ecee33113edc3393ac51083a8944

Manifest gates:

complete = true
errors = []
non_authoritative = true
accepted_authority_mutated = false
production_objects_mutated = false
production_application_performed = false
method_a_promoted = false

Scientific gate:

f2_scientific_payload_match = true
f3_scientific_payload_match = true
f4_scientific_payload_match = false
first_changed_stage = F4

### Fresh F.1 raw SHA-256

Left-lowe:
eb6f659da0511f0f6ea867420fda0ee6feb77e6698509ccd31f80a92cc541c07

Left-highe:
544ea08f71e74b6b59bc33d05458d4acc051e96fb0327f54ec092a331245c01e

Center-lowe:
2b193dec46b10aebc74dfb51634d19a7fef29362fd0a252d940f2894d18a0a16

Center-highe:
c857911396bed03f9e418ca45a609509ebf1eadf8084d4c0ff2f9cc3a9980941

Right-highe:
e03245494d54a6f812c426d189b220417fdc0888ba1f58a961fbe250ff097652

### Fresh candidate F.2

basename:
Q4p4W2p74_kaon_pion-background_hgcer-method-a-acceptance-representation-current-baseline-candidate.json

raw SHA-256:
2fa715b2b5c2d5077e38416fb9eceb814e013e72f308c43526fddff36ddd962e

representation fingerprint:
e9df6e09c31dda4d7d724792ce3c40125e27530c4c76e0dbe7f5b9861d0e6216

artifact fingerprint:
87d9a8780cb2de9f2830afd9e388064151f796654ddb5c4a8e1938443806a0c6

### Fresh candidate F.3

basename:
Q4p4W2p74_kaon_pion-background_hgcer-method-a-acceptance-map-current-baseline-candidate.json

raw SHA-256:
c5b86452b790ecbaf5b8f0df05da67efa2fa92aab157b12b153ed7e491a38228

map fingerprint:
6042e9485827d9884d2a841d5c6cb5e23401e1026c8e2d4bd2022ede15843548

algorithm fingerprint:
ba29630b2f40a87cbadbe751504ce48f23e2b17a08378a2c8a219f93131cb912

artifact fingerprint:
e1cf0dd14644034f5d21df27d29f48355cde7c8e85c6123b709a6975c152a302

candidate-construction farm-source sentinel:
0000000000000000000000000000000000000000

### Fresh candidate F.4

basename:
Q4p4W2p74_kaon_pion-background_hgcer-method-a-parent-preserving-correction-current-baseline-candidate.json

raw SHA-256:
79e7ceda7221cbeeead4ed5bc306b0e0e670741a27beaa980e22349c555e96d7

correction fingerprint:
71527eecd5828766ebcb3b9ef5e4e1931eb242fabfc1059f04e4eec104f6cc98

artifact fingerprint:
0b3899e5a3f24a34e04aa550ef162414c1d95af9a786e748612cbcb018dbcad4

It inherits the exact fresh F.3 identities above.

## 5. Scientific boundaries

Do not modify:

- random/dummy subtraction;
- slow-proton subtraction;
- baseline pion subtraction;
- pion templates/fits/priors;
- Method-A mathematics;
- F.2 representation mathematics;
- F.3 hgcer3 mathematics;
- F.4 parent-preserving mathematics;
- F.5 propagation mathematics;
- F.6.3 weighting formula;
- E.8.4 plotting mathematics;
- Method B;
- SIMC;
- cuts/windows/binning;
- efficiencies/acceptance/yields;
- L/T separation/cross sections;
- no_empirical_residual profile.

Method A remains detached/non-production.
Method B remains diagnostic/cross-check only and numerically excluded.

Do not modify the historical general authorities:

src/cuts/pion_hgcer_method_a_parent_preserving_correction.py
    ACCEPTED_F3_RUNTIME_AUTHORITY_BY_KINEMATIC

src/cuts/pion_hgcer_method_a_tphi_propagation.py
    ACCEPTED_F4_RUNTIME_AUTHORITY_BY_KINEMATIC

Do not modify the historical Left/lowe owner or its CANDIDATES.

## 6. Refresh F.6.3 private candidate lineage

Modify only the source-owned candidate records in:

src/cuts/pion_hgcer_method_a_parallel_full_procedure.py

Update:

F6_3_CANDIDATE_F1_SOURCE_FILE_SHA256

to the five fresh F.1 hashes in Section 4.

Update:

F6_3_CANDIDATE_F3_RECONSTRUCTION_AUTHORITY_BY_KINEMATIC

to the fresh F.3 raw SHA, map fingerprint, algorithm fingerprint and artifact
fingerprint in Section 4, retaining farm_source_head = forty zeroes.

Update:

F6_3_CANDIDATE_F4_VALIDATION_AUTHORITY_BY_KINEMATIC

to the fresh F.4 raw SHA, correction fingerprint and artifact fingerprint,
inheriting the refreshed F.3 identities and fresh F.1 hashes.

Set:

F6_3_CANDIDATE_VALIDATION_SOURCE_HEAD =
463d2657f696ecee33113edc3393ac51083a8944

This identifies the source under which the reviewed fresh-F1 materialization was
constructed. It is not production authority.

Preserve:

- exact candidate basenames;
- exact persisted-versus-recomputed F.4 equality gate;
- exact live-cache parity;
- unsupported-kinematic fail closed;
- no fallback to historical F.3/F.4.

## 7. Canonical-five owner must own candidate staging

Current canonical-five source incorrectly inherits:

CANDIDATES = accepted.CANDIDATES

from the historical Left/lowe owner and requires those files to already be
installed in the analysis artifact root.

Replace that inheritance only in:

testing/run_e8_4_fix5_canonical_five_plot_gate.py

with canonical-five-owned fresh candidate identities from Section 4.

Keep the historical accepted Left/lowe owner unchanged.

### 7.1 New required CLI input

Add a required argument:

--candidate-materialization-dir

It identifies the already-produced detached materialization directory.

Do not store a machine-specific concrete path in repository source or memory.

### 7.2 Verify the supplied materialization before installation

The owner must verify, fail closed, before full analysis:

- directory exists;
- comparison-input file exists and SHA is exactly
  4c3fdab5d05965a1b3f1dc835c9d8c529896e6eaf2da8eb5c8ac8fc26e95257a;
- materialization manifest exists and SHA is exactly
  e7c37f55b24739e8be0fb778a6c3491e17df746263e6d9d71a4f31ec65cd7908;
- F.2 candidate exists and SHA is exactly
  2fa715b2b5c2d5077e38416fb9eceb814e013e72f308c43526fddff36ddd962e;
- F.3 candidate exists and SHA is exactly
  c5b86452b790ecbaf5b8f0df05da67efa2fa92aab157b12b153ed7e491a38228;
- F.4 candidate exists and SHA is exactly
  79e7ceda7221cbeeead4ed5bc306b0e0e670741a27beaa980e22349c555e96d7.

Strict-parse the JSON files.

Require the manifest to reproduce exactly:

- schema_version = method_a_current_baseline_authority_materialization/v1;
- source_head = 463d2657f696ecee33113edc3393ac51083a8944;
- complete = true;
- errors = [];
- non_authoritative = true;
- accepted_authority_mutated = false;
- production_objects_mutated = false;
- production_application_performed = false;
- method_a_promoted = false;
- exact historical accepted F.2/F.3/F.4 input SHA values;
- exact comparison-input SHA;
- exact fresh five-setting F.1 raw SHA inventory;
- exact candidate F.2/F.3/F.4 hashes/fingerprints;
- exact F.3 diagnostic zero-head reconstruction record;
- F.2 match=true;
- F.3 match=true;
- F.4 match=false;
- first_changed_stage=F4.

Do not infer or discover a newer candidate.

### 7.3 Install only candidate F.3/F.4

After successful materialization verification, stage only the fresh F.3 and F.4
candidate files into the canonical analysis artifact root under their existing
candidate basenames.

Do not install/replace historical accepted F.2/F.3/F.4 authorities.

For each existing candidate target allow only:

- absent;
- exact previous historical-current-candidate hash from
  accepted.CANDIDATES;
- exact new fresh-candidate hash.

Any other existing hash fails closed.

Copy through owner-created temporary files in the destination directory,
fsync/close them, verify byte SHA, then os.replace to the candidate basename.

If the destination already has the exact new hash, record a no-op.

Never alter the source materialization directory.

No automatic rollback is required after later analysis failure: the installed
fresh candidate files are detached/non-authoritative reviewed artifacts, not
production state. Record exact before/after hashes and installation action in
the gate status and run summary.

Candidate installation must complete before launching run_Prod_Analysis.sh.

A failed candidate verification/staging gate must not start analysis.

## 8. Preserve the existing isolated runtime owner

Do not weaken or redesign:

- ordinary-checkout snapshot/preservation;
- copied ltsep overlay;
- three LTANAPATH caller probes;
- derived OUTPATH validation;
- external SIMC no-mutation preflight;
- controlled child environment;
- detached analysis worktree;
- unchanged command:
  ./run_Prod_Analysis.sh 4p4 2p74
- low/high completion markers;
- exact canonical-five inventory;
- five-setting freshness/page checks;
- full-analysis/ledger checks;
- generic collector;
- ZIP verification;
- final preservation checks;
- companion delivery.

The refreshed candidate files remain special non-fresh inputs exactly as the
existing CANDIDATES handling intends.

## 9. Required tests

Update:

testing/test_f6_3_parallel_full_procedure_method_a.py

to assert all fresh F.1/F.3/F.4 source-owned identities exactly.

Preserve its synthetic fail-closed reproduction, historical-authority
immutability, parent preservation and live-cache-parity tests.

Update:

testing/test_run_e8_4_fix5_canonical_five_plot_gate.py

with deterministic tests proving:

1. canonical-five CANDIDATES equals exactly the fresh F.3/F.4 hashes;
2. accepted.CANDIDATES remains the historical Left/lowe mapping and is not
   modified;
3. valid materialization verification passes;
4. bad comparison hash fails;
5. bad manifest hash fails;
6. wrong source_head fails;
7. wrong scientific gate fails;
8. wrong F.1 inventory fails;
9. wrong F.3/F.4 identity/fingerprint fails;
10. unknown pre-existing candidate target fails closed;
11. old known candidate target is replaced with fresh bytes;
12. absent target is installed;
13. already-fresh target is a no-op;
14. source materialization bytes remain unchanged;
15. staging failure prevents analysis;
16. successful owner flow records materialization/staging provenance and still
    performs the existing isolated analysis/verify/collect/ZIP/delivery order.

Preserve all existing canonical-five regressions.

No ROOT/farm claim follows from local tests.

## 10. Allowed paths

Source:

src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
testing/run_e8_4_fix5_canonical_five_plot_gate.py

Focused tests:

testing/test_f6_3_parallel_full_procedure_method_a.py
testing/test_run_e8_4_fix5_canonical_five_plot_gate.py

Contract:

docs/memory/phases/e8-4-canonical-five-fresh-f1-candidate-lineage-refresh-task-contract.md

Runtime-evidence checkpoint:

docs/memory/evidence/e8-4-canonical-five-fresh-f1-lineage-materialization-2026-10-05.md

Minimal active memory:

docs/memory/CURRENT.md
docs/memory/LEARNINGS.md
docs/memory/manifest.json

If another scientific/runtime source path is required, STOP and report why.

## 11. Frozen paths

Do not modify:

run_Prod_Analysis.sh
set_SymLinks.sh
farm_env/
background_samples/
src/cuts/pion_hgcer_method_a_parent_preserving_correction.py
src/cuts/pion_hgcer_method_a_tphi_propagation.py
src/cuts/pion_hgcer_method_a_acceptance_map.py
src/cuts/pion_hgcer_method_a_acceptance_representation.py
testing/run_e8_4_fix5_left_lowe_plot_gate.py
testing/pion_hgcer_validation_bundle_profile_e8_4_fix5.json
testing/pion_hgcer_validation_bundle_profile_e8_4_canonical_five.json
testing/collect_pion_hgcer_validation_bundle.py

No production/scientific mathematics may change.

## 12. Evidence checkpoint content

The new evidence record must state narrowly:

- the first isolated canonical-five run completed analysis but failed the
  downstream E.8.4 page gate;
- all five page manifests rendered only e8_4.unavailable;
- all five reasons were
  f6_3_f4_shared_reproduction_failed:f3_fingerprint_input_content_mismatch;
- fresh F.1 raw/stable/training/application identities differed from the
  previously pinned candidate lineage;
- detached comparison proved F.2 and F.3 scientific equality and F.4 as the
  first changed stage;
- detached materialization completed with no errors and no accepted-authority
  or production mutation;
- record the exact comparison, manifest, candidate F.2/F.3/F.4 identities from
  Section 4;
- canonical-five runtime is NOT validated by this evidence;
- Method A remains detached/non-production.

Use no personal filesystem paths.

## 13. CURRENT

Replace stale pre-farm wording.

During implementation, CURRENT must state that the first canonical-five farm
gate is BLOCKED by stale candidate-lineage pins exposed by fresh F.1.

After implementation/local checks, CURRENT may state:

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING

for the narrow fresh-F1 lineage refresh, conditional on:

ChatGPT actual-diff review
-> user commit/push
-> ChatGPT pushed-state synchronization
-> one new canonical-five farm run using the reviewed materialization directory.

Do not claim canonical-five runtime validation.

Keep CURRENT below its soft size limit.

## 14. Deterministic checks

At minimum run:

python -m py_compile \
  src/cuts/pion_hgcer_method_a_parallel_full_procedure.py \
  testing/run_e8_4_fix5_canonical_five_plot_gate.py \
  testing/test_f6_3_parallel_full_procedure_method_a.py \
  testing/test_run_e8_4_fix5_canonical_five_plot_gate.py

python -m unittest \
  testing.test_f6_3_parallel_full_procedure_method_a \
  testing.test_run_e8_4_fix5_canonical_five_plot_gate \
  testing.test_full_background_subtraction_plots \
  testing.test_run_prod_analysis_debug_left_low \
  testing.test_collect_pion_hgcer_validation_bundle \
  testing.test_pion_hgcer_validation_bundle_profile_e8_4_canonical_five

Then:

git diff --check

Regenerate/check docs/memory/manifest.json and run ordinary memory health and
bootstrap checks required by repository memory policy.

Report exact test totals and memory byte counts/warnings.

## 15. Review handoff

Generate one complete:

kaonlt_review.diff

containing every allowed tracked diff and every intended new file in full.

Stop before staging, commit, push, or farm execution.

Do not run the farm.

## 16. Mandatory user workflow correction: real lineage preflight

After reviewed materialization verification/staging and before `run_analysis`,
the canonical-five owner must run a no-ROOT/no-tree/no-full-analysis F.6.3
lineage preflight for exactly Left-lowe, Left-highe, Center-lowe, Center-highe,
and Right-highe. Load the current five F.1 and staged F.3/F.4 artifacts through
`accepted_f6_3_artifact_paths` / `load_accepted_f6_3_authority`. Require observed
F.1 hashes to equal the refreshed source-owned F.1 record and staged F.3/F.4
hashes to equal the canonical-five owner identities. Call the real unchanged
`reconstruct_transient_factor_map` for every setting, exercising F.5 authority,
shared F.4 reconstruction and exact persisted/recomputed equality. Require
nonempty finite positive factors, complete setting provenance and all five
settings passing before analysis. Record exact observed identities, per-setting
factor counts and F.4 reproduction success in status and run summary. Tests
must prove any F.1, F.3, F.4, reconstruction, authority or setting failure blocks
analysis. Passing local tests alone establishes neither farm readiness nor
runtime acceptance. Preserve every other contract boundary and stop before
staging, commit, push or farm execution.

### Separate real-farm gate (mandatory actual-diff orchestration repair)

The real-farm lineage preflight is a separate gate. Add tracked CLI mode
`--lineage-preflight-only` to the existing canonical-five owner. It follows
the same ordinary source preflight, reviewed-materialization verification,
candidate F.3/F.4 staging, detached worktree, copied-ltsep overlay/import/path
probes, external-link no-mutation gate and real five-setting reconstruction.
It then cleans up the worktree and proves final ordinary-checkout and installed
ltsep preservation before publishing success. Return the atomic owner status
receipt path on stdout; `--output` supplies the fresh attempt stem, without
creating a ZIP. The receipt records source commit, verification, candidate
before/after/action records, exact observed F.1/F.3/F.4 hashes, five setting IDs,
factor counts, F.4 reproduction success, cleanup and final preservation,
analysis_started=false, and an explicit preflight-only/non-runtime-validation
role. No analysis/launcher, completion-marker or PDF/artifact verification,
collector bundle, ZIP verification or normal companion delivery runs in this
mode. Failures remain failed before analysis. Normal mode retains its existing
lineage -> analysis -> artifact verification -> collection order.

After source review, user commit/push and pushed-state synchronization, the
sole next farm operation is one fast lineage-preflight-only gate, followed by
evidence review. The hours-long full canonical-five analysis is authorized
only after that separate real-farm gate is reviewed and accepted, with later
explicit authorization. This supersedes the full-run next-step wording in
Section 13; it does not reopen science or alter any candidate identity, profile,
historical authority or materialization evidence. Deterministic tests must
exercise five real reconstructions, forbidden full-run operations, final
cleanup/preservation, fail-closed identities/reconstruction/authority, and
normal-mode ordering. Rerun Section 14 and existing regressions, memory checks
and complete nine-file review bundle. Stop before Git staging/commit/push/farm.
