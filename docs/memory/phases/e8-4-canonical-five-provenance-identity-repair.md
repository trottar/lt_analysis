# E.8.4 / F.6.3 canonical-five provenance/identity repair: blocked source audit

## Authority and outcome

2026-10-07: `BLOCKED` at the mandatory source audit under the
[repair contract](e8-4-canonical-five-provenance-identity-repair-task-contract.md).
Start identity passed: branch `test`, HEAD/local `origin/test`
`ab3f29c9268a80cd872903003da0422195d118e8`. The user explicitly reopened
the formerly deferred repair. No scientific/runtime/test source was changed.

Contract section 8 freezes the tracked owner and requires stopping if it must
change to expose or enforce the equivalence gate; section 15 repeats that
hard stop. The reviewed candidate F.2 scientific payload has no runtime
handoff. Implementing only an F.3/F.4 comparison would violate section 5.3.
The original contract's owner-change hard stop was honored before implementation.
This is a source/interface blocker, not a request for another farm artifact.

## Supplied ChatGPT actual-diff decision — 2026-10-07

The user supplied that ChatGPT actual-diff review accepted the source-audit
conclusion and approved a narrowly broadened follow-up permitting the
canonical-five owner to hand off the already reviewed/hash-pinned candidate
F.2 alongside F.3/F.4. The reviewed bundle was `kaonlt_review.diff`,
40,439 bytes, SHA-256
`e07119b9b6c05ade2b151008c14999aa80e6a829ce201021d8eeb2e19030293a`.
This records supplied review acceptance, not a new independent review or
runtime result. The audit conclusion is accepted; repair remains `BLOCKED`
and no implementation is authorized by this documentation amendment.

After user commit/push and ChatGPT pushed-state synchronization of this
blocked-audit checkpoint, the substantive next gate is one new narrow
source-changing owner/F.2-handoff repair contract for that approved follow-up.
The original contract stays unchanged as the historical audit authority;
the new contract must explicitly own the narrow owner scope. Approval does
not permit scientific changes, weaker default lineage checks, a hash-only
equivalence proof, provenance relabeling, runtime acceptance or farm execution.
All existing `BLOCKED`/closed statuses and scientific boundaries are retained.
CURRENT alone owns the sole push-stable NEXT.

## Exact mismatch and consumer path — SOURCE VERIFIED

- `src/cuts/rand_sub.py:5538-5568` constructs/writes fresh F.1 contracts.
- `src/binning/calculate_yield.py:3920-3988` discovers and loads accepted
  F.6.3 authority, calls `reconstruct_transient_factor_map`, and converts its
  exception into an unavailable parallel source for downstream E.8.4.
- `src/cuts/pion_hgcer_method_a_parallel_full_procedure.py:111-145`
  discovers/loads five F.1 files and staged F.3/F.4 only. Its reconstruction
  passes staged F.3 with current F.1 into the shared F.4 builder and wraps
  failure as `f6_3_f4_shared_reproduction_failed:<reason>`.
- `src/cuts/pion_hgcer_method_a_parent_preserving_correction.py:251-317`
  retains strict persisted-F.3 validation. Lines 278-280 compare stored
  `fingerprint_inputs.input_content` with current setting IDs, stable F.1
  fingerprints and raw F.1 SHA-256 values, rejecting with
  `f3_fingerprint_input_content_mismatch`. The next check separately compares
  F.3 input setting/raw/stable/contract fingerprints to current F.1 authority.

The supplied [second-run diagnosis](../evidence/e8-4-canonical-five-second-full-run-provenance-blocker-2026-10-05.md)
is RUNTIME VERIFIED diagnostic evidence: all five regenerated F.1 raw/stable
identities differed, while exact scientific comparison found
`first_changed_stage = none` through F.2/F.3/F.4. No new runtime inspection
or canonical-five acceptance occurred in this task.

## Scientific equality versus provenance inequality

SOURCE VERIFIED: the comparator rebuilds current F.2/F.3/F.4 through their
public builders, then compares explicit scientific projections exactly.
F.2's stable F.1 identity includes raw records and authority fields
(`pion_hgcer_method_a_acceptance_representation.py:369-401`), whereas the
scientific representation is the frozen derived result. Thus identical
derived science need not imply identical serialized/source authority.
The precise low-level cause of the supplied run-to-run identity changes
remains NOT VERIFIED; no provenance field is relabeled as current authority.

The source audit covered every top-level exclusion in
`testing/compare_method_a_current_baseline_authority.py`:

- F.2/F.3/F.4 `input_fingerprints` describe input setting and F.1 identities;
  `fingerprint_inputs` combines those identities with copies of scientific
  fields retained in the compared body; `fingerprint` hashes that combination.
- F.3 `f2_representation_fingerprint` and `f2_source_file_sha256` identify
  the F.2 parent. F.3 retains its F.2 algorithm fingerprint and all model,
  scaler, coefficient, support and basis content in the projection.
- F.4 `f3_source_file_sha256`, `f3_map_fingerprint` and
  `f3_artifact_fingerprint` identify the F.3 parent. `f3_runtime_authority`
  is its source/lineage binding; the algorithm identity is also retained in
  the F.4 scientific body. All parent science remains compared.

The builders' fingerprint construction establishes these roles
(F.2:925-1095, F.3:568-626, F.4:431-446). The comparator removes only named
top-level fields, retains unknown fields, and reports exact recursive first
mismatch without tolerances. Its projections compare bodies, so any runtime
bridge must additionally enforce full wrapper contracts and flags rather
than assume body equality validates the wrappers. No field exclusion was
added or promoted into runtime by this audit.

## Narrow boundary and concrete owner hard stop

The narrow repair boundary is F.6.3 candidate/current-lineage reconstruction,
not the default F.4 validator or frozen scientific builders. A permissible
future gate would reconstruct current F.2/F.3/F.4 from actual five F.1 inputs,
enforce contracts/flags and the 15-parent inventory, compare all three full
scientific projections exactly, and only then consume current-lineage factors.
It must retain candidate authority and current provenance separately. Exact
comparison at each stage, retained unknown fields, wrapper validation and
negative tests are necessary to reject genuinely changed science.

SOURCE VERIFIED: that design cannot currently receive its reviewed F.2 input:

- `testing/run_e8_4_fix5_canonical_five_plot_gate.py:47-57,97-202`
  verifies F.2/F.3/F.4 in the external materialization directory, but
  `CANDIDATES` and `stage_candidates` deliberately stage only F.3/F.4.
  Verification returns identity records, not the F.2 scientific payload.
- Owner `analysis_environment:373-380`, preflight and `execute_gate`
  provide no materialization-directory/F.2 handoff to the analysis consumer.
- F.6.3 paths/load and the `calculate_yield.py` caller accept only F.1/F.3/F.4.
- F.3 uses a supported subset of F.2, not its full scientific payload
  (`pion_hgcer_method_a_acceptance_map.py:458-547`). Its F.2 hashes cannot
  reconstruct all F.2 groups/candidates/summaries/recommendation.
- Candidate materialization constructs F.2/F.3/F.4 with `input_paths={}`
  (`testing/materialize_method_a_current_baseline_authority.py:185,194,208`);
  the staged F.3 wrapper supplies no external F.2 path.
- `OwnerTests.test_reviewed_materialization_verifies_and_stages_only_f3_f4`
  explicitly enforces this staging boundary in the current owner regression.

INFERENCE: supplying the reviewed full F.2 payload through the owner requires
changing that frozen owner boundary. Filesystem guessing, a hash-only F.2
proof, an F.3/F.4-only gate or a refreshed pin would not satisfy this contract.
No such workaround was implemented. The audit also inspected all ten mandated
source/test files, including the exact-lineage F.6.3 and strict F.4 tests.

## Validation and retained boundaries

Only this audit record, CURRENT, roadmap and regenerated manifest are changed;
the supplied task contract is included unchanged as a complete new addition
in the review bundle. Unrelated tracked bytes, index and refs are preserved.
The seven contract regression modules ran 119 tests: six assertion failures
and two errors, all in owner fixture cases stopped by
`artifact_stale:Center_kaon_rand_sub_Q4p4W2p74_highe_full-background-subtraction-manifest.json`.
A local `/tmp` probe wrote a file after sampling `time.time_ns()` but observed
its filesystem mtime 14,346,210 ns earlier than that sample. SOURCE VERIFIED:
owner freshness requires `st_mtime_ns >= started_ns`; INFERENCE: local clock
discrepancy explains these fixture failures. The failures remain reported;
no owner/test changes or freshness bypass were made.
The other six contract modules plus the exact owner F.3/F.4 staging test
passed separately: 94 tests, zero failures/errors/skips. At the original audit checkpoint, manifest check,
ordinary memory health, bootstrap and `git diff --check` passed with no
health warnings. CURRENT/MEMORY/CURRENT_HANDOFF sizes are
8174/20623/323 bytes. These checks validate unchanged source/local bookkeeping,
not an implemented bridge. No new positive/negative bridge tests exist because
implementation stopped at the contract gate.

No ROOT/PyROOT, real analysis, farm, staging, commit, push or ref update.
Default strict F.3/F.4 validation and all scientific algorithms remain intact.
E.8.2/E.8.3 and existing F.6.3/E.8.4 Left/lowe closures retain their scopes;
the separate canonical-five preflight retains only its accepted narrow scope.
Full canonical-five runtime/PDF closure and final E.8/F.6.4 remain `BLOCKED`.
Absolute-SIMC units remain separately `BLOCKED`. Method A stays detached and
non-production; Method B stays diagnostic/cross-check only and numerically
excluded. CURRENT alone owns the exact next action.
