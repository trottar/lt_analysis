# F.6.3 — detached current-baseline candidate-lineage adoption

**Status:** `ACTIVE` — local implementation and deterministic checks complete;
independent ChatGPT actual-diff review pending. After that review passes, this
adoption is `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`. It does not close
F.6.3 or E.8.4 and does not promote Method A.

Implemented at exact `test` starting HEAD
`b349967c0d4210a78b144ce6134d3c1f15970245` under the
[adoption contract](f6-3-current-baseline-candidate-lineage-adoption-task-contract.md).

## Accepted upstream evidence

The supplied contract records independent ChatGPT acceptance of
`KaonLT_F4_Refresh2_materialization_Q4p4W2p74_20261001-120157.zip`, SHA-256
`4cbdbdd2c0403b961614a8ccbb91104f13fef59194cd48e5e3570aded0954898`,
farm HEAD `b349967c0d4210a78b144ce6134d3c1f15970245`, materializer source
`141a3d04f9e5d07be21dba14e0e63212c3990bf1`. The complete/error-free detached
gate is `CLOSED / RUNTIME VALIDATED`: F.2/F.3 scientific equality, F.4 first
changed stage, no accepted-authority mutation or production application.
Codex did not rerun the farm or independently inspect that ZIP in this task.

Candidate F.3 raw SHA-256:
`eeb480d9b7ddffab1e97c7f9f27c3099d0780c2bbb4d42a0ee68c077a5752b3d`;
map `3a9787fc58d26cc0816012bd1b637ad0c8f201b448d54cb1625a841131154728`;
algorithm `ba29630b2f40a87cbadbe751504ce48f23e2b17a08378a2c8a219f93131cb912`;
artifact `8d2d12068dfa98922d01dfafedbaa4b994bce0bf1e39ce23b1333872313ef121`.
Candidate F.4 raw SHA-256:
`1d545924eba89c7f9ffa28028e307aca9b434a89beec06863cf2893887b6b902`;
correction `bce0cc12ef283b0c361818b436c4e190fc4ea5165436f6c36912cfbb6de53368`;
artifact `4c935271a0b2723b58b02cc34d82a1e7d5757f7cd36b7707f400c2b28893668a`.
The contract and F.6.3 source pin all five current F.1 raw hashes too.

## Source and runtime path

Only `src/cuts/pion_hgcer_method_a_parallel_full_procedure.py` and its focused
test change. For Q4p4W2p74, exact canonical current F.1 files and exact
`-current-baseline-candidate.json` F.3/F.4 basenames are loaded without search
or fallback. Unsupported kinematics fail closed. Existing F.5 validation gets
the explicit F.6.3-only F.4 record, including inherited F.3/F.1 pins and farm
validation HEAD. The existing F.4 builder gets the explicit candidate-F.3
record with zero farm head `0000000000000000000000000000000000000000`.
That sentinel reproduces candidate construction only; it is not accepted farm
F.3 authority. Historical F.4/F.5 global authorities remain untouched.

Exact persisted/recomputed F.4 equality and fingerprint equality remain.
An end-to-end synthetic test exposed that F.4's sanitized calculation rows
omit `analysis_MM` and `analysis_t`. F.6.3 now joins those coordinates by event
identity from the same fingerprint-validated raw F.1 population, solely for
the unchanged live-cache parity gate. No shared calculator changes.

The existing yield producer still invokes parity before private filling,
uses transient positive/finite factors for `w0 -> w0*C`, and retains only
aggregate candidate provenance. Baseline public output, parent-t mathematics,
child normalization policy, yields, Method B exclusion, and production remain
unchanged. No event correction is persisted or production promotion performed.

## Deterministic local validation

Python 3.12.10: required four-file `py_compile` passed.
Focused F.6.3 tests: 28 passed; F.4: 8 passed; F.5: 8 passed;
tracked debug launcher: 9 passed; plotting: 84 passed, 18 skipped.
Plotting skips: four PyROOT-dependent tests and fourteen explicitly retired
cumulative presentation contracts. Exact test names/reasons accompany the
actual-diff review bundle. Synthetic tests exercise real F.4/F.5 validation and
reconstruction, identity tampering, the reconstruction sentinel, parent closure,
raw-F.1 coordinate parity, and unchanged public/private branch regressions.
They do not prove reproduction of the supplied farm bytes or farm integration.

Manifest write/check, ordinary memory health, frozen-file identity verification,
and actual-diff audit passed. No farm, local `main.py`, commit, or push occurred.

## Required next gate

Following independent actual-diff review, user commit/push, and pushed-state
synchronization, use the existing tracked Q4p4W2p74 / Left / lowe `-d` debug
path. Preserve paired low/high preflight and stop before high-epsilon full
processing. Inspect baseline/reweighted missing-mass spectra, per-t and
per-(t,phi) yields/deltas, parent-t preservation, procedure-PDF representation,
and absence of baseline production mutation. CURRENT owns the sole NEXT.
