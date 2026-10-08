# Detached F.1 baseline reproducibility diagnostic

2026-10-08: `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`, conditional on
ChatGPT actual-diff review, under the
[contract](e8-4-f1-baseline-reproducibility-diagnostic-task-contract.md).
Start branch `test`, HEAD/local `origin/test`
`ebf3a0be2f50c788c22defe3601fd5b7fee48644` matched the required source gate.

## SOURCE VERIFIED — implementation

Only two new executable/test files are added. The read-only standard-library
CLI requires explicit current F.1/output paths and raw SHA-256 pins. Reference
F.1 and reviewed F.4 require independent pins. Duplicate keys/nonfinite JSON,
wrong settings/schema/flags, fingerprints, duplicate identities, incomplete
eligible t inventory, invalid geometry and signed contributions fail closed.
The F.1 producer's pure-Python helpers supply setting/schema/hash semantics.
Actual producer Q2/W tokens and explicit numeric fixture values are supported.

The single-setting forensic parser is explicitly not the canonical F.1/F.4
validator: those frozen validators require five settings and scientific
dependencies. It verifies the serialized F.1 v2 fingerprint projection and
F.4's persisted `b=c*w0` integrity formula/default 1e-12 scaled tolerance;
reference/scalar comparisons use exact equality, without acceptance tolerance.
It imports no ROOT, NumPy/SciPy, launcher or owner and constructs no factors.
Reviewed F.4 parsing validates wrapper/parent fingerprint binding, fifteen
parent identities, full settings and current t geometry, then consumes only
three Left-lowe baseline scalar aggregates. It grants no canonical acceptance.

Inventory distinguishes all application rows from F.4's inside-phi rows,
with signed/absolute `math.fsum`, coefficient/weight/MM summaries and t/source/
phi groups. Reference rows join by `(source_label, entry_index)`; changed
variables and missing/added identities have at most ten sorted samples each.
Matched-row parent/source deltas use reference parent ownership. Alternative
raw F.1 identities are historical comparisons, never reviewed authority.
Output has no full event/factor arrays, timestamps or ambient filesystem paths.
Publication uses a flushed temporary file and atomic no-clobber hard link;
preexisting outputs, symlinks, path/inode aliases and concurrent targets remain
untouched. Inputs are hash-rechecked before publication; temporary files are
removed. No partial diagnostic is published on validation failure.

## Deterministic local validation

Linux Python 3.12.3: 74 distinct tests passed, zero failures/errors/skips:
19 diagnostic tests plus 55 unchanged F.1 acceptance-contract, detached
comparator, F.4 parent-preserving and F.6.3 parallel-procedure regressions.
Regression dependencies came from an existing isolated temporary directory;
the CLI and its fresh-process import test require only the standard library.
Synthetic contract-shaped fixtures exercise cancellation/zero/outside-phi,
row reordering, coefficient-only versus weight-only changes, MM/phi/t changes,
missing/added IDs, exact ~2.6627e-11 absolute-parent mismatch, missing reviewed
F.1 limits, byte-identical clean runs, input immutability and rejection paths.
No scientific acceptance function is patched. Syntax compilation uses temporary
bytecode outside the repository. Manifest/health/bootstrap, whitespace and
preservation checks are reported with the final review bundle.

## Evidence boundary and next gate

The contract owns the supplied newer failed pre-analysis F.4 comparison at
the start source; it is diagnostic only. Original reviewed Left-lowe F.1 was
not found in the searched locations. Full availability, event-level cause,
upstream fit/template numerical variation and original MM weight-bin assignment
remain NOT VERIFIED. No local supplied-artifact diagnosis, ROOT/farm/PDF
validation, tolerance change, pin refresh or production promotion is claimed.

CURRENT is consolidated and replaces the stale owner retry with the detached
diagnostic evidence gate conditional on actual-diff acceptance, user publication
and pushed-state synchronization. A diagnostic invocation/readiness gate must
be independently reviewed before a farm command. Canonical-five full runtime,
final E.8/F.6.4 and absolute-SIMC interpretation remain `BLOCKED`; existing
closed phases keep their scopes. Method A remains detached/non-production;
Method B remains diagnostic and numerically excluded. All frozen source,
unrelated tracked bytes, supplied contract, index and source refs are preserved.
The complete temporary review bundle includes the unchanged supplied contract.
No staging, commit, push, ref update or farm operation occurred.

## SOURCE VERIFIED — reviewed-F.4 authority pin repair

The pure comparison builder and CLI now require reviewed F.4 SHA-256
`79e7ceda7221cbeeead4ed5bc306b0e0e670741a27beaa980e22349c555e96d7`.
The CLI also retains independent raw-byte verification. A self-consistent
alternative wrapper cannot acquire reviewed-authority status using its own
hash. Calculations, parser semantics and exact comparisons are unchanged.
The combined suite passes 78 tests (23 diagnostic plus the same 55 unchanged
regressions), with zero failures/errors/skips, including fixed-pin acceptance,
alternative-identity rejection and CLI fixed-pin/raw-byte checks. Positive
builder/CLI tests temporarily substitute only the authority constant with an
actual synthetic fixture digest and restore the production pin afterward.
Different self-consistent parent payloads fail both their own-digest authority
check and, in the CLI, a falsely declared fixture pin's raw-byte check. Parser
negative tests also isolate fixture pins so schema/geometry checks still run.
The genuine reviewed F.4 bytes are not a local fixture; no successful real-
artifact CLI or farm validation is claimed. No acceptance function is patched.

The four new files retain DrvFS executable-mode warnings: the two diagnostic
Python files, this phase record and the supplied contract. Per explicit user
authorization, Git mode normalization is `DEFERRED` to user publication;
these warnings are hygiene, not implementation-test failures. No filesystem
permission change or Git-index operation is performed by this repair.
