# Work-ZIP code audit before modular testing infrastructure

Status: SOURCE REVIEWED for this targeted audit; narrow repair required before
reusable adoption. Infrastructure was implemented locally November 5 under the
[bounded contract](../decisions/modular-offline-observer-tests-task-contract.md).
No testing implementation or production source changed during this preparation.

## User request and sequence

November 4: update memory first, then audit actual code from generated Work ZIPs,
then develop a modular testing environment that accepts new inputs instead of
rewriting scripts. Opening checkpoint preserved five pending memory edits and
two untracked reviews at test HEAD `302bccddcc490237aa7128879b4b04ec979260c6`.
Scientific NEXT and original six-step gates remain distinct from this preparation.

## Actual scope

Inventory: 49 workspace KaonLT ZIPs, 260 non-snapshot Python members and 231
distinct hashes. This is not full numerical review of every program. Targeted
source/evidence review: 13 shipped scripts in the five-setting observer,
superseded low-only blocker, actual Combined54, Part3 evidence and five-setting
Method-A-sensitivity ZIPs. Those predecessor studies were not rerun or promoted.
Latest observer ZIP SHA-256:
`83ee368214dedd5fd040c5f0f7768e2856c1e7a527ecf84091c93115c4fb7e2b`.
The actual test directory is `testing/`, with 94 test_*.py files observed.

## Findings and measured checks

The exact shipped observer verifier returned PASS for four isolated altered
copies: a reported centroid changed to 9999 MeV, false signal-scope labels,
SIMC covariance multiplied by ten, and exchanged pi_n/pi_sidis component labels.
The unaltered control also passed. This establishes verification-coverage gaps,
not a wrong original scientific number. Original archives/results were preserved.

A separate direct-original recomputation checked individual components, signal
scope, principal scalar metrics/CSV and expanded covariance: PASS, 5,986 checks
over 150 records / 112 computed shapes. The current 97 child + 15 parent observer
values remain numerically supported under their conditional model. No calibrated
Q interpretation, physical-purity, matched-cell SIMC or producer/farm closure.

Work code also has import-time IO/calculation/rendering, run-specific paths and
3x9 dimensions, fixed-count verifier assertions, overwritten fixed outputs,
hard-coded narrative, split provenance gates and an eager histogram-error fallback.
The variance-only reader probe raises KeyError even with stored Sumw2; current
farm records contain both fields, so this is a generalization limitation.

Two existing ROOT-independent suites actually ran under WSL: transfer-validity
and SIMC-parent-support, 73 tests passing, zero skips (unittest 1.902 s). No full
suite, CI, ROOT, farm, Python3.9-runtime or existing-plot closure is inferred.
The repository already has argparse/pure-analyze, strict input checks, explicit
unavailable states, provenance and exclusive report-writing patterns. Extend
them; do not replace or silently broaden pinned historical diagnostic modules.

## Bounded implementation NEXT

Implement one input-driven offline observer entry point and small import-safe
package in testing/: evidence adapters, pure numerical kernel and optional
rendering. Infer valid settings/t/phi/grid from metadata; an optional profile
owns archived hashes/counts. Inspect, observer-shapes and verify consume explicit
paths and new output destinations. Low child metadata must bind to the matching
parent/manifest; never filename-only guessing. Reject all four demonstrated
mutations and malformed/mismatched inputs. Add alternate-layout and covariance/
unit-scaling tests. Preserve negative/empty/nonpositive states and fallback labels.
No source-aware fit, Method A, Outcome B, production correction or farm owner
change. Source allowlist, positive/negative/regression tests and hard stop are
in the single contract. Current production and all 338 snapshotted source/test
files are unchanged; pending memory and review files are preserved.

## Deliverable and limits

`KaonLT_Work_ZIP_Modular_Testing_Audit_20261104.zip` contains actual extracted
code, inventories, probe and independent-recheck receipts, audit, contract and
complete memory diff/additions. Large original science ZIPs are referenced by
identity rather than duplicated. Infrastructure has not yet been implemented.
No commit, push, source edits, farm execution, CI or Codex subprocess occurred.

## November 5 local implementation and review

Status: SOURCE REVIEWED for the detached local tool, not farm/producer physics.
Eight new allowlisted files implement import-safe evidence, observer and reporting
modules; inspect/observer-shapes/verify CLI; an external regression profile;
focused tests and usage documentation. New inputs, paths, settings and t/phi/MM
layouts do not require code edits. Older child metadata binds to its parent and
manifest. No existing source, test, renderer or farm owner was modified.

Thirty new tests passed with the original Work evidence explicitly supplied;
all 73 existing adjacent validator tests passed under WSL. Actual archived
97 child + 15 parent observers were reproduced, including original arrays,
principal metrics, covariance and statuses. All four prior verifier blind spots
are rejected. Tests also cover alternate two-t/four-phi kinematics/grid,
source/axis/schema/variance failures, unavailable states, CSV corruption,
partial numerical output, exact integer identities and equivalent-unit behavior
including detection at tiny input units. Repeated source/cell counts are profile
data, not kernel constants. Optional rendering produced a local five-page smoke
test; it is not ROOT procedure-PDF validation.

Numerical tests used existing Windows Python 3.12.10 / NumPy 2.3.0. Default WSL
Python lacks NumPy: metadata-only inspection worked there for both original ZIPs
and all five settings, and the existing 73 standard-library tests ran there.
No dependency was installed or changed. Native WSL numerical execution and
Python 3.9 runtime remain NOT VERIFIED; syntax parsing alone is not runtime proof.

Inputs, production source, existing tests, prior memory/accountability changes,
untracked reviews, Git index and HEAD were preserved. Manifest/diff/health and
review-package checks were completed. Current scientific gates remain unchanged.
Deliverable: `KaonLT_Modular_Offline_Testing_Implementation_20261105.zip` with
complete files/diffs, actual receipts, profile, source review and dynamic usage.
No Codex subprocess, hosted CI, ROOT, farm, staging, commit or push occurred.
The user queued the immediate normalized-MM handoff after this task; CURRENT
records that subsequent operation without promoting observer evidence.
