# F.4.Refresh.2.Validation.2 — tracked execution owner

## Status

`SOURCE REVIEWED` — independent ChatGPT actual-diff/source-runtime-path
review of `kaonlt_review(20261001-092421).diff` passed for the repaired
candidate based on
committed `test` HEAD `c406c138285727503b12115177dfb8bc7efcb7fe`.
Remote `test` was observed at that exact HEAD during the independent review;
this is supplied review evidence, not a new Codex network observation.
That review did not itself establish commit, push, or farm execution.

## Scope and delegation

`testing/run_f4_refresh2_materialize_verify_package.py` owns preflight,
ordered invocation, output completion/provenance verification, and exact
returned-ZIP verification. Explicit F.1 canonical-five, accepted F.2/F.3/F.4,
and reviewed-comparison paths are caller supplied. The new owner invokes the
existing reviewed F.4.Refresh.2 materializer CLI for all candidate scientific
construction and delegates detached ZIP collection to the existing bundle-only
wrapper and generic collector. It supplies immutable SHA-256 values for the
five profile-declared outputs. The original candidate added the owner and
focused test to the profile's committed-range allowlist. Validation.2.Fix.1
also admits exactly the two already-reviewed/pushed memory-health paths;
its frozen required analysis commit and artifact inventory remain unchanged.

The owner rejects dirty or unexpected source, unsupported kinematics, missing
inputs, stale/partial outputs, an existing ZIP, mismatched completion manifest,
or a returned ZIP with missing or substituted artifacts. It does not overwrite
or clean failed evidence. Verification checks recorded science-gate facts and
published bytes; it does not recalculate F.2/F.3/F.4 science.

## Evidence boundary and next gate

The exact source-reviewed chain is:

```text
explicit validated inputs
  -> tracked F.4.Refresh.2 execution owner
  -> unchanged reviewed materializer
  -> explicit manifest/hash completion verification
  -> unchanged reviewed package wrapper
  -> unchanged generic collector + F.4.Refresh.2 profile
  -> exact returned ZIP verification
```

Codex-reported local Python compilation passed. Focused owner tests: 11 passed, 0 skipped;
profile: 4 passed; unchanged materializer: 13 passed; collector: 28 passed;
wrapper: 3 run, 1 skipped because local `tcsh` is unavailable; memory health:
36 passed, 0 skipped. Manifest write/check, memory bootstrap, strict
warning-free health, and `git diff --check` passed. ChatGPT did not run the unit
suites; it reviewed the actual cumulative diff, source/provenance path, pushed
identities, and supplied check output. These local tests are
synthetic and do not establish JLab filesystem,
ROOT/PyROOT, materialization, packaging, candidate authority, F.6.3/E.8.4
runtime, production, or Method-A promotion acceptance. Accepted F.2/F.3/F.4
authority and all scientific/runtime source remain unchanged. F.4.Refresh.2
farm execution stays `BLOCKED` pending user commit/push and independent
pushed-state review. Final pre-push memory/status reconciliation was performed
under the [reconciliation contract](f4-refresh2-validation2-final-pre-push-source-review-reconciliation-task-contract.md)
and reviewed in `kaonlt_review(20261001-101444).diff`; that review found only
the [push-stability wording defect](f4-refresh2-validation2-final-pre-push-push-stability-repair.md).
The four substantive
candidate files and implementation contracts remain byte-identical.
After user push and pushed-state review, prepare only the single narrow
`Q4p4W2p74` F.4.Refresh.2 materialize -> verify -> package farm gate with
explicit validated inputs. F.6.3/E.8.4 cannot begin until returned evidence
is reviewed. No farm materialization, packaging, ROOT/PyROOT, full `main.py`,
F.5 authority refresh, or candidate-authority acceptance has occurred.

## Validation.2.Fix.1 repair

Independent ChatGPT review of `kaonlt_review(20260930-211453).diff` found
one source-provenance blocker: omitted post-materializer hardening paths.
The owner architecture itself was not rejected. The narrow repair retains
that architecture. Independent review of the repaired
`kaonlt_review(20261001-092421).diff` then passed for Validation.2 and Fix.1.
See the [Fix.1 record](f4-refresh2-validation2-fix1-post-hardening-source-allowlist.md)
for the exact range/blob audit and deterministic repair checks.
