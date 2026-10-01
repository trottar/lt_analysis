# F.4.Refresh.2.Validation.2.Fix.1 — post-hardening source-allowlist repair

## Status and review finding

`SOURCE REVIEWED` — independent ChatGPT actual-diff/source-runtime-path
review of repaired `kaonlt_review(20261001-092421).diff` passed for the
candidate based on
committed `test` HEAD `c406c138285727503b12115177dfb8bc7efcb7fe`.
That review did not itself establish commit, push, or farm execution.
Independent ChatGPT actual-diff review of `kaonlt_review(20260930-211453).diff`
found one narrow source-provenance blocker. The owner architecture itself was
not rejected: its post-materializer committed-range rule omitted two
already-reviewed/pushed repository-memory operational paths. This repair
implements the [Fix.1 contract](f4-refresh2-validation2-fix1-post-hardening-source-allowlist-task-contract.md).

## Provenance audit and repair boundary

The exact non-memory committed inventory for
`141a3d04f9e5d07be21dba14e0e63212c3990bf1..c406c138285727503b12115177dfb8bc7efcb7fe`
is:

```text
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
testing/test_memory_health.py
testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py
tools/check_memory_health.py
```

The hardening paths retain these pushed Git blob identities, also verified
against local worktree bytes:

```text
testing/test_memory_health.py  8f7454dea48c00148dd379a8c141a005c2367d59
tools/check_memory_health.py  2736916c9d1ae63741c1da920ad6be7c14aa0bab
```

Both the profile and owner now admit exactly six post-materializer paths:
the original profile/test, owner/test, and those two historical operational
paths. Only `docs/memory/` remains a permitted prefix.
`required_analysis_commit` remains
`141a3d04f9e5d07be21dba14e0e63212c3990bf1`.
The only owner-source edit adds the two exact allowlist entries. Preflight ->
materializer subprocess -> completion/provenance verification -> existing
package wrapper -> returned-ZIP verification remains unchanged. No scientific
or runtime file, hardening implementation, accepted authority, or production
physics changed. The original Validation.2 contract remains byte-identical.

## Deterministic local checks

Python compilation passed for the owner and focused test. Owner: 11 tests;
profile: 4; unchanged materializer: 13; collector: 28; memory health: 36;
all passed without skips. Wrapper: 3 tests run, 2 passed, with exactly one
existing skip: `tcsh is unavailable on this local host`.
Focused provenance checks accept the explicit six paths plus CURRENT and
reject arbitrary testing/tools helpers, the materializer, and scientific
source. They do not modify or monkeypatch the real hardening source files.

Manifest write/check, bootstrap, the post-update memory-health suite,
strict `--fail-on-warning` health with zero warnings, and
`git -c core.safecrlf=false diff --check` passed locally.
Exact final outputs and the cumulative diff are retained in the fresh
untracked repository-root review bundle returned with this repair.

## Next gate and evidence boundary

The supplied independent review observed remote `test` at
`c406c138285727503b12115177dfb8bc7efcb7fe`. ChatGPT reviewed the actual
cumulative diff and source/provenance path; it did not run the Codex-reported
unit suites. The collector suite's printed synthetic ZIP path is local test
output, not farm evidence. Review confirmed the exact six-path allowlist,
sole `docs/memory/` prefix, frozen source pin and hardening blobs, and
fail-closed rejection of arbitrary testing/tools/materializer/scientific paths.

Final pre-push memory/status reconciliation was performed under the
[reconciliation contract](f4-refresh2-validation2-final-pre-push-source-review-reconciliation-task-contract.md),
with the reviewed profile, owner, focused tests, and both implementation
contracts byte-identical. Independent review of `kaonlt_review(20261001-101444).diff`
found only the [push-stability wording defect](f4-refresh2-validation2-final-pre-push-push-stability-repair.md).
Farm execution remains `BLOCKED` pending user commit/push
and independent pushed-state review. Only after those gates may the single
narrow `Q4p4W2p74` F.4.Refresh.2 operation be prepared; F.6.3/E.8.4 must wait
for review of returned F.4.Refresh.2 evidence.
No farm operation, ROOT/PyROOT, full analysis, authority update, production
promotion, or Method-A promotion occurred. Local synthetic checks do not
establish candidate materialization or packaging on the JLab farm.
