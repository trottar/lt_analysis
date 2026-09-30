# F.4.Refresh.1 — current-baseline Method-A authority comparison

## Status

`CLOSED / RUNTIME VALIDATED` — the [supplied farm comparator evidence](../evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md) closes only the detached current-baseline comparison gate. F.2/F.3 scientific payloads matched exactly; F.4 did not; the first changed stage was `F4`. The earlier independent ChatGPT actual-diff/source-runtime-path review of repaired cumulative `kaonlt_review(20260930-090154).diff` PASSED at fixed `test` HEAD `a8bd4dc4e0990dd45cfade0079dc1ce40bdd36ce`. No accepted authority changed, and no refreshed F.4, F.5, F.6.3, E.8.4, or production runtime acceptance follows.

## Diagnostic boundary

`testing/compare_method_a_current_baseline_authority.py` takes explicit current canonical-five F.1 and accepted F.2/F.3/F.4 JSON paths. It computes raw input hashes, invokes the existing public F.2, F.3, and F.4 artifact builders, and hashes candidate F.2/F.3 bytes in the repository writer format. A candidate F.3 authority record is passed only in memory to the F.4 builder. It compares provenance-excluded scientific payloads exactly, reports the first changed stage, and reports each F.4 parent metric and setting/global change maxima. Scientific mismatch is diagnostic output, not a nonzero exit or promotion verdict. Only the explicit output JSON may be written.

The [direct Fix.4 evidence](../evidence/e8-4-fix4-left-lowe-runtime-closure-and-f4-baseline-divergence.md) closes the narrow cache semantics gate and establishes unchanged Left/lowe F.3 projections but changed F.4 baseline numerics. The later [comparator evidence](../evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md) resolves the diagnostic stage question. F.6.3 and E.8.4 remain `SOURCE REVIEWED`; the fresh Left/lowe runtime gate remains `BLOCKED`.

## Local checks and next action

The focused synthetic tests cover exact canonical inventory, strict JSON, raw bytes, public-builder orchestration, writer-format hashes, in-memory override with an all-zero non-farm-head sentinel, scientific projections, parent changes, stage precedence, and output/input write safety. Local suites passed: comparator 6, F.2 14, F.3 10, F.4 8, F.6.3 21, and memory-health 35 tests. Python compilation, memory manifest write/check, bootstrap, and `git diff --check` passed. Memory health returned its existing soft-size warning for `CURRENT.md` (over 8,192 bytes). No ROOT/PyROOT, farm comparator, full analysis, or production change occurred.

## Fix.1 — nested F.4 diagnostic delta completeness

Independent ChatGPT actual-diff/source-runtime-path review of `kaonlt_review(20260930-000747).diff` found one narrow defect: real `source_diagnostics` and `canonical_phi_diagnostics` are lists of dictionaries, but the candidate emitted each list as one atomic equality result. Nested numerical accepted/candidate/absolute/relative differences were therefore absent. Fix.1 adds index-preserving list recursion, explicit missing-side entries for unequal lengths, and list traversal for numerical leaves. Focused fixtures now use the real list shapes and assert source and canonical-phi nested deltas. This changes diagnostic report completeness only; the builder chain, authority boundary, scientific payload comparison, and production science remain unchanged.

F.4.Refresh.1.Fix.1 is `CLOSED / RUNTIME VALIDATED` only as part of the detached farm comparison gate that used its list-recursive output. It was `SOURCE REVIEWED` in the passed `kaonlt_review(20260930-090154).diff` cumulative review. ChatGPT did not run the Codex-reported tests: focused comparator 8, F.4 8, F.6.3 21, and memory-health 35 tests passed locally with no required-suite skips. The unchanged original F.2/F.3 candidate suites previously reported 14/10 passing tests. No current-baseline authority acceptance follows.

F.4.Refresh.2 current-baseline candidate materialization is the next `ACTIVE` development gate; its candidate awaits independent source review. Any accepted-authority decision remains later.
