# F.4.Refresh.1 — current-baseline Method-A authority comparison

## Status

`SOURCE REVIEWED` — independent ChatGPT actual-diff/source-runtime-path review of the repaired cumulative `kaonlt_review(20260930-090154).diff` PASSED at fixed `test` HEAD `a8bd4dc4e0990dd45cfade0079dc1ce40bdd36ce`. The reviewed comparator requires the exact canonical-five current F.1 inventory and raw-byte hashes, rebuilds F.2/F.3/F.4 through their public builders, hashes writer-format candidate F.2/F.3 bytes, and passes candidate-F.3 authority only in memory for diagnostic F.4 construction. It compares scientific payloads with declared provenance exclusions and reports recursive F.4 parent/source/phi numerical deltas. It changes no accepted authority, scientific/runtime source, production physics, or Method-A promotion. The comparator has not run on farm artifacts; no candidate scientific outcome or F.6.3/E.8.4 runtime acceptance is claimed.

## Diagnostic boundary

`testing/compare_method_a_current_baseline_authority.py` takes explicit current canonical-five F.1 and accepted F.2/F.3/F.4 JSON paths. It computes raw input hashes, invokes the existing public F.2, F.3, and F.4 artifact builders, and hashes candidate F.2/F.3 bytes in the repository writer format. A candidate F.3 authority record is passed only in memory to the F.4 builder. It compares provenance-excluded scientific payloads exactly, reports the first changed stage, and reports each F.4 parent metric and setting/global change maxima. Scientific mismatch is diagnostic output, not a nonzero exit or promotion verdict. Only the explicit output JSON may be written.

The [direct farm evidence](../evidence/e8-4-fix4-left-lowe-runtime-closure-and-f4-baseline-divergence.md) closes the narrow Fix.4 cache semantics gate and establishes unchanged Left/lowe F.3 projections but changed F.4 baseline numerics. It does not establish the comparator's outcome. F.6.3 and E.8.4 remain `SOURCE REVIEWED`; the fresh Left/lowe gate remains `BLOCKED`.

## Local checks and next action

The focused synthetic tests cover exact canonical inventory, strict JSON, raw bytes, public-builder orchestration, writer-format hashes, in-memory override with an all-zero non-farm-head sentinel, scientific projections, parent changes, stage precedence, and output/input write safety. Local suites passed: comparator 6, F.2 14, F.3 10, F.4 8, F.6.3 21, and memory-health 35 tests. Python compilation, memory manifest write/check, bootstrap, and `git diff --check` passed. Memory health returned its existing soft-size warning for `CURRENT.md` (over 8,192 bytes). No ROOT/PyROOT, farm comparator, full analysis, or production change occurred.

## Fix.1 — nested F.4 diagnostic delta completeness

Independent ChatGPT actual-diff/source-runtime-path review of `kaonlt_review(20260930-000747).diff` found one narrow defect: real `source_diagnostics` and `canonical_phi_diagnostics` are lists of dictionaries, but the candidate emitted each list as one atomic equality result. Nested numerical accepted/candidate/absolute/relative differences were therefore absent. Fix.1 adds index-preserving list recursion, explicit missing-side entries for unequal lengths, and list traversal for numerical leaves. Focused fixtures now use the real list shapes and assert source and canonical-phi nested deltas. This changes diagnostic report completeness only; the builder chain, authority boundary, scientific payload comparison, and production science remain unchanged.

F.4.Refresh.1.Fix.1 is `SOURCE REVIEWED` in the passed `kaonlt_review(20260930-090154).diff` cumulative review. ChatGPT did not run the Codex-reported tests: focused comparator 8, F.4 8, F.6.3 21, and memory-health 35 tests passed locally with no required-suite skips. The unchanged original F.2/F.3 candidate suites previously reported 14/10 passing tests. This review is source/runtime-path inspection, not farm comparator execution or current-baseline authority acceptance.

The next gate is user-controlled commit/push after independent final pre-push reconciliation review; farm comparison and any accepted-authority decision remain later gates.
