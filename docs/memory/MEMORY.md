# Durable KaonLT project knowledge

## Authority and status

docs/memory is the authoritative project development record; native Codex
memory is supplemental only. Start substantial work by reading CURRENT.md,
this file, and handoffs/CURRENT_HANDOFF.md; inspect relevant deep records; then
compare their source identity with live test.

For implementation: current source > exact newer source/diff > newest
authoritative handoff > older history > inference. For runtime: fresh farm
evidence > reviewed authoritative farm handoff > source evidence > history >
inference. Use only CLOSED / RUNTIME VALIDATED, SOURCE REVIEWED, DEVELOPMENT
COMPLETE, FARM VALIDATION PENDING, ACTIVE, DEFERRED, BLOCKED, and NEXT for
work-state claims.

## Source-identity semantics

- Record the reviewed analysis-source commit, the repository HEAD observed at
  the last reconciliation, and the farm-evaluated commit as separate values.
  A source review may be anchored to an analysis-source commit older than the
  observed repository HEAD.
- An observed repository HEAD is a timestamped checkout observation, never a
  permanent assertion that the same commit remains live. Every substantial task
  must establish the actual live `test` HEAD and working-tree state before work.
- Record a farm-evaluated commit only when the applicable farm artifact/evidence
  identifies it; otherwise state that no farm-evaluated commit is recorded.
- Documentation/memory-only commits after a reviewed analysis-source commit do
  not by themselves invalidate that source review. Reconcile the intervening
  diff and re-review only if the relevant analysis source or test scope changed.
- F.1.Fix.5 validation profile/bundle schema v4 pins
  dc4fc6283001739a487ec80068f951b0e388cae6. Its source-identity exception is
  deliberately narrow: three exact collector/profile/test files and the
  docs/memory/ prefix only. It must continue to reject tracked AGENTS.md and
  all other src/ and testing/ paths. The v4 collector validates the producer's v2 training and
  application populations separately; local collector checks are not farm
  evidence.

## Production and scientific boundaries

- Preserve the separate scientific owners for random subtraction, slow-proton
  PID contamination, pion-production background, HGCer diagnostics, SIMC,
  yields, and cross sections.
- Preserve random/dummy -> frozen binning -> slow proton -> pion subtraction.
  A setting-wide K Lambda gate controls whether proposed proton weights become
  applied; do not partially commit per-t results.
- Pion alignment comparisons use fixed evaluation envelopes. HGCer diagnostics
  consume the frozen baseline rather than redefine it.
- Preserve cuts, templates, priors, component definitions, normalizations,
  binning, efficiencies, acceptance, L/T separation, and uncertainties unless
  a narrow approved contract owns the change.

## HGCer boundaries

- Method A and Method B remain independent. Method B uses frozen Phase-A
  records and same-canonical-t relative closure: no Method-A numbers, cross-t
  pooling/interpolation, or absolute neutron normalization.
- Adaptive Method B is DO NOT PROMOTE: legacy Method B remains Phase-D context;
  neither changes production pion subtraction.
- Persisted diagnostics require producer -> serializer/checkpoint ->
  checkpoint-first payload -> consumer -> renderer review.
- Phase D compares frozen A/B states. Phase E is presentation-only and cannot
  recompute a diagnostic or construct a correction.
- Method B is diagnostic/cross-check/historical comparison only. Method A,
  which excludes NPE=0, is merely the future candidate numerical HGCer input
  and remains a positive-response relative diagnostic.
- Method-A response training and downstream application are separate
  populations. Train from prompt/noRF/nommcuts NPE>0 records; application uses
  the authoritative physical NPE>2 cache. Keep separate provenance and
  fingerprints. Future response application is event-level and parent-t
  normalized; never renormalize (t,phi) children separately.
- F.1-F.5 are detached. No production promotion before a separately validated
  and explicitly approved F.6 decision.

## Current evidence boundary

Phase C five-setting closure and E.3.Fix.2 Left-low are the only recovered
phase-level CLOSED / RUNTIME VALIDATED results. F.1 Fix.3/Fix.4 evidence closes
mechanical regressions only. Live F.1.Fix.5 is SOURCE REVIEWED, not farm
validated. See evidence/VALIDATION_HISTORY.md and investigations/KNOWN_GAPS.md.

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING — F.1.Fix.5 validation-gate
reconciliation updated only detached collector/profile/test infrastructure.
The immediate user-run evidence step is one targeted F.1.Fix.5 v2 farm gate
for Q4p4W2p74 Left lowe, using the completed v4 collector/profile
infrastructure. Collect fresh checker/validation JSON and bundle provenance,
inspect the F.1 v2 acceptance artifact and relevant rendered F.1 pages, inspect
traceback/log if anything fails, then establish PASS or one coherent repair.
Broaden to Left lowe, Left highe, Center lowe, Center highe, and Right highe
only after that targeted gate passes.
