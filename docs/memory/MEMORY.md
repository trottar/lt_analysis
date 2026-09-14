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

Phase C five-setting closure, E.3.Fix.2 Left-low, and F.1.Fix.5 are the
recovered phase-level CLOSED / RUNTIME VALIDATED results.

CLOSED / RUNTIME VALIDATED — F.1.Fix.5 validates the detached v2
dual-population artifacts for all five canonical Q4p4W2p74 settings. Preserve
four identities: reviewed analysis source
`dc4fc6283001739a487ec80068f951b0e388cae6`; runtime-evaluated analysis HEAD
`126fa22c19bd29b9952f55b33ab43d59f9727ef6`; collector reconciliation commit
`31dd034d8404e317863bf0933c51253e1d3deeb8`; and this docs/memory-only
closure commit, which is the F.2 baseline. The owner accepted the C1 clean
collector gate and continuity of the existing JSON/PDF/page-manifest evidence
with the supplied v4 runtime bundle. The C1 range whitespace check covers only
the profile's three explicit validation files; the global worktree check and
independent committed-file identity audit remain strict. F.1 remains detached;
F.2-F.5 remain detached and F.6 is the only possible production promotion.
