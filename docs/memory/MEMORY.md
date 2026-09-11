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
