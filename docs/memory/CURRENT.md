# Current KaonLT development state

Last source-identity reconciliation: 2026-09-11, historical-memory migration
and Fix.5 source review. Repository HEAD values in this record are timestamped
observations, not permanent claims about the live checkout.

## Repository identity at last reconciliation

- Repository: trottar/lt_analysis; branch: test.
- Repository HEAD observed at the last reconciliation:
  23bcbd2f8db05f0a2e6cf399557c70df44326121.
- Reviewed F.1.Fix.5 analysis-source commit:
  dc4fc6283001739a487ec80068f951b0e388cae6.
- F.1.Fix.5 farm-evaluated commit: none recorded; no F.1.Fix.5 farm artifact
  is present.
- The import was prepared at 7cdb7847d72501cd3dc504565fafc604e22a6132.
  It is an ancestor of the observed repository HEAD, not the reviewed analysis
  source commit.
- At this reconciliation, commits after the reviewed analysis-source commit
  through the observed HEAD are documentation/memory-only for the F.1 review:
  they contain repository guidance or memory/archive material, not F.1 analysis
  source or test changes. Such commits do not by themselves invalidate the
  source review; every substantial task must still establish the actual live
  test HEAD before work.

## Supported state

CLOSED / RUNTIME VALIDATED — Phase-C Method-B diagnostic closure for Left
lowe, Left highe, Center lowe, Center highe, and Right highe, at C.Fix.2.3
9a66bc62d20a99172e326e915866877b65ae1e5d and later accepted pre-E.3
e3853655db0809923cbf2326e2f779219128eda9. Adaptive Method B is DO NOT
PROMOTE. See evidence/phase-c-five-setting-closure.md.

CLOSED / RUNTIME VALIDATED — E.3.Fix.2 independent Method-A presentation for
Q4p4W2p74 Left-low only: implementation
eb1710f4739ba6ef14f51419806e9fc5bd53c175 and runtime/bundle HEAD
bf53dac84e1396cfd7e3f4e0234749426bfbdcf4. See
evidence/e3-fix2-left-low-runtime.md.

SOURCE REVIEWED — F.1.Fix.5 reviewed analysis source uses separate
prompt/noRF/nommcuts NPE>0 Method-A training records and authoritative
physical NPE>2 application records, with independent summaries/fingerprints.
Pages 1-4 consume training; page 5 consumes application. The path is detached,
has no Method-B numerical dependency, and makes no correction, estimator, or
production mutation. See phases/phase-f1-method-a-acceptance-farm-gate.md.

ACTIVE — F.1.Fix.5 requires source-identity reconciliation before farm review.
The checked-in profile pins d656e15761970d7d612bb028d2746d077795e9ad and
rejects later Fix.5 analysis files. No F.1.Fix.5 farm artifact is present.

## Frozen architecture

- Random subtraction, slow-proton treatment, pion subtraction, HGCer Method A,
  HGCer Method B, SIMC comparison, yields, cross sections, diagnostics/checkers,
  and presentation retain separate ownership.
- Proton cleaning precedes pion subtraction; K Lambda preservation remains a
  setting-wide production gate. Proposed proton quantities never become applied
  quantities without gate acceptance.
- Method A and Method B are independent. Method B is same-canonical-t relative
  closure and remains diagnostic/cross-check only; it never adjusts pion
  weights. Method A is positive-response relative leakage only, not absolute
  zero-NPE probability.
- Phase E presentation uses frozen upstream records and has no physics
  ownership. F.1 through F.5 remain detached; only an explicit validated F.6
  may promote a production change. Parent-level normalization forbids
  independent (t,phi)-child renormalization.

## Explicit gaps

DEFERRED — final Phase-D farm provenance; final four-setting E.3 closure;
detailed older proton farm artifacts; exact raw paths/hashes for recovered
historical bundles. F.1.Fix.5 farm status remains unresolved. See
investigations/KNOWN_GAPS.md.

## NEXT

NEXT — Create and review one narrow F.1.Fix.5 source-identity reconciliation
contract that selects the reviewed Fix.5 analysis commit and makes the
collector/profile identity rule agree with it. Only then have the user run one
fresh targeted F.1 v2 farm gate with the declared five settings, fresh
artifacts, provenance/checker inspection, and review of all five F.1 pages.
