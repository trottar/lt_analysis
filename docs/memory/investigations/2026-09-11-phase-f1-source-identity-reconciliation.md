# Phase-F.1 source-identity reconciliation

## Question

Can the detached F.1 validation collector validly evaluate the reviewed
F.1.Fix.5 analysis-source commit
dc4fc6283001739a487ec80068f951b0e388cae6 while allowing only later
non-analysis maintenance? The live test HEAD observed at the start of this
reconciliation was 7f8fcc02cde86f93e4e99042ba10c621f89dd159; it is a
timestamped observation, not a permanent live-HEAD claim.

## Source evidence reviewed

- Before reconciliation, profile v3 required
  d656e15761970d7d612bb028d2746d077795e9ad and its v1 artifact validator
  rejected a genuine Fix.5 v2 artifact.
- The range from that stale commit through reviewed Fix.5 includes the expected
  F.1 analysis changes in src/cuts/full_background_subtraction_plots.py,
  src/cuts/pion_hgcer_method_a_acceptance_contract.py, and src/cuts/rand_sub.py.
- The final tracked range from reviewed Fix.5 is limited to docs/memory/ plus
  the three exact validation files. Root AGENTS.md is local Codex guidance,
  ignored and removed from tracking; it is not an identity exception.
  Profile/bundle v4 therefore rejects it along with all other paths outside
  that narrow rule. It does not allow broad src/ or testing/ prefixes.
- The v4 detached validator mirrors the serialized v2 contract: it checks
  separate NPE>0 prompt/noRF/nommcuts training and physical NPE>2 application
  records, their summaries, population fingerprints, feature metadata,
  child-assignment projection, full fingerprint, geometry, and fixed detached
  authority/provenance flags. Local py_compile and collector tests passed.

## Conclusion

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING — The v4 profile can now provide
a valid review bundle for reviewed F.1.Fix.5 source, while retaining narrow
identity protection and rejecting v1 artifacts. This is local validation
infrastructure evidence, not an attempted farm run or a claim about runtime
analysis behavior.

## Preserved boundary and next action

The collector remains detached: it collects declared artifacts and
source/checker evidence; it does not create physics or certify runtime
acceptance.

NEXT — The user first runs one fresh targeted F.1.Fix.5 v2 gate for
Q4p4W2p74 Left lowe with the completed v4 collector/profile infrastructure.
Collect fresh checker/validation JSON and bundle provenance; inspect the F.1 v2
acceptance artifact, relevant rendered F.1 pages, and traceback/log if anything
fails; then establish PASS or one coherent repair. Only after a targeted
Left-lowe PASS may validation broaden to Left lowe, Left highe, Center lowe,
Center highe, and Right highe.
