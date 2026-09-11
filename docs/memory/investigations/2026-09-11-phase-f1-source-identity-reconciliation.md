# Phase-F.1 source-identity reconciliation

## Question

Can the detached F.1 validation collector validly evaluate the reviewed
F.1.Fix.5 analysis-source commit
dc4fc6283001739a487ec80068f951b0e388cae6 while allowing only later
non-analysis maintenance? The live test HEAD observed at reconciliation was
47723e02af03b3844dfb0323560ce39f7dfebe1a; it is a timestamped observation,
not a permanent live-HEAD claim.

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

NEXT — The user runs one fresh targeted F.1.Fix.5 v2 gate for all five declared
settings with fresh artifacts, then inspects the collector manifest, provenance,
checker output, and five F.1 pages before recording a farm conclusion.
