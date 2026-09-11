# Phase-F.1 source-identity reconciliation

## Question

Can the checked-in F.1 validation collector validly evaluate live test at
b02316f83bf8ef18641fa7217b90c04bcdca10e3?

## Source evidence reviewed

- The profile requires analysis commit
  d656e15761970d7d612bb028d2746d077795e9ad.
- The collector allows only its detached collector/profile test paths after that
  commit and reports non-allowlisted files as an identity error.
- The range from the required commit to live head contains later F.1 analysis
  changes in src/cuts/full_background_subtraction_plots.py,
  src/cuts/pion_hgcer_method_a_acceptance_contract.py, and src/cuts/rand_sub.py,
  introduced by the Fix.5 source change
  dc4fc6283001739a487ec80068f951b0e388cae6, as well as validation/memory
  material.

## Conclusion

BLOCKED — The current profile cannot provide a valid farm-review bundle for
live F.1.Fix.5 because it intentionally rejects that source identity. This is
a static source review, not an attempted farm run and not a claim that the
current analysis behavior is incorrect.

## Preserved boundary and next action

The collector remains detached: it collects declared artifacts and
source/checker evidence; it does not create physics or certify runtime
acceptance.

NEXT — Establish a narrow reviewed source-identity contract selecting the
Fix.5 analysis commit and revise/reissue the identity rule for exactly that
source. Review that diff before the user runs a fresh targeted farm gate.
