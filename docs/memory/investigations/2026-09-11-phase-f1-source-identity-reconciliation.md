# Phase-F.1 source-identity reconciliation

## Question

Can the checked-in Phase-F.1 validation collector validly evaluate the live
`test` checkout at `cf9c804b55a501af6839b8fb3b0358c56853c276`?

## Evidence reviewed

- `testing/pion_hgcer_validation_bundle_profile.json` names required analysis
  commit `d656e15761970d7d612bb028d2746d077795e9ad` and permits only
  `testing/collect_pion_hgcer_validation_bundle.py`,
  `testing/test_collect_pion_hgcer_validation_bundle.py`, and the profile
  itself after that commit.
- `testing/collect_pion_hgcer_validation_bundle.py` checks that the required
  commit is an ancestor, lists committed files in `required_commit..HEAD`, and
  reports non-allowlisted paths as an identity error.
- The live range `d656e15761970d7d612bb028d2746d077795e9ad..HEAD` contains
  later changes to `src/cuts/full_background_subtraction_plots.py`,
  `src/cuts/pion_hgcer_method_a_acceptance_contract.py`, and
  `src/cuts/rand_sub.py`, in addition to collector/profile and durable-memory
  files.

## Conclusion

`BLOCKED` — The current profile/collector will reject the live committed range
as containing unexpected files after its frozen analysis commit. The profile
therefore cannot provide a valid Phase-F.1 farm-review bundle for the current
source identity. This conclusion is a static source review; it is not an
attempted farm run and it does not say that the current analysis behavior is
incorrect.

## Constraints preserved

The collector stays detached from the analysis runtime. Its purpose is to
collect declared artifacts and source/checker evidence, not to create physics,
modify the production workflow, or assert runtime acceptance.

## Exact next action

`NEXT` — Before farm validation, establish a narrow reviewed contract that
selects the analysis commit to evaluate and updates/reissues the identity rule
for exactly that source. Review its diff and only then run one fresh targeted
gate. Keep import of older farm evidence separate from this repair.
