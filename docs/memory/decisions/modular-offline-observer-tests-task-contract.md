# Modular offline observer diagnostics — implementation contract

Audit record: [Work-ZIP review](../investigations/work-zip-modular-testing-audit-2026-11-04.md).

## Starting state and purpose

Start from `test` HEAD `302bccddcc490237aa7128879b4b04ec979260c6`. Establish live HEAD/worktree, read the repository startup core in order and this contract. Preserve all pending memory edits and both untracked reviews. If scientific source moved, inspect the gate-relevant diff before implementation; do not reset, stash or overwrite it.

User intent: stop rewriting per-run Work programs. Add a reusable, input-driven offline diagnostic in the existing `testing/` area after the completed Work-code audit. This is one testing change, not a production correction or new scientific phase.

The audited Work scripts are immutable references. Original numerical results were independently rechecked; their verifier has four demonstrated blind spots. Implement against explicit invariants rather than importing their top-level workflows.

## Allowed implementation files

New files only:

- `testing/offline_diagnostics/__init__.py`
- `testing/offline_diagnostics/evidence.py`
- `testing/offline_diagnostics/observer.py`
- `testing/offline_diagnostics/reporting.py`
- `testing/offline_diagnostics/profiles/observer_shapes_20261104.json`
- `testing/run_offline_diagnostics.py`
- `testing/test_offline_diagnostics.py`
- `testing/OFFLINE_DIAGNOSTICS.md`

Warranted updates: `docs/memory/CURRENT.md`, this audit's investigation record,
the existing baseline-pion investigation if its measured limits change, and
`docs/memory/manifest.json`. Preserve other pending memory bytes and history.
No existing source/test/renderer/owner file is in the implementation allowlist.
Stop if a broader source repair or a dependency change is necessary.

## Interface and behavior

One documented CLI provides `inspect`, `observer-shapes`, and `verify` operations.
Accept repeatable `--input` paths, an explicit new `--output` destination where
required, optional setting/epsilon selections and a diagnostic `--window LOW HIGH`.
Window defaults come from compatible retained metadata. A selected subwindow
must lie inside the retained allcuts selection; it cannot change that selection.
Support a named optional regression profile; source code must not contain this
cohort's filenames, hashes, expected 97/15/135 counts, Q/W or 3-by-9 layout.
An optional renderer must not be imported by inspect/verify or numerical imports.

Adapters initially support `pion_t_template_comparison/v1` and
`pion_t_phi_template_comparison/v1` as distinct roles. Discover input members by
supported schema and validate metadata, not only a name substring. Bind matching
parent/child identities; the older low child lacks top-level setting metadata.
Use manifest/pairing evidence for epsilon/Q/W, preserving declared versus proven
identity. Do not invent missing metadata, trust filename-only identity, or mix
cohorts because settings or source HEAD happen to match. Record raw archive,
member and manifest hashes plus producer/source/coordinate roles separately.

Use strict JSON (no duplicate keys or nonfinite constants), validated increasing
edges, finite equal-length arrays, nonnegative retained variances, explicit
flows, source/coordinate consistency and unique canonical cell identities.
If both errors and Sumw2 are present, verify consistency; explicitly choose a
supported available field rather than evaluating an eager fallback. Reject
unsupported schema/policy. Missing optional evidence is an unavailable result,
never a zero measurement or guessed production acceptance.

Validate inputs before calculating or creating outputs. Inputs and existing
results are read-only. No implicit workstation/farm discovery, archive extraction,
network operation, subprocess analysis launch or source import from `src/`.
Use explicit exceptions for input gates. All modules must be import-safe with
no reads, directories, logging, plotting or writes triggered by import.
Outputs must be exclusively created, cannot alias evidence, and cannot overwrite
the checkout or an existing result. Define structured failure/partial-result
semantics so an incomplete directory cannot masquerade as a successful report.
Missing optional plotting dependencies must not break inspect/verify; record a
render failure separately from numerical status. Distinct nonzero CLI statuses
must distinguish invalid input/usage from failed numerical verification.

## Frozen scientific arithmetic

Observer = retained allcuts kaon histogram minus its retained Local-t transferred
pion background. Check the retained application and quadratic Sumw2 relationship.
Do not refit, reconstruct a final producer, change weights or import ROOT.
Use only bins wholly contained in the selected diagnostic window. Derive their
number and effective boundaries from edges; report unresolved boundary strips.
No fractional bins, interpolation, arbitrary rebinning or hidden boundary tuning.

Shape normalization is diagnostic unit-sum only for positive signed sums.
Preserve negative bins; distinguish empty, zero, invalid, unsupported and
nonpositive-normalization states. Conditional diagonal counting covariance and
its normalization-induced bin covariance must match the audited formulas;
do not pretend fitted/training/systematic covariance is supplied. A supported
full covariance must be consumed explicitly or rejected as unsupported, not
silently reduced. No new fitted statistic, p-value or physics tolerance.

Carry the actual signal scope from the producer and verify it. Current references
are setting fallback, not matched-cell acceptance. Canonical parent data sum
children with outside-phi excluded. Use one shared identical setting template
once; distinct normalized child templates need genuine combination weights,
otherwise parent SIMC is unavailable while valid children remain processable.
Component attribution, when supported by the recorded final numerator policy,
must reproduce each component and close to retained total; if policy metadata
is insufficient, keep total observer available and attribution unavailable.

Reports consume the validated result payload. Centroid, signed RMS, peak-bin,
tail partitions, L1 and any retained descriptive Q must be recomputed where
claimed and labeled with effective domain and uncertainty scope. No hard-coded
finding prose or promotion based on a smooth plot. Keep historical baseline,
detached observer and final producer roles distinct. No ratios, Method A/B,
source-aware refitting, Outcome B, SIMC iteration, L/T or production adoption.

## Required tests and acceptance

Use existing unittest conventions and Linux/WSL development; no new test framework,
environment manager, service, CI workflow or dependency upgrade. NumPy is the
numerical dependency; Matplotlib is optional rendering only. Synthetic fixtures
test engineering behavior and must never be presented as farm/scientific data.

Positive/regression cases:

- The original two farm ZIPs and the reviewed Work package reproduce the existing
  observer arrays, principal metrics, covariance and statuses within documented
  floating-point tolerances. Original archives/source/results remain unchanged.
  Cohort-specific expectations live only in the profile or fixtures.
- Different archive names/locations and a validated alternate setting/kinematic,
  two-t/four-phi geometry and different MM grid work without source changes.
- Settings may be processed independently; an unavailable cell or parent does
  not erase supported independent comparisons. Missing requested settings are
  explicit and cannot be reported as completed.
- Plotting is absent from import/inspect/verify and optional for numerical work.

Negative cases must reject each demonstrated mutation: changed scalar metric,
false signal scope, scaled SIMC covariance, swapped component labels. Also cover
duplicate/missing cells, mismatched edges or cohorts, malformed/duplicate-key/
nonfinite JSON, negative variance, unsupported policies, missing variance,
input/output aliases, existing outputs, incompatible windows and incomplete runs.
Required structured unavailable cases include empty cells, nonpositive sums,
missing signal support and unweighted aggregation of different unit templates.

Metamorphic/oracle checks: equivalent-unit scaling of contents and variances
leaves normalized shapes/covariance invariant; canonical child order does not
change parents; normalized covariance has its unit-sum null direction and is
positive semidefinite within arithmetic tolerance. Independently check covariance
with an expanded expression, and centroid errors from original centered-bin
variances. Do not rely only on a producer/verifier pair sharing the same bug.

Run focused new tests plus the existing transfer-validity and SIMC-parent-support
suites; existing baseline observation is 73 passing tests. Use an explicit external
artifact path for actual-evidence regression; absent private inputs must be an
explicit skip, never a fabricated PASS. Test the supplied actual archives in the
authorized local environment. ROOT/farm and Python 3.9 compatibility claims need
direct evidence; report the interpreter actually tested.

## Review and hard stop

Inspect actual diff and complete new files; verify all out-of-scope source,
existing tests, pending memory changes, original ZIPs, Git index and refs remain
preserved. Regenerate/check memory manifest and health. Supply a fresh verified
ZIP containing complete actual diff/additions, test receipts, input identities
and reproduction commands with dynamic input arguments.

Stop after local implementation and review. Do not commit, push, stage, run Codex
as another process, dispatch CI or run JLab. Local tests establish local tool
behavior only. No production physics, scientific gate or farm closure changes.
