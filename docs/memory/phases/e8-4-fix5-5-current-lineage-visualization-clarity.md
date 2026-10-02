# E.8.4 Fix.5.5 — current-lineage visualization clarity

**Status:** `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`.

## Source and authority

Started on exact `test` HEAD `761fbb6c03d2d7a10bb911cf84e9ba898496fab6`,
parent `71a912bf14c59eee6d52956462aa0aa39c7f7c84`; HEAD remained unchanged
during that local implementation. The [task contract](e8-4-fix5-5-current-lineage-visualization-clarity-task-contract.md)
supplies independent ChatGPT actual-diff and pushed-state review PASS for
Fix.5.4 numerical source. That source is `SOURCE REVIEWED`; its phase remains
`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`. Fix.5.5 was `ACTIVE` during
implementation. Independent ChatGPT actual-diff and pushed-state review have
now passed for pushed `test` source
`df957a6414fc9c515d1f82228517cb801dc90350`, exactly one commit over the
starting identity with the 12 reviewed paths, per the supplied Fix.5.6 contract.
This is `SOURCE REVIEWED`, not actual ROOT/PDF or farm acceptance.

Before editing, traced the existing private producer -> validated/cloned E.8.4
payload -> renderers: child `identity_audit.signed_support` and payload
`current_lineage_aggregate_audit` already own the required values. Numerical
builders, validators, runtime finalizer, production and SIMC ownership are frozen.

## Presentation implementation

One E.8.4-specific helper styles detached displays: baseline blue solid width 2,
open marker 24 size 0.6; Method A magenta dashed width 4; common pion input black
dotted width 1; authorized SIMC black solid step width 2. The generic helper is
unchanged. Pion consequence draws input -> Method A -> baseline. Final MM draws
Method A -> baseline. Compact legends retain curve identity, existing ranges,
Lambda markers and stored yield notes. Synthetic source-authorized SIMC draws
SIMC -> Method A -> baseline; the source-unproven guard is unchanged and draws
only the explicit unavailable placeholder.

Existing yield-summary pages retain three stored-yield graphs and add a bottom
text region displaying every child's stored MM_0/MM_A/delta signed integral
and absolute support. Minimum required support fields are used for legibility;
positive/negative components remain in the stored audit. Matching stored t-parent
B_pi_0/B_pi_A/delta aggregate integrals are labeled current F.6.3 candidate
lineage, Lambda/allcut MM-template aggregate, with the stored non-comparability
reminder for broader F.4 application-population sums. No histogram integration,
child aggregation or new uncertainty is performed there. Setting summary states
stored audit ownership, diagnostic support, separate historical lineage and SIMC
absolute authority requirements. E.8.3 authority, MM, t/phi and cross-reference
visible labels now explicitly distinguish historical accepted F.6.1 from current
F.6.3/E.8.4; persisted arrays and authority/status are unchanged.

## Local verification

Python 3.12.10: six suites, **166 tests run; 148 passed; 18 retained skips**
(4 PyROOT unavailable; 14 superseded historical procedure-tail assertions):

- `testing.test_e8_4_fix5_5_visualization` (8 new deterministic tests);
- `testing.test_full_background_subtraction_plots`;
- `testing.test_e8_4_production_impact_audit`;
- `testing.test_e8_4_fix5_4_identity_audit` (unchanged);
- `testing.test_e8_3_detached_method_a_reweighting_audit`;
- `testing.test_e8_4_fix5_main_order`.

Tests prove baseline-last order and redundant styles for identical curves,
common-input order, authorized SIMC order, exact payload contents/errors/scalars/
child and aggregate identity/SIMC audits before and after rendering, all 24
unchanged page records, all six blocked SIMC statuses/literal reasons with no
SIMC draw, stored-value diagnostic ownership with histogram reads prohibited,
historical labeling, and fail-closed missing-support rendering.

In-memory syntax checks pass for all three changed/new Python files. AST audit
confirms only the ten designated presentation functions changed and three narrow
presentation helpers were added. All numerical builders/validators, generic
style, signed-difference arithmetic, page-record/orchestration functions and
runtime finalizer are unchanged. SHA-256 checks confirm `calculate_yield.py`,
`main.py` and `run_Prod_Analysis.sh` bytes unchanged. `git diff --check` passes.

## Memory and review scope

Manifest write -> check -> ordinary health runs serially with unchanged tools in
a temporary intended-candidate view: tracked memory plus this phase and its
contract. Its regenerated manifest is copied back and every indexed hash/byte
count checked against the working tree. Unrelated user-owned untracked
`workflow-continuity-hardening-task-contract.md` is untouched/untracked and
excluded from the manifest and review additions; preservation SHA-256:
`0ee396c059be03e0a528bef18a1fba1128be12e75b63ef725e9cf9a02dffe421`.
Scoped manifest: 197 indexed memory files; manifest and ordinary health PASS,
zero warnings/hard failures. CURRENT/MEMORY/CURRENT_HANDOFF byte counts: 6994/11873/323.
No index, ignore rule or memory tool changes. Complete `kaonlt_review.diff`
contains tracked changes and all three intended additions relative to starting
HEAD, preserving raw Git diff content.

## Boundaries and NEXT

Local tests do not establish ROOT/PyROOT rendering, actual PDF legibility,
observed cancellation, full `main.py` execution, numerical/visual farm closure,
SIMC absolute units, owner ZIP repair or Method-A promotion. Method A remains
non-production, Method B numerically absent. No factors, contents/errors, yields,
normalization, cuts/binning, page IDs, owner/profile/collector or farm behavior
change. No commit, push, packaging repair or farm execution occurs.

Fix.5.5 scientific/presentation source remains frozen during the separately
contracted [Fix.5.6 owner repair](e8-4-fix5-6-owner-farm-readiness-and-failure-provenance.md).
CURRENT owns NEXT: after Fix.5.6 actual-diff review, user commit/push and
pushed-state review, one narrow tracked-owner Left/lowe farm gate returns fresh
ZIP/status/log/PDF/manifest evidence. No farm execution is authorized here.
