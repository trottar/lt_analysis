# KaonLT E.8.2 Left/lowe isolated runtime scientific-audit owner — task contract

## 1. Task class and authority

This is a **source-changing validation-infrastructure task** for the active
E.8.2 science audit.

It does **not** change KaonLT production physics. It creates one narrow,
tracked, independently reviewable owner/profile/test path that can obtain
admissible Jefferson Lab farm evidence for the already source-reviewed E.8.2
baseline missing-mass/stage-yield chain.

Follow repository-root `AGENTS.md` and the exact five-file startup core before
doing any work.

Authority order remains:

```text
current source/diff
-> newest applicable farm/runtime evidence
-> tracked repository memory
-> older chat/history
```

Do not use this contract to repair the deferred canonical-five
provenance/identity blocker.

### 1.1 Revision after the failed profile gate

This contract **supersedes the prior unimplemented version at the same path**.

The prior contract incorrectly required a v4 validation profile whose declared
`settings` inventory contained only `Left / lowe`. The frozen collector rejects
that shape before collection because `load_validation_profile(...)` requires
the profile's declared settings to equal its fixed canonical-five
`_REQUIRED_SETTINGS` tuple.

`SOURCE VERIFIED` correction:

```text
profile declaration / authorization inventory:
  Left / lowe
  Left / highe
  Center / lowe
  Center / highe
  Right / highe

requested collection scope for this owner:
  Left / lowe only
```

The collector already separates those two concepts:

```text
load_validation_profile(profile_path)
-> profile must declare canonical five

resolve_settings(phi="Left", epsilon="lowe", profile)
-> returns only (("Left", "lowe"),)

collect_validation_bundle(..., phi="Left", epsilon="lowe", ...)
-> packages only the selected Left/lowe setting plus declared global artifacts
```

Therefore **do not modify the frozen collector** and do not weaken its
canonical-five profile invariant. The repair is entirely in the new profile and
owner contract: declare all five required settings, then request only Left/lowe
through the collector's existing `phi`/`epsilon` filter.

---

## 2. Exact starting identity and start gate

Required repository:

```text
https://github.com/trottar/lt_analysis/tree/test
```

Required branch:

```text
test
```

Required starting HEAD and local `origin/test`:

```text
7b8eb3cb231d289d18de9c168960ddddcaa39254
```

Commit:

```text
Fix memory manifest line-ending integrity
```

Before editing, establish:

```text
git branch --show-current
git rev-parse HEAD
git rev-parse refs/remotes/origin/test
git status --short --untracked-files=all
git diff --check
```

The only pre-existing untracked task file permitted is:

```text
docs/memory/phases/e8-2-left-lowe-runtime-scientific-audit-owner-task-contract.md
```

Temporary `kaonlt_review*.diff` files are allowed only as review artifacts.

If branch, HEAD, local `origin/test`, or unrelated worktree state differs,
**STOP as `BLOCKED`**. Do not reset, stash, clean, checkout over, overwrite, or
repair unrelated user state.

---

## 3. Established pushed-state result

The pushed-state synchronization gate at the starting HEAD is accepted for
source synchronization:

```text
SOURCE VERIFIED:
  - remote test HEAD is 7b8eb3cb231d289d18de9c168960ddddcaa39254;
  - the preceding repair commit changed exactly:
      .gitattributes
      docs/memory/manifest.json
      docs/memory/phases/memory-manifest-line-ending-integrity-repair-task-contract.md
  - .gitattributes contains:
      docs/memory/** text eol=lf
  - the regenerated memory manifest matches the pushed LF Git blobs for the
    repaired memory set;
  - CURRENT scientific/status text was unchanged by that push.

RUNTIME VERIFIED:
  - no new runtime claim follows from the memory-integrity repair.
```

Do not reopen that repair unless a new manifest/EOL integrity failure is found.

---

## 4. Source/science audit that motivates this task

ChatGPT re-inspected the live E.8.2 source at the required starting HEAD.

### 4.1 Existing production path

`SOURCE VERIFIED`:

The current active baseline chain is already represented by existing source as:

```text
same authoritative event traversal
-> random subtraction
-> dummy subtraction
-> slow-proton-cleaned production state
-> existing production prune_hist treatment
-> baseline pion subtraction B_pi^0
-> final baseline clean-kaon MM_0(t,phi)
-> existing Y_0(t,phi)
```

The slow-proton correction is applied event-by-event during the initial
production filling. It is **not** a later histogram-level subtraction inserted
between dummy and pion subtraction. E.8.2 therefore uses same-traversal detached
pre-/post-proton bookkeeping to display the conceptual slow-proton stage without
changing the production ordering.

### 4.2 Existing E.8.2 observation path

`SOURCE VERIFIED`:

Current source already captures and validates:

```text
prompt_pre_proton
random_component_pre_proton
after_random_pre_proton
dummy_component_pre_proton
after_dummy_pre_proton
proton_component_removed
after_proton_pre_prune
after_proton_post_prune
pion_input
pion_component_removed
after_pion_final
```

The producer checks:

```text
prompt - random = after_random
after_random - dummy = after_dummy
after_dummy - proton_removed = after_proton_pre_prune
after_proton_post_prune == pion_input
pion_input - B_pi^0 = after_pion_final
```

The final `Y_0`, statistical uncertainty, and existing total uncertainty are
producer-owned results; the renderer does not recalculate them.

### 4.3 Existing production prune boundary

`SOURCE VERIFIED`:

`prune_hist(...)` is unchanged production logic between:

```text
after_proton_pre_prune
and
after_proton_post_prune
```

Its existing behavior can zero high-fractional-uncertainty bins and can reset
an unusable histogram when entries/integral fail its existing criteria.

Therefore the actual Q4p4W2p74 numerical impact of this treatment is a runtime
question; it must not be assumed from source alone.

### 4.4 Dormant empirical residuals

`SOURCE VERIFIED`:

The active profile is literally:

```text
no_empirical_residual
```

and both empirical residual fit scales are forced to zero. Fit 1/Fit 2 are not
an active E.8.2 production stage.

### 4.5 Current scientific unknown

`NOT VERIFIED` until applicable farm evidence is returned:

For Q4p4W2p74 / Left / lowe, the actual per-cell numerical changes from:

```text
random subtraction
dummy subtraction
slow-proton removal
production prune_hist
baseline pion subtraction
```

have not been accepted as a current E.8.2 runtime science result.

No new scientific algorithm is justified by the source audit. What is missing
is a narrow admissible runtime artifact.

The failed implementation attempt adds one additional source-level constraint:
the existing v4 collector must retain its canonical-five profile declaration.
That is a validation-infrastructure invariant, not a reason to broaden the farm
run beyond Left/lowe.

---

## 5. Why the existing farm owners are not the correct E.8.2 gate

Do not run an existing E.8.4 owner as a substitute.

`SOURCE VERIFIED`:

- `testing/run_e8_4_fix5_left_lowe_plot_gate.py` is coupled to Method-A
  candidate identities and E.8.3/E.8.4 page requirements.
- `testing/run_e8_4_fix5_canonical_five_plot_gate.py` is coupled to fresh-F1
  candidate materialization, F.6.3 lineage reproduction, and canonical-five
  E.8.4 acceptance.
- the canonical-five full runtime gate is currently `BLOCKED`;
- its provenance/identity repair is currently `DEFERRED`.

A failed E.8.4/canonical-five gate must not be repurposed as accepted E.8.2
science evidence.

The required new owner is therefore **E.8.2-only in its acceptance contract**,
even though the unchanged analysis may still render later E.8.3/E.8.4
presentation pages.

---

## 6. Objective

Implement one tracked isolated **Left/lowe E.8.2 runtime scientific-audit
owner** that:

1. runs the existing Q4p4W2p74 Left/lowe debug analysis in an owner-created
   detached disposable analysis worktree;
2. uses the established copied-ltsep runtime-isolation/preservation machinery;
3. does not mutate the ordinary farm checkout;
4. does not require Method-A candidate materialization or F.6.3 lineage gates;
5. verifies only the E.8.2 baseline scientific/presentation invariants plus
   baseline analysis identity;
6. packages the minimum artifacts needed for ChatGPT to inspect the E.8.2
   stage spectra, stage-yield text, final MM, and final Y0;
7. fails closed on missing/stale/invalid E.8.2 evidence;
8. leaves E.8.3/E.8.4 acceptance explicitly out of scope.

The implementation must not issue or execute a farm run locally. This task ends
at source readiness and actual-diff review.

---

## 7. Exact allowed paths

Versioned changes are restricted to:

```text
testing/pion_hgcer_validation_bundle_profile_e8_2_left_lowe.json
testing/run_e8_2_left_lowe_scientific_audit_gate.py
testing/test_run_e8_2_left_lowe_scientific_audit_gate.py
docs/memory/phases/e8-2-left-lowe-runtime-scientific-audit-owner-task-contract.md
docs/memory/phases/e8-2-left-lowe-runtime-scientific-audit-owner.md
docs/memory/CURRENT.md
docs/memory/manifest.json
```

No other path may change.

The supplied task contract is an intended tracked file.

---

## 8. Frozen paths and scientific ownership

Everything outside Section 7 is frozen.

In particular, do not modify:

```text
src/
tools/
farm_env/
background_samples/
run_Prod_Analysis.sh
set_SymLinks.sh
testing/collect_pion_hgcer_validation_bundle.py
testing/run_e8_4_fix5_left_lowe_plot_gate.py
testing/run_e8_4_fix5_canonical_five_plot_gate.py
testing/pion_hgcer_validation_bundle_profile_e8_4_fix5.json
testing/pion_hgcer_validation_bundle_profile_e8_4_canonical_five.json
testing/test_e8_2_baseline_stage_audit.py
docs/memory/MEMORY.md
docs/memory/USER.md
docs/memory/handoffs/CURRENT_HANDOFF.md
docs/memory/roadmap/STATUS.md
docs/memory/evidence/
docs/memory/decisions/
docs/memory/investigations/
```

Preserve exactly:

```text
random subtraction
dummy subtraction
slow-proton factors and PID ownership
production prune_hist behavior
baseline pion w0
pion components/fits/windows/templates/priors
no_empirical_residual
SIMC
yield formulas
uncertainty propagation
cuts
canonical t/phi binning
efficiencies
acceptance
L/T separation
cross sections
```

Method A remains detached/non-production.

Method B remains diagnostic/cross-check only and numerically excluded.

No Method-A production promotion is permitted.

No empirical Fit 1/Fit 2 may be reintroduced.

---

## 9. Runtime-isolation requirement

The ordinary full-analysis launcher mutates repository-local runtime state.
Therefore the new owner must **not** run the analysis in the ordinary checkout.

Reuse the already tracked/reviewed isolation machinery from:

```text
testing/run_e8_4_fix5_canonical_five_plot_gate.py
```

where applicable, including its established concepts/functions for:

```text
ordinary-checkout snapshot/preservation
owner-created detached analysis/source worktrees
copied ltsep runtime overlay
ltsep import identity
path probes
external symlink preflight
installed-ltsep preservation
bounded owner-created cleanup
detached collector execution
```

Do not invent a second incompatible isolation architecture.

The new owner may import/reuse those helpers, but it must **not** execute or
inherit:

```text
fresh-F1 materialization
F.1/F.2/F.3/F.4 candidate staging
F.6.3 lineage reconstruction
E.8.4 page acceptance
canonical-five setting acceptance
```

If the existing isolation helpers cannot support the narrow debug run without
pulling those scientific dependencies into the E.8.2 gate, hard stop locally
and report the concrete blocker rather than weakening isolation.

---

## 10. Analysis execution owned by the new gate

The tracked owner shall own exactly one analysis command:

```text
./run_Prod_Analysis.sh -d 4p4 2p74
```

executed inside the owner-created isolated analysis worktree with the validated
copied-ltsep overlay/environment.

The owner must verify the existing debug completion markers, including:

```text
Left / lowe debug analysis completed.
Full high-epsilon processing is intentionally skipped in -d debug mode.
```

No second analysis invocation is allowed.

No canonical-five full run is allowed.

No direct `main.py` invocation is allowed.

---

## 11. New validation profile

Create:

```text
testing/pion_hgcer_validation_bundle_profile_e8_2_left_lowe.json
```

Schema:

```text
pion_hgcer_validation_bundle_profile/v4
```

Validation profile ID:

```text
phase_e8_2_left_lowe_baseline_scientific_audit/v1
```

Collection mode:

```text
generic_artifacts
```

### 11.1 Required declared profile settings

The profile's declared `settings` array must satisfy the frozen collector's
existing v4 invariant **exactly** and in this exact order:

```text
Left / lowe
Left / highe
Center / lowe
Center / highe
Right / highe
```

This is the profile's **authorization inventory**, not the scope of this farm
run.

Do not add Right/lowe.

Do not reduce the profile declaration to Left/lowe.

Do not edit `testing/collect_pion_hgcer_validation_bundle.py` to permit a
one-setting profile.

The owner must explicitly verify that the loaded profile has the canonical-five
declaration above and that:

```python
collector.resolve_settings(
    "Left", "lowe", profile
) == (("Left", "lowe"),)
```

### 11.2 Required requested collection scope

The owner must invoke the unchanged collector with the explicit setting filter:

```python
collect_validation_bundle(
    ...,
    phi="Left",
    epsilon="lowe",
    profile_path=<resolved E.8.2 profile>,
    ...
)
```

The resulting bundle manifest must therefore report exactly:

```text
requested_settings:
  Left / lowe
```

and its setting-artifact inventory must contain exactly one setting directory:

```text
Left_lowe
```

No highe, Center, or Right setting artifact may be packaged by this owner.

This distinction is mandatory:

```text
profile declares/authorizes canonical five
owner requests Left/lowe only
collector packages Left/lowe only
```

### 11.3 Artifact declaration

The profile must declare exactly these setting artifacts:

```text
procedure_pdf
page_manifest
full_analysis
correction_ledger_json
correction_ledger_csv
```

and one global artifact:

```text
run_summary
```

Use the existing established basenames:

```text
{phi}_kaon_rand_sub_{kinematic}_{epsilon}_full-background-subtraction.pdf
{phi}_kaon_rand_sub_{kinematic}_{epsilon}_full-background-subtraction-manifest.json
kaon_FullAnalysis_{kinematic}_{epsilon}.json
kaon_FullAnalysis_{kinematic}_{epsilon}_correction_ledger_no_empirical_residual.json
kaon_FullAnalysis_{kinematic}_{epsilon}_correction_ledger_no_empirical_residual.csv
```

These setting templates must remain valid for all five profile-declared
settings even though this owner requests only Left/lowe.

The run-summary basename may be newly defined, but must be unique and explicit,
for example:

```text
Q4p4W2p74_e8_2_left_lowe_scientific-audit-run-summary.json
```

### 11.4 Source identity

The tracked profile template shall pin the starting HEAD in
`required_analysis_commit`; the runtime owner shall resolve it to the exact
reviewed/pushed `--source-commit` using the existing fail-closed profile pattern.

`allowed_committed_files` must remain empty for the effective pushed source.

`allowed_non_analysis_path_prefixes` remains:

```text
docs/memory/
```

The owner/profile tests must prove that a temporary profile containing only
Left/lowe is rejected by the unchanged collector with
`validation_bundle_profile_invalid`, while the canonical-five declaration plus
the explicit Left/lowe filter resolves to one requested setting.

---

## 12. E.8.2 page verification

For each canonical parent `t1`, `t2`, and `t3`, require exactly one of each:

```text
full_background.e8_2.random.tN
full_background.e8_2.dummy.tN
full_background.e8_2.proton.tN
full_background.e8_2.pion.tN
full_background.e8_2.final_mm.tN
full_background.e8_2.stage_yields.tN
```

Total required E.8.2 pages:

```text
18
```

For every required page verify:

```text
scope == tN
t_index == N-1
represented phi inventory == phi_index 0..8 in order
invalid_unavailable_children == []
authoritative == false
presentation_only == true
```

Verify semantic-stage identity exactly:

```text
random       -> random_subtraction
dummy        -> dummy_subtraction
proton       -> slow_proton_cleaning
pion         -> baseline_pion_subtraction
final_mm     -> final_baseline_mm
stage_yields -> baseline_stage_yields
```

Require the terminal:

```text
full_background.e8.handoff
```

to remain after the E.8.2 pages.

Require the page-manifest setting identity:

```text
kinematic_token = Q4p4W2p74
epsilon_setting = low
epsilon_filename_token = lowe
phi_setting = Left
particle_type = kaon
```

Require:

```text
renderer_failures == []
```

An E.8.3 or E.8.4 unavailable page is **not** an E.8.2 failure by itself and
must neither be interpreted nor promoted by this owner.

Do not require any E.8.3 page.

Do not require any E.8.4 page.

Do not require Method-A F.1/F.3/F.4 artifacts.

---

## 13. Baseline analysis-identity verification

The new owner must verify the refreshed full-analysis JSON is current and
identifies:

```text
ParticleType = kaon
EPSSET = low
Q2 = 4p4
W = 2p74
OutFilename = FullAnalysis_Q4p4W2p74_lowe
```

Verify the debug analysis contains the intended Left low-epsilon setting and no
unexpected high-epsilon execution.

The correction-ledger JSON must verify:

```text
active_profile = no_empirical_residual
particle_type = kaon
epsset = low
q2 = 4p4
w = 2p74
```

The ledger CSV must contain the expected Left setting-total row.

Do not add a new interpretation of the correction ledger.

---

## 14. Freshness and provenance gates

Before analysis:

- verify exact branch/HEAD/local `origin/test`;
- verify source commit is a full 40-character SHA;
- verify no unreviewed source state;
- run collector source preflight from the exact detached source commit;
- record ordinary-checkout state;
- record installed-ltsep identities;
- establish isolated worktree/path agreement.

For every collected runtime artifact require:

```text
exists
nonzero size
fresh for this attempt
correct expected basename
valid JSON where applicable
```

The procedure PDF must start with `%PDF-`.

After the analysis child, including failure paths where applicable, verify:

```text
ordinary checkout unchanged
installed ltsep unchanged
owner-created worktree cleanup bounded to owner-created paths
```

Do not restore unexpected mutations. Fail instead.

---

## 15. Run summary and gate status

Create one atomic per-attempt gate-status record distinct from the ZIP.

The status schema may be new but must record at least:

```text
source_commit
kinematic = Q4p4W2p74
phi = Left
epsilon = lowe
status
stage
failure_reason
analysis_started
analysis_completed
artifact_verification_completed
collector_source_preflight_completed
collection_completed
zip_verification_completed
ordinary checkout preservation
installed ltsep preservation
```

The run summary must record:

```text
source_commit
analysis command
analysis return code
log hash
page count
18-page E.8.2 verification result
artifact hashes/sizes
isolation/preservation results
active profile = no_empirical_residual
```

and explicit boundaries:

```text
E.8.2 baseline audit only
Method A not required for gate acceptance
Method B numerically excluded
production promotion false
E.8.3 acceptance false
E.8.4 acceptance false
canonical-five acceptance false
absolute-SIMC amplitude claim false
```

The owner must not summarize E.8.4 failure as an E.8.2 scientific result.

---

## 16. ZIP verification and delivery

Use the unchanged generic collector through the new reviewed profile.

The profile loaded by that collector must still declare all five required
canonical settings. The owner must pass `phi="Left"` and `epsilon="lowe"`
to the collector so only Left/lowe is selected and packaged. Do not create a
temporary one-setting profile.

After collection, verify the ZIP:

```text
zip integrity passes
manifest complete == true
manifest errors == []
git_head == source_commit
required_analysis_commit == source_commit
profile declaration inventory == canonical five in frozen collector order
requested setting inventory == Left/lowe only
settings inventory == one Left/lowe record only
validation_profile == phase_e8_2_left_lowe_baseline_scientific_audit/v1
archive artifact inventory == one global run summary plus the five Left/lowe setting artifacts only
archive hashes/sizes match pre-collection records
```

On successful completion, return/deliver the ZIP plus the minimal companion
evidence needed for diagnosis/provenance:

```text
analysis log
run summary
gate-status receipt
```

Use the established configured farm/transfer roles. Do not hardcode a new
personal path or a guessed environment path.

Failure must not be reported as success merely because some artifacts exist.

---

## 17. Deterministic tests

Create:

```text
testing/test_run_e8_2_left_lowe_scientific_audit_gate.py
```

Cover at minimum:

1. exact profile schema/ID and **canonical-five declared setting inventory**;
2. unchanged collector rejects a one-setting v4 profile with
   `validation_bundle_profile_invalid`;
3. `collector.resolve_settings("Left", "lowe", profile)` returns only
   `(("Left", "lowe"),)`;
4. owner collection calls the unchanged collector with explicit
   `phi="Left", epsilon="lowe"`;
5. resulting collector manifest has exactly one requested/packaged setting,
   Left/lowe, despite the five-setting declaration;
6. exact E.8.2 18-page inventory;
7. page-manifest setting identity rejection;
8. duplicate/missing E.8.2 page rejection;
9. wrong t/scope/phi inventory rejection;
10. nonempty invalid-child inventory rejection;
11. renderer-failure rejection;
12. E.8.4 unavailable page does not fail E.8.2 acceptance;
13. Method-A candidate artifacts are not required;
14. source-commit/profile resolution is fail closed;
15. stale/missing artifact rejection;
16. wrong full-analysis identity rejection;
17. wrong active-profile/ledger identity rejection;
18. ZIP inventory/hash/source-identity verification, including absence of
    non-Left/lowe setting artifacts;
19. owner uses isolated worktree/runtime-overlay path rather than ordinary
    checkout execution;
20. ordinary checkout/installed-ltsep preservation failures fail the gate;
21. no second analysis command;
22. no canonical-five/F.6.3 materialization or lineage dependency.

Tests must not execute the farm.

Where practical, use pure-Python fixtures and mocks around the already-reviewed
isolation helpers.

## 18. Required local checks

Discover `<PYTHON>` using `docs/memory/TOOLS.md`.

Run at minimum:

```text
<PYTHON> -B -m py_compile \
  testing/run_e8_2_left_lowe_scientific_audit_gate.py \
  testing/test_run_e8_2_left_lowe_scientific_audit_gate.py

<PYTHON> -B -m unittest \
  testing.test_run_e8_2_left_lowe_scientific_audit_gate -v

<PYTHON> -B -m unittest \
  testing.test_e8_2_baseline_stage_audit -v
```

Also run any existing collector/profile regression that is directly needed to
show the generic collector was not changed or reinterpreted.

No ROOT/PyROOT or farm execution is permitted in local validation.

---

## 19. Memory update

Create:

```text
docs/memory/phases/e8-2-left-lowe-runtime-scientific-audit-owner.md
```

It must distinguish:

```text
SOURCE VERIFIED:
  current E.8.2 science/source audit and validation-owner mechanics.

RUNTIME VERIFIED:
  none from this implementation task.

NOT VERIFIED:
  actual Left/lowe stage magnitudes, pruning impact, PDF legibility and farm
  integration until the later user-run farm gate.
```

Update `docs/memory/CURRENT.md` concisely.

Preserve:

```text
E.8 ACTIVE
canonical-five full runtime BLOCKED
canonical-five provenance/identity repair DEFERRED
final E.8 BLOCKED
F.6.4 BLOCKED
absolute-SIMC blocker unchanged
Method A detached/non-production
Method B diagnostic/cross-check only and numerically excluded
```

Record the source-audit conclusion:

```text
no new E.8.2 scientific implementation is justified;
actual stage magnitudes require a narrow runtime gate.
```

The sole ordinary NEXT after a successful local implementation must remain
push-stable and substantive, equivalent to:

```text
NEXT — after ChatGPT actual-diff review, user commit/push and pushed-state
synchronization, run one isolated Q4p4W2p74 Left/lowe E.8.2 scientific-audit
farm gate and return the resulting ZIP/companions for scientific review.
```

Do not make commit/push itself the sole NEXT.

`CURRENT.md` is already near its soft 8-KiB threshold. Consolidate obsolete
wording rather than appending. Keep it below the soft limit if practical and
never exceed the hard limit.

Regenerate the memory manifest only after final intended contents.

---

## 20. Memory health

After all changes:

```text
<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/check_memory_health.py --root .
<PYTHON> -B tools/memory_bootstrap.py --root . --json
git diff --check
```

Do not use `--fail-on-warning` unless an ordinary warning is materially blocking
under `MAINTENANCE.md`.

Report:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
warning classification: <none or each warning: blocking | nonblocking>
manifest check: PASS | FAIL
bootstrap: PASS | FAIL
git diff --check: PASS | FAIL
```

---

## 21. Negative/regression checks

The final diff must contain **no** change to:

```text
src/
tools/
farm_env/
background_samples/
run_Prod_Analysis.sh
set_SymLinks.sh
testing/collect_pion_hgcer_validation_bundle.py
existing E.8.4 owners/profiles
testing/test_e8_2_baseline_stage_audit.py
accepted evidence
MEMORY.md
roadmap/STATUS.md
```

Search the new owner/profile/record for forbidden implications equivalent to:

```text
Method A is production
Method B contributes numerically
E.8.4 must pass for E.8.2
canonical-five runtime is accepted
canonical-five provenance repair is complete
absolute SIMC is resolved
Fit 1/Fit 2 are active
E.8.2 is runtime validated
farm run completed
```

None may be present as current claims.

---

## 22. Review bundle

Before stopping, run:

```text
git status --short --untracked-files=all
git diff --stat
git diff --check
```

Create one complete temporary root review bundle:

```text
kaonlt_review.diff
```

It must include:

- complete tracked diff for every modified tracked path;
- complete `git diff --no-index /dev/null ...` additions for every intended new
  file;
- no unrelated user state;
- no staging merely for review.

At minimum the complete additions must include:

```text
testing/pion_hgcer_validation_bundle_profile_e8_2_left_lowe.json
testing/run_e8_2_left_lowe_scientific_audit_gate.py
testing/test_run_e8_2_left_lowe_scientific_audit_gate.py
docs/memory/phases/e8-2-left-lowe-runtime-scientific-audit-owner-task-contract.md
docs/memory/phases/e8-2-left-lowe-runtime-scientific-audit-owner.md
```

---

## 23. Farm boundary

Codex must not run the farm.

Codex must not run ROOT/PyROOT full integration.

Codex must not issue a user farm command.

This task establishes source readiness only.

ChatGPT must review the actual diff, then the user controls commit/push, then
ChatGPT performs pushed-state synchronization and a full farm-readiness audit.

Only after that may ChatGPT provide one exact user-run farm command.

---

## 24. Acceptance criteria

PASS only if all are true:

1. exact starting branch/HEAD/origin/worktree gate passes;
2. only Section-7 paths change;
3. no production/scientific source changes;
4. the owner is isolated from the ordinary checkout;
5. established copied-ltsep/path/preservation machinery is reused or matched
   without importing F.6.3/E.8.4 scientific dependencies;
6. one and only one debug analysis command is owned;
7. profile declaration is exactly the collector-required canonical five, while the owner-requested collection scope is exactly Left/lowe;
8. package inventory is minimal and exact for the selected Left/lowe request;
9. all 18 E.8.2 pages are verified fail-closed;
10. all nine phi children are required on each E.8.2 parent page;
11. invalid/unavailable E.8.2 children fail the gate;
12. E.8.3/E.8.4 acceptance is not required;
13. Method-A candidates/materialization are not required;
14. no empirical residual is admitted;
15. full-analysis and ledger identity checks are explicit;
16. stale artifacts fail;
17. ZIP source/inventory/hash verification passes deterministic tests;
18. preservation failure fails the owner;
19. focused deterministic tests pass;
20. E.8.2 regression suite passes;
21. memory manifest/health/bootstrap pass;
22. `git diff --check` passes;
23. complete `kaonlt_review.diff` is prepared;
24. no commit, push, ref update, staging-for-review, or farm execution occurs.

---

## 25. Hard stop

Hard stop and report `BLOCKED` rather than broadening scope if implementation
appears to require changing
`testing/collect_pion_hgcer_validation_bundle.py`, weakening its canonical-five
profile invariant, or creating a one-setting v4 profile.

After implementation, deterministic checks, memory update, and complete review
bundle creation:

**STOP.**

Do not commit.

Do not push.

Do not run the farm.

Return:

- exact starting/ending branch, HEAD, and local `origin/test`;
- exact changed/untracked paths;
- concise implementation summary;
- exact local test results;
- exact memory-health report;
- `git diff --check` result;
- path to `kaonlt_review.diff`;
- explicit confirmation that no farm command was executed.

ChatGPT must review the actual diff before any user-controlled commit/push.
