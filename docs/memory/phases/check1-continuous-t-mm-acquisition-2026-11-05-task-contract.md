# Check 1 continuous-t/MM acquisition — Codex task contract

## Authority, starting state and hard boundary

Audited starting `test` HEAD: **6663874fae232729093a7dd51daffd6c710621a8**.
Task class: narrowly scoped additive diagnostic export and extension of the
existing isolated acquisition owner. This contract is proposed implementation
authority after independent contract review, not farm readiness or execution
authority. Work has implemented no source changes.

Read tracked root AGENTS, CURRENT, MEMORY, CURRENT_HANDOFF and USER in that
order, then MAINTENANCE, CODEX, the existing
`investigations/full-t-acceptance-edge-check1-2026-11-05.md`, the existing high
acquisition contract and the exact source owners below. Verify live HEAD, refs,
index and all tracked/untracked pending state before editing. If HEAD moved,
reconcile only the intervening relevant source diff; stop for a scientific or
ownership change, never silently repin. Preserve all pending work. In particular,
do not absorb the untracked offline-testing implementation files or their
unrelated contract into this implementation.

The prior Check 1 result, its two input dates and its five original figures are
completed evidence. Its acceptance-edge conclusion is insufficient. Do not rerun
that offline investigation, invent sparse thresholds or claim an acceptance
effect. Six roadmap statuses, Outcome B, F.1-F.4/F.6 scopes and Method A/B roles
are frozen. The Lambda/Sigma0 tmin observation is recorded context only.

## Decision and exact evidence-producing outcome

**Path B applies.** The selected detached Local-t observer projects away
continuous t; the SIMC average cache lacks MM and is allcuts-only; the normal
support serializer writes no joint t-MM. The existing owner also collects only
high artifacts in high mode. No existing controlled artifact supplies the whole
requested paired evidence. This is a source-backed persistence/collection gap,
not a claim of global farm-file absence. Path A is not established by the current
saved-object evidence. Do not leave that architecture choice to implementation.

Add one default-off Check 1 export, using the existing selected Local-t
application and existing SIMC loop, and one new mode of the **existing owner**.
One later owner-controlled attempt processes paired canonical preflight and
both epsilons once in its established private environment, validates all five
settings, then returns the complete numerical evidence before private cleanup.
No SIMC generation, new fitting procedure, production promotion or new pipeline.
The necessary existing analysis path may execute its existing fits to reproduce
the pinned cohort; the exporter itself calls no fitter and never changes a fit
or weight. Failure to reproduce the reference stops the attempt.

## Exact implementation allowlist

Existing source files, bounded edits only:

- `src/cuts/pion_t_phi_template_comparison.py`: guarded data-export call after
  its unchanged successful observation/serialization; copy scalar cache fields
  without changing legacy rows, fingerprints, fits, transfer vectors or payloads.
- `src/binning/calculate_yield.py`: optional SIMC diagnostic sink in
  `_process_yield_simc_tree`, its construction/finalization in `process_hist_simc`,
  and no other scientific algorithm edits. Keep the legacy average-cache schema.
- `testing/run_pion_t_template_comparison_gate.py`: explicit Check 1 mode,
  reference/preflight binding, combined verification/profile/collection and
  failure-safe publication; retain the existing execution/isolation machinery.
- `testing/test_pion_t_phi_template_comparison.py` and
  `testing/test_run_pion_t_template_comparison_gate.py`: relevant additions and
  regressions; do not replace existing fixtures or loosen their assertions.

New files, first confirm absent:

- `src/utility/check1_acceptance_evidence.py`: narrowly scoped ROOT-independent
  scalar/array packet builder, validator and atomic sidecar writer. No fitter,
  ROOT traversal, generic artifact framework or replacement yield producer.
- `testing/test_check1_acceptance_evidence.py`: deterministic packet and real
  seam tests, including actual writer/JSON/NPZ round trips.
- `testing/pion_check1_acceptance_evidence_profile_v1.json`: existing generic
  collector schema v4; exact combined artifact declaration described below.

Allowed memory: CURRENT, this contract, the existing Check 1 investigation and
regenerated manifest.json. Update only warranted implementation/status facts.
Everything else is frozen, especially `main.py`, `xsect_support.py`, launcher,
cuts/coordinates/canonical helpers, normalization, fitters/templates/config,
proton application, both existing profiles, the local-t event-validation capture
contract, generic collector, installed dependency/isolation helpers and their
reviewed blob pins. Stop if another file or scientific change is required.

## Source seams and before/after behavior

1. `pion_t_phi_template_comparison.observe` already has `rows`, projections,
   `result.cells[].applications.Local_t` and live experimental caches. Its
   Local-t means fresh B raw-weighted parent-t templates with fixed A alignment,
   not production A or the child-fitted Local_t_phi. After the normal successful
   result is serialized, call the new exporter **only** when
   `LT_ANALYSIS_CHECK1_EVIDENCE=1`. It copies continuous fields by the existing
   population/source/tree-entry identities and uses each actual saved Local-t
   transfer vector. No second model evaluation/refit, Global substitution or
   independent event selection. Failure propagates to the analysis/owner.
2. Existing SIMC `_process_yield_simc_tree` iterates h10 and discards non-allcuts
   events; `_init_ave_simc_event_cache` has no missmass. Add an optional sink,
   default None. In the **same** event pass, observe `missmass`, `-evt.t`, phipq,
   entry index, iter_weight and the actual allcuts flag. When the sink exists,
   also evaluate existing `apply_simc_sub_cuts` with the same hole decision for
   the existing no-MM population. Do not redefine that predicate, change any
   existing continue/fill or add another ROOT event pass. Disabled execution
   must perform exactly its old calls/fills/cache writes.
3. Finalize the SIMC packet in `process_hist_simc` after its existing histogram
   scaling, using those actual histograms as closure references and recording
   `normfac_simc` from the existing owner. The legacy `iter_weight` is unscaled;
   diagnostic event yield is iter_weight times normfac_simc **once**, as verified
   against the already-scaled producer. Never scale existing histograms again.
4. Export the whole currently accepted analysis interval, not just t3 or a
   new edge selection. Read its range from frozen configuration and its fine-t
   edges from the existing 100-bin t axis. Record canonical membership and
   preserve the full no-MM population with MM flows. Assert that canonical
   coverage spans the configured accepted interval for this frozen profile;
   if it does not, stop rather than silently truncating or extending weights
   outside their defined parents. Full here means the existing analysis-cut
   accepted interval, not generator/geometrical acceptance outside those cuts.
5. Opt-in guards must exclude prepasses/optimizer, wrong particle/kinematic,
   iteration >0 and ordinary modes. Disabled/prepass paths import no new
   diagnostic module, allocate no sink, hash no new file and emit no artifact.
   Do not change source coefficients, flags, proton factors, HGCer cuts or
   cleanup timing. Ordinary histograms, NPZs, yields and default owner modes
   keep their existing semantics.

## External reference authority and cohort gate

The Work-supplied `KaonLT_Check1_Reference_Inputs_20261105.zip` is external,
read-only reference evidence. SHA-256:
`02ff71998fd51d5e56f07484bb1c71185f78ed51d6ce409b0d34154bda7a6d7a`.
It contains `reference_manifest.json`, `current_lowe.zip`, `current_highe.zip`,
two supplied original run summaries and `SIMC_REFERENCE_lowe/highe.npz`.
Only simulation fields and axes were copied into these NPZs; no historical
experimental NPZ fields or additional normalization were introduced.

New CLI mode: `--check1-evidence --check1-reference-bundle <external-file>` on
the existing owner; Check 1 implies its already established high continuation.
Retain `--source-commit`, `--repo`, `--outdir`, `--output` and preflight-only.
Reject a missing reference, conflicting mode or injected flags. Do not execute
these operations during Codex implementation.

Validate the outer hash, safe inventory, member bytes/hashes, cohort/schema and
all five identities before scientific copies or child execution. Bind each
ordinary input through the supplied original run-summary identities to the
correct low/high role; retain distinct pre-copy and post-derived-branch runtime
identities. Do not compare a known private mutated container's hash to the
ordinary input as if they were the same role. No date-only SIMC restriction.
No silent permission to substitute another experimental cohort. If original
input identities changed, stop before analysis; do not rebuild inputs.

The child authority includes the reference digest, expected settings, source
SHA, fixed configuration/coordinate/PID identities and private output role.
Pass it through owner-controlled, hash-pinned Check 1 environment variables.
Extend child-environment whitelists only in this mode; refuse user-injected
unknown LT_ANALYSIS variables. Child metadata maps runtime files to the owner's
actual copied/consumed-file identities. Reuse already-computed hashes for large
inputs; do not hash whole output trees or skim contents.

The acquisition must reproduce the reference experimental observations,
Local-t transfer vectors and applications for every canonical cell, as well as
the unchanged SIMC native Lambda-window values/errors. Compare typed identities,
counts and axes exactly; use only existing owner's numerical closure tolerance
`rel_tol=1e-8, abs_tol=1e-12` for floating accumulation. Compare a closed scientific
projection, excluding legitimate new source SHA/runtime path fields, never an
arbitrary subset chosen after failure. No threshold relaxation, input mutation,
fit tuning or automatic retry to satisfy these gates.

## Output schema and physical ownership

For each of the five settings, write four unique basenames:
`<Setting>_kaon_Q4p4W2p74_<lowe|highe>_check1_<data|simc>.json` and `.npz`.
Metadata schema `kaonlt_check1_evidence/v1`; load arrays with allow_pickle=False.
The JSON names its NPZ, byte length/hash, column dtypes, counts, selection stages,
axes, source-file role binding, source/inventory hashes and validation status.
Use deterministic setting/epsilon, source, tree-entry, parent/child order.
Reject duplicate physical identities within a population/source/tree; kaon and
pion PID populations and separate trees remain distinct identity domains.

Required metadata:

- Q2/W, setting, epsilon, MC iteration zero, branch `detached_Local_t`,
  diagnostic_only=true, production_applied=false, scientifically_accepted=false.
- Exact source SHA, all executable exporter/observer/producer/owner/profile
  identities, reference digest, original and runtime file roles, tree names,
  configuration and coordinate/PID/cut fingerprints, parent/fit/transfer-vector
  fingerprints. Do not collapse reviewed, observed and farm-executed identities.
- Positive experimental analysis -t = -MandelT + existing kaon t_shift, MM with
  existing kaon MM shift; SIMC -t=-evt.t and MM=evt.missmass. Units GeV^2 and GeV;
  phi in degrees. Record actual shifts/conventions, never apply an extra shift.
- Source coefficients exactly as used: data normfac, -data normfac/nWindows,
  -dummy normfac, +dummy normfac/nWindows; preserve their normalization state.
  Record raw SIMC iter_weight, existing normfac and final once-scaled weight.
- Exact canonical t/phi edges and full fine-t histogram edges; MM display axis
  from the existing parent/control histogram (current reference 0.7-1.3 GeV,
  6-MeV bins), plus exact [1.10,1.16) allcuts and the native fine SIMC MM axis.
  Do not invent a new fit window or change production binning. Event rows beyond
  the broad display axis are retained in flows, not deleted.
- 15 parent/135 child slot inventory with explicit selected/empty/invalid/
  unavailable states and reasons, including separate kaon, pion and MC support;
  outside-phi contributions remain separate. No forthcoming exclusion mask.
- Full fit/covariance information is retained through references to the unchanged
  original parent/child JSON artifacts. Export conditional sumw2 honestly;
  complete fitted-background/training/shared-systematic covariance is explicitly
  unavailable where not actually retained. Do not label counting errors as its
  substitute. Preserve naturally emitted Y0/MM_0/validity if actually available
  with its real branch label; detached B is not automatically producer Y0/MM_0.

Experimental NPZ has separate kaon and pion row tables containing continuous
t/MM/phi, population/source/tree/entry identity, canonical t/phi indices, both
existing flags, source coefficient and zero/negative contributions. Kaon rows
retain proton factor and RF accept, existing pre-pion signed weight and weight^2.
Pion rows retain the exact Local-t transfer-bin index/vector fingerprint and
transfer weight; background contribution is existing source coefficient times
that weight and its squared contribution. After-subtraction contributions are
kaon pre-pion rows and negative background rows; do not pair unlike physical
events or recompute a kaon/pion classifier.

SIMC rows retain original h10 entry, t/MM/phi, canonical membership, both flags,
raw iter_weight, normfac, once-normalized weight and squared weight. Do not
substitute unit templates; do not generate new events or change cross sections.

Also supply per-setting/source count, sumw and sumw2 1D fine-t and joint t-MM
arrays for existing allcuts/nommcuts. Joint axis order is explicitly
`[source, t_with_flows, MM_with_flows]` per selection and stage; MC has one
simulation source. Values are integrated bin weights, **not densities**.
Density plots later divide by actual widths only. Store counts separately from
signed yields, with all flow bins and the unchanged native axes. Canonical/phi
identities stay in row tables rather than enormous redundant child tensors.
Conditional after sumw2 is the actual sum of squared contributing kaon and pion
weights; it excludes the separately declared unavailable fitted covariance.

Honor each producer's actual membership convention: canonical helper includes
the final endpoint while yield-loop half-open range cuts can exclude it.
Histogram upper axis edges go to ROOT overflow. Never unify these by changing
cuts or silently migrating a row. Retain original indices/flags and validate
them against their declared owner; account for any endpoint difference.

## Independent numerical and serialization gates

The new pure helper validates rows, source bindings, dtypes, finiteness and
nonnegative squared weights/counts. Compute raw/support counts, signed and
absolute weight sums separately. Verify row-to-1D and row-to-2D closure;
sum along each joint axis including flows must equal its independently formed
1D projection and actual per-source sumw2. Verify parent-child-outside-phi and
source-sum closure for both existing selections.

Project new data rows onto the **old MM axis**, per all 135 slots, to reproduce
the observer kaon/pion inputs and exact Local-t applications, before/after
contents and fixed-weight sumw2. Compare data to matching **detached B**, never
force B to equal ordinary production A. Verify SIMC allcuts on its own native
axis against actual once-scaled histograms and the supplied simulation-only
reference. Fine t can also project existing matching simulation t-vs-tmin values,
but those 2D saved arrays do not supply lost sumw2 or t-MM coupling.

Available empty cells remain explicit zero-support states, not positive evidence.
Unavailable Local-t transfer with populated pion support is a hard acquisition
failure. If the observed pion population is exactly empty, background zero is
supported by that empty population; retain missing-fit status without inventing
a fitted transfer. Invalid/nonfinite inputs, wrong stage/source or missing
required quantities fail, never become empty/zero or trigger a fallback.

Create each NPZ/JSON pair under a unique private temporary name; validate and
hash before publishing its unique metadata last. Refuse existing/symlinked
destinations and leave no apparently complete pair after failure. Packet status
means captured/pending owner verification, not scientific acceptance. Do not
alter an original source file or ordinary production artifact.

## Complete existing-owner extension

Reuse `execute_gate` and its private code worktree, copied ltsep/PYTHONPATH
overlay, volatile private data placement, existence-only skim checks, bounded
copy/hash cost inventory, input/source preservation and cleanup. The only child
analysis is the already established launcher operation, with paired canonical
preflight and low then high continuation. Do not add an acquisition wrapper,
parallel owner, new farm job or manual shell sequence.

Existing old modes/profiles remain byte-semantic compatible. Check 1 mode uses
the new profile: all **15 existing low+high parent JSON/PDF/child JSON artifacts**
plus **20 new sidecars**, exactly 35 global artifacts and the five setting
identities. Keep the ordinary collector schema and logic frozen. Retain all
available fit/uncertainty information in those original artifacts. Expected
original sidecar absence/freshness is checked before the child.

The current high committed-path guard rejects later memory changes at this
HEAD; the new mode must explicitly permit `docs/memory/` as the existing generic
owner does, plus only the exact new/modified implementation paths above. This
must not become unrestricted descendant acceptance. Preserve ancestor checks,
strict scientific allowlists, reviewed cleanup/publication blob pins and all
installed caller/dependency/namespace/alias guards. Extend committed source and
sidecar inventories to cover the actual exporter/validator/profile; do not claim
the old observer source tuple covers a new executable dependency.

High mode currently checks/collects high artifacts only. In Check 1 mode verify
both epsilon parent/child payloads through the existing verification functions,
with an explicit all-five expected inventory; then verify all new sidecars and
their common reference/axes/cohort bindings. Refuse cross-attempt low artifacts,
stale low outputs, dropped low evidence or a high-only success. Extend status,
summary, source receipts and ZIP readback to the full five-setting operation.

Use the unchanged generic collector with the combined profile. Stage the ZIP
under a visibly pending, unique destination until its readback and the existing
ordinary-input/checkout/installed-package/external-link/cleanup preservation
checks pass. Publish to the requested absent destination without clobbering a
colliding file only in the success path; retain failed status outside private
cleanup and no accepted-looking final ZIP after failure. Preserve the established
companion gate-status/run-summary/log delivery, joining them to the released
ZIP hash. No circular self-hash claim is required inside the ZIP.

## Deterministic local verification — required before review

Use the existing unittest framework and actual pure helper/writer. Do not
implement tests that only restate expected dictionaries or fake successful
validators. Required cases:

1. Controlled event arrays with independently hand-computed t/MM projections,
   exact source identities, both flags and all five logical settings.
2. Positive, negative, zero and cancelling signed weights; variance remains
   quadratic, negative data is retained and counts differ from sumw.
3. Actual NPZ/JSON round trip, dtypes/units, deterministic row/schema order and
   rejection of pickle/object arrays, NaN/Inf, duplicate or colliding identities.
4. Every canonical boundary and nextafter neighbor, final endpoint, phi-outside,
   exact MM/axis upper endpoints, underflow/overflow and wrong/transposed axes.
5. Existing data/PID/proton/RF ownership and existing Local-t weights; zero and
   unavailable transfers; populated-invalid fails, empty remains distinct.
6. New data joint projections reproduce existing real `project` and
   `apply_transfer` outputs for both flags; use their actual serializer and
   include failure/rejected/empty child evidence.
7. Real `_process_yield_simc_tree` seam with controlled fake events and the
   injectable sink: capture MM before allcuts rejection, broader no-MM support,
   correct hole and indices, once-scaled weight, native histogram closure.
   Include a non-unit normfac and a deliberate double-normalization failure.
8. Sink disabled/prepass/wrong-mode must leave every old loop/call/cache/output
   unchanged. Export failure changes no production objects. No call to any
   fitter, set_val, proton evaluator, ROOT traversal or generator by exporter.
9. Reference original/runtime identity distinction; mismatched data, MC,
   iteration, branch, coordinates, templates or transfer fails before admission.
10. Full-five owner fixture using actual profile/validator/collector round trips:
    reject missing low/high setting, one missing sidecar, native-array mismatch,
    stale artifact, mixed attempt, tampered member, manifest/hash mismatch and
    unapproved scientific-source descendant; permit a memory-only descendant
    only through the explicit Check 1 memory-prefix and ancestor checks.
11. Atomic pair/final-ZIP collision and injected write/verification/preservation/
    cleanup failures: no successful final output, nonzero result and retained
    first-invariant diagnostic; never auto-retry.
12. Preserve existing parent/child, owner high/low/preflight, optional versus
    mandatory SIMC, private path/alias/caller/dependency-closure, old serialized
    missing-fit and cleanup/source pin regressions.

Local interpreter observed by Work: Windows Python 3.12.10 / NumPy 2.3.0 usable;
WSL default python3 exists but lacks NumPy. Do not install/reconfigure packages
or fake scientific dependencies under this contract. Windows Python 3.12.10 /
NumPy 2.3.0 may execute supported deterministic tests. Missing WSL NumPy is not
a scientific failure and does not block implementation. A Linux-specific test
that genuinely cannot execute in available workstation environments is
NOT VERIFIED: report its exact purpose and environment requirements. Such
unexecuted tests prohibit farm-readiness approval, not implementation. Farm
Python, ROOT, dependency, path/alias, executable-owner, provenance and manifest
checks remain mandatory through the reviewed farm preflight before acquisition.
Do not weaken assertions, fabricate passes, substitute inputs or bypass failures.
Compile/AST-check changed/new Python without ROOT and check Python 3.9 grammar;
grammar compatibility is not execution on farm Python 3.9.

Run new unittest suite, existing child/parent and owner suites, and affected
serializer/collector profile regressions. Log exact interpreter, versions,
counts/skips/reasons. No ROOT/PyROOT/main/farm run is authorized. Work did not
run those implementation suites while preparing this contract.

## Gated later farm acquisition plan — not an execution command

**Preflight-only:** after independent actual-diff acceptance, user publication
and pushed-state identity review, user invokes the reviewed owner with Check 1
and preflight-only. It must validate the exact pushed SHA, all-five source and
profile inventories, pinned external reference, original current-cohort inputs,
existing compatible MC inputs/.hist/model, canonical pair and frozen configuration,
actual farm interpreter/module/caller/namespace identities, paths/aliases/overlay,
destinations and bounded I/O inventory. No scientific copies, child or ZIP in
this mode. Reuse existing environment roles; no personal paths in tracked memory.
Report unavailable inputs or capabilities before a lengthy analysis.

**Execution:** only after that readiness evidence is reviewed, the user runs
one owner-controlled acquisition at the same SHA/reference. Existing launcher,
paired preflight, private low/high processing and opt-in capture only. No second
manual extraction or new SIMC generation. Do not hash existence-only skim files.
The analysis is necessary because the selected Local-t application/transient
joint coordinates are not retained by the established collector; this is not
permission to rebuild producers or to rerun unaffected separate investigations.

**Validation:** existing owner payload gates plus new per-event/joint/native-MM
and fine-t closure, all source labels/units/flows, 15 parents/135 slots and five
paired settings, unchanged reference science, preserved conditional/covariance
states. Full-t rows must have genuine non-canonical coordinates; no splitting
or replication of three totals. Presence/complete=true alone never suffices.

**Collection:** 35 declared numerical/metadata/original evidence artifacts in
the existing collector, source/profile/reference identities, checks/manifest
and hash receipt, plus established status/summary/log companions. Collect and
read back before cleanup; only release after preservation succeeds. No new farm
figures are required by this export contract. Existing legitimate PDFs may be
collected unchanged; later Work plots use these actual numerical rows.

**Failure:** first failed owner/source/identity/projection/artifact/manifest/
preservation gate stops everything downstream. Record first invariant, stage and
earliest valid repair; no automatic rerun, repackaging or reinterpretation as
scientific success. Failed artifacts stay explicitly pending/failed, not
accepted. No farm invocation is given at this contract-preparation gate.

## Memory, actual-diff review, acceptance and stop

Capture full pre-edit pending state and preserve frozen bytes/index/refs. Deliver
the complete actual implementation diff and full additions for every new file,
without staging. Trace default-off and enabled paths, actual row/weight ownership,
producer -> NPZ/JSON -> owner checker -> profile/collector -> ZIP -> cleanup ->
delivery. Inspect failure semantics and the old high committed-path guard.
Codex's summary is not independent evidence or readiness approval.

Update allowlisted memory only when warranted. Regenerate/check memory manifest
and run ordinary memory health; report CURRENT/MEMORY/CURRENT_HANDOFF bytes and
classify warnings. Initial soft warnings are nonblocking, scheduled at the next
checkpoint/milestone; do not expand MEMORY near its hard limit. Leave CURRENT's
single NEXT as independent actual-diff review followed conditionally by user
publication/pushed review and the reviewed preflight/acquisition gate, never
automatic farm execution or only 'push'.

Acceptance: exact allowlist only, complete actual diff, meaningful deterministic
checks with explicit unrun cases, source and disabled-path preservation, pinned
five-setting reference/sidecar/profile contract, a complete tracked invocation
chain and no unresolved implementation/scientific invariant. Result is at most
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING, not acceptance-edge physics or
ROOT/runtime validation. User alone commits/pushes and runs the farm.

Hard stop: moved scientific source; mismatched current-cohort/reference/MC;
unavailable actual Local-t application; need to change a fit, cut, normalization,
physics model, defined analysis range, uncertainty convention or a frozen helper;
missing complete orchestration; a substantive failed numerical, integrity,
source-preservation or executable required test; or a needed out-of-allowlist
file. Unexecuted environment-specific tests are NOT VERIFIED readiness blockers
as specified above. Report the precise blocker.
Do not implement an alternative architecture, relax gates, run farm analysis,
commit/push, execute Check 2/Method A, reopen Outcome B, study Sigma0 or begin
final yield/SIMC-iteration/LT work. Stop after local implementation/checks and
actual review material; the subsequent independent review is a separate gate.
