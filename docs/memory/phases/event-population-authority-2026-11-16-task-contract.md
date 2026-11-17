# D2.a — event-level authority across physical distribution construction

## Authority and starting state

User's November 16 emergency D2.a request supersedes the preceding pion-only
recommendation. One coherent correction; no chain of replacement contracts.
Live GitHub/local test base: bd811bb3ce666afa58e654f0ded8cdb8550104c6.
Worktree is dirty: preserve the separately reviewed, default-off D1.a candidate,
its fixture/tests and unrelated work. Initial complete files/index/refs/diff and
ordered startup/health receipts are in the implementation-review bundle.
This contract is established before source editing. User has expressly authorized
implementation after this audit/contract; no additional permission gate is inferred.

## Confirmed failures and ownership

1. `src/utility/pion_component_shapes.py:load_pion_simc_component_shape`
   fills setting/native-phi populations directly but retains no direct t parent.
   Each child is independently unit-normalized. Cache/clone payloads omit parent
   population support and do not fully distinguish physical selection semantics.
2. `src/cuts/pion_component_fits.py:_resolve_component_scope_hist` and
   `_resolve_single_scope_hist` reconstruct t parents from normalized children;
   sparse children silently inherit another t/setting population. The shared
   `build_particle_subtraction_component_result` entry must reject missing scoped
   pion populations before fitting, including alternate non-parent callers.
   This affects
   pi_n/pi_delta/pi_sidis and KLambda/KSigma0, not only pion diagnostics.
3. `src/cuts/rand_sub.py:rand_sub` passes initial setting-only SIMC payloads to
   canonical parent construction. The canonical edges are resolved after main's
   initial setting-template load. Scoped payloads must be acquired from original
   SIMC events at the existing frozen-edge boundary, before parent fitting.
4. `src/cuts/pion_t_bin_parents.py:build_setting_t_bin_pion_parents` reuse must
   distinguish the population-construction/source identity, not just experimental
   inputs. Preserve experimental parent ownership and D1.a metadata handoff.
5. `src/utility/utility.py:prune_hist` arbitrarily zeros signed spectra/high-error
   bins or resets nonpositive/low-entry histograms used by yield/average paths.
   Make its existing interface observational: report candidate count, preserve
   physical contents, errors, entries and flows. Do not replace the adopted
   four-criterion Lambda mask or introduce another exclusion rule.
6. `src/utility/utility.py:TH1D_to_TH2D` and the artificial-joint section of
   `src/plotting/data_vs_simc.py:plot_data_vs_simc` multiply unrelated marginals and
   present the result as an accepted joint detector population. Require original
   selected paired (x,y,weight) records for joint filling. Never synthesize a
   joint distribution from marginals. Where the normal display payload has no
   such paired population, render an explicit unavailable notice; do not crash
   the analysis or silently substitute independent marginals as measurements.
7. `src/cuts/pion_t_phi_template_comparison.py:scope_evidence` must report strict
   unavailable/direct-event scope semantics; remove its fallback interpretation.

## Allowed files and functions

Only the seven source files named above (loader/cache/clone/empty/map helpers;
two scope resolvers, a shared strict event-scope helper and a scoped-input guard
at the shared component-result entry; rand_sub frozen-edge
SIMC handoff; parent reuse/record identity; pruning/joint utility functions;
plot joint section; detached scope diagnostics), new focused tests
`testing/test_event_population_authority.py`, and the existing directly affected
resolver expectations in `testing/test_pion_t_phi_template_comparison.py`.
Memory: CURRENT.md, this standalone contract at
`docs/memory/phases/event-population-authority-2026-11-16-task-contract.md`, one
phase record `docs/memory/phases/event-population-authority-2026-11-16.md`, and
manifest.json only. No MEMORY.md growth or roadmap milestone edits.

`src/main.py`, apply_cuts/binning owners, SIMC iteration/normalization owners and
all original event-loop owners are audited and frozen unless a newly confirmed
violation requires an explicit revision of THIS contract before that edit.

## Intended event ownership and rules

Use the same original MC tree, iter_weight/default-one weight, existing cuts,
HGCer-hole selection, -t and phi coordinate, half-open native bin conventions,
MM underflow/overflow and setting normalizations. Fill setting, t parent and
native (t,phi) directly during the original accepted MC event loop. A parent
includes its declared t population, with outside-phi events accounted explicitly;
no child is used to manufacture it. Retain raw count/signed/absolute/Sumw2 support
and provenance. Sumw2 exists before filling. Normalize only after accumulation.

Unscoped requests retain valid setting-wide purpose. Scoped requests select only
their event-originated parent/child. Missing, empty or nonpositive-normal-integral
support is explicitly unavailable, never replaced by another accepted population.
Fit/model shifts, model amplitudes, residuals and plot clones remain derived model
operations; they are not new accepted event records. No new fit estimator, physical
response calibration, amplitude floor, mask or model rescaling is introduced.

Cache/clone identity includes exact edges, source file identity, selection keys,
weight/normalization/population schema and source provenance. Raw support survives
cloning. Existing authoritative experimental parents are not reread from averages
or original experimental trees. Their frozen producer ownership stays unchanged.

## Compatibility and tests

Exercise actual loader and resolver functions with detached ROOT adapters. Cover
all five template names, unequal phi occupancy, direct parent/child sumw/Sumw2
closure, outside-phi and final boundary events, MM flows, signed cancellation,
zero/nonfinite weights, missing sources/branches/provenance, no setting substitution,
partial scopes, invalid edges, cache key/clone isolation and parent reuse guard.
Verify unscoped setting values unchanged and no normalization-factor invention.

Use all 15 retained pion source row members for 45 parent/405 channel-cell closures
(135 physical cells), preserving source hashes and cohort labels. These rows are
already accepted: they validate population arithmetic, not a fresh cut/farm replay.
KLambda/KSigma0 use deterministic implementation-path tests; genuine source-bound
event rows for these signals are not available in the retained D2 projection.
Do not fabricate them or claim empirical all-setting signal closure. Their source
and ROOT integration remain a specific outstanding farm-validation requirement.

Test prune_hist with observed signed cancellation, negative/zero integrals,
high fractional errors and low entry counts: no accepted content/error changes.
Test joint filling from retained pairs; reject absent pairs and prove the former
outer product cannot run. Render missing joint input explicitly. Test D1.a original
source-support/empty observations/cancellations/default-off/unchanged bins using
its existing focused suite. Broader relevant tests must distinguish inherited
failures from new regressions. No ROOT/main.py execution locally.

## Evidence, memory and stop

Generate source-bound before/after native-MM pion figures as individual high-
resolution PNG and vector SVG files, direct-event versus former-binned closure,
all-parent tables and detailed Center/high t3 existing-fit conditional comparison.
Existing saved raw-pion B fits freeze historical signal inputs/alignment: they
are NOT full newly executed D2.a fits, since signal-scope behavior also changes.
Report that distinction and unavailable corrected-fit quantities explicitly.

Update CURRENT narrowly to ACTIVE D2.a and exact next independent actual-diff
review; preserve D1/D1.a/D2 history, D3-D5 DEFERRED and Outcome B BLOCKED. No closed
runtime phase or six-step roadmap changes. Preserve Method A detached and Method B
diagnostic-only, random/dummy/proton subtraction, HGCer holes, binning, original
SIMC event normalization, cross-section iteration and L/T machinery.

Deliver one complete implementation-review ZIP with actual full tracked diff,
separate D2.a-only delta against the captured dirty baseline, all new files,
before/after source, operation classification, contract, tests/logs, numeric arrays,
individual PNG/SVG figures, provenance, SHA256 manifest/receipt, initial/final full
health and exact farm requirements. Validate committed-LF memory in an isolated
snapshot without changing the preserved CRLF D1.a contract checkout.

Hard stops: ambiguous physical selection/weights that would require invention;
unrelated producer/physics changes; nonempty/mutated index; new unresolved regression;
memory hard failure; missing provenance promoted to evidence. Do not hide unavailable
signal rows, ROOT validation or complete corrected fits. If a hard stop remains,
label the source gate BLOCKED and deliver only a clearly scoped diagnostic packet.
No stage, commit, push, hosted CI, farm, production-wide rerun or publication.
