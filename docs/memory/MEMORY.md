# Durable KaonLT project knowledge

## Authority and evidence

`docs/memory/` is the repository-owned durable continuity layer. Native
assistant/Codex memory is supplemental only and cannot override current source,
repository memory, or supplied farm evidence. `CURRENT.md` is the sole
ordinary active-state authority. The handoff records exceptional transfer state
only and cannot override CURRENT.

At substantial-work start, establish the actual branch, HEAD, and worktree.
For implementation, current source and relevant diff outrank stale summaries;
for runtime acceptance, fresh farm/runtime evidence outranks source review.
Source review never establishes farm validation. Historical and chat material
is lower authority than canonical current, evidence, phase, and decision
records.

## Source and provenance semantics

Reviewed scientific source, repository observation, farm-evaluated source,
renderer source, bundle/profile source, and artifact provenance are distinct
roles. A source review may be anchored to an analysis source older than the
repository observation. A stored observed HEAD is a timestamped observation,
not a permanent statement of live identity.

Record a farm-evaluated identity only when direct applicable farm evidence
supports it. Documentation or memory-only work does not automatically
invalidate an older reviewed analysis source; reconcile the intervening diff
and re-review only when the relevant analysis source or test scope changed.
Never collapse provenance roles merely for convenience.

## Scientific ownership and production ordering

Preserve separate scientific owners for random subtraction, slow-proton PID
contamination, pion-production background, HGCer diagnostics, SIMC, yields,
and cross sections. The production ordering remains random/dummy -> frozen
binning -> slow proton -> pion subtraction. A setting-wide K-Lambda gate
controls whether proposed proton weights become applied; never partially
commit per-t proton results.

Legacy empirical residual Fit 1/Fit 2 are dormant historical machinery. The
accepted active profile is `no_empirical_residual`, which forces both empirical
residual-background scales to zero; those fits are not part of the current
production or presentation chain and must not be reintroduced by E.8, Method A,
diagnostics, or presentation work.

Pion-alignment comparisons use fixed evaluation envelopes. HGCer diagnostics
consume the frozen baseline rather than redefine it. Cuts, templates, priors,
component definitions, normalizations, binning, efficiencies, acceptance, L/T
separation, and uncertainty propagation remain frozen unless a narrow approved
contract explicitly owns a change.

## HGCer and Method-A/Method-B boundaries

Method A and Method B are independent. Method B uses frozen upstream records
and same-canonical-t relative closure; it uses no Method-A numbers, no cross-t
pooling or interpolation, and no absolute neutron normalization. Adaptive
Method B remains DO NOT PROMOTE. Method B is diagnostic/cross-check/historical
comparison only and never changes production pion weights.

Presentation cannot recompute a diagnostic or construct a correction. Method A
remains a candidate numerical HGCer path only after an explicit, validated F.6
production-promotion decision. Its positive-response training and physical
application populations remain distinct: training uses the established
prompt/noRF/nommcuts NPE>0 population, while physical application uses the
authoritative NPE>2 cache. Keep their provenance and fingerprints separate.
Future response application is event-level and parent-t normalized; never
renormalize `(t,phi)` children independently. Detached F-stage diagnostics do
not become production merely because they are accepted.

## Diagnostics and runtime validation

When applicable, trace persisted diagnostics through producer ->
serializer/checkpoint -> checkpoint-first payload -> consumer -> renderer.
Local deterministic tests do not establish farm integration. Farm, ROOT/PyROOT,
and full `main.py` behavior require direct farm evidence.

E.8 presentation is the complete visual audit of the kaon missing-mass and
signal-region yield chain: every substantive stage shows its authoritative
input spectrum, applied component/treatment, and authoritative output spectrum,
culminating in canonical `(t,phi)` missing-mass spectra and extracted yields.
It only consumes authoritative runtime objects, persisted snapshots, accepted
detached Method-A artifacts, or later branch outputs. It never recomputes a
fit, factor, correction, normalization, or yield for plotting. Method-A
reweighting remains parent-preserving (`w0 -> w0*C`) with no independent child
normalization; Method B remains numerically excluded.

The current-baseline F.6.3 Method-A full-analysis path has accepted direct
runtime evidence for Q4p4W2p74 / Left / lowe only: real canonical child yields
change while signed parent-t normalization is preserved. Parent equality is a
constraint, not the magnitude of the redistribution effect. See the
[Left/lowe evidence](evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md).
No canonical-five closure or Method-A promotion follows.

E.8 shareable pages consume the existing normalized per-child SIMC MM support
from `hist["_xsect_support_simc"]["mm"]` after the existing SIMC-yield producer.
Clone only for display; never reload, refill, renormalize, or substitute a
setting-wide SIMC shape. Yield summaries consume stored branch values, and
parent closure consumes hash/fingerprint-pinned candidate F.4 records.

An explanatory aggregate page and an applied downstream branch cannot be
interpreted as the same correction lineage without exact source, artifact and
fingerprint identity. Historical accepted E.8.3 F.6.1 aggregates and the
current-baseline F.6.3/E.8.4 candidate are distinct; individual page correctness
does not establish cross-page linkage. Renderer success does not establish
scientific linkage consistency.

Current-branch yield-impact presentations require explicit histogram-scalar
closure against exact producer-owned integration semantics. When fractional
signed-yield changes appear disproportionate to spectrum overlays, report
signed integral, positive-bin support, negative-bin support and absolute
support separately; cancellation is a hypothesis until those numbers exist.
For the present Left/lowe work item, fail-closed identity checks, prior arithmetic
checks and accepted narrow evidence support treating bookkeeping as internally
consistent. Stale scalar/wrong histogram is no longer the primary unresolved
question. Arithmetic consistency does not validate the scientific origin or
production correctness of the large low-t Method-A redistribution; detector-
response validity remains unresolved. The general histogram/scalar rule remains.
SIMC/data amplitude interpretation requires proven same-cell, same-object,
same-normalization and same-unit provenance, including the relevant
integration-window integral. Cloning existing SIMC without renderer scaling
alone does not prove those identities. See the [post-farm investigation](investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md).

Fix.5.4 source audits use normal-bin, unweighted-by-width integration with
flows excluded; exact content/error ownership fingerprints include flow bins.
Current pion-template aggregates use current F.6.3 children only. F.4 parent
application sums cover broader MM support and cannot be equated to cut-window
template integrals. SIMC `iter_weight * normfac / Ncontribute` has no proven
luminosity/effective-charge unit authority in the traced source. Preserve that
absolute-comparison blocker without inventing a scale; data/SIMC agreement or
disagreement cannot be inferred from the current overlay amplitude.
Unresolved absolute SIMC units block the two absolute-SIMC page families and
their interpretation only. They do not make a valid current-F.6.3 E.8.4
data-only/identity/yield/parent-closure payload unavailable.

Fix.5.4 numerical source `761fbb6` (full identity in CURRENT and the phase record) passed
independent ChatGPT actual-diff and pushed-state review per the supplied
Fix.5.5 contract; this is source acceptance, not new farm closure.
Fix.5.5 presentation consumes stored child signed/absolute support and stored
current-F.6.3 aggregates without histogram reintegration or child reaggregation.
Baseline-last solid/open-marker blue and wider dashed magenta encode identity
beyond color. Historical E.8.3 must visibly identify its accepted F.6.1 lineage
as separate from current F.6.3/E.8.4. Fake ROOT checks establish draw commands and
immutability. The later Fix.5.8 farm gate validates the narrow Left/lowe
owner/presentation result and actual PDF legibility; canonical-five closure
and absolute-SIMC interpretation remain unresolved. See the
[accepted closure](evidence/e8-4-fix5-8-left-lowe-runtime-closure-2026-10-03.md).

Fix.5.5 presentation source `df957a6` passed independent actual-diff and
pushed-state review per the Fix.5.6 contract; full identity is in CURRENT and its
phase record. The later Fix.5.8 Left/lowe structural/visual farm gate is accepted;
it does not supply broader numerical interpretation, other-setting acceptance
or production promotion.
An analysis subprocess log contains only its child stream, not later owner
artifact/collection/ZIP failures. Absence of a gate-failure message there is
not evidence of owner success. Fix.5.6 persists atomic per-attempt owner status
outside the unchanged bundle and runs the unchanged detached collector's source
checks before analysis, retaining final collection/rechecks and ZIP verification.

A completed bundle is insufficient by itself: inspect provenance, checker
gates, structured payloads or logs, and rendered pages as applicable. When a
checker disagrees with raw evidence or implementation behavior, inspect the raw
evidence and implementation rather than forcing data to fit the checker. Use
one narrow gate -> one targeted farm run -> fresh evidence -> inspect -> one
coherent repair. See [LEARNINGS.md](LEARNINGS.md) for generalized failure and
review lessons.

Coherent implementation steps may proceed through deterministic local checks
and source review without a farm run after every edit. Farm validation remains
mandatory for ROOT/PyROOT or full-runtime claims and is scheduled at explicit
milestones; inspect fresh artifacts before broadening a narrow runtime gate.

Before a farm handoff, audit input authority, producer (if any), checker (if
any), collector, invocation owner, and returned artifact. Every required
executable step must be tracked, deterministic locally where possible,
independently source reviewed, pushed, and pushed-state reviewed. Unowned
multi-step orchestration blocks farm readiness; one reviewed CLI may own a
complete direct operation. After substantial implementation, review, closure,
or pushed-state handoff, Codex and ChatGPT report memory health. Unresolved
warnings block progress only for material active-state/provenance ambiguity
or a hard integrity violation. Record nonblocking warnings and batch their
repair at the next checkpoint or milestone audit; hard failures remain blocking.

The default loop is memory checkpoint -> substantive science/implementation ->
deterministic local test -> ChatGPT actual-diff review -> user commit/push ->
ChatGPT pushed-state synchronization review -> narrow farm/runtime validation
when required -> fresh evidence review -> memory checkpoint. Commit/push and
pushed-state review are required synchronization stages, not independent
scientific phases. A matching pushed candidate with materially accurate CURRENT
proceeds directly to the substantive gate, without a new reconciliation contract
merely to restate the push.

A memory inconsistency blocks only when it creates concrete ambiguity about
current source identity, accepted evidence, frozen scientific interfaces, active
scientific ownership, or the exact next operation. Otherwise record it, continue
the active milestone, and batch correction at a meaningful checkpoint. Every
substantive scientific loop must produce a plot, yield table, validated numerical
comparison, accepted runtime artifact, or one directly evidenced scientific/
runtime blocker with one coherent repair. Documentation alone is not scientific
completion. See the [workflow decision](decisions/memory-bracketed-scientific-throughput.md).

## Canonical record ownership

- [CURRENT.md](CURRENT.md) owns the sole active objective, blockers, and exact
  next action; the handoff is exceptional transfer state only.
- `evidence/` owns accepted farm/runtime gate facts and artifacts; `phases/`
  owns implementation and fix chronology; `decisions/` owns scientific and
  architectural contracts.
- `roadmap/STATUS.md` owns approved dependency/status structure, not the
  exact active action. `LEARNINGS.md` owns generalized reusable lessons.
- `TOOLS.md`, `COMMUNICATION.md`, and `CODEX.md` own operational references,
  farm-delivery communication, and source-changing workflow respectively.

Use the authoritative [farm validation-bundle procedure](decisions/farm-validation-bundle-procedure.md)
for Jefferson Lab `tcsh` packaging rather than duplicating its procedure here.
