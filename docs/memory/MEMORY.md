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
