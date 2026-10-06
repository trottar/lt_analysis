---
memory_schema: 3
---
# Current KaonLT development state

## Active Objective

E.8 is `ACTIVE`: audit the kaon missing-mass and signal-region yield chain
from authoritative objects. Preserve baseline production, accepted upstream
authorities and detached Method-A/Method-B boundaries.

## Current Work Item

The fresh-F1 canonical-five lineage preflight is
`CLOSED / RUNTIME VALIDATED` for pre-analysis lineage/isolation/preservation only.
The [accepted preflight](evidence/e8-4-canonical-five-lineage-preflight-runtime-closure-2026-10-05.md)
records reviewed materialization/staging, five real F.4/F.5 reconstructions,
cleanup and final preservation with analysis_started=false. Acceptance comes from the supplied receipt and prior review; Codex ran no farm gate.

The second full canonical-five runtime gate is `BLOCKED`: analysis completed
(returncode 0), but artifact/page verification failed. All five fresh manifests
contained only `full_background.e8_4.unavailable` for E.8.4, with
`f6_3_f4_shared_reproduction_failed:f3_fingerprint_input_content_mismatch`.
RUNTIME VERIFIED (supplied diagnosis): analysis regenerated different F.1
raw/stable identities after successful pre-analysis reproduction; the tracked
post-run comparator reconstructed matching F.2/F.3/F.4 scientific payloads,
`first_changed_stage = none`. INFERENCE: this gate exposes provenance/identity
coupling, not a newly established Method-A scientific defect. See the
[second-run diagnostic evidence](evidence/e8-4-canonical-five-second-full-run-provenance-blocker-2026-10-05.md).

Canonical-five provenance/identity repair is `DEFERRED` by user decision.
Full runtime/PDF closure and the low-level cause remain NOT VERIFIED. The
[first-run failure/materialization](evidence/e8-4-canonical-five-fresh-f1-lineage-materialization-2026-10-05.md)
is historical diagnostic evidence with a distinct F.4 comparison result.

The detached Left/lowe diagnostic remains `CLOSED / RUNTIME VALIDATED` for
parents t0/t1/t2, primary t1. Its [accepted checkpoint](evidence/e8-4-left-lowe-method-a-current-lineage-diagnostic-runtime-closure-2026-10-04.md)
records prior ChatGPT acceptance of supplied JSON/PDF, exact shared F.4
reproduction and eight-page review. Control identities match
(10716/18749/22356), t1 normalization contrast is -1.02%, parent residual 0.0,
OOD 21/19809 (~0.106%). INFERENCE: with accepted F.6.2 acceptance/MM/bootstrap
validation and comparator t1 continuity, no correction redesign is warranted.
Absolute mis-ID, unique zero-response proxy validity, hardware cause and
promotion remain NOT VERIFIED. Gates 1–4, external research and diagnostic interpretation are consumed.

## Verified State

SOURCE VERIFIED (Gate-3 ownership retained): NPE=0 is kaon PID, NPE>0 the
pion tree, NPE>2 physical pion control. F.3 `hgcer3` is relative response, not
absolute leakage probability. F.4 preserves signed canonical-t parents, not
children/MM subregions; baseline `w0` owns control-to-background transfer.
The HGCer hole is excluded; zero-photoelectron transfer is diagnostic-only,
slow-proton architecture analogue-only. Details remain in the Gate-3 reference.

RUNTIME VERIFIED (prior accepted reviews, not new artifact inspection): the
current-baseline comparator reproduced F.2/F.3 exactly, with F.4 first changed.
Left/lowe t1 F.4 differs at floating-point scale; substantive t0/t2 changes do
not inherit that continuity claim. F.6.2 remains acceptance-correlated/MM
validation, not absolute HGC calibration. F.4.Refresh.2 remains
`CLOSED / RUNTIME VALIDATED`; F.6.3/E.8.4 retain that status only for Left/lowe
execution, live-cache parity, real child changes and signed parent preservation.
See [branch evidence](evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md).

Fix.5.7/Fix.5.8 remain `CLOSED / RUNTIME VALIDATED` only for Left/lowe owner
provenance and presentation legibility respectively. [Fix.5.8 evidence](evidence/e8-4-fix5-8-left-lowe-runtime-closure-2026-10-03.md)
retains its 97-page structural/parent/visual scope. Earlier F-stage closures
retain their scopes in the [roadmap](roadmap/STATUS.md). E.8.2/E.8.3 and workflow
hardening remain `SOURCE REVIEWED`; historical E.8.3 F.6.1 and current F.6.3/E.8.4
are distinct lineages. Literature adds no runtime validation; no authority is
replaced.

## Source / Evidence Identity

- Second full-run farm source and this checkpoint's startup observation,
  2026-10-05, branch `test`, HEAD/local `origin/test`:
  `ace8688a71431d13b40ed19713a27746f3da6a8e`.
  Failed gate-status SHA-256:
  `ddfadaf9946e262c35795689ab8d2d52a306ede9e775d98fb05df4e15924e91f`.
  Diagnostic evidence only; no full-run acceptance.
- Accepted separate preflight source: `b8bafd3f2853523ca9deb7aa7d6b572357584d6f`.
- Accepted Left/lowe diagnostic source: `aad27a4d1639eef188dc61563fd3615682835835`.
- Accepted Fix.5.8 source: `2ddeab47d55edb57d2f313022a948c4376730c19`.
- Current F.6.3 F.4 SHA-256:
  `79e7ceda7221cbeeead4ed5bc306b0e0e670741a27beaa980e22349c555e96d7`.

Linked evidence owns receipt/artifact identities, prior timestamped observations
and dirty farm provenance; no global farm cleanliness is claimed.

## Blockers

Canonical-five full runtime gate is `BLOCKED`; provenance/identity repair is
`DEFERRED`. Final E.8 and F.6.4 remain `BLOCKED`. The accepted diagnostic
consumes prior unmeasured t1 normalization/support/population explanations;
absolute mis-ID, hardware cause and production correctness remain unproven.
RF was not performed; it is optional at low epsilon, not universal.

Absolute-SIMC interpretation remains separately `BLOCKED`:
`SIMC_normfac_luminosity_and_charge_units_not_source_proven`.
`iter_weight * normfac / Ncontribute` lacks proven luminosity/effective-charge
units. Only two absolute-SIMC page families/claims are unavailable; valid
current-F.6.3 data/identity/yield/parent-closure payloads remain available.
No normalization error, conversion or amplitude conclusion is established.

## Next Action

NEXT — after actual-diff review, user commit/push and pushed-state
synchronization, resume the E.8.2 scientific audit of the authoritative
Q4p4W2p74 baseline kaon missing-mass and stage-yield chain. Trace prompt/random
-> dummy -> slow-proton cleaning -> baseline pion subtraction -> final
canonical clean-kaon missing mass and Y0(t,phi); identify where actual
spectrum/yield changes occur before deciding whether new implementation or a
farm gate is required. This is science/audit work, not a new maintenance phase.
Canonical-five provenance repair stays `DEFERRED`; this checkpoint authorizes
no farm command or new Method-A design. Subsequent E.8.3/E.8.4 interpretation
follows the approved roadmap unless fresh evidence establishes a blocker.

## Success Criteria

Manifest, ordinary memory health, bootstrap, allowlist/byte-preservation and
complete-diff checks must pass. Memory-only checks add no runtime acceptance.
No staging, commit, push or farm execution is authorized for Codex.

## Do Not Reopen Without New Evidence

Keep `no_empirical_residual` and zero legacy residual scales. Method A stays
detached/non-production; Method B diagnostic/cross-check only, numerically
excluded. Freeze random/dummy/slow-proton/baseline-pion subtraction, SIMC,
weights, yields, uncertainties, cuts, templates, priors, binning, efficiencies,
acceptance, L/T and cross sections. No independent child normalization,
authority replacement or automatic promotion. E.8 consumes; F.6.3 owns the
private branch. Failed artifacts are inadmissible for closure; owner success
alone does not establish scientific/visual acceptance.

## Relevant References

- [Gate-3 source/science investigation](investigations/e8-4-left-lowe-method-a-detector-response-source-science-audit-2026-10-03.md)
- [Current-baseline comparator](evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md)
- [Prior Fix.5.7 evidence](evidence/e8-4-fix5-7-left-lowe-runtime-and-fix5-visual-gate-2026-10-03.md)
- [E.8 procedure roadmap](decisions/e8-full-analysis-procedure-roadmap.md)
