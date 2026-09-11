# KaonLT project history

This is the canonical chronological history reconstructed during the
2026-09-11 import. The import package remains at
docs/memory/import/chatgpt_project_history_2026-09-11 as source provenance.
Current implementation and runtime state are deliberately kept in active
memory and evidence records rather than inferred from chronology.

## Analysis foundation and production background work

[SCIENTIFIC_REFERENCE] [CHAT_HISTORY] [GIT_HISTORY]

KaonLT is the Hall-C E12-09-011 coincidence analysis. Its enduring layers
include reconstruction/PID, random subtraction, fixed acceptance/binning, SIMC
comparison, yield and average-kinematics extraction, cross-section work, L/T
separation, and uncertainties. The Trotta and Usman dissertations provide
scientific context; they do not override current source.

The June 2026 pion-background history developed staged pi_n -> pi_delta ->
pi_sidis modeling and SIMC-background inputs. Representative commits
395451f..., 6cfdd47..., 2a19bde..., 9db45f2..., 64ea04e..., 7a758491...,
e1640390..., ab43cac0..., 8284f087..., and 277f269f... are chronology only.
Dynamic pion/SIMC alignment later adopted a fixed comparison envelope, explicit
baseline candidate, common-bin ranking, support/localization/boundary checks,
and staged scans. It remains pion-model work, not later HGCer diagnostics.

The durable production ordering is:

    random/dummy correction
    -> freeze bin-counting spectra
    -> slow-proton cleaning
    -> pion component subtraction
    -> optional residual diagnostics

Slow protons are PID contamination inside kaon-selected data; pion background
is a separate pion-production population mapped from pion-control data. A
rejected broad/global proton proof of concept improved low-MM regions but could
remove real K Lambda. It was superseded by local timing/PID event treatment and
a setting-wide Lambda-preservation gate. See decisions/DECISION_HISTORY.md.

## HGCer diagnostic program

[HANDOFF] [CHAT_HISTORY] [FARM_EVIDENCE]

Phase A exposed frozen baseline pion prediction and host state without
production mutation. A.1/A.2/A.3 repaired noRF/frozen-weight provenance and
identity-host closure. Phase B/B.1 established the independent Method-A
positive-response diagnostic; its pre-C source anchor is
18b06ec1aa8ba42859dd0980705ce8374d4720ab.

Phase C implemented independent Method-B local missing-mass closure from frozen
Phase-A records. C.Fix.1 demonstrated that a live diagnostic is not enough
when its serializer omits fields consumed by checkpoint-first pages. C.Fix.2
investigated an adaptive support partition. The five-setting farm closure
recorded 30/30 rendered pages and found the adaptive prototype stable where
measurable but too sparse for promotion. Its evidence is in
evidence/phase-c-five-setting-closure.md.

Phase D compares frozen, independent A/B outputs. It retains availability
states rather than inventing ratios for zero/unavailable values and uses legacy
Method B, not adaptive Method B. Final Phase-D farm provenance was not
recovered. Phase E is presentation-only over frozen records; the one recovered
E.3.Fix.2 Left-low gate is runtime validated, but does not close the other four
settings. See evidence/e3-fix2-left-low-runtime.md.

## Phase F lineage and current source checkpoint

[HANDOFF] [FARM_EVIDENCE] [CURRENT_SOURCE]

Phase F is a staged, detached Method-A acceptance program:

    F.1 event contract
    F.2 representation freeze
    F.3 detached Method-A map
    F.4 detached parent-preserving A-only correction
    F.5 detached event propagation
    F.6 explicit production promotion only after F.5 validation

F.1 lineage is b0f6868c... (pre-F.1), b49dd1cd... (F.1), 4fa198e4...
(Fix.1), d656e157... (Fix.2), 8b1ad5b7... (Fix.3), and eb253046...
(Fix.4). Farm testing repaired NumPy edge-array truth evaluation and NumPy-bool
parity; a subsequent five-setting mechanical gate revealed rather than hid the
F.1 v1 population-ownership error. See evidence/f1-fix3-fix4-mechanics.md and
evidence/f1-v1-population-ownership-failure.md.

F.1.Fix.5 separates uncensored Method-A response training from the physically
censored downstream application cache. The actual source change is
dc4fc6283001739a487ec80068f951b0e388cae6; later commits through live
b02316f83bf8ef18641fa7217b90c04bcdca10e3 add durable-memory/import material.
The source review is in phases/phase-f1-method-a-acceptance-farm-gate.md. No
F.1.Fix.5 farm artifact was supplied in this migration.

