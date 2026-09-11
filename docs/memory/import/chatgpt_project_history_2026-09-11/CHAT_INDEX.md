# KaonLT substantive project-chat index

This is a structured map of recoverable KaonLT project conversations and handoffs, not a verbatim transcript archive.

## 2026-01-22 — early phase-space/plotting work

**Evidence:** `[CHAT_HISTORY]`

Early KaonLT plotting/phase-space work around `t`/`phi` organization and the geometry of multiple `phi` settings.

Durable value: historical context for downstream coordinates only; not a current background-validation milestone.

## 2026-06-11 onward — pion-background subtraction architecture

**Evidence:** `[CHAT_HISTORY] [GIT_HISTORY]`

Major topics:

- `pi_n`, `pi_delta`, `pi_sidis`;
- staged fitting and later joint refinement;
- pion-control to kaon mapping;
- missing-mass windows;
- SIMC component alignment;
- closure and oversubtraction concerns.

Durable conclusions:

- pion subtraction is a distinct production path;
- later HGCer work consumes the frozen baseline instead of rebuilding it;
- the pion component model must not be casually changed to make later diagnostics agree.

## 2026-07 — slow-proton contamination development

**Evidence:** `[CHAT_HISTORY] [HANDOFF]`

Development from an initial broad/global proof-of-concept to local timing/PID event weights.

Topics:

- low-MM slow-proton structure;
- RF/CT/beta timing hierarchy;
- aerogel support;
- momentum/`delta` dependence;
- protected `K Lambda`;
- fit fallback hierarchy;
- proposed versus applied semantics.

Rejected approach: a simple/global model was not trusted because it could improve broad low-MM agreement while removing protected kaon signal.

## 2026-07-15 — K-Lambda reference and offset-gating repair

**Evidence:** `[HANDOFF]`

Repairs false `K Lambda` unavailable state and skipped supported timing cells.

Durable rules:

```text
nominal offset fit
-> expanded offset fit
-> stable-center fallback
```

No hard aerogel 5-NPE production cut.

## 2026-07-16 — dynamic pion/SIMC alignment

**Evidence:** `[HANDOFF]`

Durable design:

- fit region distinct from fixed evaluation envelope;
- explicit baseline candidate;
- common-bin candidate ranking;
- support/localization/boundary/lost-integral diagnostics;
- staged component order;
- shared shifted-template implementation.

## 2026-07-28/29 — particle-background procedure consolidation

**Evidence:** `[HANDOFF] [SCIENTIFIC_REFERENCE]`

Formalized:

```text
random/dummy
-> freeze binning spectra
-> slow proton
-> pion subtraction
```

and the primary code owners for proton, pion, runtime orchestration, and downstream yields.

## 2026-08-11 — Lambda-preservation gate and cleanup

**Evidence:** `[HANDOFF]`

Key decisions:

- setting-wide Lambda safety gate;
- 10% default maximum removal;
- no partial per-`t` application;
- proposed probabilities retained after rejection;
- `applied` means committed post-gate state only;
- fix ROOT/presentation failure modes rather than accepting blank diagnostics.

## 2026-08-19 — review project status / HGCer direction

**Evidence:** `[CHAT_HISTORY]`

Transition from general background work to a structured local HGCer pion-refinement program.

Durable themes:

- proton first, pion second;
- diagnose local acceptance dependence before changing production;
- prefer independent Method-A/Method-B measurements;
- source/diff review over commit-message trust.

## 2026-08-30 — Phase C Method-B plan

**Evidence:** `[HANDOFF]`

Primary record:

```text
KaonLT_Codex_Plan_Phase_C_Method_B_Local_Pion_Background_Closure.md
```

Predecessor state:

```text
Phase A
Phase A.1
Phase A.2
Phase A.3
Phase B
Phase B.1
```

Approved B/B.1 anchor:

```text
18b06ec1aa8ba42859dd0980705ce8374d4720ab
```

Method B must use frozen Phase-A records, remain independent of Method A, protect kaon signals, use same-`t` parent-relative closure, and remain diagnostic.

## 2026-09-01/02 — Phase C review / adaptive Method B

**Evidence:** `[CHAT_HISTORY] [HANDOFF]`

Adaptive baseline-support partition developed as a diagnostic prototype.

It was forbidden from changing legacy Method B or production before data-driven validation.

## 2026-09-04 — persistence/round-trip review discipline

**Evidence:** `[CHAT_HISTORY] [HANDOFF]`

C.Fix.1 revealed that live Method-B fields were omitted by checkpoint serialization.

Permanent review rule:

```text
producer
-> serializer/checkpoint
-> checkpoint-first payload
-> consumer
-> renderer
```

## 2026-09-08 — C.Fix.2.2 handoff / Left-low gate

**Evidence:** `[HANDOFF]`

Accepted state:

```text
c7af5ada2ac6cb342c36cf351b99166ac6b1fbf1
Phase C.fix.2.2
```

Immediate next action then: one `Q4p4W2p74 Left-low` farm validation before C.Fix.2.3.

## 2026-09-09 — Phase-C closure, E.3.Fix.2 Left-low PASS

**Evidence:** `[HANDOFF] [FARM_EVIDENCE]`

Major results:

```text
Phase C: PASS
adaptive Method B: DO NOT PROMOTE
```

E.3 remained independent Method-A presentation only.

Validation workflow moved to:

```text
generic collector + JSON gate profile
```

Corrected Left-low E.3.Fix.2 farm bundle passed source/provenance checks and rendered-page review.

## 2026-09-11 — Phase F.1.Fix.5 / Method-A acceptance mapping

**Evidence:** `[HANDOFF] [CURRENT_SOURCE]`

Scientific decision:

```text
Method B:
retain as diagnostic/cross-check/historical comparison
never use numerically for future pion event-weight adjustment

Method A:
sole candidate future numerical HGCer leakage input
```

Phase-F roadmap:

```text
F.1 event contract
F.2 representation freeze
F.3 detached Method-A map
F.4 detached parent-preserving A-only correction
F.5 detached event propagation
F.6 production promotion
```

Farm-successful F.1 mechanics exposed a scientific bug: F.1 v1 trained on the already-censored `NPE>2` application cache.

Fix.5 separates uncensored Method-A response training from downstream application.

## 2026-09-11 — durable-memory bootstrap

**Evidence:** `[CURRENT_SOURCE]`

Repository now has active memory through `AGENTS.md` and `docs/memory/`.

The bootstrap intentionally refused to guess historical runtime state. This import supplies the deep history needed for the next reconciliation.
