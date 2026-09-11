# KaonLT durable decision history

## D001 — slow-proton and pion backgrounds are distinct

**Evidence:** `[SCIENTIFIC_REFERENCE] [CHAT_HISTORY]`

Do not combine slow-proton PID contamination and pion-production background into one empirical MM correction.

## D002 — proton correction precedes pion subtraction

**Evidence:** `[SCIENTIFIC_REFERENCE] [HANDOFF]`

```text
random/dummy
-> freeze binning spectra
-> slow-proton cleaning
-> pion subtraction
```

## D003 — K-Lambda preservation is a production safety gate

**Evidence:** `[HANDOFF] [FARM_EVIDENCE]`

A broad low-MM improvement is insufficient if genuine `K Lambda` is removed. Use a setting-wide gate with default maximum removal:

```text
10%
```

Gate failure bypasses the whole-setting proton correction while retaining proposed diagnostics.

## D004 — proposed versus applied proton quantities must stay distinct

**Evidence:** `[HANDOFF]`

`applied` means committed post-gate production state only.

## D005 — dynamic pion-alignment candidates use a fixed evaluation envelope

**Evidence:** `[HANDOFF]`

Moving fit windows may determine amplitudes; candidate comparison must use common evaluation bins.

## D006 — HGCer diagnostics consume the existing pion baseline

**Evidence:** `[CHAT_HISTORY] [HANDOFF]`

Do not redesign the pion component model merely to improve later HGCer agreement.

## D007 — Method A and Method B remain independent through comparison

**Evidence:** `[HANDOFF]`

Method A is detector-response driven. Method B is local MM-closure driven. No Method-A numerical values enter Method B.

## D008 — Method B is parent-relative within canonical t

**Evidence:** `[HANDOFF] [FARM_EVIDENCE]`

Use same-`t` parent normalization. Do not normalize across `t`, interpolate missing cells, or use the neutron peak as an absolute pion normalization anchor.

## D009 — adaptive Method B is not promoted

**Evidence:** `[FARM_EVIDENCE]`

```text
DO NOT PROMOTE
```

Reason: stable where measurable but too sparse and setting dependent.

Consequences:

```text
legacy Method B remains Phase-D reference
adaptive B remains diagnostic context
no C.Fix.2.4
no D.Fix.1
no production pion change
```

## D010 — persisted diagnostics require full round-trip review

**Evidence:** `[FARM_EVIDENCE]`

```text
producer
-> serializer/checkpoint
-> checkpoint-first payload
-> consumer
-> renderer
```

## D011 — Phase-E presentation does not own physics

**Evidence:** `[HANDOFF]`

Presentation consumes frozen upstream values; it cannot recompute Method A, insert Method-B numerical input, or construct corrections.

## D012 — generic collector + JSON profile is the validation architecture

**Evidence:** `[FARM_EVIDENCE]`

Change collector Python only for new capabilities; use profiles for new gates.

## D013 — validation bundles package their own review artifacts

**Evidence:** `[FARM_EVIDENCE]`

Do not make the user manually extract pages/PDFs/checkpoints when the collector can package them.

## D014 — stale checkout/source identity is a tooling responsibility

**Evidence:** `[FARM_EVIDENCE]`

Validation tooling must expose and machine-check source provenance.

## D015 — the original A+B correction plan is superseded

**Evidence:** `[HANDOFF] [SUPERSEDED]`

Earlier plan: A and B might both constrain a future correction.

Later authoritative decision:

```text
Method B:
diagnostic/cross-check/historical comparison only

Method A:
sole candidate numerical future HGCer leakage input
```

## D016 — Method A is not yet an absolute leakage probability

**Evidence:** `[HANDOFF]`

Current Method A excludes `NPE=0`. It is a positive-response relative leakage diagnostic until a dedicated zero/nonpositive-response checkpoint justifies stronger language.

## D017 — Method-A training and future application populations are distinct

**Evidence:** `[FARM_EVIDENCE] [HANDOFF]`

Training:

```text
prompt, noRF, nommcuts, NPE>0
```

Application:

```text
authoritative physical pion-control NPE>2 cache
```

F.1 v1 proved the censored application cache cannot own response training.

## D018 — response coordinate and yield coordinate intentionally differ

**Evidence:** `[HANDOFF]`

Response may be learned in:

```text
canonical t parent + acceptance coordinates
```

Application occurs event-by-event before filling:

```text
canonical (t,phi) children
```

## D019 — future normalization is parent-level, not child-level

**Evidence:** `[HANDOFF]`

Preserve baseline normalization across each canonical `t` parent. Do not separately renormalize `(t,phi)` children.

## D020 — no production promotion before F.6

**Evidence:** `[HANDOFF]`

F.1–F.5 remain detached. Production promotion is an explicit F.6 decision after F.5 validation.
