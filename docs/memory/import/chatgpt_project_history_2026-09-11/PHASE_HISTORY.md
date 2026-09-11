# KaonLT phase/fix history

This record focuses on the HGCer refinement sequence and the directly related production-background infrastructure.

## Pre-HGCer production background program

### Slow-proton subtraction

**Evidence:** `[CHAT_HISTORY] [HANDOFF]`

Scientific owner: PID contamination inside kaon-selected data.

Frozen principles:

- event-level timing/PID treatment;
- preserve protected `K Lambda`;
- setting-wide commitment gate;
- separate proposed and applied quantities;
- bypass whole setting when safety gate fails;
- remain scientifically separate from pion subtraction.

### Pion component subtraction

**Evidence:** `[GIT_HISTORY] [CHAT_HISTORY] [SCIENTIFIC_REFERENCE]`

Scientific owner: pion-production events mapped from pion-control data.

Frozen principles:

- staged `pi_n`, `pi_delta`, `pi_sidis` model;
- dynamic alignment/joint refinement belongs to this path;
- pion treatment follows slow-proton treatment;
- HGCer diagnostics consume the frozen baseline instead of redefining it.

# HGCer phased program

## Phase A — exact event-level baseline contract

**Evidence:** `[HANDOFF]`

Purpose: expose frozen pion baseline and actual host state without changing production.

Recorded as completed before Phase C.

## Phase A.1 — noRF/frozen-weight/provenance repairs

**Evidence:** `[HANDOFF]`

Recorded completed before Phase C.

## Phase A.2 — identity-host runtime wiring

**Evidence:** `[HANDOFF]`

Recorded completed before Phase C.

## Phase A.3 — independent upstream noRF identity-host closure

**Evidence:** `[HANDOFF]`

Recorded completed before Phase C.

## Phase B — Method-A HGCer-response diagnostics

**Evidence:** `[HANDOFF]`

Independent detector-response path.

## Phase B.1 — statistical/diagnostic hardening

**Evidence:** `[HANDOFF]`

Approved pre-C source:

```text
18b06ec1aa8ba42859dd0980705ce8374d4720ab
```

## Phase C — Method-B local pion-background closure

**Evidence:** `[HANDOFF]`

Purpose: independent local MM closure against actual noRF kaon host.

Forbidden:

- Method-A numerical inputs;
- tree reopening;
- rebuilt pion weights;
- rerun proton cleaning;
- absolute neutron normalization;
- correction application.

### C.Fix.1 family — persistence repair

**Evidence:** `[HANDOFF] [FARM_EVIDENCE]`

Problem: live region definitions were not serialized, so checkpoint-first rendering lost them.

Permanent result: persisted-diagnostic round-trip review rule.

### C.Fix.2 — adaptive support partition

**Evidence:** `[HANDOFF]`

Diagnostic prototype only.

### C.Fix.2.2

**Evidence:** `[HANDOFF]`

```text
c7af5ada2ac6cb342c36cf351b99166ac6b1fbf1
```

Stopped for one Left-low farm gate.

### C.Fix.2.3 / final Phase-C closure

**Evidence:** `[HANDOFF] [FARM_EVIDENCE]`

```text
implementation:
9a66bc62d20a99172e326e915866877b65ae1e5d

later accepted pre-E.3:
e3853655db0809923cbf2326e2f779219128eda9
```

Five settings:

```text
Left lowe
Left highe
Center lowe
Center highe
Right highe
```

Result:

```text
30/30 pages rendered
Phase C PASS
adaptive Method B DO NOT PROMOTE
```

Canonical migration should treat Phase C as closed after validating the imported evidence record.

## Phase D — A/B comparison

**Evidence:** `[HANDOFF]`

Uses legacy Method B, never adaptive B.

Frozen availability states:

```text
Both available
Both present; ratio undefined
Method A only
Method B only
Neither
```

When positive/comparable:

```text
r_BA = B/A
Delta_AB = ln(B)-ln(A)
```

Later project state treats D as frozen upstream comparison/closure infrastructure.

**Gap:** no single final Phase-D farm manifest is present in this import. Do not invent exact runtime provenance.

## Phase E — presentation over frozen upstream results

**Evidence:** `[HANDOFF]`

Presentation-only ownership.

## Phase E.3 — independent local Method-A presentation

**Evidence:** `[HANDOFF]`

No Method-B, adaptive-B, Phase-D numerical input, or correction.

### E.3 initial implementation

```text
517990b5a150b3b2ad2d52f50f3d01de62f4d7d0
```

### E.3.Fix.1

```text
e5f1ca4c15ae3562ed06eaf8aa0b9382cfe7acc6
```

### E.3.Fix.2

**Evidence:** `[HANDOFF] [FARM_EVIDENCE]`

```text
eb1710f4739ba6ef14f51419806e9fc5bd53c175
```

Repairs header overlap and final-page ordering.

Later validation-infrastructure runtime HEAD:

```text
bf53dac84e1396cfd7e3f4e0234749426bfbdcf4
```

Specific validated gate:

```text
Q4p4W2p74 Left-low
CLOSED / RUNTIME VALIDATED
```

**Gap:** recovered evidence does not independently establish final remaining-four-setting E.3 closure.

# Phase F — Method-A acceptance mapping

## F.1 — detached Method-A acceptance event contract

No estimator, correction, or production mutation.

Recent lineage:

```text
b0f6868c17b0aae1d97293d4a103a3736fd50fd0  pre-F.1
b49dd1cdaf579407c5375d116e048dbf7681a4f7  F.1
4fa198e4e43b776f18ef2db7fc918df18857da27  F.1.fix.1
d656e15761970d7d612bb028d2746d077795e9ad  F.1.fix.2
81b25890ff5e600617db4040cd99e3c875cb9f7c  validation bundle infra
8b1ad5b735e1b5cc95a5d876526d321296cec16b  F.1.fix.3
eb253046ff6e23ec94c8638315c4bc712fa4f992  F.1.fix.4
```

### F.1.Fix.1

Repairs Phase-A fingerprint compatibility, noRF provenance, and true `(t,phi)` page-5 mapping.

### F.1.Fix.2

Presentation-only page-5 matrix.

### F.1.Fix.3

**Evidence:** `[FARM_EVIDENCE]`

Farm found NumPy array truth-value error. Safe edge materialization repaired it. Later farm gate closed this regression.

### F.1.Fix.4

**Evidence:** `[FARM_EVIDENCE]`

Farm found Python-bool versus NumPy-bool parity mismatch. Native-scalar detachment repaired it. Later five-setting gate closed this regression.

### Post-Fix.4 mechanical farm gate

**Evidence:** `[FARM_EVIDENCE]`

All five settings produced mechanically valid F.1 output and focused test suites passed.

This did **not** scientifically validate F.1 v1; it exposed the training-population error.

### F.1.Fix.5 — dual populations

**Evidence:** `[HANDOFF] [CURRENT_SOURCE]`

Training population:

```text
Part-1 prompt/noRF/nommcuts NPE>0 response records
joined with acceptance records by (source_label, entry_index)
```

Application population:

```text
authoritative physical pion-control NPE>2 parent/child cache
```

Training owns pages 1–4; application owns page 5 and future event application.

Approved six-file implementation scope:

```text
src/cuts/pion_hgcer_method_a_acceptance_contract.py
src/cuts/rand_sub.py
src/cuts/full_background_subtraction_plots.py
testing/test_pion_hgcer_method_a_acceptance_contract.py
testing/test_pion_hgcer_phase_f_runtime_contract.py
testing/test_full_background_subtraction_plots.py
```

Expected v2 schemas:

```text
pion_hgcer_method_a_acceptance_event_contract/v2
pion_hgcer_method_a_acceptance_event_contract_fingerprint/v2
pion_hgcer_method_a_acceptance_event_contract_artifact/v2
full_background_subtraction_f1/v2
```

Current observed `test` HEAD:

```text
7cdb7847d72501cd3dc504565fafc604e22a6132
```

**Status caution:** source presence is not farm validation. Migration must inspect actual source and retain a source-level/farm-pending status unless newer runtime evidence exists.

## F.2 — freeze probability-map representation

Future after F.1 support review. Candidate reduced bases may use `(delta,xptar,yptar)` or `(delta,HGCer x,HGCer y)` depending on support.

Do not blindly select a 5D histogram.

## F.3 — detached Method-A map

Future. Requires support/OOD handling. No Method-B numerical input.

## F.4 — detached parent-preserving A-only correction

Future. Construct `C_A`; remain detached.

## F.5 — detached event-level `(t,phi)` propagation

Future. Build parallel baseline and A-adjusted templates; no production promotion.

## F.6 — production promotion

Future only after F.5 scientific validation and explicit promotion decision.
