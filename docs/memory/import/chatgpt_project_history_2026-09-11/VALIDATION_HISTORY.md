# KaonLT validation and farm-evidence history

This file separates runtime/farm evidence from source/test evidence.

## Canonical farm workflow

**Evidence:** `[HANDOFF] [FARM_EVIDENCE]`

```text
one narrow gate
-> targeted farm run
-> fresh artifacts
-> inspect provenance/checker/payload/rendered pages
-> PASS or one coherent repair
```

The user performs farm execution.

Preferred tooling:

```text
stable generic collector
+
gate-specific JSON profile
```

`complete=true` is not sufficient. Inspect evaluated commit/provenance, checker gates, required payloads/logs, and actual rendered pages.

## C.Fix.1 — farm-discovered persistence failure

**Evidence:** `[FARM_EVIDENCE]`

Live Method B had `mm_regions` / `protected_regions`; checkpoint serialization omitted them; checkpoint-first rendering therefore lost them.

Conclusion: persisted diagnostics require review of:

```text
producer
-> serializer/checkpoint
-> checkpoint-first payload
-> consumer
-> renderer
```

## C.Fix.2.2 narrow Left-low gate

**Evidence:** `[HANDOFF]`

Accepted source:

```text
c7af5ada2ac6cb342c36cf351b99166ac6b1fbf1
```

Required next gate at that point:

```text
Q4p4W2p74 Left-low only
```

before C.Fix.2.3.

## Phase-C final five-setting validation

**Evidence:** `[FARM_EVIDENCE]`

Accepted C.Fix.2.3:

```text
9a66bc62d20a99172e326e915866877b65ae1e5d
```

Later accepted pre-E.3:

```text
e3853655db0809923cbf2326e2f779219128eda9
```

Settings:

```text
Left lowe
Left highe
Center lowe
Center highe
Right highe
```

No Right-low.

Observed:

```text
30/30 C.Fix.2 pages rendered correctly
Phase C scientific verdict: PASS
adaptive Method B: DO NOT PROMOTE
```

Adaptive coverage:

```text
36/150 candidate-bearing
21/150 multi-slice
114/150 unavailable
```

This closes the Phase-C diagnostic question; it does not promote adaptive B or change production pion subtraction.

## E.3.Fix.2 source review

**Evidence:** `[HANDOFF]`

Implementation:

```text
eb1710f4739ba6ef14f51419806e9fc5bd53c175
```

Source review verified expected scope, fixed canvas/header/grid, final `t1,t2,t3` ordering, local failure behavior, and no scientific producer change.

This is source review, not farm proof.

## E.3.Fix.2 Left-low runtime PASS

**Evidence:** `[FARM_EVIDENCE]`

Runtime/bundle HEAD:

```text
bf53dac84e1396cfd7e3f4e0234749426bfbdcf4
```

Required analysis implementation:

```text
eb1710f4739ba6ef14f51419806e9fc5bd53c175
```

Collector verified implementation in runtime ancestry.

Bundle contained only:

```text
manifest.json
source_state.txt
source_checks.txt
Left_lowe/
  checkpoint JSON
  E3-validation PDF
```

Recorded source checks:

```text
required analysis ancestry: PASS
py_compile: PASS
full-background-subtraction tests: 61 OK
pion HGCer refinement tests: 63 OK
bundle collector tests: 16 OK
git diff --check: PASS
```

Source PDF had 61 pages. Collector extracted 59–61 as:

```text
E.3 t1
E.3 t2
E.3 t3
```

Visual review passed:

- complete 5x2 grids;
- all 10 delta cells/page;
- no clipping/header overlap;
- `NPE=2` visible;
- low/control labels visible;
- negative signed content retained;
- unavailable cells explicit;
- no Method-B/correction language.

Checkpoint provenance spot checks matched.

Conclusion for this exact gate:

```text
E.3.Fix.2 Q4p4W2p74 Left-low
CLOSED / RUNTIME VALIDATED
```

Do not extrapolate to the remaining four E.3 settings without evidence.

## F.1.Fix.3 farm failure

**Evidence:** `[FARM_EVIDENCE]`

Observed:

```text
ValueError:
The truth value of an array with more than one element is ambiguous
```

Cause: NumPy-like edge arrays in truth-value fallback expressions.

Repair: safe edge materialization.

Later farm output closed this mechanical regression.

## F.1.Fix.4 farm failure

**Evidence:** `[FARM_EVIDENCE]`

Observed:

```text
parent_child_parity_mismatch:allcuts
```

Cause: native Python `bool` versus NumPy boolean scalar with identity-style parity logic.

Repair: detach NumPy-like child scalar to native Python scalar at F.1 read boundary.

Later farm output closed this mechanical regression.

## Post-Fix.4 five-setting mechanical gate

**Evidence:** `[FARM_EVIDENCE]`

All five settings produced mechanically valid F.1 artifacts:

```text
status=available
diagnostic_stage=complete
records nonempty
fingerprints populated
five F.1 pages rendered
```

Focused checks:

```text
full-background tests: 90 PASS
F.1 contract tests: 10 PASS
Phase-F runtime tests: 3 PASS
bundle collector tests: 14 PASS
py_compile: PASS
git diff --check: PASS
```

Interpretation: Fix.3/Fix.4 mechanics were validated.

This did **not** scientifically validate F.1 v1. Successful output exposed the population-ownership bug.

## F.1 population-ownership failure

**Evidence:** `[FARM_EVIDENCE]`

Method-A counts:

```text
Left low:    52397 positive, 576 low, 51821 control
Left high:   14862 positive, 285 low, 14577 control
Center low:  37035 positive, 568 low, 36467 control
Center high: 22429 positive, 479 low, 21950 control
Right high:  23963 positive, 736 low, 23227 control
```

F.1 v1:

```text
low = 0
control = exactly Method-A control
```

for every setting.

Interpretation: F.1 used the downstream application cache already censored at `NPE>2`. The application cache remained correct; it was the wrong scientific owner for Method-A training.

This motivated F.1.Fix.5 dual populations.

## F.1.Fix.5 current validation state

**Evidence:** `[CURRENT_SOURCE]`

Observed current `test` HEAD:

```text
7cdb7847d72501cd3dc504565fafc604e22a6132
```

This proves Fix.5-era source exists on `test`.

It does not prove:

- a farm run at this HEAD;
- ROOT rendering PASS;
- F.1 v2 five-setting closure;
- readiness for F.2.

Migration should inspect live Fix.5 source and use a source-level/farm-pending status unless newer farm evidence is present.
