# F.1 v1 population-ownership failure

## Evidence

[FARM_EVIDENCE] [HANDOFF]

The post-Fix.4 five-setting mechanical output compared Method-A counts to F.1
v1 results:

| Setting | positive | low | control |
|---|---:|---:|---:|
| Left low | 52397 | 576 | 51821 |
| Left high | 14862 | 285 | 14577 |
| Center low | 37035 | 568 | 36467 |
| Center high | 22429 | 479 | 21950 |
| Right high | 23963 | 736 | 23227 |

For every setting, F.1 v1 reported low=0 and control exactly equal to the
Method-A control count.

## Scientific diagnosis

F.1 v1 had consumed Phase-A pion_records from the authoritative downstream
physical pion-control cache, which correctly already restricts events to NPE>2.
That cache was correct for downstream application, but it cannot train a
Method-A response measurement requiring 0<NPE<=2 low observations.

## Corrective decision and current boundary

This evidence motivated the F.1.Fix.5 dual-population architecture: response
training uses prompt/noRF/nommcuts NPE>0 Part-1 records joined to acceptance
records by source_label and entry_index; application retains the authoritative
physical NPE>2 parent/child cache. The live source is SOURCE REVIEWED for this
architecture. F.1.Fix.5 farm evidence is absent, so its runtime status remains
ACTIVE pending the identity-reconciled targeted farm gate.

