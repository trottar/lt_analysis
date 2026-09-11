# Approved KaonLT roadmap

This roadmap separates implemented scope from acceptance status. It is based
on the supplied project history and the reviewed current source; it does not
promote a phase because a commit or test exists.

## Implemented scope recorded by history/source

The established background program includes random subtraction, pion-background
handling, slow-proton contamination treatment, HGCer diagnostics, Method-A /
Method-B comparison infrastructure, validation/checker infrastructure, and
later presentation of frozen diagnostics. These are architecture facts, not
runtime-status claims. Slow-proton and pion-background work remain separate.

The source also identifies the HGCer sequence: independent Method A and Method
B diagnostics; Phase-D comparison/closure infrastructure; Phase-E
presentation-only work; and a Phase-F.1 Method-A acceptance diagnostic/farm
review path. Earlier phase acceptance is intentionally unclassified pending
farm evidence.

## `CLOSED / RUNTIME VALIDATED`

No phase is recorded in this category. No farm evidence was supplied or found
in repository-owned durable memory during initialization.

## `SOURCE REVIEWED`

Phase F.1 Method-A acceptance diagnostic and its renderer/collector path are
source reviewed at live `test` head
`cf9c804b55a501af6839b8fb3b0358c56853c276`. It is non-authoritative and does
not license a production correction. See
`../phases/phase-f1-method-a-acceptance-farm-gate.md`.

## `ACTIVE`

Phase-F.1 source-identity reconciliation for the farm-review profile. Its
present profile pins `d656e15761970d7d612bb028d2746d077795e9ad`, while the
live source contains later non-allowlisted analysis changes. The current
collector therefore cannot constitute a valid gate for live `test` until that
identity rule is reconciled.

## `DEFERRED`

Import authoritative farm history before assigning runtime status to the broad
background program, Method-A/Method-B work, Phase D, or Phase E. This record
does not invent future scientific phases.

## Future scope

No future phase is approved in this record. Add one only through an explicit
project decision or authoritative handoff; do not infer it from source history.

## `NEXT`

Establish one narrow Phase-F.1 source-identity reconciliation contract, then
run one targeted farm gate with fresh artifacts and rendered-page inspection.
