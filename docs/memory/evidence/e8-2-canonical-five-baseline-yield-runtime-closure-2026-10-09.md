# E.8.2 canonical-five baseline yield/MM runtime closure — 2026-10-09

Status: `CLOSED / RUNTIME VALIDATED` only for the independent Q4p4W2p74
canonical-five baseline extracted-yield and matching MM numerical/procedure gate.

## Acceptance and exact identities

RUNTIME VERIFIED (attributed): user-supplied Jefferson Lab farm artifacts and
independent prior ChatGPT archive/status/receipt/hash review own this acceptance.
The [closure contract](../phases/e8-2-canonical-five-baseline-yield-runtime-closure-task-contract.md)
transfers that completed acceptance to repository memory. Codex records it
without rerunning, reinterpreting or independently validating external artifacts.

Farm source: `8e68d2acae91fd97622793df563aba0296403a48`.
Attempt: `Q4p4W2p74_E8_2_baseline_five_8e68d2ac_20261009T223015Z`.

| Supplied artifact | Bytes | SHA-256 |
| --- | ---: | --- |
| Validation ZIP | 111668823 | `df5a64563c7322416f0895c52557bea69f49d73ac0cb47836b0635a27e7660c4` |
| Companion gate-status JSON | 20273 | `a8cdfa6ee35a461f3313f63c1b8faa7c4db5416c9dfffafff5f87b4f986fad29` |
| Companion run-summary JSON | 965570 | `856128f6c41f27e5a09a52aff01c37825c8d9de7279cda36d11009f3e0d2f583` |
| Companion child log | 292966650 | `f630c33b5d759c3c7541fe552f83524a63edb923ce7c9e51c05badeaeba6c282` |

The accepted owner receipt reports `status=success`, `stage=complete`,
`failure_reason=null`, analysis subprocess return code 0, completed analysis,
artifact checks, collection, ZIP verification and disposable-worktree cleanup.
Primary checkout identity/ref/index/content and original installed-ltsep
preservation passed, including final checks. The full launcher ran only in
the isolated disposable worktree. This does not assert a currently clean
ordinary checkout.

## Accepted inventory and presentation scope

Exactly five settings: `Left/lowe`, `Left/highe`, `Center/lowe`, `Center/highe`,
`Right/highe`; no Right/lowe. Both epsilon analyses used `no_empirical_residual`
and existing 0th-iteration cut-selected MM support.

ZIP CRC, manifest/source/provenance and all 29 required artifact hashes passed.
Required inventory: five procedure PDFs, five page manifests, two epsilon
full-analysis JSON files, two JSON and two CSV correction ledgers, two source
NPZ MM support files, ten byte-for-byte producer yield `.dat` tables and one
canonical owner receipt. The ZIP has 38 entries including directories and
three collector control entries. Packaged owner receipt and companion summary
are byte-identical.

Five PDFs contain 97 pages each (485 total). Matching manifests report zero
renderer failures and all 18 E.8.2 stage/final pages per setting (90 records).
Representative final-MM/stage-yield pages for each setting were rendered and
inspected. Structural inventory and spot-checked readability are accepted;
independent pixel-by-pixel visual acceptance of all 485 pages is not claimed.

## Accepted numerical support

Exactly 135 `(t_index,phi_index)` cells: three canonical t parents × nine phi
children × five settings. Producer-to-cut-selected-NPZ gates passed for
135/135 data and 135/135 SIMC cells, using one-based producer `(phi_bin,t_bin)`
to zero-based NPZ child mapping, normal-bin integration and exact existing
four-decimal text formatting. Prior independent receipt reintegration found
zero discrepancies. Source table uncertainties were neither changed nor
rederived as a scientific replacement.

The cut-selected MM window is 1.10–1.16 GeV, with 100 normal bins per child.
These integrated numerical-support objects are distinct from wider stage-display
histograms. Diagnostic Lambda-window stage integrals, including after-pion
`PI`, are not authoritative extracted baseline `Y0` or another subtraction stage.

Four signed-negative producer data yields are retained (one-based producer
`(t_bin, phi_bin)` values):

| Setting | t_bin (one-based) | phi_bin (one-based) | Data yield |
| --- | ---: | ---: | ---: |
| Left/lowe | 2 | 7 | -0.0001 |
| Left/highe | 2 | 7 | -0.0017 |
| Right/highe | 1 | 2 | -0.0003 |
| Right/highe | 1 | 8 | -0.0001 |

The existing nonnegative SIMC-yield convention yields 36 zero SIMC cells.
All 135 recorded parent policies are `valid`; no accepted skip/zero parent
policy was needed. Negative/zero cells are observed producer results, not
grounds for invented normalization or an established physics defect.

## Preserved limits and ownership

This expands accepted baseline `Y0/MM_0` runtime coverage beyond the prior
[Left/lowe E.8.2 closure](e8-2-left-lowe-runtime-scientific-closure-2026-10-06.md)
without changing its historical scientific/visual scope or any other closed phase.
The accepted yield/support producers, page producers, frozen profile and
`testing/run_e8_2_canonical_five_baseline_yield_gate.py` remain unchanged.
No physics, yield, uncertainty, cut, normalization or production source changes.

- Private Method-A canonical-five F.4 scientific equivalence remains `BLOCKED`
  at `$.parents[0].absolute_baseline_parent_sum` (owner failure:
  `f6_3_scientific_equivalence_mismatch:f4:$.parents[0].absolute_baseline_parent_sum`).
  No original-owner retry, reviewed/current F.4 acceptance, private `YA/MM_A`
  comparison, numerical Method-A improvement or promotion follows.
- Absolute-SIMC amplitude/units remain independently `BLOCKED` by
  `SIMC_normfac_luminosity_and_charge_units_not_source_proven`. Matching
  numerical support proves no absolute-scale agreement, rate or luminosity.
- Final E.8 and F.6.4 remain `BLOCKED`. Method A stays private/nonproduction;
  Method B stays diagnostic/cross-check only and numerically excluded.
  E.8.4 unavailable-page records establish no private-branch availability.
- Earlier E.8.2/E.8.3/F.6.3/E.8.4 Left/lowe and other closed scopes remain intact.
  No cross-section accuracy, total-systematics, detector-response causality or
  physical explanation of negative/zero yields is accepted here.

[CURRENT](../CURRENT.md) owns the sole push-stable substantive next gate:
baseline-five extracted-yield scientific QA after actual-diff acceptance,
user publication and pushed-state synchronization. This closure authorizes
no execution of that successor gate.
