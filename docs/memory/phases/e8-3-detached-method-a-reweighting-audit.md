# E.8.3 — detached Method-A reweighting audit

## Status

`SOURCE REVIEWED` — independent ChatGPT actual-diff review passed the complete
cumulative E.8.3 + Fix.1 candidate. This record does not claim ROOT/PyROOT,
full-analysis, farm, or runtime validation.

## Starting identity and scope

- Starting branch/HEAD: `test` at
  `91bb7809d27d84d7709a6600dbb0dc9ab514a458`.
- E.8.3 adds a fail-closed procedure-PDF appendix reader/renderer only. It
  changes `src/cuts/full_background_subtraction_plots.py`, the existing
  `src/cuts/rand_sub.py` procedure route, focused E.8.3 test coverage, and the
  warranted active-state memory records.
- It creates no Method-A producer, event-factor table, checkpoint, independent
  PDF lifecycle, F.6.3 branch, or production-pion template.

## Frozen accepted authority

E.8.3 consumes only the following byte-pinned persisted detached evidence for
`Q4p4W2p74` and selects one exact canonical setting from `Left-lowe`,
`Left-highe`, `Center-lowe`, `Center-highe`, and `Right-highe`:

- F.4 JSON SHA-256
  `adcc01900b8adc211e296aaedaa314fdb86207d51483ebfd4b3edf69f558f188`;
  correction fingerprint
  `362241005c02f2149e260c391b5c3d35793287573128b42cf5ed693419d9d2f3`;
  artifact fingerprint
  `c4b9f513d5918ca77179bbe5fab28c5d9d14e63a440338545e73961fa50b9d67`.
- F.5 JSON SHA-256
  `143e3af6b1c69560e5bf351155b570d05f19e7e0ceb14665f5448c351ad501be`;
  propagation fingerprint
  `d11b728d1089301a12c29e7f8b1798c6e0b6021ac47bd5d2afc2b62b47a1effa`;
  artifact fingerprint
  `261968ee7d9590d7d0afe0cd15155ef745a4169f95f52ffffd392aadd320de63`.
- F.6.1 JSON SHA-256
  `62bad2dae0f65f1fff87eeffd071c29dbbd4cdc15a43f37d6d9853cb58b25bd6`;
  artifact fingerprint
  `377872a218a780481347402e4410c81bd568a4fb7682cd49f0f3352b499bad41`;
  validation fingerprint
  `d992d789b3434897d76691df27c190a0f51479f67c0d5b496e84835e517c77f8`.
- The existing F.6.2 reader remains the sole F.6.2 authority: JSON SHA-256
  `5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1`,
  artifact fingerprint
  `ee713b70de898bad8fa61164cbc8a54886af1ea1eb1e712ed95de17df8890cd0`,
  and validation fingerprint
  `7edc73fce20ad7dc7622b8c23a7ba7e8986e367ace605884e010945595370b3b`.

Missing, stale, malformed, nonfinite, noncanonical, or wrong-authority inputs
produce an explicit E.8.3 unavailable procedure page. They do not trigger a
fallback, reconstruction, production failure, or a replacement calculation.

## Presentation and ownership

The appendix states the detached comparison
`b_j^0 = s_j*w0_j` and `b_j^A = s_j*w0_j*C_j`. F.4 parent normalization remains
frozen and no canonical `(t,phi)` child is independently normalized.

For the selected setting, E.8.3 copies the accepted F.6.1 signed parent-level
`analysis_MM` baseline, Method-A, and delta arrays; it draws a ratio only at
finite nonzero persisted baseline bins. It copies F.5's persisted 3 x 9 signed
baseline/Method-A/delta children and persisted parent closure, preserving explicit
zero/empty children. Existing E.8/F.6.2 pages remain authoritative for the
L/B/A, acceptance, support/OOD, effective-statistics, and kaon-window context.

No clipping, cap, smoothing, interpolation, normalization, correction
recalculation, event-level factor persistence, yield construction, cross section,
Method-B numerical input, empirical Fit 1/Fit 2 activation, or F.6.3 work is
introduced. Baseline production, `w0`, cuts, binning, normalizations, templates,
SIMC, efficiencies, and yields remain unchanged.

## Runtime path and local checks

`rand_sub.py` locates the three deterministic accepted F-stage basenames only
from established `OUTPATH`, builds the E.8.3 payload, passes it into the existing
procedure renderer, and retains it in the existing E.8.2 pair-safe post-yield
rerender state. It does not create another PDF lifecycle.

Codex ran deterministic `py_compile`, focused E.8.3, procedure-page/runtime
contract, and consumed-artifact F.4/F.5/F.6.1/F.6.2 regression checks after
implementation. They do not establish ROOT/PyROOT or farm behavior; no farm
action was run.

## Fix.1 — gap-safe ratio and persisted empty-child identity

Fix.1 repairs only source-review findings in the local E.8.3 candidate. The
persisted F.5 `event_counts` matrix is now copied unchanged into each selected
parent display payload. The `(t,phi)` page prints `EMPTY` only when that copied
count is zero; a populated child with a zero signed aggregate remains
`POPULATED`. No aggregate, closure, weight, bin, or child normalization is
recomputed or changed.

The persisted F.6.1 ratio remains defined only for finite, nonzero baseline
denominators, but its ROOT representation is now points-only (`AP`). It cannot
draw a line across an omitted denominator bin. No undefined value, cap, clip,
smoothing, interpolation, or replacement scientific array is introduced.

The focused pure-Python module now explicitly covers all original
authority/fail-closed identities; missing/duplicate canonical F.4/F.5/F.6.1
inventory; malformed F.5/F.6.1 geometry; exact persisted arrays, counts,
closure, and all nine children; the gap-safe renderer; existing `OUTPATH` and
pair-safe procedure lifecycle; unavailable presentation; and the detached
non-mutation/Method-B/empirical-residual/yield/F.6.3 boundaries.

### Codex-run deterministic local checks

Interpreter: `Python 3.12.10` at
`C:\Users\trott\AppData\Local\Programs\Python\Python312\python.exe`.

All results below are **Codex-reported checks; NOT RUN by ChatGPT.**

```text
PASS  python -B -m py_compile src/cuts/full_background_subtraction_plots.py src/cuts/rand_sub.py testing/test_e8_3_detached_method_a_reweighting_audit.py
      exit 0
PASS  python -B -m unittest testing.test_e8_3_detached_method_a_reweighting_audit -v
      Ran 21 tests in 11.222s; OK
PASS  python -B -m unittest testing.test_full_background_subtraction_plots -v
      Ran 102 tests in 3.178s; OK (skipped=18: expected PyROOT-unavailable or retired-tail cases)
PASS  python -B -m unittest testing.test_pion_hgcer_phase_e_runtime_contract -v
      Ran 3 tests in 0.004s; OK
PASS  python -B -m unittest testing.test_pion_hgcer_method_a_parent_preserving_correction -v
      Ran 8 tests in 32.081s; OK
PASS  python -B -m unittest testing.test_pion_hgcer_method_a_tphi_propagation -v
      Ran 8 tests in 10.779s; OK
PASS  python -B -m unittest testing.test_pion_hgcer_method_a_reweighting_validation -v
      Ran 9 tests in 47.647s; OK
PASS  python -B -m unittest testing.test_pion_hgcer_method_a_acceptance_refinement_validation -v
      Ran 15 tests in 55.381s; OK
```

The first newly-added focused fake-ROOT run reported two test-double failures
because `_FakeHistogram` did not implement the renderer's pre-existing `Draw`
method. This was repaired in the test double only; the immediate rerun above
passed all 21 focused tests. An earlier combined F.4/F.5 terminal invocation
ended before its summary was captured, so it is not used as evidence; the
separate complete reruns above are the recorded results.

```text
PASS  python -B tools/update_memory_manifest.py --root . --write
      MANIFEST: WROTE docs/memory/manifest.json
PASS  python -B tools/update_memory_manifest.py --root . --check
      MANIFEST: PASS
PASS  python -B tools/check_memory_health.py --root .
      MEMORY HEALTH: PASS
PASS  python -B tools/memory_bootstrap.py --root . --json
      exit 0; reports test at 91bb7809d27d84d7709a6600dbb0dc9ab514a458 and the intentional dirty worktree
PASS  python -B -m unittest testing.test_memory_health -v
      Ran 35 tests in 4.142s; OK
PASS  git -c core.safecrlf=false diff --check
      exit 0; no whitespace errors
```

No pre-existing unrelated failure occurred. No ROOT/PyROOT, full `main.py`,
farm, or runtime validation was run or is claimed. `src/cuts/rand_sub.py` is
frozen for Fix.1 and remained byte-unchanged during this repair (SHA-256 before
the Fix.1 edit: `cc5c4b51acf472990536c39ebe2b79e1f47d2b8cd791b291600e0d63af94a739`).

## Independent source-review closure

Independent ChatGPT actual-diff review of the complete cumulative
`kaonlt_review(20260924-064150).diff` passed. It accepted the detached,
byte-pinned F.4/F.5/F.6.1/F.6.2 authority path; the Fix.1 points-only (`AP`)
gap-safe ratio; persisted F.5 `event_counts` empty-child identity; and the
focused contract coverage. No source or test changed after that review.

The Codex-reported checks above were **NOT RUN by ChatGPT**. This is source/diff
review only: no ROOT/PyROOT, full `main.py`, Jefferson Lab farm, procedure-PDF
farm rendering, runtime acceptance, Method-A production promotion, or F.6.3
implementation is claimed.

## Next

`NEXT` — user-controlled commit/push of the reviewed cumulative
post-E.8.2-reconciliation + E.8.3/Fix.1 source/test/memory set
-> ChatGPT pushed-state review
-> F.6.3 source/runtime-path audit and standalone implementation contract.

E.8.3 is `SOURCE REVIEWED`; F.6.3 is dependency `NEXT` but does not begin
until that pushed-state review passes. E.8.4 remains `BLOCKED` pending F.6.3,
final E.8 remains `BLOCKED` pending E.8.4 and its runtime/visual gate, and
F.6.4 remains `BLOCKED` pending full production-impact evidence.
