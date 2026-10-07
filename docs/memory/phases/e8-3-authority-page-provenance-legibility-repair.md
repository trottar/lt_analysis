# E.8.3 authority-page provenance-legibility repair

Status: `CLOSED / RUNTIME VALIDATED` only for the Q4p4W2p74 / Left / lowe
authority-page presentation repair and narrow E.8.3 runtime/presentation gate.
The implementation history and local checks below precede this acceptance.

## Starting identity and accepted continuity

Branch `test`; starting HEAD/local `origin/test`:
`0ac6134c6f76efa731f94390b986b9533b53a7f6`.
The [contract](e8-3-authority-page-provenance-legibility-repair-task-contract.md)
and [fresh supplied evidence](../evidence/e8-3-left-lowe-runtime-readiness-visual-blocker-2026-10-07.md)
own the narrow scope and accepted runtime/visual facts.

SOURCE VERIFIED: the farm `706ae708ae69be0585472982a4df3d87c63f186f` to current
source compare is documentation-only (seven memory paths). Original E.8.3
`c9ed0d6b4d0013f7475eedbfca57a096730ea840` reader/validator/payload/ratio/histogram/
unavailable functions retain source continuity; `rand_sub.py` and the starting
focused test share the original blobs. Current rendering already labels
historical accepted F.6.1 distinctly from current F.6.3/E.8.4.

RUNTIME VERIFIED by supplied artifacts/prior ChatGPT review, not a Codex farm
run: all eight E.8.3 pages 65–72 are present with no renderer failures; pages
66–72 are legible. Page 65 clips the F.5/F.6.2 SHA-256 values because it puts
two full hashes on each long line. That is the sole fresh visual blocker.

## Presentation-only implementation

Adjacent pure `_e8_3_authority_lines(payload)` formats four input SHA-256 values
and three accepted fingerprints without shortening any value. Separate
SHA-256/fingerprint headings and short owner labels place each complete value
on its own line. All lines are at most 90 characters. No calculation or
mutation is performed.

`_e8_3_render_authority_page()` passes that tuple into the same text renderer.
Canvas/name/title, text size 0.023 and single-page inventory remain unchanged.
The seven already-legible E.8.3 pages, page IDs, schemas, payload/validation
logic, authority constants, ratio definition, child identity, parent closure,
OUTPATH route and pair-safe PDF/manifest lifecycle are unchanged.

Four new deterministic tests retain all 21 existing E.8.3 tests. They check
all seven full values exactly once, bounded lines with no combined hashes,
lineage/ownership/setting statements, the exact one-page renderer call and
payload immutability. The bounded-line test fails against the original
renderer's two-hash construction; after repair all 25 focused tests pass.

## Exact changed paths

- `src/cuts/full_background_subtraction_plots.py`
- `testing/test_e8_3_detached_method_a_reweighting_audit.py`
- `docs/memory/CURRENT.md`
- `docs/memory/manifest.json`
- `docs/memory/evidence/e8-3-left-lowe-runtime-readiness-visual-blocker-2026-10-07.md`
- `docs/memory/phases/e8-3-authority-page-provenance-legibility-repair-task-contract.md`
- `docs/memory/phases/e8-3-authority-page-provenance-legibility-repair.md`

All other tracked starting files, including four unrelated stat-only memory
files, are byte-preserved. All existing module nodes/functions except the
authority renderer are unchanged; only one adjacent pure helper is added.
`rand_sub.py`, owner/profile/collector and production remain byte-preserved.

## Codex-run deterministic local checks

WSL Ubuntu, Linux 5.15.167.4-microsoft-standard-WSL2; /usr/bin/python3 3.12.3.
System Python lacked test dependencies. NumPy 2.5.3, SciPy 1.18.1 and matplotlib
3.11.2 were installed only into task-specific temporary WSL storage; no
repository dependency or system-package configuration changed. Initial import
failures preceded those dependency additions. An initial command-transport
string error was corrected before final validation; it changed no scientific
logic.

All results below are Codex-reported checks; NOT RUN by ChatGPT.

```sh
PYTHONPYCACHEPREFIX=/tmp/e8-3-provenance-pycache python3 -B -m py_compile src/cuts/full_background_subtraction_plots.py testing/test_e8_3_detached_method_a_reweighting_audit.py
# PASS; compile cache outside repository
PYTHONPATH=/tmp/e8-3-provenance-deps python3 -B -m unittest testing.test_e8_3_detached_method_a_reweighting_audit -v
# PASS; 25 tests, no skips
PYTHONPATH=/tmp/e8-3-provenance-deps MPLCONFIGDIR=/tmp/e8-3-provenance-matplotlib python3 -B -m unittest testing.test_full_background_subtraction_plots -v
# PASS; 102 tests, 18 retained presentation-policy/PyROOT-unavailable skips
PYTHONPATH=/tmp/e8-3-provenance-deps python3 -B -m unittest testing.test_e8_2_baseline_stage_audit -v
# PASS; 26 tests, no skips
PYTHONPATH=/tmp/e8-3-provenance-deps python3 -B -m unittest testing.test_pion_hgcer_phase_e_runtime_contract -v
# PASS; 3 tests, no skips
python3 -B tools/update_memory_manifest.py --root . --write
python3 -B tools/update_memory_manifest.py --root . --check
python3 -B tools/check_memory_health.py --root .
python3 -B tools/memory_bootstrap.py --root . --json
git diff --check
```

Manifest write/check, ordinary memory health, bootstrap JSON, sole NEXT,
function/byte/index preservation and git diff --check: PASS, no health warnings.
CURRENT/MEMORY/CURRENT_HANDOFF sizes: 8190 / 20085 / 323 bytes.
These deterministic checks established SOURCE VERIFIED behavior only; physical
PDF legibility was then NOT VERIFIED pending fresh farm validation.

## Accepted later farm closure — 2026-10-07

RUNTIME VERIFIED by supplied artifacts/prior independent ChatGPT review at
farm source `63eb2bf66d98baaffa39e88cf1db7ac37397d364`; see the
[closure evidence](../evidence/e8-3-left-lowe-runtime-closure-2026-10-07.md).
Page 65 now passes full-provenance legibility, with four complete input hashes
and three complete fingerprints visible without clipping or overlap. Pages
66–72 are pixel-identical to the prior accepted render; page inventory and
manifest are unchanged, with no renderer failures. Owner/provenance/artifact/
ZIP and final preservation checks passed. Codex records this acceptance without
rerunning or independently validating external farm artifacts.
This closure changes documentation only; no scientific or production source
changed. The accepted repair changed only authority-page presentation and
does not establish whole-analysis JSON byte identity.

## Boundary and successor

Method A stays detached/non-production; Method B diagnostic/cross-check only
and numerically excluded. No payload/scientific redesign, owner/profile/
collector change, production mutation, ROOT/PyROOT, analysis, farm execution,
staging, commit, push or ref update occurs. Canonical-five provenance repair
remains `DEFERRED`; final E.8/F.6.4 remain `BLOCKED`.

The implementation's then-required successor was: after ChatGPT actual-diff review, user commit/push and pushed-state
synchronization/farm-readiness review, run one fresh isolated Q4p4W2p74 / Left /
lowe procedure owner gate. Review E.8.3 pages 65–72, all four full hashes and
three full fingerprints, no clipping/overlap and unchanged later-page values
before any narrow E.8.3 runtime closure. That gate is now consumed by the
accepted later farm evidence. CURRENT owns the sole ordinary NEXT: user decision
on whether to explicitly reopen deferred canonical-five provenance/identity
repair. No further automatic Left/lowe E.8.3 work or farm command is authorized.
