# KaonLT E.8.3 Left/lowe authority-page provenance-legibility repair — task contract

## 1. Task class and authority

This is a **narrow presentation-only source repair** for the existing detached
E.8.3 Method-A reweighting audit.

It is triggered by a fresh farm procedure PDF that was generated successfully
at source:

```text
706ae708ae69be0585472982a4df3d87c63f186f
```

and independently inspected after the E.8.2 closure. The E.8.3 runtime route,
payload construction, page inventory, and seven substantive comparison pages
rendered, but the E.8.3 authority/provenance page has clipped SHA-256 text.

The current pushed `test` branch is one documentation-only closure commit later:

```text
0ac6134c6f76efa731f94390b986b9533b53a7f6
Record narrow E8.2 Left lowe runtime scientific closure
```

The `706ae708... -> 0ac6134...` compare changes only seven
`docs/memory/` paths. Therefore the farm-rendered E.8.3 source/runtime evidence
applies directly to the current E.8.3 scientific/presentation source.

This task must not redesign E.8.3, Method A, F.6.3, or any production path.

Before editing, read repository-root `AGENTS.md`, then the startup core in this
exact order:

1. `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read only task-relevant records:

- `docs/memory/MAINTENANCE.md`
- `docs/memory/CODEX.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- `docs/memory/phases/e8-3-detached-method-a-reweighting-audit.md`
- `docs/memory/phases/e8-3-detached-method-a-reweighting-audit-fix1-task-contract.md`
- `docs/memory/evidence/e8-2-left-lowe-runtime-scientific-closure-2026-10-06.md`

Authority remains:

```text
current source/diff
-> newest applicable farm/runtime evidence
-> newest applicable validation/handoff artifact
-> tracked durable memory
-> older chat/history
```

---

## 2. Exact starting identity and startup gate

Repository:

```text
https://github.com/trottar/lt_analysis/tree/test
```

Required branch:

```text
test
```

Required starting HEAD and local `origin/test`:

```text
0ac6134c6f76efa731f94390b986b9533b53a7f6
```

Before editing, establish in the normal **WSL/POSIX shell**:

```bash
git branch --show-current
git rev-parse HEAD
git rev-parse refs/remotes/origin/test
git status --short --untracked-files=all
git diff --check
```

Do not use PowerShell commands for this workstation workflow.

The only expected new task file before implementation is:

```text
docs/memory/phases/e8-3-authority-page-provenance-legibility-repair-task-contract.md
```

Temporary root-level `kaonlt_review*.diff` files are permitted only as review
artifacts.

If branch, HEAD, local `origin/test`, or unrelated workstation state differs,
stop as `BLOCKED`. Do not reset, stash, clean, restore, or overwrite unrelated
user state.

If the previously observed `/mnt/c` stat-only entries recur, do not treat
porcelain alone as content evidence. Reuse the established procedure: prove
HEAD = index = worktree with Git object hashes before classifying them. Do not
modify content-identical unrelated files.

No farm command is authorized in this implementation task.

---

## 3. Accepted fresh E.8.3 runtime-readiness evidence

The fresh procedure PDF belongs to package stem:

```text
KaonLT_E8_2_Left_lowe_scientific_audit_20261006-233426
```

Farm source:

```text
706ae708ae69be0585472982a4df3d87c63f186f
```

Scope:

```text
Q4p4W2p74 / Left / lowe
```

Procedure PDF:

```text
Left_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction.pdf
SHA-256:
aacac5daa935eb1f2af6dacdf0ddfdcdd4f5b202e0d54fc11b8edd254576e61e
```

Page manifest:

```text
Left_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction-manifest.json
SHA-256:
e51fe03cff97ce00e05c850b031206e97a9cb9d15b0c71279618fda1a9577f02
```

The existing accepted E.8.2 closure evidence owns the complete package,
owner/provenance/artifact, and companion identities.

RUNTIME VERIFIED from the fresh supplied artifacts and ChatGPT inspection:

- the procedure PDF has 74 pages;
- page-manifest `renderer_failures` is empty;
- setting identity is exactly Left / lowe / Q4p4W2p74 / kaon;
- all eight E.8.3 page records are present:
  - page 65 `full_background.e8_3.authority`;
  - page 66 `full_background.e8_3.mm`, t1;
  - page 67 `full_background.e8_3.tphi`, t1;
  - page 68 `full_background.e8_3.mm`, t2;
  - page 69 `full_background.e8_3.tphi`, t2;
  - page 70 `full_background.e8_3.mm`, t3;
  - page 71 `full_background.e8_3.tphi`, t3;
  - page 72 `full_background.e8_3.f6_2_cross_reference`;
- pages 66–72 are visually complete and legible;
- all three MM pages show persisted baseline, Method-A, delta, and gap-safe
  defined-bin-only ratio points;
- all three `(t,phi)` pages show all nine children, persisted EMPTY/POPULATED
  identity, baseline/Method-A/delta values, and persisted parent closure;
- the F.6.2 cross-reference page is legible;
- no E.8.3 renderer failure is present.

This evidence establishes that the existing E.8.3 runtime route is active in
the real farm procedure path. It does **not** close E.8.3 because page 65 fails
visual provenance legibility.

---

## 4. First failed invariant: authority-page SHA-256 clipping

The sole observed E.8.3 visual blocker is page 65:

```text
full_background.e8_3.authority
```

The page is otherwise legible, including:

- historical accepted F.6.1 lineage label;
- detached/non-production statement;
- Method-B exclusion;
- dormant empirical residual statement;
- F.6.3/E.8.4/F.6.4 ownership text;
- setting identity;
- the three fingerprint lines.

But these two current lines place two 64-character SHA-256 values on one line:

```text
Accepted input SHA-256: F.4=<64 hex> F.5=<64 hex>
F.6.1=<64 hex> F.6.2=<64 hex>
```

The physical A4 farm PDF clips the right-hand hashes:

- the F.5 SHA-256 is visibly truncated;
- the F.6.2 SHA-256 is visibly truncated.

PDF text extraction confirms the same right-edge loss.

SOURCE VERIFIED cause in:

```text
src/cuts/full_background_subtraction_plots.py
```

`_e8_3_render_authority_page()` passes those long strings directly into the
non-wrapping `_e8_text_page()` / `_e8_add_text()` path.

This is a presentation-only defect. The persisted authority values themselves
are correct and already pass source/runtime authority checks.

Therefore the current E.8.3 readiness state is:

```text
E.8.3:
SOURCE REVIEWED

E.8.3 Left/lowe runtime-readiness gate:
BLOCKED at rendered authority-page provenance legibility
```

Do not reinterpret or recompute the E.8.3 scientific arrays.

---

## 5. Source continuity already established

The original E.8.3 implementation commit is:

```text
c9ed0d6b4d0013f7475eedbfca57a096730ea840
Add E8.3 detached Method-A reweighting audit
```

Current-source audit established:

- `src/cuts/rand_sub.py` has the **same Git blob** at the original E.8.3 commit
  and current `0ac6134...`:
  `375a6f167962f3ab0e5bd0fe3731f8c5424e1d07`;
- `testing/test_e8_3_detached_method_a_reweighting_audit.py` has the **same Git
  blob** at both:
  `204b718234d60a2d85e7ca49d8a84573d3d0b320`;
- the E.8.3 reader/validator/payload functions in
  `full_background_subtraction_plots.py` are byte-identical between the
  original E.8.3 implementation and current source:
  - `_e8_3_unavailable`;
  - `_e8_3_read_json`;
  - `_e8_3_hash`;
  - `_e8_3_setting_id`;
  - `_e8_3_matrix`;
  - `_e8_3_validate_f4`;
  - `_e8_3_validate_f5`;
  - `_e8_3_validate_f6_1`;
  - `build_full_background_subtraction_e8_3_payload`;
  - `_e8_3_ratio_points`;
  - `_e8_3_signed_histogram`;
  - `_render_full_background_subtraction_e8_3_unavailable_page`.

Later changes to the four E.8.3 page-rendering functions only clarified the
historical F.6.1 versus current F.6.3/E.8.4 lineage in displayed text. They did
not change the persisted arrays, authority checks, ratio definition, child
identity, parent closure, or runtime input path.

Preserve this continuity.

---

## 6. Frozen scientific ownership

E.8.3 remains detached and presentation-only.

It consumes only byte-pinned accepted authorities:

```text
F.4 JSON SHA-256:
adcc01900b8adc211e296aaedaa314fdb86207d51483ebfd4b3edf69f558f188

F.5 JSON SHA-256:
143e3af6b1c69560e5bf351155b570d05f19e7e0ceb14665f5448c351ad501be

F.6.1 JSON SHA-256:
62bad2dae0f65f1fff87eeffd071c29dbbd4cdc15a43f37d6d9853cb58b25bd6

F.6.2 JSON SHA-256:
5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1
```

Freeze:

- all E.8.3 input SHA/fingerprint constants;
- all F.4/F.5/F.6.1/F.6.2 validation;
- `b_j^0 = s_j * w0_j`;
- `b_j^A = s_j * w0_j * C_j`;
- parent-preserving normalization;
- persisted F.5 `event_counts`;
- persisted F.5 `3 x 9` child values;
- persisted F.6.1 MM arrays;
- ratio definition and points-only `AP` rendering;
- all eight E.8.3 page IDs;
- `E8_3_PRESENTATION_SCHEMA_VERSION`;
- existing `OUTPATH` input route;
- existing E.8.2 pair-safe PDF/manifest lifecycle;
- Method B numerical exclusion;
- dormant empirical residual Fit 1/Fit 2;
- baseline production;
- F.6.3/E.8.4 ownership;
- F.6.4 as the only promotion decision.

Do not change any pion/proton/random/dummy subtraction, weights, cuts,
normalization, binning, templates, SIMC, yields, efficiencies, acceptance,
L/T separation, or cross sections.

---

## 7. Exact allowed versioned paths

Only these seven paths may change:

```text
src/cuts/full_background_subtraction_plots.py
testing/test_e8_3_detached_method_a_reweighting_audit.py
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/evidence/e8-3-left-lowe-runtime-readiness-visual-blocker-2026-10-07.md
docs/memory/phases/e8-3-authority-page-provenance-legibility-repair-task-contract.md
docs/memory/phases/e8-3-authority-page-provenance-legibility-repair.md
```

No other versioned path may change.

In particular freeze:

```text
src/cuts/rand_sub.py
testing/run_e8_2_left_lowe_scientific_audit_gate.py
testing/collect_pion_hgcer_validation_bundle.py
testing/pion_hgcer_validation_bundle_profile_e8_2_left_lowe.json
run_Prod_Analysis.sh
set_SymLinks.sh
```

and all other `src/`, `testing/`, `tools/`, `farm_env/`,
`background_samples/`, and unrelated memory paths.

---

## 8. Exact repair objective

Repair only the E.8.3 authority-page provenance text layout.

Preferred narrow implementation:

1. introduce a small pure deterministic helper, e.g.
   `_e8_3_authority_lines(payload)`, or equivalent narrowly local formatting;
2. keep each full 64-character hash/fingerprint intact;
3. display one authority value per bounded line rather than combining two
   SHA-256 values on one line;
4. preserve the existing historical-lineage and ownership statements;
5. keep the same authority page title and page ID;
6. keep the same single-page inventory; do not add continuation pages.

A suitable displayed structure is equivalent to:

```text
Accepted input SHA-256:
F.4:   <full 64 hex>
F.5:   <full 64 hex>
F.6.1: <full 64 hex>
F.6.2: <full 64 hex>

Accepted fingerprints:
F.4 correction:   <full 64 hex>
F.5 propagation:  <full 64 hex>
F.6.1 validation: <full 64 hex>
```

The exact typography may differ, but every complete value must fit visibly
inside the page.

Do not abbreviate, truncate, wrap by dropping characters, or replace a hash with
a shortened prefix.

If modestly reducing the text size is needed to accommodate the additional
lines, keep it readable and scoped only to this authority page.

Do not modify the seven already-legible E.8.3 pages.

---

## 9. Deterministic regression coverage

Extend only:

```text
testing/test_e8_3_detached_method_a_reweighting_audit.py
```

No PyROOT dependency is permitted.

At minimum add tests proving:

### 9.1 Full provenance completeness

For a valid payload, capture the authority-page line tuple before rendering and
assert that:

- full F.4 SHA appears exactly once;
- full F.5 SHA appears exactly once;
- full F.6.1 SHA appears exactly once;
- full F.6.2 SHA appears exactly once;
- full F.4 correction fingerprint appears exactly once;
- full F.5 propagation fingerprint appears exactly once;
- full F.6.1 validation fingerprint appears exactly once.

### 9.2 Bounded-line policy

Assert no authority-page line contains two distinct 64-character provenance
values.

The test must fail on the current two-hashes-per-line implementation.

Do not use a test that accepts truncated prefixes.

### 9.3 Semantic preservation

Assert the authority-page lines still visibly state:

- historical accepted F.6.1 lineage;
- not current F.6.3/E.8.4 lineage;
- detached presentation only;
- Method B has no numerical input;
- empirical residual Fit 1/Fit 2 are inactive;
- F.6.3 owns the private parallel branch;
- Method A is not production-promoted;
- F.6.4 remains the promotion decision;
- current setting identity.

### 9.4 Renderer contract

Patch/mock `_e8_text_page` and prove `_e8_3_render_authority_page()`:

- uses the existing authority-page canvas/name/title identity;
- passes the complete bounded lines;
- performs no calculation or mutation;
- does not create a second page.

Retain all existing 21 E.8.3 tests.

---

## 10. Source immutability gate

Within `src/cuts/full_background_subtraction_plots.py`, only these may change:

```text
_e8_3_render_authority_page
```

plus at most one new pure E.8.3 authority-line formatting helper located
adjacent to it.

All other functions in the file must remain byte-identical in function body
where practical to verify, especially:

```text
_e8_3_validate_f4
_e8_3_validate_f5
_e8_3_validate_f6_1
build_full_background_subtraction_e8_3_payload
_e8_3_ratio_points
_e8_3_signed_histogram
_e8_3_render_mm_page
_e8_3_render_tphi_page
_render_full_background_subtraction_e8_3_pages
_render_full_background_subtraction_e8_3_unavailable_page
```

`src/cuts/rand_sub.py` must remain byte-identical.

No page IDs or schema constants may change.

---

## 11. Required evidence record

Create:

```text
docs/memory/evidence/e8-3-left-lowe-runtime-readiness-visual-blocker-2026-10-07.md
```

Record:

- farm source `706ae708...`;
- current source `0ac6134...`;
- the compare between them is documentation-only;
- procedure PDF and manifest identities;
- all eight E.8.3 page IDs and page numbers 65–72;
- `renderer_failures=[]`;
- pages 66–72 visually accepted;
- page 65 visual blocker exactly as described above;
- E.8.3 source/runtime route is exercised on the farm;
- E.8.3 remains `SOURCE REVIEWED`, not `CLOSED / RUNTIME VALIDATED`;
- runtime-readiness is `BLOCKED` only at authority-page provenance legibility;
- no Method-A promotion or production correctness follows.

Do not record new scientific magnitude conclusions in this repair task.

---

## 12. Required phase record and CURRENT update

Create:

```text
docs/memory/phases/e8-3-authority-page-provenance-legibility-repair.md
```

Record:

- starting identity;
- source/runtime continuity audit;
- fresh farm page inventory;
- first failed invariant;
- exact presentation-only repair;
- changed paths;
- local checks;
- no farm execution;
- fresh PDF validation still pending.

Update `docs/memory/CURRENT.md` compactly.

Required meaning after successful local implementation:

```text
E.8 remains ACTIVE.
E.8.2 Left/lowe remains CLOSED / RUNTIME VALIDATED.
E.8.3 remains SOURCE REVIEWED.

Fresh E.8.3 farm evidence at source 706ae708... exercises the real Left/lowe
runtime path and produces all eight E.8.3 pages with no renderer failures.

Pages 66–72 pass visual review.
Page 65 is BLOCKED only by clipped provenance SHA-256 text.

The narrow authority-page provenance-legibility repair is:
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING.

No E.8.3 scientific redesign is justified.
Canonical-five provenance/identity repair remains DEFERRED.
Final E.8 and F.6.4 remain BLOCKED.
```

The sole ordinary NEXT after successful implementation must be equivalent to:

```text
NEXT — after ChatGPT actual-diff review, user commit/push and pushed-state
synchronization/farm-readiness review, run exactly one fresh isolated
Q4p4W2p74 / Left / lowe procedure owner gate and visually inspect E.8.3 pages
65–72, especially the full provenance values on the authority page, before any
E.8.3 runtime closure.
```

Do not make the commit/push itself the scientific NEXT.

Keep CURRENT within the tracked size policy by consolidating consumed wording
rather than appending a long block.

Regenerate/check `docs/memory/manifest.json`.

---

## 13. Deterministic local checks

Discover the repository-appropriate interpreter according to tracked
`docs/memory/TOOLS.md`.

Do not run ROOT/PyROOT, production analysis, or the farm.

Run at minimum:

```text
py_compile:
  src/cuts/full_background_subtraction_plots.py
  testing/test_e8_3_detached_method_a_reweighting_audit.py

testing.test_e8_3_detached_method_a_reweighting_audit
testing.test_full_background_subtraction_plots
testing.test_e8_2_baseline_stage_audit
testing.test_pion_hgcer_phase_e_runtime_contract

manifest write/check
ordinary memory health
bootstrap JSON
git diff --check
```

Report exact commands, interpreter/OS, counts, skips, and results.

Fake-ROOT/local checks are `SOURCE VERIFIED` only. Physical PDF legibility
remains farm-only.

---

## 14. Farm-validation boundary

No farm run occurs in Codex implementation.

After:

```text
implementation
-> deterministic local checks
-> ChatGPT actual-diff review
-> user commit/push
-> ChatGPT pushed-state synchronization
-> farm-readiness review
```

the user will run exactly one fresh isolated Q4p4W2p74 / Left / lowe procedure
owner gate using the existing reviewed owner/collector chain.

Do not change the owner/profile/collector merely to validate this presentation
repair.

Fresh E.8.3 visual acceptance requires:

- the same eight page IDs;
- no renderer failures;
- all four complete input SHA-256 values visible on page 65;
- all three complete fingerprint values visible on page 65;
- no clipping or overlap;
- pages 66–72 remain legible;
- scientific arrays and child/parent values unchanged.

Only after this gate may E.8.3 be considered for narrow Left/lowe runtime
closure.

---

## 15. Final integrity and review bundle

At the end:

1. regenerate/check the memory manifest;
2. run ordinary memory health;
3. run bootstrap JSON;
4. run `git diff --check`;
5. confirm exactly the seven allowed versioned paths changed;
6. confirm frozen source/runtime paths remain byte-identical;
7. create one complete root-level `kaonlt_review.diff` containing:
   - tracked diffs;
   - complete no-index additions for each new contract/evidence/phase file.

Do not stage merely to create the review bundle.

Report:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
warning classification: <none or blocking/nonblocking>
manifest check: PASS | FAIL
bootstrap: PASS | FAIL
git diff --check: PASS | FAIL
```

---

## 16. Hard stops

Stop as `BLOCKED` rather than broadening scope if any of the following appears
necessary:

- changing E.8.3 authority values;
- changing E.8.3 payload or validation logic;
- changing ratio definition;
- changing child identity or parent closure;
- changing page IDs or schema;
- adding continuation pages;
- changing `rand_sub.py`;
- changing owner/profile/collector;
- changing E.8.2 or E.8.4 scientific/presentation behavior;
- changing F.6.3/F.6.4 ownership;
- activating Method A in production;
- activating Method B numerically;
- activating empirical residual fits;
- changing SIMC/normalization/yields;
- running ROOT/PyROOT locally;
- running the farm.

---

## 17. Acceptance state and Codex stop point

After successful local implementation:

```text
E.8.3:
SOURCE REVIEWED

E.8.3 authority-page repair:
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

Codex stops after:

- narrow renderer implementation;
- deterministic local tests;
- evidence/phase/CURRENT update;
- manifest/health/bootstrap checks;
- complete `kaonlt_review.diff`.

Codex must not:

- stage merely for review;
- commit;
- push;
- update refs;
- run ROOT/PyROOT;
- run production analysis;
- run the farm;
- begin F.6.3/F.6.4 work.

Return:

- starting/ending branch;
- starting/ending HEAD and local `origin/test`;
- exact changed paths;
- concise repair summary;
- exact local checks/results;
- frozen-path preservation result;
- memory-health report;
- `kaonlt_review.diff` path/size;
- explicit confirmation that no ROOT/PyROOT, production analysis, farm command,
  commit, push, or ref update occurred.
