# KaonLT — E.8.4 Fix.5.4 post-farm identity-audit memory checkpoint

## 1. Purpose

Record the direct `Q4p4W2p74 / Left / lowe` Fix.5 farm-render findings before any further scientific-source change.

This is a **memory-only checkpoint**. It must not modify analysis, tests, farm drivers, profiles, collectors, plotting source, Method-A mathematics, SIMC normalization, yield extraction, production physics, or accepted authorities.

The checkpoint establishes the exact next development sequence:

1. this memory checkpoint -> ChatGPT PASS -> user commit/push -> pushed-state review;
2. separate Fix.5.4 numerical audit/fix contract;
3. Codex numerical implementation;
4. ChatGPT actual-diff review;
5. user commit/push;
6. pushed-state review;
7. only then write/run the visualization-only contract from that reviewed pushed numerical source, after numerical closure;
8. Codex visualization implementation;
9. ChatGPT actual-diff review;
10. user commit/push;
11. pushed-state review;
12. one narrow Q4p4W2p74 / Left / lowe farm run;
13. fresh scientific and visual evidence review.

---

## 2. Exact starting state

Required branch:

```text
test
```

Required starting HEAD:

```text
f9d70732290ea461096374ca1270b47452644991
```

Commit subject:

```text
E8.4 Fix.5: complete Left/lowe shareable-page farm gate
```

Before changing files:

```bash
git branch --show-current
git rev-parse HEAD
git status --short --untracked-files=all
```

Requirements:

- branch must be `test`;
- HEAD must equal the required starting HEAD;
- inspect and preserve the existing worktree;
- root `AGENTS.md` and `.codex/` are local-only/untracked and must remain so;
- do not reset, clean, stash, commit, push, or run the farm;
- unrelated tracked changes are a blocker.

---

## 3. Required startup reading

Read in this order:

1. root `AGENTS.md`;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read:

- `docs/memory/LEARNINGS.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/e8-4-fix5-shareable-method-a-impact-pages.md`
- `docs/memory/evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md`
- `docs/memory/decisions/e8-full-analysis-procedure-roadmap.md`
- only source excerpts required to state the observed lineage distinction accurately; no source modifications are allowed.

---

## 4. Direct new farm evidence to record

A fresh tracked Fix.5 Left/lowe analysis was run from pushed source:

```text
f9d70732290ea461096374ca1270b47452644991
```

The tracked owner did not return its final ZIP, so this is not accepted bundle closure. However, the expensive full-analysis/render stage completed and produced fresh artifacts in the canonical OUTDIR:

```text
Left_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction.pdf
Left_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction-manifest.json
```

Observed farm properties:

- PDF timestamp: `2026-10-01 23:34` local farm-visible time;
- PDF size approximately `2.2 MB`;
- manifest timestamp: `2026-10-01 23:34`;
- `page_count = 97`;
- `renderer_failures = []`;
- all required Fix.5 E.8.4 pages were present:
  - `full_background.e8_4.method_a_vs_simc.t1/t2/t3`
  - `full_background.e8_4.baseline_method_a_simc.t1/t2/t3`
  - `full_background.e8_4.yield_summary.t1/t2/t3`
  - `full_background.e8_4.parent_closure`
- prior E.8.4 pages were retained.

The user copied inspection-only copies to Globus:

```text
KaonLT_E8_4_Fix5_Left_lowe_20261001-234439.pdf
KaonLT_E8_4_Fix5_Left_lowe_20261001-234439-manifest.json
```

The absence of the owner ZIP is an operational post-render failure still to be diagnosed separately. Do not require an analysis rerun merely to explain that packaging failure.

---

## 5. Confirmed scientific/presentation finding: lineage mismatch

Record as **CONFIRMED**:

The sequential E.8.3 and E.8.4 pages do not describe one identical Method-A lineage.

E.8.3 explicitly presents the historical accepted persisted F.6.1 aggregate lineage. Its page records the historical accepted input identities, including:

```text
F.4:
adcc01900b8adc211e296aaedaa314fdb86207d51483ebfd4b3edf69f558f188

F.5:
143e3af6b1c69560e5bf351155b570d05f19e7e0ceb14665f5448c351ad501be

F.6.1:
62bad2dae0f65f1fff87eeffd071c29dbbd4cdc15a43f37d6d9853cb58b25bd6
```

By contrast, the active F.6.3 branch consumes the separately validated current-baseline candidate lineage. Its candidate F.4 raw input is:

```text
1d545924eba89c7f9ffa28028e307aca9b434a89beec06863cf2893887b6b902
```

Therefore E.8.3 must not be presented or interpreted as the aggregate explanation of the current F.6.3/E.8.4 branch unless exact lineage identity is proved.

This does not invalidate the historical E.8.3 `SOURCE REVIEWED` work or the accepted historical F.6.1 authority. It is a cross-page provenance/presentation mismatch exposed by the current Fix.5 review.

---

## 6. Unresolved numerical identities — scientific interpretation BLOCKED

Record each item separately and do not convert a missing check into a claimed failure.

### 6.1 Pion-template to final-MM algebra

For every populated canonical child, the current branch must prove bin by bin:

```text
MM_A - MM_0 = -(B_pi_A - B_pi_0)
```

apart from defined floating-point tolerance.

The current presentation plots both sides but does not fail closed on their numerical equality.

Status: **NOT YET VALIDATED**.

### 6.2 Histogram to scalar-yield closure

For every populated canonical child, prove:

```text
Integral(MM_0) = Y0
Integral(MM_A) = YA
Integral(MM_A - MM_0) = YA - Y0
```

using the exact branch histograms and the exact producer-owned integration semantics.

The fresh `t1` pages show visually near-overlapping baseline/Method-A MM spectra while stored `Y0/YA` values differ substantially in several populated children. This is an observed inconsistency requiring numerical audit, not yet proof of which producer or presentation object is wrong.

Representative displayed values include:

```text
t1 phi [-180,-140):
Y0 = 0.033906
YA = 0.025165

t1 phi [140,180):
Y0 = 0.029873
YA = 0.016411
```

Status: **NOT YET VALIDATED**.

### 6.3 Signed-cancellation diagnostics

Because the relevant spectra contain signed random/dummy/background-subtracted content, each child audit must report separately:

- signed integral;
- positive-bin support;
- negative-bin support;
- absolute support.

This is required to distinguish a genuine object mismatch from a large fractional change caused by cancellation of positive and negative content.

Do not infer cancellation without these numbers.

### 6.4 SIMC source/normalization/unit identity

The Fix.5 renderer clones existing per-child SIMC MM support from:

```python
hist["_xsect_support_simc"]["mm"]
```

without renderer-side renormalization.

Nevertheless, before interpreting data/SIMC amplitude disagreement, prove that the plotted SIMC child is exactly the authoritative same-cell yield-chain SIMC object in the same normalization and units as `MM_0/MM_A`.

Record its relevant integration-window integral and provenance.

The fresh plots visually show a large data/SIMC amplitude difference in some populated cells. This observation alone does not establish a physics discrepancy until normalization identity closes.

Status: **NOT YET VALIDATED**.

---

## 7. Visualization finding

Record as **CONFIRMED presentation deficiency**:

In many baseline / Method-A / SIMC panels, the baseline blue curve is almost completely obscured by the magenta Method-A curve.

The plot therefore makes near-equality or small differences difficult to evaluate visually.

This must be repaired only **AFTER** the numerical identity audit closes and the numerical implementation passes ChatGPT actual-diff review, user commit/push and pushed-state review; only then write/run its separate visualization-only contract from that reviewed pushed source.

Future visualization work may change only presentation properties such as:

- draw order;
- line style;
- line width;
- marker style;
- marker frequency/size;
- legend clarity;
- difference or ratio support panel only if it consumes already-validated stored objects without numerical transformation beyond explicitly presentation-owned arithmetic.

It must not change histogram contents, normalization, cuts, binning, fit results, Method-A factors, SIMC scaling, yield extraction, or production objects merely to improve visibility.

---

## 8. Status discipline

Preserve existing closed upstream phases unless the new evidence directly invalidates what they claimed.

Specifically:

- F.4.Refresh.2 remains `CLOSED / RUNTIME VALIDATED`.
- F.6.3 current-baseline Left/lowe runtime mechanics remain `CLOSED / RUNTIME VALIDATED` for what the accepted evidence established:
  - branch execution;
  - live-cache parity;
  - real child changes;
  - signed parent preservation.
- Historical F.6.1/F.6.2 closures remain unchanged.
- E.8.3 remains `SOURCE REVIEWED` as the historical accepted-lineage detached presentation; it is not the current F.6.3 lineage explanation.
- Method A remains detached/non-production.
- Method B remains diagnostic-only and numerically absent.
- Final E.8 remains `BLOCKED`.
- F.6.4 remains `BLOCKED`.

Create the new active work item:

```text
E.8.4 Fix.5.4 — current-lineage numerical identity and SIMC normalization audit
```

Status:

```text
ACTIVE
```

The scientific interpretation of the fresh Fix.5 Method-A/SIMC pages remains `BLOCKED` pending Fix.5.4.

The visualization-improvement stage follows Fix.5.4 numerical closure, ChatGPT actual-diff review, user commit/push and pushed-state review; only then write/run the visualization-only contract from that reviewed pushed numerical source. It remains dependency-blocked until numerical identity closes.

No further farm run is authorized by this checkpoint.

---

## 9. Required CURRENT.md state

`CURRENT.md` must become materially accurate for pushed HEAD:

```text
f9d70732290ea461096374ca1270b47452644991
```

and the fresh farm observation.

It must state:

```text
ACTIVE:
E.8.4 Fix.5.4 current-lineage numerical identity/SIMC-normalization audit.
```

```text
BLOCKED:
scientific interpretation of the current Fix.5 comparison pages until the
four numerical/provenance checks above close.
```

```text
NEXT:
audit the exact current-lineage producer -> sidecar -> payload -> renderer path
and implement only the numerical invariants or narrow source repairs warranted
by that audit.
```

Do not make visualization work or another farm run the sole ordinary NEXT yet.

Record the planned dependency order:

1. this memory checkpoint -> ChatGPT PASS -> user commit/push -> pushed-state review;
2. separate Fix.5.4 numerical audit/fix contract;
3. Codex numerical implementation;
4. ChatGPT actual-diff review;
5. user commit/push;
6. pushed-state review;
7. only then write/run the visualization-only contract from that reviewed pushed numerical source, after numerical closure;
8. Codex visualization implementation;
9. ChatGPT actual-diff review;
10. user commit/push;
11. pushed-state review;
12. one narrow Q4p4W2p74 / Left / lowe farm run;
13. fresh scientific and visual evidence review.

---

## 10. Required durable-memory changes

Update only warranted durable knowledge.

### CURRENT.md

Replace the stale pre-farm active state and NEXT as described above.

### MEMORY.md

Add durable rules:

1. an explanatory aggregate page and a downstream applied branch must not be treated as the same correction lineage without exact source/artifact/fingerprint identity;
2. E.8 current-branch yield-impact presentations require explicit histogram-scalar closure;
3. signed spectra require positive/negative/absolute support diagnostics when fractional signed-yield changes appear visually disproportionate;
4. SIMC/data amplitude interpretation requires proven same-object, same-normalization, same-unit provenance;
5. renderer success does not establish scientific linkage consistency.

Do not encode transient numerical plot values as universal constants.

### LEARNINGS.md

Add generalized lessons:

- sequential procedure pages can be individually correct yet scientifically misleading when they silently cross lineage boundaries;
- every displayed scalar yield attached to a histogram should be auditable against that histogram's producer-owned integral;
- signed-integral percentages can be misleading without positive/negative/absolute support;
- overlapping comparison curves must remain independently visible;
- visualization fixes must never alter physics or normalization.

### roadmap/STATUS.md

Add Fix.5.4 as `ACTIVE` and state that the fresh Fix.5 scientific interpretation is blocked on numerical/provenance closure.

Do not downgrade historical closures outside the exact affected interpretation.

### Existing Fix.5 phase record

Append a dated post-farm observation section. Preserve the original chronology.

Do not rewrite the historical local-development account as though the new findings were known earlier.

### New investigation record

Create:

```text
docs/memory/investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md
```

It must separate:

- direct observations;
- confirmed lineage finding;
- unvalidated hypotheses/checks;
- scientific boundaries;
- exact planned audit;
- visualization issue;
- farm/package caveat;
- exact next sequence.

### New phase record

Create:

```text
docs/memory/phases/e8-4-fix5-4-current-lineage-identity-audit.md
```

Status `ACTIVE`.

This record owns the planned identity audit but does not authorize scientific source edits by itself. A separate implementation contract will be written after this memory checkpoint is independently reviewed and pushed.

### USER.md

Do not change merely for this task-specific visualization issue.

### CURRENT_HANDOFF.md

Do not change; no exceptional handoff exists.

### manifest.json

Regenerate whenever versioned memory changes.

---

## 11. Allowed files

- `docs/memory/CURRENT.md`
- `docs/memory/MEMORY.md`
- `docs/memory/LEARNINGS.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/e8-4-fix5-shareable-method-a-impact-pages.md`
- `docs/memory/phases/e8-4-fix5-4-post-farm-identity-audit-memory-checkpoint-task-contract.md`
- `docs/memory/phases/e8-4-fix5-4-current-lineage-identity-audit.md`
- `docs/memory/investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md`
- `docs/memory/manifest.json`

No other tracked file may change.

---

## 12. Frozen files and ownership

All analysis/source/runtime/test files are frozen.

In particular do not modify:

- `src/`
- `testing/`
- `tools/`
- `run_Prod_Analysis.sh`
- background configuration
- Method-A/F.3/F.4/F.5/F.6 source
- E.8 plotting code
- SIMC/yield producers
- collectors/profiles/owners.

No scientific calculation is changed by this checkpoint.

---

## 13. Checks

Positive:

- CURRENT identifies `f9d707...` accurately.
- CURRENT has exactly one substantive NEXT: Fix.5.4 audit/fix.
- confirmed and unvalidated findings are clearly separated.
- old closed phases remain closed for their actual accepted scope.
- investigation and phase records exist.
- manifest includes every changed/new versioned-memory file.

Negative:

- no claim that the histogram/yield mismatch is already proven;
- no claim that signed cancellation is already the explanation;
- no claim that SIMC is incorrectly normalized;
- no claim that Method A is production-ready;
- no farm authorization;
- no source/test change.

Regression:

- Method B remains diagnostic-only.
- Method A remains detached.
- F.6.4 remains blocked.
- canonical-five E.8 remains unclosed.
- no accepted authority is replaced.

---

## 14. Local validation

Use the repository-selected Python interpreter.

Run the repository memory manifest writer/checker as established by current memory tooling, then:

```bash
<PYTHON> -B tools/check_memory_health.py --root .
```

Also run:

```bash
git diff --check
```

Report:

- `CURRENT.md` byte count;
- `MEMORY.md` byte count;
- `CURRENT_HANDOFF.md` byte count;
- memory-health hard failures;
- warnings classified blocking/nonblocking;
- manifest regeneration/check result.

Do not use `--fail-on-warning` unless required by an existing blocking warning.

---

## 15. Diff audit

Show:

```bash
git status --short --untracked-files=all
git diff --stat
git diff -- \
  docs/memory/CURRENT.md \
  docs/memory/MEMORY.md \
  docs/memory/LEARNINGS.md \
  docs/memory/roadmap/STATUS.md \
  docs/memory/phases/e8-4-fix5-shareable-method-a-impact-pages.md \
  docs/memory/manifest.json
```

For each new file include the complete addition with:

```bash
git diff --no-index -- /dev/null <new-file> || true
```

If the result is too large for terminal review, write one temporary `kaonlt_review.diff` in the repository root containing the complete tracked diff and complete additions for all new files.

Do not stage merely to make the diff reviewable.

---

## 16. Acceptance criteria

PASS only if:

- exact starting HEAD was respected;
- only allowlisted memory files changed;
- the fresh farm observations are recorded accurately;
- confirmed lineage mismatch is distinguished from unvalidated numerical hypotheses;
- CURRENT/NEXT sequence checkpoint review/push/synchronization -> numerical contract/implementation -> actual-diff review/push/pushed-state review -> visualization-only contract/implementation from that reviewed pushed numerical source -> actual-diff review/push/pushed-state review -> narrow farm run -> fresh evidence review;
- historical closed work is not unnecessarily reopened;
- manifest is regenerated;
- memory-health hard checks pass;
- `git diff --check` passes;
- no analysis source/test/runtime file changed.

---

## 17. Hard stop

After the memory-only checkpoint and its deterministic checks:

**STOP.**

Do not implement Fix.5.4 scientific/source changes.

Do not implement visualization changes.

Do not run the farm.

Do not commit.

Do not push.

Return the actual diff/review bundle to ChatGPT for independent review.
