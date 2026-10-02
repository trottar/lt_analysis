# KaonLT — E.8.4 Fix.5.4 Fix.1 SIMC-availability granularity and manifest-scope repair

## Status

`ACTIVE`

This is a narrow repair to the local E.8.4 Fix.5.4 implementation candidate
built on committed `test` HEAD:

```text
71a912bf14c59eee6d52956462aa0aa39c7f7c84
```

The current working tree already contains the reviewed-but-not-accepted
Fix.5.4 implementation. Preserve that candidate and repair only the two
defects identified by independent ChatGPT actual-diff review.

Do not commit, push, or run the farm.

---

## 1. Starting state and preservation

Required committed branch/HEAD:

```text
branch: test
HEAD:   71a912bf14c59eee6d52956462aa0aa39c7f7c84
```

Before editing, report:

```bash
git branch --show-current
git rev-parse HEAD
git status --short --untracked-files=all
```

The working tree is expected to contain the current Fix.5.4 candidate and may
also contain unrelated pre-existing user-owned untracked files.

Do not:

- reset;
- clean;
- stash;
- discard;
- overwrite;
- commit;
- push;
- run the farm.

Preserve unrelated user-owned work exactly.

Read the normal startup core and the existing implementation contract:

```text
AGENTS.md
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/handoffs/CURRENT_HANDOFF.md
docs/memory/USER.md
docs/memory/phases/e8-4-fix5-4-current-lineage-identity-audit-and-fix-task-contract.md
```

This Fix.1 contract narrows only the two review findings below. All scientific
ownership, frozen interfaces, tests, and hard-stop requirements from the parent
Fix.5.4 contract remain in force.

---

## 2. Review finding A — SIMC unavailability is currently too coarse

### Defect

The local candidate currently makes the whole E.8.4 payload availability equal
to:

```python
simc_audit["absolute_comparison_available"]
```

The source audit correctly found that the current repository does **not** prove
an absolute same-unit luminosity/charge contract for the SIMC MM support, so
that flag is currently false.

The existing renderer treats `e8_4_payload["available"] == False` as meaning
that the entire E.8.4 audit is unavailable. Therefore the candidate would
suppress not only the absolute SIMC overlays but also already-valid current
F.6.3 data-only evidence such as:

- authority/lineage information;
- pion consequence;
- final `MM_0` versus `MM_A`;
- signed differences;
- `Y0/YA` impact;
- current-lineage identity/aggregate audit;
- parent closure.

That is broader than the parent contract.

### Required behavior

Separate these two concepts:

```text
E.8.4 current-F.6.3 audit availability
```

and:

```text
absolute SIMC amplitude-comparison availability
```

If the current F.6.3 identity audit and all non-SIMC payload validation pass:

```python
payload["available"] is True
```

even when:

```python
payload["simc_absolute_comparison_available"] is False
```

The SIMC audit must remain attached with its exact unavailable reason and
source provenance.

Only presentation whose scientific meaning requires an **absolute**
data-versus-SIMC amplitude comparison may be unavailable/provenance-blocked.

At minimum this applies to:

```text
full_background.e8_4.method_a_vs_simc.t1/t2/t3
full_background.e8_4.baseline_method_a_simc.t1/t2/t3
```

Do not render the black SIMC curve as an authoritative absolute-amplitude
comparison when the audit says that comparability is unproven.

The rest of E.8.4 must remain renderable and auditable.

Implement a narrow explicit SIMC-unavailable page/status for those affected
SIMC comparison pages, or another equally explicit page-level mechanism that:

- preserves the expected page-level provenance;
- states the literal normalization/provenance blocker;
- does not claim absolute agreement/disagreement;
- does not rescale SIMC;
- does not suppress unrelated E.8.4 pages.

Do not perform visualization-style changes in this repair.

### Required tests

Update/add deterministic tests proving:

1. valid current-F.6.3 identity + unresolved absolute SIMC units:
   - overall E.8.4 payload remains available;
   - `simc_absolute_comparison_available` is false;
   - literal reason is preserved;
   - data-only/current-lineage E.8.4 content remains available;
   - absolute-SIMC pages are explicitly unavailable/provenance-blocked.

2. invalid current-F.6.3 identity:
   - overall E.8.4 remains unavailable as before.

3. source-proven comparable SIMC synthetic case:
   - normal SIMC comparison path remains available;
   - no display normalization is added.

4. no production, Method-A, yield, SIMC normalization, or child contents are
   mutated.

Do not weaken any existing current-lineage identity checks.

---

## 3. Review finding B — unrelated untracked contract contaminated manifest/review scope

### Defect

The current review bundle includes the unrelated pre-existing user-owned file:

```text
docs/memory/phases/workflow-continuity-hardening-task-contract.md
```

and `docs/memory/manifest.json` now indexes it.

That file is outside Fix.5.4 scope. It must neither be silently committed as
part of Fix.5.4 nor be deleted/discarded.

A candidate manifest that references a file which will remain untracked after
the Fix.5.4 commit is also invalid.

### Required behavior

Preserve this unrelated file exactly as user-owned untracked work:

```text
docs/memory/phases/workflow-continuity-hardening-task-contract.md
```

but exclude it from:

- the Fix.5.4 candidate diff;
- the Fix.5.4 review bundle;
- the versioned `docs/memory/manifest.json` produced for this candidate.

Do not delete or rewrite it.

The intended Fix.5.4 implementation contract:

```text
docs/memory/phases/e8-4-fix5-4-current-lineage-identity-audit-and-fix-task-contract.md
```

is part of this Fix.5.4 work and may remain in the candidate/manifest.

This Fix.1 repair contract itself is also intended to become a versioned phase
contract and therefore may be included in the candidate/manifest.

Regenerate and check `docs/memory/manifest.json` against the **intended
versioned candidate state**, not against unrelated user-owned untracked
memory files.

Use a safe temporary/detached candidate view or another non-destructive method
if needed. Do not use `git clean`, `git reset`, or `git stash`.

Update any task memory wording that currently says the manifest intentionally
indexes both pre-existing contracts. That statement must be removed/corrected.

---

## 4. Scientific boundaries

Do not change:

- Method-A factors or correction mathematics;
- `B_pi_0`, `B_pi_A`, `MM_0`, `MM_A`, `Y0`, or `YA` arithmetic;
- the new binwise/yield/fingerprint identity audit;
- current-lineage aggregate mathematics;
- the scientific conclusion that absolute SIMC amplitude comparability is not
  source-proven unless a concrete source contradiction is discovered;
- SIMC `normfac`, event weighting, model weighting, or histogram contents;
- random/dummy/proton/pion subtraction;
- cuts, windows, binning, fits, priors, templates, efficiencies, acceptance,
  yields, cross sections, or L/T separation;
- Method-B diagnostic-only status;
- production objects.

Do not add any SIMC scale or child renormalization.

Do not do the later blue/magenta visualization work here.

The owner/collector/profile packaging failure remains a separate issue.

---

## 5. Allowed substantive paths

Repair only as required within the current Fix.5.4 candidate:

```text
src/cuts/full_background_subtraction_plots.py
testing/test_e8_4_production_impact_audit.py
testing/test_e8_4_fix5_4_identity_audit.py
```

`src/binning/calculate_yield.py` should remain unchanged unless a deterministic
test exposes a direct need to propagate an already-existing field; no numerical
producer behavior may change.

Warranted memory paths:

```text
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md
docs/memory/phases/e8-4-fix5-4-current-lineage-identity-audit.md
docs/memory/roadmap/STATUS.md
docs/memory/manifest.json
docs/memory/phases/e8-4-fix5-4-current-lineage-identity-audit-and-fix-task-contract.md
docs/memory/phases/e8-4-fix5-4-fix1-simc-availability-and-manifest-scope-task-contract.md
```

The parent implementation contract should not be rewritten unless necessary to
record a factual Fix.1 cross-reference.

Everything else is frozen.

In particular, the unrelated:

```text
docs/memory/phases/workflow-continuity-hardening-task-contract.md
```

must remain untouched and untracked by this candidate.

---

## 6. Memory/status requirements

After successful local repair, Fix.5.4 may remain:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

only if all deterministic checks still pass.

Memory must say precisely:

- current F.6.3 numerical identity audits are implemented locally;
- absolute SIMC amplitude comparability remains source-unproven;
- this blocks only authoritative absolute SIMC comparison/interpretation, not
  the entire current-F.6.3 E.8.4 audit;
- visualization work remains the next source-changing stage only after
  ChatGPT actual-diff review, user commit/push, and pushed-state review;
- no farm run is authorized yet.

Do not upgrade to runtime validation.

---

## 7. Local validation

Rerun the focused Fix.5.4 tests and relevant regressions, including at minimum:

```text
testing/test_e8_4_fix5_4_identity_audit.py
testing/test_e8_4_production_impact_audit.py
testing/test_full_background_subtraction_plots.py
```

plus the same F.6.3/yield/SIMC regressions used by the parent task.

Run syntax checks for all changed/new Python files.

Run:

```bash
git diff --check
```

Regenerate/check the memory manifest using only the intended candidate state,
then run:

```bash
<PYTHON> -B tools/check_memory_health.py --root .
```

Report all tests, skips, warnings, and hard failures accurately.

Local checks do not establish ROOT/PyROOT, full `main.py`, procedure-PDF, or
farm behavior.

---

## 8. Diff/review requirements

The returned `kaonlt_review.diff` must contain the complete intended Fix.5.4
candidate relative to committed HEAD
`71a912bf14c59eee6d52956462aa0aa39c7f7c84`.

It must include all intended new files in full.

It must **not** include:

```text
docs/memory/phases/workflow-continuity-hardening-task-contract.md
```

because that file is unrelated preserved user work.

Before stopping, report:

```bash
git status --short --untracked-files=all
git diff --stat
git diff --check
```

The status may still show the unrelated contract as an untracked user file.
That is acceptable and must be reported, not removed.

---

## 9. Acceptance criteria

PASS only if:

1. the overall E.8.4 current-F.6.3 payload remains available when its identity
   audit passes, independent of unresolved absolute SIMC units;
2. absolute SIMC comparison availability is represented separately;
3. SIMC comparison pages fail closed/are provenance-blocked without suppressing
   data-only E.8.4 pages;
4. no SIMC rescaling or physics change is introduced;
5. all Fix.5.4 numerical identity work remains intact;
6. the unrelated workflow-continuity contract remains untouched/untracked;
7. the versioned manifest does not index that unrelated file;
8. the review bundle does not contain that unrelated file;
9. deterministic tests and memory-health checks pass;
10. no visualization-style changes, commit, push, or farm run occur.

---

## 10. Hard stop

After the narrow repair, deterministic checks, memory reconciliation, manifest
regeneration, and complete review bundle preparation:

**STOP.**

Do not implement the later visualization stage.

Do not repair the farm owner/packaging path.

Do not commit.

Do not push.

Do not run the farm.

Return the complete repaired `kaonlt_review.diff` to ChatGPT.
