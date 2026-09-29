# KaonLT E.8.4 Bundle-Profile Re-pin — Final Pre-Push Source-Review Reconciliation

## Purpose

Record the independent ChatGPT review result for the complete bundle-profile
re-pin candidate after Fix.1.

Reviewed artifact:

```text
kaonlt_review(20260929-100154).diff
```

Independent review result:

```text
PASS
```

The substantive bundle-profile re-pin is accepted for source/provenance scope:

- `testing/pion_hgcer_validation_bundle_profile_e8_1.json` requires pushed
  E.8.4 analysis source `1aa1fd4184a6f8b20043e00ebb1ed3e9505a4935`;
- `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py` has the matching
  exact source expectation;
- the profile/test-only committed-range allowlist is unchanged;
- canonical-five setting and artifact inventories are unchanged;
- collector, tcsh wrapper, analysis/runtime source, launcher, Method A/B
  ownership, and production physics are unchanged;
- `CURRENT.md` preserves the existing detailed E.8.3/F.6.3/E.8.4 source-review
  provenance;
- `docs/memory/USER.md` contains the durable one-gate-at-a-time rule;
- Codex-reported focused and memory checks are recorded in the phase record;
- `check_memory_health.py` exited 0 with a soft-size warning for CURRENT
  (`8349 > 8192` bytes), not a hard failure;
- no ROOT/PyROOT, procedure-PDF, Jefferson Lab farm, or runtime validation is
  claimed.

This task is **memory/status reconciliation only**.

---

# 1. Exact starting state

Required branch:

```text
test
```

Required committed HEAD:

```text
1aa1fd4184a6f8b20043e00ebb1ed3e9505a4935
Add E8.4 Method-A production-impact audit
```

The worktree must contain the exact cumulative candidate reviewed in:

```text
kaonlt_review(20260929-100154).diff
```

plus this reconciliation contract after the user copies it into
`docs/memory/phases/`.

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git log -1 --oneline
git status --short
```

Hard stop if any substantive validation file, USER record, collector/wrapper,
analysis/runtime source, or launcher differs from the reviewed candidate.

Do not reset, stash, clean, commit, push, or run the farm.

---

# 2. Required startup reading

Read in order:

1. local root `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read only:

- `docs/memory/phases/e8-4-bundle-profile-repin-task-contract.md`
- `docs/memory/phases/e8-4-bundle-profile-repin-fix1-memory-evidence-task-contract.md`
- `docs/memory/phases/e8-4-bundle-profile-repin.md`
- `docs/memory/roadmap/STATUS.md`
- this reconciliation contract

---

# 3. Frozen files

Remain byte-identical to `kaonlt_review(20260929-100154).diff`:

```text
testing/pion_hgcer_validation_bundle_profile_e8_1.json
testing/test_pion_hgcer_validation_bundle_profile_e8_1.py
docs/memory/USER.md
```

Remain byte-identical to committed HEAD:

```text
testing/collect_pion_hgcer_validation_bundle.py
testing/package_pion_hgcer_validation_bundle.tcsh
src/cuts/full_background_subtraction_plots.py
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
src/binning/calculate_yield.py
src/cuts/rand_sub.py
src/cuts/pion_component_subtraction.py
src/main.py
run_Prod_Analysis.sh
```

Do not change scientific or runtime behavior.

---

# 4. Required status reconciliation

Update only warranted active-state/history records so that the bundle-profile
re-pin is now:

```text
SOURCE REVIEWED
```

Specifically:

## `docs/memory/CURRENT.md`

Preserve all existing detailed E.8.3/F.6.3/E.8.4 source-review provenance.

Change only the bundle-profile re-pin state from `ACTIVE` to
`SOURCE REVIEWED`, recording:

- independent ChatGPT review of
  `kaonlt_review(20260929-100154).diff`;
- review result `PASS`;
- source/provenance scope only;
- no ROOT/PyROOT, procedure-PDF, farm, or runtime claim;
- Codex checks were not run by ChatGPT;
- the memory-health soft-size warning is non-fatal and remains explicitly
  distinct from a runtime/source failure.

The exact current NEXT must be only:

```text
user-controlled commit/push of the reviewed bundle-profile re-pin
```

Do not provide later-gate farm or packaging commands in CURRENT.

## `docs/memory/phases/e8-4-bundle-profile-repin.md`

Append a chronology-preserving final independent-review section recording:

```text
SOURCE REVIEWED
artifact: kaonlt_review(20260929-100154).diff
result: PASS
```

Record that:

- substantive profile/test content passed;
- Fix.1 restored provenance and recorded required checks;
- the `check_memory_health.py` result was exit 0 with the soft CURRENT-size
  warning, not a hard error;
- no farm/runtime claim exists;
- no production promotion exists.

Do not erase the prior `ACTIVE` chronology or Fix.1 history.

## `docs/memory/roadmap/STATUS.md`

Change only the bundle-profile re-pin status from `ACTIVE` to
`SOURCE REVIEWED`.

Keep:

```text
E.8.4 — SOURCE REVIEWED
final E.8 — BLOCKED
F.6.4 — BLOCKED
```

---

# 5. Allowed changes

Only:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/e8-4-bundle-profile-repin.md
docs/memory/phases/e8-4-bundle-profile-repin-final-source-review-reconciliation-task-contract.md
docs/memory/roadmap/STATUS.md
```

No other file may change.

---

# 6. Memory integrity

Run:

```bash
python -B tools/update_memory_manifest.py --root . --write
python -B tools/update_memory_manifest.py --root . --check
python -B tools/check_memory_health.py --root .
python -B tools/memory_bootstrap.py --root . --json
python -B -m unittest testing.test_memory_health -v
git -c core.safecrlf=false diff --check
```

Record actual outcomes. A repeat of the known CURRENT soft-size warning is
acceptable only if the command exits 0 and no hard health error appears.

Do not claim ChatGPT ran these commands.

---

# 7. Diff audit

Run:

```bash
git status --short
git diff --stat
git diff -- docs/memory/CURRENT.md
git diff -- docs/memory/phases/e8-4-bundle-profile-repin.md
git diff -- docs/memory/roadmap/STATUS.md
git diff -- testing/pion_hgcer_validation_bundle_profile_e8_1.json
git diff -- testing/test_pion_hgcer_validation_bundle_profile_e8_1.py
git -c core.safecrlf=false diff --check
```

Confirm the substantive profile/test diff is identical to
`kaonlt_review(20260929-100154).diff`.

Create a fresh complete cumulative review bundle including no-index sections
for intended new/untracked memory files.

---

# 8. Hard stop

Stop after memory/status reconciliation, manifest regeneration/checks, diff
audit, and review-bundle creation.

Do not commit.
Do not push.
Do not run the farm.
Do not provide or prepare farm/package commands.

Exact NEXT:

```text
independent ChatGPT review of the final reconciliation diff
```
