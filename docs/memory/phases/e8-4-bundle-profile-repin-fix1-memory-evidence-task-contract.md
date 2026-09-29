# KaonLT E.8.4 Bundle-Profile Re-pin Fix.1 — Preserve CURRENT Provenance and Record Required Memory Checks

## Purpose

Repair two narrow memory/evidence issues found by independent ChatGPT review of:

```text
kaonlt_review(20260929-010558).diff
```

The substantive bundle-profile re-pin itself is correct and must remain byte-identical:

- `testing/pion_hgcer_validation_bundle_profile_e8_1.json`
- `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`

This Fix.1 is **memory/evidence only**.

Independent review found:

1. `docs/memory/CURRENT.md` unnecessarily compressed away the existing exact
   E.8.3 / F.6.3 / E.8.4 source-review provenance that was present at pushed
   HEAD `1aa1fd...`. The re-pin status should be appended without deleting that
   established provenance.
2. The phase record reports the focused profile/test checks, but does not record
   the required memory-integrity checks or `git diff --check` result from the
   task contract. Those checks must be run now and their actual results recorded.
   Do not claim PASS for anything that did not pass.

No scientific, analysis, runtime, collector, wrapper, profile, or focused-test
change is authorized.

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

The worktree must contain the existing uncommitted cumulative bundle-profile
re-pin candidate corresponding to:

```text
kaonlt_review(20260929-010558).diff
```

plus this Fix.1 task contract after the user copies it into repository memory.

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git log -1 --oneline
git status --short
```

Hard stop if:

- committed HEAD differs from `1aa1fd...`;
- either substantive validation file differs from the reviewed candidate;
- collector/wrapper or any analysis/runtime file changed;
- unrelated work is present.

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
- `docs/memory/phases/e8-4-bundle-profile-repin.md`
- this Fix.1 contract
- the pushed `1aa1fd...` version of `docs/memory/CURRENT.md` for the exact
  pre-repin provenance text
- `docs/memory/roadmap/STATUS.md`

Do not reopen scientific architecture.

---

# 3. Freeze the accepted substantive re-pin

These files must remain byte-identical to the candidate reviewed in
`kaonlt_review(20260929-010558).diff`:

```text
testing/pion_hgcer_validation_bundle_profile_e8_1.json
testing/test_pion_hgcer_validation_bundle_profile_e8_1.py
```

In particular preserve:

```text
required_analysis_commit =
1aa1fd4184a6f8b20043e00ebb1ed3e9505a4935
```

and preserve the existing profile/test-only committed-range allowlist, exact
canonical-five setting inventory, and artifact inventory.

Also remain byte-unchanged:

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

---

# 4. CURRENT.md repair

In:

```text
docs/memory/CURRENT.md
```

restore the exact established E.8.3 / F.6.3 / E.8.4 source-review provenance
text that existed at pushed HEAD `1aa1fd...`.

Specifically, do not replace the detailed records naming:

```text
kaonlt_review(20260924-064150).diff
kaonlt_review(20260928-161516).diff
kaonlt_review(20260928-233040).diff
```

with a compressed summary.

Retain the new separate bundle-profile re-pin bullet:

```text
E.8.4 validation-bundle profile provenance re-pin — ACTIVE
required analysis source 1aa1fd...
pending independent ChatGPT review
no farm/runtime claim
```

Retain the current immediate NEXT:

```text
independent ChatGPT review of the refreshed complete bundle-profile re-pin diff
```

and retain the durable one-gate-at-a-time / no-premature-farm-command rule.

This is a restoration of existing provenance plus the new re-pin state, not a
redesign of CURRENT.

---

# 5. Run and record the required memory checks

Run exactly:

```bash
python -B tools/update_memory_manifest.py --root . --write
python -B tools/update_memory_manifest.py --root . --check
python -B tools/check_memory_health.py --root .
python -B tools/memory_bootstrap.py --root . --json
python -B -m unittest testing.test_memory_health -v
git -c core.safecrlf=false diff --check
```

If any check fails, stop and report the actual failure. Do not record PASS.

In:

```text
docs/memory/phases/e8-4-bundle-profile-repin.md
```

append a concise Fix.1 review-repair section that records the actual outcome of
each command above.

Also preserve the already-recorded focused checks:

```text
py_compile profile test
json.tool profile
5 focused profile tests
28 collector tests
```

Do not claim ChatGPT ran any of these local commands.

---

# 6. USER.md

The newly added durable preference is correct and must remain unchanged:

```text
Do not provide or authorize farm-run commands until the exact validation
bundle/profile required for that gate is ready, independently source-reviewed,
pushed, and pushed-state reviewed. Advance one workflow gate at a time; do not
provide commands for later gates before the current gate has passed.
```

No further USER.md edit is needed.

---

# 7. Allowed changes

Only:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/e8-4-bundle-profile-repin.md
docs/memory/phases/e8-4-bundle-profile-repin-fix1-memory-evidence-task-contract.md
```

`docs/memory/roadmap/STATUS.md`, `docs/memory/USER.md`, the profile, and its test
must remain byte-identical to the reviewed candidate unless a concrete blocker
requires stopping.

---

# 8. Status

Keep:

```text
E.8.4 — SOURCE REVIEWED
E.8.4 bundle-profile re-pin — ACTIVE
final E.8 — BLOCKED
F.6.4 — BLOCKED
```

Do not promote the bundle-profile re-pin to `SOURCE REVIEWED` yourself. That
status is reserved for the subsequent independent ChatGPT review.

No farm/runtime status is changed.

---

# 9. Diff audit

Run:

```bash
git status --short
git diff --stat
git diff -- docs/memory/CURRENT.md
git diff -- docs/memory/phases/e8-4-bundle-profile-repin.md
git diff -- testing/pion_hgcer_validation_bundle_profile_e8_1.json
git diff -- testing/test_pion_hgcer_validation_bundle_profile_e8_1.py
git -c core.safecrlf=false diff --check
```

Verify the two substantive validation files are unchanged from the reviewed
`20260929-010558` candidate.

Create a fresh complete cumulative review bundle containing all tracked changes
and no-index sections for all intended new/untracked memory files.

Do not stage merely to create the bundle.

---

# 10. Hard stop

Stop after:

- restoring the exact prior CURRENT source-review provenance;
- running and recording the required memory checks;
- regenerating/checking the manifest;
- auditing the cumulative diff;
- creating the refreshed review bundle.

Do not commit.
Do not push.
Do not run the farm.
Do not provide or prepare farm/package commands.

Exact NEXT:

```text
independent ChatGPT review of the refreshed complete bundle-profile re-pin diff
```
