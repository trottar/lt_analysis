# E.8.1.Fix.5 — actual-diff review repair

## Purpose

Repair two narrow issues found by independent ChatGPT inspection of the complete
E.8.1.Fix.5 `kaonlt_review.diff`.

The Fix.5 renderer implementation itself is not being reopened by this repair.
The reviewed renderer patch remains the intended candidate:

- overlay canvas `1800 x 3600` -> `3600 x 3600`;
- draw/populate the canonical `3 x 9` grid first;
- draw the reserved header pad as the final top-level sibling;
- keep the L/B/A legend in that header;
- call `Modified()` / `Update()` before PDF `Print()`.

This repair only corrects a non-representative real-ROOT regression/evidence
label and restores two CURRENT continuity facts that were accidentally removed.

## Exact committed base

Committed `test` HEAD must remain:

`fdd368f2084c9b9508e7ce679b8f7391f5b556f1`

The existing uncommitted E.8.1.Fix.5 change set is intentional.

Before editing:

1. establish branch/HEAD/worktree;
2. read the normal five-file repository-memory startup sequence;
3. read the Fix.5 task contract and phase record;
4. confirm no unrelated changes exist.

## Independent review result

ChatGPT inspected the complete actual Fix.5 diff.

The presentation implementation is narrowly scoped and otherwise matches the
Fix.5 contract. However, the new conditional real-ROOT PDF regression changes
the frozen fixture's `phi_index` values from the canonical persisted `0..8` to
synthetic `1..9` values:

```python
for phi_index, child in enumerate(parent["children"], 1):
    child["phi_index"] = phi_index
```

and then hard-codes:

```text
Left-lowe phi1 [-180, -140)
```

That is not the actual persisted E.8 child-label contract.

The authoritative reader requires each child `phi_index` to equal the
zero-based `child_position` from `enumerate(children)`. `_e8_child_label()`
prints that stored `phi_index` directly. Therefore the actual first canonical
Left/lowe child is:

```text
Left-lowe phi0 [-180, -140)
```

The supplied fresh farm PDF confirms that `phi0 [-180, -140)` is absent from
page 38. `phi1 [-140, -100)` is also absent, while later child titles such as
`phi2 [-100, -60)` are extractable.

The regression must test the real persisted contract; it must not mutate frozen
child identity merely to satisfy a mistaken example string in the original
Fix.5 contract.

ChatGPT also found that the Fix.5 CURRENT update removed two still-valid
continuity facts unrelated to this repair:

- lifecycle-hook dispatch remains non-required and `BLOCKED / DEFERRED`;
- E.8.1.Debug.1.Fix.1 remains `SOURCE REVIEWED`, and the earlier incomplete
  output still does not prove launcher provenance or the intentional
  full-high-epsilon skip.

The relevant wrapper/debug links were also removed from CURRENT even though the
next farm gate still uses both the debug launcher and bundle wrapper. Restore
those links rather than performing unrelated CURRENT cleanup.

## Status

Preserve:

- E.8.1 — `ACTIVE`
- E.8.1.Fix.1 — `CLOSED / RUNTIME VALIDATED`
- E.8.1.Fix.2 — `CLOSED / RUNTIME VALIDATED`
- E.8.1.Fix.3 — `SOURCE REVIEWED`
- E.8.1.Fix.4 — `CLOSED / RUNTIME VALIDATED` only for its narrow provenance
  re-pin
- E.8.1.Fix.5 — `ACTIVE` pending a new independent actual-diff review
- E.8.1.Debug.1.Fix.1 — `SOURCE REVIEWED`
- F.6.3 — `BLOCKED`
- F.6.4 — `BLOCKED`
- lifecycle-hook dispatch — `BLOCKED / DEFERRED`

Do not mark Fix.5 `SOURCE REVIEWED` in this repair pass.

## Allowed modifications in this repair

Only:

- `testing/test_full_background_subtraction_plots.py`
- `docs/memory/CURRENT.md`
- `docs/memory/evidence/e8-1-fix4-left-lowe-overlay-blocker.md`
- `docs/memory/phases/e8-1-fix5-overlay-pdf-geometry-task-contract.md`
- `docs/memory/manifest.json`

This repair contract is also an intended tracked file:

- `docs/memory/phases/e8-1-fix5-actual-diff-review-repair-task-contract.md`

## Frozen during this repair

Do not make any further change to:

- `src/cuts/full_background_subtraction_plots.py`
- `docs/memory/phases/e8-1-fix3-procedure-pdf-layout.md`
- `docs/memory/phases/e8-1-fix4-bundle-profile-repin.md`
- `docs/memory/phases/e8-1-fix5-overlay-pdf-geometry.md`
- any profile;
- any collector or wrapper;
- `run_Prod_Analysis.sh`;
- `src/main.py`;
- `src/cuts/rand_sub.py`;
- accepted F.6.2 artifacts;
- production/scientific source.

The current uncommitted Fix.5 renderer patch must remain byte-for-byte
unchanged during this repair.

## Required real-ROOT regression correction

In:

`testing/test_full_background_subtraction_plots.py`

repair only the new
`test_e8_real_root_overlay_pdf_emits_header_legend_and_first_child_title`.

### Do not mutate frozen child identity

Delete the synthetic one-based rewrite:

```python
for phi_index, child in enumerate(parent["children"], 1):
    child["phi_index"] = phi_index
```

Do not replace it with any other mutation of `phi_index`, `phi_low`,
`phi_high`, `setting_id`, or canonical t identity.

### Use a validated canonical parent

Prefer to derive the test parent through the actual E.8 validation contract:

```python
artifact = deepcopy(_e8_frozen_artifact())
validation = plots._e8_validate_artifact(artifact)
parent = deepcopy(next(
    parent
    for parent in validation["parents"]
    if parent["setting_id"] == "Left-lowe"
    and parent["canonical_t_index"] == 0
))
```

An equivalent approach is acceptable only if it proves the parent satisfies the
same zero-based persisted child contract before rendering.

Add an explicit contract assertion:

```python
self.assertEqual(
    [child["phi_index"] for child in parent["children"]],
    list(range(9)),
)
```

Derive the first child display label from production code:

```python
first_child_title = plots._e8_child_label(parent["children"][0])
self.assertEqual(
    first_child_title,
    "Left-lowe phi0 [-180, -140)",
)
```

Then require `first_child_title` in the extracted overlay page instead of the
incorrect hard-coded `phi1 [-180, -140)` string.

Retain the header and all three L/B/A legend assertions.

Retain the multipage-PDF context.

Do not weaken the conditional PyROOT/`pdftotext` behavior.

## Required evidence correction

In:

`docs/memory/evidence/e8-1-fix4-left-lowe-overlay-blocker.md`

replace the incorrect first-child statement with the actual persisted labels.

Record:

- actual first canonical title
  `Left-lowe phi0 [-180, -140)` is absent on page 38;
- `Left-lowe phi1 [-140, -100)` is also absent;
- later titles, including `Left-lowe phi2 [-100, -60)`, are extractable.

Do not change the bundle/provenance/hash/test facts.

## Required original-contract correction

In:

`docs/memory/phases/e8-1-fix5-overlay-pdf-geometry-task-contract.md`

correct the mistaken real-ROOT regression example from:

```text
Left-lowe phi1 [-180, -140)
```

to:

```text
Left-lowe phi0 [-180, -140)
```

Add a short note that the canonical persisted `phi_index` contract is zero
based (`0..8`) and regression fixtures must not renumber it.

Do not otherwise rewrite the original contract.

## Required CURRENT restoration

In `docs/memory/CURRENT.md` preserve the new Fix.4/Fix.5 state, but restore:

```text
E.8.1.Debug.1.Fix.1 remains SOURCE REVIEWED. Its earlier incomplete output does
not prove launcher provenance or the intentional full-high-epsilon skip.
```

and in Blockers restore the still-current statement:

```text
Lifecycle-hook dispatch is non-required and BLOCKED / DEFERRED.
```

Under Relevant References restore:

- `[Wrapper reconciliation](phases/farm-validation-bundle-wrapper-fix1-state-reconciliation.md)`
- `[E.8.1.Debug.1 local debug record](phases/e8-1-debug-1-left-low-debug-mode.md)`
- `[E.8.1.Debug.1.Fix.1 local repair record](phases/e8-1-debug-1-fix1-preserve-diamond-cut.md)`

Do not restore the superseded Fix.2 blocker as the primary current blocker; the
fresh Fix.4 overlay evidence remains the active E.8.1 blocker reference.

## Manifest and validation

Regenerate:

`docs/memory/manifest.json`

Run at minimum:

```text
python -B -m py_compile \
  src/cuts/full_background_subtraction_plots.py \
  testing/test_full_background_subtraction_plots.py

python -B -m unittest \
  testing.test_full_background_subtraction_plots.FullBackgroundSubtractionE8Tests \
  -v

python -B -m unittest testing.test_full_background_subtraction_plots -v
```

The conditional real-PDF test may still skip locally when PyROOT is
unavailable. Report exact tests and skips.

Run the normal memory manifest/integrity/health/bootstrap checks and:

```text
git diff --check
```

## Required diff audit

Before stopping:

1. confirm committed HEAD remains
   `fdd368f2084c9b9508e7ce679b8f7391f5b556f1`;
2. confirm `src/cuts/full_background_subtraction_plots.py` has no additional
   change relative to the already reviewed Fix.5 patch;
3. confirm no profile/collector/wrapper/production source changed;
4. confirm the real-ROOT regression now consumes zero-based canonical child
   identity without mutation;
5. confirm evidence and original contract use `phi0 [-180, -140)` for the first
   canonical child;
6. confirm CURRENT again records lifecycle `BLOCKED / DEFERRED` and the debug
   launcher provenance caveat;
7. regenerate one complete root-level `kaonlt_review.diff`, including this new
   untracked contract with `git diff --no-index /dev/null ...`.

Do not stage merely to create the review file.

## Hard stop

Do not commit, push, alter the validation profile, or run the Jefferson Lab
farm.

Return the refreshed root-level `kaonlt_review.diff` for independent ChatGPT
review.
