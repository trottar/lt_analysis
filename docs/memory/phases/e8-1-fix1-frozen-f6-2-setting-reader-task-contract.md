# E.8.1.Fix.1 — align the ordinary E.8 frozen-F.6.2 reader with the accepted artifact

## Objective

Repair the ordinary E.8 procedure-PDF reader so it accepts the already accepted,
immutable F.6.2 aggregate artifact using that artifact's actual persisted setting
schema, while preserving every scientific, production, provenance, and
presentation boundary.

This is a narrow presentation-reader/test repair discovered by fresh Jefferson
Lab farm evidence. It does not reopen F.6.2 science and must not modify the
accepted F.6.2 JSON.

## Exact starting state

Start only from committed `test` HEAD:

`1c1a0be5da997043378fb49263b1de6ffc0d473c`

Before editing:

1. establish branch/HEAD/worktree;
2. read the normal five-file repository-memory startup sequence;
3. confirm the task-contract file itself is the only expected user-created
   uncommitted file;
4. stop if the starting HEAD differs or unrelated worktree changes exist.

The farm bundle that exposed this repair is:

`KaonLT_E8_1_Debug1_Fix1_Q4p4W2p74_Left_lowe_20260923-012755.zip`

Its relevant evidence is:

- bundle SHA-256:
  `490036e8d6c53d6e025da99a5c6c4ddcd529b6573007534efb7aaceac1acd002`
- bundle manifest: `complete=false`
- frozen F.6.2 JSON SHA-256 remained exactly:
  `5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1`
- Left/lowe procedure-PDF SHA-256:
  `d1bb7f945aee7b04b4c9baa0e25ba9b9e0d38956ad42086737b938a6b97ade01`
- Left/lowe page-manifest SHA-256:
  `6bebd1b8b0382b946c7e9e4ef5e0eca5613dc2871d5e31f5c9020ddc09574828`
- the procedure page manifest ended with:
  `full_background.e8.unavailable`
- the literal E.8 reason was:
  `frozen_f6_2_authority_rejected:parent_setting_identity_invalid`
- the bundle source-check failure was the stale
  `testing.test_full_background_subtraction_plots`
  `test_root_rendering_preserves_source_and_page_order` expectation.

The accepted F.6.2 JSON is immutable scientific evidence. Its parent setting
metadata uses, for example, Left/lowe:

```json
{
  "Q2": "4p4",
  "W": "2p74",
  "epsilon_filename_token": "lowe",
  "epsilon_setting": "low",
  "kinematic_token": "Q4p4W2p74",
  "particle_type": "kaon",
  "phi_setting": "Left"
}
```

High-epsilon parents analogously use
`epsilon_filename_token="highe"` and `epsilon_setting="high"`.
The accepted artifact does not persist an `ordinal` field.

## Scientific ownership

This repair is presentation-reader and regression-test work only.

F.6.2 remains `CLOSED / RUNTIME VALIDATED`. Its accepted JSON, artifact
fingerprint, validation fingerprint, numerical content, populations, geometry,
bootstrap results, Method-A quantities, and all scientific interpretation are
frozen.

Method A remains detached/non-production. Method B remains diagnostic only.
This task must not change production weights, corrected yields, cuts, binning,
normalization, templates, priors, subtraction formulas, SIMC, efficiencies,
cross sections, L/T separation, or uncertainties.

## Allowed implementation files

Only:

- `src/cuts/full_background_subtraction_plots.py`
- `testing/test_full_background_subtraction_plots.py`

## Allowed durable-memory files

Update only as warranted:

- `docs/memory/CURRENT.md`
- `docs/memory/phases/e8-1-fix1-frozen-f6-2-setting-reader.md`
- `docs/memory/manifest.json`

This task contract itself is also an intended tracked file:

- `docs/memory/phases/e8-1-fix1-frozen-f6-2-setting-reader-task-contract.md`

Regenerate `docs/memory/manifest.json` after all versioned memory changes.

Do not change `MEMORY.md`, `USER.md`, `CODEX.md`, `TOOLS.md`,
`CURRENT_HANDOFF.md`, decisions, roadmap, or evidence records unless an actual
blocker proves that their owned durable knowledge changed. Stop and report such
a blocker rather than expanding scope automatically.

## Frozen files and interfaces

Do not modify:

- the accepted F.6.2 JSON or any farm artifact;
- `run_Prod_Analysis.sh`;
- `src/main.py`;
- `src/cuts/diamond.py`;
- `src/cuts/rand_sub.py`;
- pion/proton subtraction, HGCer scientific calculation, SIMC, binning, yield,
  efficiency, cross-section, or uncertainty code;
- `testing/render_pion_hgcer_method_a_acceptance_refinement_figure_library.py`;
- E.8.1 collector/profile files;
- `testing/collect_pion_hgcer_validation_bundle.py`;
- `testing/package_pion_hgcer_validation_bundle.tcsh`;
- any F.6.2 analyzer/renderer/profile;
- the ordinary D.6-D.9 retained procedure-page semantics.

Do not reintroduce D.10, D.11, E.2-E.7.2 into the ordinary E.8 procedure PDF.
The current public renderer intentionally retains D.6-D.9 and appends E.8.

## Required source repair

In `src/cuts/full_background_subtraction_plots.py`, keep the existing E.8
reader fail-closed and keep all accepted SHA/fingerprint/schema checks.

Repair only the persisted parent-setting validation.

The current invalid assumption is equivalent to:

```python
"{}-{}".format(
    setting.get("phi_setting"),
    setting.get("epsilon_setting"),
) == setting_id
```

plus a required synthetic `ordinal`.

That is incompatible with the accepted F.6.2 artifact because
`epsilon_setting` is semantic (`low` / `high`) while the canonical setting ID
uses filename tokens (`lowe` / `highe`), and `ordinal` is not persisted.

Validate the accepted setting contract explicitly instead:

- `Q2 == "4p4"`
- `W == "2p74"`
- `kinematic_token == "Q4p4W2p74"`
- `particle_type == "kaon"`
- `phi_setting` equals the phi component of `setting_id`
- `epsilon_filename_token` equals the epsilon-token component of `setting_id`
- semantic epsilon is consistent:
  - `lowe -> low`
  - `highe -> high`

Do not infer or manufacture an `ordinal`.
Do not mutate the loaded artifact to add one.
Do not weaken the canonical 15-parent / 135-child inventory checks.
Do not replace exact accepted SHA/fingerprint gates with schema-only acceptance.

Prefer a small dedicated helper for this validation if it makes the contract
clearer, but do not refactor adjacent E.8 logic.

## Required regression coverage

Update `testing/test_full_background_subtraction_plots.py` so its E.8 authority
fixture represents the actual accepted persisted setting schema rather than the
synthetic `epsilon_setting="lowe"/"highe"` plus `ordinal` shape.

Positive coverage must prove that an otherwise valid frozen authority using the
actual persisted setting metadata is accepted and selects exactly the expected
three parents for the requested setting.

Negative coverage must fail closed for at least:

- mismatched `phi_setting`;
- mismatched `epsilon_filename_token`;
- inconsistent semantic `epsilon_setting`;
- wrong `Q2`;
- wrong `W`;
- wrong `kinematic_token`;
- wrong `particle_type`;
- missing/non-mapping parent `setting`.

Do not weaken unrelated malformed-payload tests.

## Required stale ROOT-test repair

The farm bundle also exposed one stale ROOT-aware assertion in
`test_root_rendering_preserves_source_and_page_order`.

Update that test to match the current public-renderer contract:

- retained D.6-D.9 pages remain in canonical per-t order;
- D.10/D.11 and E.2-E.7.2 are not ordinary pages;
- when no E.8 payload is supplied, the final page is the explicit E.8
  unavailable setting page;
- when valid E.8 is supplied in focused coverage, the final E.8 section follows
  the current E.8 page contract;
- source ROOT objects remain unmodified.

Do not change production rendering merely to make the stale test pass.

## Before/after behavior

Before:

`accepted F.6.2 JSON -> exact SHA/fingerprint checks -> incorrect parent setting
sub-dict check -> parent_setting_identity_invalid -> E.8 unavailable page`

After:

`accepted F.6.2 JSON -> exact SHA/fingerprint checks -> exact persisted setting
schema check -> full canonical inventory validation -> select current setting's
three parents -> ordinary E.8 final pages`

All non-E.8 science and D.6-D.9 page generation are unchanged.

## Deterministic local validation

Run, at minimum:

```text
python -m py_compile \
  src/cuts/full_background_subtraction_plots.py \
  testing/test_full_background_subtraction_plots.py
```

Run the focused E.8 tests from `testing.test_full_background_subtraction_plots`
that do not require unavailable local ROOT first.

Then run the complete module:

```text
python -m unittest testing.test_full_background_subtraction_plots -v
```

If local ROOT is unavailable, report the actual skip behavior; do not claim the
farm-only ROOT path passed locally.

Also run the normal repository memory manifest/health/bootstrap checks required
by current repository instructions and:

```text
git diff --check
```

Codex must report every command and exact result. Local checks are not farm
validation.

## Required diff audit

Before stopping:

1. show exact starting HEAD and current HEAD;
2. show `git status --short`;
3. show `git diff --stat`;
4. show `git diff --name-only`;
5. confirm every changed path is allowlisted;
6. inspect the actual diff for scientific/provenance boundary violations;
7. confirm frozen files are unchanged.

If the diff is too large for practical review, write one complete temporary
root-level review file such as `kaonlt_review.diff` containing the tracked diff
plus complete `git diff --no-index /dev/null ...` additions for every intended
new file. Do not stage files merely to create the review.

## Durable-memory update

`CURRENT.md` should record the fresh farm finding without overstating it:

- E.8.1 is `ACTIVE` because farm validation exposed a concrete E.8 presentation
  reader/test blocker.
- E.8.1.Fix.1 is `ACTIVE` until ChatGPT independently reviews the actual diff.
- F.6.2 and F.6.2.Fix.5 remain `CLOSED / RUNTIME VALIDATED`.
- E.8.1.Debug.1.Fix.1 remains `SOURCE REVIEWED`; the fresh Left/lowe output is
  runtime evidence that the targeted execution reached downstream procedure
  rendering, but the supplied bundle alone does not establish full launcher
  provenance or prove the intentional skip of full high epsilon.
- F.6.3 remains `BLOCKED` pending E.8.1.
- the next source action is this reader/test repair.
- after this repair is reviewed and user-pushed, a separate narrow E.8.1 bundle
  profile re-pin is required before the next farm gate so the collector
  provenance recognizes the repaired E.8 source commit.
- after that profile re-pin, rerun `Q4p4W2p74 / Left / lowe` first because it is
  the setting that exposed this blocker. Only after that narrow gate passes
  return to the planned detailed `Left / highe` E.8.1 review and then broaden.

Create
`docs/memory/phases/e8-1-fix1-frozen-f6-2-setting-reader.md`
as the permanent phase/fix chronology for this repair. Do not create an
accepted runtime evidence record yet.

## Forbidden shortcuts

Do not:

- edit or regenerate the accepted F.6.2 JSON;
- change its expected SHA-256 or fingerprints;
- accept arbitrary same-schema input;
- add an `ordinal` to frozen data;
- reinterpret `low` as `lowe` by mutating persisted content;
- suppress the unavailable/error path;
- delete the failing test;
- restore obsolete D.10/D.11/E.2-E.7.2 ordinary pages merely to satisfy the
  stale assertion;
- modify production physics;
- change the E.8.1 profile in this task;
- guess the future repaired source commit in the profile;
- commit, push, or run the farm.

## Acceptance criteria

Source-review candidate PASS requires all of the following:

- the E.8 reader validates the exact accepted F.6.2 parent-setting schema;
- exact accepted SHA/fingerprint/schema and canonical-inventory checks remain;
- actual persisted setting metadata is covered positively;
- malformed/mismatched setting metadata fails closed;
- stale public-renderer test is aligned to current D.6-D.9 + E.8 behavior;
- no frozen scientific or production source changed;
- deterministic local checks pass to local capability;
- durable memory accurately records the farm-discovered blocker and next
  profile-repin dependency;
- actual diff is narrow and allowlisted.

Farm/runtime acceptance is explicitly out of scope for this task.

## Hard stop

Do not commit, push, alter the E.8.1 profile, or run Jefferson Lab farm
validation.

Return to ChatGPT with:

- implementation summary;
- exact changed files;
- exact local test/check results;
- runtime/source trace;
- actual diff review artifact if needed;
- proposed user Git handoff commands, not executed.

ChatGPT must inspect the actual diff before the user commits or pushes.
