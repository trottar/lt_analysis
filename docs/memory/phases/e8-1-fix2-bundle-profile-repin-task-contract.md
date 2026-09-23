# E.8.1.Fix.2 — re-pin the E.8.1 validation bundle profile to the repaired reader source

## Objective

Re-pin the existing E.8.1 canonical-five-setting validation bundle profile to the
now-pushed, independently source-reviewed E.8.1.Fix.1 reader/test commit.

This is a narrow validation/provenance packaging change only. It must not alter
analysis physics, the ordinary procedure-PDF renderer, the generic collector, the
bundle wrapper, artifact declarations, canonical settings, or accepted F.6.2
evidence.

## Exact starting state

Start only from committed `test` HEAD:

`8985d9a212799c021c4ad1a759a689ea3826e0ea`

Commit message:

`Fix E8 frozen F6.2 setting reader`

Its parent is:

`1c1a0be5da997043378fb49263b1de6ffc0d473c`

Before editing:

1. establish branch/HEAD/worktree;
2. read the normal five-file repository-memory startup sequence;
3. confirm this task-contract file is the only expected user-created uncommitted
   file;
4. stop if the starting HEAD differs or unrelated worktree changes exist.

The pushed E.8.1.Fix.1 source commit is now the analysis/procedure source that
the E.8.1 collector profile must require:

`8985d9a212799c021c4ad1a759a689ea3826e0ea`

The current E.8.1 profile still pins the older reviewed procedure source:

`0bf1281ebf44868b68ff808090436653eea6ed60`

That old pin must not be broadened by allowlisting later analysis source. It
must be replaced by the exact repaired source commit.

## Scientific and architectural ownership

E.8.1 remains presentation-only and `ACTIVE` pending fresh farm validation.

F.6.2 and F.6.2.Fix.5 remain `CLOSED / RUNTIME VALIDATED`.

E.8.1.Fix.1 remains `SOURCE REVIEWED`; this task does not upgrade it to farm
validation.

Method A remains detached/non-production. Method B remains diagnostic only.

This task changes validation provenance only. It must not change:

- accepted F.6.2 JSON or fingerprints;
- cuts, binning, normalizations, templates, priors, fits, subtraction formulas,
  production weights, SIMC, yields, efficiencies, cross sections, L/T
  separation, or uncertainties;
- ordinary D.6-D.9 + E.8 procedure-PDF semantics;
- canonical-five-setting membership.

## Allowed implementation files

Only:

- `testing/pion_hgcer_validation_bundle_profile_e8_1.json`
- `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`

## Allowed durable-memory files

Update only as warranted:

- `docs/memory/CURRENT.md`
- `docs/memory/phases/e8-1-fix1-frozen-f6-2-setting-reader.md`
- `docs/memory/phases/e8-1-fix2-bundle-profile-repin.md`
- `docs/memory/manifest.json`

This task contract is also an intended tracked file:

- `docs/memory/phases/e8-1-fix2-bundle-profile-repin-task-contract.md`

Regenerate `docs/memory/manifest.json` after all versioned memory changes.

Do not change `MEMORY.md`, `USER.md`, `CODEX.md`, `TOOLS.md`,
`CURRENT_HANDOFF.md`, roadmap, decisions, or evidence records unless an actual
blocker proves that their owned durable knowledge changed. Stop and report such
a blocker rather than expanding scope automatically.

## Frozen files and interfaces

Do not modify:

- anything under `src/`;
- `run_Prod_Analysis.sh`;
- `testing/collect_pion_hgcer_validation_bundle.py`;
- `testing/package_pion_hgcer_validation_bundle.tcsh`;
- any other validation profile;
- any F.6.2 analyzer/renderer/profile;
- any farm artifact;
- accepted F.6.2 evidence;
- E.8.1 artifact declarations;
- E.8.1 canonical setting declarations.

The generic collector and wrapper already implement the required fail-closed
provenance semantics. Do not edit them.

## Required profile change

In:

`testing/pion_hgcer_validation_bundle_profile_e8_1.json`

change only:

```json
"required_analysis_commit": "0bf1281ebf44868b68ff808090436653eea6ed60"
```

to:

```json
"required_analysis_commit": "8985d9a212799c021c4ad1a759a689ea3826e0ea"
```

Preserve exactly:

- `schema_version`:
  `pion_hgcer_validation_bundle_profile/v4`
- `validation_profile`:
  `phase_e8_1_full_background_procedure_pdf_farm_review/v1`
- `collection_mode`:
  `generic_artifacts`
- the canonical five settings, in their existing order:
  - Left / lowe
  - Left / highe
  - Center / lowe
  - Center / highe
  - Right / highe
- the one required global frozen F.6.2 JSON artifact;
- the two required per-setting artifacts:
  - full-background-subtraction PDF
  - full-background-subtraction page-manifest JSON
- `allowed_committed_files`, exactly:
  - `testing/pion_hgcer_validation_bundle_profile_e8_1.json`
  - `testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`
- `allowed_non_analysis_path_prefixes`, exactly:
  - `docs/memory/`

Do not bump the profile schema or `validation_profile` identifier. The evidence
package contract is unchanged; only its exact required analysis-source identity
is advancing to the repaired source.

Do not add `src/cuts/full_background_subtraction_plots.py` or any other analysis
file to `allowed_committed_files`. The correct provenance model is a new exact
required source commit, not a broader exception list.

## Required test change

In:

`testing/test_pion_hgcer_validation_bundle_profile_e8_1.py`

update the exact `REVIEWED_SOURCE` constant from the old E.8.1 procedure source
to:

`8985d9a212799c021c4ad1a759a689ea3826e0ea`

Preserve the current source-identity tests and their fail-closed behavior.

The tests must continue to prove that:

- the profile declares exactly the canonical five settings;
- artifact declarations are unchanged;
- the profile's `required_analysis_commit` equals `REVIEWED_SOURCE`;
- only the profile/test paths are allowed committed files after the required
  source;
- `docs/memory/` remains the only allowed non-analysis prefix;
- a complete synthetic package can be collected;
- missing/invalid required artifacts remain incomplete;
- missing required source ancestry remains incomplete;
- analysis-source changes after the required source, including
  `src/cuts/full_background_subtraction_plots.py`, remain rejected;
- output collisions still fail without overwrite.

Do not weaken or remove a negative provenance test merely because the repaired
reader itself now resides at the new required commit.

## Before/after provenance behavior

Before:

`profile required source = 0bf1281...`

so a bundle collected from the current reviewed state sees the later E.8 reader
repair as an unexpected committed analysis change.

After:

`profile required source = 8985d9a...`

and a future profile/test commit may differ from that source only by the two
explicitly allowed profile/test files plus `docs/memory/`.

This preserves fail-closed provenance while admitting the exact reviewed
E.8.1.Fix.1 source.

## Deterministic local validation

Run at minimum:

```text
python -B -m py_compile testing/test_pion_hgcer_validation_bundle_profile_e8_1.py
python -m json.tool testing/pion_hgcer_validation_bundle_profile_e8_1.json
python -B -m unittest testing.test_pion_hgcer_validation_bundle_profile_e8_1 -v
```

Also run the generic collector regression module if it is locally deterministic:

```text
python -B -m unittest testing.test_collect_pion_hgcer_validation_bundle -v
```

If that module has environment-specific skips, report them exactly.

Run the normal memory manifest/integrity/health/bootstrap checks required by
current repository instructions and:

```text
git diff --check
```

Codex must report every command and exact result. These are source/provenance
checks only, not farm validation.

## Required diff/runtime-path audit

Before stopping:

1. confirm starting HEAD and current HEAD remain
   `8985d9a212799c021c4ad1a759a689ea3826e0ea`;
2. show `git status --short`;
3. show `git diff --stat`;
4. show `git diff --name-only`;
5. confirm every changed path is allowlisted by this contract;
6. verify no `src/` file changed;
7. verify profile settings/artifacts/allowlists are byte-for-byte unchanged
   except the exact required commit field;
8. verify the test changes only the expected source identity unless a narrowly
   necessary assertion update is required;
9. trace the intended farm path:
   future profile commit -> wrapper `--bundle-commit` -> detached worktree ->
   E.8.1 profile -> generic collector -> source-ancestry/range checks ->
   existing artifacts only.

If a diff is too large for practical terminal review, create one complete
temporary root-level `kaonlt_review.diff` containing the tracked diff plus full
no-index additions for every intended new file. Do not stage merely for review.

## Durable-memory update

Update `CURRENT.md` to record:

- live/pushed E.8.1.Fix.1 source commit:
  `8985d9a212799c021c4ad1a759a689ea3826e0ea`;
- E.8.1 remains `ACTIVE`;
- E.8.1.Fix.1 remains `SOURCE REVIEWED`, now user-pushed;
- E.8.1.Fix.2 bundle-profile re-pin is `ACTIVE` until ChatGPT independently
  reviews its actual diff;
- F.6.2 and F.6.2.Fix.5 remain `CLOSED / RUNTIME VALIDATED`;
- F.6.3 remains `BLOCKED`;
- no new farm/runtime acceptance is claimed;
- NEXT after source review is user commit/push of the profile re-pin, then a
  fresh `Q4p4W2p74 / Left / lowe` farm bundle using the newly pushed profile
  commit as `--bundle-commit`;
- only after Left/lowe passes should the gate return to detailed Left/highe and
  then broaden.

Update the existing E.8.1.Fix.1 phase record only to append the actual pushed
source commit and the fact that profile re-pin is now the next provenance step;
do not rewrite its accepted source-review history.

Create:

`docs/memory/phases/e8-1-fix2-bundle-profile-repin.md`

as the permanent chronology for this profile-only repair. Do not create a farm
evidence record yet.

## Farm boundary and future bundle command

Do not run the farm in this task.

Do not hard-code or guess the future profile-repin commit into any file. After
the user commits and pushes this task, ChatGPT will inspect the pushed commit.
That future pushed commit becomes the wrapper `--bundle-commit` value for the
next Left/lowe bundle-only collection.

The profile's `required_analysis_commit` remains the already known repaired
analysis source:

`8985d9a212799c021c4ad1a759a689ea3826e0ea`

These two identities must remain distinct.

## Forbidden shortcuts

Do not:

- broaden `allowed_committed_files`;
- allow `src/` paths after the required analysis commit;
- change artifact declarations;
- change canonical settings;
- change the collector or wrapper;
- change analysis source;
- modify accepted F.6.2 artifacts;
- bump the profile/version merely to make the provenance change look larger;
- insert the future profile commit as `required_analysis_commit`;
- commit, push, or run the farm.

## Acceptance criteria

Source-review candidate PASS requires:

- profile required source is exactly
  `8985d9a212799c021c4ad1a759a689ea3826e0ea`;
- profile/test are the only allowed committed files;
- canonical settings and artifact declarations are unchanged;
- source identity remains fail-closed;
- analysis-source changes after the required source remain rejected;
- generic collector and wrapper are untouched;
- no `src/` or production file changed;
- deterministic local tests pass to local capability;
- memory accurately distinguishes repaired analysis source from future
  profile/bundle commit;
- actual diff is narrow and allowlisted.

## Hard stop

Do not commit, push, or run Jefferson Lab farm validation.

Return to ChatGPT with:

- implementation summary;
- exact changed files;
- exact local test/check results;
- provenance/runtime-path trace;
- actual review diff artifact if needed;
- proposed user Git handoff commands, not executed.

ChatGPT must inspect the actual diff before the user commits or pushes.
