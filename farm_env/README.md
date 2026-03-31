# Farm Workflow Helpers

This directory contains the Python helpers wrapped by the repo-root script
`run_farm.sh`.

## What `run_farm.sh` Does

From the base `lt_analysis` directory, `run_farm.sh` manages two workflow
families for a `Q2/W` setting:

- replay submission
- applyCuts submission

It can also rebalance either workflow after jobs have been submitted.

The wrapper scans all matching manifests under `input/kaon` for the requested
family, so all supported phi/epsilon/target variants are included
automatically.

For isolated testing, you can point the wrapper at a separate manifest tree with
`-m`, for example `input/kaon_test`.

## Basic Usage

From the repo root:

```bash
./run_farm.sh Q2 W
```

Example:

```bash
./run_farm.sh 3p0 3p14
```

That prints a replay dry-run plan. Nothing is submitted unless `-s` is used.

## Flags

`-h`

- Print the help message and exit.

`-s`

- Actually submit jobs for the selected mode.
- Without `-s`, submit mode is dry-run only.

`-r`

- Switch to rebalance mode.
- In this mode the wrapper calls `farm_env/rebalance_swif.py` for the selected
  workflow name.

`-a`

- Only valid together with `-r`.
- Actually apply the `swif2 modify-jobs` changes.
- Without `-a`, rebalance mode is dry-run only.

`-n`

- Do not call `swif2 run` after submit or rebalance.

`-c`

- Use applyCuts mode instead of replay mode.
- Replay is still the default if `-c` is omitted.

`-w <workflow-name>`

- Override the auto-generated workflow name.

`-m <manifest-dir>`

- Override the manifest directory scanned by replay/applyCuts submission.
- Default is `input/kaon`.
- Useful for small test manifests stored separately from production.

## Replay Mode

Default mode:

```bash
./run_farm.sh 3p0 3p14
./run_farm.sh -s 3p0 3p14
./run_farm.sh -m input/kaon_test -s 3p0 3p14
```

This calls `farm_env/submit_replay.py`, which:

- finds all matching `Q{Q2}W{W}*.json` manifests under the selected manifest directory
- merges runs across those manifests
- submits one replay job per unique run

Default workflow name:

```text
kaonlt_Q{Q2}W{W}_${USER}
```

## applyCuts Mode

Use `-c`:

```bash
./run_farm.sh -c 3p0 3p14
./run_farm.sh -c -s 3p0 3p14
./run_farm.sh -c -m input/kaon_test -s 3p0 3p14
```

This calls `farm_env/submit_applycuts.py`, which:

- finds all matching manifests under the selected manifest directory
- keeps each manifest variant separate
- submits one job per manifest variant + run
- only plans a job if the full replay ROOT file already exists in the ltsep
  replay ROOT area
- skips the job if both skim outputs already exist in the ltsep skim ROOT area
- runs `applyCuts_Prod.sh` once per planned job, letting that script process
  both kaon and pion for the run

Default workflow name:

```text
kaonlt_Q{Q2}W{W}_applycuts_${USER}
```

## Rebalancing

Replay workflow:

```bash
./run_farm.sh -r 3p0 3p14
./run_farm.sh -r -a -n 3p0 3p14
```

applyCuts workflow:

```bash
./run_farm.sh -r -c 3p0 3p14
./run_farm.sh -r -c -a -n 3p0 3p14
```

## Jasmine Uploads

`farm_env/jasmine_put_from_manifest.py` now resolves volatile source roots from
ltsep at runtime instead of storing them in the manifests.

Replay upload:

```bash
python farm_env/jasmine_put_from_manifest.py center high 3p0 3p14 lh2 --product-kind replay
python farm_env/jasmine_put_from_manifest.py center high 3p0 3p14 lh2 --manifest-dir input/kaon_test --product-kind replay
```

Skim upload:

```bash
python farm_env/jasmine_put_from_manifest.py center high 3p0 3p14 lh2 --product-kind skim
python farm_env/jasmine_put_from_manifest.py center high 3p0 3p14 lh2 --manifest-dir input/kaon_test --product-kind skim
```

What it uses:

- replay source: the ltsep-resolved replay ROOT directory
- skim source: the ltsep-resolved skim ROOT directory
- replay tape destination: `destination` from the manifest
- skim tape destination: `skim_destination` from the manifest

Small files are still tar-grouped before `jput` submission.

## Helper Scripts

`submit_replay.py`

- Plans or submits one replay job per unique run.

`submit_applycuts.py`

- Plans or submits one applyCuts job per manifest variant + run.

`rebalance_swif.py`

- Inspects an existing SWIF2 workflow and proposes or applies resource bumps.

`jasmine_put_from_manifest.py`

- Packages replay or skim products for tape using the manifest plus ltsep path
  discovery.

## Supported Kinematics

- `Q2=5p5`, `W=3p02`
- `Q2=4p4`, `W=2p74`
- `Q2=3p0`, `W=3p14`
- `Q2=3p0`, `W=2p32`
- `Q2=2p1`, `W=2p95`
- `Q2=0p5`, `W=2p40`

## Notes

- Run everything from the repo root so `run_farm.sh` can find `farm_env/`,
  `input/kaon/`, and `applyCuts_Prod.sh`.
- Full submission and real file checks still depend on the farm environment,
  ltsep, and visible replay/skim storage trees.
