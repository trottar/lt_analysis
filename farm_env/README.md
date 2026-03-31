# Farm Workflow Helpers

This directory contains the Python helpers used by the repo-root wrapper
`run_farm.sh`.

`run_farm.sh` lives in the base `lt_analysis` directory and is the main entry
point for farm replay workflow work. It wraps:

- `farm_env/submit_replay.py`
- `farm_env/rebalance_swif.py`

## What `run_farm.sh` Does

Given `Q2` and `W`, the wrapper can do one of two things:

1. Submit replay jobs for that kinematic family.
2. Inspect and rebalance an existing SWIF2 workflow for that family.

In submit mode, it automatically scans all matching JSON manifests under
`input/kaon` for the requested `Q2/W` family. That means all matching settings
for that family are included automatically, including dummy variants when those
JSONs exist.

Dry-run is the default in both modes.

## Basic Usage

From the repo root:

```bash
./run_farm.sh Q2 W
```

Example:

```bash
./run_farm.sh 3p0 3p14
```

That prints the replay submission plan without actually submitting anything.

## Flags

`-h`

- Print the help message and exit.

`-s`

- Actually submit jobs via `farm_env/submit_replay.py`.
- Without `-s`, submit mode is dry-run only.

`-r`

- Switch to rebalance mode.
- In this mode the wrapper calls `farm_env/rebalance_swif.py` for the workflow
  derived from `Q2`, `W`, and `USER`.

`-a`

- Only valid together with `-r`.
- Actually apply the `swif2 modify-jobs` changes.
- Without `-a`, rebalance mode is dry-run only.

`-n`

- Do not call `swif2 run` after submit or rebalance.

`-w <workflow-name>`

- Override the auto-generated workflow name.

## Common Examples

Dry-run submit:

```bash
./run_farm.sh 3p0 3p14
```

Actually submit jobs:

```bash
./run_farm.sh -s 3p0 3p14
```

Dry-run rebalance of an existing workflow:

```bash
./run_farm.sh -r 3p0 3p14
```

Apply rebalance changes but do not rerun the workflow:

```bash
./run_farm.sh -r -a -n 3p0 3p14
```

Use a custom workflow name:

```bash
./run_farm.sh -s -w my_custom_workflow 3p0 3p14
```

## Supported Kinematics

- `Q2=5p5`, `W=3p02`
- `Q2=4p4`, `W=2p74`
- `Q2=3p0`, `W=3p14`
- `Q2=3p0`, `W=2p32`
- `Q2=2p1`, `W=2p95`
- `Q2=0p5`, `W=2p40`

## Helper Scripts

`submit_replay.py`

- Finds all JSON manifests matching `Q{Q2}W{W}*.json` under `input/kaon`.
- Merges run lists across those variants.
- Plans or submits one SWIF2 job per unique run.

`rebalance_swif.py`

- Looks for problem jobs in an existing SWIF2 workflow.
- Detects likely memory, disk, or time-limit failures.
- Plans or applies `swif2 modify-jobs` resource bumps.

## Notes

- The wrapper expects to be run from the repo root so that `run_farm.sh` can
  find `farm_env/` and `input/kaon/`.
- `submit_replay.py` and `rebalance_swif.py` can still be run directly, but
  `run_farm.sh` is the intended convenience entrypoint.
- Full end-to-end execution still depends on the farm environment, available
  run lists, and `swif2`.
