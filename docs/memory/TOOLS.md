# KaonLT tools and operational commands

This concise reference owns generic Linux/POSIX repository operations and
durable Jefferson Lab environment facts. Command output never replaces reading
the actual source or diff.

## Repository identity

```sh
git status --short --branch
git rev-parse HEAD
git branch --show-current
```

## Diff and source inspection

```sh
git diff
git diff --check
git diff --name-only
git diff --stat
git show --stat --oneline HEAD
```

## Memory controls

Discover a working interpreter and represent it as `<PYTHON>`; do not assume a
workstation-specific executable.

```sh
<PYTHON> -B tools/check_memory_health.py --root .
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/memory_bootstrap.py --root . --json
```

## JLab environment

- Farm shell: `tcsh`.
- Repository and analysis root:
  `/group/c-kaonlt/USERS/trottar/lt_analysis`.
- Analysis artifact root:
  `/group/c-kaonlt/USERS/trottar/lt_analysis/OUTPUT/Analysis/KaonLT`.
- Transferable bundle root: `/volatile/hallc/c-kaonlt/trottar/globus`.

Farm runtime remains farm-only; do not represent ROOT/PyROOT or full-analysis
behavior as locally validated. A project-wide farm Python executable is not
assumed here.

## Provenance inspection

```sh
git status --short
git rev-parse HEAD
git show --stat --oneline <commit>
git diff <base>..<head> -- <path>
```

## Procedure ownership

- [farm-validation-bundle-procedure.md](decisions/farm-validation-bundle-procedure.md)
  owns detailed farm packaging.
- [COMMUNICATION.md](COMMUNICATION.md) owns farm request/return format.
- [CODEX.md](CODEX.md) owns source-changing workflow.
- [MAINTENANCE.md](MAINTENANCE.md) owns memory maintenance.

This record does not own scientific status, an active NEXT, delivery style, or
phase-specific farm commands.
