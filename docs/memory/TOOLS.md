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
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/check_memory_health.py --root .
<PYTHON> -B tools/memory_bootstrap.py --root . --json
```

Run manifest write then check in sequence after memory changes, followed by
ordinary task-final health as shown above. Reserve `--fail-on-warning` for
explicit milestone/zero-warning audits or materially blocking-warning cases.
Hard failures block; record nonblocking warnings and batch their correction.
Report exact CURRENT/MEMORY/CURRENT_HANDOFF byte counts and warnings after
substantial work.

## JLab environment

- Farm shell: `tcsh`.
- Repository and analysis root: `<KAONLT_REPO_ROOT>`.
- Analysis artifact root: `<KAONLT_ARTIFACT_ROOT>`.
- Transferable bundle root: `<KAONLT_BUNDLE_ROOT>`.

Concrete values are supplied by the active external environment/ChatGPT Project
configuration, not stored in public workflow memory. These are symbolic roles,
not guessed replacement paths. Never invent or substitute a generic path;
request a necessary value only when unavailable from active configuration.

Farm runtime remains farm-only; do not represent ROOT/PyROOT or full-analysis
behavior as locally validated. A project-wide farm Python executable is not
assumed here.

## E.8.4 tracked-owner diagnostics

The Left/lowe owner uses one timestamped output stem for the intended ZIP,
child analysis `.log`, and `-gate-status.json`. The log and status live in the
`<KAONLT_ARTIFACT_ROOT>`; the ZIP lives in `<KAONLT_BUNDLE_ROOT>`. Status schema
`e8_4_fix5_owner_gate_status/v1` records source, setting, expected ZIP/log paths,
UTC timestamps, stage, running/failed/success state, literal failure reason and
stage completion flags. It is an owner diagnostic outside the frozen bundle.
Owner failure prints the status path to stderr; success stdout is one ZIP path.

The child log contains only launcher output and is unchanged after summary
hashing. Inspect gate-status for later owner failures; its absence or a running
state does not establish successful delivery. Early detached source checks and
final detached collection both use the unchanged collector. This reference
supplies no farm execution authorization or command.

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
