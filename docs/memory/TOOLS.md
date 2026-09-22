# KaonLT tools and operational commands

## Repository identity

At substantial-work start, establish the actual checkout state:

```sh
git status --short --branch
git rev-parse HEAD
git branch --show-current
```

For a diff audit:

```sh
git diff --check
git diff --name-only
git diff --stat
```

These commands do not replace inspection of the actual diff.

## Memory controls

Discover a working Python interpreter and represent it as `<PYTHON>`. On the
intended Linux/JLab environment this will normally be an available interpreter
such as `python3`, but the executable must be discovered rather than assumed.

```sh
<PYTHON> -B tools/check_memory_health.py --root .
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/memory_bootstrap.py --root . --json
```

## Local source review

Use lightweight repository inspection as appropriate:

```sh
git diff
git show --stat --oneline HEAD
git status --short
```

Do not represent ROOT or farm runtime commands as locally validated behavior.

## JLab farm boundary

The established farm shell is `tcsh`; actual analysis runtime is farm-only.
The detailed packaging procedure belongs to
[farm-validation-bundle-procedure.md](decisions/farm-validation-bundle-procedure.md),
and farm delivery communication belongs to [COMMUNICATION.md](COMMUNICATION.md).

## Ownership pointers

- [AGENTS.md](AGENTS.md) — behavioral and scientific boundaries.
- [CURRENT.md](CURRENT.md) — active state.
- [USER.md](USER.md) — collaboration preferences.
- [COMMUNICATION.md](COMMUNICATION.md) — farm-delivery communication.
- [CODEX.md](CODEX.md) — source-changing Codex workflow.
- [MAINTENANCE.md](MAINTENANCE.md) — memory maintenance.

This is a concise operations reference, not active state, scientific evidence,
a farm validation-bundle procedure, or Windows transport infrastructure.
