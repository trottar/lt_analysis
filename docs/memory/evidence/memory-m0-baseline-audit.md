# Memory M0 baseline audit

## Identity

- Audit date: 2026-09-22.
- Repository and branch: `trottar/lt_analysis`, `test`.
- Exact starting HEAD: `d186943ce5aab19dd3de2428c76ba459efb9af9f`.
- Initial `git status --short --branch`: `## test...origin/test`; there were
  no tracked or untracked worktree entries.
- `git branch --show-current`: `test`.

The expected M0 starting identity was confirmed before any repository file was
created or changed.  This audit is source/structural memory work only; it
contains no farm, ROOT/PyROOT, full-analysis, or runtime-validation claim.

## Existing-tool baseline

The contract's literal `python3` invocations were attempted first.  On this
Windows checkout, `python3.exe` resolved to the inaccessible App Execution
Alias at `C:\Users\trott\AppData\Local\Microsoft\WindowsApps\python3.exe`.
Each command returned 1 with `ResourceUnavailable` before its script started.
The literal command result is therefore `FAIL` at the launcher; it is recorded
as a launcher-environment warning for M0 because the memory tools themselves
did not start.  M0 does not alter the launcher or the tools.

The installed interpreter
`C:\Users\trott\AppData\Local\Programs\Python\Python312\python.exe`
then executed the unmodified scripts with `-B` (to avoid creating bytecode).
The pre-existing `tools/__pycache__/` directory was timestamped 2026-09-18 and
was not changed or removed by M0.

| Required command | Literal `python3` result | Executed unchanged-tool equivalent | Return code / result | Relevant observation | Validation class |
| --- | --- | --- | --- | --- | --- |
| `python3 tools/check_memory_health.py --root .` | 1, FAIL (launcher unavailable) | `...python.exe -B tools/check_memory_health.py --root .` | 0, PASS | `MEMORY HEALTH: PASS` | Structural/source only |
| `python3 tools/update_memory_manifest.py --root . --check` | 1, FAIL (launcher unavailable) | `...python.exe -B tools/update_memory_manifest.py --root . --check` | 0, PASS | `MANIFEST: PASS` | Integrity/source only |
| `python3 tools/memory_bootstrap.py --root . --json` | 1, FAIL (launcher unavailable) | `...python.exe -B tools/memory_bootstrap.py --root . --json` | 0, PASS | Reported `test`, `d186943ce5aab19dd3de2428c76ba459efb9af9f`, three active records, five CURRENT references, and health `pass` | Startup summary/source only |
| `python3 tools/check_memory_health.py --self-test` | 1, FAIL (launcher unavailable) | `...python.exe -B tools/check_memory_health.py --self-test` | 0, PASS | `SELF-TEST: PASS` | Temporary-fixture structural test only |
| `python3 tools/update_memory_manifest.py --self-test` | 1, FAIL (launcher unavailable) | `...python.exe -B tools/update_memory_manifest.py --self-test` | 0, PASS | `SELF-TEST: PASS` | Temporary-fixture integrity test only |
| `python3 tools/memory_bootstrap.py --self-test` | 1, FAIL (launcher unavailable) | `...python.exe -B tools/memory_bootstrap.py --self-test` | 0, PASS | `SELF-TEST: PASS` | Temporary-fixture startup-summary test only |

The successful bootstrap reported active-state metadata from `CURRENT.md`:
E.8.1 is `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`, with the
canonical-five-setting gate and `Left / highe` first.  It also reported the
live branch/HEAD dynamically.  No tool run is scientific, farm, or runtime
validation.

## Existing architecture inventory

| Existing record or directory | Current role observed in M0 |
| --- | --- |
| `AGENTS.md` | Repository-memory authority, startup, scientific boundaries, and actor limits. |
| `CODEX.md` | Codex workflow, task classes, source/farm evidence limits. |
| `COMMUNICATION.md` | Concise JLab farm delivery and communication rules. |
| `CURRENT.md` | Concise active objective, active work, identities, blockers, exact next action, and direct references. |
| `MAINTENANCE.md` | Record roles, consolidation procedure, size limits, and manifest/health procedure. |
| `MEMORY.md` | Durable rules plus substantial detailed phase, commit, artifact, and current-gate chronology. |
| `README.md` | Memory navigation and control-layer entry point. |
| `manifest.json` | Generated SHA-256/byte-count inventory with copied active-state metadata, generation date, and observed-Git-HEAD field. |
| `chats/` | Distilled chat-era navigation. |
| `decisions/` | Scientific, architectural, and implementation-contract decisions. |
| `evidence/` | Farm/runtime evidence records, plus a validation-history navigation summary. |
| `handoffs/` | Current continuation/transfer record. |
| `history/` | Historical and chat navigation material. |
| `import/` | Immutable 2026-09-11 imported history/provenance package and checksums. |
| `investigations/` | Durable active, closed, and deferred investigation records and known gaps. |
| `phases/` | Phase/fix implementation and validation summaries. |
| `roadmap/` | Approved Phase-F dependency structure, but also current frontier prose. |
| `sources/` | Source and artifact indexes with implementation-path and provenance guidance. |
| `templates/` | New-work contract, memory-update, and runtime-evidence templates. |

## Ownership-overlap findings

### CURRENT versus CURRENT_HANDOFF

`CURRENT.md` and `handoffs/CURRENT_HANDOFF.md` contain the same schema-2
frontmatter: active objective, work item, status, exact next action, baseline,
source, and bundle/profile commits.  Both also narrate the E.8.1 work item,
source identities, farm-only boundary, and downstream blockers.  The handoff
adds resume and transfer wording, but presently duplicates ordinary active
state.  M3's target is `CURRENT.md` as sole active-state authority and an
exceptional-only handoff; M0 makes no change.

### CURRENT versus roadmap/CURRENT

`roadmap/CURRENT.md` duplicates the same schema-2 frontmatter and repeats the
active E.8.1 objective, state, canonical-five-setting next action, blockers,
and reviewed source/bundle identities found in `CURRENT.md`.  It also contains
the legitimate Phase-F dependency structure.  M6's target is roadmap ownership
of dependencies/status structure, with today's authoritative frontier and
exact action left to `CURRENT.md`.

### MEMORY versus phase/evidence chronology

`MEMORY.md` contains durable cross-phase rules: evidence precedence, distinct
source/runtime identities, scientific ownership, Method-A/Method-B boundaries,
and farm-bundle procedure.  It also carries detailed F.1 through F.6.2
chronology, individual commits, bundle hashes, Fix.5 presentation history, and
the E.8.1/F.6.3 current gate.  The canonical deeper owners already exist:
F.1-F.6 implementation summaries in `phases/phase-f*-*.md`; accepted gate
facts in `evidence/*runtime-closure.md`; scientific/architectural policy in
`decisions/`; and current blockers in `CURRENT.md`.  M4 will normalize this
without losing unique provenance.

### Static versus live identity

The reviewed scientific, renderer, bundle/profile, and farm-evaluated commits
in phase/evidence records are historical/provenance identities and must remain
fixed with their roles.  In contrast, `sources/SOURCE_INDEX.md` labels the old
`b02316f83bf8ef18641fa7217b90c04bcdca10e3` as a "Live migration head", while
the actual M0 HEAD is `d186943ce5aab19dd3de2428c76ba459efb9af9f`.  The current
manifest also stores `observed_git_head` as
`c25f9d8248f2acfb8b0443a482fc76781eb405ca`.  The manifest check passed because
it validates that field's type, not equality with the dynamic HEAD.  M6/M8 will
make dynamic repository identity explicit without rewriting scientific
provenance.

### Control-tool representation

`check_memory_health.py` requires the control files and checks links, exact
CURRENT section order, runtime-claim evidence links, top-level role uniqueness,
authority markers, active-file size limits, resolvable commit syntax, and a
fresh manifest.  It parses all three active surfaces and requires their
frontmatter dictionaries to be identical; it therefore currently enforces
mirrored active-state metadata.

`update_memory_manifest.py` builds a non-authoritative schema-2 inventory of
versionable memory files excluding itself, with byte counts and SHA-256 values.
It copies `CURRENT.md` active state and records generation date plus observed
HEAD.  Its check compares files and copied active state, but does not assert
that `observed_git_head` equals the live checkout.

`memory_bootstrap.py` lists `CURRENT.md`, `MEMORY.md`, and the handoff as
active records; obtains active state from `CURRENT.md`; reports direct resolved
CURRENT references and live Git branch/HEAD; and runs health.  It is a narrow
summary, not a substitute for reading its records.

## Migration map

| Existing information | Existing owner(s) | Problem | Future canonical owner | Phase | Evidence/provenance to preserve |
| --- | --- | --- | --- | --- | --- |
| Ordinary active objective, work item, status, exact next action | CURRENT, handoff, roadmap, manifest frontmatter | Multiple mirrored authorities | CURRENT | M3/M8 | Exact active state and historical transition context |
| Exceptional transfer/resume context | CURRENT, handoff | Handoff mixes transfer facts with ordinary state | handoffs/CURRENT_HANDOFF | M3 | Resume-specific facts and transfer rationale |
| Stable rules and scientific boundaries | MEMORY, decisions, AGENTS | Some duplication is useful but roles are not yet sharply separated | MEMORY and decisions/AGENTS by rule type | M4/M5 | Method boundaries, evidence precedence, actor boundary |
| Detailed phase implementation history | MEMORY, phases, roadmap | Chronology inflates durable memory and repeats roadmap/current state | phases | M4/M6 | Contracts, source-review identities, phase scope |
| Farm/runtime closure facts and hashes | MEMORY, evidence, phases, roadmap | Acceptance evidence is repeated outside its canonical gate record | evidence | M4/M6 | All supplied artifacts, hashes, farm identities, limitations |
| Architectural/scientific decisions | decisions, MEMORY, phases | Decision rationale can be duplicated as phase chronology | decisions | M4 | Decision IDs, rationale, supersession and boundaries |
| Dynamic branch/HEAD/worktree observations | source index, manifest, bootstrap, handoff | Stored values can look permanently live | Dynamic startup query | M6/M8 | Historical source/evaluated/provenance SHAs with explicit roles |
| Source/artifact navigation | sources, evidence, import | Index labels mix navigation and live-state language | sources | M6 | Immutable import package and artifact limitations |
| Dependency/status structure | roadmap, CURRENT, phases | Roadmap repeats exact current action and source identity | roadmap | M6 | F.6 dependency ordering and non-promotion boundary |
| Operational execution guidance | AGENTS, CODEX, COMMUNICATION, MAINTENANCE, templates | Boundaries overlap but no separate user/tools/learnings roles yet | AGENTS, USER, TOOLS, COMMUNICATION, CODEX, MAINTENANCE | M2/M5 | Linux/JLab farm procedure and Codex/user/farm sequence |
| Historical/chat-era material | chats, history, import | Can be mistaken for active startup control | history/import/chats navigation | M6 | Imported source, checksums, unresolved-history limits |
| Integrity metadata and startup summary | manifest, health checker, bootstrap | Manifest copies semantic state and stale observed identity; bootstrap lists only three active files | manifest integrity only; bootstrap dynamic facts | M7/M8 | SHA-256 inventory and structural test coverage |

## Scientific frontier preservation

M0 did not change E.8.1; the F.6.2 or F.6.2.Fix.5 closures; F.6.3/F.6.4
ordering; Method A; Method B; analysis source; validation profiles; random,
proton, or pion subtraction; SIMC treatment; yields; or cross sections.  Method
B remains diagnostic/cross-check only.  Method A remains outside production
unless a later explicitly validated F.6 promotion authorizes it.

No farm command was run, no commit/push or remote update occurred, and M0 makes
no farm/runtime acceptance assertion.

## Unknowns and deferred questions

- M0 does not decide the schema-v3 field set or the exact legacy
  compatibility period; that is M1 work.
- M0 does not decide which current handoff facts are genuinely exceptional
  transfer state; M3 must classify them before removing any mirrored content.
- M0 does not collapse MEMORY chronology or move unique provenance; M4 must
  trace each fact to its canonical deeper record first.
- M0 does not decide final role-specific size limits, duplicate-role rules, or
  rendering checks; those require M7's contract and deterministic tests.
- M0 records the launcher mismatch but does not add Windows-specific
  infrastructure; Linux/JLab operational behavior belongs to M5.

## M0 conclusion

`SOURCE REVIEWED` — the M0 baseline audit and migration map are complete from
the unmodified starting state, successful local structural tool executions,
and the actual repository inventory.  The `python3` App Execution Alias
failure is preserved as an environment warning; the tools themselves passed
through the installed interpreter.  This status is not farm/runtime validation.

## NEXT

NEXT — Memory M1: schema-v3 contract and compatibility health checker.
