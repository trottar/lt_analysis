# Memory M9 fresh-session integration audit

## Scope and evidence class

This is repository-memory integration evidence only. It records a read-only
fresh-session recovery audit and deterministic local memory checks. It is not
Jefferson Lab farm, ROOT/PyROOT, full-main, scientific, or runtime validation.

## Starting observation

- Branch: test.
- HEAD: d22b61228841d39cd5d815a3b094cf115b388eb5.
- Worktree: clean (## test...origin/test, no short-status entries).
- The observed HEAD equalled the M9 audit baseline. No reset, stash, clean, or
  unrelated modification was performed.

## Fresh-session boundary

After the M9 audit prompt was read, no prior-chat transcript, external
continuation, or external KaonLT summary was deliberately consulted or used as
audit evidence. A brief native-memory lookup occurred before the M9 prompt was
read; it was not audit evidence and is explicitly non-disqualifying under the
audit contract. Every substantive result below was independently recovered
from the repository records listed here after the audit clock began.

## Startup recovery

docs/memory/AGENTS.md supplied the tracked startup contract. After dynamic Git
identity was established, its five-file core was read completely in this exact
order:

1. docs/memory/AGENTS.md
2. docs/memory/CURRENT.md
3. docs/memory/MEMORY.md
4. docs/memory/handoffs/CURRENT_HANDOFF.md
5. docs/memory/USER.md

The contract assigns schema-3 CURRENT.md sole ordinary active-state authority,
the handoff exceptional-transfer-only authority, and roadmap/STATUS.md
dependency/status ownership rather than exact active action.

## Stage-A recovery table

| Recovery item | Recovered result | Canonical source(s) | Confidence / ambiguity |
| --- | --- | --- | --- |
| A. Repository identity semantics | Establish branch, HEAD, and worktree dynamically. Current source/diff and direct applicable farm evidence outrank durable memory; historical/chat material is lower authority. Stored SHAs are provenance roles, never permanent live HEAD. | AGENTS.md; MEMORY.md; decisions/memory-schema-v3-contract.md; phases/memory-m8-manifest-bootstrap-semantics.md | High; no ambiguity. |
| B. Active-state ownership | CURRENT.md solely owns active objective, work item, blockers, and exact NEXT. Roadmap owns approved dependency/status structure; handoff owns exceptional transfer only. Before closure the active objective/work item/NEXT were M0-M9 refinement, Memory M9, and this audit. | CURRENT.md; MEMORY.md; handoffs/CURRENT_HANDOFF.md; roadmap/STATUS.md | High; no ambiguity. |
| C. Scientific frontier | F.1 through F.6.2, including F.6.2.Fix.5 presentation closure, are CLOSED / RUNTIME VALIDATED under their direct evidence owners. E.8.1 is DEVELOPMENT COMPLETE, FARM VALIDATION PENDING; F.6.3 is BLOCKED pending E.8.1; F.6.4 is BLOCKED pending F.6.3 evidence. | roadmap/STATUS.md; evidence/f6-2-scientific-runtime-closure.md; evidence/f6-2-fix5-presentation-runtime-closure.md; phases/phase-f6-method-a-production-promotion.md | High. Older E.8 S1 phase wording is superseded for present state by CURRENT and roadmap ownership. |
| D. Scientific ownership boundaries | Production ordering is random/dummy -> frozen binning -> slow proton -> pion subtraction. Method B remains diagnostic/cross-check only and never changes production pion weights. Method A remains detached/non-production; only explicit validated F.6.4 may decide promotion. Presentation cannot recompute diagnostics or construct corrections. | MEMORY.md; roadmap/STATUS.md; phases/phase-f6-method-a-production-promotion.md | High; no ambiguity. |
| E. Source review versus farm evidence | Local/source review and deterministic checks do not establish farm, ROOT/PyROOT, full-main, or procedure-PDF runtime behavior. A completed bundle alone is insufficient; inspect provenance, checker gates, payload/log evidence, and rendered pages. If a checker conflicts with raw evidence, inspect raw evidence and implementation. | MEMORY.md; USER.md; COMMUNICATION.md; LEARNINGS.md | High; no ambiguity. |
| F. Actor sequence and workflow | Codex may make allowlisted local changes and deterministic checks. ChatGPT audits the actual diff and later fresh evidence. The user alone commits/pushes accepted work and runs farm validation. Codex may not commit, push, update remote refs, or initiate farm execution. | AGENTS.md; CODEX.md; COMMUNICATION.md | High; no ambiguity. |
| G. Linux/JLab operational context | The runtime target is the Jefferson Lab Linux farm using tcsh. Repository/artifact/bundle roots are /group/c-kaonlt/USERS/trottar/lt_analysis, /group/c-kaonlt/USERS/trottar/lt_analysis/OUTPUT/Analysis/KaonLT, and /volatile/hallc/c-kaonlt/trottar/globus. Detailed packaging belongs to the farm validation-bundle procedure. | TOOLS.md; decisions/farm-validation-bundle-procedure.md; USER.md | High; no ambiguity. |
| H. Handoff state | No exceptional transfer state is present. The handoff cannot override CURRENT. A new session resumes from CURRENT, follows the five-file core, and then opens task-relevant canonical records. | handoffs/CURRENT_HANDOFF.md; AGENTS.md | High; no ambiguity. |
| I. Selective expansion | The five-file core plus only schema/M8, F.6.2 evidence, F.6 phase, roadmap, operational/workflow records, farm procedure, E.8 contract, and learnings were sufficient; no whole-tree content load was needed. | Records listed below | High; no ambiguity. |

## Selective records opened

After the five-file core, the following additional records were opened because
they directly answer a recovery item: decisions/memory-schema-v3-contract.md
and phases/memory-m8-manifest-bootstrap-semantics.md (identity/bootstrap);
evidence/f6-2-scientific-runtime-closure.md,
evidence/f6-2-fix5-presentation-runtime-closure.md,
phases/phase-f6-method-a-production-promotion.md, and roadmap/STATUS.md
(frontier, dependency, and Method-A boundary); TOOLS.md, CODEX.md,
COMMUNICATION.md, and decisions/farm-validation-bundle-procedure.md
(actors and Linux/JLab operation);
decisions/e8-f6-2-figure-library-implementation-contract.md (historical E.8
scope); and LEARNINGS.md (checker/raw-evidence principle).

The complete canonical repository set consulted was the five-file core plus
exactly the records named above. No imported chat history, native-memory
record, or external summary is an evidentiary source for this audit.

## Bootstrap cross-check

Using C:\Users\trott\AppData\Local\Programs\Python\Python312\python.exe:

    <PYTHON> -B tools/memory_bootstrap.py --root . --json

returned dynamic branch test, HEAD d22b61228841d39cd5d815a3b094cf115b388eb5,
a clean worktree, ordered core records and byte sizes, CURRENT direct
references, memory-health pass (return code 0), and
exceptional_handoff.present = false. It did not emit active-state or
scientific-state synthesis, as required. This agreed with the manual
repository-led recovery.

## Verdict

Stage A passed. The repository supplied every required A-I recovery item,
canonical ownership resolved the older phase chronology without ambiguity, and
task-directed expansion was sufficient. Stage B records the audit without
changing science, runtime state, production behavior, source history, farm
state, profiles, collectors, or lifecycle hooks.

NO CHAT CONTINUATION REQUIRED
