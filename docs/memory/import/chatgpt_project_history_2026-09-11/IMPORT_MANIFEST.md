# KaonLT historical-memory import manifest

Import ID:

```text
kaonlt_chatgpt_project_history_2026-09-11/v1
```

Target:

```text
repository: trottar/lt_analysis
branch: test
```

Observed source identity when this package was prepared:

```text
test HEAD:
7cdb7847d72501cd3dc504565fafc604e22a6132

subject:
Pion subtraction hgcer checker| Phase F.1.fix.5 || Added working memory
```

This is a source-identity observation, not a runtime-validation claim.

## Files

```text
README.md
IMPORT_MANIFEST.md
PROJECT_HISTORY.md
CHAT_INDEX.md
PHASE_HISTORY.md
VALIDATION_HISTORY.md
DECISION_HISTORY.md
SOURCE_INDEX.md
ARTIFACT_INDEX_SEED.md
KNOWN_GAPS.md
CODEX_MIGRATION_CONTRACT.md
CHECKSUMS.sha256
```

## Provenance classes

- `[CURRENT_SOURCE]` — observed directly from current repository/source.
- `[GIT_HISTORY]` — Git commit/ancestry metadata. Commit subjects are chronology, not validation.
- `[CHAT_HISTORY]` — substantive KaonLT project discussion.
- `[HANDOFF]` — dedicated continuation/handoff preserving accepted state.
- `[FARM_EVIDENCE]` — supplied/reviewed farm output, checker/bundle, rendered-PDF evidence, or an authoritative handoff explicitly recording that review.
- `[SCIENTIFIC_REFERENCE]` — thesis/analysis-note scientific context.
- `[INFERENCE]` — derived conclusion not directly established by one authoritative record.
- `[SUPERSEDED]` — historically real plan/decision replaced by newer specific evidence.

For implementation state:

```text
current source
> newer exact source/diff evidence
> newest authoritative handoff
> older history
> inference
```

For runtime state:

```text
fresh farm evidence for the evaluated commit/setting
> reviewed authoritative farm handoff
> source evidence
> historical summary
> inference
```

Never assign `CLOSED / RUNTIME VALIDATED` from a commit subject, source existence, unit tests, or a collector alone.

## Migration policy

Do not paste this package wholesale into `MEMORY.md`.

Distribute it into deep memory:

```text
docs/memory/history/
docs/memory/phases/
docs/memory/decisions/
docs/memory/evidence/
docs/memory/sources/
docs/memory/investigations/
```

Then regenerate/reconcile only the concise active layer:

```text
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/handoffs/CURRENT_HANDOFF.md
docs/memory/roadmap/CURRENT.md
```

## Non-goals

This import does not:

- alter analysis code or production physics;
- promote Method A or Method B;
- validate Phase F.1.Fix.5 on the farm;
- infer missing farm results from Git;
- replace the JSON-driven validation workflow;
- claim literal preservation of every old ChatGPT utterance.
