# KaonLT source and evidence index

## Authority order

For implementation, current source and an exact newer diff outrank handoffs,
history, and inference. For runtime, fresh farm evidence and a reviewed farm
handoff outrank source. Git subjects, local tests, and collectors alone never
establish runtime validation.

## Live repository identity

- Repository: trottar/lt_analysis; branch: test.
- Live migration head: b02316f83bf8ef18641fa7217b90c04bcdca10e3.
- Fix.5 analysis-source commit: dc4fc6283001739a487ec80068f951b0e388cae6.
- Import-preparation observation: 7cdb7847d72501cd3dc504565fafc604e22a6132.

The current head must be rechecked before new work. The latter two identities
are source-history observations, not runtime claims.

## Repository-owned records

- AGENTS.md — memory, ownership, and runtime-evidence protocol.
- CURRENT.md, MEMORY.md, handoffs/CURRENT_HANDOFF.md — concise active state.
- import/chatgpt_project_history_2026-09-11/ — immutable imported historical
  source provenance for this migration.
- decisions/DECISION_HISTORY.md, phases/PHASE_HISTORY.md, and
  evidence/VALIDATION_HISTORY.md — reconciled durable records.

## Relevant implementation paths

Production owners include src/cuts/rand_sub.py,
src/cuts/proton_contamination_weights.py, src/cuts/pion_component_fits.py,
src/cuts/pion_component_subtraction.py, src/binning/calculate_yield.py,
src/binning/ave_per_bin.py, and src/utility/background_config.py.

The HGCer/frozen-record path includes pion_hgcer_refinement_method_a.py,
pion_hgcer_refinement_method_b.py, pion_hgcer_refinement_checkpoint.py,
pion_hgcer_refinement_comparison.py,
pion_hgcer_method_a_acceptance_contract.py, rand_sub.py, and
full_background_subtraction_plots.py under src/cuts/. The collector/profile
remain validation tooling and are not farm evidence.

## Historical external references

The Trotta and Usman dissertations, particle-background procedure documents,
Phase-C and F.1 continuation handoffs, and Lambda/timing repair contracts are
scientific or historical authorities for their matching decision/evidence
records. They are not current source authority. The import SOURCE_INDEX.md
preserves original titles and identifiers.

