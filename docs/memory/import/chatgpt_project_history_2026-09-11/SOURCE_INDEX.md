# KaonLT source and reference index

## S-GIT-001 — current repository

```text
https://github.com/trottar/lt_analysis/tree/test
```

Observed HEAD when import prepared:

```text
7cdb7847d72501cd3dc504565fafc604e22a6132
```

Authority: current implementation state, not automatic runtime authority.

## S-MEM-001 — `AGENTS.md`

Role: active Codex memory protocol. Requires loading current memory, relevant deep records, live-HEAD reconciliation, and memory updates as part of substantive work.

## S-MEM-002 — active memory

```text
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/handoffs/CURRENT_HANDOFF.md
docs/memory/roadmap/CURRENT.md
```

Caution: initial versions were intentionally conservative and predate this full historical migration.

## S-PLAN-C-001 — Phase-C Method-B contract

```text
KaonLT_Codex_Plan_Phase_C_Method_B_Local_Pion_Background_Closure.md
```

Key authority for:

- predecessor phases A/A.1/A.2/A.3/B/B.1;
- approved B/B.1 HEAD `18b06ec...`;
- Method-A/Method-B independence;
- Phase-A-record-only Method B;
- protected signal and parent-relative closure rules;
- no correction application.

## S-HANDOFF-C-001 — C.Fix.2.2 continuation

Document begins:

```text
KaonLT Continuation — Phase C.Fix.2.2 Complete;
Next Step Is Single Left-Low Farm Validation
```

Key facts: accepted `c7af5ada...`, one Left-low gate before C.Fix.2.3, persistence/serializer lesson, Codex phase/fix discipline.

## S-HANDOFF-E3-001 — E.3.Fix.2 Left-low PASS continuation

```text
KaonLT_Continuation_Prompt_After_E3_Fix2_LeftLow_PASS.md
```

Key facts:

- Phase-C five-setting PASS;
- adaptive B `DO NOT PROMOTE`;
- C.Fix.2.3 `9a66bc...`;
- pre-E.3 `e385365...`;
- E.3 contract;
- E.3.Fix.2 `eb1710f...`;
- JSON-driven collector/profile;
- Left-low farm PASS at `bf53dac...`.

## S-HANDOFF-F1-001 — Phase F.1.Fix.5 continuation

Title:

```text
KaonLT Continuation Prompt — Phase F.1.Fix.5 and Method-A Acceptance Mapping
```

Key facts:

- Method B diagnostic-only future role;
- Method A sole future numerical HGCer input;
- F.1→F.6 roadmap;
- F.1 lineage;
- Fix.3/Fix.4 farm failures;
- F.1 v1 population bug;
- exact training/application separation;
- six-file Fix.5 allowlist;
- v2 schema requirements.

## S-PROC-001 — production particle-background procedure

Documents:

```text
KaonLT_particle_background_subtraction_updated.tex
KaonLT Particle-Background Subtraction Procedure
```

Key authority for production ordering, slow-proton versus pion ownership, code owners, fallback semantics, and reproducibility records.

## S-PROC-002 — proton/pion/HGCer physics procedure

```text
KaonLT_Proton_Pion_HGCer_Refinement_Physics_Procedure.md
```

Scientific explanation of background separation and Method A/B.

Caution: its older long-term A+B correction statement is superseded by the later explicit Method-B diagnostic-only decision.

## S-PROT-001 — K-Lambda reference / offset-gating repair

```text
codex_prompt_repair_lambda_reference_and_offset_gating.md
```

Historical implementation contract for K-Lambda provenance, timing fallback hierarchy, stable-center behavior, and no hard 5-NPE cut.

## S-PROT-002 — Lambda-preservation gate

```text
codex_prompt_lambda_preservation_gate.md
```

Historical contract for the setting-wide 10% protected-signal gate.

## S-PROT-003 — timing/Lambda cleanup

```text
codex_prompt_final_timing_t_lambda_gate_cleanup.md
```

Historical contract for proposed/applied semantics, ROOT page-lifetime repair, gate readability, and regression reporting.

## S-PION-001 — dynamic pion/SIMC alignment

```text
codex_prompt_dynamic_pion_simc_alignment.md
```

Historical design contract for fixed-envelope candidate comparison, staged scans, and alignment payload.

## S-THESIS-001 — Trotta PhD dissertation

```text
PhDThesis.pdf
```

Scientific/experimental reference for KaonLT apparatus, kinematics, ltsep/Python analysis, PID, offsets, SIMC/data treatment, L/T cross sections, and HGCer inefficiency context.

Not current 2026 code authority.

## S-THESIS-002 — Usman PhD dissertation

```text
Usman,Ali_PhD_PHYS_Thesis_2025Fall.pdf
```

Complementary KaonLT/Hall-C reference for pion electroproduction, HMS/SHMS, HGCer/aerogel PID, target variables/`delta`, SIMC generation, detector-response limitations, SIMC resolution, and pion backgrounds.

Not current KaonLT kaon-production code authority.

# Current production/source paths to inspect when relevant

```text
src/cuts/rand_sub.py
src/cuts/proton_contamination_weights.py
src/cuts/pion_component_fits.py
src/cuts/pion_component_subtraction.py
src/binning/calculate_yield.py
src/binning/ave_per_bin.py
src/utility/background_config.py

src/cuts/pion_hgcer_refinement_method_a.py
src/cuts/pion_hgcer_refinement_method_b.py
src/cuts/pion_hgcer_refinement_checkpoint.py
src/cuts/pion_hgcer_refinement_comparison.py
src/cuts/full_background_subtraction_plots.py
src/cuts/pion_hgcer_method_a_acceptance_contract.py

testing/collect_pion_hgcer_validation_bundle.py
testing/pion_hgcer_validation_bundle_profile.json
```

# Git chronology anchors

Early pion work:

```text
395451f...
6cfdd47...
2a19bde...
9db45f2...
64ea04e...
```

SIMC/background implementation:

```text
7a758491...
e1640390...
ab43cac0...
8284f087...
277f269f...
```

HGCer/Phase C onward:

```text
18b06ec1...
c7af5ada...
9a66bc62...
e3853655...
eb1710f4...
bf53dac8...
b0f6868c...
b49dd1cd...
4fa198e4...
d656e157...
81b25890...
8b1ad5b7...
eb253046...
7cdb7847...
```

Use exact handoffs/source diffs to interpret these commits; subjects alone are not proof.
