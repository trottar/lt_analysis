# Production Analysis Workflow

## Overview

The production LT workflow keeps the existing split:

- `src/main.py` is the 0th-iteration, data-side authority.
- `src/main_iter.py` and `src/main_auto.py` reuse frozen data-side products and only update later SIMC-weight iterations.

The nominal physics procedure is unchanged when `BG_OPT_ACTIVE_PROFILE = "nominal_weighted"`.

## 0th Iteration

The 0th iteration is run through `src/main.py` for low epsilon and then high epsilon.

During that pass the workflow now:

- applies the usual cuts and subtraction chain;
- optimizes the shared `t` and `phi` binning at low epsilon;
- optimizes empirical background scales through the active background profile;
- freezes the selected bin edges and background prescription;
- writes a correction ledger for each epsilon pass;
- writes a 0th-iteration input bundle that records the exact `main.py` arguments;
- finalizes a frozen manifest after the high-epsilon pass;
- writes a high/low empirical-correction comparison;
- writes a review-oriented final analysis summary.

## What Is Optimized

The active profile in `src/utility/background_config.py` controls:

- selection mode;
- optimizer metric weights;
- forced Fit 1 and Fit 2 scales;
- whether high and low epsilon use independent or common scales;
- future profile-specific constraints.

Nominal mode remains:

- weighted optimizer;
- default metric weights;
- independent low/high epsilon scale optimization;
- current sidebands and fit models.

## What Is Frozen

After the high-epsilon 0th pass, the workflow writes:

- `OUTPATH/<ParticleType>_Q<Q2>W<W>_frozen_manifest.json`
- `OUTPATH/<ParticleType>_Q<Q2>W<W>_0th_iteration_input_bundle.json`

The manifest records:

- particle, polarity, `Q2`, `W`, and low/high epsilon values;
- MM cut window;
- frozen `t` and `phi` bin edges;
- background model keys, sidebands, and fit expressions;
- selected Fit 1 and Fit 2 scale maps;
- active optimizer profile and resolved settings;
- timestamp, git hash if available, and key file hashes.

## Later Iterations

`src/main_iter.py` and `src/main_auto.py` now validate the frozen contract before reuse.

They require:

- frozen manifest;
- previous JSON;
- previous ROOT file;
- correction ledgers;
- support NPZ files.

They also validate:

- current `Q2`/`W`;
- current `t` and `phi` bin edges;
- bin counts;
- tracked config/code hashes.

Config drift is controlled by `ALLOW_CONFIG_DRIFT` in `src/utility/background_config.py`.

## Correction Ledgers

Each epsilon pass writes:

- `OUTPATH/<ParticleType>_<OutFilename>_correction_ledger.json`
- `OUTPATH/<ParticleType>_<OutFilename>_correction_ledger.csv`

The ledger records Lambda-window yields after:

- raw prompt;
- random subtraction;
- dummy subtraction;
- pion subtraction;
- Fit 1;
- Fit 2;
- final subtraction.

It also records:

- fractional empirical corrections;
- stage-to-stage percent changes;
- empirical-fit uncertainty estimates;
- failed-fit and zero-background fallback counts;
- over-subtraction diagnostics;
- per-`(t,phi)` stage yields for later comparison tools.

## Epsilon Comparison

After the high-epsilon pass the workflow writes:

- `OUTPATH/<ParticleType>_Q<Q2>W<W>_epsilon_empirical_compare.json`
- `OUTPATH/<ParticleType>_Q<Q2>W<W>_epsilon_empirical_compare.csv`

This compares the empirical correction fraction between low and high epsilon for each `phi` setting and `t-phi` bin.

Thresholds are configured in `src/utility/background_config.py`:

- `EMP_EPS_DIFF_WARN_THRESHOLD`
- `EMP_EPS_DIFF_FAIL_THRESHOLD`

## Background Profiles

Profiles are selected only in `src/utility/background_config.py` through:

- `BG_OPT_ACTIVE_PROFILE`

Implemented profiles:

- `nominal_weighted`
- `ratio_first`
- `kinematic_diagnostic_only`
- `no_empirical_residual`
- `fit1_only`
- `fit2_only`
- `common_epsilon_scales`

Systematic replay lists are configured through:

- `BG_SYSTEMATIC_PROFILES`

## Systematic Replays

`src/utility/run_bg_systematics.py` replays the stored 0th-iteration inputs through each configured profile without permanently editing the tracked config file.

It uses profile snapshot overrides via environment-driven config snapshots and writes:

- `OUTPATH/<ParticleType>_Q<Q2>W<W>_bg_systematics_summary.json`
- `OUTPATH/<ParticleType>_Q<Q2>W<W>_bg_systematics_summary.csv`

Failed profiles are kept in the summary with an error message.

## Non-KLambda Cross-Checks

`src/utility/nonklambda_crosscheck.py` is an independent post-production interface.

It is not part of the production yield extraction.

It compares empirical residual corrections against externally supplied non-KLambda expectations and writes:

- `OUTPATH/<ParticleType>_Q<Q2>W<W>_nonklambda_crosscheck.json`
- `OUTPATH/<ParticleType>_Q<Q2>W<W>_nonklambda_crosscheck.csv`
- `OUTPATH/<ParticleType>_Q<Q2>W<W>_nonklambda_crosscheck.pdf`

## Final Summary

After the high-epsilon pass, and after later iteration high-epsilon updates, the workflow writes:

- `OUTPATH/<ParticleType>_Q<Q2>W<W>_final_analysis_summary.md`
- `OUTPATH/<ParticleType>_Q<Q2>W<W>_final_analysis_summary.csv`
- `OUTPATH/<ParticleType>_Q<Q2>W<W>_final_analysis_summary.json`

This summary is intended for note-ready review and collects:

- frozen configuration and profile metadata;
- selected scale maps;
- low/high empirical correction summaries;
- epsilon-difference diagnostics;
- cross-section outputs when present;
- optional systematic-envelope and non-KLambda status.

## Recommended Review Artifacts

For production review, inspect:

- frozen manifest JSON;
- low and high correction ledgers;
- BG optimization JSON/CSV/PDF;
- epsilon empirical compare JSON/CSV;
- final analysis summary Markdown;
- optional background-systematics summary;
- optional non-KLambda cross-check outputs.
