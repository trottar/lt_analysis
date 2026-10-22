# KaonLT — Matched-cohort A/B/C event evidence, one science-producing gate

**Task class:** One coherent, detached, science-producing recovery from the October 21 accepted-input mismatch. **Not** a new phase, source-authority reset, or production change.

**Audit anchor (verify live):** `test` at `532b661f73a827034edc269f85dc776b4170f574`.

**Scientific question:** In each Q4p4W2p74 Center/Left low-ε canonical t1/t2/t3 population, does transfer-variance regularization actually reduce *signed-event-applied* pion-background variance without worsening the physically relevant pion-control MM shape or causing unacceptable changes beneath K+Lambda?

## 1. Existing evidence and the one real blocker

The historical October 19–20 A/B JSON and 42-fit A/B/C feasibility result are accepted only for their **original input cohort**. The current four farm `Center/Left × Data/Dummy` ROOT files (all `Q4p4W2p74_lowe`) were independently observed with SHA-256 **and size** different from all four original accepted A/B pins. There is no established event-population equivalence. The first October 21 farm capture at `04f1ec403…` failed in a scalar optimizer prepass; the prepass guard was reviewed and published at `532b661f…`. Neither event capture nor C scientific validation succeeded. Do not repin an old cohort, waive strict historical validation, pretend a matching filename proves identical events, or retry the failed owner unchanged.

The original scientifically accepted A/B/C results and their exact archived hashes must remain immutable and attributable to their original cohort. Historic C is **not** a valid fitted candidate for the presently different input files without a newly matched fit.

## 2. Startup, audit and hard stops

Read tracked repository-root `AGENTS.md`, then `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, `docs/memory/USER.md` in order; then `docs/memory/MAINTENANCE.md`, `CODEX.md`, `TOOLS.md`, `COMMUNICATION.md`, the October 19 A/B template and owner contracts, October 20 C contract, both October 21 event/capture-fix contracts, and only their directly relevant source. Establish live `test`/`origin/test`, clean/dirty staged/unstaged/untracked state and current source identity; emit full health receipt. Stop on unexpected source, incompatible worktree, missing authoritative data source, or needed physics behavior change.

Audit **actual call order**, source identity and output destinations of:

- `src/cuts/rand_sub.py::observe_template_comparison` (already-opt-in A/B comparator);
- `src/binning/calculate_yield.py` and `testing/pion_local_t_event_validation_capture.py` (already-opt-in, post-cache event capture; prepass excluded);
- `testing/run_pion_t_template_comparison_gate.py` and `testing/run_pion_local_t_event_validation_gate.py` (existing detached owners and isolation/collector helpers);
- `testing/pion_local_t_transfer_variance_comparison.py` and `testing/pion_local_t_event_validation_replay.py` (offline C fit and event replay, both currently bound to historic accepted cohort).

Before implementing, prove that A/B frozen parents and signed source events can be captured from **one unmodified scientific run and the exact same four physical inputs** with an owner-managed, unique attempt output, preserving unrelated outputs and raw ROOT files. If safe orchestration is impossible within these boundaries, stop with the **first concrete invariant**; do not build auxiliary owner chains.

## 3. One source-coherent scientific acquisition, not a hash workaround

Build **one** tracked user-run owner using existing isolated disposable worktree and ltsep-overlay helpers. It may jointly enable the existing A/B template and event observer **only in that owned child**, with strict environment and path checks. Record the SHA-256 and byte size of the four actual backing ROOT inputs before and after the child; reject any in-flight identity change or mismatch between recorded A/B inputs and event packet source identities. The current cohort is a **new cohort**, with independently recorded identity, not a replacement of accepted October 19–20 A/B authority. Preserve the historical `ACCEPTED_INPUT_IDENTITIES` and default historical-mode checks unchanged. If a separate explicitly named matched-cohort capture mode is necessary, it must require the tracked owner, same-run A/B identity joins and exactly the reviewed source provenance; fail closed when these links are missing, inconsistent or incomplete. Do not add an unrestricted pin-bypass option.

Use a new attempt-owned destination namespace, never overwrite or delete original accepted A/B or failed October 21 receipts, ROOT inputs, model files, or unrelated output. Reuse existing collector where possible; one owned analysis run, verified Center+Left A/B JSON and PDFs, two Center+Left signed-event JSON packets, six canonical-t parent/source joins, source/selection and histogram closure, verified ZIP/manifest plus log/receipts. Distinguish the captured independent cohort from scientific acceptance. Respect farm `tcsh`, and do not run the farm or publish from Codex/ChatGPT.

Once new A/B and signed-event packets are collected, recompute **all seven** prospective C strengths using the same reviewed objective/optimizer on **these new A/B inputs**, not historic C values or fitted priors. Any new cohort-specific CLI flag must be explicit, validated against the exact source/manifest/parent content and leave original CLI/default/historical accepted validations intact. Do not modify the C objective, physical weights, protected kaon model, templates or priors. Preserve fixed `K_B` for B/C and distinct A convention. Use the newly generated C JSON and captured packets together for the detached event replay with exact hashes and full source/coordinate/parent/model join checks.

## 4. The physics result must be the deliverable

Produce an independently reviewable concise **six-parent scientific comparison** (structured JSON and a small detached PDF, not a many-page gate dump) containing, for A, B and each of seven C strengths:

1. Observed signed pion-control MM and A/B/C fitted MM with standardized residuals and meaningful fit diagnostics, for each Center/Left t1–t3.
2. Actual event-applied pion background MM compared to the fitted kaon-side pion background, with full-MM and magnified K+Lambda displays, using the same normalization/cuts and separate `nommcuts` vs `allcuts` populations.
3. Explicit signed prompt/random/dummy/dummy-random sums and quadratic Sumw2; support, large-weight tails, denominator zeros and positive/negative cancellation, by physical MM domain (0.7–1.3, 1.10–1.16, 1.105–1.135, 0.94–1.03 GeV), plus full normal-MM.
4. Demonstrated parent/child closure and actual Center/Left phi support; no inferred cross-setting equivalence. Display a signal-region effect with errors, and distinguish fitted model, applied subtraction and authoritative yield.
5. Report which lambda values are supported vs dominated by fit degradation, variance amplification or denominator gaps **without selecting/adopting any strength**. Mark held-out fitted-estimator evaluation, genuine bias truth and MSE selection `unavailable` if no training-only refits/truth cohort exists. A deterministic record fold alone is not independent validation.

Do not infer a cross-section, bias, physical purity, acceptance or final Lambda yield from histogram-only C risk values. An off-Lambda catastrophic weight must not be described as an observed large Lambda contamination. Present fresh numerical/visual evidence, including negative results; preserve the broader full-MM and eventual yield-science objective.

## 5. Minimal source and tests, no architectural detour

Prefer reusing existing source and helpers. Audit existence before adding any path. Conditional allowlist for the *narrowest coherent implementation*:

- `testing/run_pion_local_t_matched_cohort_gate.py` and `testing/test_run_pion_local_t_matched_cohort_gate.py` (only if an adapted existing owner cannot safely provide the one operation);
- `testing/pion_local_t_matched_cohort_gate_profile.json` (if required by the unchanged collector);
- `testing/pion_local_t_event_validation_capture.py` and `testing/test_pion_local_t_event_validation.py` (strict separate new-cohort authority and tests, never historical pin weakening);
- `testing/pion_local_t_transfer_variance_comparison.py`, `testing/test_pion_local_t_transfer_variance_comparison.py`, `testing/pion_local_t_event_validation_replay.py` (only minimal explicitly new-cohort CLI/provenance support and tests);
- `src/binning/calculate_yield.py` or `src/cuts/rand_sub.py` **only if a missing read-only dispatch prevents simultaneously using the existing observers**, with no physics/fit/selection changes;
- `docs/memory/CURRENT.md`, `docs/memory/investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md`, `docs/memory/manifest.json`, and this single contract.

No changes to production physics, optimizer, parent fitter, particle subtraction, SIMC, model/rate normalization, HGCer Method A/B, generic collector, trusted isolation helpers, or existing accepted/frozen archival artifacts. If the allowlist cannot support the requirement, return BLOCKED with the exact reason, not an improvised expansion.

Deterministic tests must cover failed old-cohort pins in historical mode; correct four same-run identity joins in explicitly new mode; rejection of mixed old/new cohort; wrong paths/parents/source/phi; optimizer prepass no import/I/O; real authority after cache population; no duplicate or overwritten outputs; positive/negative weights and Sumw2; owner failure order and cleanup; unchanged direct old A/B/C behavior; full normal-MM and Lambda-domain source closure. Compile affected modules and rerun existing relevant tests, `git diff --check`, memory manifest health. Report skips accurately. No fake ROOT/farm acceptance.

## 6. Workflow and stop conditions

This is one science-producing task, not a sequence of memory/provenance feature phases. Codex returns a complete actual diff including all new files, test receipts, and one exact input→owned child→verification→ZIP→detached science-result chain. ChatGPT reviews the **actual diff**. User alone commits and pushes, then ChatGPT synchronizes published HEAD and verifies complete farm-readiness. User alone runs one farm attempt after review, with fresh output names. Any first failed gate stops; no automatic retries, repinning or acceptance of partial outputs. If the same-run coherent A/B and event capture cannot be obtained, report that exact impediment and propose no auxiliary support loop.

**NEXT:** The next admissible science evidence is six matched-source parent A/B/C event-applied background/variance comparisons, with realistic Lambda-region plots—not another unmatched-input hash inventory or memory-only closure.
