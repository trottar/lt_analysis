# KaonLT — October 18 low-ε baseline-science memory checkpoint

**Task:** one bounded **documentation-only** update. **Status:** `ACTIVE` baseline-pion scientific investigation. **Source reference:** `test` at `f1d1984b5e0fca9ba3c5666477996af783423db7` when this contract was authored; verify current local HEAD and the entire worktree before editing.

## Objective and scope

Correct tracked `CURRENT.md`, which still directs a review of the already-published debug patch and a farm run that has already completed. Record the October 17 farm evidence, first scientific diagnosis, and the agreed next read-only science gate. **Do not rebuild or modify the initially useful staged `pi_n`/`pi_sidis`/`pi_delta` model.** Method A stays detached/non-production; Method B stays independent/diagnostic-only. This is not a new phase, full farm closure, analysis rerun, source implementation, or template/fit correction.

Read tracked root `AGENTS.md`, then in exact order `docs/memory/CURRENT.md`, `docs/memory/MEMORY.md`, `docs/memory/handoffs/CURRENT_HANDOFF.md`, `docs/memory/USER.md`; then `docs/memory/MAINTENANCE.md`, `docs/memory/CODEX.md`, the investigation below, the manifest and needed direct CURRENT references. Inspect local branch/HEAD, tracked/staged/untracked files and review bundles. Preserve all unrelated changes. If a true overlap/conflict exists, stop `BLOCKED` at the first precise invariant; do not reset/clean/stash or invent another contract.

**Only these repository paths may change:**

- `docs/memory/CURRENT.md` — narrow state and single NEXT update; correct obsolete pre-run assertions.
- `docs/memory/investigations/baseline-pion-control-center-lowe-regression-2026-10-14.md` — append one concise dated evidence section, retaining historical text.
- `docs/memory/manifest.json` — regenerate with existing tools.
- `docs/memory/phases/kaonlt-low-epsilon-baseline-science-memory-checkpoint-2026-10-18-task-contract.md` — this contract once transferred into the repo.

Do **not** edit `MEMORY.md` unless a truly absent reusable rule is demonstrated and explicitly escalated. No changes under `src/`, `testing/`, launcher, scripts, symlinks, models or other records. Codex does not commit/push or run farm/ROOT. Keep unrelated historical normal-launcher/attempt-isolation blockers separate from the successful plain `-d` case.

## Evidence to preserve, with strict labels

- **SOURCE VERIFIED:** the reviewed `-d` launcher at `f1d1984b5` uses ordinary available low-ε settings and exits before full high ε; plain `-d` does not invoke `-c` merge. That source change did not alter analysis physics. GitHub `test` HEAD was `f1d1984b5` at this checkpoint; farm checkout HEAD was also explicitly verified before the run.
- **RUNTIME VERIFIED:** farm invocation `./run_Prod_Analysis.sh -d 4p4 2p74` ended with `DEBUG_EXIT_STATUS=0`, reported all available low-ε settings complete and skipped high ε. Center and Left low-ε produced analysis and full-subtraction PDFs. There is no established populated Right/lowe result.
- **RUNTIME VERIFIED:** one Center/lowe canonical t1/t2/t3 seven-stage JSON/PDF pair was captured, PDF 19 pages. Producer status is **`captured_pending_independent_validation`**, completeness `captured`, PDF `available`; its `observed_source_HEAD` remains **`not_provided`**. Do **not** imply full independent provenance/checker acceptance or broader production scientific validation. The reviewed `KaonLT_loweps_20261018_093137.zip` has eight intact members: six Center/Left analysis, full-subtraction and pion-fit-debug PDFs plus that Center JSON/PDF pair. The user's desired future handoff is one ZIP per farm run, with original evidence and checksums.
- **RUNTIME VERIFIED, Center/t1:** the staged SIDIS normalization describes its local anchor well, but joint full-range refinement heavily suppresses SIDIS and underpredicts pion control around Λ; the staged fit itself overpredicts portions of the high-MM tail. This is a *transition-specific fit/model tension*, not proof that the optimizer is broken.
- **RUNTIME VERIFIED, Center/t2:** one high-leverage negative signed random contribution near 0.910–0.916 GeV contributes to a negative neutron least-squares numerator (~−61.674). Initial staged `B_n=0`, refined `B_n≈0.000202`, and the ratio reaches `w0≈182.82`. Center/t2 `K_model` pion integral ≈0.2529 versus event-applied background ≈6.5562. Model algebra/after-subtraction arithmetic close, but signed-event/model agreement fails in each Center t bin. Do not propose dropping the negative event or assert sole root cause.
- **RUNTIME VERIFIED, Left:** setting-wide neutron amplitude remains positive but fit shape quality is poor; its production canonical-t plots show concerns. Setting-wide model closure is **not** canonical-t event/model closure; matching source-resolved Left canonical-t diagnostics are **NOT VERIFIED**.
- **SOURCE VERIFIED:** `src/main.py` calls the pion-shape loader without `t_bins`/`phi_bins`, so canonical t scopes may fall back to setting-wide component shapes; the Center artifact reports no constructed t-binned SIMC shapes. This is a *possible limitation*, not a proven cause or permission to activate new shapes.
- **INFERENCE:** a statistical-support-gated fit per t followed by conditional parent-anchored `(t,φ)` refinement could be useful. It is **not** an accepted/implemented physics change. Parent-child event covariance, Sumw2, signed cancellations and available SIMC support must be examined first. Preserve the current staging, normalization, weights, cuts, priors and templates.

Do not mislabel full-range integrals as Λ cross sections or extracted yields; do not infer a historical first-bad commit. Accepted earlier phase closures and their original scopes remain untouched.

## State and NEXT

Replace the obsolete, multi-kilobyte pre-run NEXT in `CURRENT.md` with the **single substantive NEXT**:

> `NEXT`: Read-only audit of canonical-t and `(t,φ)` statistical/SIMC support for Q4p4W2p74 Center/Left low ε, using existing evidence first. Quantify signed/absolute counts, Sumw2 and cancellation, effective statistical support where meaningful, per-child SIMC occupancy, parent fit quality, component support and under-Λ residuals. Identify any actual missing evidence. Test the *justification* for hierarchical refinement; do not implement or rerun the farm from this memory checkpoint.

Make current-work-item text reflect *farm observed and evidence reviewed*, not “no capture/no farm.” Keep distinct producer capture versus checker/provenance acceptance. Preserve existing CLOSED statuses and external normal-launcher isolation blockers. `CURRENT.md` was ≈13.8 KiB (above 8 KiB soft warning; below 16 KiB hard limit); compact obsolete detail by referencing historical files. Do not increase above hard limit or spend a separate cycle on aesthetic memory compaction. The investigation remains `ACTIVE`.

## Deterministic checks and handoff

With a discovered suitable Python interpreter, run `tools/update_memory_manifest.py --root . --write` followed by `--check`, `tools/check_memory_health.py --root .`, and `tools/memory_bootstrap.py --root . --json`. Then `git diff --check`, allowlist/staged/untracked review, exact CURRENT/MEMORY/CURRENT_HANDOFF byte counts, sole-NEXT and frozen-status checks. Report nonblocking warnings; stop on a hard failure.

Prepare **one complete actual review diff**, including full untracked contract contents and the other changed memory files. Keep that diff untracked and outside publication. Do **not** stage, commit, push, run farm or alter science. After independent ChatGPT diff review and user publication, immediately return to the read-only statistical-support science audit. No additional task contract or farm-readiness project follows from this documentation checkpoint.
