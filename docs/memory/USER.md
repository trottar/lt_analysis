# KaonLT user collaboration context

## Project working context

KaonLT is a mature Jefferson Lab Hall C analysis. Do not restart its
architecture or reopen closed phases without concrete new evidence. Its
operational/runtime target is the Jefferson Lab Linux farm; ROOT/PyROOT, full
`main.py` behavior, and production procedure PDFs require farm evidence rather
than local-source assumptions. Discover the actual local interpreter and OS
for each task rather than treating workstation configuration as project truth.

## Collaboration preferences

Use Work as the preferred environment for substantial offline KaonLT scientific development, numerical experiments, statistical validation, prototyping, and genuine-data diagnostic plots. Establish scientific feasibility and test failure modes before requesting farm-oriented integration. Use Codex for narrowly scoped tracked implementation changes when Jefferson Lab farm integration becomes necessary. Preserve independent source review, user-only Git publication, user-only farm execution, and applicable runtime validation. This separation is intended to minimize repeated costly farm attempts and speculative infrastructure development; offline tests never substitute for actual farm evidence.

For forthcoming E.8 quick-debug science, use **Left/lowe + Center/lowe**
together by default: adjacent low-epsilon runs are quick to analyze. Left/lowe
has restricted phi support; Center/lowe covers a wider phi range. Report actual
populated phi coverage and within-setting results separately; compare overlap
with acceptance and normalization limits explicit. Preserve earlier Left-only
conclusions as setting-local. This preference implies neither complete physical
acceptance, equal cross-setting normalization, pooled event populations nor
other-setting/high-epsilon analysis; preserve existing cuts and sources.

Use one narrow phase or fix at a time. Inspect actual source and diffs rather
than trusting summaries. Preserve frozen scientific interfaces unless the task
contract explicitly owns a change. Do not request information already present
in repository memory, supplied artifacts, or current source. Prefer the
smallest complete validation-artifact package needed for the question, and keep
source review explicitly distinct from farm validation.

Do not provide or authorize farm-run commands until the exact validation
bundle/profile required for that gate is ready, independently source-reviewed,
pushed, and pushed-state reviewed. Advance one workflow gate at a time; do not
provide commands for later gates before the current gate has passed.
Before a farm command, verify the complete tracked, reviewed, pushed, and
pushed-state-reviewed path from input authority through production/checking,
packaging, invocation, and returned artifact. Multi-step operations require a
tracked orchestration owner; otherwise farm readiness is `BLOCKED` and no
interactive shell sequence is provided. A single reviewed CLI may own a
complete direct operation. Codex and ChatGPT each report memory health after
substantial gates. Hard failures and material active-state/provenance ambiguity
block progress; record nonblocking warnings and batch correction at the next
checkpoint or milestone audit.

Repository memory exists to improve scientific throughput and prevent
regression; maintenance must not displace substantive physics work. The normal
loop brackets science and validation with memory checkpoints. ChatGPT
actual-diff review -> user commit/push -> ChatGPT pushed-state synchronization
review remains required so Codex, ChatGPT, and the farm use the same source.
Pushed-state review verifies identity, gate-relevant changed paths/blobs, and
CURRENT/NEXT continuity. Once synchronization passes and active memory is
materially accurate, proceed directly to the substantive gate; a push alone
requires no new source-change contract or broad memory reconciliation. Favor
early narrow runtime tests over prolonged speculative hardening. Record
nonblocking drift while continuing the active milestone, then batch cleanup
at meaningful milestones; standalone repairs require concrete active-state
ambiguity. See the [workflow decision](decisions/memory-bracketed-scientific-throughput.md).

After Codex completes a local source-changing task, ChatGPT must inspect the
actual diff before the user commits or pushes. If the diff is too large to
review comfortably in terminal output, Codex must provide an exact command that
writes one complete review bundle directly in the `lt_analysis` repository
root, using a clearly temporary filename such as `kaonlt_review.diff`, for the
user to upload. The review file must be removed before commit/push unless it is
explicitly intended to be tracked. Do not require manual terminal copy/paste or
staging merely to make a diff reviewable. A review bundle must include both the
tracked diff and complete `git diff --no-index /dev/null ...` representations
for every intended new/untracked file.

Tracked KaonLT analysis/procedure plots are implemented through the normal
source-changing Codex workflow and farm-rendered from repository source.
ChatGPT-generated/extracted plots are inspection aids only when explicitly
requested, never substitutes for tracked E.8 deliverables.

## Scientific output and correction requirements

The user's explicit October 22 correction identifies repeated assistant
farm-readiness approvals without representative environmental tests as an
assistant judgment failure, not a user farm-operation failure. Prioritize fresh
scientific throughput and concrete accountability. Do not repeat synthetic-only
farm-ready claims, reflexively rerun a failed owner or create recursive support
contracts. Record the first exact blocker, cover the supplied failure classes
and adjacent guards in the existing bounded repair, and state the unknown
installed environment honestly. See the
[incident](investigations/kaonlt-matched-cohort-farm-readiness-failure-2026-10-22.md).

Scientific tasks must produce the requested evidence. The user does not accept
multi-day support-contract spirals, historical PDFs or illustrations in place of
requested fresh farm results. After two consecutive non-evidence support gates,
stop new support work and apply the MAINTENANCE science-gate displacement alarm.
Preserve pending local work; inspect available source and contracts before asking
for genuinely missing information. A user correction requires an immediate
professional re-anchor and the shortest justified action, not rationalization.
A comprehensive durable memory correction must update the appropriate governing
records as well as the investigation; a standalone incident write-up is incomplete.
For source-changing tasks, supply the exact direct task-contract move command
using active external configuration and one short inline Codex launch prompt;
keep personal paths out of tracked memory. Authorized Work memory-only updates
need neither a Codex launch nor a farm operation.

### Work owns authorized memory updates too

Work can and should apply authorized repository-memory updates directly,
including the governing scientific objective, CURRENT, durable knowledge,
workflow preferences and accountability history. Do not hand memory-only work
back to the user as a Codex task or leave it proposed solely because it edits
tracked Markdown. No separate Codex contract/launch or farm run is required.
Read the required context, preserve pending work, make the bounded update,
review its actual diff, regenerate/check the manifest, run memory health and
deliver the requested ZIP/receipt. Source-changing implementation retains its
established Codex workflow. The user alone publishes or runs the farm.
The user explicitly clarified this on November 4 after the assistant wrongly
deferred the scientific-objective correction and supplied Codex instructions.

### Scientific objective controls recommendations

The user's established objective is clean background-subtracted K+Lambda
data yields and corresponding combined data/SIMC ratios flat versus phi for
every t at both epsilons, with later SIMC-model iterations approaching unity.
Iteration-zero pion subtraction stays fixed. Setting-local E.8 debugging and
historical parent sums are not the final combined experimental observable.
Before substantial Work, identify exactly how its new result advances this
objective under the approved six steps, whether the needed producer/cohort
inputs actually exist and whether the calculation can answer that physics
question. If the first input is known missing, identify it and stop that
comparison; do not authorize another report/study to restate the boundary.
Method A must use matched selected Local-t branches and the same SIMC
denominator; assess ratio shape, normalization, background validity and
uncertainties separately. Ratios cannot justify forcing background or SIMC
normalization. Existing safeguards must be applied before spending resources;
new documentation is not scientific throughput. See the existing November 4
accountability continuation. No numerical Work follows from this correction.

### Complete scientific requirements and distinct investigations

Keep these requirements together in scientific recommendations: physically
supported pion subtraction; close cleaned-data/K+Lambda-SIMC missing-mass
shape agreement in the actual Lambda window (peak and relevant tail); then
corresponding combined data/SIMC ratios flat versus phi at each t/both epsilons
and eventual unity through existing SIMC-model iteration with experimental
subtraction fixed at iteration zero, followed by existing L/T extraction.
The experimental yield itself need not be flat. Integrated-yield agreement is
insufficient. Match selections, bins, resolution/acceptance and documented
normalization; keep supported shape-normalized diagnostics separate from
absolute yields. Shape/ratio agreement cannot alone prove purity or justify
tuning background or scaling SIMC. Full-MM plots do not replace signal-window QA.

Global / Local-t / Local-(t,phi) are the three model scopes; Local-t is selected,
not automatically production-adopted. Source-aware statistical fitting is a
separate promising investigation in the user's assessment, not Local-t's
definition/automatic extension or a fourth scope. Consult the detailed November
2-3 record before interpreting it: complete staged component refits changed
applied estimates and studied conditional coverage, variance, bias/MSE and
covariance; this did not establish universal superiority, truth or adoption.
Do not repeat source-aware fits or reopen model selection through memory work.

The primary pair is selected Local-t versus the SAME Local-t + Method A,
matched in source/cohort/fit/cuts/window/normalization and using the same initial
SIMC denominator. Assess both Lambda-window MM shape and ratio-versus-phi,
plus uncertainty/background validity. Historical Left/lowe yield shifts, parent
stability or attractive errors cannot substitute. Retain six presentations:
the primary Local-t pair and four historical Global/Local-(t,phi) before/after-A
comparisons. Method A remains detached/non-production; B diagnostic/excluded.

### November 3 explicit user correction: judgment before commissioning Work

ChatGPT repeatedly reopened the selected Local-t method, commissioned redundant
Outcome B and SIMC-normalization studies, and displaced the approved roadmap.
These were avoidable ChatGPT planning/review errors, not user failures or a
consequence of Work lacking farm access. Existing instructions were read but
not applied. See the [accountability record](investigations/kaonlt-assistant-roadmap-displacement-2026-11-03.md).

The user-approved roadmap governs priorities. Preserve completed decisions,
accepted scoped runtime evidence and mature analysis functionality. Before any
substantial Work assignment, inspect relevant source and accepted runtime records:
is it implemented, has its applicable gate passed, has the investigation already
concluded, what specific new observable/result/invariant can this task evaluate,
and does it directly advance the approved roadmap? Do not launch duplicate work
or work with no demonstrated path to new evidence. A report/receipt is not itself
scientific progress. Apply the existing displacement alarm before spending resources.

A blocked gate or technically investigable independent topic does not authorize
a priority change. A proposed alteration requires concrete new evidence and user
authorization before becoming an execution task. Missing files in an offline
review package do not establish missing production functionality or justify
normalization archaeology, redesign, refactoring or repeated validation. Preserve
the exact absolute-unit evidence limitation without inventing a correction or
expanding it into an analysis-wide prerequisite. No new SIMC normalization
correction was demonstrated necessary. Preserve accepted scopes without either
expanding or downgrading them. Use accumulated source/runtime knowledge for
professional approval judgment; do not restart validation in each new chat.

### November 4 correction: inspect both epsilons before declaring missing data

For the Lambda-window investigation, use the already supplied low AND high
farm evidence and all five requested settings. Starting with low epsilon is
an ordering choice, not permission to omit high epsilon. Before asking for
data or issuing a task-wide blocker, inventory original packets/receipts and
retained data/SIMC arrays; distinguish an archive-local omission from exact
selected-Local-t lineage and from globally unavailable evidence. Assess each
setting independently and continue supported comparisons. Do not confuse
135 canonical diagnostic slots with 135 populated/accepted measurements.
Do not use hard-coded absence fields or integrity check counts as proof of
scientific unavailability. Preserve historical-versus-selected branch identity,
available errors and scope; never fabricate comparisons to satisfy the request.
Do not substitute another blocker ZIP, memory loop or promise for actual work.
These are execution requirements following repeated assistant failures recorded
in the existing accountability record, not a new scientific phase. No absolute
promise of error-free future behavior is warranted; demonstrate compliance by
the evidence actually inspected and the comparisons actually completed.

### Reusable testing tools instead of per-run rewrites

The user requests a modular testing environment in the existing test area:
accept new artifact paths, cohorts, settings and supported analysis options
as validated runtime inputs instead of rewriting scripts for each run. Begin
with the usual memory checkpoint and actual audit of generated Work-ZIP code;
use that audit and existing repository tools to scope infrastructure afterward.
Keep adapters, pure numerical checks and reporting separate; retain provenance,
signed/invalid/unavailable states, and diagnostic versus producer authority.
Do not bake workstation paths or change production physics to accommodate tests.

### November 5 meeting decisions and immediate priorities

The user reports colleagues' general Lambda peak/relevant-tail shape agreement
for latest Local-t data versus unchanged iteration-zero SIMC at parent and
populated child level. Retain Local-t as selected working baseline; this is
qualitative, not final yield/covariance/purity acceptance. Keep BOTH unit-normalized
shape and existing yield-normalized magnitude diagnostics. Unit normalization
belongs only to the shape view and never rescales authoritative yield spectra.
The old iteration-zero cross-section parameterization need not give D/M unity.

Test Center/low and Center/high t3 restricted high-minus-t coverage first:
acceptance-edge/sparse/nonuniform kinematics is a plausible explanation, not
established cause or evidence of pion-model failure. Review effects beneath
Lambda even when unusual pion shape extends beyond it. Do not reopen model
selection or repeat normalization/source-aware investigations without new evidence.

Approved order: all-five full-t/acceptance-edge diagnostic with canonical
boundaries; colleague's actual exclusion procedure when received; remaining
beneath-Lambda background discrepancies; Method A only if a physical target is
justified; existing final-yield gates. Empty/inadequately populated cells should
be excluded by the colleague's pending procedure, but no numerical criterion or
final mask is agreed. Review support/errors/acceptance/angular coverage separately
for individual and combined settings; no invented threshold or source-data deletion.
No production adoption before relevant scientific/source gates. Do not force A.

Preserve SAME Local-t versus Local-t+A comparison, detached A, signed parent
normalization and diagnostic/excluded B. A slow-proton probability analogy does
not validate a pion probability model. Outcome B and canonical-five equivalence
remain unresolved. Final combined ratios still need producer-owned records,
errors/validity/lineage and covariance evidence; available offline acceptance
diagnostics do not require that later gate to be complete. After acceptance,
sum complementary setting yields, inspect six phi ratios, continue existing SIMC
iterations with experimental subtraction fixed at iteration zero, then existing L/T.
Preserve all six roadmap statuses/dependencies. This task updates memory only
and stops after verification; Check 1 is a separate task, with no farm authority.

Earlier November 5 comparison requirements below retain their yield-calculation
scope; the meeting explicitly preserves a separate unit-normalized shape view.

The November 5 clarification keeps the immediate MM investigation at full
setting scope across all five settings and all canonical t/phi cells. Do not
replace it with one Center/low t3 cell or describe an input inventory as the
completed scientific comparison. The single next evidence gate must retain
complete-setting scope; priorities within it do not redefine the deliverable.

The user's November 5 repeated-failure correction explicitly requires using
the implemented get_eff_charge normalization and existing histogram/yield chain.
Do not commission further normalization investigations, add scale factors, or
turn historical unit-provenance caveats into a new analysis prerequisite.
Inspect actual saved numerical content before claiming it unavailable; JSON
omissions or NPZ provenance questions are not missing normalization functionality.
Record responsibility under ChatGPT without blaming user, Codex, Work or JLab,
then resume all-five-setting Lambda-window physics rather than another general
audit or documentation-only substitute when applicable spectra are available.

For the November 5 normalized Lambda-window calculation, preserve existing
experimental and SIMC yield normalization in both coarse/fine views; do not
unit-area normalize, multiply normfac_simc twice, fit a scale or peak match.
Use the user-supplied All Five Settings Observer Shape Atlas as the visual
reference only: five-setting overview, consistent t columns, black data points
with errors, blue SIMC steps, light grids, clear units/legends, and aligned
residual/component panels. Gray boundary shading denotes actual unavailable
support only. Update labels to yield units; retain signed negative observations
and explicit empty/invalid states. The reference's observer/template content
and unit-area normalization do not authorize numerical substitutions.

November 5 correction: inspect retained ordinary standard analysis results,
including compatible other-run producer spectra where scientifically appropriate;
temporary diagnostic cleanup does not establish their absence. Use genuine
paired data/SIMC outputs immediately, preserving yields and recording the exact
cohort/subtraction scope. Do not silently relabel an older baseline as later
selected-Local-t validation. A fixed per-setting SIMC normalization constant is
not a statement of equal histogram contents. Do not repeat normalization work,
generalized inventories, producer rebuilding, or requests for the same run files.
Only a genuine accessible-output limitation may stop the affected calculation;
the user alone authorizes any subsequent bounded farm acquisition.

November 5 final correction: the experimental cohort for this task is October
28-30 Local-t, unchanged. The October 9 comparison is explicitly rejected.
Compatible other-run SIMC alone may be used when its actual selection, iteration,
weights and bin contents match; this never licenses older experimental data.
Do not confuse unit-integral SIMC fit/reference shapes with absolute yield
histograms or restore their amplitude with normfac_simc or fitted signal scales.
Do not split nonaligned bins by assuming an unretained within-bin distribution.
Recover matching yield information directly; if it is genuinely not retained in
accessible evidence, state the exact limitation and stop. No broad investigation,
new producer, automatic Method-A adoption, source modification or farm operation.

The corrective latest-Local-t/existing-SIMC handoff explicitly permits the
retained October 9 KLambda SIMC yield arrays because that simulation is unchanged.
Keep October 28-30 experimental data; never replace them with October 9 data.
Different run dates alone do not invalidate unchanged simulation. The factor is
already applied to retained SIMC yields: no second factor or unit-area scaling.
Use native bin edges and supported bounds where exact rebinning is unavailable;
do not turn an axis mismatch into another normalization or acquisition study.

## Delivery preferences

For a source-changing task, ChatGPT provides exactly one standalone Markdown
**task contract**, placed in the repository, followed by one short **inline
Codex launch prompt** telling Codex to read repository memory and that contract.
The contract is the durable scoped implementation authority; the inline prompt
is only the launch instruction. Do not create a second standalone Codex-prompt
Markdown file unless the user explicitly asks. The contract file, not a shell
heredoc, is the normal transfer artifact. After Codex completes, ChatGPT reviews
the actual diff rather than accepting the Codex summary.

When ChatGPT generates the single task contract as a downloadable file, provide
the exact direct move command from the configured Downloads location into the
exact repository contract path, using active external environment values. Do
not ask the user to place it manually or insert unnecessary file-existence
checks when a direct scoped move suffices. Concrete personal paths remain
external and must not be stored in public repository memory.

For farm work, provide concise, exact instructions appropriate to the established
JLab environment, with reproducible commands and supplied paths. Preserve
cumulative context in durable repository records rather than giant
chat-continuation prompts.

Personal workstation, Downloads and machine-specific path values are external
environment/session configuration, not durable public repository memory. Use
the supplied active environment/ChatGPT Project configuration; never substitute
a guessed generic path. Ask only when a needed value is genuinely unavailable.

When the user identifies workflow, source-state, pathing, procedure, evidence or
task misalignment, stop extending that approach and follow root AGENTS'
professional re-anchor: re-establish identity, reread the startup core and
relevant CURRENT records, identify the exact error and resume from the earliest
valid gate. Preserve accepted work and closed phases; do not rationalize the
previous approach or invent a replacement without a concrete blocker.

Treat farm requests as safety-critical for an ordinary JLab user account, not
an administrator account. Every farm command must be valid `tcsh`, scoped to
the user's designated paths and the requested gate, and minimize changed
state. Do not prescribe privileged, scheduler, service, shared-filesystem, or
ordinary-checkout cleanup actions; any temporary-worktree cleanup must be
explicitly path-bounded and limited to the worktree created by the request.

For every completed local change set, include an exact, non-executed,
user-controlled Git handoff block with the scoped `git add`, `git commit`, and
`git push origin test` commands. State the exact changed paths and make clear
that the user reviews and runs those commands; Codex never commits or pushes.

## Boundaries

This record is stable collaboration context only; it is not a command
reference, execution-authority policy, farm procedure, active state,
scientific evidence, phase chronology, or source-identity ledger. See
[root AGENTS.md](../../AGENTS.md) for behavior and scientific boundaries,
[CURRENT.md](CURRENT.md) for active state, [TOOLS.md](TOOLS.md) for canonical
operations, [COMMUNICATION.md](COMMUNICATION.md) for farm delivery, and
[CODEX.md](CODEX.md) for source-changing workflow.
