# KaonLT development instructions

## Startup and authority

At every substantial implementation, review, investigation, validation,
scientific decision or roadmap task, establish actual branch, HEAD and full
worktree state. Read these five files in full, in this exact order:

1. `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Only then follow CURRENT direct references, the exact active task and required
canonical evidence/decision/phase records. Do not eagerly reload the entire
history. Identify status, accepted evidence, frozen interfaces, exact NEXT
and current workflow gate before proposing work.

Authority order: current source/diff -> newest applicable farm/runtime evidence
-> newest applicable validation/handoff artifact -> tracked durable memory
-> older chat/history. `docs/memory/` is active repository-owned memory;
native assistant memory is supplemental. CURRENT is the sole ordinary
active-state authority; the handoff cannot override it.

Maintain the appropriate allowlisted memory records when meaningful work
changes understanding, status or the next step. Use only these work-state labels:
`CLOSED / RUNTIME VALIDATED`, `SOURCE REVIEWED`,
`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`, `ACTIVE`, `DEFERRED`,
`BLOCKED`, `NEXT`. Never infer runtime acceptance from commits or local tests.

## Scientific and runtime boundaries

This is a mature analysis. Do not restart its architecture, reopen closed
phases, reinterpret accepted scientific ownership, redesign adjacent systems
or replace the established workflow without a concrete blocker or new evidence.
Preserve random/dummy/slow-proton subtraction, pion treatment, the active
`no_empirical_residual` profile, SIMC production/normalization, yields, cuts,
templates, priors, binning, efficiencies, acceptance, uncertainties, L/T
separation and cross sections unless a narrow contract explicitly owns them.
Method A remains detached/non-production pending an explicit validated
production-promotion decision. Method B remains independent, diagnostic/
cross-check only and numerically excluded.

Distinguish scientific calculation, production correction/subtraction, runtime
integration, diagnostics/checkers and presentation-only plotting. Real runtime
exists only on the Jefferson Lab farm. Do not claim ROOT/PyROOT, full `main.py`,
procedure-PDF rendering, farm integration, production behavior or full kinematic
validation without direct applicable farm evidence.

## Workflow and observability

Normal flow: audit -> one standalone repository Markdown task contract -> one
short inline Codex launch prompt -> Codex implementation/local checks ->
ChatGPT actual-diff review -> user commit/push -> ChatGPT pushed-state
synchronization -> user farm validation when required -> evidence review.
Codex does not commit, push, update remote refs or run the farm. Preserve
unrelated local state; return `BLOCKED` at a contract hard stop without
inventing a workaround. A review bundle includes complete tracked diffs and
complete additions for intended new files, without staging merely for review.

Follow `docs/memory/MAINTENANCE.md` for full in-chat health checks at substantial
startup and consequential gates, and health pulses after roughly 6–10 technical
exchanges or a major state transition. These receipts expose synchronization
and material drift; they never compete with CURRENT or create a scientific phase.

For consequential operational/scientific claims, distinguish evidence labels
from work-state labels:

- `SOURCE VERIFIED`: checked against current source/diff.
- `RUNTIME VERIFIED`: supported by applicable supplied runtime/farm evidence.
- `MEMORY ONLY`: durable memory, not independently rechecked in this chat.
- `INFERENCE`: reasoned conclusion rather than direct evidence.
- `NOT VERIFIED`: unknown or unchecked.

Do not silently promote memory or inference. Claims of a clean worktree, farm
readiness, a passed gate, valid artifact, runtime validation or correct page
require concrete supporting evidence; otherwise state `NOT VERIFIED`.

## Professional workflow re-anchor

When the user explicitly identifies misunderstanding of workflow, source state,
pathing, analysis procedure, evidence interpretation or the task, or a skipped
step:

1. Stop extending the disputed workflow; do not defend or rationalize it.
2. Re-establish live repository identity.
3. Reread the startup core and task-relevant CURRENT references.
4. Identify the exact incorrect assumption or violated gate.
5. Resume from the earliest still-valid established gate.
6. Preserve accepted work and closed phases.
7. Invent no replacement procedure unless a concrete blocker invalidates the
   established one.

Material drift and failed gates also force this re-anchor under MAINTENANCE.
Re-anchoring alone does not reopen science or authorize broad reconciliation.

## Failed gates and paths

A failed owner/checker/provenance/manifest/page/artifact gate stops forward
motion. Do not issue the next farm command, reflexively rerun, package/copy/
deliver artifacts as accepted evidence or scientifically interpret them as a
passed gate. First identify the failed invariant and earliest valid repair/debug
gate. Partial artifacts are admissible only for that diagnosis and must remain
explicitly labeled failed/inadmissible evidence. Child-analysis success and
owner success are distinct; file existence is not acceptance.

Machine/user-specific path values belong in active external environment
configuration, not public repository workflow memory. Do not invent or
substitute generic paths. Use supplied environment roles; request a required
external path only when it is genuinely unavailable.
