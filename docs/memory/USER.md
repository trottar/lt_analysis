# KaonLT user collaboration context

## Project working context

KaonLT is a mature Jefferson Lab Hall C analysis. Do not restart its
architecture or reopen closed phases without concrete new evidence. Its
operational/runtime target is the Jefferson Lab Linux farm; ROOT/PyROOT, full
`main.py` behavior, and production procedure PDFs require farm evidence rather
than local-source assumptions. Discover the actual local interpreter and OS
for each task rather than treating workstation configuration as project truth.

## Collaboration preferences

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
