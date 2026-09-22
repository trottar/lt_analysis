# KaonLT user collaboration context

## Project working context

KaonLT is a mature Jefferson Lab Hall C analysis. Do not restart its
architecture or reopen closed phases without concrete new evidence. The
operational/runtime target is the Jefferson Lab Linux farm; ROOT/PyROOT, full
`main.py` behavior, and production procedure PDFs require farm evidence rather
than local-source assumptions. Discover the actual local interpreter and OS
for each task rather than encoding a workstation path as project truth.

## Collaboration preferences

Use one narrow phase or fix at a time. Inspect actual source and diffs rather
than trusting summaries. Preserve frozen scientific interfaces unless the task
contract explicitly owns a change. Do not request information already present
in repository memory, supplied artifacts, or current source. Prefer the
smallest complete validation-artifact package needed for the question, and keep
source review explicitly distinct from farm validation.

## Delivery preferences

Deliver Codex plans and prompts as standalone Markdown files. For farm work,
provide concise exact commands appropriate to the established JLab environment,
with reproducible commands and explicit paths. Preserve cumulative context in
durable repository records rather than relying on giant chat-continuation
prompts once the memory migration is complete.

## Boundaries

See [AGENTS.md](AGENTS.md) for execution and scientific operating rules,
[CURRENT.md](CURRENT.md) for active state, [TOOLS.md](TOOLS.md) for operational
references, and [COMMUNICATION.md](COMMUNICATION.md) for farm-delivery
conventions. This record is stable collaboration context only; it is not active
state, scientific evidence, phase chronology, or a source-identity ledger.
