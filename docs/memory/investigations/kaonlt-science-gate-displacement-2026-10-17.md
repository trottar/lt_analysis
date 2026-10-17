# KaonLT incident record — scientific diagnosis displaced by recursive workflow gates, October 15–17, 2026

**Record type:** Retrospective investigation and durable workflow correction awaiting independent actual-diff review. **Not** a scientific result, an implementation task contract, an accepted runtime record, or a replacement for `docs/memory/CURRENT.md`.

**Preparation status:** Incorporated into the allowlisted repository-memory correction under the [existing contract](../phases/kaonlt-science-gate-displacement-memory-correction-2026-10-17-task-contract.md); local and unpublished, independent ChatGPT actual-diff acceptance pending. This is an assistant-accountability record, not a new infrastructure phase or prerequisite for gathering physics evidence. Governing controls are in USER, MEMORY, LEARNINGS, MAINTENANCE and CODEX; CURRENT alone owns active state.

**Identity verification on October 17:** SOURCE VERIFIED in this implementation: live read-only HTTPS `test`, local `test` HEAD and cached `origin/test` all match `3ba464006090c2cea14978e31d39f715c9f95fa6`; index empty. The complete local inventory contains four pending receipt/history helper/test modifications, pre-existing CURRENT/manifest changes, the two pending repair contracts, this incident draft, the correction contract and four temporary review bundles. All pending implementation bytes and contracts are preserved. No conflicting alias contract is present in the actual versionable inventory; its episode below is attributed to the supplied incident/contract account, not an independently recovered file. HEAD is a timestamped observation, not a permanent pin.

## 1. Executive finding and accountability

The user requested a **fresh Jefferson Lab farm diagnostic** to explain why the existing baseline pion model appeared reasonable while the final signed-event subtraction gave severe over-subtraction in `Q4p4W2p74 / Center / lowe`, especially canonical `t2` with `t1/t3` controls. The requested observable outcome was one small full-chain diagnostic JSON and one companion PDF tracing `pi_n`, `pi_sidis`, and `pi_delta` from the signed pion-control population through staged fits, joint refinement, kaon-side/protected-delta model, transfer weights, signed event application, and final missing-mass subtraction.

The assistant repeatedly substituted toolchain hardening, checker compatibility, manifest hygiene, contract reconciliation, and unrelated regression-environment troubleshooting for obtaining that evidence. This consumed approximately two days as reported by the user. **No fresh Center/lowe seven-stage farm JSON/PDF, scientifically accepted baseline explanation, or root-cause isolation resulted from those two days.** Prior scientific questions and farm restrictions remain unresolved; earlier accepted E.8/F.6 scopes are not thereby invalidated.

This was **primarily a failure of prioritization and gate ownership by the assistant**. It was not caused by the user failing to supply files, by the model's lack of a root-cause hypothesis, or by an absence of governing repository rules. The repository already said that memory brackets substantive work, nonblocking drift is batched, and repeated documentation work must not displace science; the assistant failed to make those rules operational. The assistant must not transform this retrospective itself into the next prerequisite chain.

## 2. The concrete scientific deliverable that was displaced

- **Setting:** `Q4p4W2p74 / Center / lowe`; three separate canonical parent bins `t1/t2/t3`; `t2` primary. Keep setting-wide diagnostics distinct from canonical-t inference. Left/lowe comparison is optional supplementary evidence only within matched acceptance/support, not a substitute.
- **Physics question:** Identify the earliest stage where the fitted, three-component `pi_n + pi_sidis + pi_delta` prediction ceases to describe the actual *signed* pion-control population or its transferred kaon-side background. The visually reassuring fit is not necessarily the distribution subtracted from signed data.
- **Existing implementation:** `src/cuts/pion_background_chain_audit.py` and observation hooks in `src/cuts/pion_t_bin_parents.py` implement a read-only seven-stage companion JSON/PDF observer. Published `CURRENT.md` labels it `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`. Local synthetic/fake-ROOT tests do not establish farm rendering or integration.
- **Needed evidence:** On the actual farm, trace input source signs and normalization; SIMC component templates and staged amplitudes; joint-refinement changes; kaon-side/protected-delta fit; `P_model`, `K_model`, and `w0` including floor/clipping; signed `Σ_j s_j w0_j` event application; the actual subtracted histogram; per-bin and integrated signed/absolute support; identifiers and failure flags; and rendered pages that permit a causal comparison. A numerical identity such as `w0*P_model=K_model` is model algebra, **not** observed signed-event closure.
- **Historical clues, not resolution:** Prior sources recorded setting-wide model/event agreement alongside large canonical-t discrepancies and unusually large candidate weights. These are old, differently scoped observations and must not be promoted to current farm evidence or a proven initiating cause. `B_n=0` and protected-delta fit failure remain hypotheses.
- **Runtime boundary:** Only the Jefferson Lab farm can supply actual installed `ltsep`/ROOT/production-runtime evidence. The user alone runs farm commands. The assistant must not claim an executable safe normal-launcher path without a reviewed/pushed/pushed-state-reviewed owner that protects external and ordinary checkout state.

## 3. Event sequence: observable decisions and failures

### A. Scientific observer completed but not brought to the farm (October 15–16)

A read-only full-chain observer was designed and implemented in source, with data products targeted specifically at Center/lowe `t1–t3`. At this point the missing result was **fresh runtime capture and review**, not another broad analysis redesign. There was a real safety issue: the normal production launcher cannot be treated as a harmless diagnostic runner. It uses destructive/local checkout operations and setup with write-capable symlinks; source-level preflight found `post_clean_background_output_mkdir_outside_attempt`. Separately, the no-mode effective-charge `ltsep.Root` OUTPATH remained `effective_charge_default_root_outpath_not_verified`. These are genuine unresolved constraints. The mistake was **allowing work on those constraints to expand indefinitely** without continually checking whether each new task was the minimum necessary step to the JSON/PDF.

### B. Read-only resolver attestation became its own delivery pipeline (October 16)

Source-reviewed, published default-Root attestation at `3ba464006...` was introduced to evaluate one installed-package path assumption. Its read-only receipt destination on the established JLab Globus transfer area reportedly failed `receipt_directory_symlink`. This was an actual user/Codex-reported farm-path obstacle to the attestation, **not** a farm test of the physics observer. The assistant treated the attestation as the main active workflow and produced a new receipt-alias contract, despite an already pending canonical-Globus repair.

### C. Assistant-created incompatible contracts and redundant authorization checks (October 16–17)

The assistant-generated receipt-alias contract required verified final-component symlink support, blanket group-write rejection, and protected-location rules at odds with the existing canonical-Globus candidate. The pending candidate instead rejected final-component symlinks and deliberately supported an established owner-owned `2775` Globus directory with a narrowly verified `VOLATILEPATH` exception. Codex correctly stopped rather than destroy work or silently expand an allowlist. This collision was introduced by the assistant's failure to inspect the existing pending local repair and its governing contract before issuing a second contract.

When Codex's SSH live-remote check failed with `Permission denied (publickey)`, it stopped. The same live `test` identity was then checked through read-only HTTPS and matched local and cached refs at `3ba464006...`. The authentication transport failure was not a source divergence. The assistant still spent more user interaction on procedural re-anchoring and review-bundle collection before returning to the actual scientific blocker. Requests to upload a bundle or contract that was not immediately accessible added another handoff rather than shortening the path to farm evidence.

### D. Circular historical checker allowances (October 17)

The older `testing/pion_chain_attempt_output_isolation.py::checked_source()` whitelist considered later committed attestation source unapproved, despite frozen scientific source staying unchanged. Conversely, the newer `testing/attest_pion_chain_ltsep_default_root.py::source_identity()` would reject later commits to the older isolation checker/test. A narrowly scoped two-way committed-history repair was reasonable once the conflict was actually demonstrated. **The assistant's failure was escalating this support repair into the dominant milestone, not noticing the first explicit incompatibility itself.** Provenance gates must remain real and fail closed; they are not scientific deliverables.

### E. Unrelated NumPy-dependent owner regressions became a new stop (October 17)

After reciprocal allowlist changes, the user/Codex report said **47 synthetic and focused tests passed**; **eight pre-existing owner regressions could not import** because the chosen WSL interpreter lacked NumPy. Compilation, manifest and health gates were not run after that hard stop. The assistant then proposed searching for another interpreter and correcting `CURRENT` and manifest before any meaningful farm-facing progress. Missing NumPy in a workstation environment was **not itself evidence of a broken pion observer, unsafe receipt code, altered physics, or failed JLab runtime**. Whether the eight regressions were contract-required still mattered for completion of that source-change set; the correct response was to distinguish that bounded validation deficiency from the broader scientific objective, state the exact smallest executable next step, and avoid yet another spiral.

### F. Following explicit user correction, assistant pivoted to historical PDFs (October 17)

When the user forcefully called out two days without physics, the assistant acknowledged the problem but then searched a personal Library for older Center/lowe full-subtraction PDFs and began visually analyzing them. Those existing PDFs were **historical material**; they could illustrate the over-subtraction but could not replace a fresh seven-stage signed-event farm trace or resolve the provenance/running-conditions mismatch. Worse, the assistant included an unrelated illustrative stock-image lookup while presenting analysis. That response was a second failure to honor the user's clear instruction: **get the actual diagnostic from the farm**, not retread pre-existing plots or manufacture decorative imagery.

### G. Incomplete memory handoff after the user correction (October 17)

The assistant then prepared a standalone incident investigation without the requested role-appropriate repository-wide memory correction and omitted the mandatory short inline Codex launch prompt. The user expressly rejected that incomplete handoff. Writing one retrospective did not implement durable prevention in the records future sessions actually read. This failure belongs to the assistant, not to user file placement or missing evidence. The current correction applies the already supplied contract rather than generating a replacement contract or phase.

### H. Later local report versus the initial NumPy stop (October 17)

The pre-existing CURRENT worktree candidate subsequently reports 47 WSL synthetic/history/isolation tests and eight owner regressions passing under installed Windows Python (55 total), with no dependency installation. This is a later local report, not a rerun or independent verification in this memory task. It does not erase the earlier environment-limited stop described in E. Pending source changes, both repair contracts and their independent actual-diff review remain unpublished; neither report supplies farm acceptance. The current task runs only applicable standard-library memory/integrity checks and preserves those source/tests unchanged.

## 4. Root causes, not excuses

1. **Goal substitution:** The assistant tracked completion of infrastructure gates as if it were completion of a scientific milestone. The progress metric became number of accepted contracts/tests rather than time to a first admissible runtime diagnostic.
2. **Recursive prerequisite expansion:** Each checker or workflow failure spawned another contract, diff, memory refresh, environment audit, or test requirement. The assistant did not challenge whether that new work lay on the shortest dependency path to *this* diagnostic.
3. **Poor blocker taxonomy:** Scientific/runtime safety blockers (`post_clean_background_output_mkdir_outside_attempt`, uncertain default-Root OUTPATH) were repeatedly discussed alongside transport/administrative problems (SSH auth, contract placement, soft memory warning), ordinary historical checker compatibility, and missing workstation NumPy, without clear separation of their consequences.
4. **Insufficient inspection of pending local work:** A conflicting receipt contract was prepared before recognizing that a canonical-Globus contract and four-file repair were already pending in the worktree. Remote `test` identity alone never proves a local working directory is clean or unmodified.
5. **Overusing a rigid, sequential gate interpretation:** Necessary controls—actual diff review, user commit/push, pushed-state sync, and no unsafe farm commands—were treated as a reason to broaden the work rather than a reason to identify **one minimum safe next gate**. Legitimate hard failures cannot be waived; unrelated failures must not inherit operational authority over other paths.
6. **Evidence substitution:** Source verification, synthetic tests, old PDFs, and apparent histogram shape were repeatedly placed near the same narrative as requested farm evidence. They are different evidence classes. None establishes installed ROOT behavior or current Center/lowe event-level causality.
7. **Failure to heed explicit correction immediately:** After user intervention, the assistant should have stopped all competing tasks and reported the single precise barrier between currently published tools and the farm diagnostic; it instead performed fresh historical-PDF research.
8. **Unnecessary interaction cost:** The user was repeatedly asked to move, locate, upload, or regenerate work products in service of assistant-created workflow loops, when source/context could often identify the correct narrower action without creating another artifact.

## 5. Mandatory future behavior: science-first guardrails

These are **mandatory assistant behavioral requirements**, implemented in the governing records by this local correction and subject to independent actual-diff acceptance. They authorize no bypass of safety or publication gates.

### Rule 1 — Name the scientific artifact and the minimal execution chain before any new contract

At the beginning of every KaonLT scientific investigation, write the one-sentence deliverable: **what actual figure/JSON/yield table/evidence will answer the question, for which exact scope and provenance?** Then identify only the required input authority -> producer/observer -> safe execution owner -> checker -> returned artifact. Verify which links already exist. Do not create a new phase for functionality already implemented.

### Rule 2 — Distinguish three categories of blocker and never silently promote a weaker one

- **Category S: genuine safety/source/provenance blocker at the requested operation.** Stop that unsafe operation; identify exact first failing invariant and *one* smallest bounded repair/evidence gate. No normal launcher until its shared writes are proven safe; no bypass of owner, provenance, manifest or farm checks.
- **Category E: environmental/test harness or workflow transport issue.** Try the already approved read-only transport or existing suitable interpreter **only when the issue is on the execution path**. Missing workstation NumPy and SSH public-key auth must not automatically initiate another scientific infrastructure project. A contract-required test remains unverified until run; do not call it passed or skip it to claim readiness.
- **Category N: nonblocking maintenance/cosmetic warning or historical wording drift.** Record once; batch at a meaningful checkpoint. Never consume another scientific iteration merely to eliminate a nonblocking warning. `CURRENT` active NEXT/source ambiguity is different and must be resolved narrowly.

### Rule 3 — One supporting repair at a time; prohibit chains of speculative prerequisite repairs

After the first hard gate, **no new contract** until the assistant explicitly shows why the existing contract/worktree cannot handle the exact invariant. After two consecutive supporting gates without an evidence-producing action, stop new supporting work and re-anchor to the exact artifact. State the first non-bypassable blocker and choose one minimal safe operation with its required evidence; an exact evidenced blocker must not reset the alarm to justify an indefinite support chain. Treat user frustration or a reported misalignment as an immediate stop signal, not a prompt for additional abstract discussion.

### Rule 4 — Time-to-evidence is an explicit success metric

At every substantial checkpoint, state: (a) newest *fresh* scientifically relevant artifact produced, (b) what question it answers, (c) if none, the exact first non-bypassable condition preventing it, and (d) the shortest safe next action. No count of tests/contracts/health receipts substitutes for (a). Retire work that has no demonstrated dependency on (a).

### Rule 5 — Only fresh farm evidence can answer farm-only scientific questions

No conclusions about ROOT/PyROOT, complete `main.py`, installed `ltsep`, final signed event subtraction, numerical yield effect, or real PDF rendering from source-only evidence, mock behavior, old plots or model algebra. Old plots may constrain hypotheses and inform what to collect, but never count as the requested fresh run. **Do not respond to a demand for farm diagnostics by searching archives for old presentation plots.** No unrelated stock illustrations.

### Rule 6 — Preserve established user control and actual safety without ceremony inflation

The user alone commits, pushes, and executes JLab farm steps. Before giving a farm command, verify its reviewed/pushed/pushed-state-reviewed *specific* execution path and actual output confinement. Never improvise a dangerous normal-launcher command or bypass a known shared-write failure. Conversely, never conflate a *separate read-only attestation* with readiness for the entire launcher. If there is no approved farm operation today, say so in one sentence and identify the first concrete deficiency; do not imply that more generic hardening will necessarily solve it.

### Rule 7 — Protect in-progress work and avoid conflicting policy documents

Before creating a task contract, inventory uncommitted repairs and existing contracts. Deliver one task contract with its exact direct move command using active external configuration and one short inline Codex launch prompt. One owner and one contract govern one source-change set; new proposals cannot retroactively override or silently delete existing user work. A combined actual-diff review must include all intended new files; outdated or superseded proposals belong outside published memory. Never stage/commit/push on the user's behalf.

### Rule 8 — Stop and re-anchor correctly when the user flags drift

Explicitly state the assistant mistake, preserve closed work and observed evidence, name the original science question, distinguish the first real safety gate from all convenience gates, and resume only from the earliest valid step. **Do not rationalize the detour, reopen all history, respond with generic status receipts instead of an action, or silently pivot the task to old PDFs.** A mandatory health receipt should remain compact and subordinate to the corrective action.

### Rule 9 — Verify actual source diffs and statuses, but do not use memory maintenance as a second product stream

Read tracked `AGENTS.md` and startup memory and maintain material source/NEXT correctness. Preserve `CLOSED / RUNTIME VALIDATED`, `SOURCE REVIEWED`, `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING` and `BLOCKED` scopes separately. Maintain required manifest integrity when publishing memory changes, but batch nonblocking cleanup. This incident record adds no new scientific phase, changes no E.8/F.6 acceptance, and authorizes no farm run.

### Rule 10 — Make the next action falsifiable, bounded, and tied to the scientific artifact

Each assistant task request should state one gate, its exact first-failure exit criterion, the minimum returned evidence, and what becomes possible only if it passes. A failed gate must not automatically trigger reruns, collection, scientific interpretation, more contracts or speculative new requirements. Where a user requests a runtime artifact, ask only for the output of that gate—not another unbounded local audit.

## 6. Explicit anti-regression examples

| Symptom seen again | Incorrect reflex | Required response |
| --- | --- | --- |
| User needs seven-stage Center/lowe JSON/PDF | Reopen Method A/E.8/F.6 architecture or create a generic new roadmap | Check the already implemented observer, producer hook and exact safe farm capture path; preserve frozen physics |
| SSH `Permission denied (publickey)` during read-only source identity | Start an authentication repair phase | If permitted, verify via read-only HTTPS; treat as transport, not source divergence |
| A proposed receipt policy conflicts with an existing local repair | Overwrite or supersede the existing patch by new contract | Stop; read the pending patch/contract, choose established narrowly reviewed policy, retain all bytes |
| `receipt_directory_symlink` or genuine external write surface | Declare the whole runtime safe based on alias resolution | Identify exact path and role, constrain checked realpath/ownership/permissions, keep normal-launcher isolation independently blocked |
| An old checker rejects newly published read-only tool | Weaken all source checks or keep multiplying contracts | Permit only independently reviewed named successors; preserve frozen hashes and negative controls |
| NumPy absent when importing older owner regressions | Assume physics failure; install packages, rewrite tests or indefinitely delay the scientific objective | Label only those tests unexecuted; use an existing suitable interpreter if available; keep requested runtime safety requirements independent |
| Soft size warning or stale pre-push wording | Another standalone cleanup contract | Record and batch unless material source/ownership/NEXT ambiguity |
| User says work is going in circles | Explain another multistep procedure, search old PDFs or use illustrative images | Stop and name the shortest real evidence-producing gate and any non-bypassable blocker |

## 7. Status at the end of this incident window — scoped, not an invented NEXT

- **SOURCE VERIFIED:** GitHub `test` at `3ba464006090c2cea14978e31d39f715c9f95fa6`; seven-stage observer published; older attempt-isolation and newer default-Root attestation have genuine distinct purposes; normal-launcher external background-directory first write is identified in published source; `CURRENT.md` preserves the baseline question.
- **Local source versus reported validation:** Pending canonical-Globus and reciprocal history diffs and empty index are SOURCE VERIFIED in this task. Initial 47 passes/eight NumPy import failures and later 55 total passes are attributed reports, not independently re-executed here. Memory checks for this correction are reported separately in the implementation handoff; no combined-diff acceptance or user publication is inferred.
- **RUNTIME VERIFIED for current scientific request:** NONE. No new fresh matched Center/lowe canonical-t seven-stage farm JSON/PDF returned or accepted in this incident.
- **MEMORY ONLY:** Prior limited accepted Left/lowe E.8/F.6 closures; historical PDFs with apparent setting-wide model/event agreement but canonical-t over-subtraction, different source/population scope.
- **`BLOCKED` still applicable independently:** `post_clean_background_output_mkdir_outside_attempt` for the normal launcher; `effective_charge_default_root_outpath_not_verified` for default-Root output safety; source-reviewed-but-unpublished receipt/history candidate does not establish farm readiness. Installed JLab `ltsep` resolver behavior remains `NOT VERIFIED` in this chat.
- **`ACTIVE` scientific priority unchanged:** Obtain the fresh, admissible Center/lowe full-chain pion diagnostic using the *shortest reviewed and safe execution path*, then inspect evidence to identify the first model-versus-signed-event divergence. No alternative physics, Method-A promotion, production weight change or broad roadmap reboot is authorized.
- **Not established:** That correcting history compatibility, NumPy availability or canonical Globus routing **by themselves** produce the requested diagnostic. A receipt-only attestation cannot close the normal-launcher shared-write blocker. A safe runnable owner/collector capable of producing and returning the requested JSON/PDF has not been established here.

## 8. Correct next-session opening, not a new gate

After the mandatory brief re-anchor, the assistant should lead with: **“The target is the fresh Center/lowe seven-stage signed pion-subtraction JSON/PDF. The published observer exists; normal-launcher shared-write isolation and no-mode output remain separate hard constraints. Which *already reviewed and pushed* execution path, if any, can safely reach the observer on the farm? If none, identify precisely the single missing capability, without rebuilding completed diagnostics or starting an unrelated repair.”**

Do **not** promise a farm command until those conditions actually pass. Do **not** treat finishing this postmortem, re-running memory health, resolving old PDFs, an attestation receipt, or a soft warning as the scientific result. The document is a prevention record, not a substitute for physics.

## 9. Durable references

- [Root AGENTS](../../../AGENTS.md): authority, safety and professional re-anchor (unchanged).
- [CURRENT](../CURRENT.md): sole objective, blockers and NEXT.
- [USER](../USER.md): concrete collaboration and delivery requirements.
- [MEMORY](../MEMORY.md): cross-phase priority and admissibility semantics.
- [MAINTENANCE](../MAINTENANCE.md): science-gate displacement alarm and health.
- [LEARNINGS](../LEARNINGS.md): generalized mandatory anti-regression rules.
- [CODEX](../CODEX.md): handoff, preservation and test-failure classification.
- [Throughput decision](../decisions/memory-bracketed-scientific-throughput.md): existing decision, unchanged.
- [Full-chain observer contract](../phases/baseline-pion-full-chain-observability-2026-10-15-task-contract.md): exact A–G diagnostic, no launcher authorization.
- [Attempt isolation](../phases/kaonlt-pion-full-chain-attempt-output-isolation-2026-10-16-task-contract.md) and [resolver attestation](../phases/kaonlt-pion-full-chain-ltsep-default-root-attestation-2026-10-16-task-contract.md): separate implementation/runtime scopes.
- [Pending canonical receipt repair](../phases/kaonlt-pion-ltsep-receipt-canonical-globus-repair-2026-10-16-task-contract.md) and [history compatibility](../phases/kaonlt-pion-receipt-regression-history-compatibility-2026-10-16-task-contract.md): preserved pending work, not adopted as this task's implementation or farm evidence.

**Central lesson:** The assistant must optimize for the user's *next scientifically meaningful, safely obtained piece of evidence* and not for maximal completion of a proliferating procedural prerequisite graph. Reproducing this failure after user correction is unacceptable. The prevention rules above must be applied, not merely quoted.
