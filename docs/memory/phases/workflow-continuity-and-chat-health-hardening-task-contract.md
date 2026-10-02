# KaonLT workflow-continuity and in-chat health hardening — task contract

## Status

`ACTIVE`

## Independent actual-diff review repair

Independent ChatGPT actual-diff review of the first hardening candidate found
two contract-level blockers that must be repaired before commit:

1. the new public root `AGENTS.md` is not yet reconciled with the existing
   repository-memory startup surfaces and bootstrap tooling, which still treat
   `docs/memory/AGENTS.md` as the first startup record; and
2. this contract itself contains literal examples of personal filesystem paths,
   contradicting the task's public-path externalization requirement.

The same review also found one push-stability issue in the hardening phase
record: it labels itself `SOURCE REVIEWED` while simultaneously describing
independent actual-diff review as still pending. The committed form must state
one unambiguous post-review status.

This revised contract owns only those narrow repairs in addition to the
already-implemented hardening candidate. Preserve all first-pass work unless
this contract explicitly requires a repair.

A second independent actual-diff review found two wording contradictions in
this contract itself: the opening froze all tests even though the dedicated
memory-health test is explicitly allowlisted, and the frozen-area wording
prohibited every checker even though the memory-health checker is explicitly
allowlisted. This revision narrows those prohibitions to scientific/runtime
tests and scientific/runtime checkers. No implementation scope is expanded.


This is a repository workflow/memory hardening task only.

The user has explicitly paused scientific analysis until the workflow hardening
discussed on 2026-10-02 is implemented in the public repository and the matching
ChatGPT Project configuration is updated separately.

This task must not alter KaonLT scientific calculations, production behavior,
scientific diagnostic/validation mathematics, plotting implementation, runtime
owners, farm profiles, farm wrappers, scientific/runtime tests, or physics
outputs. The memory-health/bootstrap tools and their dedicated deterministic
memory test are allowed only as explicitly listed below.

## Required starting source and repair worktree

Required branch:

`test`

Required committed HEAD and `origin/test`:

`ccaf15358efc205cd601aeccecdd1ab1b1360dff`

Expected commit subject:

`E8.4 Fix.5.6: harden owner failure provenance`

Codex must establish and report:

```text
git status --short --branch --untracked-files=all
git rev-parse HEAD
git branch --show-current
git rev-parse origin/test
git log --oneline -1
```

Required:

- branch is exactly `test`;
- local HEAD is exactly the required committed HEAD;
- `origin/test` is exactly the required committed HEAD.

This is a review-repair continuation, so the worktree is intentionally dirty
with the first-pass hardening candidate. Preserve that candidate and repair it
in place.

At repair start, changes are permitted only on:

- the already-proposed hardening paths from the first review bundle;
- the additional startup-consistency paths newly allowlisted below;
- the current hardening review bundle, which may be regenerated;
- the three explicitly preserved pre-existing untracked files described below.

Do not pull, reset, clean, stash, commit, push, move, delete, rename, or discard
pre-existing state.

The following pre-existing files remain frozen in place and excluded from this
task's commit/review content:

- `docs/memory/phases/workflow-continuity-hardening-task-contract.md`
- `kaonlt_review.diff`

The current contract remains the governing task file:

- `docs/memory/phases/workflow-continuity-and-chat-health-hardening-task-contract.md`

The task review bundle remains:

- `kaonlt_hardening_review.diff`

It may be regenerated for this repaired candidate but must never overwrite the
older `kaonlt_review.diff`.

Any change outside the complete revised allowlist is a hard stop.

## Authoritative starting facts

At the required starting commit:

- public root `AGENTS.md` is absent;
- `docs/memory/CURRENT.md` still contains consumed pre-push wording for Fix.5.6;
- Fix.5.6 actual-diff review and pushed-state review were completed before the
  latest supplied farm attempt;
- the latest supplied Fix.5.6 farm attempt ran at pushed source
  `ccaf15358efc205cd601aeccecdd1ab1b1360dff`;
- that attempt completed the analysis child but failed the tracked owner at
  `verify_artifacts` with literal failure
  `page_manifest_setting_invalid`;
- `analysis_started=true` and `analysis_completed=true`;
- collection and ZIP creation did not start;
- a fresh PDF and manifest existed, but because the owner gate failed they are
  not accepted validation artifacts and must not be promoted, packaged, or
  scientifically interpreted as a successful gate;
- the failed gate does not invalidate earlier accepted closures and does not
  promote Method A;
- no rerun is authorized by this hardening task.

These are continuity facts only. This task does not diagnose or repair the
scientific/runtime failure.

## Mandatory startup reading

Because root `AGENTS.md` is not yet public at the starting commit, use this
one-time bootstrap order:

1. inspect any existing local root `AGENTS.md` without assuming it is safe to
   publish;
2. read `docs/memory/CURRENT.md`;
3. read `docs/memory/MEMORY.md`;
4. read `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. read `docs/memory/USER.md`;
6. read `docs/memory/MAINTENANCE.md`;
7. read `docs/memory/CODEX.md`;
8. read `docs/memory/COMMUNICATION.md`;
9. read `docs/memory/TOOLS.md`;
10. read `docs/memory/LEARNINGS.md`;
11. read this contract in full.

After this task, tracked public root `AGENTS.md` becomes the first file in the
ordinary startup core and this bootstrap exception disappears.

## Objectives

Implement all agreed repository-side workflow protections from the 2026-10-02
hardening audit.

The task has seven objectives:

1. make root `AGENTS.md` a tracked public stable authority;
2. add explicit in-chat synchronization, health-check, evidence-label, drift,
   and re-anchor behavior;
3. make failed-gate handling stop forward motion instead of encouraging reruns,
   packaging, artifact consumption, or scientific interpretation;
4. make the ChatGPT -> Codex handoff exact: one standalone Markdown task
   contract plus one short inline Codex launch prompt;
5. remove personal/machine-specific filesystem path values from public workflow
   memory and replace them with symbolic environment roles;
6. reconcile `CURRENT.md` to the already-pushed Fix.5.6 state and the supplied
   failed owner gate without performing scientific diagnosis;
7. regenerate repository memory integrity metadata and produce one complete
   review bundle for independent ChatGPT actual-diff review.

Do not redesign the scientific architecture or repository memory hierarchy.

## Allowed paths

Only these paths may change:

- `AGENTS.md`
- `docs/memory/AGENTS.md`
- `docs/memory/README.md`
- `docs/memory/CURRENT.md`
- `docs/memory/MAINTENANCE.md`
- `docs/memory/USER.md`
- `docs/memory/CODEX.md`
- `docs/memory/COMMUNICATION.md`
- `docs/memory/TOOLS.md`
- `docs/memory/LEARNINGS.md`
- `docs/memory/phases/workflow-continuity-and-chat-health-hardening.md`
- `docs/memory/phases/workflow-continuity-and-chat-health-hardening-task-contract.md`
- `docs/memory/manifest.json`
- `tools/check_memory_health.py`
- `tools/memory_bootstrap.py`
- `testing/test_memory_health.py`

The task-contract file is an input. Do not rewrite it unless an objective
contradiction is discovered. If the contract itself is inconsistent, stop and
report the blocker.

## Frozen paths

Everything not explicitly allowlisted is frozen. The two explicitly
preserved pre-existing files below are frozen-in-place even though they are not
part of the task allowlist:

- `docs/memory/phases/workflow-continuity-hardening-task-contract.md`
- `kaonlt_review.diff`

Do not edit, move, delete, stage, include, or otherwise normalize either one.

In particular, do not modify:

- `src/**`
- `testing/**` except `testing/test_memory_health.py`
- `tools/**` except `tools/check_memory_health.py` and `tools/memory_bootstrap.py`
- `.codex/**`
- `.gitignore`
- `run_Prod_Analysis.sh`
- any farm profile
- any collector
- any scientific/runtime checker
- any runtime owner
- any plotting implementation
- any scientific producer/serializer/consumer
- `docs/memory/MEMORY.md`
- `docs/memory/handoffs/CURRENT_HANDOFF.md`
- `docs/memory/roadmap/**`
- `docs/memory/evidence/**`
- `docs/memory/decisions/**`
- `docs/memory/investigations/**`
- all existing scientific phase records, including Fix.5.4/Fix.5.5/Fix.5.6

No Project-level ChatGPT file or setting is edited by Codex.

## Scientific ownership and frozen interfaces

This task owns no scientific quantity.

Preserve exactly:

- random subtraction;
- dummy subtraction;
- slow-proton subtraction;
- pion-background treatment;
- the active `no_empirical_residual` profile;
- HGCer Method-A mathematics and authority boundaries;
- HGCer Method-B diagnostic-only role;
- SIMC production and normalization;
- yield extraction;
- cuts;
- templates;
- priors;
- binning;
- efficiencies;
- acceptance;
- uncertainty propagation;
- L/T separation;
- cross sections;
- F.6.2 accepted science;
- F.6.3 branch behavior;
- E.8 scientific/presentation source;
- Fix.5.4 numerical source;
- Fix.5.5 presentation source;
- Fix.5.6 owner/runtime source.

Method A remains detached/non-production pending an explicit validated F.6.4
decision. Method B remains diagnostic/cross-check only and numerically excluded.

## 1. Public root `AGENTS.md`

Create a concise tracked root `AGENTS.md`.

It is a stable repository-wide behavior/scientific-boundary authority. It must
not become a mutable state ledger.

It must not contain:

- a live/current HEAD;
- a current phase or fix;
- a current blocker;
- a current `NEXT`;
- temporary validation state;
- personal Windows/WSL paths;
- personal JLab filesystem paths;
- Downloads paths;
- ChatGPT Project-file paths;
- profanity or quoted correction markers.

### Required startup contract

At substantial-work start:

1. establish actual branch, HEAD, and worktree state;
2. read, in this exact order:
   - `AGENTS.md`
   - `docs/memory/CURRENT.md`
   - `docs/memory/MEMORY.md`
   - `docs/memory/handoffs/CURRENT_HANDOFF.md`
   - `docs/memory/USER.md`
3. follow CURRENT direct references and only task-relevant canonical records;
4. identify current status, accepted evidence, frozen interfaces, exact NEXT,
   and current workflow gate before proposing work.

Do not eagerly reload the entire history.

### Required authority order

Preserve:

current source/diff
-> newest applicable farm/runtime evidence
-> newest applicable validation/handoff artifact
-> tracked durable memory
-> older chat/history material.

`CURRENT.md` remains the sole ordinary active-state authority.

### Mature-analysis rule

Do not restart the analysis architecture, reinterpret accepted ownership,
reopen closed phases, redesign adjacent systems, or replace the established
workflow without a concrete blocker or new evidence.

### Runtime boundary

Do not claim ROOT/PyROOT, full `main.py`, procedure-PDF rendering, farm
integration, production behavior, or full kinematic validation without direct
applicable farm evidence.

### Source-changing workflow

Preserve the normal flow:

audit
-> one standalone repository task contract
-> one short inline Codex launch prompt
-> Codex implementation/local checks
-> ChatGPT actual-diff review
-> user commit/push
-> ChatGPT pushed-state synchronization
-> user farm validation when required
-> evidence review.

Codex does not commit, push, or run the farm.

### Professional workflow re-anchor

If the user explicitly states that the established workflow, source state,
pathing, analysis procedure, evidence interpretation, or requested task has
been misunderstood or that an established step was skipped:

1. stop extending the disputed workflow;
2. do not defend or rationalize the preceding approach;
3. re-establish live repository identity;
4. reread the startup core and task-relevant CURRENT references;
5. identify the exact incorrect assumption or violated gate;
6. resume from the earliest still-valid established gate;
7. preserve accepted work and closed phases;
8. do not invent a replacement procedure unless a concrete blocker invalidates
   the established one.

A re-anchor does not itself reopen science or authorize broad reconciliation.

### Evidence discipline

For consequential operational/scientific claims, distinguish:

- `SOURCE VERIFIED` — checked against current source/diff;
- `RUNTIME VERIFIED` — supported by applicable supplied runtime/farm evidence;
- `MEMORY ONLY` — supported by durable repository memory but not independently
  rechecked in the current chat;
- `INFERENCE` — reasoned conclusion rather than direct evidence;
- `NOT VERIFIED` — unknown or not checked.

Do not silently promote memory or inference into verified fact.

Claims such as `worktree clean`, `farm ready`, `gate passed`, `artifact valid`,
`runtime validated`, or `page correct` require concrete supporting evidence.
Otherwise state `NOT VERIFIED`.

### Failed-gate rule

A failed owner/checker/provenance/manifest/page/artifact gate stops forward
motion.

After a failed gate:

- do not issue the next farm command;
- do not reflexively rerun;
- do not package/copy/deliver failed-gate artifacts as accepted evidence;
- do not scientifically interpret failed-gate artifacts as though the gate
  passed;
- first identify the failed invariant and the earliest valid repair/debug gate.

Partial artifacts may be inspected only to diagnose the failed invariant and
must remain explicitly labeled failed/inadmissible evidence.

### Path discipline

Do not invent or substitute generic filesystem paths.

Machine/user-specific path values belong in the active external environment
configuration, not public repository workflow memory. If a required external
path is not available, request it rather than guessing.

## 1A. Reconcile the universal startup core

The repository currently has two different AGENTS surfaces:

- new public root `AGENTS.md`, intended as the stable repository-wide authority;
- existing `docs/memory/AGENTS.md`, which is a memory-specific record.

The first hardening candidate did not fully reconcile them. Existing
`docs/memory/AGENTS.md`, `docs/memory/README.md`, `MAINTENANCE.md`,
`tools/check_memory_health.py`, and `tools/memory_bootstrap.py` still encode a
startup core beginning with the memory-local AGENTS record.

Repair this so there is one unambiguous universal startup core:

1. repository-root `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Required behavior:

- root `AGENTS.md` is the only repository-wide startup/behavior authority;
- `docs/memory/AGENTS.md` remains tracked as memory-specific operating guidance
  but explicitly states that it is supplemental and not part of the universal
  five-file startup core;
- `docs/memory/README.md` and `MAINTENANCE.md` must identify the root AGENTS file
  explicitly rather than ambiguously referring to the memory-local file;
- `tools/memory_bootstrap.py` must report the universal five-file core with root
  `AGENTS.md` first;
- `tools/check_memory_health.py` must require root `AGENTS.md`, verify that the
  universal startup surfaces agree on the root-first sequence, and preserve all
  existing schema-3 integrity checks;
- `testing/test_memory_health.py` must cover the new root-first contract,
  including failure when root AGENTS is missing or when startup order drifts;
- existing `docs/memory/AGENTS.md` remains available for task-directed memory
  maintenance but is no longer a competing startup authority.

Do not alter scientific/runtime tooling. These two tools and one test are
memory-health/bootstrap infrastructure only.

## 2. `MAINTENANCE.md` — in-chat health protocol

Preserve existing memory roles, byte limits, manifest policy, and
memory-bracketed scientific-throughput rules.

Add an explicit assistant health protocol.

### Full health check triggers

Require a compact full health check:

- at the first substantial KaonLT task in a new chat/session;
- before a source-changing Codex handoff;
- before a farm-readiness declaration or farm command;
- before accepting/interpreting a returned farm artifact after a state change;
- before changing a phase/status/NEXT;
- after any failed owner/checker/provenance/manifest/page/artifact gate;
- after an explicit user workflow correction/re-anchor.

The health check must expose at least:

```text
KaonLT health check

repository:
  branch: <value>
  live HEAD: <sha>
  worktree: clean | dirty | NOT VERIFIED
  source checked live: YES | NO

state:
  active item: <item>
  authoritative NEXT: <gate>
  current gate: <gate>
  gate status: PASS | BLOCKED | NOT YET EVALUATED

evidence:
  newest accepted runtime evidence: <record/artifact/NONE>
  current source/diff reviewed: YES | NO
  runtime/farm validated in this chat: YES | NO
  unresolved contradiction: NONE | <exact issue>

ownership:
  scientific source frozen: YES | NO
  production logic touched: YES | NO
  Method A role: <role>
  Method B role: <role>

verification:
  SOURCE VERIFIED: <claims or NONE>
  RUNTIME VERIFIED: <claims or NONE>
  MEMORY ONLY: <claims or NONE>
  INFERENCE: <claims or NONE>
  NOT VERIFIED: <claims or NONE>

next action:
  <single next action>
```

The health check is an observability aid, not a competing active-state record.

### Periodic health pulse

For a long technical chat, require a shorter health pulse after roughly 6–10
substantive technical exchanges or at a major workflow-state transition,
whichever comes first:

```text
KaonLT health pulse
HEAD: <sha>
active item: <item>
NEXT: <gate>
current gate: <gate>
source/memory consistency: PASS | DRIFT
unverified assumptions: NONE | <list>
```

Do not emit the pulse for every ordinary explanatory turn.

### Automatic drift/re-anchor triggers

Force a re-anchor when any of these occurs:

- live HEAD differs from the identity currently being used;
- CURRENT and live source materially disagree;
- a gate/checker/owner returns failure;
- the worktree becomes unexpectedly dirty;
- an unrecorded filesystem path is proposed;
- source changes are proposed without the required contract;
- artifact interpretation begins before provenance/checker/page validation;
- the assistant contradicts an accepted status or NEXT;
- the user explicitly states workflow/state/path/task misalignment;
- the assistant cannot identify evidence supporting an operational claim.

Report `DRIFT` when appropriate and stop forward motion until the material
contradiction is resolved.

## 3. `USER.md` — exact delivery semantics

Preserve existing collaboration preferences.

Replace ambiguous wording that says plans and prompts are both standalone
Markdown files.

The stable rule must be:

- ChatGPT provides exactly one standalone Markdown **task contract** for a
  source-changing task;
- that contract is placed in the repository;
- ChatGPT then provides one short **inline Codex launch prompt** telling Codex
  to read repository memory and the contract;
- do not create a second standalone Codex-prompt Markdown file unless the user
  explicitly asks for one;
- for normal use the contract file, not a shell heredoc, is the transfer
  artifact;
- after Codex completes, ChatGPT reviews the actual diff rather than accepting
  the Codex summary.

Do not put the user's Downloads or workstation paths in public repository
memory.

Add the professional correction/re-anchor preference consistent with
`AGENTS.md`.

## 4. `CODEX.md` — contract and review handoff

Preserve existing source-changing ownership and hard boundaries.

Clarify:

- one standalone Markdown task contract;
- one short inline Codex launch prompt;
- no duplicate prompt artifact by default;
- never use a giant shell heredoc as the normal task-contract delivery method;
- contract placement is external/session-specific and must not be guessed in
  public repo memory;
- actual-diff review remains mandatory;
- new/untracked files must be represented in the review bundle;
- a failed hard stop returns `BLOCKED` and Codex must not invent a workaround.

Do not alter the scientific contract template or tool code in this task.

## 5. `COMMUNICATION.md` — failed-gate and artifact-consumption discipline

Preserve the existing farm-readiness audit and `tcsh` requirement.

Strengthen it with the following rules.

### Farm readiness is gated

Before a farm command, state either:

`Farm readiness: PASS`

or

`Farm readiness: BLOCKED`

and identify the evidence/source owners supporting that result.

Do not infer readiness merely because a previous run worked or a PDF/ZIP path
exists.

### Failed gate stops progression

If an owner/checker/provenance/manifest/page/artifact verification step fails:

1. stop the operation at that gate;
2. request/inspect the exact failure evidence;
3. determine whether produced artifacts are admissible for diagnosis only;
4. do not move to collection, packaging, copying, scientific interpretation, or
   another run until the failed invariant is understood.

### No rerun reflex

A failed farm attempt does not automatically authorize another run.

The next action must follow the failure class:

- source/state failure -> source/state repair;
- artifact-schema/provenance failure -> artifact/provenance diagnosis;
- renderer/page failure -> renderer/presentation diagnosis;
- scientific closure failure -> producer/payload/consumer closure audit.

### Verification before consumption

Encode the intended order:

```text
run
-> owner/checker PASS
-> provenance/freshness PASS
-> structured payload/page-manifest PASS
-> rendered-page inspection
-> packaging/handoff
-> scientific interpretation
```

A file existing on disk is not gate acceptance.

### Rendered content is part of acceptance

For procedure-PDF work, `complete=true`, analysis completion, or file existence
is insufficient. Required page IDs, represented children, renderer failures,
provenance, required comparison objects, and actual visual legibility must be
inspected as applicable.

### Suspicious physics output

If scalar-yield changes and displayed spectra appear inconsistent, or a required
comparison such as SIMC is absent, do not explain the physics first.

Trace:

producer
-> serializer/sidecar/checkpoint
-> payload
-> consumer
-> renderer

and test producer-owned closure/integration semantics before interpretation.

### Tracked owner remains authoritative

Once a reviewed owner/wrapper exists, do not improvise an alternate manual
analysis/package workflow around it unless a concrete source-level blocker
requires an explicit repair contract.

## 6. `TOOLS.md` — externalize actual path values

Keep durable generic operations and the durable fact that the JLab farm shell
is `tcsh`.

Remove personal/user-specific filesystem values from the public workflow
reference, including exact workstation and JLab user paths.

Replace exact path values with symbolic roles such as:

```text
<KAONLT_REPO_ROOT>
<KAONLT_ARTIFACT_ROOT>
<KAONLT_BUNDLE_ROOT>
```

State that their concrete values are supplied by the active external
environment/ChatGPT Project configuration.

Do not invent replacement values.

Retain semantic ownership of where owner status/logs/ZIPs live by referring to
the symbolic roots, not personal literal paths.

## 7. `LEARNINGS.md` — durable lessons from the failure

Add concise reusable lessons, not chronology:

- a failed downstream gate invalidates forward progression even when the child
  analysis completed;
- child-process success and owner success are distinct;
- artifact existence is not artifact admissibility;
- never rerun by reflex after a failed gate;
- investigate the first failed invariant;
- suspicious scalar/histogram disagreement requires closure before physics
  interpretation;
- generated PDF existence does not justify packaging or handoff after a failed
  manifest/page gate;
- long-chat health checks and explicit evidence labels reduce silent drift.

Do not quote profanity or private chat text.

## 8. `CURRENT.md` reconciliation

Make only the minimum current-state correction needed for accuracy.

Preserve the scientific objective and all accepted closures.

Record that:

- Fix.5.6 is still `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`;
- its source at
  `ccaf15358efc205cd601aeccecdd1ab1b1360dff`
  had independent actual-diff review and pushed-state review before the supplied
  farm attempt;
- the supplied farm attempt completed the analysis child but failed the owner at
  `verify_artifacts` with `page_manifest_setting_invalid`;
- collection/ZIP did not begin;
- generated PDF/manifest are failed-gate diagnostic artifacts, not accepted
  validation evidence;
- no scientific interpretation, Method-A promotion, canonical-five closure, or
  production change follows.

Replace the consumed old NEXT with one push-stable scientific NEXT:

after workflow hardening is complete and the user resumes scientific work,
diagnose the exact `page_manifest_setting_invalid` failure and the corresponding
page/payload provenance before any rerun, packaging, or scientific
interpretation.

Do not diagnose the failure in this task.
Do not authorize a farm command.

## 9. Hardening phase record

Create:

`docs/memory/phases/workflow-continuity-and-chat-health-hardening.md`

It must record:

- purpose and scope;
- starting source identity;
- public `AGENTS.md` addition;
- in-chat health protocol;
- evidence labels;
- re-anchor protocol;
- failed-gate stop rule;
- exact task-contract/inline-prompt delivery convention;
- path externalization;
- current-state reconciliation;
- frozen scientific source;
- deterministic validation;
- no farm/runtime claim.

Write it in intended final form for a candidate that may be committed only
after independent ChatGPT actual-diff review passes.

Do not claim farm validation.

## 9A. Hardening phase-record push stability

Repair the phase-record opening so it does not simultaneously claim
`SOURCE REVIEWED` and describe independent actual-diff review as still pending.

The intended committed status is:

`SOURCE REVIEWED`

That status is warranted only if the repaired candidate receives independent
ChatGPT actual-diff PASS before the user commits it. Write the phase record in
push-stable form: no consumed sentence may imply that the independent review is
still pending after publication.

Do not claim pushed-state synchronization before it occurs.

## Before/after behavior

### Before

A new or long-running chat can:

- start without public `AGENTS.md`;
- silently use stale mutable Project checkpoint information;
- drift without an explicit health pulse;
- promote memory/inference into operational fact;
- guess a familiar-looking path;
- continue after a failed farm gate;
- rerun reflexively;
- consume/copy/interpret artifacts before gate acceptance;
- confuse child-analysis completion with owner success;
- deliver a task contract through a fragile nested heredoc;
- interpret ambiguous `USER.md` wording as requiring multiple Markdown files.

### After

The public repository provides:

- one tracked stable root startup authority;
- explicit health checks and periodic health pulses;
- explicit evidence labels;
- professional workflow re-anchor behavior;
- automatic material-drift triggers;
- failed-gate stop/no-rerun discipline;
- verification-before-consumption ordering;
- exact one-contract-plus-inline-prompt delivery semantics;
- no-guess path discipline with concrete path values externalized;
- accurate current Fix.5.6 failed-gate continuity;
- unchanged scientific/runtime implementation.

## Positive checks

Verify all of the following:

1. root `AGENTS.md` is an intended tracked public file;
2. `AGENTS.md` contains the five-file startup order;
3. `AGENTS.md` contains evidence labels;
4. `AGENTS.md` contains the professional re-anchor rule;
5. `AGENTS.md` contains the failed-gate stop rule;
6. `AGENTS.md` contains no mutable current HEAD/phase/NEXT;
7. `MAINTENANCE.md` contains full health-check triggers and the periodic pulse;
8. `MAINTENANCE.md` defines material drift/re-anchor triggers;
9. `USER.md` specifies one standalone task contract plus one short inline Codex
   prompt;
10. `CODEX.md` matches that delivery convention;
11. `COMMUNICATION.md` states PASS/BLOCKED farm readiness and failed-gate stop;
12. `COMMUNICATION.md` contains no-rerun and verification-before-consumption
    behavior;
13. `TOOLS.md` uses symbolic environment roots rather than personal literal
    path values;
14. `LEARNINGS.md` records the reusable failed-gate/health lessons;
15. CURRENT records the pushed Fix.5.6 source and supplied failed
    `verify_artifacts/page_manifest_setting_invalid` attempt accurately;
16. CURRENT authorizes no rerun or artifact promotion;
17. all scientific/runtime source remains unchanged.

## Negative checks

Changed public workflow/memory files, including this contract, must not contain
literal personal or user-specific filesystem values. Prohibit concrete values
for all of these categories:

- Windows or WSL user-home paths;
- personal Downloads paths;
- workstation-specific repository clone paths;
- user-specific Jefferson Lab repository/output roots;
- user-specific transferable-bundle roots;
- guessed home-relative repository paths.

Use symbolic roles or generic category descriptions only.

Do not require removing historical/scientific path literals from frozen
scientific/evidence records outside this task.

The changed public workflow files must not contain profanity or quoted private
correction markers.

`AGENTS.md` must not contain a 40-character live source SHA as current state,
an active phase/fix, or substantive current NEXT.

## Regression checks

Confirm no diff under:

```text
src/
docs/memory/evidence/
docs/memory/decisions/
docs/memory/investigations/
docs/memory/roadmap/
```

Under `tools/`, only these may change:

```text
tools/check_memory_health.py
tools/memory_bootstrap.py
```

Under `testing/`, only this may change:

```text
testing/test_memory_health.py
```

Confirm no diff to any existing scientific Fix.5.4/Fix.5.5/Fix.5.6 phase record.

Confirm no farm profile, collector, wrapper, owner, checker, renderer,
scientific source, or scientific/runtime test changed.

Confirm:

- no status upgrade to `CLOSED / RUNTIME VALIDATED` for Fix.5.6;
- no Method-A production promotion;
- no Method-B numerical role;
- no normalization/cut/template/prior/binning/yield change;
- no new farm authorization.

## Local validation

Discover a working Python interpreter according to repository policy.

After versioned memory changes:

```text
<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/check_memory_health.py --root .
<PYTHON> -B tools/memory_bootstrap.py --self-test
<PYTHON> -B testing/test_memory_health.py
git diff --check
```

Also run deterministic textual checks sufficient to verify:

- all required health/re-anchor/failure rules are present;
- no changed public workflow record, including this contract, contains literal
  personal filesystem values;
- root `AGENTS.md` has no mutable current state;
- the universal startup core is root-AGENTS-first in docs, health checks and
  bootstrap output;
- `docs/memory/AGENTS.md` is explicitly supplemental rather than a competing
  startup authority;
- changed paths are exactly within the revised allowlist.

Report:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
warning classification: <blocking | nonblocking>
manifest check: PASS | FAIL
```

Do not use `--fail-on-warning` unless an actual warning is materially blocking
under existing policy.

## Diff audit and review bundle

Report:

```text
git status --short --untracked-files=all
git diff --stat
git diff --name-only
git diff --check
```

Do not stage files merely to make them reviewable.

Create exactly one temporary review bundle in the repository root:

`kaonlt_hardening_review.diff`

It must contain:

1. the complete tracked diff for every modified allowlisted tracked file;
2. complete `git diff --no-index /dev/null ...` additions for every intended
   new/untracked allowlisted file, including root `AGENTS.md`, the hardening
   phase record, and this contract when untracked.

Use `|| true` only for `git diff --no-index`, where exit status 1 means a
difference was emitted.

Do not overwrite any unrelated pre-existing review bundle.

The hardening review bundle must exclude:

- `docs/memory/phases/workflow-continuity-hardening-task-contract.md`
- `kaonlt_review.diff`

The hardening review bundle is temporary and must not be committed.

## User-controlled Git handoff

After deterministic checks pass, provide exact non-executed user commands for:

```text
git add <exact intended paths only>
git commit -m "Harden KaonLT workflow continuity and chat health"
git push origin test
```

The exact `git add` command must exclude:

- `docs/memory/phases/workflow-continuity-hardening-task-contract.md`
- `kaonlt_review.diff`
- `kaonlt_hardening_review.diff`

Do not execute them.

The user alone commits/pushes.

## Farm boundary

No farm execution belongs to this task.

Do not provide a farm command.

No ROOT/PyROOT, full `main.py`, procedure-PDF, ZIP, page-manifest, numerical
closure, visual-quality, or production-validation claim is established here.

The existing supplied failed Fix.5.6 attempt remains failed/inadmissible for
validation. Hardening does not convert it into accepted evidence.

## Acceptance criteria

PASS only if:

- public root `AGENTS.md` is ready to track;
- the universal startup order is root AGENTS -> CURRENT -> MEMORY -> handoff -> USER;
- `docs/memory/AGENTS.md` is supplemental and no longer a competing startup authority;
- memory health/bootstrap tooling enforces/reports the root-first startup core;
- deterministic memory-health tests cover the repaired startup contract;
- full chat health checks and periodic pulses are durable;
- evidence labels are durable;
- professional re-anchor behavior is durable;
- material drift automatically stops forward motion;
- failed gates stop progression and rerun reflexes;
- verification precedes artifact consumption/interpretation;
- procedure-PDF acceptance includes rendered-content checks;
- suspicious physics output requires closure before interpretation;
- Codex delivery is exactly one standalone task contract plus one inline launch
  prompt;
- public workflow memory contains no personal concrete path values owned by the
  external environment;
- CURRENT reflects the pushed Fix.5.6 source and failed owner gate accurately;
- no scientific/runtime/test source changed;
- the earlier superseded hardening contract remains untouched;
- the pre-existing `kaonlt_review.diff` remains untouched;
- manifest check passes;
- memory health has no hard failure;
- diff audit is allowlist-clean;
- one complete `kaonlt_hardening_review.diff` exists.

## Hard stop

Stop and report `BLOCKED` without broadening scope if:

- branch is not `test`;
- local HEAD or `origin/test` is not
  `ccaf15358efc205cd601aeccecdd1ab1b1360dff`;
- any worktree change exists outside the revised allowlist, the current review
  bundle, or the explicitly preserved pre-existing files;
- either preserved pre-existing file changes during the task;
- repository evidence contradicts the supplied failed-gate continuity facts;
- a required change would touch a frozen file;
- scientific/runtime/test changes appear;
- a hard memory-health failure cannot be repaired inside the allowlist;
- making a local `AGENTS.md` public would expose private/machine-specific
  content that cannot be cleanly separated.

Do not reset, stash, clean, commit, push, run the farm, or invent a workaround.
