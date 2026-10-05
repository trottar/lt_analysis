# E.8.4 Left/lowe current-lineage Method-A diagnostic runtime-closure and workflow-correction task

## 1. Task class and objective

This is a **closure task plus narrow workflow-correction checkpoint** under
`docs/memory/CODEX.md`.

It replaces an earlier local, uncommitted draft of this same contract that
incorrectly required Codex to locate or independently revalidate farm artifacts
that ChatGPT had already reviewed and accepted.

Objectives:

1. record the accepted Jefferson Lab runtime evidence for the detached
   `Q4p4W2p74 / Left / lowe` current-lineage Method-A diagnostic;
2. update active/durable repository memory after the completed scientific
   interpretation gate;
3. record the workflow failures identified by the user so they are not repeated;
4. leave one accurate, push-stable substantive NEXT;
5. change no scientific/analysis source and authorize no farm execution.

This task must not alter Method A mathematics, Method B, baseline pion
subtraction, random/dummy subtraction, slow-proton treatment, SIMC, yields,
cuts, templates, priors, binning, efficiencies, acceptance, uncertainties,
L/T separation, or cross sections.

---

## 2. Exact starting identity and local-state rule

Required branch:

```text
test
```

Required current local `HEAD` and local `origin/test`:

```text
aad27a4d1639eef188dc61563fd3615682835835
```

Commit subject:

```text
Add detached current-lineage Method-A diagnostic
```

Before editing, report:

```bash
git branch --show-current
git rev-parse HEAD
git rev-parse origin/test
git status --short --untracked-files=all
```

The local worktree may already contain this contract at:

```text
docs/memory/phases/e8-4-left-lowe-method-a-current-lineage-diagnostic-runtime-closure-task-contract.md
```

because the user moved the earlier blocked draft there. That pre-existing local
contract is expected and must be replaced by this corrected version.

Hard stop if:

- branch is not `test`;
- local `HEAD` differs from the required SHA;
- local `origin/test` differs from the required SHA;
- unrelated local state cannot be preserved safely;
- any path outside the allowlist below would need modification;
- current source or accepted repository evidence materially contradicts the
  runtime facts or scientific boundaries below.

Do not reset, clean, stash, overwrite unrelated files, update refs, stage,
commit, push, or run the Jefferson Lab farm.

---

## 3. Mandatory startup reads

Read in exact repository order:

1. `AGENTS.md`
2. `docs/memory/CURRENT.md`
3. `docs/memory/MEMORY.md`
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`
5. `docs/memory/USER.md`

Then read:

- `docs/memory/MAINTENANCE.md`
- `docs/memory/CODEX.md`
- `docs/memory/COMMUNICATION.md`
- `docs/memory/LEARNINGS.md`
- `docs/memory/TOOLS.md`
- `docs/memory/templates/CODEX_CONTRACT.md`

Then read only task-relevant records:

- `docs/memory/phases/e8-4-left-lowe-method-a-current-lineage-diagnostic-measurement-task-contract.md`
- `docs/memory/evidence/f6-2-scientific-runtime-closure.md`
- `docs/memory/evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md`
- `docs/memory/evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/phases/phase-f6-method-a-production-promotion.md`

Do not reopen closed F stages.

---

## 4. Evidence handoff rule for this closure

The returned farm artifacts were already supplied to ChatGPT and independently
reviewed before this closure task.

For this task:

- Codex is **not responsible for independently locating, opening, hashing, or
  revalidating the raw returned JSON/PDF**.
- Codex must **not** block because those raw artifacts are absent from the local
  workstation/repository.
- Codex must **not** require a second external "ChatGPT acceptance review" file.
- The exact accepted runtime facts and artifact identities are provided below so
  Codex can record them in repository memory/evidence.
- The new repository evidence record must state transparently that the runtime
  artifacts were user-supplied and independently reviewed by ChatGPT; Codex is
  recording that accepted evidence and did not itself rerun or revalidate the
  farm artifacts.

The task contract itself is not runtime evidence. The runtime evidence is the
user-supplied artifact pair identified by exact SHA-256 below, whose acceptance
was completed by ChatGPT before this closure task.

This is the established closure workflow:

```text
user farm run
-> artifacts returned to ChatGPT
-> ChatGPT artifact/provenance/page/science review PASS
-> one tracked closure contract
-> Codex repository-memory/evidence reconciliation
-> ChatGPT actual-diff review
-> user commit/push
-> ChatGPT pushed-state synchronization
```

Do not invent an additional external review artifact between ChatGPT evidence
acceptance and the tracked closure contract.

---

## 5. Accepted returned artifact identity

Accepted artifacts:

```text
Q4p4W2p74_Left_lowe_kaon_pion-background_hgcer-method-a-current-lineage-diagnostic.json
Q4p4W2p74_Left_lowe_kaon_pion-background_hgcer-method-a-current-lineage-diagnostic.pdf
```

Accepted SHA-256:

```text
JSON:
ecd547de15aa19a50b73bfd8d9dd1c159027c6438ca38b79bc45ff69a8d0e8c7

PDF:
0f3e0d016ce895c6ff1aa390f28c3f3e85c99e2cc8733c7fa0a176aec48f0a73
```

Accepted JSON size:

```text
690916 bytes
```

Accepted PDF size:

```text
71818 bytes
```

Runtime provenance recorded by the JSON:

```text
git_head:
aad27a4d1639eef188dc61563fd3615682835835

git_status_short:
 M src/models/xmodel_kaon_pl.f
?? src/kaon/functions/Q4p4W2p74.model

runtime_claim:
detached diagnostic only; no production promotion
```

Do **not** rewrite the farm checkout as clean.

The unrelated dirty paths above were not consumed by the detached JSON
diagnostic, which consumed exact F.1/F.3/F.4 artifact inputs.

---

## 6. Accepted current-lineage authority

The accepted exact current-lineage identities are:

```text
candidate validation source head:
b349967c0d4210a78b144ce6134d3c1f15970245

Left/lowe F.1 SHA-256:
10086d7d4c42389c9fdd16471c49d59cd189a30980767bcdb87b85914169ef95

candidate F.3 source SHA-256:
eeb480d9b7ddffab1e97c7f9f27c3099d0780c2bbb4d42a0ee68c077a5752b3d

candidate F.3 map fingerprint:
3a9787fc58d26cc0816012bd1b637ad0c8f201b448d54cb1625a841131154728

candidate F.4 source SHA-256:
1d545924eba89c7f9ffa28028e307aca9b434a89beec06863cf2893887b6b902

candidate F.4 correction fingerprint:
bce0cc12ef283b0c361818b436c4e190fc4ea5165436f6c36912cfbb6de53368

candidate F.4 artifact fingerprint:
4c935271a0b2723b58b02cc34d82a1e7d5757f7cd36b7707f400c2b28893668a

diagnostic fingerprint:
2c68804d60aaea423cf3dffcd55c3e324c812b97298149df8599ba193ee23616

artifact fingerprint:
7b6c1b244eb15d5d262391ca4fb288d896b9776989c9f658f2877b05dd02205e
```

The F.4 reproduction reported:

```text
payload_identical = true
fingerprint = bce0cc12ef283b0c361818b436c4e190fc4ea5165436f6c36912cfbb6de53368
```

The input-authority gate reported:

```text
accepted_authority_match = true
```

Important correction from the blocked draft:

- the exact Left/lowe F.1 hash contains `...bcdb87...`;
- the exact untracked farm-local path is
  `src/kaon/functions/Q4p4W2p74.model`.

Do not reproduce the prior transcription mistakes.

---

## 7. Accepted artifact gates

Record as `RUNTIME VERIFIED` at the exact detached
`Q4p4W2p74 / Left / lowe` diagnostic scope:

### JSON

The returned JSON passed independent review for:

- exact SHA-256 match;
- schema/version validity;
- `non_authoritative = true`;
- artifact fingerprint recomputation;
- diagnostic fingerprint recomputation;
- exact scope `Q4p4W2p74 / Left / lowe`;
- parents `t0`, `t1`, `t2`;
- primary parent `t1`;
- exact current-lineage F.1/F.3/F.4 authority;
- exact accepted-authority match;
- exact shared F.4 reproduction and payload identity;
- no Method-B numerical dependency;
- no production-object mutation;
- no alternative correction construction;
- aggregate-only persisted diagnostics;
- no forbidden persisted per-event factor/identity table;
- finite numerical payload;
- source-decomposition closure;
- child/MM/coordinate aggregate reconstruction;
- signed canonical-parent closure.

### PDF

The returned PDF passed independent review for:

- exact SHA-256 match;
- exactly eight pages;
- no missing page;
- no clipping;
- no overlap;
- no broken plot/glyph;
- legible page content.

Page roles:

1. authority/scope/limitations;
2. weak-positive/control response;
3. `hgcer3` topology;
4. signed/absolute parent normalization;
5. t1 source decomposition/tails;
6. t1 child/MM sensitivity;
7. t1 coordinate sensitivity;
8. interpretation boundary.

The smaller weak-positive curves reflect lower statistics and are not a
rendering failure.

---

## 8. Accepted Measurement-A facts

Required limitation flags from the accepted JSON:

```text
operational_kaon_pid_npe_zero = true
direct_hgc_free_true_pion_tag_present_in_consumed_artifacts = false
direct_pion_to_kaon_misid_calibration_performed = false
absolute_misid_probability_constructed = false
hgc_free_pion_tag_available_in_consumed_artifacts = false
proxy_validity_established = false
proxy_validity_evidence_label = NOT VERIFIED
```

These flags limit the **absolute mis-identification interpretation**. They do
not invalidate Method A's accepted relative acceptance-correlated refinement.

Training counts:

```text
t0:
  control-response = 10716
  weak-positive    = 71

t1:
  control-response = 18749
  weak-positive    = 283

t2:
  control-response = 22356
  weak-positive    = 222
```

Training-control versus physical-control identity matching is exact:

```text
t0: 10716 matched; 0 training-only; 0 physical-only
t1: 18749 matched; 0 training-only; 0 physical-only
t2: 22356 matched; 0 training-only; 0 physical-only
```

Persisted matched-population coordinate/response quantile differences are zero.

For t1, the accepted relative-response summaries include approximately:

```text
weak-positive:
  median = 0.872
  p90    = 1.629
  p99    = 3.295

control-response:
  median = 0.729
  p90    = 1.373
  p99    = 3.022
```

These are relative response/topology measurements, not absolute
`P(NPE=0 | true pion, x)`.

---

## 9. Accepted Measurement-B facts

For t1:

```text
B_signed      = 0.3374726818838243
A_abs         = 0.3845152767785784
cancellation  = 0.8776574098983241

U_signed      = 0.29922739636212625
V_abs         = 0.34446789697425667

N_signed      = 0.8866714623885792
N_abs         = 0.8958497042306518

N_signed / N_abs
              = 0.9897547079619178

relative normalization difference
              = -0.010245292038082221
              ≈ -1.02%
```

Accepted t1 signed-parent closure residual:

```text
0.0
```

Same-setting context:

```text
t0 relative signed-vs-absolute normalization difference:
+0.0027001080906590147 ≈ +0.27%

t2 relative signed-vs-absolute normalization difference:
+0.03129771404630625 ≈ +3.13%
```

Therefore t1 is not the parent with the largest signed-versus-positive-support
normalization contrast.

### t1 source support

```text
prompt:
  count = 18749
  absolute-support fraction = 0.9383929586334123 ≈ 93.84%

dummy:
  count = 93
  absolute-support fraction = 0.04915800899474474 ≈ 4.92%

random:
  count = 963
  absolute-support fraction = 0.012013286056093246 ≈ 1.20%

dummy-random:
  count = 4
  absolute-support fraction = 0.0004357463157497171
```

Each source class is internally sign-pure in this decomposition; cancellation
occurs between source classes.

### t1 support/OOD

```text
in_support_count = 19788
ood_count        = 21
ood_fraction     = 0.001060124185976119 ≈ 0.106%
```

A broad F.3 out-of-support application failure is not supported.

### t1 response/correction tail

```text
raw r:
  p50 = 0.7291584557619877
  p90 = 1.3661220249087134
  p99 = 2.988399963863005
  max = 53.331269127265315

accepted C_signed:
  p50 = 0.8223547127565471
  p90 = 1.5407307924725082
  p99 = 3.3703576698100095
  max = 60.14772256636942
```

The finite long tail remains real. Do not clip, cap, winsorize, hide, or
reinterpret it as a different correction.

---

## 10. Accepted F.6.2 / current-lineage continuity

Preserve the accepted F.6.2 status:

```text
CLOSED / RUNTIME VALIDATED
```

F.6.2 is the accepted **acceptance-correlated Method-A refinement validation**.
Its scientific evidence includes:

- normalized low-response / baseline / Method-A shapes;
- F.3 construction variables;
- independent acceptance variables including `SHMS_xptar` and `SHMS_yptar`;
- phi and missing mass;
- MM x acceptance localization;
- support/OOD;
- diagnostic variance/effective-statistics information;
- fixed kaon-window pion-change metrics;
- paired bootstrap intervals.

Do not reduce the scientific interpretation to "HGC response only".

The accepted current-baseline comparator established:

```text
F.2 scientific payload match = true
F.3 scientific payload match = true
first scientifically changed stage = F.4
```

For the three Left/lowe parents:

- current-lineage t1 differs from historical accepted F.4 only at
  floating-point scale;
- t0 and t2 contain the substantive current-baseline F.4 changes.

Therefore accepted historical F.6.2 acceptance/MM evidence remains directly
relevant to present t1 interpretation.

Do not automatically extend that exact t1 continuity claim to current-lineage
t0/t2.

RF is additional corroborating PID/timing information **where experimentally
available at low epsilon only**. It is not present at high epsilon and therefore
must not become a universal Method-A requirement. The returned diagnostic did
not itself perform an RF corroboration.

---

## 11. Scientific interpretation and status decision

Record carefully by evidence class.

### SOURCE VERIFIED

- Method A is an acceptance-correlated PID refinement, not merely an HGCer
  response curve.
- F.3 construction uses
  `SHMS_delta`, `P_hgcer_xAtCer`, `P_hgcer_yAtCer`.
- accepted F.6.2 validation uses independent acceptance and MM information.
- baseline `w0` retains pion-control -> kaon-background transfer ownership.
- F.4 preserves the signed canonical-t parent.
- Method A remains detached/non-production.
- Method B remains diagnostic/cross-check only and numerically excluded.

### RUNTIME VERIFIED

At the exact returned Left/lowe diagnostic scope:

- JSON/PDF acceptance gates passed;
- current t1 application is overwhelmingly in F.3 support;
- training-control and physical-control populations match exactly;
- t1 signed-versus-positive-support normalization differs by only about 1.02%;
- t1 parent closure is exact at recorded precision;
- accepted F.6.2 remains closed at its recorded scope;
- comparator F.2/F.3 scientific equality remains accepted;
- comparator t1 F.4 change is floating-point scale.

### INFERENCE / scientific decision

Record:

```text
The current Q4p4W2p74 / Left / lowe / t1 diagnostic does not support strong
signed-normalization pathology, broad F.3 OOD application, or
training-control -> physical-control population mismatch as explanations of the
large redistribution. Together with the accepted F.6.2 acceptance/MM validation
and t1 current-lineage continuity, no concrete evidence from this diagnostic
warrants a Method-A redesign or replacement correction for the t1 issue.
```

Do not promote this inference to a measured detector hardware cause.

The narrow current-lineage diagnostic itself may be recorded as:

```text
CLOSED / RUNTIME VALIDATED
```

at this exact measurement scope only.

### NOT VERIFIED

Preserve:

- direct absolute `P(NPE=0 | true pion, x)` calibration;
- absolute pion-to-kaon HGC mis-ID probability from F.1/F.3/F.4;
- direct proof that weak-positive response uniquely determines latent zero
  response;
- a specific PMT/mirror/optical/hardware cause;
- RF corroboration in the returned diagnostic;
- Method-A production correctness outside accepted detached scopes;
- current-lineage t0/t2 inheritance of the historical t1 F.4 continuity claim;
- canonical-five current-lineage closure;
- production promotion;
- absolute-SIMC amplitude interpretation.

---

## 12. User-identified workflow failures that must be recorded

This task must update durable workflow memory so these failures do not recur.

### Failure A — wrong post-farm closure dependency

The blocked draft incorrectly made Codex require local access to the raw
returned JSON/PDF after ChatGPT had already independently accepted those
artifacts.

Correct rule:

```text
For a closure task after ChatGPT has accepted supplied farm evidence, the
tracked closure contract may carry the exact accepted artifact identities,
measurements, evidence labels, and closure decision. Codex records/reconciles
that accepted evidence and must not invent a new local-artifact prerequisite
unless the contract explicitly requires Codex to perform a separate local
artifact check.
```

The contract is not runtime evidence; it is the scoped authority for recording
runtime evidence that ChatGPT already accepted.

### Failure B — invented second external review artifact

After Codex correctly blocked on the faulty contract, ChatGPT invented a
separate "ChatGPT acceptance review" Markdown file in Downloads and tried to
make that another Codex dependency.

Correct rule:

```text
Do not create a second standalone acceptance-review artifact between completed
ChatGPT farm-evidence review and the normal tracked closure contract unless the
user explicitly requests such an artifact or a concrete repository workflow
requires it.
```

The usual flow is farm artifacts -> ChatGPT review -> tracked closure contract.

### Failure C — broken generated-contract delivery handoff

ChatGPT initially told the user to "place" a generated task contract in the
repository, then added an unnecessary intermediate `ls` step, instead of using
the established environment-resolved download-to-repository move handoff.

Correct collaboration rule:

```text
When ChatGPT generates the one task-contract file as a downloadable artifact,
give the exact active-environment move command from the configured Downloads
location to the exact repository contract path. Do not make the user manually
place it or perform an unnecessary existence-check step when a direct scoped
move is sufficient.
```

Do not store concrete workstation/Downloads paths in public repository memory;
store only this behavioral rule. Concrete paths remain in active external
environment configuration.

### Failure D — manual provenance transcription defects

The blocked/replacement draft introduced transcription errors:

- incorrect Left/lowe F.1 hash text by dropping the `d` in `...bcdb87...`;
- shortened dirty-path text instead of preserving
  `src/kaon/functions/Q4p4W2p74.model`.

Correct rule:

```text
When carrying accepted hashes, fingerprints, or provenance paths into a closure
contract/evidence record, copy them exactly from the accepted artifact/source
and independently cross-check them before handoff. Do not manually normalize,
shorten, or retype provenance strings.
```

These are workflow failures, not scientific evidence changes.

---

## 13. Required workflow-memory edits

### `docs/memory/CODEX.md`

Add a concise post-farm closure-evidence handoff rule:

- after ChatGPT independently accepts supplied farm artifacts, a closure
  contract may carry the exact accepted facts to Codex for repository
  reconciliation;
- Codex must not require raw artifacts locally unless the contract explicitly
  assigns Codex a local artifact-validation task;
- do not invent a second external acceptance-review file;
- preserve transparent attribution that Codex records previously accepted
  runtime evidence rather than claiming it independently validated the farm.

Keep this general and concise.

### `docs/memory/USER.md`

Under delivery/collaboration preferences, add the stable rule:

- when ChatGPT generates the single downloadable task contract, use the active
  external environment to provide the exact direct move command from configured
  Downloads into the exact repository path;
- do not tell the user to place the file manually;
- do not add unnecessary intermediate file-existence checks when the direct
  scoped move is sufficient;
- concrete personal paths remain external and must not be written into public
  repository memory.

### `docs/memory/LEARNINGS.md`

Add concise reusable lessons covering:

- post-farm closure evidence must not be converted into a redundant Codex
  local-artifact prerequisite after ChatGPT acceptance;
- avoid invented second acceptance-review artifacts;
- exact hashes/fingerprints/provenance paths must be copied and cross-checked,
  never casually retyped or shortened;
- generated contract delivery follows the established direct move handoff using
  external environment values.

Do not add phase chronology to LEARNINGS.

---

## 14. Required scientific-memory/evidence edits

### New evidence record

Create:

```text
docs/memory/evidence/e8-4-left-lowe-method-a-current-lineage-diagnostic-runtime-closure-2026-10-04.md
```

It must compactly include:

- artifact filenames, exact SHA-256 and sizes;
- exact farm provenance including the dirty worktree paths;
- exact F.1/F.3/F.4 authority identities;
- JSON/PDF acceptance gates;
- key Measurement-A/B results;
- t1 source support, OOD and tail facts;
- limitation flags;
- accepted F.6.2 acceptance/MM continuity;
- comparator t1 continuity versus t0/t2 distinction;
- evidence labels and scientific inference;
- explicit no-production/no-promotion boundary;
- remaining `NOT VERIFIED` items;
- attribution that artifacts were supplied by the user and independently
  accepted by ChatGPT before this Codex closure task.

Do not claim Codex independently revalidated the artifacts.

### `docs/memory/CURRENT.md`

Compact the consumed farm-pending wording.

It must state:

- current-lineage Left/lowe diagnostic:
  `CLOSED / RUNTIME VALIDATED` at the exact narrow scope;
- link to the new evidence record;
- t1 signed normalization, broad OOD, and training/application mismatch are not
  supported as explanations by the returned measurements;
- accepted F.6.2 acceptance/MM evidence and comparator t1 continuity are part of
  the interpretation;
- no correction redesign is warranted by this diagnostic;
- direct absolute mis-ID calibration and hardware cause remain `NOT VERIFIED`;
- Method A remains detached/non-production;
- Method B remains diagnostic only;
- canonical-five remains `DEFERRED`;
- final E.8/F.6.4 remain `BLOCKED`;
- absolute-SIMC interpretation remains separately `BLOCKED`.

Keep CURRENT compact and below repository size limits.

### `docs/memory/MEMORY.md`

Update durable statements made stale by the accepted diagnostic:

- do not continue saying current-lineage t1 signed/absolute support, source
  decomposition, response/correction tails, child/MM redistribution,
  coordinate dependence, or support/OOD are unmeasured;
- preserve relative-response != absolute leakage probability;
- preserve baseline `w0` transfer ownership;
- preserve acceptance/MM as part of Method-A scientific validation;
- preserve direct absolute proxy validity and hardware cause as `NOT VERIFIED`;
- preserve the exact current-lineage t1 continuity boundary and do not extend it
  to t0/t2.

Add no mutable NEXT or permanent live-HEAD ledger.

### `docs/memory/roadmap/STATUS.md`

Consume:

```text
narrow Left/lowe farm diagnostic measurement
scientific interpretation / method-selection decision
```

Record that this diagnostic did not warrant a correction-design implementation.

Preserve:

```text
later canonical-five reconsideration
-> final E.8
-> F.6.4 explicit production-promotion decision
```

Canonical-five remains `DEFERRED` unless the user explicitly lifts that status
after this closure.

Do not downgrade or reopen historical closed phases.

---

## 15. Sole push-stable NEXT

Leave one ordinary CURRENT NEXT equivalent to:

```text
NEXT — after this closure is independently diff-reviewed, user-committed/pushed,
and pushed-state synchronized, return to the user-deferred current-lineage
canonical-five reconsideration gate. Do not start canonical-five source/farm
work unless the user explicitly authorizes lifting DEFERRED. If authorized,
first audit the exact five-setting current-lineage authority/validation scope
before writing any implementation or farm contract.
```

Commit/push itself must not be CURRENT's substantive NEXT.

---

## 16. Allowed versioned paths

Only these versioned paths may change:

```text
docs/memory/phases/e8-4-left-lowe-method-a-current-lineage-diagnostic-runtime-closure-task-contract.md
docs/memory/evidence/e8-4-left-lowe-method-a-current-lineage-diagnostic-runtime-closure-2026-10-04.md
docs/memory/CURRENT.md
docs/memory/MEMORY.md
docs/memory/roadmap/STATUS.md
docs/memory/CODEX.md
docs/memory/USER.md
docs/memory/LEARNINGS.md
docs/memory/manifest.json
```

Everything else is frozen.

The contract and evidence record are intended tracked files. The corrected
contract replaces the earlier local uncommitted draft at the same path.

---

## 17. Frozen paths and scientific ownership

Do not modify any file under:

```text
src/
testing/
farm_env/
```

Also freeze:

```text
AGENTS.md
docs/memory/handoffs/CURRENT_HANDOFF.md
docs/memory/MAINTENANCE.md
docs/memory/COMMUNICATION.md
docs/memory/TOOLS.md
docs/memory/decisions/
docs/memory/investigations/
```

unless explicitly allowlisted above.

No changes to:

- F.1/F.2/F.3/F.4/F.5/F.6.1/F.6.2/F.6.3 calculations;
- baseline `w0`;
- random/dummy/slow-proton subtraction;
- Method B;
- production profiles;
- SIMC;
- yields;
- uncertainties;
- cuts/templates/priors/binning;
- efficiencies/acceptance;
- L/T separation;
- cross sections.

---

## 18. Before / after behavior

### Before

- diagnostic source is pushed at `aad27a4d...`;
- farm JSON/PDF were returned and accepted by ChatGPT;
- tracked CURRENT still says the farm diagnostic is pending;
- durable memory still carries pre-measurement unresolved wording;
- no repository evidence record owns the accepted 2026-10-04 diagnostic;
- workflow memory does not explicitly prevent the four user-identified
  handoff/transcription failures above;
- the local contract path may contain the earlier blocked draft.

### After

- one evidence record owns the accepted diagnostic result;
- CURRENT accurately reflects the consumed farm/scientific gate;
- MEMORY accurately reflects measured versus still-unverified t1 questions;
- roadmap advances to the canonical-five reconsideration dependency while
  preserving DEFERRED;
- CODEX/USER/LEARNINGS record the workflow corrections;
- no science source changes;
- no production behavior changes;
- no farm run is required;
- Method A remains detached/non-production.

---

## 19. Local validation and memory health

Discover `<PYTHON>` using repository conventions.

After memory edits:

```bash
<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/update_memory_manifest.py --root . --check
<PYTHON> -B tools/check_memory_health.py --root .
<PYTHON> -B tools/memory_bootstrap.py --root . --json
git diff --check
```

Report:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
warning classification: <none or each warning>
manifest check: PASS | FAIL
git diff --check: PASS | FAIL
```

Hard failures and blocking warnings stop the task.

---

## 20. Positive checks

Confirm:

```text
branch/head/origin exact
only allowlisted files changed
corrected contract replaced blocked draft
no raw-artifact local prerequisite
no second external acceptance-review dependency
workflow failures recorded in CODEX/USER/LEARNINGS
exact JSON/PDF hashes preserved
exact Left/lowe F.1 hash includes ...bcdb87...
exact farm dirty path includes src/kaon/functions/Q4p4W2p74.model
diagnostic CLOSED / RUNTIME VALIDATED only at exact Left/lowe scope
F.6.2 acceptance/MM validation preserved
t1 comparator continuity preserved without extending to t0/t2
signed-normalization/support/population results stated quantitatively
absolute proxy validity remains NOT VERIFIED
hardware cause remains NOT VERIFIED
RF is low-epsilon-only corroboration, not universal requirement
no Method-A redesign authorized
no production promotion
canonical-five remains DEFERRED
final E.8/F.6.4 remain BLOCKED
absolute-SIMC blocker retained
one push-stable NEXT
manifest and memory-health checks pass
```

---

## 21. Negative checks / forbidden shortcuts

Confirm:

```text
no src/ edit
no testing/ edit
no farm_env/ edit
no production mutation
no new correction
no new normalization
no child renormalization
no clipping/winsorization
no Method-B numerical use
no claim of absolute pion mis-ID calibration
no claim of specific hardware cause
no claim RF validated the returned artifact
no automatic canonical-five expansion
no automatic F.6.4 promotion
no SIMC amplitude interpretation
no requirement to locate raw farm JSON/PDF locally
no new external ChatGPT acceptance-review file
no farm command/run
no staging
no commit/push by Codex
no reset/clean/stash
```

Do not describe the farm worktree as clean.

Do not convert `INFERENCE` into `RUNTIME VERIFIED`.

Do not downgrade closed phases.

---

## 22. Farm-validation boundary

This closure task does **not** lead directly to a farm run.

No farm command is authorized.

The next possible farm-bound task exists only after:

1. this closure diff passes ChatGPT review;
2. the user commits/pushes;
3. ChatGPT performs pushed-state synchronization;
4. the user explicitly authorizes lifting canonical-five `DEFERRED`;
5. a later audit/contract establishes the exact current-lineage five-setting
   validation chain.

---

## 23. Complete actual-diff review bundle

Create temporary repository-root:

```text
kaonlt_review.diff
```

It must contain the complete diffs for all intended changed/new files,
including complete additions for the new evidence record and this corrected
contract.

A suitable pattern is:

```bash
git diff -- \
  docs/memory/CURRENT.md \
  docs/memory/MEMORY.md \
  docs/memory/roadmap/STATUS.md \
  docs/memory/CODEX.md \
  docs/memory/USER.md \
  docs/memory/LEARNINGS.md \
  docs/memory/manifest.json \
  > kaonlt_review.diff

git diff --no-index -- /dev/null \
  docs/memory/phases/e8-4-left-lowe-method-a-current-lineage-diagnostic-runtime-closure-task-contract.md \
  >> kaonlt_review.diff || true

git diff --no-index -- /dev/null \
  docs/memory/evidence/e8-4-left-lowe-method-a-current-lineage-diagnostic-runtime-closure-2026-10-04.md \
  >> kaonlt_review.diff || true
```

If the corrected contract path is tracked in local state rather than untracked,
use the appropriate complete tracked diff instead; the final review bundle must
contain its entire effective change.

Do not stage merely for review.

Do not include unrelated local state.

Stop before commit/push.

---

## 24. Acceptance criteria

PASS only if:

- exact starting identity passes;
- only allowlisted paths change;
- the prior false local-artifact dependency is removed;
- no second external acceptance-review artifact is introduced;
- all four user-identified workflow failures are durably recorded in the
  appropriate canonical memory roles;
- exact artifact/provenance strings are correct;
- runtime facts and scientific inference remain distinct;
- accepted F.6.2 acceptance/MM evidence is preserved;
- t1 continuity is not generalized to t0/t2;
- no Method-A redesign is invented;
- Method A remains detached/non-production;
- canonical-five remains deferred pending explicit user authorization;
- all manifest/memory/diff checks pass;
- complete `kaonlt_review.diff` is produced;
- no farm execution, staging, commit, or push occurred.

If all pass, return a completed **closure candidate** ready for independent
ChatGPT actual-diff review.

---

## 25. Hard stop

Return `BLOCKED` and stop if:

- branch/HEAD/origin differs from the exact required identity;
- unrelated local state cannot be preserved;
- repository source/evidence contradicts the accepted facts above;
- closure would require analysis-source changes;
- closure would require inventing an absolute pion mis-ID calibration;
- closure would require claiming a hardware cause;
- closure would require lifting canonical-five deferral without explicit user
  authorization;
- a memory edit creates material ambiguity in source identity, accepted
  evidence/status, frozen interfaces, active ownership, blocker state, or NEXT;
- manifest or memory-health hard checks fail.

Do not invent a workaround.
