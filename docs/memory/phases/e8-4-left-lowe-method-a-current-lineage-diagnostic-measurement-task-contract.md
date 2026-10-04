# KaonLT — Left/lowe current-lineage Method-A diagnostic-measurement task

## 1. Task class and objective

This is a **detached diagnostic implementation task** under `docs/memory/CODEX.md`.

It follows the completed Gate-4 source/science audit and the completed,
user-authorized external PID/background-methods review. The next scientific
question is deliberately split into two independent measurements:

1. **F.3 response/proxy question**
   Determine what the accepted current-lineage F.3 weak-positive HGC response
   actually demonstrates about relative pion→kaon HGC mis-identification
   topology across `hgcer3 = (SHMS_delta, P_hgcer_xAtCer, P_hgcer_yAtCer)`,
   without inventing a truth label or treating `NPE=0` as a pion sample.

2. **F.4 signed-normalization question**
   Quantify whether the large `Q4p4W2p74 / Left / lowe / t1` redistribution is
   associated with strong signed-weight cancellation or a large difference
   between the existing signed F.4 normalization measure and a descriptive
   positive-support comparator.

This task creates a **new detached, non-authoritative, current-lineage
diagnostic only**. It does not alter F.1, F.3, F.4, F.6.3, baseline pion
subtraction, production yields, or any production weight.

The implementation must answer with measurements and explicit limitations,
not with a replacement correction or a promotion decision.

---

## 2. Exact starting identity and hard precondition

Required branch:

```text
test
```

Required local `HEAD` and local `origin/test` at task start:

```text
6542bba01946e370842566f24e2c8861071ecd09
```

Before editing, report:

```bash
git branch --show-current
git rev-parse HEAD
git rev-parse origin/test
git status --short --untracked-files=all
```

Hard stop if:

- branch is not `test`;
- `HEAD` differs from the required SHA;
- local `origin/test` differs from the required SHA;
- unexpected unrelated worktree state cannot be preserved safely;
- this contract is absent from the exact path specified below;
- current source materially contradicts the scientific/source facts in this
  contract.

Do not reset, clean, stash, overwrite unrelated files, update refs, commit,
push, or run the Jefferson Lab farm.

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

Then read only task-relevant records/source:

- `docs/memory/investigations/e8-4-left-lowe-method-a-detector-response-source-science-audit-2026-10-03.md`
- `docs/memory/investigations/e8-4-left-lowe-method-a-external-pid-background-methods-review-2026-10-03.md`
- `docs/memory/roadmap/STATUS.md`
- `docs/memory/evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md`
- `docs/memory/evidence/f6-3-e8-4-left-lowe-runtime-closure-2026-10-01.md`
- `src/utility/background_config.py`
- `src/cuts/pion_hgcer_method_a_acceptance_contract.py`
- `src/cuts/pion_hgcer_method_a_acceptance_map.py`
- `src/cuts/pion_hgcer_method_a_parent_preserving_correction.py`
- `src/cuts/pion_hgcer_method_a_parallel_full_procedure.py`
- `src/cuts/pion_hgcer_method_a_reweighting_validation.py`
- `testing/analyze_pion_hgcer_method_a_parent_preserving_correction.py`
- `testing/test_f6_3_parallel_full_procedure_method_a.py`
- `testing/test_pion_hgcer_method_a_parent_preserving_correction.py`

Authority remains:

```text
current source/diff
-> newest applicable accepted farm/runtime evidence
-> tracked durable repository memory
-> external literature
-> inference
```

---

## 4. SOURCE VERIFIED scientific/source baseline

### 4.1 Operational HGC PID semantics

Current source defines:

```text
kaon PID category:
  P_hgcer_npeSum == 0

pion tree:
  P_hgcer_npeSum > 0

physical pion control:
  P_hgcer_npeSum > 2
```

`NPE=0` is therefore the **kaon-selected PID category**, not an observed pion
population.

Conditional on true species,

```text
P(NPE=0 | true pion, x)
```

means the pion→kaon HGC mis-identification probability into the kaon-selected
sample.

Do not use the phrase “NPE=0 pion sample/class/population” unless explicitly
qualified as a latent true-pion component inside the kaon-selected PID category.

### 4.2 F.1 populations

The accepted dual-population F.1 structure is:

```text
training:
  prompt / noRF / nommcuts / P_hgcer_npeSum > 0

weak-positive training class:
  0 < P_hgcer_npeSum <= 2

control-response training class:
  P_hgcer_npeSum > 2

physical application:
  authoritative physical-pion-control population with P_hgcer_npeSum > 2
```

Training and physical application are distinct.

The current F.1 sidecar contains HGC response and trajectory/kinematic
coordinates but does **not** itself provide an HGC-independent true-species
pion tag for the kaon-selected `NPE=0` category.

Therefore this task must not claim direct measurement of
`P(NPE=0 | true pion, x)` from F.1/F.3/F.4 alone.

### 4.3 F.3 meaning

The accepted current-lineage F.3 basis is:

```text
hgcer3 =
  SHMS_delta
  P_hgcer_xAtCer
  P_hgcer_yAtCer
```

F.3 is a support-aware **relative weak-positive detector-response
representation**.

It is not:

- an absolute pion→kaon mis-ID probability;
- a true-pion efficiency calibration;
- a production correction;
- a replacement for `w0`.

### 4.4 F.4 mathematics

For each canonical-`t` parent:

```text
b_j = signed_source_coefficient_j * w0_j
r_j = F.3 relative response

B = sum_j b_j
U = sum_j b_j * r_j

N_signed = U / B
C_signed,j = r_j / N_signed
```

and by construction:

```text
sum_j b_j * C_signed,j = B
```

This preserves the signed canonical-`t` parent only.

It does not independently preserve a `(t,phi)` child, MM region, Lambda-window
yield, or individual event contribution.

### 4.5 Baseline `w0`

`w0` already owns the accepted pion-control → kaon-background transfer from the
baseline pion component model.

This task must not create another transfer normalization or interpret an
alternative response normalization as a new physical pion-background
normalization.

---

## 5. Exact current-lineage authority

This task is **not** allowed to consume historical F.6.1 F.3/F.4 authority as
though it were the present current-baseline lineage.

Use the current F.6.3 candidate authority exposed by:

```text
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
```

At the required starting source, the current Q4p4W2p74 candidate identity is:

```text
candidate validation source head:
  b349967c0d4210a78b144ce6134d3c1f15970245

Left-lowe F.1 source SHA-256:
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
```

The implementation must obtain/validate authority through existing
current-lineage source helpers and must fail closed on any mismatch. Do not
duplicate these hashes into a second independent authority system unless a
test fixture needs an exact expected value.

Use the existing current-lineage helpers where appropriate:

```text
accepted_f6_3_artifact_paths(...)
load_accepted_f6_3_authority(...)
reconstruct_transient_factor_map(...)
```

and the accepted shared F.4 calculator/review-data path.

No filename fallback, nearest match, stale historical artifact, or
“best available” artifact is permitted.

---

## 6. Scope

Target scope:

```text
kinematic: Q4p4W2p74
phi setting: Left
epsilon token: lowe
canonical t parents: t0, t1, t2
primary scientific target: t1
```

`t0` and `t2` are required as same-setting context/control parents.

Do not expand to the canonical five settings in this task.

Canonical-five expansion remains `DEFERRED`.

---

## 7. Required new detached diagnostic module

Create:

```text
src/cuts/pion_hgcer_method_a_current_lineage_diagnostic.py
```

The module must:

- be pure Python/NumPy only;
- import no ROOT/PyROOT;
- traverse no ROOT tree;
- consume exact F.1/F.3/F.4 current-lineage JSON artifacts;
- validate current-lineage authority before producing scientific output;
- reuse existing accepted F.3/F.4 evaluators rather than reimplementing their
  scientific calculations;
- reproduce the accepted F.4 payload exactly before measuring it;
- expose only aggregate/descriptive diagnostic output;
- never modify F.3 or F.4;
- never construct a production correction;
- never alter baseline `w0`;
- never use Method B numerically.

Suggested schema names:

```text
pion_hgcer_method_a_current_lineage_diagnostic/v1
pion_hgcer_method_a_current_lineage_diagnostic_artifact/v1
```

Use repository conventions for canonical JSON fingerprints.

---

## 8. Required measurement A — weak-positive F.3 proxy evidence

This is a **proxy-evidence diagnostic**, not direct pion→kaon mis-ID
calibration.

### 8.1 Required limitation status

Persist explicit flags equivalent to:

```text
operational_kaon_pid_npe_zero: true
direct_hgc_free_true_pion_tag_present_in_consumed_artifacts: false
direct_pion_to_kaon_misid_calibration_performed: false
absolute_misid_probability_constructed: false
proxy_validity_established: false
```

The exact field names may follow repository naming conventions, but the meaning
must be unambiguous and tested.

The diagnostic must state that direct proxy validity remains `NOT VERIFIED`
from the consumed F.1/F.3/F.4 artifacts alone.

### 8.2 Required descriptive populations

For each Left/lowe canonical-`t` parent, summarize the prompt training
populations:

```text
weak-positive:
  0 < NPE <= 2

control-response:
  NPE > 2
```

Required aggregate counts and feature summaries:

```text
count
NPE distribution summary
SHMS_delta summary
P_hgcer_xAtCer summary
P_hgcer_yAtCer summary
F.3 raw relative-response summary
F.3 support / OOD counts where defined
```

At minimum use robust distribution summaries sufficient to reveal tails:

```text
min
p01
p10
p25
median
p75
p90
p99
max
```

Do not add clipping/winsorization.

### 8.3 Positive-NPE diagnostic bands

For descriptive visualization only, summarize prompt positive-HGC training in
fixed diagnostic bands:

```text
0 < NPE <= 1
1 < NPE <= 2
2 < NPE <= 4
NPE > 4
```

These are **diagnostic bins**, not new PID cuts, production categories, or
optimization boundaries.

For each band and canonical-`t` parent, summarize:

```text
event count
F.3 raw relative-response distribution
hgcer3 coordinate distributions
```

This may reveal whether weaker observed pion response occupies the same
acceptance regions assigned larger F.3 relative response, but it must not be
described as direct mis-ID validation.

### 8.4 Training-control vs physical-control population shift

For prompt events that can be identity-matched between the F.1 training
control population and physical application control population, report:

```text
matched count
training-control-only count
physical-control-only count
```

and aggregate `hgcer3` / F.3-response distribution differences.

This is a population-transport/support diagnostic only.

Do not create a calibration weight.

### 8.5 HGC-independent tag availability

Within the **consumed artifacts and fields only**, explicitly audit whether an
HGC-independent true-pion tag is already represented.

Expected current result is that F.1/F.3/F.4 do not provide such a truth/species
tag.

Do not search raw ROOT branches or invent a tag in this task.

Persist:

```text
hgc_free_pion_tag_available_in_consumed_artifacts: false
```

if source confirms that result.

A broader search for a scientifically defensible independent pion tag may be
authorized later; absence here is not proof that no such information exists
anywhere in KaonLT.

### 8.6 Forbidden proxy claims

The module/report must not claim:

- F.3 is an absolute pion→kaon mis-ID map;
- weak-positive events are kaons;
- `NPE=0` is a pion sample;
- weak-positive/control separation proves true-pion leakage;
- positive-NPE data uniquely determine the zero-response probability;
- a particular Poisson/truncated response model is required.

---

## 9. Required measurement B — F.4 signed-normalization sensitivity

This measurement uses the **existing accepted F.4 current-lineage raw response
and signed baseline contributions**.

For each Left/lowe canonical-`t` parent compute and persist:

```text
B_signed      = sum_j b_j
A_abs         = sum_j |b_j|
cancellation  = |B_signed| / A_abs

U_signed      = sum_j b_j * r_j
V_abs         = sum_j |b_j| * r_j

N_signed      = U_signed / B_signed
N_abs         = V_abs / A_abs
```

`N_abs` is a **descriptive positive-support comparator only**.

It must never be used to construct or persist an alternative correction applied
to data, templates, yields, or production.

Persist descriptive contrasts such as:

```text
normalization_ratio = N_signed / N_abs
relative_normalization_difference = (N_signed - N_abs) / N_abs
```

only when denominators are finite and strictly positive. Fail closed rather
than silently substituting a value.

### 9.1 Exact F.4 reproduction

Before computing the comparison, reproduce the persisted current-lineage F.4
payload exactly through the shared calculator.

Require exact fingerprint/payload identity.

No diagnostic result is valid if the current F.4 artifact cannot be reproduced
exactly.

### 9.2 Source-separated decomposition

For each canonical-`t` parent and each authoritative source label present,
including the expected:

```text
prompt
rand
dummy
dummy_rand
```

persist:

```text
event count
sum b_j
sum |b_j|
sum b_j*r_j
sum |b_j|*r_j
signed fraction of B_signed
absolute-support fraction of A_abs
source cancellation ratio = |sum b_j| / sum |b_j|
```

Do not relabel these source classes as physical pion species.

They are authoritative statistical/source classes.

### 9.3 Event-distribution summaries

Persist aggregate summaries, not event identities, for:

```text
r_j                           raw F.3 relative response
C_signed,j                    accepted F.4 correction factor
w0_j                          baseline pion transfer weight
b_j                           signed baseline contribution
|b_j|                         absolute support
```

Use robust tails including p01/p10/p25/p50/p75/p90/p99/min/max.

Also persist:

```text
in-support count
OOD count
OOD fraction
```

for the F.3 application population.

### 9.4 Child / MM / coordinate sensitivity

For each parent, and with t1 emphasized, persist aggregate descriptive
measurements of the accepted signed baseline and accepted F.4-adjusted
contributions versus:

```text
phi child
analysis_MM
SHMS_delta
P_hgcer_xAtCer
P_hgcer_yAtCer
```

Requirements:

- canonical phi edges must be the existing frozen child geometry;
- MM and coordinate display bins are diagnostic-only;
- no adaptive optimization;
- no child renormalization;
- no fitting to make curves agree;
- no clipping/winsorization;
- no new physics cut.

For fixed diagnostic binning, choose deterministic source-owned bins in the
new module and document them in the artifact `display_policy`.

For each bin, record at minimum:

```text
count
signed baseline sum
signed adjusted sum
signed delta
absolute baseline support
```

where meaningful.

### 9.5 No verdict metric

Do not create:

- an overall Method-A score;
- a pass/fail “scientific validity” classifier;
- a winner between signed and absolute normalization;
- a promotion recommendation;
- a tuned threshold selected from these measurements.

The output is descriptive evidence for later scientific interpretation.

---

## 10. Persisted artifact contract

The diagnostic JSON must be aggregate-only.

It must persist:

```text
schema/version
status / available / reason
scope identity
current-lineage input identities and hashes/fingerprints
F.4 exact-reproduction result
PID semantic statement
measurement-A limitation flags
proxy-evidence aggregate tables
measurement-B signed/absolute parent tables
source decompositions
distribution summaries
binned child/MM/hgcer3 summaries
display policy
provenance
diagnostic fingerprint
```

It must **not** persist:

```text
event entry_index
(source_label, entry_index) identity tables
per-event correction factors
per-event alternative factors
raw event feature rows
a production weight
a new response matrix
an absolute pion mis-ID probability
a replacement F.4 normalization
```

Transient event-level objects may exist only inside the process long enough to
construct aggregate diagnostics.

Add an explicit recursive persistence guard/test that rejects forbidden
event-level identity/factor payloads if practical under existing repository
patterns.

---

## 11. Deterministic analyzer and review PDF

Create:

```text
testing/analyze_pion_hgcer_method_a_current_lineage_diagnostic.py
```

The analyzer must:

- be deterministic;
- require an existing output/artifact directory and `--kinematic`;
- support only the exact intended current-lineage candidate identity;
- require/select `Left / lowe` explicitly, either by frozen arguments or exact
  CLI choices;
- load artifacts through current-lineage authority helpers;
- fail closed on stale/missing/wrong artifacts;
- refuse to overwrite existing output;
- write exactly one diagnostic JSON and one review PDF;
- record Git provenance using existing testing/analyzer conventions;
- import no ROOT.

Freeze deterministic output names:

```text
Q4p4W2p74_Left_lowe_kaon_pion-background_hgcer-method-a-current-lineage-diagnostic.json
Q4p4W2p74_Left_lowe_kaon_pion-background_hgcer-method-a-current-lineage-diagnostic.pdf
```

If repository naming conventions require hyphen normalization, keep the exact
two filenames source-owned and test them. Do not discover names dynamically
from directory contents.

### 11.1 Required PDF content

The review PDF should be compact and scientific, approximately 6–8 pages:

1. **Authority / scope / limitations**
   - exact current F.1/F.3/F.4 identities;
   - `NPE=0` kaon PID semantics;
   - no HGC-free pion truth tag in consumed artifacts;
   - no direct absolute mis-ID calibration.

2. **Weak-positive vs control response**
   - counts and F.3 raw-response distributions;
   - positive-NPE band comparison.

3. **`hgcer3` proxy topology**
   - weak-positive/control distributions or ratio-style descriptive overlays
     versus `delta`, HGC x, HGC y;
   - support/OOD.

4. **Parent normalization table**
   - t0/t1/t2:
     `B_signed`, `A_abs`, cancellation ratio, `U_signed`,
     `V_abs`, `N_signed`, `N_abs`, normalization contrast.

5. **t1 source decomposition**
   - prompt/rand/dummy/dummy_rand contributions;
   - raw `r` and accepted `C_signed` tails/ECDFs.

6. **t1 child/MM sensitivity**
   - baseline vs accepted F.4-adjusted signed contribution by canonical phi
     child and MM.

7. **t1 coordinate sensitivity**
   - baseline/adjusted contribution versus `delta`, HGC x, HGC y.

8. **Interpretation boundary**
   - measured facts;
   - `INFERENCE`;
   - `NOT VERIFIED`;
   - no correction/promotion decision.

Do not include Method B.

Do not include an absolute-SIMC amplitude comparison.

---

## 12. Required tests

Create:

```text
testing/test_pion_hgcer_method_a_current_lineage_diagnostic.py
```

At minimum test:

### Authority / scope

- wrong kinematic rejected;
- wrong setting rejected;
- stale/wrong F.1/F.3/F.4 hash or fingerprint rejected;
- exact current-lineage F.4 reproduction required;
- historical F.6.1 artifact is not silently accepted.

### PID semantics

- artifact says `NPE=0` is kaon PID category;
- no standalone observed NPE=0 pion class is emitted;
- direct HGC-free pion tag is absent from consumed artifacts;
- no absolute mis-ID probability is constructed;
- weak-positive/control definitions are exact.

### Measurement A

- correct weak-positive/control counts on synthetic fixtures;
- deterministic NPE-band boundaries;
- response/coordinate summaries deterministic;
- no calibration weight generated;
- no proxy-validity PASS verdict generated.

### Measurement B

Use synthetic rows where signed cancellation is controlled exactly.

Test at least:

1. no-cancellation case:
   `|B|/A_abs = 1`;

2. cancellation case:
   `|B|/A_abs < 1`;

3. hand-calculable:
   `B`, `A_abs`, `U`, `V_abs`, `N_signed`, `N_abs`,
   source decomposition and contrasts;

4. zero/invalid denominators fail closed;

5. descriptive `N_abs` never becomes an applied event factor;

6. t1 can differ strongly while exact accepted signed parent closure remains
   unchanged.

### Persistence / frozen boundaries

- no event identity table persisted;
- no per-event correction factor persisted;
- no Method-B numerical field/dependency;
- no production object mutation;
- no child renormalization;
- no yield/cross-section construction;
- deterministic filename contract;
- deterministic JSON fingerprint for same inputs.

### Analyzer

Add focused analyzer tests if needed to validate:

- deterministic paths;
- overwrite refusal;
- one JSON + one PDF;
- failure on wrong current-lineage authority;
- successful fake/pure-Python rendering without ROOT.

---

## 13. Required regressions

Run the new tests plus existing current-lineage/F.4 regressions at minimum:

```bash
<PYTHON> -B testing/test_pion_hgcer_method_a_current_lineage_diagnostic.py
<PYTHON> -B testing/test_f6_3_parallel_full_procedure_method_a.py
<PYTHON> -B testing/test_pion_hgcer_method_a_parent_preserving_correction.py
<PYTHON> -B testing/test_analyze_pion_hgcer_method_a_parent_preserving_correction.py
```

Run `py_compile` on every new/modified Python source.

If these tests use repository-specific invocation conventions, follow the
existing conventions exactly and report the actual commands.

Do not run ROOT/PyROOT or full `main.py` locally merely to satisfy this task.

---

## 14. Memory updates

This task owns only compact active-state/roadmap updates.

### 14.1 `docs/memory/CURRENT.md`

Allow update, but **compact rather than append**.

`CURRENT.md` is already near the 8-KiB soft threshold. After this task it must
remain below the threshold.

Replace the consumed research-checkpoint live wording with:

- research checkpoint consumed;
- detached current-lineage diagnostic implementation is the active work item;
- exact two-question scope:
  1. F.3 weak-positive relative proxy evidence;
  2. F.4 signed-normalization sensitivity;
- no direct HGC-free true-pion tag exists in consumed F.1/F.3/F.4 artifacts;
- direct proxy validity therefore remains unverified until independent tagging
  or other defensible evidence exists;
- current task is diagnostic/non-production only;
- canonical-five remains `DEFERRED`;
- final E.8/F.6.4 and absolute-SIMC remain separately `BLOCKED`.

If local implementation and deterministic tests pass, use the work-state:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

for this new diagnostic implementation.

The sole push-stable substantive NEXT should then be:

```text
NEXT — after independent actual-diff review, user commit/push, and pushed-state
synchronization: run one narrow Q4p4W2p74 / Left / lowe farm diagnostic gate to
produce the current-lineage diagnostic JSON/PDF, then review the measurements
before any method-selection or correction-design decision.
```

Do not place the actual farm command in CURRENT.

### 14.2 `docs/memory/roadmap/STATUS.md`

Update only the active dependency line minimally:

```text
external methods review consumed
-> detached current-lineage diagnostic implementation
-> narrow Left/lowe farm diagnostic measurement
-> scientific interpretation / method-selection decision
-> implementation contract only if warranted
-> later canonical-five reconsideration
-> final E.8
-> F.6.4 promotion decision
```

Do not change accepted historical closure scopes.

### 14.3 `docs/memory/MEMORY.md`

Frozen for this task.

The durable research lessons are already recorded.

### 14.4 Manifest

Regenerate `docs/memory/manifest.json` after all intended repository-memory
changes.

---

## 15. Allowed versioned paths

Only these versioned paths may change:

```text
src/cuts/pion_hgcer_method_a_current_lineage_diagnostic.py
testing/analyze_pion_hgcer_method_a_current_lineage_diagnostic.py
testing/test_pion_hgcer_method_a_current_lineage_diagnostic.py
docs/memory/CURRENT.md
docs/memory/roadmap/STATUS.md
docs/memory/manifest.json
docs/memory/phases/e8-4-left-lowe-method-a-current-lineage-diagnostic-measurement-task-contract.md
```

The first three and the contract are intended new tracked files.

Everything else is frozen.

---

## 16. Explicit frozen paths / scientific ownership

Do not modify:

```text
src/utility/background_config.py
src/cuts/pion_hgcer_method_a_acceptance_contract.py
src/cuts/pion_hgcer_method_a_acceptance_map.py
src/cuts/pion_hgcer_method_a_parent_preserving_correction.py
src/cuts/pion_hgcer_method_a_tphi_propagation.py
src/cuts/pion_hgcer_method_a_reweighting_validation.py
src/cuts/pion_hgcer_method_a_parallel_full_procedure.py
src/cuts/pion_hgcer_transfer.py
src/cuts/rand_sub.py
src/cuts/particle_subtraction.py
src/cuts/pion_component_subtraction.py
src/cuts/calculate_yield.py
src/cuts/ave_per_bin.py
src/cuts/full_background_subtraction_plots.py
main analysis orchestration
testing existing analyzers/tests except the one new analyzer/test
testing validation-bundle collectors/profiles/packagers
run_Prod_Analysis.sh
farm_env/
docs/memory/MEMORY.md
docs/memory/evidence/
docs/memory/decisions/
docs/memory/sources/SOURCE_INDEX.md
docs/memory/handoffs/CURRENT_HANDOFF.md
docs/memory/USER.md
docs/memory/MAINTENANCE.md
docs/memory/CODEX.md
AGENTS.md
```

Preserve:

```text
baseline production authoritative
no_empirical_residual
zero legacy empirical residual scales
random/dummy subtraction frozen
slow-proton subtraction frozen
baseline pion subtraction frozen
F.3 frozen
F.4 frozen
F.6.3 private branch frozen
Method A detached/non-production
Method B diagnostic/cross-check only and numerically excluded
SIMC production/normalization frozen
yield formulas frozen
uncertainties frozen
cuts/templates/priors/binning/efficiencies/acceptance frozen
L/T separation frozen
cross sections frozen
canonical-five expansion DEFERRED
final E.8 BLOCKED
F.6.4 BLOCKED
absolute-SIMC interpretation separately BLOCKED
```

---

## 17. Before / after behavior

### Before

At starting HEAD:

- current-lineage F.3/F.4 artifacts exist and are authority-pinned;
- accepted Left/lowe F.6.3 evidence proves branch execution, real child changes,
  and signed parent closure;
- no detached current-lineage diagnostic decomposes signed versus absolute
  support at t1;
- no current-lineage diagnostic explicitly separates weak-positive proxy
  evidence from direct pion→kaon mis-ID calibration;
- no HGC-free true-pion tag is represented in consumed F.1/F.3/F.4 artifacts.

### After

After this task:

- one new pure-Python detached diagnostic can consume exact current-lineage
  F.1/F.3/F.4 artifacts;
- it reproduces F.4 exactly before analysis;
- it measures weak-positive/control response topology without claiming direct
  mis-ID calibration;
- it measures signed cancellation and signed-vs-positive-support normalization
  sensitivity;
- it emits aggregate JSON + PDF only;
- no production object, correction, yield, or Method-B quantity changes;
- local deterministic source is ready for one narrow farm measurement gate;
- scientific interpretation remains pending returned farm evidence.

---

## 18. Positive checks

Confirm all of the following:

```text
exact start HEAD/origin/test
only allowlisted paths changed
current F6.3 candidate authority fail-closed
F4 exact reproduction required
NPE=0 correctly described as kaon PID category
weak-positive = 0<NPE<=2
physical pion control = NPE>2
no HGC-free pion tag claimed from F1/F3/F4
direct absolute pion->kaon mis-ID calibration = false
proxy validity not declared established
signed B/A/U quantities measured
positive-support comparator measured descriptively
source decomposition preserved
t1 highlighted with t0/t2 context
no alternative correction applied
aggregate-only persistence
no Method B numerical dependency
no ROOT import
no production mutation
CURRENT below 8 KiB
manifest regenerated/checks pass
```

---

## 19. Negative checks / forbidden shortcuts

Confirm all of the following:

```text
no edit to F1/F3/F4/F6.3 calculators
no edit to production pion subtraction
no edit to w0
no new child normalization
no fitted/tuned correction
no absolute mis-ID probability
no claim NPE=0 is an observed pion sample
no invented truth label
no mandatory zero-truncated/Poisson model
no imported historical F6.1 authority as current authority
no Method B numerical use
no canonical-five expansion
no SIMC amplitude interpretation
no yield/cross-section construction
no production promotion
no farm run
no commit/push
```

No “helpful fallback” to a stale file, nearest filename, other setting, or
historical authority is allowed.

---

## 20. Local validation

Discover `<PYTHON>` through `docs/memory/TOOLS.md`.

Required:

```text
py_compile: PASS
new diagnostic tests: PASS
F6.3 regression: PASS
F4 regression: PASS
F4 analyzer regression: PASS
manifest write/check: PASS
memory health: PASS
git diff --check: PASS
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
```

A nonblocking memory-size warning follows `MAINTENANCE.md`; however this contract
explicitly requires `CURRENT.md` to stay below its 8-KiB soft threshold.

No local test establishes ROOT/PyROOT, full `main.py`, farm rendering, or
scientific validity of the measured quantities.

---

## 21. Farm-validation boundary

**Do not run the farm in this task.**

If the source candidate later passes:

```text
Codex local implementation/checks
-> ChatGPT actual-diff review
-> user commit/push
-> ChatGPT pushed-state synchronization
```

then and only then may a separate farm-readiness review authorize one narrow
`Q4p4W2p74 / Left / lowe` farm diagnostic execution.

The future farm gate must return the diagnostic JSON and PDF with exact source
and input provenance.

Do not provide or execute the farm command now.

---

## 22. Complete actual-diff review bundle

Create temporary repository-root:

```text
kaonlt_review.diff
```

It must include:

- complete tracked diffs for all modified tracked files;
- complete `git diff --no-index /dev/null ...` additions for:
  - the new source module;
  - the new analyzer;
  - the new test;
  - this task contract.

Use the POSIX/Git-Bash procedure in `docs/memory/CODEX.md`.

Do not stage merely for review.

Do not include unrelated local state.

Stop before staging, commit, push, or farm execution.

---

## 23. Acceptance criteria

The local candidate passes only if:

- repository identity/preconditions pass;
- only allowlisted paths change;
- exact current-lineage authority is consumed and verified;
- accepted F.4 reproduces exactly;
- proxy evidence and signed-normalization sensitivity remain separate;
- no direct mis-ID calibration is claimed without an HGC-independent tag;
- no production or correction object is constructed;
- aggregate JSON/PDF generation is deterministic;
- all specified local tests/regressions pass;
- memory/manifest/diff checks pass;
- complete actual-diff bundle is produced;
- no farm run, staging, commit, or push occurred.

If all pass, local status is:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

not runtime validated.

---

## 24. Hard stop

Return `BLOCKED` and stop if:

- branch/HEAD/origin identity differs;
- unexpected worktree state cannot be preserved;
- exact current-lineage F.1/F.3/F.4 authority cannot be established;
- F.4 exact reproduction fails;
- the requested diagnostic would require changing F.3/F.4 or production;
- direct proxy validity would require inventing a truth label;
- event-level identities/factors would have to be persisted;
- Method B would become numerically required;
- a required regression fails;
- memory health has a hard or materially blocking failure;
- CURRENT cannot remain below its required size bound without losing necessary
  active-state information.

Do not invent a workaround.

The purpose of this task is to obtain current-lineage diagnostic evidence,
not to force Method A to pass or fail.
