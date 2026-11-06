# Check 1: full-t acceptance-edge assessment - November 5, 2026

`SOURCE REVIEWED` for detached retained-input arithmetic and plots only.
Scientific conclusion: **evidence remains insufficient**, separately for
Center/low and Center/high. No acceptance-edge cause is established or excluded.
Source observation: local/origin/live HTTPS test 6663874fae232729093a7dd51daffd6c710621a8; scientific
production source remains 9172aba853fd5b5b744c19add992fa4f27466180.

## Measurements and their meaning

Latest October 28 low / October 30 high Local-t evidence and unchanged,
already-normalized iteration-zero KLambda SIMC only. Positive analysis -t is
-MandelT plus the existing kaon-derived t_shift. Shared canonical boundaries
are 0.4, 0.5666666666666667, 0.7333333333333334 and 0.9 GeV^2.
Do not equate shared nominal intervals to shared occupied phase space.

| Setting | Lambda source records t1/t2/t3 | t3 record share | t3 MC yield share |
| --- | --- | --- | --- |
| Left/low | 501/3957/3708 | 45.4% | 52.1% |
| Center/low | 2907/3045/647 | 9.8% | 7.4% |
| Left/high | 122/898/807 | 44.2% | 57.1% |
| Center/high | 1065/1035/336 | 13.8% | 12.4% |
| Right/high | 960/1124/1023 | 32.9% | 22.5% |

Counts include prompt/random/dummy streams and potentially zero proton/RF
weights; they are not net kaon counts, significance or a sparse-bin criterion.
Center broad-MM t3 records are 7779/2624 (low/high), 41.6%/40.2% of their
setting totals: no corresponding broad-population depletion. The Lambda-only
drop is compatible with restricted KLambda support but does not measure its
accepted span or establish why it occurs. Cross sections/backgrounds/weights
also affect population; no overlap length or geometrical fraction is inferred.

Center/low t3 applied B/pre-pion K is 2.65 +/- 0.18% over Lambda and 1.17%
in the whole peak-containing coarse bins [1.102,1.126) GeV. Cleaned high-tail
fraction is 36.47 +/- 2.53% versus MC bounds [24.02,24.80]%.
Center/high t3 B/K is 18.29 +/- 1.47%, peak-containing fraction 6.18%; tail
18.55 +/- 3.47% versus MC [17.91,18.62]%. Tail is [1.132,1.16) GeV.
Errors are conditional counting only, weights fixed; no fitted/shared-systematic
covariance or calibrated significance. B/K is an estimator ratio, not true-pion
probability. Low has a whole-t3 tail difference; high's tail summaries overlap.
Neither identifies edge/interior behavior or validates the background transfer.

## Exact retained-input limitation

The named current pion-t/pion-t-phi JSONs retain MM histograms and source
accounting by canonical t index, not continuous analysis_t or fine t histograms.
The observer row/projection retains ti while discarding continuous t. SIMC-only
NPZs retain canonical t/phi-by-MM yield/error arrays, not fine t or joint t-MM.
Thus full accepted-t distributions, Center's occupied t3 fraction, within-t3
MM/background dependence and matched overlapping-interior tests are not
identifiable from these marginals. No full-t plot or edge explanation is claimed.
This is not a claim that ordinary farm files or implemented functionality are
absent. The precise needed existing evidence is full accepted-t spectra and
joint t-MM/source/selection/applied-weight/Sumw2 information for current data,
the unchanged Local-t application and compatible KLambda MC. No farm operation
or acquisition is authorized here. CURRENT owns the sole follow-up action.

## Deliverable and checks

KaonLT_Check1_Acceptance_Edge_20261105.zip contains a six-page report with five
original figures, embedded HTML/Markdown, all-15-parent occupancy table and
numerical arrays, named original inputs, reproducible scripts and review receipts.
PDF SHA-256: 418e83be531f1e7ad0d4662b3464193b7013ddf553802100032a9f2887786d58.
125 extraction/closure checks; 262 independent direct-parent/serialized-aggregate
and shape-error-Jacobian checks passed. Packaged-input reproduction gives four
byte-identical numerical result files. All six PDF pages rendered and inspected.
These are local detached checks, not new farm/ROOT/production acceptance.

Separate yield-normalized and display-only unit-normalized shape views are
preserved. No new factor, peak match, fit amplitude, cut, binning, subtraction,
model or source change. Native coarse SIMC edges retain their 0.2-MeV offset
from data; no fractional refill. No exclusion threshold/mask; no deleted data.
Method A remains detached, Method B diagnostic/excluded. Outcome B, canonical-
five equivalence, later producer-evidence/final-yield gates and all six roadmap
statuses remain unchanged. Stop after this Check 1 assessment: no Check 2,
Method-A assessment, final acceptance, farm run, commit or push performed.

## Check 1 continuation: existing-output access gate

Live/local/origin test remains 6663874fae232729093a7dd51daffd6c710621a8.
Completed Check 1 arithmetic/figures were reused unchanged, not regenerated.
SOURCE VERIFIED: calculate_yield.py builds continuous adj_t/adj_MM/source/selection/
proton-factor caches (lines 521-618), 100-bin H_t_DATA per canonical cell (1944),
and t support histograms (2914, 4418); it exposes data/sub/SIMC caches (4067-4068,
4582). ave_per_bin.py fills t from adj_t/minus_t and consumes these caches.
main.py writes eligible top-level histograms to OUTPATH/kaon_<OutFilename>.root
(294, 1449-1481) and records analysed data/dummy ROOT outputs (983). This proves
existing producer capability, not that every nested cache/histogram is persisted
or that an ordinary saved histogram has selected-Local-t lineage.

Supplied October 28/30 run summaries identify ordinary output and SIMC roles;
their personal paths remain in the external gate receipt. The local analysis
checkout has no OUTPUT or simc_gfortran link. Inspected local simulation output
locations and Downloads have no ROOT event files; the current Work projections
are already documented. A direct retained-JSON continuous-coordinate field
search found only a historical Center F.1 contract, with no established current
selected-Local-t data/MC pair; it is not substituted. No broad archive study,
new extraction, farm login or execution was initiated.

First precise missing requirement: a locally readable current-cohort numerical
source containing continuous analysis_t/MM and source/selection/unchanged-weight
information, starting with Center/low and Center/high and retaining all-five
scope. The owners above already exist. This is an access/availability limit in
this session, not evidence that ordinary farm outputs were deleted or never
produced. Their current remote existence is NOT VERIFIED. Single next action:
make the existing matching outputs or their existing numerical exports locally
readable, then resume the full-t/common-interior test. No rerun, new producer,
normalization investigation or source change follows from this gate.

User-provided context, record only: K+Sigma0 and K+Lambda have different tmin;
the user notes a possible combined larger-t effect of reduced Lambda population
and relatively greater Sigma0 contribution. This is not independently tested,
not an additional Check 1 hypothesis, not an attribution of the discrepancies,
and not authority for a Sigma0 study or changed pion decomposition.

Both Center conclusions remain insufficient. Original numerical results/plots,
six-step roadmap, accepted scopes, Outcome B and Method A/B boundaries remain.
KaonLT_Environment.md was explicitly waived by the user as ChatGPT-specific.
Stop at this input access requirement; no Check 2, Method A or farm operation.

## Check 1 farm-evidence acquisition: pre-Codex contract

User subsequently authorized the bounded source audit and contract preparation,
not implementation or execution. At test 6663874fae232729093a7dd51daffd6c710621a8,
Path B is established for the existing controlled acquisition: the detached
Local-t observer omits continuous t from its serialized MM projections, SIMC's
average cache has t but no MM and is allcuts-only, and normal xsect support has
no t-MM. The retained simulation t-vs-tmin arrays do support a fine-t projection,
but not t-MM or its lost 2D Sumw2. Ordinary ROOT contents/remote existence remain
unverified; no global absence claim follows.

The [single proposed contract](../phases/check1-continuous-t-mm-acquisition-2026-11-05-task-contract.md)
specifies default-off scalar exports consuming the actual selected Local-t
application and the unchanged SIMC loop, with no new fit, cut or normalization.
Its existing-owner extension must verify and collect both epsilons: 15 original
artifacts plus 20 new sidecars, with five-setting provenance/closure, preserved
conditional Sumw2 and explicit unavailable fitted covariance. A pinned external
reference contains October 28/30 experimental evidence and unchanged simulation-
only arrays; no October 9 experimental NPZ fields are included.

Gate: independent contract review; no production source changed. Later gates
are scoped Codex implementation/checks, independent actual-diff review, user
publication/pushed review, then reviewed user preflight/acquisition. Existing
source suites NOT RUN by Work for contract preparation; usable Windows Python/
NumPy observed, default WSL Python lacks NumPy. Manifest/ordinary health checked;
soft size warnings remain nonblocking and scheduled for the next checkpoint.
Both Center physics conclusions, six roadmap statuses and the user-provided
untested Lambda/Sigma0 observation remain unchanged. No Check 2, Method A,
farm execution, commit or push occurred.

## Codex implementation startup: required-test hard stop

`BLOCKED`, November 5: independent contract review passed per the user's
implementation request. SOURCE VERIFIED: live branch test and HEAD
6663874fae232729093a7dd51daffd6c710621a8 match the contract; index empty,
with pre-existing pending memory/offline-testing/review files. No implementation
source was changed and the three proposed new implementation files are absent.

The bounded existing-environment inventory found Ubuntu system Python 3.12.3
and one existing virtual environment, also Python 3.12.3; both fail NumPy import.
No other candidate was found in the inspected Linux home/opt/usr-local roots.
Windows Python 3.12.10 / NumPy 2.3.0 imports successfully, but does not establish
POSIX-sensitive owner/alias/permission regression coverage.

Actual Linux command, using the existing virtual-environment interpreter:
`python -B -m unittest testing.test_pion_t_phi_template_comparison`.
Result: one unittest loader error, zero substantive test cases executed;
`test_pion_component_dynamic_alignment` imports `pion_component_fits`, which
raises `ModuleNotFoundError: No module named 'numpy'`. The system interpreter
separately fails the direct NumPy import probe. No dependency was installed,
reconfigured or faked. This is the contract's explicit "required test cannot be
run" hard stop, not evidence of an implementation or scientific defect.

Unexecuted: new packet/round-trip/actual-source projection/default-off/SIMC-seam/
five-setting/negative/failure suites (implementation not begun); remaining
parent/owner/profile/collector regressions and changed-Python grammar checks.
ROOT/PyROOT, main.py and farm execution remain unexecuted and NOT VERIFIED.
Resume at the same scoped implementation/local-test gate when an existing
suitable Linux environment is available, then independent actual-diff review;
user publication/pushed review and reviewed user preflight/acquisition remain
conditional later gates. Prior Check 1 science and all frozen scopes remain.

## Resumed implementation under the reviewed testing correction

Interim checkpoint, superseded by the completion record below: `ACTIVE`, partial
implementation. The user explicitly approved supported Windows
Python 3.12.10 / NumPy 2.3.0 deterministic execution and corrected the contract:
unexecuted Linux-specific tests are NOT VERIFIED readiness conditions, not an
implementation hard stop. Failed numerical/integrity/preservation/executable
tests remain hard stops. No package installation or environment repair follows.
The prior blocked-attempt record above is preserved as history, not current NEXT.

Local additions: ROOT-independent packet/array validation, signed/conditional
Sumw2 and flows, existing-cache identity join and existing Local-t projection
closure, atomic NPZ/JSON writer/readback, and the combined schema-v4 profile
(15 original artifacts plus 20 sidecars). calculate_yield has an optional sink
at the actual SIMC loop before its allcuts continue; existing callers still pass
no sink. No observer/export invocation, process_hist_simc construction/finalization,
Check 1 owner mode, reference admission or final publication path is yet wired.
The source-through-collection/cleanup route is therefore NOT VERIFIED/incomplete.
Do not label these changes DEVELOPMENT COMPLETE or authorize acquisition.

Windows new suite: 14 tests passed, covering real project/apply_transfer and
serializer closure, exact cache joins, signed cancellation/quadratic Sumw2,
normalization-once, actual SIMC loop disabled-path fills/cache parity, native
histogram closure, NPZ roundtrip/object rejection, collisions and injected writer
failures. An initial NameError caused by a misplaced test block was diagnosed
and corrected within that test file; no assertion was removed or weakened.
Existing Windows parent/child suite: 57 tests, two symlink-capability skips.
Existing Windows collector suite: 28 tests passed. Existing Linux owner suite
does not need NumPy: 95 tests passed with one unavailable-Poppler skip, using the
existing Python 3.12.3 virtual environment. These component results do not prove
full-five owner/sidecar/reference admission, final-ZIP release or cleanup safety
for a Check 1 mode that has not yet been implemented.

Remaining input requirement: the actual external
KaonLT_Check1_Reference_Inputs_20261105.zip at the contract's pinned SHA-256.
It was not found in the checked session Downloads location; an accessible path
was requested while independent component implementation/tests continued.
Its actual member schema and original/current-runtime identities must be read
before implementing exact cohort/reference/MC admission; do not invent fields
or replace it with another cohort. This is an input-access limitation, not the
superseded WSL NumPy blocker or a global claim of evidence absence.

NOT VERIFIED: NumPy-dependent Linux parent/child numerical suites (requires an
existing Linux NumPy environment), the two Windows symlink cases (requires
symlink capability), Linux owner's real PDF-page regression (requires Poppler
pdfinfo), and all farm Python/ROOT/dependency/caller/path/provenance conditions.
The full-five owner/collector negative/failure tests remain unexecuted pending
completion of that owner path. Mandatory reviewed farm preflight remains before
any scientific acquisition. No farm/ROOT/main execution, commit or push occurred.

## Check 1 implementation completed locally; actual-diff review pending

`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`. Live branch test/HEAD remain
6663874fae232729093a7dd51daffd6c710621a8; changes are an uncommitted local candidate.
The external reference subsequently became accessible and was verified read-only
at SHA-256 02ff71998fd51d5e56f07484bb1c71185f78ed51d6ce409b0d34154bda7a6d7a.
Its seven-member manifest, byte hashes, nested original bundles, original run
summaries and simulation-only arrays were checked. All 135 reference Local-t
applications are available. No reference content was copied into the repository.

SOURCE VERIFIED route: default-off observer call after its unchanged successful
serialization -> exact cache/entry/tree/flag join and actual Local-t vectors ->
validated row/projection NPZ plus metadata-last JSON; same existing SIMC loop
captures both existing flags before its allcuts continue, with the same hole
decision -> finalize against actual once-scaled histograms. Original caches,
fit/transfer/PID/proton machinery, cuts, normalizations and algorithms are unchanged.
Rows keep signed/zero contributions, raw and nonzero-support counts, conditional
quadratic Sumw2, flows, parent/child/outside-phi identities and unavailable fit
states. Fitted/training/shared-systematic covariance is not manufactured.

The existing execute_gate now owns the explicit Check 1 mode: strict descendant
allowlist plus memory prefix and unchanged historical blob pins -> pinned
reference and original source/cohort/template admission -> existing private
code/data/ltsep isolation and paired low/high processing -> existing parent/child
checks for both epsilons plus row/native-MM/fine-t/reference checks -> unchanged
schema-v4 collector with exactly 35 artifacts -> pending ZIP readback -> private
cleanup and all original/input/installed-package/external-link/35 ordinary-target
preservation checks -> hash-bound companions, completed staging cleanup, then
atomic no-overwrite final-ZIP rename as the last release operation.
Missing verification/cleanup/preservation receipts or any collision/failure
prevent final release. First-invariant status survives private cleanup. Old
owner modes and both old profiles retain their default semantics.

Coordinate fingerprints include diagnostic plot paths in the frozen source.
Reference and current coordinate self-identities are validated independently;
both fingerprints are retained with their distinct roles. A new regression
exposed the parent builder's kaon_data_coordinate_fingerprint alias. The repair
validates all three source-declared coordinate identity aliases against the
complete self-identity and retains their fields under the full path-independent
transform/provenance identity; it does not drop the failing alias. Corrupt aliases
and changed MM transforms must still fail. The closed comparison checks all
remaining scientific values/provenance, excludes explicitly declared execution/
lineage fields and compares artifact/root filenames by basename.
No transform, template contents, transfer, threshold or numerical tolerance is
changed; comparison remains rel_tol=1e-8, abs_tol=1e-12, with exact typed counts
and identities. Old source/runtime-container hashes are never interchanged.

Local checks: final consolidated Windows Python 3.12.10 / NumPy 2.3.0 run:
116 tests, no failures, four explicit skips (two existing symlink-capability
cases, new writer symlink capability, real-Poppler PDF-inclusive owner). This
includes the actual read-only pinned external reference and 31 new
source/observer/numerical/export/collection cases alongside 85 existing cases. An inherited mocked PDF page count was
removed from the new fixture: actual JSON/source/numerical validators and the
unchanged collector run directly; the entire PDF-inclusive owner requires the
real existing Poppler executable and is NOT VERIFIED here. No unavailable page
gate is reported passed. The source-traced alias regression was repaired and
its MM-shift negative requires the precise earlier physical-identity failure.
Existing Ubuntu virtual-environment Python 3.12.3 owner regressions: 100 tests,
one missing-pdfinfo skip, no failures. Changed/new Python compiles without ROOT
and passes Python 3.9 grammar checks; this is not execution on farm Python 3.9.
Direct kernel no-overwrite rename/collision checks passed on Windows and
Linux; the actual farm output filesystems remain unverified until preflight.
Failure tests include receipt and post-link cleanup errors, rollback and
preservation failures, nonfinite values and even tiny negative moment values.
The final-ZIP rename occurs only after staging cleanup, with no fallible
post-publication unlink/cleanup step.
Complete actual diff includes all additions; nothing was staged, committed or
pushed. CURRENT alone owns the independent actual-diff-review NEXT.

The executable inventory explicitly adds the exporter, observer, producer,
owner, profile, immutable snapshot helper and prompt-tree resolver. Original
fit/covariance evidence is linked by the actual serialized child JSON hash, not
just a basename. Supplied summaries show all ten pion containers unchanged
from pre-copy originals; those hashes are reused and independently checked
against the actual post-child containers. Kaon post-derived observer hashes and
new SIMC consumed-file hashes retain their distinct original/runtime roles.

NOT VERIFIED: NumPy-dependent Linux execution of the new numerical/export/five-
setting suites; actual symlink packet/collision cases on Windows; the real
PDF-inclusive five-setting owner and Linux six-page pdfinfo case; farm
interpreter/ROOT/root_numpy/dependency identities,
actual caller/namespace/path/alias/overlay behavior and real producer integration.
These are readiness conditions, not fabricated passes. The new reviewed farm
preflight must execute dependency/interpreter/ROOT and no-clobber capability
probes, complete source/profile/reference/cohort/manifest/caller/I/O checks and
preservation, without scientific copies, analysis or ZIP. Only reviewed successful
preflight evidence at the independently accepted/pushed SHA can precede a
user-authorized acquisition. No farm readiness, physics acceptance, Check 2,
Method-A adoption, final yields, ROOT/main/farm run, commit or push is claimed.

## Independent actual-diff review: two narrow repairs

`ACTIVE`: local Check 1 candidate requires a narrow repair; source-review gate
is `BLOCKED`. At unchanged test 6663874fae232729093a7dd51daffd6c710621a8,
the supplied 224585-byte complete 12-file implementation diff matched SHA-256
ade4d945c198ec6968e7f3c907f55f381aab828354c83d861649778ab0f6924d and passed
reverse apply --check against the actual files without applying it.

SOURCE VERIFIED P1: main.py creates Setting_kaon_Simc_Q4p4W2p74_epsilon.root,
compare_simc opens it, and calculate_yield_simc passes it to export_simc. The
new exporter instead requires Prod_Coin_Q4p4W2p74setting_epsilon.root and raises
check1:simulation_runtime_role. The owner also compares runtime bindings only
with the archived-copy map. Preserve archived-input admission and explicitly
bind/hash the consumed derived model as a distinct runtime role through exporter,
owner and same-attempt validation. Accepting another basename alone is incomplete.
The new fixture uses the archive name and misses the live producer seam.

SOURCE VERIFIED P2: write_packet's rollback does not surround temporary-context
exit. An injected cleanup exception after both links left the complete JSON/NPZ
pair and read_packet accepted it as captured_pending_owner_verification. No final
scientific ZIP release is inferred. Enclose context exit in pair rollback,
preserve the first failure and collision safety, and add the cleanup regression.

Independent reviewer checks: Windows Python 3.12.10/NumPy 2.3.0 numerical/export/
parent-child/collector suites ran 111 tests, no failures, 3 skips (one real-Poppler
gate and two symlink-capability cases), with the real pinned reference read-only.
Existing Ubuntu environment owner suite ran 100 tests, no failures, 1 pdfinfo
skip. The two focused local probes reproduced the findings. Linux NumPy/export
parity, real PDF-inclusive checks and all farm/ROOT conditions remain NOT VERIFIED.
Manifest/ordinary health pass with existing nonblocking soft warnings; batch at
the next checkpoint. No source patch, staging, commit, push or farm run by Work.

Single NEXT: repair these two defects and their regressions under the existing
contract/allowlist, then independent actual-diff review. Do not commission another
support architecture/normalization investigation or issue a farm command. Both
Center physics conclusions, accepted closures, Outcome B and Method A/B roles
remain unchanged. Detailed review and candidate/source snapshots are in the
external KaonLT_Check1_Independent_Review_20261105.zip review package.

## Local repair of the two independent review findings

`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`: local repairs to exactly P1
and P2 under the same contract/allowlist; independent actual-diff acceptance
remains pending. The review findings above remain the record of the failed
candidate. No new architecture, contract, environment repair or farm task.

SOURCE VERIFIED P1: original Prod_Coin ROOT/.hist admission is unchanged.
The owner authority declares each exact archive-to-Setting_kaon_Simc model
mapping and existing copy/iteration-zero operation. Export hashes the actual
consumed model ROOT and copied model .hist; metadata retains distinct archive
original, staged archive, consumed model and histogram roles. The owner hashes
all five model ROOTs and five model .hist files independently after the child,
checks unchanged archives against the pinned original summaries and admitted
staging identities, and checks packet bindings against those actual runtime
hashes. Archive/model hashes are never equated; model .hist must equal its
unchanged archived .hist. Missing or substituted roles fail. Cost/status/summary
records include the derived-file scope and hashes. Main, iter_weight and
compare_simc are now explicit immutable executable provenance dependencies;
none of those files, the event weights or normalization was modified.

The five-setting fixture now uses the real producer runtime basenames. It
executes the actual source-extracted _prepare_iter_weighted_kaon_simc copy/name/
iteration-call function with controlled bytes and a callback that marks only
the copied model, proving the original and model hashes differ. This is a
ROOT-independent wiring test, not execution of the scientific ROOT operation.
Positive full-five numerical/collector closure and strict negative tests cover
missing archive/model/.hist roles, changed hashes, archive-for-model substitution,
changed derivation authority, actual same-basename model replacement, archive
basename passed to export, and model .hist substitution. Original archive/.hist
reference admission remains mandatory and is tested independently.

SOURCE VERIFIED P2: the rollback boundary now surrounds TemporaryDirectory
context exit, not just writes/links. Rollback checks device/inode ownership and
never removes an existing competitor replacement. Original exceptions survive
subsequent cleanup or rollback errors; secondary rollback errors are attached
to that exception. The exact regression patches TemporaryDirectory.cleanup
after both published sidecars have been read successfully, raises the original
cleanup exception, and requires both attempt-owned sidecars removed. Adjacent
cases preserve unrelated bytes/replacement files and retain a prior link failure
if context cleanup also fails. Hard-link no-clobber semantics remain unchanged.

Executed in order: five focused P1/P2 tests, then the focused archive/.hist
admission test: six passes, no skips. Final Windows Python 3.12.10 / NumPy 2.3.0
consolidated suite: 122 tests, no failures, four skips (two existing symlink
capability cases, new sidecar symlink capability, real-Poppler PDF-inclusive
five-setting owner). The pinned external reference was read-only. Existing
Ubuntu virtual-environment Python 3.12.3 owner suite: 100 tests, no failures,
one real-pdfinfo skip. All seven changed/new implementation Python files compile
and pass Python 3.9 grammar checks; this does not execute farm Python 3.9.
Existing numerical, default-off, source/provenance, five-setting, collector,
negative/failure and original owner regressions remain enforced.

NOT VERIFIED: Linux NumPy-dependent numerical/export/five-setting parity (the
existing Linux interpreters lack NumPy); real PDF-inclusive owner and six-page
pdfinfo checks (requires Poppler); skipped Windows symlink cases (requires
symlink capability); actual farm Python/ROOT/root_numpy/dependencies, namespace,
caller/path/alias, source/profile/reference/cohort/provenance/manifest and real
producer integration. Reviewed farm preflight remains mandatory before any
scientific acquisition is authorized. No packages installed or assertions
weakened; no ROOT/main/farm execution or readiness approval.

HEAD remains test 6663874fae232729093a7dd51daffd6c710621a8. Pending files and
reviewer memory changes are preserved. Only the exporter, existing owner and
Check 1 test file changed for these repairs, plus CURRENT/this record/manifest.
The existing contract, other implementation files and prior review bundles
remain unchanged. NEXT is independent review of the complete updated actual
diff; no staging, commit, push or farm execution was performed.

## Independent P1/P2 repair review: source acceptance

`SOURCE REVIEWED`: independent actual-diff review passes the repaired local
candidate at unchanged test 6663874fae232729093a7dd51daffd6c710621a8. Submitted
complete 12-file diff SHA-256:
4d9133db1f514467dfd9ac841945980cad266afe16d57bf8461f799be32d7e9f.
Hash and reverse apply --check matched the actual files before review-memory
updates. Only exporter/owner/Check 1 tests changed in source since the failed
candidate; other producer hooks, profile, tests and task contract are unchanged.

P1 is resolved through distinct admitted archive ROOT/.hist and consumed model
ROOT/.hist roles, exact existing derivation authority and main/iter_weight/
compare_simc executable identities. The owner independently hashes the actual
derived runtime files and compares the packets against those identities before
cleanup. Archive reference/cohort admission remains strict. The real producer
copy/name/iteration-call seam is covered with controlled scalar bytes; that
test does not execute ROOT or establish scientific model runtime validation.

P2 is resolved: rollback surrounds temporary-context exit, preserves the first
exception, tracks attempt-owned paths by device/inode and leaves replacement
or unrelated files intact. The actual post-publication cleanup regression and
adjacent double-failure cases passed. Existing no-clobber behavior is preserved.

Reviewer executed the six focused new regressions first: 6 passes, no skips.
Windows Python 3.12.10/NumPy 2.3.0 then ran the named numerical/export/parent-
child/collector subset: 117 tests, no failures, 3 skips, with the pinned real
reference read-only. This is distinct from Codex's reported 122-test run.
Existing Ubuntu Python 3.12.3 owner suite: 100 tests, no failures, 1 pdfinfo
skip. Seven-file compilation and Python 3.9 grammar pass without bytecode;
grammar is not execution on farm Python 3.9. Manifest and ordinary health pass;
existing CURRENT/MEMORY soft warnings remain nonblocking, batched at checkpoint.

Linux NumPy/export parity, skipped Windows capabilities, real PDF-inclusive/page
checks and actual farm Python/ROOT/dependency/caller/path/alias/source/profile/
reference/cohort/provenance/manifest conditions remain NOT VERIFIED. Source
acceptance does not authorize acquisition or clear the scientific Check 1 gate.
No implementation source was changed by the reviewer; unrelated pending bytes,
previous review bundles, raw index and refs remain preserved. No staging,
commit, push, farm or ROOT/main execution by Work.

Single NEXT: user publication, then independent pushed-state synchronization,
then reviewed preflight-only readiness evidence before a separately authorized
user acquisition. No new contract, normalization study or farm command follows
from this acceptance. Both Center conclusions, Outcome B, Method A/B roles and
all accepted runtime closures remain unchanged. Review materials are preserved
in the external KaonLT_Check1_P1P2_Independent_Review_20261105.zip package.

## November 6 requested contract revalidation

The repeated implementation request was checked against actual test HEAD
6663874fae232729093a7dd51daffd6c710621a8, the full dirty worktree and empty
index. All eight implementation/test/profile diffs exactly match the source-
accepted P1/P2 candidate. The older full review diff no longer reverse-checks
against the three subsequently updated review-memory records; its source-only
reverse check passes, and exact per-file implementation diff comparisons pass.
No additional source change is warranted or made. Preserve the independent
SOURCE REVIEWED acceptance above and CURRENT's conditional publication/pushed-
state/preflight sequence; this receipt neither reopens source nor authorizes it.

Fresh Windows Python 3.12.10 / NumPy 2.3.0 run: 122 tests, no failures, four
skips (real PDF-inclusive five-setting owner requires Poppler; sidecar dangling
symlink and two existing publication collision subcases require Windows symlink
capability). It used the pinned real reference ZIP read-only. Existing-source
project/apply_transfer/serializer closure, exact identities, signed contributions/
Sumw2, normalization-once, default-off, five-setting numerical/profile/collector,
P1/P2 and strict negative/failure regressions passed where executable. Ubuntu
virtual-environment Python 3.12.3: 100 owner tests, no failures, one unavailable
real-pdfinfo six-page gate. Seven Python files compile and pass Python 3.9 grammar.

NOT VERIFIED: Linux NumPy-dependent numerical/export parity (requires an
available NumPy-enabled Linux interpreter), real PDF-inclusive/page gates
(requires Poppler), the skipped Windows symlink-capability cases, farm Python/
ROOT/root_numpy/dependency/caller/path/alias/provenance/manifest conditions and
real producer integration. No environment repair, input substitution, package
installation, ROOT/main/farm execution, readiness claim or Git publication.
Only this allowlisted validation receipt and its manifest change in this pass;
all implementation, CURRENT/reviewer state, unrelated pending work and prior
review bundles are preserved. Return the fresh complete actual diff for the
requested independent review and stop before staging/commit/push/farm work.
