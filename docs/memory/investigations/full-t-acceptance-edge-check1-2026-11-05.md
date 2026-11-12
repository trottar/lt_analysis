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

## Lustre publication-capability failure and local repair

`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`. User-supplied farm observation
at test 374bed263e077c356ea06ae7eea6416886b61d4a: check1_runtime_preflight
failed with OSError errno 22 (EINVAL) from renameat2(RENAME_NOREPLACE) on the
private Lustre temporary renamed target. Preceding hard-link/no-clobber probes
succeeded. Source/reference admission, environment setup and preservation had
no reported failures; analysis never started. This is attributed supplied runtime
information, not an independently executed farm check here. The local rename
pass did not establish Lustre support. The existing contract was amended with
the exact failure and bounded invariants before any source edit; no new contract.

SOURCE VERIFIED route: both private-output and delivery-parent capability probes
now execute the same selected hard-link no-overwrite publication helper used by
release, verify complete target bytes and pending-source removal, and require an
existing-target collision to fail. Unsupported renameat2 is removed entirely;
EINVAL is not caught or ignored and no overwriting rename is substituted.
The final ZIP is linked atomically from an already verified same-directory
pending ZIP. os.link refuses existing targets, including concurrent creation and
symlinks. Device/inode checks protect source and target identities. Source cleanup
failure rolls back only this attempt's link; its identity remains in the outer
release ledger even if the first rollback fails.

Release retains rollback through pending unlink, status/receipt publication and
TemporaryDirectory context exit. Success requires completing that entire path.
On failure, publication_complete is false, original status becomes failed, an
owned success receipt is invalidated before removal, and rollback deletes only
owned regular inodes. Independently replaced files/symlinks are preserved. If
rollback deletion itself fails, residual files remain failed/inadmissible evidence
and secondary errors are recorded; neither their existence nor an old success
receipt establishes acceptance. First errors survive later context cleanup,
receipt-write and rollback errors. The actual owner's failure-finally path also
retains the original error, protects a replaced pending file, and executes
preservation even when pending cleanup fails. Ordinary modes retain their old
failure behavior. No scientific/helper P1/P2 algorithms or collector changed.

Focused tests executed first: eight Linux tests passed; eight Windows tests ran,
no failures, one unavailable-symlink skip. The exact old committed Linux function
was executed against an injected renameat2 errno-22 response after real successful
hard-link/no-clobber checks. This reproduces the supplied syscall failure class,
not the unavailable Lustre filesystem. Positive local filesystem operations,
file/concurrent/symlink collisions, source replacement, post-final-link source and
context cleanup failures, rollback errors, receipt replacement/invalidation,
first-error/status correctness and both role probes passed where executable.
The actual owner finally AST is replayed with real files/status writes to prove
that source-cleanup failure still runs preservation and keeps the first error.
The dependency subprocess is mocked only in capability-unit tests; no installed
ROOT/dependency pass is claimed from those tests.

Full Windows Python 3.12.10 / NumPy 2.3.0 regression: 129 tests, no failures,
five skips (real PDF-inclusive owner; existing sidecar symlink; two existing
publication symlink subcases; new publication symlink test). It includes the
pinned external reference read-only, existing-source projection/transfer/serializer
closure, signed weights/Sumw2, exact identities, normalization-once, P1/P2, default-
off, provenance, all-five numerical/profile/collector and negative/failure tests.
Existing Ubuntu virtual-environment Python 3.12.3 owner suite: 107 tests, no
failures, one real-pdfinfo skip; actual regular/dangling symlink publication
collisions passed there. Seven Python files compile and pass Python 3.9 grammar;
that is not farm Python 3.9 execution. Source scope is the owner and owner tests;
memory scope is CURRENT, existing contract, this investigation and manifest.

NOT VERIFIED: actual private Lustre and delivery-parent behavior; Linux NumPy
numerical/export parity; real Poppler PDF-inclusive/six-page checks; unavailable
Windows symlink cases; full farm Python/ROOT/root_numpy/dependency/caller/path/
alias/source/profile/reference/cohort/provenance/manifest readiness and producer
integration. No package installation, environment repair, assertion weakening,
scientific change, farm execution, staging, commit, push or readiness approval.
Preserve all pending unrelated bytes/index/refs and prior review bundles. After
independent actual-diff acceptance, user publication and pushed-state review,
the next farm gate remains the SAME preflight-only operation. Full acquisition
is forbidden until its successful readiness evidence is reviewed and the user
separately authorizes acquisition. Prior P1/P2 review and science scopes stand.

## Independent Lustre repair review: pending-file ownership gap

`ACTIVE`; independent source-review gate `BLOCKED`, narrow repair required.
Reviewed six-file actual diff at local/origin test
374bed263e077c356ea06ae7eea6416886b61d4a: 46302 bytes, SHA-256
081b36f4265e89e9a58c1212be2690cac0349332599a9f40046fea4f339e4329.
Hash, inventory and read-only reverse check matched actual files. The hard-link
helper, dual filesystem probes, receipt invalidation, first-failure retention
and context-exit rollback were inspected; no scientific/P1/P2 algorithm changed.

SOURCE VERIFIED P2-priority finding: execute_gate captures pending_identity only
immediately before release (line 1968), after collection/ZIP verification/private
cleanup/preservation/summary writing. Its finally branch (1992-1994) blindly
unlinks a regular non-symlink pending path when that identity is None. A real-file
replay of the actual finally AST under a simulated summary failure, with an
independently replaced pending file and no ownership ID, deleted that replacement.
The first exception was preserved and preservation ran; ownership protection
failed. Existing new finally coverage always supplies a known identity.

Single NEXT: fix deletion ownership across the entire pending-file lifecycle;
without trustworthy ownership, retain the visibly pending/failed file rather
than delete it. Add failure-before-identity-capture/replacement coverage while
preserving known-owned cleanup, first-error and preservation checks. Use the
same contract/allowlist; no new architecture or normalization investigation.

Reviewer Windows Python 3.12.10/NumPy 2.3.0 full numerical/publication cohort:
129 tests, no failures, 4 skips, pinned reference read-only. WSL Python 3.12.3
full owner suite: 107 tests, no failures, 1 pdfinfo skip. Initial focused command
included a misnamed extra method (reviewer invocation error); seven publication
tests passed there and the correct full owner suite subsequently passed. The
additional actual-finally reproduction exposes a branch beyond those passes.
Compilation/3.9 grammar is implementer-reported, not newly rerun in this review.
Manifest/ordinary health pass with existing nonblocking soft-size warnings;
batch at checkpoint. Actual Lustre, Linux NumPy/export parity, real PDF/page,
skipped capabilities and farm environment/provenance conditions remain NOT VERIFIED.

No source edit, staging, commit, push or farm operation by Work. Unrelated bytes,
previous bundles, raw index and refs remain preserved. After later acceptance,
user publication and pushed-state review, the next farm action remains the SAME
preflight-only gate. No acquisition authority follows. Prior P1/P2 source review,
scientific Check 1 insufficiency, Outcome B, Method A/B and accepted closures stand.
External review package: KaonLT_Check1_Lustre_Independent_Review_20261106.zip.

## Pending ZIP ownership correction after review

`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`, independent actual-diff review
pending. The reproduced finding above remains the history of the failed local
candidate; existing contract scope and all accepted scientific/P1/P2 results stand.
A regression replayed the actual execute_gate finally AST before source repair:
both unknown-identity subcases deleted the pending file, including an independent
replacement. Known-identity cases were already protected. No scientific failure
or farm execution was inferred from that expected pre-fix regression failure.

SOURCE VERIFIED correction: unknown ownership no longer has any path-inspection/
unlink fallback. It retains visibly pending/failed evidence and still executes
preservation and propagates the original exception. Normal successful ZIP
verification now opens the pending file, validates its actual CRC/inventory/
manifest/source/member hashes with the existing validator, hashes that same open
file, and requires stable device/inode/size/mtime plus the same regular inode at
the pending path. Only successful validation returns the trusted cleanup identity.
Capture occurs before private/source cleanup, preservation and summary writing.
Validation failure or replacement during verification leaves identity unknown,
so failure cleanup retains the file. Pre-release identity comparison rejects a
replacement rather than adopting the newly inspected path. Existing known-owned
cleanup, first-error retention, hard-link publication, rollback and ordinary modes
remain unchanged. No collector/schema or contract change was needed.

The two added regressions cover unknown/known ownership with original/replaced
pending files through the actual finally block, retained first exceptions and
preservation, valid versus wrong-source ZIPs through the real 35-member validator,
and replacement of an open validated ZIP with identical bytes. Linux executes
that replacement and rejects the inode mismatch. Windows refuses unlink of the
open file with WinError 32; the test explicitly requires that OS behavior and
preservation of the original inode. An initial Linux-style Windows fixture exposed
that platform difference and was corrected without allowing another exception or
weakening the Linux mismatch assertion. No dependency or environment repair.

Focused tests then passed: Linux ten, no skips; Windows ten run, one unavailable-
symlink skip. Full Windows Python 3.12.10 / NumPy 2.3.0 cohort: 131 tests, no
failures, five skips (real-Poppler full owner; three existing symlink-capability
cases; publication symlink case). It includes the real pinned reference read-only,
P1/P2, existing numerical/projection/default-off/provenance/five-setting/collector
and failure regressions. Existing Ubuntu virtual-environment Python 3.12.3 owner
suite: 109 tests, no failures, one real-pdfinfo skip. Both changed Python files
compile and pass Python 3.9 grammar checks. These do not execute farm Python 3.9.

NOT VERIFIED: actual Lustre/delivery-parent and installed farm runtime/preflight;
Linux NumPy-dependent numerical/export parity; real Poppler PDF/page checks;
skipped Windows symlink capabilities. Prior test results are not farm readiness.
No scientific, P1/P2 or existing hard-link helper/release changes in this pass;
source scope is the owner ownership lifecycle and two added owner regressions.
Only CURRENT/this record/manifest are updated in memory; the existing contract,
reviewer findings and unrelated pending work remain preserved. Return the complete
revised actual diff, then stop before staging, commit, push or farm execution.
After later acceptance, user publication and pushed-state review, the next farm
operation remains the SAME preflight-only gate, never automatic acquisition.

## Independent pending ownership review: publication identity handoff

`ACTIVE`; source-review gate `BLOCKED`, one narrow repair remains. Reviewed
six-file diff at test 374bed263e077c356ea06ae7eea6416886b61d4a: 61912 bytes,
SHA-256 aada0840aafc661658df8758595d828916d974be8a85dae11ad1f13bbd95efcd.
Hash, inventory, read-only reverse application and whitespace checks pass.

The prior unknown-identity deletion is corrected: finally retains unknown-owned
files, stable open-ZIP validation returns trusted device/inode identity before
private cleanup, and a pre-release comparison rejects replacements visible there.
SOURCE VERIFIED remaining P2-priority issue: release_check1 receives no trusted
identity and _check1_publish_no_replace independently captures source.stat as
authority. An identical-byte replacement after the last guard but before that
stat is adopted. A real temporary-file release probe injected such replacement
at the helper seam: release status success, published inode differs from trusted,
replacement pending path removed. This is an ownership invariant failure, not
evidence of different scientific ZIP content. Original inode was retained under
a hard link, so the reproduction does not depend on inode reuse.

Single NEXT: pass the validated identity through release/publication and check
against it rather than adopting a new inode. Add a same-content/different-inode
replacement test between the last guard/hash and linking. Reject without deleting
or adopting that replacement. Preserve the corrected unknown retention and all
existing rollback/collision/first-error/preservation checks. Existing contract
and allowlist suffice; no new architecture or scientific source change.

Reviewer Linux full owner suite: 109 tests, no failures, one pdfinfo skip,
41.680 seconds. Windows full rerun did not complete: initial automatic approval
review timed out and the allowed retry was interrupted. The implementer's
131-test Windows and compilation/3.9 grammar results remain attributed, not
independently rerun here. Additional real-file release reproduction above exposes
a gap beyond the passing Linux suite. Manifest/ordinary health pass with existing
nonblocking soft warnings, scheduled for checkpoint. Actual Lustre/farm, Linux
NumPy/export parity, real PDF/page and unavailable capabilities remain NOT VERIFIED.

No implementation source, stage, commit, push or farm operation by Work. Previous
bundles, unrelated pending bytes, raw index and refs remain preserved. Prior
P1/P2 acceptance and scientific scopes remain unchanged. After later source
acceptance, user publication and pushed-state review, the next farm gate remains
the SAME preflight-only operation. No acquisition authority follows. External
package: KaonLT_Check1_Pending_Identity_Review_20261106.zip.

## Validated identity handoff correction

`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`: local correction only; independent
actual-diff review remains next. Same test HEAD
374bed263e077c356ea06ae7eea6416886b61d4a and existing contract. The review finding
above is preserved. A newly added regression executed the actual execute_gate
final guard and injected same-byte replacement at the publication-helper seam
(after the guard and release hash). Before repair it failed with ValueError not
raised: the replaced inode had been adopted. Its original inode was retained by
a hard link, so the result does not rely on inode-number reuse.

SOURCE VERIFIED: the validated tuple now travels from execute_gate through
release_check1 into _check1_publish_no_replace, and missing/invalid identity is
rejected. The helper no longer captures a source.stat tuple as authority. Its
lstat and pre/post-link checks compare against the supplied identity, and the
rollback ledger retains that identity. The capability probe passes an identity
captured at its own test-file creation. Existing release fixtures likewise pass
fixture-owned identities recorded before injecting failures; no assertion was
removed or weakened. The existing sidecar exporter/P1/P2 code, normal owner modes,
collector schema, scientific inputs/weights, normalization and contract bytes
are unchanged in this correction.

The real-file post-guard regression now requires failure, absent final ZIP,
replacement inode/bytes retained, failed source gate receipt, publication_complete
false and removal of attempt-owned companions. Additional tests reject missing
and stale identities at both API boundaries, and preserve an independently
replaced target after linking with failed release status. Existing source
replacement, collision/symlink, cleanup/rollback/diagnostic failure, first-error,
unknown-identity retention and actual-finally preservation regressions still pass.

Focused tests ran first: Linux 13 passed; Windows 13 run, no failures, one
symlink-capability skip. Full Windows Python 3.12.10 / NumPy 2.3.0 suite: 134 run,
no failures, five skips (real-Poppler full owner; sidecar symlink; two existing
publication symlink subcases; publication symlink test). The pinned external
reference was read-only. This includes P1/P2, existing project/transfer/serializer,
normalization-once, default-off, provenance, five-setting/profile/collector and
negative/failure regressions. Existing Ubuntu virtual-environment Python 3.12.3
full owner suite: 112 run, no failures, one real-pdfinfo skip. All three changed
Python files compile and pass Python 3.9 grammar checks; no farm Python execution.

NOT VERIFIED: actual Lustre/delivery-parent and installed farm preflight/runtime;
Linux NumPy-dependent numerical/export parity; real Poppler PDF/page checks;
unavailable Windows symlink capabilities. No packages, environment repair,
assertion weakening, staging, commit, push, farm execution or acquisition authority.
Memory health has only the existing nonblocking CURRENT/MEMORY soft-size warnings,
reserved for the next checkpoint/milestone audit. Six changed paths in this pass:
owner, owner tests, Check 1 exporter tests (required signature wiring only),
CURRENT, this investigation and manifest. The existing contract and its earlier
pending amendment are byte-preserved and supplied separately as snapshot/context;
this correction introduces no contract edit. Deliver complete six-path diff,
pre/post source snapshots, actual logs and health/preservation receipts. After
later independent acceptance, user publication and pushed-state review, the next
farm action is the SAME preflight-only gate, never automatic acquisition.

## Independent identity handoff review: source acceptance

`SOURCE REVIEWED`: bounded source review passes the cumulative hard-link /
pending-ownership repair at test 374bed263e077c356ea06ae7eea6416886b61d4a.
Submitted six-path diff: 78009 bytes, SHA-256
971f96b1da6224d99ade353680b85191b72f9b61e1f9833cf13cdd7c756d5924.
Hash, read-only reverse application and diff checks matched actual files. Supplied
review ZIP CRC, every member hash/size, six-path diff and current snapshots were
independently verified. No additional blocking source finding in this review.

The validated tuple now passes through execute_gate -> release_check1 -> hard-link
helper. Missing or mismatched authority fails; stats only compare against the
trusted tuple. Actual-final-guard identical-byte replacement now fails with no
final ZIP, replacement survival and failed owner status. New missing-identity and
post-link target-replacement tests retain existing cleanup, rollback, collisions,
first-error, receipt and preservation checks. Unknown retention and stable
open-ZIP validation remain intact. Scientific/P1/P2 code and normal modes are
unchanged. Five-setting tests only adapt to the required identity parameter.

Reviewer focused Linux: 13 passed, 1.025s. Reviewer focused Windows: 13 run,
no failures, one symlink skip, 0.843s. Reviewed supplied full logs identify the
matching current source: Linux112/1 pdfinfo skip, Windows134/5 capability skips,
no failures. Those full suites were NOT rerun by this reviewer. Independently
compiled and 3.9-grammar-checked the three changed Python files without bytecode;
this is not farm Python 3.9 execution. Manifest/ordinary health pass with existing
nonblocking soft-size warnings, scheduled for checkpoint/milestone.

Publication reconciliation: the latest correction has six paths, but the earlier
52-line Lustre contract amendment is still uncommitted at this base. Its bytes
match the previously reviewed contract and supplied context snapshot. The current
manifest hashes it. The complete accepted publication cohort is SEVEN paths:
the submitted six plus the existing contract amendment. No additional source
repair/contract cycle is needed; publish that preserved authority with the repair.

Single NEXT: user publication -> independent pushed-state synchronization ->
reviewed SAME preflight-only gate. Actual Lustre/delivery-parent, installed farm
runtime, Linux NumPy parity, real PDF/page and unavailable capability conditions
remain NOT VERIFIED. No acquisition approval or farm command follows. Source,
contract, unrelated pending bytes, previous bundles, index and refs remain
preserved; no stage/commit/push/farm/ROOT operation by Work. Prior P1/P2 acceptance,
both Center conclusions, Outcome B, Method A/B and accepted closures stand.
External package: KaonLT_Check1_Identity_Handoff_Acceptance_20261106.zip.


## November 6 optimizer-reference destination collision repair

`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; actual-diff review pending.
The user reports a November 6 farm `check1:destination_exists`: Step 4 optimizer
reference acquisition wrote a SIMC packet before legitimate Step 6 publication
at the same setting/epsilon basename. The supplied failed gate-status fields,
quoted in the user's provenance correction, identify source commit
`2ebcf8e45c4a3aa40737e37fb0700e976d6c8698` and attempt
`kaonlt-e8-4-canonical-five-analysis-a95bbrb9`. Failed stage: `completion_markers`;
analysis return code: `1`; owner invariant: `analysis_failed`; underlying exception:
`check1:destination_exists`. Scientific evidence is inadmissible and not accepted.

`SOURCE VERIFIED` at local test 2ebcf8e45c4a3aa40737e37fb0700e976d6c8698:
both `_build_simc_reference` callers use ordinary `find_yield_simc`, whose
iteration-zero Check 1 guard previously admitted the optimizer. The existing
contract was amended before source edits with an exact nine-path allowance.
The copied candidate input now sets `check1_export_context` to
`optimizer_simc_reference`; only the Check 1 SIMC sink/export guard consumes it.
Original inputs are unchanged; no prepass physics context is reused. The caller
is added to exact Check 1 source-admission and executable-hash inventories.
Ordinary owner inventories/modes, schema and the publication helper are frozen.

The real optimizer caller -> find/calculate/bin/process SIMC -> event loop ->
native closure -> authority -> sidecar writer sequence executes in new local
controlled-event/histogram tests for low/high epsilon. Before the fix both
published during Step 4 (2 expected failures); default-off and ordinary-scope
controls passed. Afterward repeated Step 4 calls publish zero, Step 6 publishes
one JSON/NPZ pair, and duplicate final calls fail `check1:destination_exists`
with both prior hashes unchanged. Per-cell yield/error remain 6/6 from raw
iteration weight 2 and normfac 3; MM/t/Q2/W histogram contents and errors agree.
ROOT histogram/geometry adapters are controlled fixtures, not farm validation;
normalization provenance unavailable in that fixture stays explicitly unavailable.
Existing real data/Local-t source projection, signed/Sumw2, P1/P2, five-setting,
source/provenance/default-off and collector tests also pass where executable.

Adjacent source audit: both optimizer SIMC reference sites route through the
new context. Repeated data candidates already use the real existing prepass
configuration and data-observer guards; a new regression evaluates those actual
guards for both epsilons. Main Step 6 has one data and one SIMC call; main_auto
and main_iter SIMC calls use iteration=True and remain excluded. High-epsilon
selected-data rerun routing and the distinct parent observer are unchanged.
No additional collision path was identified within this bounded source audit;
full main/ROOT/farm runtime remains NOT VERIFIED.

Executed local evidence: Windows Python 3.12.10 / NumPy 2.3.0, focused 5 passed;
affected suites 144 run, 5 skips (one real Poppler PDF-inclusive five-setting
check, one sidecar-symlink test, two older high-publication symlink subcases,
one owner-publication symlink test: Windows privilege unavailable). The pinned
external ZIP test passed read-only after authorized access; its fixed digest
remains 02ff71998fd51d5e56f07484bb1c71185f78ed51d6ce409b0d34154bda7a6d7a.
Linux WSL Python 3.12.3 / GCC 13.3.0 owner suite: 112 run, 1 skip (real Poppler
pdfinfo absent), including publication ownership/rollback/symlink checks.
Numerical WSL tests require NumPy, absent in the existing WSL environments:
NOT VERIFIED there; no packages installed or scientific inputs substituted.
Five changed Python files compile and parse with Python 3.9 grammar; this is
not farm Python 3.9 execution. Exact normalized-source comparisons prove only
the context assignment, Check 1 guard condition and three inventory additions
changed executable source; all prior Check 1 tests are retained unchanged.

Sole NEXT is independent actual-diff review, then conditional user publication,
pushed-state synchronization and the SAME reviewed preflight-only gate. Actual
private/delivery filesystem capabilities, farm Python/ROOT/dependencies,
path/alias and executable/source/provenance/manifest checks remain mandatory
before any separately authorized scientific acquisition. No commit, stage,
push, farm run, acquisition authorization or failed-farm artifact delivery.

## Independent optimizer context review: source acceptance

`SOURCE REVIEWED`: bounded nine-path source review passes at test
2ebcf8e45c4a3aa40737e37fb0700e976d6c8698. Submitted diff: 29288 bytes,
SHA-256 e2f7abfcb7daf3caab39472be1b63e394f3cc87e77768e9c51a5cdfab7769d24.
Hash, read-only reverse application and actual nine-path inventory match. Review
ZIP CRC, diff/current snapshots and every declared current/protected/unrelated
file hash were independently checked. No blocking source finding in this scope.

Both Step 4 paths use _build_simc_reference, which adds the explicit context only
to its existing copied input. Only the default-off Check 1 SIMC sink/export guard
consumes it. Scientific prepass context is not reused; existing SIMC numerical
loop, histograms, normalization, weights and yield extraction still run. Original
production input remains unchanged and final Step 6 publication remains eligible.
Repeated optimizer calls emit zero packets; final emits one pair; repeated final
publication raises destination_exists with both previous packet hashes intact.
The owner/helper inventories and exact Check 1 allowlist now cover the actual
new caller. Publication/rollback, writer algorithms, collector/schema and prior
scientific/P1/P2 code remain unchanged. Existing data preparent guards remain.

Reviewer Windows Python 3.12.10/NumPy 2.3.0 executed all five new sequence tests
and all five existing preparent-contract tests: 10 passes, no failures/skips.
Both epsilons, repeated references, default-off, ordinary scope, actual adjacent
guard expressions, unchanged input, yields and MM/t/Q2/W bin contents/errors are
covered with controlled ROOT/cut/geometry adapters. These establish local wiring
and numerical checks, not full ROOT/farm or physical validation. Independently
compiled and 3.9-grammar-checked all five changed Python files without bytecode.
Supplied full Windows144/5 skips and Linux owner112/1 pdfinfo skip logs were
reviewed, not independently rerun here. Manifest/ordinary health pass; existing
soft-size warnings remain nonblocking and scheduled for the next milestone.

The supplied failed gate-status identifies source
`2ebcf8e45c4a3aa40737e37fb0700e976d6c8698` and attempt
`kaonlt-e8-4-canonical-five-analysis-a95bbrb9` for the `check1:destination_exists`
failure documented above. The earlier missing-SHA claim is corrected; this failed,
inadmissible and unaccepted attempt does not establish runtime validation of the
independently SOURCE REVIEWED local repair.
Actual Lustre/delivery-parent, farm interpreter/ROOT/dependencies/path/provenance/
manifest readiness, Linux NumPy parity and unavailable PDF/capability gates remain
NOT VERIFIED. No source edit, staging, commit, push or farm operation by Work.
Source/contract, unrelated pending bytes, prior bundles, index and refs preserved.

Single NEXT: user nine-path publication -> independent pushed-state synchronization
-> SAME reviewed preflight-only gate. No new source task, environment repair or
normalization study is needed. Acquisition still requires reviewed successful
preflight and separate user authority. Prior P1/P2, both Center conclusions,
Outcome B, Method A/B and accepted closures stand. External package:
KaonLT_Check1_Optimizer_Context_Acceptance_20261106.zip.


## November 8 scientific-identity failure diagnostics

`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; independent ChatGPT actual-diff
review remains required. Local test base is
`912ca31d0afdd7b9757395d205d7c54bc49018e8`, matching the November 8 diagnostic
contract. Pre-existing memory/offline-testing/review work is preserved.

The contract supplies the November 7 attempt
`KaonLT_Check1_Acquisition_Q4p4W2p74_20261107T163552Z-3952780`, source as above,
and frozen reference SHA-256
`02ff71998fd51d5e56f07484bb1c71185f78ed51d6ce409b0d34154bda7a6d7a`.
Reported analysis return code 0 and completed analysis are distinct from the
failed owner: stage `verify_artifacts`, first invariant
`check1:scientific_identity:Left:lowe.cells[0].fit.pion_component_alignment.applied_components.pi_delta.parent_alignment_hash`.
Artifact verification, collection and ZIP verification did not complete; both
scientific acceptance fields remain false. The destination collision did not
recur. These are contract-supplied failed-evidence facts; this implementation
neither independently parsed the November 7 farm artifacts nor recovered their
missing actual hash. Failed artifacts remain diagnosis-only. Cause NOT VERIFIED.

SOURCE VERIFIED locally: the unchanged comparator's failure is intercepted at
the existing packet gate; bounded diagnostic-only values and first differences
are written through the existing owner GateStatus before private cleanup. The
same exception object is re-raised. The receipt retains source/reference roles,
cell/t/phi/component, complete differing hash strings, parent and child selection,
parent/child/inherited-hash consistency, exact parent-preimage operand summaries
and complete canonical digests. Parent hash reproduction mirrors the unchanged
producer for saved JSON-native values. Exact operand comparison has no tolerance;
no recomputed hash is substituted into acceptance or scientific output. Missing
parents and inconsistent links are explicit. Diagnostic construction/write/stderr
failures cannot replace the original scientific error. Success remains silent.

Controlled regressions serialize a child through the real JSON/coordinate and
packet gate, write the actual owner status, remove private inputs and read the
surviving failed receipt. They cover checksum and sub-tolerance shift changes,
missing/corrupt parents, inconsistent lineage, large identifiers, unrelated
float/key/list/setting/type failures, unchanged successful status, exception
identity, and diagnostic/status/stderr failures. Frozen comparator, projection
and coordinate-validator functions are unchanged; the alignment producer and
all protected scientific source bytes match startup. No fitting, configuration,
threshold, normalization, template, weight or publication/collector change.

Local checks: supported Windows Python 3.12.10 / NumPy 2.3.0 suites ran 148 tests,
no failures, five capability skips (real Poppler PDF-inclusive owner, sidecar
symlink, publication symlink, two older parent publication symlink subcases).
Linux WSL Python 3.12.3 / GCC 13.3.0 full standard-library owner suite ran 112
tests, no failures, one real-pdfinfo skip. An initial broad Windows invocation
included POSIX-only fixtures and an incorrectly named test module: 212 tests,
112 environment/loader errors, four skips; its log is retained and not called a
pass. Supported Windows suites and Linux owner coverage were then run separately.
No dependencies were installed or fixtures weakened. Linux NumPy numerical parity,
farm Python 3.9/ROOT/main, real PDF/page rendering and actual farm readiness remain
NOT VERIFIED. Three changed Python files compile and parse with Python 3.9 grammar;
that is not execution on Python 3.9. Existing memory soft-size warnings remain
nonblocking and scheduled for the next checkpoint/milestone.

Final focused rerun: 11 tests passed without skips, including all nine new
diagnostic cases and both immutable-reference tests. These additionally exercise
the real immutable reference's
Left/lowe cell-0 fallback/accepted-parent distinction and full inherited hash.
The first new reference test incorrectly assumed stored/recomputed parent hash
equality; it failed and was corrected to assert parity with the unchanged
producer while preserving both identities and the reported mismatch. The real
frozen parent preimage recomputes differently from its stored hash under that
producer definition; this is diagnostic-only metadata consistency, not a cause
for the November 7 comparison or authority to change the reference. The earlier
failed test log is retained. No comparator assertion or acceptance gate changed.

NEXT: independent actual-diff review, then conditional user publication and exact
pushed-state synchronization. Only a separately reviewed evidence-preserving
narrow diagnostic execution may then be considered. No second full acquisition,
Check 2, scientific interpretation or farm command is authorized. The original
scientific-identity invariant still blocks scientific acquisition; prior P1/P2,
both Center conclusions, Outcome B, Method A/B and accepted closed scopes stand.


### November 8 independent-review continuation: persisted bound and exact patch

`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; independent actual-diff review
is still required. The supplied review found one source blocker: the compact
64 KiB checks did not bound GateStatus._persist's indented UTF-8 text output.
Its exact fixture was reproduced locally: 19,244 compact bytes, 90,582 bytes in
the actual saved diagnostic. Earlier compact-size claims did not prove this gate.

The local owner now measures a same-depth one-key envelope using the writer's
indent=2, sort_keys/default ASCII escaping, UTF-8 and platform text newlines.
Envelope overhead is included conservatively. Oversized diagnostics summarize
inline containers with their existing complete canonical SHA-256 digests while
retaining full first-difference hashes, exact scalar operand differences,
parent/child selection and consistency. A further ancillary-context fallback
retains both first differences and complete context digests. The supplied
oversized fixture now persists 21,621 diagnostic bytes. The actual writer/readback
regression measures raw saved bytes without newline normalization, and every
controlled failed receipt checks the same bound and failed/no-publication state.
The original exception/invariant, science comparator, projection, coordinate
validator, parent-hash producer, reference and publication path remain unchanged.

The review also found that the original patch's contract addition had been
read using Windows-1252 instead of UTF-8. Replacement packaging uses Git's raw
subprocess bytes throughout, including the complete new-contract addition;
no shell redirection or default text decoding/re-encoding is used. Old review
bundles and the original contract bytes are preserved. The exact replacement
patch must pass forward/reverse checks and content reconstruction before delivery.

Continuation checks: Windows Python 3.12.10 / NumPy 2.3.0 affected suites ran
151 tests with no failures and five existing capability skips (real Poppler
owner page gate; sidecar/publication symlinks and two parent symlink subcases).
All 11 focused diagnostic cases passed. Linux WSL Python 3.12.3 standard-library
owner suite ran 112 tests with no failures and one real-pdfinfo skip. Real farm,
ROOT/main/PDF, unavailable Windows capabilities and Linux NumPy numerical parity
remain NOT VERIFIED. No dependency installation or scientific change occurred.
The initial new size test expected an empty rejection-reason list to remain
inline; compaction was narrowed to preserve empty containers, then all focused
and affected regressions passed. Existing memory soft-size warnings remain
nonblocking, scheduled for the next checkpoint/milestone.

NEXT remains independent review of the repaired full candidate and exact patch;
conditional later user publication and pushed-state synchronization authorize
no farm command or acquisition by themselves. The original scientific mismatch
and its cause remain unresolved; no science, reference repair or second full
acquisition follows. Method A/B, Outcome B and accepted closed scopes are preserved.


### November 8 v2-review continuation: preserve short first-difference windows

`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; independent actual-diff review
remains required. V2 review found that ancillary compaction also removed short
exact values in first_parent_operand_difference even though the receipt fit.
The supplied pi_delta.effective_windows fixture was reproduced: the v2 writer
saved 21,586 diagnostic bytes containing only the differing windows' digests.

SOURCE VERIFIED locally: both first_difference and first_parent_operand_difference
now remain exact while ancillary values are compacted, then ancillary context
is summarized if necessary. Only if the remaining receipt still exceeds the
actual persisted budget are inline first operands summarized, largest first,
with a size recheck after each removal. Complete canonical digests are retained;
other exact short values and full hashes survive. No scientific comparator,
projection, coordinate, producer-hash, reference or acceptance logic changed.

The real JSON/packet gate and unchanged GateStatus writer/readback now retain
expected [[1.1,1.16]] and actual [[1.1,1.16],[1.16,1.18]] at
parent_hash_operands.components.pi_delta.effective_windows. Saved diagnostic:
21,690 bytes, with the original inherited-hash failure invariant, failed stage,
failure reason, false admissibility/acceptance and no collection/ZIP/publication.
Additional saved-receipt tests cover short scientific first-difference lists
through ancillary-context fallback, and an operand whose own deeply nested
representation exceeds the budget: only that value is summarized, with its
complete canonical digest and the other exact short operand preserved.
All 14 focused diagnostic regressions passed locally. Affected Windows Python
3.12.10 / NumPy 2.3.0 suites ran 154 tests with no failures and five existing
capability skips (real-Poppler owner page gate; sidecar/publication symlinks and
two parent publication symlink subcases). Linux WSL Python 3.12.3 / GCC 13.3.0
owner suite ran 112 tests with no failures and one real-pdfinfo skip. Farm
ROOT/main/PDF, unavailable Windows capabilities and Linux NumPy numerical parity
remain NOT VERIFIED. These are deterministic fixture checks, not farm or
scientific validation. Prior bundles and the
unchanged contract are preserved. NEXT remains independent review of the complete
candidate/exact patch; no publication, farm command or acquisition authority.
The original scientific mismatch/cause, Method A/B, Outcome B and accepted
closed scopes remain unchanged. Existing nonblocking memory size warnings
remain scheduled for the next checkpoint/milestone.


## November 9 deterministic scientific-identity correction

`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`: local comparison-only candidate at test base
fde9b5002e08bf692d706a4303ee23b99edd7adf under the
[single November 9 contract](../phases/check1-deterministic-scientific-identity-fix-2026-11-09-task-contract.md).
The startup manifest omission was exactly this new contract (19170 bytes,
SHA-256 461ddafbb6f23f5631d78eaf5ef710e95642b6ccda66ea5d30aea4f6a92234fd).
The user explicitly authorized manifest reconciliation and resumption; regenerated
manifest/check/ordinary health passed with the existing nonblocking size warnings.
Unrelated worktree bytes and the empty index are preserved.

Evidence: supplied KaonLT_Check1_Independent_Root_Cause_Audit_20261109.zip,
SHA-256 15e345b7d637ded5fbb1914ebbdcc08726df34377bccc14ba90ff85c86e1af3b,
and read-only pinned reference
02ff71998fd51d5e56f07484bb1c71185f78ed51d6ce409b0d34154bda7a6d7a.
SOURCE VERIFIED: all 15 retained t-parent preimages reconstruct their stored
setting-wide ancestor hashes. Each Left/lowe t-bin yields
c03775871b5033a2d9682e174293ecefe599ebefb1a2e9a10cbc102dd223c82a.
The t-bin-0 local rehash is 92983dbde534bab2cd80aa229677658b75f6a0969997279ea1118472cafc6fce;
the child invocation-shift predicate rehash is
fc192ead77430b4adb0edecf80e8f2d241cd5e62585ffce113ace39441409f4d.
These are distinct scopes. The earlier November 8 statement of a frozen-parent
metadata mismatch used the wrong scope; it is superseded by this reconstruction.
The controlled UUID-only replay proves a sufficient mechanism, not the unique
historical cause: November 7 actual generated operands/hash remain unavailable.

The correction leaves producer/cache/output/provenance identities unchanged.
A dedicated Check 1 view validates each side's original setting-wide raw hash,
all t/child/applied links, component/name/checksum aliases and configuration
binding before deriving a deterministic cross-run hash. Separate parent-artifact
and embedded-parent alignment records must match exactly within a run.
The closed name grammar erases only the generated and clone 10-hex suffixes
of H_MM_component_shape_full pion names in the source-proven caller context.
Stable component/setting/context, axis and content/error checksum survive.
Stable t/phi names remain verbatim. Unknown generated syntax/location fails closed.
Alignment creation timestamps are validated as producer-format UTC and represented
only in this view as creation metadata. All other scientific fields remain;
the existing scientific_projection, compare_science and coordinate-validator
functions are unchanged. Both child and separate-parent owner seams use this view.
Bounded failure receipts distinguish setting-wide reconstruction from t-bin and
child-predicate rehashes and retain raw identity/selection context. No production
parent_valid, fallback, fitting, selection, normalization, yield or error change.

Final local evidence: Windows Python 3.12.10 / NumPy 2.3.0 combined suites ran
143 tests, no failures, four capability skips (real-Poppler whole owner,
sidecar symlink and two parent-publication symlink subcases). The tightened
UUID-only replay subsequently passed separately across all five child and
parent comparison paths; timestamp-only variation was tested independently.
The final diagnostic-scope/immutable-reference rerun passed 17 tests without
skips, including oversized receipts, missing-parent/operand and write-failure
cases. Six new real-reference tests cover 135 nominal cells, all 15 parents,
16 physical/schema negatives with regenerated raw lineage and 12 lineage/
identity-syntax negatives. Every negative reports its reproducible path.
Unchanged producer AST/hash reconstruction parity passes without ROOT.

Final WSL Python 3.12.3 / GCC 13.3.0 standard-library owner suite: 112 tests,
no failures, one missing-pdfinfo skip; WSL NumPy is unavailable. An interim
Linux run caught CRLF introduced by this implementation into the two edited
memory records: their raw bytes would change under Git clean conversion.
Only those allowlisted records were restored to required LF, the manifest
regenerated/checked, and the full 112-test Linux suite passed. The failed log
is retained for review, not counted as a pass. No test assertion was weakened.
All three changed Python files compile and parse with Python 3.9 grammar;
this is not farm Python 3.9 execution. Underlying projection/comparator/
coordinate-validator ASTs, seven frozen producer/source file hashes, original
reference bytes, 35 pre-existing pending files, empty index and refs are
preserved. git diff --check passes. Manifest and ordinary memory health pass
with the documented nonblocking size warnings.
Local validation is not ROOT/PyROOT/main/PDF/farm acceptance. Source review is
pending; no scientific closure or readiness is inferred. CURRENT owns the single
push-stable NEXT: review a representative bounded early runtime identity diagnostic
conditional on source acceptance, user publication and exact pushed-state review.
Such an equivalent entry point is not proven; preflight-only is insufficient.
No full acquisition, Check 2, Method A/B change or Outcome B change is authorized.
Existing size warnings remain nonblocking, scheduled for the next meaningful
checkpoint/milestone audit; no unrelated cosmetic consolidation follows.


### November 9 independent-review P1/P2 continuation

`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`: same-contract narrow repair; replacement actual-diff review remains
required. Independent review found that candidate B could change raw UUID/name
aliases while retaining baseline A's ancestor string, and that a pre-projection
raw rejection could be explained by an unrelated raw-document UUID difference.
The frozen producer explicitly requires candidate alignment equality to baseline
A (candidate_alignment_changed). The verifier now checks that equality for every
t parent and both artifacts/sides; both raw sides pass before either projection.
Generated whole-run UUID changes remain permitted only with valid raw lineage.

Raw-lineage/syntax diagnostics now preserve the actual exception invariant,
rejecting side, parent/cell location and available raw operands. No raw-document
science first difference is computed on this branch; projected first difference
and cross-run parent operand difference are unavailable rather than causal claims.
The legitimate UUID plus broken third-parent applied pi_delta regression must
identify parents[2]. Existing bounded receipt, exception identity, false
acceptance and no collection/ZIP/publication boundaries remain.

The submitted test-file raw after-snapshot had CRLF while its Git patch used LF.
Only this allowlisted test source is normalized to LF; the replacement package
must verify raw/Git-clean equality and byte-exact forward/reverse reconstruction
for every task path. Prior snapshots/bundles and unrelated files remain preserved.
Replacement validation: Windows Python 3.12.10 / NumPy 2.3.0 combined reference,
all-five/135-cell, physical/lineage-negative, packet, parent/child and collector
suites ran 146 tests without failures, with the same four capability skips.
All 19 final affected causal-diagnostic, candidate-B and immutable-reference
tests passed without skips. WSL Python 3.12.3 / GCC 13.3.0 standard-library
owner suite ran 112 tests without failures, one missing-pdfinfo skip; WSL NumPy
remains unavailable. Three changed Python files compile/parse with 3.9 grammar,
not farm Python 3.9 execution. An initial two-test focused run retained two
exception-identity assertion failures because its injection still mocked the
old later projection seam. The fixture now injects at the new first raw check;
the same assertIs assertion and all actual rejection/receipt checks pass.
The failed log is retained, not counted as a pass. No assertion was weakened.
Raw/Git-clean source equality, frozen scientific-source/function identity,
full exact patch reconstruction and manifest/ordinary health are verified in
the replacement review receipt. Existing nonblocking memory size warnings stay
scheduled for the next meaningful checkpoint/milestone. No scientific producer, reference, numeric comparator,
coordinate contract, Method A/B, Outcome B or farm authority is changed.


### November 9 v2-review embedded/separate binding continuation

`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`: same-contract error-path repair. The user reports acceptance
at source-design level for all other identity behavior. The remaining binder
correctly rejected a changed separate-parent alignment but omitted its known
operands; a missing alignment instead escaped as KeyError. Exact equality is
retained. Structural checks now identify embedded/separate artifact, parent
index, branch and available/missing operands in scientific-lineage ValueErrors;
unequal alignments attach their bounded exact first difference. Ancillary
receipt context handles malformed records without replacing that evidence.
No producer, reference, projection, tolerance, fitting or acceptance change.

Three initial focused tests pass: the real pinned third-parent hash/missing-
alignment cases pass the existing legacy parent/child JSON checks before the
binder; twenty malformed-record cases cover both artifacts. The real owner
writer preserves the same exception, failed verify_artifacts status, false
acceptance, bounded receipt, private cleanup and no ZIP/publication. Final Windows Python 3.12.10 / NumPy 2.3.0 relevant suites ran 149 tests,
no failures, four unchanged capability skips (real Poppler owner page gate,
sidecar symlink and two parent-publication symlink subcases). Exact source
hashes are recorded before/after the run. WSL Python 3.12.3 / GCC 13.3.0 owner
suite ran 112 tests, no failures, one missing-pdfinfo skip; WSL NumPy remains
unavailable. Both actual and reference binding-error receipts are covered.
Seven frozen-source hashes, accepted owner source and original comparison/
coordinate functions remain unchanged. Changed Python compiles/parses with
3.9 grammar, not farm 3.9 execution. The replacement receipt verifies raw/
Git-clean bytes, exact forward/reverse patch reconstruction, unchanged
reference/unrelated state and manifest/ordinary health. Existing two size
warnings remain nonblocking, scheduled for the next meaningful checkpoint.
The only next gate remains CURRENT's conditional narrow readiness review;
no publication, farm execution, new acquisition or scientific closure follows.


## November 9 continuous-t evidence and PID follow-up

The published comparison-only correction at
`6fb1a94ed8f6b7af61e130c9b42b9762164444b9` matches the seven independently
reviewed v3 blobs. The corresponding acquisition
`KaonLT_Check1_Acquisition_Q4p4W2p74_20261109T181613Z-98329` reports owner
success through collection, cleanup, preservation and publication. Direct
review of the primary ZIP, log and summary matched the owner hashes, including
ZIP SHA-256 `4fcc463bae457fdd49f5fef5cb814143df9ce7313be73623a7d674d123554213`.
Independent readback checked all 35 artifacts, both scientific comparison
paths, ten numerical packets and the October data/dummy and unchanged SIMC
bindings. This verifies the retained diagnostic payload; the owner's stored
admissibility and scientific-acceptance flags remain false.

`KaonLT_Check1_Full_t_Matched_Interior_20261109.zip` retains the continuous -t,
occupancy and matched-interior results. Center prompt and SIMC have broad
observed t3 support across the retained canonical interval. With the Lambda
window, empirical common coverage is 97% at low epsilon and 100% at high
epsilon at the 0.005 GeV^2 resolution. These are occupancy measurements and
do not isolate geometric acceptance. Center/low broad-MM t2 has a very large
negative signed subtraction that needs reconciliation against the source
weights, application and plotting. Center/low t3 has weak and poorly
resolved Lambda signal despite its measured -t occupancy. The prior
lower-t3 tail physics interpretation and the resulting conclusion against
an upper-edge explanation are withdrawn. Center/high t3 has a potentially
substantial applied pion contribution whose physical accuracy is unverified.

The November 9 user decision adds an acceptance-dependent PID investigation
within Check 1, using established detector-response observables before any
Method A consideration. Separate the operational zero, low-positive and
strong-positive HGCer categories and the training and physical populations.
Compare all five settings where supported; prioritize the Center t3 anomalies.
CURRENT owns the active action. Earlier phase closures and the six roadmap
statuses, unresolved Outcome B and canonical-five equivalence, detached
Method A and diagnostic/excluded Method B retain their recorded boundaries.

## November 9 Check 1.B results and Check 1.C decision

The two November 9 analytical packages passed member-hash readback. The
supplied primary acquisition archive, log, summary and delivered status
also match the preceding readback receipt. This confirms those exact bytes;
the reported owner outcome and the stored false evidence-admissibility and
scientific-acceptance flags retain their separate meanings.

The retained Center/low t2 broad-MM arithmetic reproduces +1.405651 before
pion subtraction, +148.831062 applied pion contribution and -147.425411 after.
The Lambda-window result remains positive (+0.715728). The maximum transfer,
about 3435.63, is associated with a small fitted pion-control denominator
below the Lambda window. See the [October investigation](baseline-pion-control-center-lowe-regression-2026-10-14.md)
for the staged/joint-fit, source-aware statistics and component history.
Different cohorts do not establish a regression from their transfer maxima.

Five-setting joint pion-control records show SHMS delta and HGCer-position
changes with continuous -t. The conditional low-positive training fraction
does not establish a universal Center upper-t3 increase. About 90-97% of
Center/high t3 signed beneath-Lambda pion background falls in the retained
negative-delta grouping (-10 to -1 percent) across the three subintervals.
Prompt, random, dummy and dummy-random identify experimental sources. The
total-background tables do not separate pi_n, pi_Delta and pi_SIDIS. The
analytical packet lacks matched detector rows for the cleaned kaon-selected
events; this is a package-specific limitation. Weak Center/low t3 Lambda
support and unresolved Center/high beneath-Lambda background remain.

The user selects the five-setting matched event-level kaon/pion PID and
three-channel SIMC acceptance comparison within Check 1, including the
hypothesized SIDIS onset near zero SHMS delta. CURRENT alone owns NEXT.
The six milestones, historical phase closures, Outcome B and canonical-five
F.4 blockers, detached Method A and excluded Method B retain their boundaries.

## November 10 ten-bin channel constraint

Check 1.C verified 188253 joint physical pion-control rows and reconstructed
all 15 retained Local-t transfer arrays. On shared joint detector support, the
descriptive low-positive fraction rose from middle to upper t3 at Center/low
(0.664 to 0.809 percent) and fell at Center/high (1.356 to 1.046 percent).
Neither result calibrates true-pion leakage. Matched kaon detector observables
and channel-resolved SIMC delta rows were unavailable in the reviewed analytical
evidence. Center/high t3 fitted kaon-side SIDIS amplitude was zero; its nonzero
pion-side amplitude contributes to the control denominator. The applied
channel allocation is algebra on the same signed experimental rows. Weak
Center/low t3 Lambda support and the earlier Center/low t2 denominator finding
retain their respective scopes.

The user directs a five-setting comparison of the three existing fitted
pion-control components to signed experimental pion control, preserving
joint delta, -t, MM and phi. Every delta diagnostic uses the established
ten proton-cleaning bins with edges [-10,-7,-4,-1,2,5,8,11,14,17,20] percent.
Keep the empty cells, native MM edges and original parent normalization.
Use the pion-side B amplitudes without refitting. Quantify the SIDIS
neighboring-t stability and differing Delta shapes as hypotheses; the
[-1,2) percent cell cannot resolve an onset exactly at zero. Use readable,
visually reviewed figures and retain signed values. Do not repeat the
Center/low t2 denominator investigation or treat loader fallback details
as a reason to reject shape comparisons. No experimental event is labeled
as a species from the fitted allocation.

This continues the existing six-step roadmap. Outcome B and the canonical-five
equivalence gate remain blocked, Method A stays detached, and Method B stays
diagnostic/excluded. CURRENT owns the sole NEXT.

## November 10 joint pion-component acceptance constraint

The ten-bin archive passed hash and numerical-receipt readback. It preserves
188253 observed pion-control rows, 300 setting/t/delta/selection cells and
45 channel parent closures on phi and MM. Center/high t3 has 90.6 percent of
signed Lambda-window pion control below delta=-1 percent. No channel-specific
delta templates passed admission in that archive. These are diagnostic results.

The user authorizes a detached joint (MM,delta) candidate across all 15 parents:
one constrained amplitude per channel across the exact ten delta bins within
a parent, with existing normalization and no cross-setting equality. First
reconstruct original parent templates from source-bound delta slices plus
explicit flow. Compare fixed-model residuals before the candidate fit. Reusing
the same pion-control data is appropriate for development if fit-data dependence
and source covariance are accounted for. Evaluate finite-SIMC, alignment and
acceptance uncertainties separately from conditional counting precision. A
zero fitted kaon-side SIDIS term is a physical-model concern; pion-control
identification does not calibrate the kaon-side leakage. Preserve Outcome B,
all six roadmap states, historical denominator evidence and the existing
Method A/B boundaries.

### Closing offline checkpoint

The all-five-setting joint (MM,delta) candidate remains BLOCKED at template
admission. No exact channel h10 event file or equivalent same-event ssdelta
projection was available in the inspected retained evidence. The observed
ten-bin spectra remain usable, with empty and unsupported states kept separate.
No ten-bin model residual, fit, uncertainty reduction or kaon-side correction
is inferred.

A separate phi/MM local sensitivity calculation covers all 15 parents and
uses common observed support with the same fixed pion-side amplitudes.
It propagates original signed counting errors and conditional finite-SIMC
variance through the shared parent normalization and fixed alignment.
For Center/high t3 complete Lambda bins, the least-constrained MM contrast
has a conditional standard error of 3.168 in MM
and 0.812 in phi/MM; including conditional finite
SIMC statistics gives 3.366 and
0.892. This measures phi information only.
Data-to-fit covariance for independent validation and acceptance/alignment
model covariance remain unavailable. Outcome B still requires independent
kaon-side anchor/gauge information and applicable covariance.

## November 11 regenerated pion SIMC projection compatibility

The user intentionally reran all 15 neutron, Delta and SIDIS channels.
Different hashes, entry counts and accepted histogram contents are expected.
The November 9 owner provenance remains historical. Its resolved originals
were shared background-sample files; the E.8.4-named worktree was the successful
November 9 attempt's temporary analysis alias. The subsequent handoff used
fit-packet alias paths as its primary lookup and shared storage as a fallback.

The uploaded numbered (5) handoff is byte-identical to the original 214916-byte
ZIP, SHA-256 584b7a1cfb2627a515df18be17d77b26fdfbec0d7d33fdad52a8236d59e54233.
Its extractor, reference builder and all six members were inspected. The
builder continues to construct the historical reference. Compatibility is
implemented in the reviewed candidate testing extractor with an explicit
reference option and regenerated input root. Default historical hash,
entry-count and histogram checks remain. Regenerated mode records fresh
file identities, tree/branch metadata and before/after checks; it constructs
new direct parents with the existing scalar integrator, preserves fixed
alignment and parent normalization, and closes the ten delta slices and
flows against those same-cohort parents. Historical differences are separate
receipt fields. This establishes no historical-template equivalence or
fitted-mixture acceptance.

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING, subject to independent actual-
diff review. Windows Python 3.12.10 / NumPy 2.3.0 ran 26 focused tests and
30 existing integrator/observer regressions, 56 total, no failures and two
existing symlink-capability skips. Real historical arrays replayed 45
content/Sumw2 closures. Python 3.9 grammar/compilation passed. Tests cover
changed hashes/entries/weights, signed contributions, flows and endpoints,
all-15 mocked execution, changed inputs/reference, alias duplication,
source pins, exclusive publication and rollback through context exit.
Actual regenerated ROOT reads, PyROOT and farm behavior remain NOT VERIFIED.
The production source and baseline, all six roadmap states, Outcome B,
detached Method A and diagnostic/excluded Method B retain their boundaries.
No Git publication, farm operation or new scientific fit was executed.


## November 11 fixed-model MM-delta comparison

The user requested fixed-model residuals before any detached joint fit.
Local and public test HEAD and supplied projection runtime HEAD agree at
f346579713168bb9c36b48128df51117c83059ae. Six published snapshots match the
reviewed corrected candidate. The farm projection ZIP SHA-256 is
17244335766882fceaf7c42e74e6f2c3509f759b0b78c774ae3e596fb17e709c.
Independent readback verifies 15 source-bound roles, 45 raw/aligned closures,
canonical-phi and explicit flows; 188253 observed source/entry records bind
to the November 9 acquisition. New hashes were admitted as regenerated inputs.
Measured historical raw-parent differences are zero, with aligned differences
at most 3.42e-15; historical equality was not required.

The fixed-B comparison covers all five settings and 15 parents, preserving
ten delta bins, signed sources, existing shifts and parent normalization.
Center/high t3 broad native MM [0.7,1.3) gives data 5.094172 +/- 0.083070,
model 4.997829 and residual +0.096342 +/- 0.083070. Delta [-10,-7) gives
+0.350696 +/- 0.046359; [-4,-1) gives -0.174789 +/- 0.037512.
All five t3 parents have positive broad-MM residuals in [-10,-7).
Complete native Lambda bins [1.102,1.156) give data 0.417555,
model 0.443693 and residual -0.026137 +/- 0.026175 at Center/high t3.
Exact [1.10,1.16) data are 0.469769 +/- 0.025871; model boundary-bin content
bounds are [0.443693,0.549561], without fractional fills.

Center/high t3 has 186/1000 occupied native broad MM/delta cells, six negative
signed cells, and three occupied cells without retained MC support. Lambda
has 28/90 occupied cells; four empty cells carry model yield totaling 0.002487.
Empty data, finite-MC absence and physical acceptance are distinct.
All observed/model delta-flow counts are zero. Native-MM flows remain explicit.

Errors are conditional signed counting plus finite-MC propagation through the
shared stochastic parent normalization and fixed interpolation. Retained B
covariance matches the applied amplitudes, but data-to-B, MC-to-B, upstream
cleaning, acceptance/alignment and shared systematic covariance remain absent.
Residual quadratics are descriptive, without independent-fit p-values.
Joint designs have numerical rank three in all parents/regions, with weak
physical separation: Center/high t3 broad Delta/SIDIS local correlation is
-0.920; Lambda joint local errors (n,Delta,SIDIS) are (0.727,1.242,1.455).
No fitted values or uncertainty reduction have been adopted.

The retained scientific bundle is
KaonLT_Fixed_LocalT_Pion_MM_Delta_Comparison_20261111.zip, with original plots,
all native-cell residuals/support/errors, numerical inputs, source identities,
replay scripts and independent checks. Checks include 30600 exact observed
native-cell readbacks, 15 original fit-input histograms, 45 derivative checks,
90 signed-source aggregate checks and 30 parent residual sums.
Independent fixed-model review is NEXT. No detached fit, source change,
Git publication or farm execution was performed by the reviewer. Outcome B
remains unresolved, Method A detached and Method B diagnostic/excluded.


## November 11 detached joint MM-delta fit

The user accepted the regenerated 15-source/45-closure farm projection and the
all-five/15-parent fixed-model comparison as diagnostic evidence. Their joint-fit
handoff authorizes revised relative amplitudes as the primary goal; possible
acceptance/shape mismatch is a model-uncertainty concern, not a deferral gate.
The implementation follows the single joint-fit contract at starting test HEAD
f346579713168bb9c36b48128df51117c83059ae. Production and all six roadmap states
are unchanged. The new testing CLI is ROOT-independent and detached.

One nonnegative amplitude per channel is shared across ten delta bins within
each parent; settings/parents are independent. Original staged targets,
regularization widths, MM fit domain, positive finite signed Sumw2 rule,
alignment and single parent normalization remain. The primary quadratic uses
counting plus correlated finite-MC precision fixed at original B. All 15 original
marginal joint fits replay within 1.02e-8 before joint fitting. Empirical empty
and unsupported cells are retained separately without fictitious variance.

Center/high t3 fixed (n,Delta,SIDIS)=(0.235833,1.777489,2.984507) becomes
(0.189216,1.395621,3.208202). Conditional estimator errors are
(0.044257,0.237827,0.256694); Delta/SIDIS correlation is -0.93855.
Counting errors are (0.034146,0.183664,0.197950), finite-MC errors
(0.028155,0.151093,0.163425). MC propagation includes shared parent normalization,
fixed interpolation and the design-transpose/residual term in the differentiated
normal equations. Within-candidate regional residual errors include data/fit and
direct-template/estimator cross terms, conditional on stored priors and precision.
Old-versus-new fit covariance, prior/data reuse covariance, upstream cleaning and
physical acceptance/alignment covariance remain unavailable; no independent
validation or calibrated confidence interval is claimed.

The same-precision data residual quadratic changes 306.2545 -> 296.3345 (3.24%).
Native complete Lambda model integral changes 0.443693 -> 0.439954, with signed
data 0.417555 and fitted residual -0.022399 +/- 0.023905 conditional in-sample.
Its fitted (n,Delta,SIDIS) contributions are (0.027441,0.074970,0.337543).
The broad model integral changes to 4.793039 against data 5.094172: objective
reduction is not integrated-yield closure. The applied production numerator,
denominator, background and event weights are unmodified.

For [-10,-7) delta at Center/high t3 the signed broad-MM fraction is
0.287128 +/- 0.007357 conditional counting, above every retained channel's
fraction. The maximum is SIDIS 0.231521 +/- 0.003762 conditional MC. Thus a
positive three-weight mixture alone cannot exactly reproduce the observed
fraction and total with these central shapes; the empirical gap 0.055606 is
a marginal total-variation discrepancy lower bound, not calibrated physical
acceptance covariance or an adopted correction.

Leave-one-occupied-delta refits span n=[0.136914,0.247227],
Delta=[1.037405,1.545062], SIDIS=[2.718317,3.572385]. Retained non-rejected
alignment/baseline shifts span n=[0.166520,0.189216],
Delta=[1.252947,1.395621], SIDIS=[3.208202,3.367674]. These are empirical
influence/scan envelopes with no assigned probability law or quadrature sum.
Counting-only Center/high t3 gives (0,2.113750,2.480697); covariance evaluated at
the candidate gives (0.183024,1.368450,3.240491). No-prior and Lambda-excluded
sensitivities are also retained. Center/low t3 neutron is at zero; Center/low
t2 neutron is prior-dominated (retained width 0.001 around zero), with central
0.000129 versus unregularized 0.278196. Tiny conditional errors are not species
absence. All 15 constrained solutions, correlations, contours and residuals
are retained in the scientific bundle.

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING, subject to independent review.
Recommendation: further offline model/weight-robustness validation is warranted
after actual-diff/scientific review; no production adoption or farm integration
run follows. Full outputs, embedded project-data plots, exact source diffs,
before/after snapshots and deterministic/local numerical checks are packaged in
KaonLT_Joint_MM_Delta_Fit_20261111_review_bundle.zip. Outcome B remains unresolved,
Method A detached and Method B diagnostic/excluded. No commit, push or farm run.


## November 11 source-aware bootstrap and conditional transfer

The user accepted the detached joint-fit source-level implementation and requested
source-aware data/SIMC bootstrap, correlated model uncertainty and propagation
through existing pion transfer without waiting for externally owned PID results.
Public/local test remained f346579713168bb9c36b48128df51117c83059ae. The accepted
joint source/test/contract bytes match its review ZIP and remain unpublished.
The single bootstrap contract adds a detached CLI/tests; production, prior
pending source, all six roadmap milestones, index and refs are preserved.

1024 paired data-only, MC-only and joint replicas cover all five settings and
15 parents, with 46080 fits and 46080 parent/channel/phi closure checks.
One Poisson(1) multiplicity per retained cohort/source/PID-tree/entry is reused
across selections, parent labels, MM/delta/phi cells, fitting and application.
Independent MC channel/file-entry draws preserve raw parent/child relations.
Replica Sumw2 is sum(k*w^2), not sum(k^2*w^2). Fixed alignment and common
stochastic parent normalization are rebuilt. Counting and correlated finite-MC
precision are recomputed per replica at the original-B convention; retained
staged targets/widths and upstream calibration remain conditional/fixed.
No failed fits or unsupported broad/Lambda transfers occurred in the final run.

Center/high t3 amplitude bootstrap means (n,Delta,SIDIS) are
(0.152453,1.678388,2.826742), with SD (0.046726,0.283777,0.303005), compared
with unchanged detached central (0.189216,1.395621,3.208202). Empirical
Delta/SIDIS correlation is -0.956133. Bootstrap mean displacement is not silently
adopted as physical bias correction. Complete-native-Lambda control mean/SD is
0.409727/0.019403 versus central 0.439954. Exact-allcuts Lambda background
mean/SD is 0.065674/0.004331 versus detached central 0.061554 and unchanged
original applied 0.060786. Conditional background variance components are
data-only 1.3422952e-5, MC-only 4.4857441e-6, nonlinear interaction 6.2050247e-7,
paired cross 2.3114461e-7, total 1.8760343e-5. Their paired identity closes;
the nonlinear output covariance is not an assumed independent quadrature sum.

Protected kaon A coefficients are fixed, but the SAME bootstrapped SIMC shapes
enter numerator and denominator. SAME pion draws enter B fitting and event
application. Original transfer ratios and signed cached application replay for
all parents with the original floor/clipping/zero-flow convention. The explicit
clip_max sentinel is traced to owner-emitted None through snapshot(None), not
treated as an invented default. Post-cleaning kaon rows do exist and reproduce
the frozen kaon MM input; raw learned populations and common physical run/event
IDs across PID trees are not retained for the full learned/protected-numerator
refit here. Complete numerator-fit, training/prior/PID, shared-systematic and
physical-model covariance remain unavailable. Conditional transfer is not full
physical subtraction uncertainty; no independence or zero missing error is used.

Full MM/phi covariance is stored losslessly as centered replica factors with
explicit axes/common-support indices; exact Lambda/phi covariance is separate
because partial native-MM boundary bins cannot be fractionally filled. An
independent MC-only numerical replay checks all 15360 numerator/denominator
vectors (maximum denominator difference 2.23e-16). All 45 amplitude and 45
Lambda-background covariance readbacks and factor reconstruction pass.

Source-defined alignment/leave-delta response vectors and measured residual
bias directions preserve the negative-delta mismatch. Their outer products are
correlated response kernels without an assigned probability law; physical model
covariance remains null pending justified nuisance covariance. External PID
covariance/cross terms are explicitly null, not zero. The interface requires
matched axes/units/populations/lineage and overlap checks before combination.

Center/high t3 first-half/full amplitude SD differences are about 1.1/2.1/3.8
percent. Center/low t3 neutron is bound in 1018/1024 replicas; six positive-tail
outcomes control its variance, with roughly 54 percent half/full SD change.
That rare-tail precision is limited, not a stable Gaussian absence claim.
Center/low t2 retains its strongly constraining prior. Conditional empirical
percentiles are not calibrated physical coverage.

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING, subject to independent actual-
diff/scientific review. KaonLT_Pion_Bootstrap_Uncertainty_20261111_review_bundle.zip
contains source diffs, exact snapshots, all replicas/covariance factors, correlated
model responses, external-PID interface, original project-data plots and checks.
No production update, commit, push, farm run, Method A/B promotion or Outcome B
closure occurred. Independent bootstrap/scientific review is NEXT.


## November 11 publication checkpoint and uncertainty propagation

The user records independent source review of the detached joint fit and
source-aware bootstrap/conditional transfer: both are SOURCE REVIEWED. User
publication and pushed-state synchronization completed at test
5daceacf43be914de7cce2868fc9ff6dfc8f2c86; this continuation verifies live public
and local identity and exact reviewed source snapshots. The regenerated h10
projection is CLOSED / RUNTIME VALIDATED only within its accepted 15-source,
45-parent diagnostic-acquisition scope, as recorded by the supplied handoff.
No new farm execution or physical-subtraction acceptance is claimed here.
The earlier unpublished/development-review statements above describe historical
checkpoints and are superseded for current status by this paragraph.

ACTIVE: complete existing statistical and model-dependence propagation using
the retained 1024 replicas and alignment/delta-influence variations. Deliver
supported MM/phi and exact-Lambda sensitivity ranges without waiting for PID
or requiring a new calibrated acceptance-nuisance phase. Broad-MM negative-delta
residuals are not automatically Lambda uncertainty. Preserve all six roadmap
states/dependencies, Outcome B BLOCKED and detached Method A/excluded Method B.


### Completed five-setting statistical/model-dependence propagation

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING, subject to independent source/
scientific review. The new detached consumer retains the published 1024-replica
NPZ byte-for-byte and performs no new fits/resampling. It evaluates 116 already
approved alignment/delta-influence cases across all five settings/15 parents.
All broad and exact-Lambda transfers are supported here. Every central/model
denominator, Lambda integral and active MM/phi response replays against the
accepted bootstrap result. Protected A is fixed, with the same varied shapes
in numerator/denominator and the original signed event application unchanged.

Center/high t3 central exact-allcuts Lambda background is 0.0615541025,
conditional bootstrap SD 0.0043313212 and mean 0.0656744772. Existing combined
sensitivity range is [0.0583032301,0.0752890737], a central shift of
[-0.0032508723,+0.0137349712]. Alignment alone spans [0.0615541025,0.0752890737];
delta influence spans [0.0583032301,0.0679239910]. Complete-native-Lambda control
denominator central 0.4399542638 has existing-model range [0.397390,0.491410].
No model probability law, calibrated one-sigma or quadrature total is invented.

For the most-negative delta bin [-10,-7) percent, the central fitted residual is
0.388750 in broad native MM, -0.001978 in complete-native Lambda and 0.389176 in
complete high-MM bins above 1.16 GeV. Boundary-straddling bins remain separate.
The broad discrepancy is predominantly higher MM; its full magnitude is not a
Lambda uncertainty. The transferred model range follows actual coherent
scenarios, not direct allocation of the broad residual beneath the signal.

The package tabulates all-45 central/bootstrap amplitudes, full conditional
covariance/correlation, model amplitudes/boundaries, denominator/numerator MM,
joint control MM/delta, transferred MM/phi and exact Lambda/phi vectors. Full
statistical MM/phi factors and common-replica masks, aggregate MM/phi/Lambda
covariance and denominator/background cross covariance are preserved. Each
model scenario remains coherent; a pointwise envelope is not a single model.
Response co-movement matrices are descriptive geometry, not probabilistic
model correlation. No Center/low alignment alternatives exist in the retained
cases; that family is unavailable, not a zero error. Supported delta ranges
are delivered. Empty/signed/support states, Center/low t3 rare neutron tail and
Center/low t2 prior limitation are preserved without reopening their study.

Independent readback validates 90 all-parent covariance/cross-covariance
aggregates, 15 amplitude covariance matrices, 15 Lambda/phi variance sums,
and 116 direct signed-event transfer integrations from the original acquisition:
maximum transfer discrepancy zero. All 1024 empirical replicas and original
source identities remain unchanged. Deterministic suites run 145 tests with
two existing Windows symlink-capability skips and no failures. Python 3.9
grammar/compilation and exact patch reconstruction are separate local checks;
Linux numerical parity and ROOT/farm integration are NOT RUN in this task.

The review artifact is KaonLT_Pion_Uncertainty_Propagation_20261111_review_bundle.zip,
including exact eight-path diff/snapshots, source/member hashes, full numerical
arrays, eight embedded original project-data plots, CSVs, validation and
conditional user-only Git commands. Missing full protected-numerator/upstream
learning/physical-model and external PID covariance/cross terms remain null.
The supported statistical and sensitivity result is complete without waiting
for PID or a new acceptance-nuisance calibration. No physical coverage or
production-subtraction acceptance follows. All six roadmap states/dependencies,
independent Outcome B BLOCKED, Method A detached and B excluded are unchanged.
No staging, commit, push or farm execution occurred. Independent source/science
review is NEXT, then return to existing physical subtraction/final-yield gates
under their unchanged dependencies and user-owned publication/execution.


## November 11 Outcome B source/evidence audit and publication synchronization

The supplied continuation handoff records completed independent review and
publication of the joint fit, source-aware 1024-replica bootstrap and five-setting/
15-parent uncertainty propagation (116 coherent model cases): all SOURCE
REVIEWED. Local and live public test match d28ad6e7f620f0c4cc18c389373560f56fd9e84a.
All three review bundles pass CRC/member hashes; frozen source/test/contracts
match earlier reviewed snapshots, and all eight latest propagation paths match
published bytes. Earlier bundles' active-memory snapshots are superseded by the
latest checkpoint, not required to match it. The regenerated projection retains
its accepted diagnostic-acquisition CLOSED / RUNTIME VALIDATED scope; there is
no new farm or physical-subtraction acceptance. Previous review-pending wording
above remains historical and is superseded for current status here.

SOURCE VERIFIED: pion_component_shapes.py normalizes selected SIMC shapes.
pion_component_fits.py::_run_staged_component_pass builds component residuals
from the same target and earlier fitted amplitudes; the protected Delta path
constructs R_preDelta from that kaon target minus early fitted terms and calls
_fit_lambda_gauge on it. The original Center/high t3 receipt confirms gauge
source_histogram=R_preDelta, window [1.105,1.125], one-component weighted least
squares and available conditional fit uncertainty. Existing coefficient
covariance is present; it is not an independent anchor/gauge calibration.

The original acquisition source packet explicitly labels covariance conditional
fixed-weight Sumw2 with fitted/training/shared-systematic covariance unavailable;
it contains 6527 processed kaon and 24967 pion rows for Center/high. Thus observed
data are present, not globally missing. Its applied kaon amplitudes are
(n,Delta,SIDIS)=(0.1822022441,0.5690006384,0), while original control SIDIS is
2.9845068988 and detached joint central SIDIS is 3.2082020108. These are fitted
conditional allocations, not physical truth or an independent zero-SIDIS bound.

SOURCE VERIFIED: bootstrap_pion_mm_delta.py imports the same protected A values
and holds their estimates fixed, sharing resampled MC shapes with the fitted
control denominator. propagate_pion_mm_delta_uncertainty.py consumes the retained
replicas/model scenarios and explicitly preserves missing full physical/PID
covariance. This provides useful conditional statistical covariance and model
sensitivity; neither estimates the missing independently calibrated kaon-side
means. The three current supplied review packages and original acquisition
record provide no qualifying independent joint-constraint/calibration record.
This is a scoped admissibility finding about inspected evidence, not a claim
that no such measurement can exist elsewhere. No historical archive search,
zero-photoelectron study, synthetic experiment, source change or refit occurred.

Outcome B / BLOCKED remains. First missing input is a jointly specified pair
of independent aggregate component means for Center/high t3: non-SIDIS mean in
the kaon SIDIS anchor [1.0601751288600272,1.1101751288600272] GeV and pion-background
mean entering the existing Lambda gauge [1.105,1.125] GeV. Supply their 2x2 mean
covariance and cross covariance to evaluated signed prompt/random/dummy/
dummy-random observations, with exact component basis, PID/source normalization,
cuts, coordinate/native-bin integration convention and training/calibration
identities. Overlapping intervals require joint treatment. If disjointness
justifies zero cross covariance, document it; do not assume independence.
An independently calibrated aggregate constraint suffices without eventwise
species labels. Same-data fits, observed total gauge counts and copied control
fractions are not substitutes. External PID uncertainty remains separately
owned and unavailable, not a supplied independent constraint by default.

Stop at the independent-input gate; no implementation contract/source change
is warranted without that input. Independent review of qualifying new evidence
precedes the existing blind predictive closure. All six roadmap dependencies,
Local-t selection, no_empirical_residual, Method A detached/B diagnostic-excluded,
frozen normalization/production/yield logic, unrelated bytes, index and refs
are preserved. The review ZIP contains the four-file memory diff/snapshots,
unchanged source snapshots, actual receipt/result extracts, input hashes and
memory-integrity/preservation checks. No new scientific covariance/mean or
combined one-sigma physical error was invented. No tests of changed analysis
code were needed because no analysis code changed; farm/runtime NOT RUN.
