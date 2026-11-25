# KaonLT Step 2 — isolated genuine K⁺Σ⁰ generation

## Authority and start

Start from `test` at `e98889fb7aadb54005c34a0b3c02d81f05b4f8ff`. Verify the live branch, HEAD, index, full worktree and repository startup core before editing. A moved HEAD needs an independently reviewed re-anchor. Preserve existing unstaged memory, unrelated work, the LF repair, and the parked external-diamond contract. This is the separately authorized Step 2 infrastructure task. The Q4p4 offline diamond comparison remains blocked on its existing pre-diamond data/dummy/KΛ evidence and is outside this contract.

Use GPT-6.2 Sol Medium at Standard speed for shell, configuration and input handling. Fallback GPT-6.1 Sol Medium. Consider High only if a verified Fortran, detector or reconstruction incompatibility seems to require a separately scoped scientific-source repair.

Use the externally held `KaonLT_Colleague_KSigma0_AllSettings_20261123T231712Z-2526928.zip` with SHA-256 `bc68d67dd1e2501492827f9a855347a3ca006f28a149522454fece6edba7af26`. The companion `generation_inputs_manifest.json` (SHA-256 `4b795a56666b96986f135f54917e4f4c377f76740204d0afa5c016140c15e5ef`) and `settings25.csv` pin the exact input members, settings and historical report identities. Copy or generate only allowlisted provenance metadata in the public repository. Keep the archive, its Fortran/C++/table/input bytes, original script files, and private physical paths out of tracked material. Do not import or execute the supplied launchers or submission scripts.

## Allowlist and objective

Implement one opt-in, isolated generation path for the verified 25 Q2p1_W2p95, Q3p0_W2p32, Q3p0_W3p14, Q4p4_W2p74 and Q5p5_W3p02 input combinations. Reuse the existing background launcher and generation → raw ROOT conversion → reconstruction runner. One branch in the launcher dispatches the Sigma-only route before legacy symlink changes.

Allow only:

- `background_samples/set_BackgroundInput.sh` and `background_samples/background_samples.conf`, for explicit dispatch and empty/default-off external configuration.
- `background_samples/ksigma0_inputs_manifest.json`, reviewed metadata only, and `background_samples/prepare_ksigma0.py`, for verified extraction, validation and dedicated workspace preparation. Use the supplied manifest as the sole 25-setting lookup.
- `src/setup/run_simc_tree`, for the opt-in isolated adapter while preserving the existing caller contract.
- `src/setup/validate_ksigma0_output.py`, for generated-output verification and the detached consumer-ready copy.
- `testing/test_ksigma0_generation.py`, for deterministic fixtures.
- `docs/memory/CURRENT.md`, the existing event-population-authority phase record and the manifest, only for the warranted source-review closeout.

All other tracked files are frozen. In particular, preserve all current generation channels, source resolvers, `background_sample_manifest.py`, the `LT_BG_SIGMA0_*_ROOT` interface, component loaders/fits, normalization, cuts and background ownership. Do not edit the local reference SIMC repository or the colleague's tree. If an actual model or reconstruction change is required, stop on the earliest proved incompatibility and preserve evidence for separate review.

## Source and input admission

1. Use an explicit Sigma-only selector and one exact setting or one explicit valid batch of settings. Unset selector preserves the old mode and arguments. Reject unknown kinematic/ε combinations, Right/low, ambiguous mixes with ordinary channels, duplicates and malformed values before writes or launches. The Sigma branch runs no neutron, Delta or SIDIS work. Default background configuration remains unchanged. Document the new local CLI in launcher help once implemented and reviewed.
2. Require an approved, externally configured dedicated base identity (commit plus tracked source and deployed build/runtime manifest). Compare the available reference tree at `641b8c07182fd045fd8a5c007cdc1b71a418f3ce` only for interface evidence. Its local `physics_iterate.f` is modified and must not be incorporated from that worktree. The supplied archive does not contain the complete historical `simc_gfortran_corrected` tree. Do not guess it or use an arbitrary substitute. Verify the proposed base builds `physics_iterate.o`, provides the compatible `peek` model link and Fortran module/include interface, uses the documented `which_kaon=1` proton-to-Σ⁰ recoil, and provides the matching drift chamber, database, binary writer, ROOT converter and runtime libraries. A missing tree or facility toolchain is a runtime prerequisite, not a reason to invent source.
3. Admit exact archive bytes by outer SHA-256, then use unique safe member names, bounded verified reads, per-member SHA-256 and rejection of symlink, traversal and duplicate entries. Stage setting bytes unchanged at `input/<canonical-stem>.inp` where the existing naming form is `Prod_Coin_Q…W…<arm>_<lowe|highe>`. Only the basename changes. Preserve reaction flags, requested successes, decay, radiation, energy loss and smearing. Do not reuse historical logs, ROOT-like names or generated metadata as outputs of a new attempt.
4. Stage the exact `physics_iterate.f`, `par.pl`, `shms/mc_shms_hut.f` from the archive's matching `mc_shms_hut.f`, and the exact reconstruction pair under the private reconstruction directory. Preserve the archived Saghai and weight tables actually read by the approved base's startup at their expected relative locations. Verify the approved base's actual active data paths. In the local reference, `dbase_read` opens `weightc.dat`, `weightd.dat` and both Saghai files for kaon startup. The supplied `peek` fit uses `par.pl`; those startup reads do not make Saghai the active cross-section model. The supplied `default.dat` is identical to `par.pl` and is provenance unless an approved base documents a separate read. Preserve and audit supplied parameter and detector bytes. No refit, default model override, coefficient change, mass change, re-tuning or production parameter helper.

## Isolation and runner contract

Use explicit read-only external archive and base roles plus a user-owned private attempt root. Create one fresh exclusive working root per setting attempt. Resolve every effective path before use. Require build/worksim/input/OUTPUTS, reconstruction, logs and temporary outputs to stay inside that root. Fit parameters and all Fortran-relative table reads must resolve there. Stage independent files, without mutable hardlinks or symlinks back to protected paths. Detect overlaps, existing outputs, escaping link targets, stale outputs, collisions or owner mismatch before writes. Keep existing general runner invocations behaving as before, including their existing argument defaults. The Sigma branch must never call a shared rewiring or a production setup helper.

Use the approved native generator and converter in the dedicated tree through `run_simc_tree` with the existing resolved runtime/log roles. Verify private executable and compiler/ROOT/environment identity and build success. Preserve the native `input/<stem>.inp` lookup, `worksim/<stem>.bin` conversion and production of `.hist`, `.gen` and `.geni`. Reject early failures, stale files or incomplete sidecars. Archive source model copy, generator metadata and logs by digest.

For the audited macro use the exact reconstruction call shape `recon_hcana.C("<stem>", "kaon", 0.493677)` with the original generated `<stem>.root` and `<stem>.hist` co-located in its private working directory. It produces a separate `<stem>_recon.root`. The old two-argument `production` invocation and implicit file routing are incompatible. Use the exact macro and header, and keep the raw and reconstructed products distinct. No global macro replacement. No shell-evaluated arbitrary paths, stem or reaction strings.

The macro already copies `Weight` via its common branch writer. Prove it on the actual surviving entries in order, with independently checked membership for the archived tray-cut logic or equivalent reliable entry identity. Require ROOT `h10`, typed finite branches needed by the current loader (Q², W, t, phipq, missmass, HMS/SHMS coordinates and HGCer axes) and finite event weights. Verify the coordinate definitions and support, including GeV²/GeV, signed t versus −t, radians, and missing mass with the fixed kaon hypothesis. Validate recoil identity and MM at the first runtime gate. Keep native Weight unchanged. If producing a consumer-ready artifact, make a separate clone with `iter_weight` equal to the preserved per-event `Weight`, with identical surviving count/order and checked value and Sumw2 equality. Do not invoke a production reweight model or rely on the unchanged consumer's unit fallback. Validate this as a lossless generation adapter before later binding. No automatic export of `LT_BG_SIGMA0_*_ROOT` or import into working subtraction configurations.

## Provenance, tests and fail-closed conditions

- Positive coverage: all 25 exact member hashes and output stem mappings, both ε and all five physical settings per kinematic, exact kaon/Σ⁰ flags, and approved base provenance. Keep per-setting historical report values clearly identified as historical. Include build commands/versions, loaded parameter/model/table hashes, normalized roots, output and sidecar digests, reconstructed row counts/membership and weighting receipts for any future actual attempt.
- Negative coverage: wrong archive/base hash, wrong reaction, missing dependency, altered table/parameter/member, unauthorized file inclusion, hostile archive paths, unknown setting, Right/low, no opt-in, mixed ordinary and Sigma jobs, ambiguous symlinks, collisions, path escapes, existing outputs, mid-run failure, missing or incompatible artifacts, and export of consumer bindings. Assertions must occur before side effects whenever possible. An interrupted run can be diagnosed with a failed receipt and partial logs. Only checked artifacts may become accepted candidates.
- Regression: demonstrate unchanged default command routing, links/paths, inputs, expected runner argument shape and outputs for neutron, Delta and SIDIS using disposable deterministic stubs. Preserve working production invocations and test the existing Sigma external-source binding tests. Exercise shell argument validation and return codes, path and archive security checks, deterministic runner sequencing, sidecar parser, Sigma normalization fields, raw-to-reconstructed Weight preservation and exact identity-copy branches. Use mock/stub tests for ROOT behaviors that cannot be proven locally, and label those limits.
- The unchanged baseline has a WSL cleanup-guard regression for a symlink into the checkout with a missing intermediate directory. Keep the baseline result visible when comparing regressions. The dedicated Sigma route must bypass this cleanup code and prove strict private path handling independently. If it depends on the failing guard or reaches shared setup, stop. Repairing the legacy production guard is outside this contract.
- Run scoped shell syntax and Python tests, verify no external or installed tree mutation, inspect the complete actual tracked diff and all intended additions, verify the memory manifest and run ordinary non-strict memory health. Report CURRENT/MEMORY/CURRENT_HANDOFF bytes and any warnings. Do not equate passing tests with actual ROOT or farm behavior.

## November 24 user-authorized actual-diff repair

The independent review held publication. The follow-up local repair is limited
to the three identified defects. Add only a path-specific root .gitattributes
entry to preserve the exact reviewed JSON bytes in a fresh checkout. The
isolated runner must re-admit the independent external approval, physical
base/archive/private roles and staged inventory before any build.

The archived reconstruction pair remains hash-pinned. The exact original macro
is retained under private provenance. The sole permitted reconstruction
derivation is replacing the pinned macro destructor with a defaulted destructor:
Process owns the two trees through its TFiles, and the object must not delete
the dangling tree pointers. Record both hashes and keep the original invocation,
header, reaction, kinematics, cuts, normalization and branch handling. Reject
different source bytes or a different destructor. This scoped amendment
supersedes the byte-identical executed macro requirement only for that one
ownership edit. No other scientific or reconstruction behavior may change.

Require a fresh-checkout manifest test, independent runner admission tests and
an exact-byte ownership-derivation test. ROOT lifetime execution and scientific
acceptance remain unverified pending the existing review/publication sequence.

## November 24 user-authorized reference compatibility amendment

The November 24 user follow-up authorizes only the concrete base-admission and
ROOT identity compatibility repair demonstrated by reference commit
`641b8c07182fd045fd8a5c007cdc1b71a418f3ce`. Preserve all three accepted repairs
above, including the sole ownership derivation of the archived macro.

Admit the reference's exact Git-identified immutable reconstruction C/header and
launcher metadata, plus `cosy/COSY.bin` and `cosy/JLAB.bin`; these are tracked
dependencies, not generated job products. Do not execute or stage the base
reconstruction launcher for execution. The archived Sigma0 pair remains the
private reconstruction overlay. Only `shms/shms_forward.dat` and
`shms/shms_recon.dat` may be admitted as tracked relative optics symlinks, with
their exact Git link and target identities. Verify the on-disk link spelling,
regular tracked targets and containment; materialize independent regular copies
of the approved target bytes at the private alias paths. No other symlink,
untracked file, arbitrary binary output or escaping target is authorized.

Inspect the actual converter Makefile's effective ROOTSYS, ROOTDIR and absolute
root-config calls. A private command-line ROOTSYS selection may override the
hardcoded installation without changing the Makefile. Require the effective
converter configuration, reconstruction ROOT and validator PyROOT to identify
the same approved installation/version and loaded core library. Run a lightweight
PyROOT import/identity probe in the validator interpreter/environment before any
generation; missing or mismatched ROOT/PyROOT blocks the attempt. Record the
identities in the private receipt. No facility software installation or shared
environment/source mutation is authorized.

Limit edits to this amendment, the existing preparer/output validator and tests,
and the existing allowlisted memory closeout records. Regression coverage must
use the actual reference Git tree inventory and source-defined optics/converter
paths, with strict negative cases. Reference Git objects may be read into
disposable local fixtures; never import the modified reference worktree model.
Return one complete actual-diff bundle with all intended additions and test/
memory receipts. Independent review still precedes publication and runtime;
all original scientific, diamond and user-only Git/farm restrictions remain.

## November 24 user-authorized runtime-identity closeout

Independent review verified the compatibility diff and held publication on two
runtime identity gaps. This follow-up is limited to the existing preparer,
validator, deterministic fixtures and the corresponding memory/contract records.

The pre-generation probe in the private attempt and final event validation must
use the same explicit Python executable, entry point, private working directory
and approved environment, with user site disabled. Retain the early
external-base probe. Require identity equality in the private probe before any
build and in the final validator before opening event files. Preserve the
validator's existing numerical checks and result schema.

Require the GNU/Linux loader inspector `ldd` as an independently pinned tool
in the external approval manifest. Resolve the actual built native converter
under its launch environment, without running the converter's event input.
Require one resolved libCore path and digest exactly matching the admitted
ROOT/PyROOT runtime, with no missing dependencies. Record the loader output and
resolved identity, and stop before native conversion if the check fails.

Use deterministic tests for actual dispatch context and both failure orders.
The previous archive/model/optics and reconstruction ownership decisions remain
in scope as established. ROOT/farm execution and publication remain separate.

## Subsequent runtime gate and hard stop

After independent actual-diff review, user publication and pushed-state verification, the smallest runtime gate is one user-controlled Q4p4W2p74 Center/high attempt at its archived settings. Verify dedicated base/model identity and physical paths before compilation, then actual generator, native converter, macro and output-checker execution in the same owner attempt. Require histogram reaction/normalization fields including current Ngen, Ncontribute, normfac, and distinguish generated successes and converted or post-tray-cut entries. Require finite signed Weight and Sumw2 checks, event/branch identity, genuine missing-mass and accepted kinematic support. Any Sig0 diagnostic using the existing acceptance/cuts remains separately scoped. After independent scientific review expand first to the four other Q4p4 settings, then the other twenty when authorized. No synthetic MM plots.

Keep the unrelated diamond evidence gate separate. Preserve Local-t, W03/F120, L26, separate F280, D1.a OFF, the 67/135 mask, Stage 1–3 and all six milestones. Stage 4 and D3–D5 remain deferred, Outcome B blocked, Method A detached and Method B diagnostic-only/excluded. KΣ⁰ generation gives no authority to fit or adopt a template.

Stop after the narrowly allowed implementation, deterministic local checks and actual-diff handoff. No farm execution, Git publication, use of colleague directories, diamond acquisition or background refitting. No unreviewed farm command.
