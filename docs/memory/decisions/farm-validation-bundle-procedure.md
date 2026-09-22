# Farm validation-bundle procedure

This is the durable workflow for packaging KaonLT validation evidence on the
Jefferson Lab farm. It is procedural guidance, not scientific evidence and not
a phase-status decision. A successfully created ZIP is not accepted evidence
until its manifest/provenance and the requested scientific or presentation
evidence have been reviewed separately.

This procedural update was prepared at observed `test` HEAD
`69629c90929eea95ab7bfa465d44cf423521fcab`. It documents the exercised
F.6.2.Fix.4 candidate procedure only; it does not accept its candidate PDF or
ZIP, or make a phase-status decision.

## Farm shell and communication

The Jefferson Lab farm shell is `tcsh`. Write farm commands in valid `tcsh`,
not Bash. For a simple operation, communicate concisely: state the exact
artifact path or paths, state any narrow caveat, provide one exact `tcsh`
command block, then ask for the full command output on failure or the resulting
artifact/ZIP on success. Do not hide a simple bundle operation behind wrapper
scripts, large shell programs, Bash heredocs, or guessed alternate paths.

The canonical farm paths are:

- analysis artifacts: `/group/c-kaonlt/USERS/trottar/lt_analysis/OUTPUT/Analysis/KaonLT`
- transferable bundles: `/volatile/hallc/c-kaonlt/trottar/globus`

## Checkout state and provenance

An ordinary farm checkout may intentionally contain modified or untracked,
unpushed files unrelated to a detached validation. Known examples include
`src/models/xmodel_kaon_pl.f`, `src/models/xmodel_pion_pl.f`,
`src/kaon/functions/Q4p4W2p74.model`,
`src/kaon/xsects/unsep_Q44W274.csv`, and a local helper script. Do not require
that ordinary checkout to be globally clean before detached bundle creation;
do not delete, reset, stash, checkout-over, or otherwise disturb those files.
Record relevant local state in provenance, but treat it as contamination only
when it is actually on the execution or artifact-consumption path.

This allowance never weakens source-provenance checks. When the normal checkout
has advanced beyond the reviewed collector/profile state and later committed
files would be rejected, use a temporary detached Git worktree at the exact
reviewed bundle/profile commit. Run the collector and profile from that clean
worktree, point `--outdir` to the existing canonical artifact directory, write
the ZIP to the canonical Globus directory, then remove and prune the temporary
worktree. Never broaden a profile allowlist merely to admit a later unrelated
commit, and never modify the ordinary checkout to create a bundle.

## Determine the requested farm operation first

Before issuing a farm command, distinguish a full scientific rerun, a
presentation-only rerender plus bundle, and a bundle-only operation. Do not
default to the full analyzer. A presentation-source repair must not regenerate
an accepted frozen scientific payload.

For a full scientific rerun, follow that phase's scientific contract and run
its analyzer before collection. This procedure does not authorize one.

For bundle-only work, when the required JSON and PDF are already the intended
versions, run only the validation-bundle collector. Do not rerun analysis or
delete/regenerate the frozen JSON. Use a unique timestamped ZIP name, inspect
its integrity and manifest, and transfer only that exact fresh ZIP.

For a presentation-only rerender plus bundle, treat the accepted JSON as an
immutable input. Verify its accepted SHA-256 before rendering; load that JSON;
call only the phase presentation renderer directly; write a temporary PDF in
the artifact directory; require that new PDF not to have the old known PDF
identity; atomically replace only the PDF; and recheck that the JSON SHA-256 is
unchanged. Do not invoke analyzer `main()` for this operation: it rebuilds the
scientific artifact and is not the requested step. Then collect and inspect a
fresh bundle.

The collector intentionally refuses to overwrite an existing ZIP. Use a unique
timestamped ZIP name for repeatable attempts. If a phase contract instead
requires a deterministic name, remove its old ZIP only when that removal is
explicitly intended and safe. Do not create multiple user-visible helper files
with the same name.

For a temporary-worktree collection, use the phase contract's exact collector
arguments. Repeated **bundle-only** gates may instead use the transparent,
input-driven `testing/package_pion_hgcer_validation_bundle.tcsh` wrapper under
its [wrapper contract](farm-validation-bundle-wrapper-contract.md). It accepts
only explicit bundle inputs, invokes the unchanged generic collector, and
cannot run an analyzer, renderer, scheduler, or arbitrary command. Full
scientific reruns and presentation-only rerenders remain phase-specific direct
operations; do not generalize either into the wrapper.

## F.6.2.Fix.4 example only

The following identities are specific to the current F.6.2.Fix.4 packaging
workflow and must not be generalized to later phases:

- reviewed clean bundle/profile commit:
  `c3252f1ec43dcca0357aebe8c6ca75d7d198f99a`
- profile's reviewed Fix.4 presentation source:
  `b929815517643edaa75949f7e61ce46c6e8f0d63`
- accepted F.6.2 scientific JSON SHA-256:
  `5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1`
- prior accepted F.6.2 artifact fingerprint:
  `ee713b70de898bad8fa61164cbc8a54886af1ea1eb1e712ed95de17df8890cd0`
- prior accepted F.6.2 validation fingerprint:
  `7edc73fce20ad7dc7622b8c23a7ba7e8986e367ace605884e010945595370b3b`
- old pre-Fix.4 review-PDF SHA-256:
  `4d07c2989c217d18d088b181390e9478fae51e4168dba2d64faea690b40f0f85`

If the accepted JSON and intended PDF already exist, this is bundle-only: run
the collector from the reviewed clean worktree and inspect the resulting
manifest before transfer. Do not rerun the F.6.2 analyzer or regenerate the
accepted JSON.

If the accepted JSON is frozen but the PDF still has the old identity, use this
presentation-only rerender plus bundle command. It runs the renderer from a
temporary detached worktree at the clean collector/profile commit (which
contains the reviewed presentation source), points only `--outdir` at the
canonical artifact directory, and never changes the ordinary checkout.

```tcsh
cd /group/c-kaonlt/USERS/trottar/lt_analysis
set py = /u/group/c-kaonlt/USERS/trottar/replay_lt_env/bin/python
setenv PYTHONDONTWRITEBYTECODE 1
set artifact_dir = /group/c-kaonlt/USERS/trottar/lt_analysis/OUTPUT/Analysis/KaonLT
set expected_json = 5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1
set old_pdf = 4d07c2989c217d18d088b181390e9478fae51e4168dba2d64faea690b40f0f85
set timestamp = `date -u +%Y%m%d-%H%M%S`
set worktree = "/tmp/kaonlt-f6-2-fix4-$timestamp-$$"
set zip = "/volatile/hallc/c-kaonlt/trottar/globus/KaonLT_PhaseF6_2_Fix4_validation_Q4p4W2p74_$timestamp.zip"
git worktree add --detach "$worktree" c3252f1ec43dcca0357aebe8c6ca75d7d198f99a
if ($status != 0) exit 1
$py -c 'import atexit,hashlib,json,os,sys,tempfile; from pathlib import Path; repo,artifact_dir=Path(sys.argv[1]),Path(sys.argv[2]); expected_json,old_pdf=sys.argv[3],sys.argv[4]; sys.path.insert(0,str(repo/"testing")); import analyze_pion_hgcer_method_a_acceptance_refinement_validation as f6_2; json_path=artifact_dir/"Q4p4W2p74_kaon_pion-background_hgcer_method-a-acceptance-refinement-validation.json"; pdf_path=artifact_dir/"Q4p4W2p74_kaon_pion-background_hgcer_method-a-acceptance-refinement-validation.pdf"; sha256=lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest(); check=lambda ok,label: ok or (_ for _ in ()).throw(RuntimeError(label)); check(json_path.is_file(),"f6_2_json_missing"); check(pdf_path.is_file(),"f6_2_pdf_missing"); check(sha256(json_path)==expected_json,"f6_2_json_sha_mismatch"); check(sha256(pdf_path)==old_pdf,"f6_2_existing_pdf_not_old_identity"); artifact=json.loads(json_path.read_text(encoding="utf-8")); fd,name=tempfile.mkstemp(prefix=".f6_2_fix4_",suffix=".pdf",dir=artifact_dir); os.close(fd); temporary_pdf=Path(name); atexit.register(lambda: temporary_pdf.exists() and temporary_pdf.unlink()); f6_2.write_review_pdf(temporary_pdf,artifact); check(sha256(temporary_pdf)!=old_pdf,"f6_2_pdf_sha_unchanged"); check(sha256(json_path)==expected_json,"f6_2_json_changed_before_replace"); os.replace(temporary_pdf,pdf_path); check(sha256(json_path)==expected_json,"f6_2_json_changed_after_replace"); print("candidate PDF SHA-256: "+sha256(pdf_path))' "$worktree" "$artifact_dir" "$expected_json" "$old_pdf"
if ($status != 0) then
  git worktree remove --force "$worktree"
  git worktree prune
  exit 1
endif
$py "$worktree/testing/collect_pion_hgcer_validation_bundle.py" --outdir "$artifact_dir" --kinematic Q4p4W2p74 --profile "$worktree/testing/pion_hgcer_validation_bundle_profile_f6_2.json" --output "$zip"
if ($status != 0) then
  git worktree remove --force "$worktree"
  git worktree prune
  exit 1
endif
unzip -t "$zip"
if ($status != 0) then
  git worktree remove --force "$worktree"
  git worktree prune
  exit 1
endif
unzip -p "$zip" manifest.json | $py -m json.tool
set manifest_status = $status
git worktree remove --force "$worktree"
git worktree prune
if ($manifest_status != 0) exit 1
echo "candidate ZIP: $zip"
```

The observed 2026-09-18 outputs from this procedure are candidate / pending
evidence review only: JSON SHA-256 remained
`5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1`; the
candidate Fix.4 PDF SHA-256 is
`a014ba89d36a38d1387448b5176b17e932f1257ab0c8fa6b605271009758358a`; and
the candidate bundle
`KaonLT_PhaseF6_2_Fix4_validation_Q4p4W2p74_20260918-002336.zip` has SHA-256
`2118c5e360a701295eda45c02aa76d6cfa0774c6e95f431d6f369607e8f854f7`.
ZIP creation and manifest completion do not accept the PDF or close Fix.4;
request the exact ZIP for independent ChatGPT evidence review.

## Non-negotiable boundaries

Never clean, reset, stash, or otherwise disturb the ordinary farm checkout
without explicit user direction; delete local model files to make status clean;
weaken collector provenance gates; rerun expensive analysis for a bundle-only
request; regenerate a frozen artifact during a presentation-only or bundle-only
step; or claim acceptance from ZIP creation alone.
