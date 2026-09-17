# Farm validation-bundle procedure

This is the durable workflow for packaging KaonLT validation evidence on the
Jefferson Lab farm. It is procedural guidance, not scientific evidence and not
a phase-status decision. A successfully created ZIP is not accepted evidence
until its manifest/provenance and the requested scientific or presentation
evidence have been reviewed separately.

This record was updated from the observed `test` HEAD
`31be08bcd8956b4ee3593b0919ed2fcdf5101e0d`. The supplied starting identity
was `ab751c82f66eb1f16d0f9b4dbc6e7c3f6b323ea7`; the intervening changes only
removed a farm helper file and did not affect this procedure.

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

## Match the command to the requested step

Determine first whether analysis or a presentation rerender has already run.
When the required JSON/PDF artifacts exist and the request is bundle-only, run
only the validation-bundle collector: do not rerun the analyzer and do not
delete or regenerate a frozen scientific JSON. If a fresh analysis or rerender
is actually required for that phase, run that phase's analyzer first and its
collector second.

The collector intentionally refuses to overwrite an existing ZIP. Use a unique
timestamped ZIP name for repeatable attempts. If a phase contract instead
requires a deterministic name, remove its old ZIP only when that removal is
explicitly intended and safe. Do not create multiple user-visible helper files
with the same name.

For a temporary-worktree collection, use the phase contract's exact collector
arguments, with this concise `tcsh` sequence: enter the normal repository;
define timestamp, temporary-worktree, artifact-directory, and unique-ZIP
variables; add the detached worktree at the reviewed bundle/profile commit;
run its collector/profile against the existing artifacts; inspect the ZIP
manifest; remove and prune the worktree; and print the exact ZIP path for
transfer. Keep the phase-specific invocation in one short command block rather
than generalizing it into a new script.

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

For this example, if the approved JSON/PDF already exist and the requested
operation is packaging, do not rerun the F.6.2 analyzer or regenerate the
accepted JSON. Run the collector from the reviewed clean worktree and inspect
the resulting manifest before transfer.

## Non-negotiable boundaries

Never clean, reset, stash, or otherwise disturb the ordinary farm checkout
without explicit user direction; delete local model files to make status clean;
weaken collector provenance gates; rerun expensive analysis for a bundle-only
request; regenerate a frozen artifact during a presentation-only or bundle-only
step; or claim acceptance from ZIP creation alone.
