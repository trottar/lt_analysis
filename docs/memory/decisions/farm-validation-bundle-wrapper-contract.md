# Parameterized farm validation-bundle wrapper contract

## Objective and exact starting HEAD

At observed `test` HEAD `45633698d624304005afdf7a9db25cd4d67ec573`, add one
input-driven `tcsh` wrapper for repeated **bundle-only** KaonLT validation
collection. The wrapper reduces hand-written farm command blocks without
changing the generic collector, any profile, analysis source, or runtime gate.
It is local source work only until the user reviews, commits, pushes, and then
runs a fresh farm bundle.

## Allowed and frozen files

Allowed source files are:

- `testing/package_pion_hgcer_validation_bundle.tcsh`
- `testing/test_package_pion_hgcer_validation_bundle_tcsh.py`

Allowed memory records are this contract, the farm-validation procedure,
decision history, the wrapper phase record, `CURRENT.md`, `USER.md`, and the
generated manifest. All analysis source, the generic collector, every
validation profile, the ordinary farm checkout, frozen artifacts, production
weights, yields, cuts, normalizations, templates, and binning are frozen.

## Before and after behavior

Before this change, each bundle-only gate repeats a phase-specific detached
worktree/collector command block. Afterwards, a caller supplies the exact
Python executable, full bundle commit, repository-relative profile,
artifact directory, kinematic token, output ZIP, optional authorized setting,
and zero or more immutable `PATH SHA256` assertions. The wrapper creates a
detached worktree at the supplied full commit, invokes the unchanged generic
collector, rechecks every supplied immutable input, verifies the new ZIP and
its manifest, and removes only the created temporary worktree.

## Preserved boundaries and forbidden shortcuts

The wrapper packages existing profile-declared evidence only. It must not run
an analyzer, renderer, scheduler, or arbitrary caller-supplied command; accept
an abbreviated/unknown commit; overwrite a ZIP; weaken profile provenance;
regenerate a frozen artifact; clean, reset, stash, or otherwise disturb the
ordinary checkout; or use administrator, scheduler, or shared-filesystem
actions. It must be valid `tcsh` and use only the caller's designated paths.

Full scientific reruns and presentation-only rerenders remain separate,
phase-owned operations. Their source invocation is never parameterized by this
wrapper; only their already-created artifacts can be packaged after the
appropriate source/runtime procedure has completed.

## Checks, farm boundary, and hard stop

The local test statically verifies required inputs, generic-collector use,
detached-worktree confinement, immutable hashing, ZIP inspection, and the
absence of analysis/scheduler/ordinary-checkout mutation paths. It runs
`tcsh -n` only when a local `tcsh` executable exists. Local checks do not
establish farm behavior.

Codex performs only local editing and checks; the user alone commits, pushes,
and executes farm work. A fresh user-supplied ZIP remains required for evidence
review. Stop rather than add a generic analysis runner, a shell `eval`, a
fallback path, an implicit input, or a bypass for a failed collector/hash/ZIP
check.
