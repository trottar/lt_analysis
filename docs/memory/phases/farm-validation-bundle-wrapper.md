# Parameterized farm validation-bundle wrapper

## Scope

At observed `test` HEAD `45633698d624304005afdf7a9db25cd4d67ec573`, the
user authorized an input-driven reusable farm bundle wrapper. It is
`SOURCE REVIEWED` local/proposed work only until the user reviews, commits,
and pushes it. It does not alter E.8.1's procedure-PDF source/profile identity,
the frozen F.6.2 JSON, or any farm/runtime status.

## Implementation

`testing/package_pion_hgcer_validation_bundle.tcsh` accepts explicit Python,
full bundle commit, repository-relative profile, canonical artifact directory,
kinematic, ZIP, optional setting, and immutable-file SHA-256 assertions. It
uses a detached worktree at the supplied commit, invokes the unchanged generic
collector, verifies immutable inputs before and after collection, checks the
new ZIP/manifest, and deletes only its created `/tmp` worktree.

It is deliberately bundle-only: no analyzer, renderer, scheduler, arbitrary
caller command, ZIP overwrite, `git clean`, reset, stash, or ordinary-checkout
mutation path is available. The paired test statically enforces those
boundaries and parses with `tcsh -n` when that interpreter is available.

## Local review and runtime boundary

Focused wrapper tests and the E.8.1 profile suite passed locally. The local
host has no `tcsh`, so syntax/runtime execution remains farm-only and has not
been claimed. A fresh user-supplied farm ZIP is still required before any
bundle or E.8.1 runtime conclusion.

## Remaining action

The user reviews, commits, and pushes this source change. Then use the wrapper
only after the phase-owned source rerun or presentation rerender has produced
the intended artifacts; it never substitutes for that operation.
