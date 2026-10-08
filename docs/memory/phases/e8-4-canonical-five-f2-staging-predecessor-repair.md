# Canonical-five F.2 staging-predecessor migration repair

## Scope and state

2026-10-07: `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING` under the
[task contract](e8-4-canonical-five-f2-staging-predecessor-repair-task-contract.md).
Branch `test`, HEAD/local `origin/test`
`946b5512bb00d22bb446bd96d6ab9ee35f5f989b` passed the exact start gate.
The prior F.2/F.3/F.4 equivalence repair remains source-reviewed/pushed per
the supplied contract and `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`.

The [supplied failed-attempt diagnosis](../evidence/e8-4-canonical-five-f2-staging-predecessor-blocker-2026-10-07.md)
owns the fresh candidate_staging failure, analysis_started=false, observed
identities and source-owned exact scientific comparison. No local farm
artifact inspection or new farm acceptance is claimed.

## SOURCE VERIFIED — one owner staging exception

Only the canonical-five owner and its tests change executable behavior.
`RECOGNIZED_STAGING_PREDECESSORS` contains exactly the explicit F.2 basename:

`Q4p4W2p74_kaon_pion-background_hgcer-method-a-acceptance-representation-current-baseline-candidate.json`

mapped to the observed target SHA-256:
`182433ccd2d13d3f7af80b963e69190ad8f2eb4e295c8216b36585e8ee7409b2`.

This state may only be replaced by reviewed candidate SHA-256:
`2fa715b2b5c2d5077e38416fb9eceb814e013e72f308c43526fddff36ddd962e`.

The predecessor is absent from materialization, candidate, historical-input,
Left/lowe and runtime-authority pins. It is not accepted scientific authority
or a fallback. No scientific comparison is performed during staging.
Every other unknown F.2 target still raises `unknown_candidate_target:<basename>`
before replacement. Existing `accepted.CANDIDATES.get(name)` behavior for
F.3/F.4 is unchanged.

Source hashing, temporary copy/hash, destination recheck, atomic replacement
and final reviewed hash enforcement are unchanged, as are symlink/non-file/
source-target-alias checks. Installation records retain path, before_sha256,
after_sha256 and action. The recognized transition records `replaced`.
Isolation, ordinary-checkout/ltsep preservation, lineage, freshness, artifacts,
page, collection, ZIP and final delivery gates are unchanged.

## SOURCE VERIFIED — deterministic checks and limits

Six new focused tests passed. Existing tests retain absent/current F.2 and
recognized F.3/F.4 predecessor coverage. New checks cover the exact singleton
policy, successful predecessor replacement, unrelated/nearby hashes, wrong
source bytes, concurrent target change, symlink, non-file and source alias.

Observed predecessor bytes are not available locally. The positive test
explicitly mocks only the pre-existing F.2 planning/recheck hashes. Reviewed
source fixtures use their actual synthetic hashes via the existing fixture-pin
mechanism; source, temporary-copy and final hashes remain real. Final bytes
equal source bytes, status carries the before/after transition, and the
source materialization is unchanged. No bytes are invented for the farm hash.

All 119 tests passed, zero failures/errors/skips, across the five required
modules plus the complete frozen Left/lowe owner suite. Linux Python 3.12.3
used isolated dependencies outside the repository. The initial combined run
had four failures confined to the historical Left/lowe synthetic fixture's
freshness condition. A direct clock probe measured new /tmp file mtime
4,364,612 ns earlier than the pre-write time sample.

The final test process controlled mtimes only for `Path.write_bytes` calls
from the historical test module's synthetic `run` function, setting them to
time.time_ns()+1 second. Its deliberate stale cases subsequently set mtime
to 1 and still exercised strict rejection. All 35 affected writes were
test-owned. The temporary runner did not modify repository source or tests;
production freshness semantics remained unchanged. The required modules and
six new tests did not require this historical fixture control.

Manifest write/check, ordinary memory health, bootstrap, diff check and final
byte/index/ref/allowlist preservation checks accompany the complete review
bundle. CURRENT's added diagnostic links exceed the 8 KiB soft limit only;
classify that warning as nonblocking and batch focused consolidation at the
next checkpoint or milestone audit. No active-state/provenance ambiguity or
hard integrity failure is introduced. MEMORY and the handoff are unchanged.

## Retained boundaries and next gate

All scientific/runtime-analysis source, equivalence projection and strict
validators are byte-preserved. Profiles, generic collector, launcher, external
configuration and unrelated tracked files are unchanged. The supplied contract
is byte-preserved; index, HEAD and local origin/test are unchanged. No farm,
ROOT/PyROOT, real analysis, staging, commit, push or ref update occurred.

Local checks establish only source behavior. Repaired farm staging and full
canonical-five runtime/PDF closure remain NOT VERIFIED. Canonical-five full
runtime, final E.8/F.6.4 and absolute-SIMC interpretation remain `BLOCKED`.
Existing closed scopes remain intact; Method A stays detached/non-production
and Method B stays diagnostic/cross-check only and numerically excluded.
CURRENT alone owns the fresh farm-readiness NEXT conditional on actual-diff
review, user commit/push and pushed-state synchronization. Only Farm readiness:
PASS may authorize the isolated retry.
