# Local (t,phi) missing-fit serialization repair

Starting test/origin/test: 240f30c987cd9cd7cf51cc65958ee05ebf757396.
User authorizes this narrowly scoped repair after the supplied October 27
21:05:03 UTC attempt failed at verify_artifacts / child_failed_fit_evidence.
Follow root AGENTS and MAINTENANCE. No scientific phase or production promotion.

## Audit and failed invariant

The supplied gate-status receipt reports analysis_completed=true, returncode=0,
artifact verification/collection incomplete. Original inputs, accepted artifacts,
checkout, external links and installed ltsep preservation passed. The triggering
cell and fit exception are not retained in this receipt. Cleanup may remove the
private diagnostic JSON on unwinding; do not infer its continued existence.
Receipt SHA256: 376ae2d6a6c93c4332c4a64c84d67320fdf2e1fc6f84284a4f0a29e8b87f5adf.
Source-verified: observer finally writes snapshot(result); the real audit
snapshot(None) becomes exactly {"status":"unavailable","reason":"not retained by producer"}.
Owner instead demands literal None for non-returned fit/fit_acceptance. Existing
fixtures write raw None or mock snapshot with deepcopy, concealing the mismatch.

## Allowed and frozen files

Only testing/run_pion_t_template_comparison_gate.py, its existing test file,
docs/memory/CURRENT.md, this contract and regenerated docs/memory/manifest.json.
Freeze every other file, USER, refs and index. Freeze production, observer,
serializer, fitter, fit acceptance, protected models, alignment, templates,
cuts, normalizations, signed weights/Sumw2, launcher, profiles, collector,
private-output isolation and cleanup. Preserve existing failed/accepted evidence.

## Before/after and tests

Recognize only explicit null or the exact real serializer missing-value record
for both non-returned fit fields. Require field presence, false acceptance,
nonempty string reason and unchanged protected/status checks. Do not interpret
generic unavailable, malformed or populated fit evidence as absent. Returned
fits must still pass all evidence/acceptance guards. Rejected or failed child
fits must never supply Local_t_phi transfer weights.

Use the actual ROOT-independent audit module, snapshot and JSON roundtrip in
the existing owner suite; demonstrate rejection before the fix, acceptance
after it for failed and wholly empty cells, and idempotent serialization.
Test missing fields, altered/extra sentinel fields, wrong types, claimed
acceptance, mismatched protected status, false emptiness, returned missing
evidence and attempted transfer from failed fits. Retain accepted/rejected and
signed application regressions. No fake ROOT result is farm evidence.

Run WSL owner suite and affected observer/renderer/event/statistics regressions
using available local dependencies. Check Python 3.9 grammar, full diff,
manifest and ordinary memory health. Verify allowed changes only, index/refs
and frozen file hashes. Export complete tracked plus new-file diff without
staging. No commits, pushes, farm execution, automatic retries or historical
repinning. Hard stop if production/serializer/cleanup changes become necessary.

## Acceptance and NEXT

Local deterministic validation plus full actual diff prepared for independent
ChatGPT review. User publication and pushed-state verification precede a single
existing private Center/Left six-artifact gate. Farm/ROOT/scientific acceptance
remains NOT VERIFIED. No farm command belongs in this repair.

## Local outcome

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING. Both failed and wholly empty
serialized fixtures reproduced child_failed_fit_evidence before repair. The
exact unavailable records now pass child and six-artifact checks with zero
accepted fits. Explicit null compatibility remains tested by existing cases.
Two new test methods cover serialization idempotence, both populations/settings,
missing/malformed/extra-field records, false acceptance, protected status,
returned missing evidence, false emptiness and forbidden failed-fit transfer.
The scalar application fixture also uses the real serializer instead of deepcopy.
An initial new fixture reused Center's PDF identity for Left; corrected only
that fixture to preserve setting-specific PDF/runtime provenance.

Final WSL owner suite: 77 passed. Existing observer, renderer, event and
statistics compatibility suites: 98 cases, 96 passed, 2 historical-artifact
skips. Combined: 175 cases, 173 passed, 2 skipped, zero failures. These are
local deterministic results, not farm/PyROOT acceptance. Only the bounded
owner absent-evidence branch and tests changed; fit decisions, application
checks, source pins, isolation and cleanup remain unchanged.

Supplied receipt proves analysis completion and the listed preservation checks,
not fit adequacy, collection, ZIP delivery or continued private-file existence.
No farm rerun, commit, push, index or reference mutation occurred. Independent
actual-diff review remains NEXT; author inspection is not independent review.
