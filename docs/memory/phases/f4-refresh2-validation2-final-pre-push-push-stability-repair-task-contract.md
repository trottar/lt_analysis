# F.4.Refresh.2.Validation.2/Fix.1 — Final Pre-Push Push-Stability Repair

## 1. Objective

Repair one narrow repository-memory continuity defect found during independent
ChatGPT review of the final pre-push reconciliation bundle:

```text
kaonlt_review(20261001-101444).diff
```

The substantive F.4.Refresh.2.Validation.2/Fix.1 implementation review remains
**PASS**. Do not reopen or alter the reviewed execution owner, profile, focused
tests, materializer, collector, wrapper, scientific/runtime source, or accepted
authority.

The final reconciliation correctly advanced Validation.2 and Fix.1 to
`SOURCE REVIEWED`, preserved the substantive candidate byte-for-byte, retained
strict warning-free memory health, and used the correct push-stable ordinary
NEXT.

The remaining defect is only consumed pre-push wording in active/durable memory.
Several records still say the reviewed owner is currently "uncommitted" and/or
that the "final reconciliation review is pending." Those statements become
false immediately after the already-in-progress ChatGPT reconciliation review
passes and/or after the user's subsequent push. That would recreate the exact
post-push CURRENT-continuity defect the repository hardening was intended to
prevent.

Repair only that wording. Make the live memory push-stable across the upcoming
user-controlled commit/push while preserving the evidence chronology.

Required committed base:

```text
c406c138285727503b12115177dfb8bc7efcb7fe
Harden post-push CURRENT continuity
```

Remote `test` was independently observed by ChatGPT at that exact HEAD during
review of `kaonlt_review(20261001-101444).diff`.

---

## 2. Review finding that owns this repair

Independent ChatGPT review found these substantive implementation facts remain
accepted:

- F.4.Refresh.2.Validation.2 and Validation.2.Fix.1 are `SOURCE REVIEWED`;
- the reviewed implementation bundle is
  `kaonlt_review(20261001-092421).diff`;
- the final-reconciliation bundle is
  `kaonlt_review(20261001-101444).diff`;
- the four substantive candidate files in the latter are byte-identical to the
  source-reviewed implementation bundle;
- the tracked owner still follows:
  preflight -> reviewed materializer -> explicit completion/provenance
  verification -> reviewed package wrapper/collector -> returned-ZIP
  verification;
- the six-path post-materializer allowlist remains exact;
- `required_analysis_commit` remains
  `141a3d04f9e5d07be21dba14e0e63212c3990bf1`;
- frozen hardening and delegated component blobs remain unchanged;
- Codex-reported strict memory health and manifest checks passed with zero
  warnings;
- no farm, ROOT/PyROOT, full-analysis, accepted-authority, production, F.6.3,
  E.8.4, or Method-A-promotion validation occurred.

Do not change any of those conclusions.

The narrow blocker is push-stability wording only.

---

## 3. Exact starting state

Required branch:

```text
test
```

Required committed HEAD:

```text
c406c138285727503b12115177dfb8bc7efcb7fe
```

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short --untracked-files=all
git log -1 --oneline
```

The existing cumulative candidate must contain exactly the already reviewed
versionable paths:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/f4-refresh2-validation2-final-pre-push-source-review-reconciliation-task-contract.md
docs/memory/phases/f4-refresh2-validation2-fix1-post-hardening-source-allowlist-task-contract.md
docs/memory/phases/f4-refresh2-validation2-fix1-post-hardening-source-allowlist.md
docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner-task-contract.md
docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner.md
docs/memory/roadmap/STATUS.md
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
testing/run_f4_refresh2_materialize_verify_package.py
testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py
testing/test_run_f4_refresh2_materialize_verify_package.py
```

This newly placed repair contract is additionally permitted:

```text
docs/memory/phases/f4-refresh2-validation2-final-pre-push-push-stability-repair-task-contract.md
```

Temporary repository-root `kaonlt_review*.diff` files may remain untracked and
must not be staged, removed, or rewritten except for creation of the one new
review bundle required by this repair.

Root `AGENTS.md` and `.codex/` remain local-only/untracked.

If committed HEAD differs or unrelated tracked changes exist, **STOP**.

---

## 4. Mandatory startup

Read in this exact order:

1. root `AGENTS.md` if present;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read only task-relevant records:

```text
docs/memory/CODEX.md
docs/memory/MAINTENANCE.md
docs/memory/roadmap/STATUS.md
docs/memory/phases/post-push-current-continuity.md
docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner.md
docs/memory/phases/f4-refresh2-validation2-fix1-post-hardening-source-allowlist.md
docs/memory/phases/f4-refresh2-validation2-final-pre-push-source-review-reconciliation-task-contract.md
```

Do not reopen adjacent scientific architecture.

---

## 5. Frozen substantive candidate

These four source-reviewed substantive files must remain byte-identical to
`kaonlt_review(20261001-092421).diff` and
`kaonlt_review(20261001-101444).diff`:

```text
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py
testing/run_f4_refresh2_materialize_verify_package.py
testing/test_run_f4_refresh2_materialize_verify_package.py
```

Their reviewed Git object identities are:

```text
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
ae86d67dd892cfe7f71bd8f1f473c530b3cb3d09

testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py
069789422798f80fa51b1d66048980231855d874

testing/run_f4_refresh2_materialize_verify_package.py
4e523de26673653d62847fb9fa76ed00df3b42dd

testing/test_run_f4_refresh2_materialize_verify_package.py
32d06d7e299b12897d64efa001bf075b09851676
```

Verify them before and after the repair with `git hash-object`.

Also freeze both implementation contracts:

```text
docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner-task-contract.md
docs/memory/phases/f4-refresh2-validation2-fix1-post-hardening-source-allowlist-task-contract.md
```

and the final-pre-push reconciliation contract:

```text
docs/memory/phases/f4-refresh2-validation2-final-pre-push-source-review-reconciliation-task-contract.md
```

Do not edit any materializer, comparator, collector, package wrapper,
memory-health implementation, scientific/runtime/production source, accepted
authority, Method B code, or F.6.3/E.8.4 runtime logic.

---

## 6. Allowed changes

Modify only:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner.md
docs/memory/phases/f4-refresh2-validation2-fix1-post-hardening-source-allowlist.md
docs/memory/roadmap/STATUS.md
```

Create exactly:

```text
docs/memory/phases/f4-refresh2-validation2-final-pre-push-push-stability-repair-task-contract.md
docs/memory/phases/f4-refresh2-validation2-final-pre-push-push-stability-repair.md
```

No other versionable file may change.

---

## 7. Exact wording defect to remove

Remove present-tense/forward-looking wording whose truth is consumed by this
review or the immediately following user push.

At minimum, eliminate current-state formulations equivalent to:

```text
SOURCE REVIEWED as an uncommitted candidate
the source-reviewed owner remains uncommitted
The uncommitted owner ...
the repaired uncommitted candidate ...
the uncommitted candidate remains ...
Final memory/status reconciliation review is pending
Independent ChatGPT final reconciliation review remains pending
still requires independent ChatGPT final review
```

from `CURRENT.md`, `roadmap/STATUS.md`, and the live phase records listed above.

Do **not** erase history. Replace those statements with push-stable historical
phrasing, for example:

```text
Independent ChatGPT source review passed for the candidate based on committed
test HEAD c406c138... . That source review did not itself establish commit,
push, farm execution, or runtime validation.
```

and:

```text
At source-review closure, farm execution remained BLOCKED pending
user-controlled commit/push and independent pushed-state review.
```

Historical phrasing that explicitly describes what was true **at the time of a
past review** is acceptable.

Avoid any statement that claims the new candidate is already pushed. The user
has not pushed it yet.

Avoid any statement that will become false merely because the user performs the
next authorized commit/push.

---

## 8. CURRENT.md requirements

`CURRENT.md` must:

1. keep F.4.Refresh.2.Validation.2/Fix.1 `SOURCE REVIEWED`;
2. keep farm readiness `BLOCKED` pending user-controlled commit/push and
   independent pushed-state review;
3. contain no claim that the owner is currently "uncommitted";
4. contain no claim that this final reconciliation review is still pending;
5. not claim the candidate is already pushed;
6. preserve exactly one ordinary NEXT;
7. keep that NEXT push-stable.

Use this exact ordinary NEXT:

```text
NEXT — after user-controlled commit/push and pushed-state review of the SOURCE REVIEWED F.4.Refresh.2.Validation.2/Fix.1 execution owner, prepare the single narrow Q4p4W2p74 F.4.Refresh.2 materialize -> verify -> package farm gate; do not begin F.6.3/E.8.4 until the returned F.4.Refresh.2 evidence is reviewed.
```

That NEXT remains correct before the push, immediately after the push, and
until pushed-state review resolves the next gate.

Do not modify `CURRENT_HANDOFF.md`.

---

## 9. Phase/roadmap requirements

### Validation.2 phase

Keep:

```text
F.4.Refresh.2.Validation.2 — SOURCE REVIEWED
```

Preserve the reviewed chain and evidence boundary.

Phrase the commit/push history push-stably. Preferred form:

```text
Independent ChatGPT source review passed for the candidate based on committed
test HEAD c406c138... . That review did not itself establish commit, push, or
farm execution.
```

Do not say the candidate is currently uncommitted.

Do not say final pre-push review remains pending.

### Validation.2.Fix.1 phase

Keep:

```text
F.4.Refresh.2.Validation.2.Fix.1 — SOURCE REVIEWED
```

Preserve the exact provenance audit and hardening blob identities.

Do not say the candidate currently remains uncommitted.

Do not say final reconciliation review remains pending.

### Roadmap

Keep all existing phase statuses.

State that Validation.2/Fix.1 is `SOURCE REVIEWED` and that the farm gate
remains `BLOCKED` pending user-controlled commit/push and pushed-state review.

Do not describe final reconciliation review as pending.

Do not describe the owner as currently uncommitted.

### New repair phase record

Create:

```text
docs/memory/phases/f4-refresh2-validation2-final-pre-push-push-stability-repair.md
```

Record:

- independent ChatGPT review of
  `kaonlt_review(20261001-101444).diff` found one memory-continuity blocker;
- substantive Validation.2/Fix.1 implementation remained source-review PASS;
- four substantive candidate files were independently compared against
  `kaonlt_review(20261001-092421).diff` and found byte-identical;
- the blocker was consumed "uncommitted"/"final review pending" wording;
- this repair only makes CURRENT/roadmap/phase prose push-stable;
- no implementation, test, scientific/runtime source, authority, or farm action
  changed;
- farm readiness remains `BLOCKED` until user commit/push and pushed-state
  review;
- the next artifact is one repaired cumulative review bundle.

Status of this memory-only repair before independent ChatGPT review:

```text
DEVELOPMENT COMPLETE, FARM VALIDATION PENDING
```

Do not mark the repair itself `SOURCE REVIEWED` until ChatGPT reviews its new
cumulative diff.

---

## 10. Push-stability scan

After editing, perform an explicit scan of the intended live memory records:

```text
docs/memory/CURRENT.md
docs/memory/roadmap/STATUS.md
docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner.md
docs/memory/phases/f4-refresh2-validation2-fix1-post-hardening-source-allowlist.md
```

Search for at least these case-insensitive terms:

```text
uncommitted
review is pending
review remains pending
still requires independent ChatGPT final review
local candidate
local implementation
```

Any occurrence must be examined.

It may remain only if clearly historical and explicitly time-bounded such that
it will still be true after the user's push. Prefer eliminating ambiguous
present-tense occurrences entirely.

The final report must show the scan output and explain any retained match.

---

## 11. Deterministic checks

No substantive implementation suite rerun is required because all reviewed
implementation/test files are frozen byte-identically.

Run:

```bash
git hash-object \
  testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json \
  testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py \
  testing/run_f4_refresh2_materialize_verify_package.py \
  testing/test_run_f4_refresh2_materialize_verify_package.py

python -B tools/update_memory_manifest.py --root . --write
python -B tools/update_memory_manifest.py --root . --check
python -B tools/memory_bootstrap.py --root .
python -m unittest testing.test_memory_health
python -B tools/check_memory_health.py --root . --fail-on-warning
git -c core.safecrlf=false diff --check
```

The four substantive object identities must remain exactly:

```text
ae86d67dd892cfe7f71bd8f1f473c530b3cb3d09
069789422798f80fa51b1d66048980231855d874
4e523de26673653d62847fb9fa76ed00df3b42dd
32d06d7e299b12897d64efa001bf075b09851676
```

Strict memory health must have zero warnings.

Required health report:

```text
Memory health: PASS | BLOCKED
CURRENT bytes: <integer>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: <none or exact list>
manifest check: PASS | FAIL
```

---

## 12. Diff audit

Before stopping:

```bash
git status --short --untracked-files=all
git -c core.safecrlf=false diff --stat
git -c core.safecrlf=false diff --no-ext-diff
```

The cumulative intended versionable candidate after this repair must contain
exactly:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/f4-refresh2-validation2-final-pre-push-push-stability-repair-task-contract.md
docs/memory/phases/f4-refresh2-validation2-final-pre-push-push-stability-repair.md
docs/memory/phases/f4-refresh2-validation2-final-pre-push-source-review-reconciliation-task-contract.md
docs/memory/phases/f4-refresh2-validation2-fix1-post-hardening-source-allowlist-task-contract.md
docs/memory/phases/f4-refresh2-validation2-fix1-post-hardening-source-allowlist.md
docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner-task-contract.md
docs/memory/phases/f4-refresh2-validation2-tracked-execution-owner.md
docs/memory/roadmap/STATUS.md
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
testing/run_f4_refresh2_materialize_verify_package.py
testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py
testing/test_run_f4_refresh2_materialize_verify_package.py
```

Temporary `kaonlt_review*.diff` files remain outside the versionable candidate.

No `git add -A`.

---

## 13. Fresh repair review bundle

Create one fresh repository-root:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

It must contain byte-faithfully:

1. branch;
2. committed HEAD;
3. full `git status --short --untracked-files=all`;
4. cumulative tracked `git diff --stat`;
5. complete cumulative tracked `git diff --no-ext-diff`;
6. complete `git diff --no-index /dev/null ...` for every intended new/untracked
   versionable file;
7. the four substantive candidate `git hash-object` identities;
8. the explicit push-stability scan and any retained-match explanation;
9. exact manifest/memory/check command outputs;
10. memory-health report;
11. final cumulative changed-path inventory;
12. explicit statement that substantive implementation remains byte-identical
    and no farm/commit/push occurred.

Do not stage the review bundle.

---

## 14. Acceptance criteria

This repair is complete locally only if:

1. committed HEAD remains
   `c406c138285727503b12115177dfb8bc7efcb7fe`;
2. all four substantive source-reviewed candidate files retain their exact
   reviewed Git object identities;
3. both implementation contracts and the final-reconciliation contract remain
   unchanged;
4. Validation.2 and Fix.1 remain `SOURCE REVIEWED`;
5. CURRENT/roadmap/phase live prose contains no consumed present-tense claim that
   the owner is currently uncommitted;
6. CURRENT/roadmap/phase live prose contains no consumed claim that the final
   reconciliation review is still pending;
7. no claim says the new candidate is already pushed;
8. CURRENT keeps exactly the required single push-stable NEXT;
9. farm readiness remains `BLOCKED` pending user commit/push and pushed-state
   review;
10. all scientific/runtime/authority statuses and boundaries remain unchanged;
11. manifest check passes;
12. strict memory health passes with zero warnings;
13. `git diff --check` passes;
14. one fresh complete cumulative review bundle is created;
15. no commit, push, farm, ROOT/PyROOT, full analysis, authority update,
    F.6.3/E.8.4, production correction, or Method-A promotion occurs.

---

## 15. Hard stop

After the push-stability wording repair, manifest regeneration/check, strict
memory-health validation, identity audit, diff audit, and creation of the fresh
cumulative review bundle:

**STOP.**

Do not edit substantive implementation source.

Do not perform another general final-pre-push redesign.

Do not commit.

Do not push.

Do not run the farm.

Do not provide a farm command.

Return the fresh cumulative review bundle for independent ChatGPT review. If
that repair review passes, the next action is the user-controlled scoped
commit/push followed by independent pushed-state review.
