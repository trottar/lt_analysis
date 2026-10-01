# KaonLT Post-Push CURRENT Continuity — Final Pre-Push Reconciliation

## 1. Purpose

Perform the required final pre-push memory/status reconciliation after
independent ChatGPT actual-diff review **PASSED** the local post-push CURRENT
continuity / push-stable NEXT candidate.

Reviewed cumulative bundle:

```text
kaonlt_review(20260930-194901).diff
```

Required committed base:

```text
3fd4edfcd05fd4cc25fb4e2119fb5f212ee916a1
```

Independent review established that:

- remote `test` still points to exactly
  `3fd4edfcd05fd4cc25fb4e2119fb5f212ee916a1`;
- the candidate changes only the seven allowlisted memory/workflow paths;
- no scientific/runtime/production source changes;
- CURRENT no longer describes the already-pushed hardening work as local;
- CURRENT records pushed hardening source `3fd4edfcd...` and independent
  pushed-state source/provenance PASS;
- F.4.Refresh.2 farm execution remains `BLOCKED`;
- CURRENT has exactly one push-stable NEXT:
  `NEXT — after user-controlled commit/push and pushed-state review of this
  continuity reconciliation, audit and contract the missing tracked
  F.4.Refresh.2 materialize -> verify -> package execution owner.`;
- CODEX, MAINTENANCE, and the contract template now require push-stable NEXT
  semantics and an explicit consumed-pre-push-language check during
  pushed-state review;
- no missing execution owner is implemented;
- no farm command/run occurs.

The original continuity task contract in the reviewed candidate is
byte-identical to the ChatGPT-provided contract:

```text
SHA-256:
4590c953b6822da8eb2922b14c6f03adb32f8c6be50c955ca0ee805391b20b62
```

Codex-reported deterministic checks:

```text
manifest write/check                         PASS
non-strict memory health                     PASS
strict --fail-on-warning memory health       PASS
memory bootstrap                             PASS
git diff --check                             PASS

CURRENT bytes                                7045
MEMORY bytes                                 6997
CURRENT_HANDOFF bytes                        323
health warnings                              none
live-state continuity                        PASS
```

No unit suite was required because checker implementation did not change.
These checks were **NOT RUN by ChatGPT**.

This task is final memory/status reconciliation only. It must not alter the
reviewed push-stable rule implementation or any scientific/runtime source.

---

## 2. Exact starting state

Required branch:

```text
test
```

Required committed HEAD:

```text
3fd4edfcd05fd4cc25fb4e2119fb5f212ee916a1
```

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
git log -1 --oneline
```

Expected existing local candidate paths before placing this contract:

```text
docs/memory/CODEX.md
docs/memory/CURRENT.md
docs/memory/MAINTENANCE.md
docs/memory/manifest.json
docs/memory/templates/CODEX_CONTRACT.md
docs/memory/phases/post-push-current-continuity-task-contract.md
docs/memory/phases/post-push-current-continuity.md
```

Pre-existing `kaonlt_review*.diff` files may remain untracked and outside the
candidate.

If committed HEAD differs, substantive reviewed bytes differ, or unrelated
tracked changes exist, **STOP**.

Do not reset, stash, clean, discard, commit, push, or run the farm.

---

## 3. Mandatory startup

Read in exact order:

1. root `AGENTS.md`;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read only:

- `docs/memory/phases/post-push-current-continuity.md`;
- `docs/memory/phases/post-push-current-continuity-task-contract.md`;
- `docs/memory/CODEX.md`;
- `docs/memory/MAINTENANCE.md`;
- `docs/memory/templates/CODEX_CONTRACT.md`.

Do not reopen adjacent scientific phases.

---

## 4. Files that must remain byte-identical

Do not edit:

```text
docs/memory/CODEX.md
docs/memory/MAINTENANCE.md
docs/memory/templates/CODEX_CONTRACT.md
docs/memory/phases/post-push-current-continuity-task-contract.md
```

Also freeze:

```text
tools/check_memory_health.py
testing/test_memory_health.py
testing/materialize_method_a_current_baseline_authority.py
testing/compare_method_a_current_baseline_authority.py
testing/collect_pion_hgcer_validation_bundle.py
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
testing/package_pion_hgcer_validation_bundle.tcsh
```

and all scientific/runtime/production source.

If a reviewed substantive file needs repair, **STOP**; do not hide a repair
inside reconciliation.

---

## 5. Allowed reconciliation changes

Only:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/post-push-current-continuity.md
```

Create exactly:

```text
docs/memory/phases/post-push-current-continuity-final-pre-push-reconciliation-task-contract.md
```

No other path may change.

---

## 6. Required status reconciliation

Advance:

```text
Post-push CURRENT continuity / push-stable NEXT — ACTIVE
```

to:

```text
Post-push CURRENT continuity / push-stable NEXT — SOURCE REVIEWED
```

Record that independent ChatGPT actual-diff review of:

```text
kaonlt_review(20260930-194901).diff
```

passed.

Preserve:

```text
memory-health / operational-completeness hardening — SOURCE REVIEWED
F.4.Refresh.2.Validation.1                           — SOURCE REVIEWED
F.4.Refresh.2 farm execution gate                   — BLOCKED
```

No farm/runtime validation follows.

Do not implement the missing execution owner.

---

## 7. CURRENT reconciliation

Update CURRENT minimally to make the continuity candidate itself source
reviewed while preserving push-stable semantics.

Requirements:

- CURRENT remains <= 7168 bytes;
- schema 3 and required section order remain unchanged;
- exactly one ordinary NEXT remains;
- no wording says this continuity candidate has already been committed/pushed;
- no wording makes commit/push itself the sole NEXT;
- no scientific/runtime status changes.

The exact sole NEXT must remain:

```text
NEXT — after user-controlled commit/push and pushed-state review of this continuity reconciliation, audit and contract the missing tracked F.4.Refresh.2 materialize -> verify -> package execution owner.
```

This NEXT is intentionally unchanged across the later commit/push.

---

## 8. Required health / continuity checks

Use the established `<PYTHON>`.

Run:

```bash
<PYTHON> -B tools/update_memory_manifest.py --root . --write
<PYTHON> -B tools/update_memory_manifest.py --root . --check

<PYTHON> -B tools/check_memory_health.py --root .
<PYTHON> -B tools/check_memory_health.py --root . --fail-on-warning

<PYTHON> -B tools/memory_bootstrap.py --root . --json

wc -c \
  docs/memory/CURRENT.md \
  docs/memory/MEMORY.md \
  docs/memory/handoffs/CURRENT_HANDOFF.md

git -c core.safecrlf=false diff --check
```

Acceptance requires:

```text
Memory health: PASS
CURRENT bytes: <= 7168
health warnings: none
manifest check: PASS
live-state continuity: PASS
```

No maintenance exception is allowed.

---

## 9. Required cumulative diff audit

Before stopping:

```bash
git status --short
git -c core.safecrlf=false diff --stat
git -c core.safecrlf=false diff --no-ext-diff
```

The cumulative candidate must contain exactly these eight paths:

```text
docs/memory/CODEX.md
docs/memory/CURRENT.md
docs/memory/MAINTENANCE.md
docs/memory/manifest.json
docs/memory/templates/CODEX_CONTRACT.md
docs/memory/phases/post-push-current-continuity-task-contract.md
docs/memory/phases/post-push-current-continuity.md
docs/memory/phases/post-push-current-continuity-final-pre-push-reconciliation-task-contract.md
```

Temporary review bundles remain untracked and outside the candidate.

---

## 10. Required final review bundle

Create one fresh repository-root bundle:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

It must contain byte-faithfully:

1. branch and committed HEAD;
2. `git status --short`;
3. cumulative diff stat;
4. complete tracked cumulative diff;
5. complete no-index diffs for every intended untracked candidate file,
   including this reconciliation contract;
6. exact deterministic reconciliation outputs;
7. exact health/continuity report;
8. final cumulative candidate inventory.

Do not stage the review bundle.

---

## 11. Acceptance criteria

Final reconciliation passes only if:

1. committed HEAD remains
   `3fd4edfcd05fd4cc25fb4e2119fb5f212ee916a1`;
2. reviewed push-stable rule files remain byte-identical;
3. original continuity task contract remains byte-identical with SHA-256
   `4590c953b6822da8eb2922b14c6f03adb32f8c6be50c955ca0ee805391b20b62`;
4. no scientific/runtime/production source changes;
5. continuity hardening advances to `SOURCE REVIEWED`;
6. F.4.Refresh.2 farm execution remains `BLOCKED`;
7. CURRENT contains exactly the required push-stable NEXT;
8. CURRENT contains no consumed push wording;
9. CURRENT remains <= 7168 bytes;
10. strict health has zero warnings;
11. manifest/bootstrap/diff checks pass;
12. live-state continuity reports PASS;
13. no execution-owner implementation occurs;
14. no farm command/run occurs;
15. one fresh cumulative review bundle is produced.

---

## 12. Hard stop

After final reconciliation, checks, health/continuity report, diff audit, and
fresh review bundle:

**STOP.**

Do not commit or push.

Do not implement the missing F.4.Refresh.2 execution owner.

Do not provide or run a farm command.

Return the fresh cumulative review bundle for independent ChatGPT final
pre-push review.
