# KaonLT Memory-Health / Operational-Completeness Hardening — Final Pre-Push Reconciliation

## 1. Purpose

Perform the required final pre-push memory/status reconciliation after
independent ChatGPT actual-diff review **PASSED** the local
memory-health and farm-gate operational-completeness hardening candidate.

Reviewed cumulative bundle:

```text
kaonlt_review(20260930-173355).diff
```

Required committed base:

```text
86590fa655512926f2e4d0c50bf12b57d5198da5
```

Independent review established that the candidate:

- changes only the allowlisted repository-memory, memory-health tooling, and
  focused memory-health test paths;
- records the F.4.Refresh.2.Validation.1 operational-readiness failure
  accurately and without changing scientific/runtime status;
- preserves the pushed F.4.Refresh.2 materializer/profile/collector/wrapper and
  all scientific/runtime/production source unchanged;
- adds the required pre-farm complete-execution-chain audit and the rule that
  unowned multi-step orchestration blocks farm handoff;
- requires `Farm readiness: PASS` plus named tracked source owners before any
  farm command and `Farm readiness: BLOCKED` with no command otherwise;
- adds `--fail-on-warning` to `tools/check_memory_health.py` while preserving
  default diagnostic warning behavior;
- adds focused deterministic tests for healthy strict mode, soft-warning
  compatibility, strict warning failure, and hard-limit failure;
- makes unresolved health warnings workflow-stopping unless the user approves a
  recorded maintenance exception;
- requires explicit Codex and ChatGPT health updates after substantial gates;
- compacts `docs/memory/CURRENT.md` from 16,338 bytes to 7,035 bytes while
  preserving schema 3, the required section order, one ordinary NEXT, active
  E.8 state, frozen F.6.2 identity, F.4.Refresh.1 closure, and F.4.Refresh.2
  source-review boundaries;
- leaves F.4.Refresh.2.Validation.1 `SOURCE REVIEWED` while marking its farm
  execution gate `BLOCKED` pending a separately implemented/reviewed execution
  owner;
- does not implement that missing execution owner.

The supplied hardening task contract in the candidate is byte-identical to the
ChatGPT-provided contract, SHA-256:

```text
9ad4767cb66ce98796b77fc08b1c52ee6478cf1132b9f652c6cf57d560a39322
```

Codex-reported local validation in the reviewed bundle:

```text
py_compile                                            PASS
testing.test_memory_health                            36 tests OK, 0 skips
manifest write/check                                  PASS
non-strict memory health                              PASS
strict --fail-on-warning memory health                PASS
memory bootstrap                                      PASS
git diff --check                                      PASS

CURRENT bytes                                         7035
MEMORY bytes                                          6997
CURRENT_HANDOFF bytes                                 323
health warnings                                       none
```

These unit tests were **NOT RUN by ChatGPT**. ChatGPT independently inspected
the actual cumulative diff, the strict-mode implementation and focused test,
the final candidate inventory, the contract identity, and the reported health
evidence.

This task is final memory/status reconciliation only. Do not alter the reviewed
tooling, tests, investigation, operational rules, or scientific/runtime source.

---

## 2. Exact starting state

Required branch:

```text
test
```

Required committed HEAD:

```text
86590fa655512926f2e4d0c50bf12b57d5198da5
```

Before editing:

```bash
git branch --show-current
git rev-parse HEAD
git status --short
git log -1 --oneline
```

The cumulative local candidate must already contain exactly the reviewed
hardening paths from `kaonlt_review(20260930-173355).diff`, plus pre-existing
temporary `kaonlt_review*.diff` files outside the candidate.

If committed HEAD differs, if any substantive reviewed file differs from the
reviewed bundle, or if unrelated tracked changes are present, **STOP**.

Do not reset, stash, clean, discard, commit, push, or run the farm.

---

## 3. Mandatory startup

Read in this exact order:

1. root `AGENTS.md`;
2. `docs/memory/CURRENT.md`;
3. `docs/memory/MEMORY.md`;
4. `docs/memory/handoffs/CURRENT_HANDOFF.md`;
5. `docs/memory/USER.md`.

Then read only:

- `docs/memory/phases/memory-health-operational-completeness-hardening.md`
- `docs/memory/investigations/f4-refresh2-validation1-operational-readiness-failure.md`
- `docs/memory/phases/f4-refresh2-validation1-farm-materialization-bundle-profile.md`
- `docs/memory/roadmap/STATUS.md`
- this reconciliation contract.

Do not reopen adjacent phases.

---

## 4. Files that must remain byte-identical

Do not edit:

```text
docs/memory/AGENTS.md
docs/memory/CODEX.md
docs/memory/COMMUNICATION.md
docs/memory/LEARNINGS.md
docs/memory/MAINTENANCE.md
docs/memory/MEMORY.md
docs/memory/TOOLS.md
docs/memory/USER.md
docs/memory/templates/CODEX_CONTRACT.md

docs/memory/investigations/f4-refresh2-validation1-operational-readiness-failure.md
docs/memory/phases/memory-health-operational-completeness-hardening-task-contract.md

tools/check_memory_health.py
testing/test_memory_health.py
```

Also freeze:

```text
testing/materialize_method_a_current_baseline_authority.py
testing/test_materialize_method_a_current_baseline_authority.py
testing/compare_method_a_current_baseline_authority.py
testing/collect_pion_hgcer_validation_bundle.py
testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json
testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py
testing/package_pion_hgcer_validation_bundle.tcsh
```

and every scientific/runtime/production source path, including all `src/`,
`src/main.py`, and `run_Prod_Analysis.sh`.

If any reviewed substantive file needs a repair, **STOP**. Do not hide a repair
inside reconciliation.

---

## 5. Allowed reconciliation changes

Only:

```text
docs/memory/CURRENT.md
docs/memory/manifest.json
docs/memory/phases/memory-health-operational-completeness-hardening.md
docs/memory/roadmap/STATUS.md
```

Create exactly:

```text
docs/memory/phases/memory-health-operational-completeness-hardening-final-pre-push-reconciliation-task-contract.md
```

If roadmap status already expresses the required reviewed hardening state
without modification, it may remain unchanged. No other path may change.

---

## 6. Required status reconciliation

Advance:

```text
Memory health / operational-completeness hardening — ACTIVE
```

to:

```text
Memory health / operational-completeness hardening — SOURCE REVIEWED
```

Record that independent ChatGPT actual-diff review of:

```text
kaonlt_review(20260930-173355).diff
```

passed.

Preserve:

```text
F.4.Refresh.1              — CLOSED / RUNTIME VALIDATED
F.4.Refresh.1.Fix.1        — CLOSED / RUNTIME VALIDATED
                              detached comparator gate only

F.4.Refresh.2              — SOURCE REVIEWED
F.4.Refresh.2.Fix.1        — SOURCE REVIEWED
F.4.Refresh.2.Validation.1 — SOURCE REVIEWED

F.4.Refresh.2 farm execution gate — BLOCKED

F.6.3                      — SOURCE REVIEWED
E.8.4                      — SOURCE REVIEWED
E.8                        — ACTIVE
final E.8                  — BLOCKED
F.6.4                      — BLOCKED
lifecycle hook             — BLOCKED / DEFERRED
```

Do not claim farm/runtime validation for the hardening task.

Do not claim the missing F.4.Refresh.2 execution owner exists.

---

## 7. Required health boundary

The final reconciliation must itself satisfy the new strict health policy.

Required final health report:

```text
Memory health: PASS
CURRENT bytes: <integer <= 7168>
MEMORY bytes: <integer>
CURRENT_HANDOFF bytes: <integer>
health warnings: none
manifest check: PASS
```

If any warning exists, **STOP**.

Do not use a maintenance exception in this reconciliation.

---

## 8. CURRENT.md exact NEXT

Replace the review NEXT with exactly:

```text
NEXT — user-controlled commit/push of the independently reviewed memory-health and farm-gate operational-completeness hardening candidate.
```

There must be exactly one ordinary NEXT.

Keep `CURRENT.md` at or below:

```text
7168 bytes
```

Prefer retaining the current compact form; do not re-expand historical prose.

---

## 9. Required deterministic reconciliation checks

Use the already established local `<PYTHON>`.

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

No tooling/unit-suite rerun is required for this memory-only reconciliation.
If Codex reruns tests, report them accurately.

Strict health must exit 0 with no warnings.

---

## 10. Required cumulative diff audit

Before stopping:

```bash
git status --short
git -c core.safecrlf=false diff --stat
git -c core.safecrlf=false diff --no-ext-diff
```

The cumulative candidate may contain only the previously reviewed 18 paths plus
the new final-reconciliation contract:

```text
docs/memory/AGENTS.md
docs/memory/CODEX.md
docs/memory/COMMUNICATION.md
docs/memory/CURRENT.md
docs/memory/LEARNINGS.md
docs/memory/MAINTENANCE.md
docs/memory/MEMORY.md
docs/memory/TOOLS.md
docs/memory/USER.md
docs/memory/investigations/f4-refresh2-validation1-operational-readiness-failure.md
docs/memory/manifest.json
docs/memory/phases/f4-refresh2-validation1-farm-materialization-bundle-profile.md
docs/memory/phases/memory-health-operational-completeness-hardening-task-contract.md
docs/memory/phases/memory-health-operational-completeness-hardening.md
docs/memory/phases/memory-health-operational-completeness-hardening-final-pre-push-reconciliation-task-contract.md
docs/memory/roadmap/STATUS.md
docs/memory/templates/CODEX_CONTRACT.md
testing/test_memory_health.py
tools/check_memory_health.py
```

Temporary review bundles remain untracked and outside the candidate.

---

## 11. Required final review bundle

Create one fresh repository-root bundle:

```text
kaonlt_review(YYYYMMDD-HHMMSS).diff
```

It must contain:

1. branch and committed HEAD;
2. `git status --short`;
3. cumulative diff stat;
4. complete tracked cumulative diff;
5. complete no-index diffs for every intended untracked candidate file,
   including this reconciliation contract;
6. exact deterministic reconciliation commands/results;
7. the exact final health update;
8. final cumulative changed-path inventory.

The review bundle remains untracked.

---

## 12. Acceptance criteria

Final reconciliation is acceptable only if:

1. committed HEAD remains
   `86590fa655512926f2e4d0c50bf12b57d5198da5`;
2. all reviewed hardening implementation/tool/test/rule files remain unchanged;
3. the original hardening task contract remains byte-identical;
4. no scientific/runtime/production source changes;
5. hardening advances to `SOURCE REVIEWED`;
6. F.4.Refresh.2.Validation.1 remains `SOURCE REVIEWED`;
7. F.4.Refresh.2 farm execution remains `BLOCKED`;
8. no missing execution driver is implemented;
9. CURRENT has exactly the required one NEXT and remains <= 7168 bytes;
10. manifest check passes;
11. strict health exits 0 with no warnings;
12. bootstrap and diff check pass;
13. exact memory byte counts are reported;
14. no farm command or farm run occurs;
15. one fresh cumulative review bundle is produced.

---

## 13. Hard stop

After reconciliation, strict health, manifest/bootstrap checks, diff audit, and
creation of the fresh review bundle:

**STOP.**

Do not commit or push.

Do not implement the missing F.4.Refresh.2 execution driver.

Do not provide or run a farm command.

Do not run the materializer or collector.

Do not create a farm ZIP.

Do not change accepted authority.

Return the fresh cumulative review bundle for independent ChatGPT final
pre-push review.
