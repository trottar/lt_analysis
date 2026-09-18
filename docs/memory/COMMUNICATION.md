# KaonLT farm communication

Use the detailed [farm validation-bundle procedure](decisions/farm-validation-bundle-procedure.md)
as the authoritative workflow. This page records the concise delivery rules.

- The Jefferson Lab farm shell is `tcsh`.
- Analysis artifacts use
  `/group/c-kaonlt/USERS/trottar/lt_analysis/OUTPUT/Analysis/KaonLT`.
- Transferable bundles use `/volatile/hallc/c-kaonlt/trottar/globus`.
- Ordinary unrelated modified or untracked farm-local files are not
  automatically contamination. Do not clean, reset, stash, or disturb them
  merely to package detached evidence.
- For a bundle-only request, do not rerun analysis or regenerate a frozen
  scientific artifact. Use a temporary detached worktree only when source
  provenance requires it.
- State a simple farm operation as one concise, exact `tcsh` command block.
  Where a name can collide, give the user a unique filename.
- Never claim farm, ROOT/PyROOT, full-runtime, or production validation without
  direct supplied evidence.

## Execution authority

Codex may make allowlisted local changes and run deterministic local checks. It
must not commit, push, update remote refs, or initiate Jefferson Lab farm
execution. ChatGPT independently audits the actual diff. The user alone
commits/pushes accepted changes and runs farm validation. Workflow: Codex local
changes -> ChatGPT audit -> user commit/push -> user farm run -> ChatGPT
evidence review.
