"""Read-only, fail-closed diagnosis of the unchanged normal launcher's paths.

This is a bounded BLOCKED receipt, not a complete isolation implementation.
It stops before the first external write and never executes setup or analysis.
The tests reproduce that write only in an exclusively synthetic filesystem.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess


SOURCE_HEAD = "86c7570826b824b4f2bd2eabae57dffeab4afbc0"
SOURCE_PATHS = (
    "run_Prod_Analysis.sh", "set_SymLinks.sh",
    "farm_env/print_ltsep_path_fields.py", "farm_env/ltsep_paths.py",
    "testing/run_e8_4_fix5_canonical_five_plot_gate.py",
    "src/cuts/pion_background_chain_audit.py", "src/cuts/pion_t_bin_parents.py",
)
FIRST_SHARED_WRITE = "post_clean_background_output_mkdir_outside_attempt"


def git(root: Path, *args: str) -> bytes:
    return subprocess.run(["git", *args], cwd=root, check=True,
                          capture_output=True).stdout


def checked_source(root: Path) -> tuple[dict, dict, str]:
    """Bind the current commit to the frozen ancestor and protected contents.

    Only repository-memory and the isolation/attestation helpers/tests may
    advance after the pin.
    Protected committed and index bytes must match the pin exactly; worktree
    comparison permits only Git's Windows CRLF/LF checkout conversion.
    """
    head = git(root, "rev-parse", "HEAD").decode().strip()
    ancestor = subprocess.run(["git", "merge-base", "--is-ancestor", SOURCE_HEAD, head],
                              cwd=root, capture_output=True)
    if ancestor.returncode != 0:
        raise ValueError("scientific_source_not_verified_ancestor")
    texts, hashes = {}, {}
    for name in SOURCE_PATHS:
        raw = git(root, "show", f"{SOURCE_HEAD}:{name}")
        if (git(root, "ls-tree", head, "--", name) !=
                git(root, "ls-tree", SOURCE_HEAD, "--", name) or
                git(root, "show", f":{name}") != raw):
            raise ValueError("protected_source_changed:" + name)
        # WSL and Windows Git can have different autocrlf defaults. Compare
        # bytes allowing only newline conversion, not arbitrary source filters.
        if ((root / name).is_symlink() or
                (root / name).read_bytes().replace(b"\r\n", b"\n") != raw.replace(b"\r\n", b"\n")):
            raise ValueError("protected_source_changed:" + name)
        texts[name] = raw.decode("utf-8").replace("\r\n", "\n")
        hashes[name] = hashlib.sha256(raw).hexdigest()
    allowed = {"testing/pion_chain_attempt_output_isolation.py",
               "testing/test_pion_chain_attempt_output_isolation.py",
               "testing/attest_pion_chain_ltsep_default_root.py",
               "testing/test_attest_pion_chain_ltsep_default_root.py"}
    changed = git(root, "diff", "--name-only", "--no-renames", "-z", SOURCE_HEAD, head)
    for name in changed.decode("utf-8").split("\0"):
        if name and name not in allowed and not name.startswith("docs/memory/"):
            raise ValueError("unapproved_source_change:" + name)
    launcher = texts["run_Prod_Analysis.sh"]
    setup = texts["set_SymLinks.sh"]
    ordered = ('mkdir -p "${SKIM_OUTPUT_DIR}"', 'mkdir -p "${OUTPATH}"',
               'git clean -fdx\n', './set_SymLinks.sh $ParticleType')
    positions = [launcher.index(token) for token in ordered]
    if positions != sorted(positions):
        raise ValueError("normal_launcher_order_changed")
    assignment = 'BACKGROUND_OUTPUT_DIR="${VOLATILEPATH}/OUTPUT/Analysis/${ANATYPE}LT/background_samples/"'
    mkdir = 'mkdir -p "${BACKGROUND_OUTPUT_DIR}"'
    if not setup.index(assignment) < setup.index(mkdir) < setup.index('# Check symlinks'):
        raise ValueError("background_output_order_changed")
    if 'declare -a EPS=("low" "high")' not in launcher:
        raise ValueError("normal_epsilon_scope_changed")
    return texts, hashes, head


def target_record(role: str, path: Path, attempt: Path, operation: str) -> dict:
    """Resolve existing chains strictly; resolve missing leaves via real parents."""
    path = Path(path)
    if not path.is_absolute():
        raise ValueError("relative_path:" + role)
    # strict=False alone silently accepts a dangling symlink. Inspect each hop.
    for node in (path, *path.parents):
        if node.is_symlink():
            try:
                node.resolve(strict=True)
            except (OSError, RuntimeError) as error:
                raise ValueError("broken_or_cyclic_symlink:" + role) from error
    real = path.resolve()
    return {"role": role, "path": str(path), "realpath": str(real),
            "owner": "attempt" if real.is_relative_to(attempt) else "external",
            "operation": operation, "exists": path.exists()}


def preflight(source_root: Path, attempt_root: Path, roles: dict,
              command: tuple = ("./run_Prod_Analysis.sh", "4p4", "2p74")) -> dict:
    """Inspect supplied probe roles without creating, unlinking or writing paths.

    A later owner would have to establish real copied-ltsep/import identities.
    Neither supplied roles nor a local fixture authorize an actual launcher.
    No PASS route exists while downstream write coverage remains unproved.
    """
    receipt = {"schema_version": 1, "checked_source_HEAD": None,
               "scientific_source_HEAD": SOURCE_HEAD,
               "analysis_command_identity": list(command),
               "epsilon_scope": ["low", "high"], "passed": False,
               "failed_invariant": None, "farm_execution_authorized": False,
               "role_resolved_write_targets": [], "protected_source_hashes": {},
               "pre_post_clean_link_state": "not_executed_by_read_only_preflight",
               "external_immutable_identities": "requires_independent_owner_probe",
               "verified_operations": [],
               "unknown_operations": ["real installed-ltsep identity",
                   "post-setup attestation", "downstream cleanup/merge/SIMC/model/log writes",
                   "companion JSON/PDF destinations", "full-launcher isolation"]}

    def block(invariant):
        receipt["failed_invariant"] = invariant
        return receipt

    try:
        _, hashes, checked_head = checked_source(Path(source_root))
        receipt["checked_source_HEAD"] = checked_head
        receipt["protected_source_hashes"] = hashes
        receipt["verified_operations"].append("normal launcher: mkdirs -> clean -> link setup; low+high")
        if command != ("./run_Prod_Analysis.sh", "4p4", "2p74"):
            return block("normal_launcher_command_required_no_debug_shortcut")
        attempt_root = Path(attempt_root)
        if not attempt_root.is_absolute() or attempt_root.is_symlink():
            return block("attempt_root_not_owned_real_directory")
        attempt = attempt_root.resolve(strict=True)
        if not attempt.is_dir() or (hasattr(os, "getuid") and attempt.stat().st_uid != os.getuid()):
            return block("attempt_root_not_owned_real_directory")
        if any(attempt.iterdir()):
            return block("attempt_root_collision_or_stale")
        # Roles come from path probes, never ambient overrides or guessed values.
        required = ("LTANAPATH", "VOLATILEPATH", "OUTPATH", "SKIMPATH", "ANATYPE")
        if any(not isinstance(roles.get(k), str) or not roles[k] for k in required):
            return block("required_path_role_missing")
        if roles["ANATYPE"] != "Kaon":
            return block("analysis_type_not_kaon")
        worktree = target_record("LTANAPATH", Path(roles["LTANAPATH"]), attempt, "worktree")
        if worktree["owner"] != "attempt":
            return block("ltanapath_outside_attempt")
        # Launcher mkdir -p of an already existing real directory is a no-op.
        skim = Path(roles["SKIMPATH"])
        if skim.name != "KaonLT":
            skim = skim / "Analysis/KaonLT"
        operations = (
            ("SKIM_OUTPUT_DIR", skim, "pre_clean_skim_mkdir_outside_attempt"),
            ("OUTPATH", Path(roles["OUTPATH"]), "pre_clean_outpath_mkdir_outside_attempt"),
            ("BACKGROUND_OUTPUT_DIR", Path(roles["VOLATILEPATH"]) /
             "OUTPUT/Analysis/KaonLT/background_samples", FIRST_SHARED_WRITE),
        )
        for role, path, failure in operations:
            record = target_record(role, path, attempt, "mkdir -p")
            receipt["role_resolved_write_targets"].append(record)
            # Even an existing shared background directory remains a reachable
            # output surface; mkdir's no-op cannot establish downstream safety.
            if record["owner"] != "attempt":
                if role == "BACKGROUND_OUTPUT_DIR":
                    return block(failure)
                if not path.is_dir():
                    return block(failure)
            receipt["verified_operations"].append(role + ": existing directory or attempt-confined mkdir")
        return block("normal_launcher_downstream_write_coverage_unproven")
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError, KeyError) as error:
        return block("preflight_indeterminate:" + str(error))


if __name__ == "__main__":
    # No operational launcher CLI is provided at this blocked gate.
    print(json.dumps({"schema_version": 1, "passed": False,
                      "failed_invariant": "independent_path_probe_required_use_preflight_api",
                      "farm_execution_authorized": False}, sort_keys=True))
    raise SystemExit(1)
