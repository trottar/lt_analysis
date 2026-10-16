#!/usr/bin/env python3
"""Read-only resolver attestation; never a full-launcher isolation approval.

Local fixtures validate this tool only. An installed-package receipt requires
running the reviewed CLI in that package's environment after publication.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import uuid

try:
    from . import run_e8_4_fix5_canonical_five_plot_gate as primitives
except ImportError:
    import run_e8_4_fix5_canonical_five_plot_gate as primitives

BASELINE = "92fa0dd578c8911e503f4b0ae9bf2434f01864b1"
TOOL = "testing/attest_pion_chain_ltsep_default_root.py"
TEST = "testing/test_attest_pion_chain_ltsep_default_root.py"
CALLERS = (
    # role, actual resolver argument, mode, inherited cwd, frozen source owner.
    # The normal owner launches from the worktree root. Shell SCRIPT_DIR command
    # substitutions do not change that cwd, including before set_SymLinks.sh.
    ("effective_charge", "src/setup/findEffectiveCharge.py", (), "src/setup", "src/setup/findEffectiveCharge.py"),
    ("efficiency_reader", "src/setup/getEfficiencyValue.py", (), "src/setup", "src/setup/getEfficiencyValue.py"),
    ("data_table_reader", "src/setup/getDataTable.py", (), "src/setup", "src/setup/getDataTable.py"),
    ("launcher", ".", ("Plot_LTSep",), ".", "run_Prod_Analysis.sh"),
    ("analysis", "src/main.py", ("Plot_LTSep",), "src", "src/main.py"),
    ("model_setup", "src/setup/set_sig_fortran.py", ("Plot_LTSep",), "src/setup", "src/setup/set_sig_fortran.py"),
    ("symlink_setup", ".", ("Plot_LTSep",), ".", "set_SymLinks.sh"),
)
FROZEN = tuple(source for _, _, _, _, source in CALLERS) + (
    "run_Prod_Analysis.sh", "farm_env/ltsep_paths.py",
    "farm_env/print_ltsep_path_fields.py", "src/utility/utility.py",
    "src/utility/background_sample_manifest.py",
    "background_samples/background_samples.conf",
    "testing/run_e8_4_fix5_canonical_five_plot_gate.py",
    "testing/run_e8_4_fix5_left_lowe_plot_gate.py",
    "testing/collect_pion_hgcer_validation_bundle.py",
)
FIELDS = primitives.PATH_FIELDS + ("CACHEPATH",)
LIMITS = ("resolver only; complete launcher isolation NOT VERIFIED",
          "ROOT integration and pion scientific correctness NOT VERIFIED",
          "post_clean_background_output_mkdir_outside_attempt remains BLOCKED")


def require(condition, invariant):
    if not condition:
        raise ValueError(invariant)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def git(repo, *args):
    env = dict(os.environ, GIT_OPTIONAL_LOCKS="0")
    return subprocess.run(["git", *args], cwd=repo, env=env, check=True,
                          capture_output=True).stdout


def real(path):
    """Resolve missing leaves, but never accept dangling/cyclic symlinks."""
    path = Path(path)
    require(path.is_absolute() and ".." not in path.parts, "path_not_absolute_or_traversal")
    try:
        for node in (path, *path.parents):
            if node.is_symlink():
                node.resolve(strict=True)
        return path.resolve()
    except (OSError, RuntimeError) as error:
        raise ValueError("path_broken_or_cyclic") from error


def tree_bytes(root):
    """Package inventory only: no scientific input file reads or link traversal."""
    result = {}
    for path in sorted(Path(root).rglob("*")):
        if "__pycache__" in path.parts or path.suffix == ".pyc":
            continue
        require(not path.is_symlink(), "package_symlink_not_supported")
        if path.is_file():
            result[path.relative_to(root).as_posix()] = digest(path.read_bytes())
    return result


def repository_state(repo):
    names = git(repo, "ls-files", "-z", "--cached", "--others", "--exclude-standard")
    files = {}
    for name in names.decode().split("\0"):
        if not name:
            continue
        path = repo / name
        files[name] = ("link:" + os.readlink(path) if path.is_symlink() else
                       digest(path.read_bytes()) if path.is_file() else "missing")
    index = Path(git(repo, "rev-parse", "--git-path", "index").decode().strip())
    if not index.is_absolute():
        index = repo / index
    return {"files": files, "index_sha256": digest(index.read_bytes()),
            "refs": git(repo, "show-ref").decode(),
            "head": git(repo, "rev-parse", "HEAD").decode().strip(),
            "status": git(repo, "status", "--porcelain=v1", "--untracked-files=all").decode()}


def source_identity(repo, commit):
    require(re.fullmatch(r"[0-9a-f]{40}", commit) is not None, "source_commit_requires_full_sha")
    require(git(repo, "branch", "--show-current").decode().strip() == "test", "wrong_branch")
    for ref in ("HEAD", "refs/remotes/origin/test"):
        require(git(repo, "rev-parse", ref).decode().strip() == commit, "source_identity_mismatch:" + ref)
    git(repo, "merge-base", "--is-ancestor", BASELINE, commit)
    changed = git(repo, "diff", "--name-only", "--no-renames", "-z", BASELINE, commit)
    require(all(n in {TOOL, TEST} or n.startswith("docs/memory/")
                for n in changed.decode().split("\0") if n), "unapproved_source_successor")
    identities = {}
    for name in (*FROZEN, TOOL, TEST):
        raw = git(repo, "show", commit + ":" + name)
        require(git(repo, "show", ":" + name) == raw, "source_index_changed:" + name)
        path = repo / name
        require(not path.is_symlink() and path.read_bytes().replace(b"\r\n", b"\n") ==
                raw.replace(b"\r\n", b"\n"), "source_worktree_changed:" + name)
        if name in FROZEN:
            require(git(repo, "ls-tree", BASELINE, "--", name) ==
                    git(repo, "ls-tree", commit, "--", name), "protected_source_changed:" + name)
        identities[name] = digest(raw)
    require(digest(Path(__file__).read_bytes().replace(b"\r\n", b"\n")) ==
            digest((repo / TOOL).read_bytes().replace(b"\r\n", b"\n")), "running_tool_identity_mismatch")
    require(not git(repo, "diff", "HEAD", "--", "src", "farm_env", "testing",
                    "run_Prod_Analysis.sh", "set_SymLinks.sh", "background_samples"), "dirty_protected_source")
    require(not git(repo, "ls-files", "--others", "--exclude-standard", "--",
                    "src", "farm_env", "testing", "background_samples"), "untracked_protected_source")
    return identities


# Run only the resolver wrapper, never the caller modules. Audit hooks prevent
# ordinary Python filesystem/process mutation by imports/Root/SetPath as well.
PROBE = r'''
import sys, os, json, inspect, importlib.util, sysconfig
from pathlib import Path
spec=importlib.util.find_spec("ltsep")
if spec is None or not spec.origin: raise ImportError("installed_ltsep_unavailable")
allowed_reads=[Path(spec.origin).resolve().parent,Path(sys.argv[1],"farm_env").resolve(),
               Path(sysconfig.get_paths()["stdlib"]).resolve()]
def readonly(event, args):
    if event == "import" and args[0].split(".")[0] in {"ROOT","numpy","root_numpy"}:
        raise RuntimeError("scientific_module_import_forbidden")
    if event == "open":
        mode, flags = args[1], args[2]
        if (isinstance(mode,str) and any(c in mode for c in "wax+")) or (isinstance(flags,int) and flags & (os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC|os.O_APPEND)):
            raise RuntimeError("resolver_write_attempt")
        if isinstance(args[0],(str,bytes,os.PathLike)):
            path=Path(os.fsdecode(args[0])).resolve()
            # -B prevents writes but otherwise still accepts preexisting pyc.
            # Force source loading so the hashed Python source owns behavior.
            if path.suffix == ".pyc": raise FileNotFoundError("source_only_resolver_probe")
            if not any(path.is_relative_to(root) for root in allowed_reads):
                raise RuntimeError("resolver_data_read_forbidden")
    if event in {"os.mkdir","os.remove","os.rmdir","os.rename","os.link","os.symlink","os.chmod","os.chown","os.truncate","os.utime","subprocess.Popen","os.system","os.exec","os.posix_spawn"}:
        raise RuntimeError("resolver_mutation_attempt:"+event)
sys.addaudithook(readonly)
sys.path.insert(0,sys.argv[1])
import ltsep, ltsep.pathing
from farm_env.ltsep_paths import create_ltsep_root
caller=sys.argv[2]; modes=json.loads(sys.argv[3]); fields=json.loads(sys.argv[4])
lt=create_ltsep_root(caller,*modes)
paths={k:getattr(lt,k,None) for k in fields}
matches=[]
p=ltsep.pathing.SetPath(caller)
for f in sorted((Path(ltsep.__file__).parent/"PATH_TO_DIR").glob("*.path")):
    values={}
    for line in f.read_text(encoding="utf-8").splitlines():
        if not line.strip(): continue
        key,sep,value=line.partition("=")
        if not sep or not key.strip() or key.strip() in values: raise ValueError("ltsep_config_invalid")
        values[key.strip()]=value.strip()
    if "LTANAPATH" not in values: continue
    observed={k:p.getPath(k) for k in values}
    expanded={k:v.replace("${USER}",str(paths["USER"])) for k,v in values.items()}
    if expanded==observed: matches.append(str(f.resolve()))
print(json.dumps({"caller":caller,"caller_realpath":os.path.realpath(caller),"mode":modes,"cwd":os.getcwd(),
    "paths":paths,"config_matches":matches,"package_file":ltsep.__file__,
    "pathing_file":ltsep.pathing.__file__,"root_file":inspect.getsourcefile(type(lt)),
    "resolver_file":sys.modules["farm_env.ltsep_paths"].__file__,
    "python":sys.executable,"python_version":sys.version},allow_nan=False))
'''


def probe(tree, caller, modes, env, cwd="."):
    path = tree / caller
    require(path.exists() and real(path) == path, "caller_path_invalid:" + caller)
    working_directory = tree / cwd
    require(working_directory.is_dir() and real(working_directory) == working_directory,
            "caller_cwd_invalid:" + caller)
    result = subprocess.run([sys.executable, "-B", "-c", PROBE, str(tree),
                             str(path), json.dumps(modes), json.dumps(FIELDS)],
                            cwd=working_directory, env=env, capture_output=True, text=True, timeout=60)
    require(result.returncode == 0, "resolver_subprocess_failed:" + caller + ":" + result.stderr[-2000:])
    require(not result.stderr.strip(), "resolver_stderr_contamination:" + caller)
    try:
        record = json.loads(result.stdout)
    except (ValueError, TypeError) as error:
        raise ValueError("resolver_json_invalid:" + caller) from error
    require(record["caller"] == str(path) and record["caller_realpath"] == str(path) and
            record["mode"] == list(modes), "resolver_caller_mismatch")
    require(record["cwd"] == str(working_directory), "resolver_cwd_mismatch")
    require(real(Path(record["resolver_file"])) == tree / "farm_env/ltsep_paths.py",
            "resolver_wrapper_import_mismatch")
    require(real(Path(record["python"])) == real(Path(sys.executable)), "interpreter_mismatch")
    require(len(record["config_matches"]) == 1, "ltsep_matching_config_not_unique")
    fields = record["paths"]
    for name in primitives.PATH_FIELDS:
        require(isinstance(fields.get(name), str) and fields[name].strip(),
                "default_root_unresolved:" + caller + ":" + name)
    record["unavailable_fields"] = [k for k, v in fields.items() if v is None]
    return record


def module_identity(record):
    package_file = Path(record["package_file"])
    require(package_file.name == "__init__.py" and real(package_file) == package_file,
            "ltsep_package_layout_invalid")
    package = package_file.parent
    for key in ("pathing_file", "root_file", "config_matches"):
        value = record[key][0] if key == "config_matches" else record[key]
        require(value and real(Path(value)).is_relative_to(package), "ltsep_module_identity_invalid:" + key)
    return package


def overlay(repo, worktree, baseline, package, before):
    """Use the reviewed matching/value-span algorithm with strict read-only probes.

    prepare_runtime_overlay's subprocess helper accepts stderr and lacks an
    import-time write guard; reuse its config parser/environment/cleanup bounds
    here, adding per-caller identity and whole-package failure preservation.
    """
    primitives.cleanup_path(repo, worktree.parent, worktree, "analysis")
    parent = worktree.parent / "runtime-python"
    require(not parent.exists(), "ltsep_overlay_already_exists")
    parent.mkdir()
    copied = parent / "ltsep"
    shutil.copytree(package, copied, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    require(tree_bytes(copied) == before, "ltsep_copy_identity_mismatch")
    selected = Path(baseline["config_matches"][0]).relative_to(package)
    target = copied / selected
    original = target.read_bytes()
    primitives.read_path_config(original)
    matches = list(re.finditer(rb"(?m)^(LTANAPATH[ \t]*=[ \t]*)([^\r\n]*?)([ \t]*)(\r?$)", original))
    require(len(matches) == 1, "ltsep_ltanapath_assignment_not_unique")
    match = matches[0]
    patched = original[:match.start(2)] + str(worktree).encode() + original[match.end(2):]
    target.write_bytes(patched)
    expected = dict(before, **{selected.as_posix(): digest(patched)})
    require(tree_bytes(copied) == expected, "overlay_unapproved_change")
    return copied, expected, primitives.analysis_environment(parent)


def compare(records, attempt, worktree, targets):
    """Validate each default caller independently, then compare mode groups."""
    expected = worktree / "OUTPUT/Analysis/KaonLT"
    for row in records:
        base, child = row["baseline"]["paths"], row["overlay"]["paths"]
        for key in FIELDS:
            if key not in {"LTANAPATH", "OUTPATH"}:
                require(base.get(key) == child.get(key), "scientific_input_path_changed:" + row["caller"] + ":" + key)
        row["input_path_identity"] = {}
        for key in primitives.PATH_FIELDS:
            if key in {"USER", "HOST", "ANATYPE", "LTANAPATH", "OUTPATH"}:
                continue
            baseline_real = row["baseline"]["input_realpaths"][key]
            require(str(real(Path(child[key]))) == baseline_real,
                    "scientific_input_realpath_changed:" + row["caller"] + ":" + key)
            row["input_path_identity"][key] = {"path": child[key], "realpath": baseline_real,
                                                "read_only": True, "data_bytes_read": False}
        require(child["ANATYPE"] == "Kaon", "analysis_type_not_kaon")
        require(real(Path(child["LTANAPATH"])) == worktree, "ltanapath_not_isolated")
        try:
            output = real(Path(child["OUTPATH"]))
        except ValueError as error:
            raise ValueError("default_root_unresolved:" + row["caller"] + ":" + str(error)) from error
        prefix = "default_root" if not row["mode"] else "plot_root"
        targets.append({"caller": row["caller"], "mode": row["mode"],
                        "outpath": child["OUTPATH"], "realpath": str(output),
                        "owner": "attempt" if output.is_relative_to(attempt) else "external",
                        "csv_parent": str(output) if not row["mode"] else None})
        require(output.is_relative_to(attempt), prefix + "_outpath_external")
    defaults = {t["realpath"] for t in targets[:3]}
    plots = {t["realpath"] for t in targets[3:]}
    require(len(defaults) == 1, "default_root_caller_mismatch")
    require(len(plots) == 1, "plot_root_caller_mismatch")
    require(defaults == plots, "default_root_outpath_differs_from_plot")
    require(defaults == {str(real(expected))}, "analysis_output_root_mismatch")
    for key in FIELDS:
        if key in {"LTANAPATH", "OUTPATH"}:
            continue
        default_inputs = {r["overlay"]["paths"].get(key) for r in records[:3]}
        plot_inputs = {r["overlay"]["paths"].get(key) for r in records[3:]}
        require(len(default_inputs) == 1, "default_root_input_caller_mismatch:" + key)
        require(len(plot_inputs) == 1 and plot_inputs == default_inputs,
                "input_mode_or_plot_caller_mismatch:" + key)
    return targets


def attest(repo, source_commit):
    """Return resolver-only PASS/BLOCKED, including preservation on failure."""
    receipt = {"schema_version": "pion_chain_ltsep_default_root/v1", "tool_version": 1,
               "baseline_source_sha": BASELINE, "checked_commit": source_commit,
               "passed": False, "failed_invariant": None, "farm_execution_authorized": False,
               "scope_limits": list(LIMITS), "callers": [], "resolved_targets": [],
               "source_sha256": {}, "preservation": {"installed_package": None,
                   "copied_package": None, "ordinary_checkout": None}, "verified_operations": [],
               "unknown_operations": ["real installed-package behavior until directly probed",
                                      "complete launcher write confinement"]}
    before = package = package_before = None

    def block(error):
        if receipt["failed_invariant"] is None:
            receipt["failed_invariant"] = str(error)
        else:
            receipt.setdefault("additional_failures", []).append(str(error))
        receipt["passed"] = False

    try:
        repo = real(Path(repo))
        receipt["source_sha256"] = source_identity(repo, source_commit)
        before = repository_state(repo)
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", GIT_OPTIONAL_LOCKS="0")
        receipt["configured_external_inputs"] = {k: v for k, v in env.items()
            if k.startswith("LT_BG_") or k == "BACKGROUND_SIMC_PATH"}
        baseline = []
        receipt["baseline_callers"] = baseline
        # Record installed identity/inventory before any disposable worktree or
        # copied-config patch, so partial setup failures still verify originals.
        for role, caller, mode, cwd, source in CALLERS:
            record = probe(repo, caller, mode, env, cwd)
            record.update(role=role, source_file=source)
            current_package = module_identity(record)
            if package is None:
                package = current_package
                package_before = tree_bytes(package)
                receipt["installed_package"] = {"path": str(package), "files": package_before}
            require(current_package == package, "installed_package_identity_changed")
            require(record["config_matches"] == baseline[0]["config_matches"] if baseline else True,
                    "default_root_config_mismatch")
            require(real(Path(record["paths"]["LTANAPATH"])) == repo, "baseline_ltanapath_not_ordinary_repo")
            record["input_realpaths"] = {k: str(real(Path(record["paths"][k])))
                for k in primitives.PATH_FIELDS
                if k not in {"USER", "HOST", "ANATYPE", "LTANAPATH", "OUTPATH"}}
            baseline.append(record)
        with primitives.owned_worktree(repo, source_commit, "analysis") as worktree:
            copied, expected, child_env = overlay(repo, worktree, baseline[3], package, package_before)
            receipt["overlay_package"] = {"path": str(copied), "files": expected,
                                           "only_changed_field": "LTANAPATH"}
            try:
                for (role, caller, mode, cwd, source), base in zip(CALLERS, baseline):
                    record = probe(worktree, caller, mode, child_env, cwd)
                    record.update(role=role, source_file=source)
                    require(module_identity(record) == copied, "copied_import_identity_mismatch")
                    for key in ("package_file", "pathing_file", "root_file"):
                        require(Path(record[key]).relative_to(copied) == Path(base[key]).relative_to(package),
                                "copied_module_identity_mismatch:" + key)
                    require(Path(record["config_matches"][0]).relative_to(copied) ==
                            Path(base["config_matches"][0]).relative_to(package), "copied_config_selection_changed")
                    receipt["callers"].append({"role": role, "source_file": source,
                                               "caller": caller, "mode": list(mode),
                                               "baseline": base, "overlay": record})
                compare(receipt["callers"], worktree.parent, worktree, receipt["resolved_targets"])
                receipt["verified_operations"] = ["exact default and Plot resolver calls", "input-path equality",
                    "LTANAPATH-only copied config patch", "resolved CSV/Plot parent containment; no output creation"]
                receipt["unknown_operations"] = ["complete launcher write confinement", "ROOT and science acceptance"]
                receipt["passed"] = True
            except (OSError, RuntimeError, ValueError, KeyError, TypeError, subprocess.SubprocessError) as error:
                block(error)
            finally:
                receipt["preservation"]["copied_package"] = False
                require(tree_bytes(copied) == expected, "copied_package_changed_during_probe")
                receipt["preservation"]["copied_package"] = True
    except (OSError, RuntimeError, ValueError, KeyError, TypeError, subprocess.SubprocessError) as error:
        block(error)
    finally:
        if package_before is not None:
            try:
                receipt["preservation"]["installed_package"] = False
                require(tree_bytes(package) == package_before, "installed_package_changed")
                receipt["preservation"]["installed_package"] = True
            except (OSError, ValueError) as error:
                block(error)
        if before is not None:
            try:
                receipt["preservation"]["ordinary_checkout"] = False
                require(repository_state(repo) == before, "ordinary_checkout_changed")
                receipt["preservation"]["ordinary_checkout"] = True
            except (OSError, ValueError, subprocess.SubprocessError) as error:
                block(error)
    return receipt


def protected_receipt_locations(repo, receipt):
    """Establish exclusions from real probe evidence, never from output strings.

    On incomplete discovery no receipt is published. Include both observed and
    current realpaths so an input alias change cannot erase its earlier target.
    Reading source configuration here does not execute it or open science data.
    """
    package = receipt.get("installed_package", {}).get("path")
    baseline = receipt.get("baseline_callers", [])
    require(package and len(baseline) == len(CALLERS), "receipt_protected_paths_not_verified")
    locations = {real(Path(repo)), real(Path(package))}
    for record in [*baseline, *(r["overlay"] for r in receipt["callers"])]:
        for field in FIELDS:
            if field in {"USER", "HOST", "ANATYPE"}:
                continue
            value = record["paths"].get(field)
            if value is not None:
                locations.add(real(Path(value)))
            observed = record.get("input_realpaths", {}).get(field)
            if observed is not None:
                locations.add(real(Path(observed)))
    for key, value in receipt.get("configured_external_inputs", {}).items():
        if not value or key in {"LT_BG_SAMPLE_Q2", "LT_BG_SAMPLE_W", "LT_BG_SAMPLE_EPSILON"}:
            continue
        path = real(Path(value))
        # Sample BASE is a filename prefix; SAMPLE_BASE is a directory.
        locations.add(path if key in {"BACKGROUND_SIMC_PATH", "LT_BG_SAMPLE_BASE"} else path.parent)
    # set_SymLinks sources this frozen config, which may override the ambient
    # background SIMC path. Accept only its simple literal/${USER} assignment.
    config = Path(repo) / "background_samples/background_samples.conf"
    assignments = [line.strip().partition("=")[2].strip() for line in config.read_text().splitlines()
                   if re.match(r"^\s*BACKGROUND_SIMC_PATH\s*=", line)]
    require(len(assignments) <= 1, "receipt_background_config_not_verified")
    if assignments:
        value = assignments[0]
        if value[:1] in {"\"", "'"}:
            require(len(value) >= 2 and value[-1] == value[0], "receipt_background_config_not_verified")
            value = value[1:-1]
        value = value.replace("${USER}", baseline[0]["paths"]["USER"])
        require(not any(c in value for c in "$`\n\r"), "receipt_background_config_not_verified")
        if value:
            locations.add(real(Path(value)))
    return sorted(locations)


def output_directory(directory, repo, protected=()):
    path = Path(directory)
    require(path.is_absolute() and ".." not in path.parts, "receipt_directory_invalid")
    require(all(not p.is_symlink() for p in (path, *path.parents)), "receipt_directory_symlink")
    require(path.is_dir() and not real(path).is_relative_to(real(repo)), "receipt_directory_not_external")
    resolved = real(path)
    for location in protected:
        require(not resolved.is_relative_to(location), "receipt_directory_protected_input:" + str(location))
    require(hasattr(os, "getuid") and path.stat().st_uid == os.getuid(), "receipt_directory_not_owned")
    require(path.stat().st_mode & 0o022 == 0, "receipt_directory_shared_writable")
    return path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    fd = None
    receipt = None
    try:
        directory = output_directory(args.output_dir, args.repo)
        fd = os.open(directory, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
        identity = os.fstat(fd)
        attempt_id = uuid.uuid4().hex
        name = "pion-ltsep-default-root-" + attempt_id + ".json"
        require(not os.path.lexists(directory / name), "receipt_collision")
        receipt = attest(args.repo, args.source_commit)
        protected = protected_receipt_locations(args.repo, receipt)
        receipt["attempt_id"] = attempt_id
        receipt["status"] = "PASS" if receipt["passed"] else "BLOCKED"
        payload = (json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
        output_directory(directory, args.repo, protected)
        after = directory.stat()
        require((identity.st_dev, identity.st_ino) == (after.st_dev, after.st_ino), "receipt_directory_replaced")
        out = os.open(name, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600, dir_fd=fd)
        with os.fdopen(out, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        print(json.dumps({"receipt": str(directory / name), "sha256": digest(payload),
                          "bytes": len(payload), "status": receipt["status"],
                          "farm_execution_authorized": False}))
        return 0 if receipt["passed"] else 1
    except (OSError, ValueError, RuntimeError) as error:
        first = receipt.get("failed_invariant") if receipt is not None else None
        print(json.dumps({"status": "BLOCKED", "failed_invariant": first or str(error),
                          "receipt_publication_failure": str(error),
                          "farm_execution_authorized": False}), file=sys.stderr)
        return 1
    finally:
        if fd is not None:
            os.close(fd)


if __name__ == "__main__":
    raise SystemExit(main())
