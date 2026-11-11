import argparse
import hashlib
import io
import json
import math
import os
import stat
import subprocess
import sys
import tempfile
import zipfile
from copy import deepcopy
from pathlib import Path

import numpy as np

DELTA = (-10., -7., -4., -1., 2., 5., 8., 11., 14., 17., 20.)
REQUIRED = ("missmass", "t", "phipq", "ssdelta", "ssxptar", "ssyptar",
            "hsdelta", "hsxptar", "hsyptar", "Q2", "W", "phgcer_x_det", "phgcer_y_det")
OPTIONAL = ("ssxfp", "ssyfp", "ssxpfp", "ssypfp", "ssytar",
            "hsxfp", "hsyfp", "hsxpfp", "hsypfp", "hsytar")
CUT_KEYS = ("W", "Q2", "EPSSET", "POL", "ParticleType", "tmin", "tmax",
            "cut_poly", "TMIN_RESOLUTION_THRESHOLD")


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def close(a, b, label):
    np.testing.assert_allclose(a, b, rtol=1e-10, atol=1e-12, err_msg=label)


def transform(parent):
    ed = np.asarray(parent["mm_edges"])
    mid = (ed[:-1] + ed[1:]) / 2
    size = len(ed) + 1
    shift = parent["shift"]
    if abs(shift) <= 1e-12:
        L = np.eye(size)
    else:
        L = np.zeros((size, size))
        for i, x in enumerate(mid-shift, 1):
            if x < mid[0] or x > mid[-1]:
                continue
            if parent["interpolation"] == "nearest":
                L[i, int(np.argmin(abs(mid-x)))+1] = 1.
                continue
            require(parent["interpolation"] == "linear", "unsupported frozen interpolation")
            anchor = max(0, min(len(mid)-1, int(np.searchsorted(ed, x, side="right"))-1))
            if abs(x-mid[anchor]) <= 1e-12:
                L[i, anchor+1] = 1.
                continue
            left = anchor if x > mid[anchor] else anchor-1
            right = left+1
            if left >= 0 and right < len(mid):
                f = (x-mid[left])/(mid[right]-mid[left])
                L[i, left+1] = 1.-f
                L[i, right+1] = f
    raw = np.asarray(parent["raw"]["values"])
    den = (L @ raw)[1:-1].sum() if parent["renormalize"] else raw[1:-1].sum()
    require(den > 0, "nonpositive original parent normalization")
    return L/den


def check_moments(weights, bins, ref, n, label):
    require(len(weights) == ref["count"], label + " count")
    for key, moment in (("signed", weights), ("absolute", np.abs(weights)), ("sumw2", weights**2)):
        close(moment.sum(), ref[key], label + " " + key)
    for key, moment in (("values", weights), ("variances", weights**2)):
        close(np.bincount(bins, weights=moment, minlength=n), ref[key], label + " " + key)


def check_role(rows, role):
    n = len(role["parents"][0]["mm_edges"]) + 1
    w = rows["weight"]
    ti, pi, di = (rows[k] for k in ("t_bin", "phi_bin", "delta_bin"))
    mb = rows["mm_bin"]
    accepted = (ti >= 0) & (pi >= 0)
    for label, sel in (("selected", np.ones(len(w), bool)), ("accepted", accepted), ("excluded", ~accepted)):
        check_moments(w[sel], mb[sel], role[label], n, role["id"] + " " + label)
    receipts, arrays = [], {}
    for t, parent in enumerate(role["parents"]):
        sel = accepted & (ti == t)
        L = transform(parent)
        for phi, child in enumerate(parent["children"]):
            cut = sel & (pi == phi)
            check_moments(w[cut], mb[cut], child, n, f"{role['id']} t{t+1} phi{phi+1}")
        check_moments(w[sel], mb[sel], parent["raw"], n, f"{role['id']} t{t+1}")
        # Last entry explicitly preserves delta flow even though the original cuts
        # constrain all admitted rows to the closed [-10, 20] interval.
        part = np.where(di[sel] < 0, 10, di[sel])
        idx = part*n + mb[sel]
        S = np.bincount(idx, weights=w[sel], minlength=11*n).reshape(11,n)
        V = np.bincount(idx, weights=w[sel]**2, minlength=11*n).reshape(11,n)
        C = np.bincount(idx, minlength=11*n).reshape(11,n)
        aligned = S @ L.T
        aligned_var = V @ (L**2).T
        close(aligned.sum(0), parent["aligned"], "aligned parent closure")
        close(aligned_var.sum(0), parent["aligned_sumw2"], "aligned parent Sumw2 closure")
        for key, value in (("raw_sumw", S), ("raw_sumw2", V), ("count", C),
                           ("aligned", aligned), ("aligned_sumw2", aligned_var), ("alignment_operator", L)):
            arrays[f"t{t+1}_{key}"] = value
        receipts.append({"t_bin": t+1, "raw_entries": int(C.sum()), "delta_counts": C.sum(1).tolist(),
                         "delta_raw_sumw": S.sum(1).tolist(), "delta_raw_sumw2": V.sum(1).tolist(),
                         "delta_flow_count": int(C[-1].sum()), "mm_underflow": int(C[:,0].sum()),
                         "mm_overflow": int(C[:,-1].sum()),
                         "max_aligned_error": float(np.max(abs(aligned.sum(0)-parent["aligned"])))})
    return receipts, arrays


def find_input(role, roots, regenerated_root=None):
    if regenerated_root is not None:
        historical = Path(role["source"]["path"])
        candidate = regenerated_root / historical.parent.name / historical.name
        require(candidate.is_file(), "regenerated input missing: " + str(candidate))
        return candidate.absolute()
    original = Path(role["source"]["path"])
    attempts = []
    for path in dict.fromkeys([original, *(r / original.parent.name / original.name for r in roots)]):
        if path.is_file():
            size = path.stat().st_size
            if size == role["source"]["bytes"] and digest(path) == role["source"]["sha256"]:
                return path.resolve()
            attempts.append(str(path) + " [size or hash mismatch]")
        else:
            attempts.append(str(path) + " [absent]")
    raise RuntimeError(role["id"] + " exact input unresolved: " + "; ".join(attempts))


def verify_source(ref, repo):
    for name, expected in ref["source_hashes"].items():
        source = (repo/name).read_bytes().replace(bytes((13,10)),bytes((10,)))
        require(hashlib.sha256(source).hexdigest() == expected, "source bytes differ: " + name)
    cmd = ["git", "-C", str(repo), "rev-parse", "HEAD"]
    return subprocess.check_output(cmd, text=True).strip()


def open_tree(ROOT, role, path, regenerated=False):
    rootfile = ROOT.TFile.Open(str(path), "READ")
    require(bool(rootfile) and not rootfile.IsZombie(), "cannot read ROOT " + str(path))
    try:
        tree = rootfile.Get(role["source"]["tree_name"])
        require(bool(tree) and tree.GetEntries() > 0, "missing or empty tree: " + role["id"])
        if not regenerated:
            require(tree.GetEntries() == role["source"]["tree_entries"],
                    "tree entries changed: " + role["id"])
        has = lambda k: bool(tree.GetBranch(k) or tree.GetLeaf(k))
        missing = [k for k in REQUIRED if not has(k)]
        require(not missing, role["id"] + " missing branches " + str(missing))
        if not regenerated:
            require(has("iter_weight") == role["source"]["weight_branch_present"], "weight branch mismatch")
        return rootfile, tree, [k for k in OPTIONAL if has(k)]
    except BaseException:
        rootfile.Close()
        raise


def read_rows(ROOT, cuts, bin_for, role, path, admitted, regenerated=False):
    state = role["source"]["active_selection_state"]
    cuts.set_val(dict(state, cut_mode="poly", poly_points=state["cut_poly"]))
    require(all(getattr(cuts,k) == state[k] for k in CUT_KEYS if k != "cut_poly"),
            "cut-state mismatch: " + role["id"])
    require([list(p) for p in cuts.cut_poly] == state["cut_poly"], "frozen diamond mismatch")
    points = role["source"]["hgcer_hole_points"]
    require(len(points) >= 3, "missing frozen hole")
    hole = ROOT.TCutG("delta_projection_" + role["id"], len(points))
    for j, (x, y) in enumerate(points):
        hole.SetPoint(j, x, y)
    f, tree, optional = open_tree(ROOT, role, path, regenerated)
    fields = ("entry", "weight", "minus_t", "phi_degrees", "delta_bin",
              "t_bin", "phi_bin", "mm_bin") + REQUIRED + tuple(optional)
    values = {k: [] for k in fields}
    rejected = {"cuts": 0, "hgcer_hole": 0}
    try:
        for entry, evt in enumerate(tree):
            if not cuts.apply_simc_sub_cuts(evt):
                rejected["cuts"] += 1
                continue
            if hole.IsInside(evt.phgcer_x_det, evt.phgcer_y_det):
                rejected["hgcer_hole"] += 1
                continue
            x = {k: float(getattr(evt,k)) for k in REQUIRED + tuple(optional)}
            x.update(entry=entry, weight=float(getattr(evt,"iter_weight",1.0)),
                     minus_t=-x["t"], phi_degrees=x["phipq"]*(180.0/math.pi))
            require(all(math.isfinite(x[k]) for k in REQUIRED + ("weight",)), "nonfinite selected row")
            x.update(delta_bin=bin_for(x["ssdelta"], DELTA),
                     t_bin=bin_for(x["minus_t"], role["t_edges"]),
                     phi_bin=bin_for(x["phi_degrees"], role["phi_edges"]),
                     mm_bin=np.searchsorted(role["parents"][0]["mm_edges"],x["missmass"],side="right"))
            for k in fields:
                values[k].append(x[k])
    finally:
        f.Close()
    rows = {k: np.array(v,dtype=np.int64 if k=="entry" or k.endswith("_bin") else np.float64)
            for k,v in values.items()}
    require(len(rows["entry"])+sum(rejected.values()) == admitted["tree_entries"], "tree accounting")
    require(input_identity(path) == admitted["identity"], "input changed during read " + role["id"])
    return rows, rejected


def input_identity(path):
    before = path.stat()
    result = {"requested_path": str(path.absolute()), "resolved": str(path.resolve(strict=True)),
              "sha256": digest(path), "bytes": before.st_size,
              "device": before.st_dev, "inode": before.st_ino, "mtime_ns": before.st_mtime_ns}
    after = path.stat()
    require((before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) ==
            (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns),
            "input changed during hash: " + str(path))
    require(str(path.resolve(strict=True)) == result["resolved"], "input alias changed during hash")
    return result


def validate_reference(ref):
    expected = {setting+"_"+channel for setting in
                ("Center_highe", "Center_lowe", "Left_lowe", "Left_highe", "Right_highe")
                for channel in ("pi_n", "pi_delta", "pi_sidis")}
    require(tuple(ref["delta_edges"]) == DELTA and len(ref["roles"]) == 15, "reference inventory")
    require({r["id"] for r in ref["roles"]} == expected, "reference role identities")
    for role in ref["roles"]:
        require(len(role["parents"]) == 3 and len(role["t_edges"]) == 4 and
                len(role["phi_edges"]) == 10, "reference canonical inventory")
        for edges in (role["t_edges"], role["phi_edges"],
                      *(p["mm_edges"] for p in role["parents"])):
            require(all(math.isfinite(x) for x in edges) and
                    all(a < b for a, b in zip(edges, edges[1:])), "reference edges")
        require(all(p["mm_edges"] == role["parents"][0]["mm_edges"]
                    for p in role["parents"]), "parent MM axes differ")


def regenerated_role(rows, historical, admitted):
    from pion_t_template_comparison import integrate_raw_templates
    source = dict(historical["source"], status="available",
                  path=admitted["identity"]["requested_path"],
                  sha256=admitted["identity"]["sha256"], bytes=admitted["identity"]["bytes"],
                  tree_entries=admitted["tree_entries"],
                  weight_branch_present=admitted["weight_branch_present"])
    coordinates = zip(rows["missmass"], rows["minus_t"],
                      rows["phi_degrees"], rows["weight"])
    payload = integrate_raw_templates(coordinates, historical["t_edges"],
                                     historical["phi_edges"],
                                     historical["parents"][0]["mm_edges"], source)
    require(payload["invalid_records"] == 0 and payload["canonical_selection_closure"],
            "regenerated canonical closure")
    role = deepcopy(historical)
    role["source"] = source
    for key in ("selected", "accepted", "excluded"):
        role[key] = payload[key]
    differences = []
    for j, (parent, evidence) in enumerate(zip(role["parents"], payload["parents"])):
        require(evidence["status"] == "genuinely_local" and
                evidence["weighted_child_to_parent_closure"], "regenerated parent unavailable")
        parent["raw"] = evidence["raw"]
        parent["children"] = evidence["raw_children"]
        L = transform(parent)
        parent["aligned"] = (L @ np.asarray(parent["raw"]["values"])).tolist()
        parent["aligned_sumw2"] = (L**2 @ np.asarray(parent["raw"]["variances"])).tolist()
        original = historical["parents"][j]
        differences.append({"t_bin": j+1,
                            "count_difference": parent["raw"]["count"]-original["raw"]["count"],
                            "raw_values_difference": (np.asarray(parent["raw"]["values"]) -
                                                       original["raw"]["values"]).tolist(),
                            "aligned_values_difference": (np.asarray(parent["aligned"]) -
                                                           original["aligned"]).tolist(),
                            "historical_equality_required": False})
    return role, differences


def npz_bytes(arrays):
    buf = io.BytesIO()
    np.savez_compressed(buf, **arrays)
    return buf.getvalue()


def publish_projection(output, output_data, receipt):
    published = False
    pending_identity = None
    try:
        with tempfile.TemporaryDirectory(prefix=".kaonlt_delta_", dir=output.parent) as temp:
            pending = Path(temp)/"projection.zip"
            with zipfile.ZipFile(pending, "x", zipfile.ZIP_DEFLATED) as archive:
                for name, payload in output_data.items():
                    archive.writestr(name, payload)
                archive.writestr("receipt.json", json.dumps(receipt, indent=2, allow_nan=False))
            with zipfile.ZipFile(pending) as archive:
                require(archive.testzip() is None, "ZIP integrity")
            info = pending.lstat()
            pending_identity = (info.st_dev, info.st_ino)
            result_hash = digest(pending)
            os.link(pending, output)
            published = True
        final = output.lstat()
        require(stat.S_ISREG(final.st_mode) and
                (final.st_dev, final.st_ino) == pending_identity,
                "published projection replaced")
        return result_hash
    except BaseException:
        if published:
            try:
                final = output.lstat()
                if stat.S_ISREG(final.st_mode) and (final.st_dev, final.st_ino) == pending_identity:
                    output.unlink()
            except FileNotFoundError:
                pass
            except OSError as cleanup_error:
                print("projection rollback failed:", cleanup_error, file=sys.stderr)
        raise


def main():
    cli = argparse.ArgumentParser(description="Read existing h10 and verify ten-delta parent closure.")
    cli.add_argument("--repo", type=Path, required=True)
    cli.add_argument("--input-root", type=Path, action="append", default=[])
    cli.add_argument("--regenerated-input-root", type=Path,
                     help="Explicit new cohort root; neutron/delta/sidis subdirectories.")
    cli.add_argument("--reference", type=Path, default=Path(__file__).with_name("reference.json"))
    cli.add_argument("--output", type=Path, required=True)
    args = cli.parse_args()
    repo, output = args.repo.resolve(), args.output.absolute()
    require(not os.path.lexists(output), "output already exists")
    require(output.parent.is_dir(), "output parent must already exist")
    ref_path = args.reference.resolve()
    reference_bytes = ref_path.read_bytes()
    reference_hash = hashlib.sha256(reference_bytes).hexdigest()
    ref = json.loads(reference_bytes)
    validate_reference(ref)
    regenerated = args.regenerated_input_root is not None
    require(not (regenerated and args.input_root), "choose one input admission mode")
    head = verify_source(ref, repo)
    for rel in ("src/cuts", "src/binning", "src/utility"):
        sys.path.insert(0, str(repo/rel))
    import ROOT
    import apply_cuts as cuts
    from canonical_binning import find_canonical_bin
    import pion_t_template_comparison as integration
    for module, path in ((cuts, "src/cuts/apply_cuts.py"),
                         (sys.modules["canonical_binning"], "src/utility/canonical_binning.py"),
                         (sys.modules["theta_cm"], "src/binning/theta_cm.py"),
                         (sys.modules["data_coordinates"], "src/cuts/data_coordinates.py"),
                         (integration, "src/cuts/pion_t_template_comparison.py")):
        require(Path(module.__file__).resolve()==(repo/path).resolve(), "source import alias")
    ROOT.gROOT.SetBatch(True)
    resolved = []
    admitted_inputs = []
    physical_inputs = set()
    for role in ref["roles"]:
        p = find_input(role, args.input_root, args.regenerated_input_root)
        identity = input_identity(p)
        require(Path(identity["resolved"]).name == p.name, "input role filename alias")
        physical = (identity["device"], identity["inode"])
        require(physical not in physical_inputs, "duplicate physical input")
        physical_inputs.add(physical)
        f, tree, optional = open_tree(ROOT, role, p, regenerated)
        try:
            item = {"identity": identity, "tree_entries": int(tree.GetEntries()),
                    "tree_name": role["source"]["tree_name"],
                    "weight_branch_present": bool(tree.GetBranch("iter_weight") or tree.GetLeaf("iter_weight")),
                    "required_branches": list(REQUIRED), "optional_branches": optional,
                    "available_branches": [branch.GetName() for branch in tree.GetListOfBranches()]}
        finally:
            f.Close()
        require(input_identity(p) == identity, "input changed during admission")
        admitted_inputs.append(item)
        resolved.append(p)
        print("INPUT VERIFIED", role["id"], p, flush=True)
    require(len(set(resolved))==15, "unexpected duplicate input path")
    output_data = {}
    receipt = {"status":"verified_projection", "runtime_head":head, "root":ROOT.gROOT.GetVersion(),
               "python":sys.version.split()[0], "numpy":np.__version__, "reference_sha256":reference_hash,
               "primary_reference_sha256":ref["primary_zip_sha256"], "delta_edges":DELTA,
               "delta_flow_index":10, "roles":[],
               "cohort_kind": "regenerated" if regenerated else "historical",
               "parent_closure_reference": "current_inputs" if regenerated else "historical_reference",
               "historical_template_equivalent": False if regenerated else True,
               "scientifically_accepted": False}
    for role, path, admitted in zip(ref["roles"], resolved, admitted_inputs):
        rows, rejected = read_rows(ROOT, cuts, find_canonical_bin, role, path, admitted, regenerated)
        current, comparison = (regenerated_role(rows, role, admitted) if regenerated
                               else (role, []))
        checks, parts = check_role(rows, current)
        receipt["roles"].append({"id":role["id"],"actual_path":str(path),"sha256":admitted["identity"]["sha256"],
                                 "input":admitted, "historical_source":role["source"],
                                 "historical_comparison":comparison, "closure_reference":current,
                                 "tree":role["source"]["tree_name"],"rejected":rejected,"parents":checks})
        output_data[role["id"] + ".npz"] = npz_bytes(dict(rows,**parts))
        print("PARENT CLOSURE", role["id"], len(checks), flush=True)
    for path, admitted in zip(resolved, admitted_inputs):
        require(input_identity(path) == admitted["identity"], "input changed before publication")
    require(digest(ref_path) == reference_hash, "historical reference changed during extraction")
    require(verify_source(ref, repo) == head, "repository HEAD changed during extraction")
    receipt["parent_channel_closures"] = sum(len(r["parents"]) for r in receipt["roles"])
    require(receipt["parent_channel_closures"]==45, "incomplete parent closure")
    receipt["members"] = {n:{"bytes":len(b),"sha256":hashlib.sha256(b).hexdigest()} for n,b in output_data.items()}
    output_hash = publish_projection(output, output_data, receipt)
    print("PASS 45 parent-channel closures", output, output_hash, flush=True)


if __name__ == "__main__":
    main()
