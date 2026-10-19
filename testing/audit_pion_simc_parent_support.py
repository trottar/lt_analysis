"""Detached existing-SIMC support measurement; never a production correction.

CLI: --config PATH --config-sha256 SHA --settings Center Left --report NEW.json
The explicit config is the resolved production inpDict *after* DiamondPlot,
not the early zeroth-iteration input bundle (which lacks the polygon). It must
carry analysis_runtime_config/hash, source_head, and the original input roles.
No config, interval, cache, or ROOT input is written. Missing config provenance
is a prerequisite failure, never an invitation to reconstruct cuts by fitting.
Optional --center-evidence is the existing seven-stage JSON, not a template.

Exit 0 means arithmetic/production-loader parity only, 1 a failed measurement,
2 invalid invocation. Neither occupancy nor closure proves scientific adequacy.
All ROOT/farm imports are inside scan(); the numerical API is standard library.
"""
from __future__ import annotations

import argparse
import bisect
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

COMPONENTS = ("pi_n", "pi_sidis", "pi_delta")
MOMENTS = ("count", "sum_w", "sum_abs_w", "sum_w2")
SCHEMA = "pion_simc_parent_support/v1"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def encoded(value):
    return json.dumps(value, sort_keys=True, allow_nan=False, separators=(",", ":")).encode()


def file_identity(path):
    path = Path(path).resolve(strict=True)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def edges(values):
    values = list(map(float, values))
    require(len(values) >= 2 and all(map(math.isfinite, values))
            and all(a < b for a, b in zip(values, values[1:])), "invalid axis")
    return values


def bin_index(value, axis):
    """Same left-closed/right-open rule as binning_helpers.find_bin_index."""
    if not math.isfinite(value) or value < axis[0] or value >= axis[-1]:
        return None
    return bisect.bisect_right(axis, value) - 1


def empty_moments():
    return dict.fromkeys(MOMENTS, 0)


def add(moment, weight):
    moment["count"] += 1
    moment["sum_w"] += weight
    moment["sum_abs_w"] += abs(weight)
    moment["sum_w2"] += weight * weight
    require(all(math.isfinite(moment[k]) for k in MOMENTS), "moment overflow")


def support(moment):
    result = dict(moment)
    w, a, q = (moment[k] for k in ("sum_w", "sum_abs_w", "sum_w2"))
    result.update(neff_signed=(w / math.sqrt(q)) ** 2 if q > 0 else None,
                  neff_absolute=(a / math.sqrt(q)) ** 2 if q > 0 else None,
                  cancellation_fraction=max(0., 1. - abs(w) / a) if a > 0 else None)
    return result


class Histogram:
    """Raw weighted moments, including flows and event-defined regions."""
    def __init__(self, axis, regions=None):
        self.edges = edges(axis)
        self.bins = [empty_moments() for _ in range(len(self.edges) + 1)]
        self.regions = regions or {}
        for windows in self.regions.values():
            for lo, hi in windows:
                require(math.isfinite(lo) and math.isfinite(hi) and lo <= hi,
                        "invalid region window")
        self.region_moments = {name: empty_moments() for name in self.regions}
        self.region_bins = {name: [empty_moments() for _ in self.bins] for name in self.regions}

    def fill(self, mm, weight):
        require(math.isfinite(mm) and math.isfinite(weight)
                and math.isfinite(weight * weight), "invalid fill")
        index = (0 if mm < self.edges[0] else len(self.edges)
                 if mm >= self.edges[-1] else bin_index(mm, self.edges) + 1)
        add(self.bins[index], weight)
        for name, windows in self.regions.items():
            if any(lo <= mm <= hi for lo, hi in windows):
                add(self.region_moments[name], weight)
                add(self.region_bins[name][index], weight)

    def report(self):
        total = {k: math.fsum(b[k] for b in self.bins) for k in MOMENTS}
        integral = math.fsum(b["sum_w"] for b in self.bins[1:-1])
        meaningful = math.isfinite(integral) and integral > 0
        return {"edges": self.edges, "total": support(total),
                "normal_bins": [support(b) for b in self.bins[1:-1]],
                "underflow": support(self.bins[0]), "overflow": support(self.bins[-1]),
                "normalization_integral": integral, "normalization_meaningful": meaningful,
                "normalized_content": [b["sum_w"] / integral for b in self.bins[1:-1]]
                if meaningful else None,
                "regions": {k: dict(support(v), normal_bins=[support(b) for b in self.region_bins[k][1:-1]],
                                    underflow=support(self.region_bins[k][0]),
                                    overflow=support(self.region_bins[k][-1]))
                            for k, v in self.region_moments.items()},
                "region_rule": "event MM inclusive bounds; raw coordinates; no data-MM shift",
                "status": "NOT VERIFIED", "scientific_adequacy": "NOT VERIFIED"}


def closure(whole, pieces):
    require(all(p.edges == whole.edges and p.regions == whole.regions for p in pieces),
            "closure axes/regions mismatch")
    residuals = []
    rows = list(zip(whole.bins, zip(*(p.bins for p in pieces))))
    rows += [(whole.region_moments[k], [p.region_moments[k] for p in pieces])
             for k in whole.regions]
    for name in whole.regions:
        rows += list(zip(whole.region_bins[name], zip(*(p.region_bins[name] for p in pieces))))
    for target, terms in rows:
        row = {}
        for key in MOMENTS:
            value = math.fsum(p[key] for p in terms)
            row[key] = target[key] - value
            if key == "count":
                require(row[key] == 0, "bookkeeping closure failed: count")
            scale = max(1., target["sum_abs_w"], value if key != "sum_w" else
                        math.fsum(abs(p[key]) for p in terms))
            require(abs(row[key]) <= 1e-10 * scale, "bookkeeping closure failed: " + key)
        residuals.append(row)
    return {"status": "SOURCE VERIFIED", "meaning": "arithmetic only, flows and regions included",
            "bin_then_region_residuals": residuals}


def shape_comparison(left, right, *, shared=False):
    """Delta-method weighted variances include normalization covariance.

    Inclusive-setting comparisons share the parent's events. Return that
    covariance explicitly. Disjoint complement variances may be added; no
    p-values or adequacy thresholds are inferred from these sparse shapes.
    """
    require(left.edges == right.edges, "shape axes mismatch")
    a, b = left.report(), right.report()
    if not a["normalization_meaningful"] or not b["normalization_meaningful"]:
        return {"status": "NOT VERIFIED", "reason": "nonpositive normal-bin integral"}
    sa, sb = a["normalization_integral"], b["normalization_integral"]
    qa = math.fsum(x["sum_w2"] for x in left.bins[1:-1])
    qb = math.fsum(x["sum_w2"] for x in right.bins[1:-1])
    rows = []
    for x, y in zip(left.bins[1:-1], right.bins[1:-1]):
        p, q = x["sum_w"] / sa, y["sum_w"] / sb
        va = ((1 - p) ** 2 * x["sum_w2"] + p * p * (qa - x["sum_w2"])) / sa ** 2
        vb = ((1 - q) ** 2 * y["sum_w2"] + q * q * (qb - y["sum_w2"])) / sb ** 2
        cov = (((1 - p) * (1 - q) * x["sum_w2"] + p * q * (qa - x["sum_w2"]))
               / (sa * sb)) if shared else 0.
        rows.append({"difference": p - q, "left_variance": va, "right_variance": vb,
                     "shared_covariance": cov, "difference_variance": max(0., va + vb - 2 * cov)})
    return {"status": "SOURCE VERIFIED", "bins": rows, "p_value": None,
            "caveat": "parent is a subset of inclusive setting; correlated bins and shared events"
            if shared else "disjoint events; normalized bins remain correlated",
            "variance_model": "independent event sumw2, delta method; template systematics unresolved",
            "scientific_adequacy": "NOT VERIFIED"}


class Population:
    def __init__(self, mm_edges, t_edges, phi_edges=None, regions=None):
        self.t_edges = edges(t_edges)
        self.phi_edges = edges(phi_edges) if phi_edges is not None else None
        self.make = lambda: Histogram(mm_edges, regions)
        self.setting = self.make()
        self.parents = [self.make() for _ in self.t_edges[:-1]]
        self.unassigned = self.make()
        self.complements = [self.make() for _ in self.parents]
        self.children = ([[self.make() for _ in self.phi_edges[:-1]] for _ in self.parents]
                         if self.phi_edges is not None else [])
        self.phi_unassigned = [self.make() for _ in self.parents]
        self.coverage = {key: None for key in ("t_min", "t_max", "phi_min", "phi_max")}
        self.invalid_t = self.invalid_phi = 0

    def fill(self, mm, t, weight, phi=None, finder=None):
        self.setting.fill(mm, weight)
        ti = bin_index(t, self.t_edges)
        pi = None
        if self.phi_edges is not None:
            require(phi is not None, "phi missing in requested child population")
            if finder is not None:
                ti, pi = finder(t, phi, self.t_edges, self.phi_edges)
            elif ti is not None:
                pi = bin_index(phi, self.phi_edges)
            self.invalid_phi += not math.isfinite(phi)
        self.invalid_t += not math.isfinite(t)
        for label, value in (("t", t), ("phi", phi)):
            if value is not None and math.isfinite(value):
                for suffix, fn in (("min", min), ("max", max)):
                    key = label + "_" + suffix
                    self.coverage[key] = value if self.coverage[key] is None else fn(self.coverage[key], value)
        (self.unassigned if ti is None else self.parents[ti]).fill(mm, weight)
        for i, complement in enumerate(self.complements):
            if i != ti:
                complement.fill(mm, weight)
        if ti is not None and self.phi_edges is not None:
            (self.phi_unassigned[ti] if pi is None else self.children[ti][pi]).fill(mm, weight)

    def report(self):
        parents = []
        for i, parent in enumerate(self.parents):
            item = {"t_index": i, "t_edges": self.t_edges[i:i + 2], "histogram": parent.report(),
                    "disjoint_complement": self.complements[i].report(),
                    "parent_vs_setting": shape_comparison(parent, self.setting, shared=True),
                    "parent_vs_complement": shape_comparison(parent, self.complements[i]),
                    "complement_closure": closure(self.setting, [parent, self.complements[i]])}
            if self.children:
                item.update(children=[{"phi_index": j, "phi_edges": self.phi_edges[j:j + 2],
                                       "histogram": h.report()} for j, h in enumerate(self.children[i])],
                            phi_unassigned=self.phi_unassigned[i].report(),
                            child_closure=closure(parent, self.children[i] + [self.phi_unassigned[i]]))
            parents.append(item)
        return {"setting": self.setting.report(), "parents": parents,
                "unassigned_t": self.unassigned.report(),
                "closure": closure(self.setting, self.parents + [self.unassigned]),
                "coverage": self.coverage, "invalid_t": self.invalid_t, "invalid_phi": self.invalid_phi,
                "phi_status": "SOURCE VERIFIED" if self.children else "NOT VERIFIED",
                "coverage_limit": "rightmost t/phi edge excluded; unassigned support retained"}


def validate_config(config, source_head):
    required = ("Q2", "W", "EPSSET", "POL", "ParticleType", "tmin", "tmax", "mm_min", "mm_max",
                "NumtBins", "NumPhiBins", "cut_mode", "poly_points", "analysis_runtime_config",
                "analysis_runtime_config_hash", "source_head")
    require(all(k in config for k in required), "missing resolved production configuration/provenance")
    require(config["source_head"] == source_head, "configuration source HEAD mismatch")
    require(config["Q2"] == "4p4" and config["W"] == "2p74" and config["EPSSET"] == "low"
            and config["ParticleType"] == "kaon", "outside contract setting scope")
    require(config["cut_mode"] == "poly" and len(config["poly_points"]) >= 3, "missing production polygon")
    for point in config["poly_points"]:
        require(len(point) == 2 and all(math.isfinite(float(v)) for v in point), "invalid polygon")
    runtime = config["analysis_runtime_config"]
    require(runtime.get("config_hash") == config["analysis_runtime_config_hash"]
            and bool(config["analysis_runtime_config_hash"]), "runtime config identity mismatch")
    for lo, hi in (("tmin", "tmax"), ("mm_min", "mm_max")):
        edges([config[lo], config[hi]])


def validate_branches(available, *, phi):
    required = {"missmass", "t", "ssdelta", "ssxptar", "ssyptar", "hsdelta", "hsxptar",
                "hsyptar", "Q2", "W", "phgcer_x_det", "phgcer_y_det"}
    if phi:
        required.add("phipq")
    require(not (required - set(available)), "missing branches: " + ", ".join(sorted(required - set(available))))


def production_parity(raw, histogram, diagnostic):
    """Compare against the actual source loader, including flows and Sumw2."""
    require(not diagnostic.get("fallback_used"), "production loader fallback; parity unavailable")
    require(len(raw.edges) == histogram.GetNbinsX() + 1, "production axis size mismatch")
    for i, edge in enumerate(raw.edges, 1):
        require(math.isclose(edge, histogram.GetXaxis().GetBinLowEdge(i), abs_tol=1e-12, rel_tol=1e-12),
                "production axis mismatch")
    integral = raw.report()["normalization_integral"]
    require(integral > 0 and diagnostic.get("normalized"), "production shape not normalizable")
    require(math.isclose(integral, diagnostic["weighted_integral_before_norm"], rel_tol=1e-10, abs_tol=1e-12),
            "production raw integral mismatch")
    for i, moment in enumerate(raw.bins):
        # Production TH1D and independent raw sums differ only in summation
        # rounding. Also perform exact ROOT-against-ROOT parity below.
        require(math.isclose(moment["sum_w"] / integral, histogram.GetBinContent(i),
                             rel_tol=1e-10, abs_tol=1e-12), "production shape mismatch")
        require(math.isclose(moment["sum_w2"] / integral ** 2, histogram.GetBinError(i) ** 2,
                             rel_tol=1e-10, abs_tol=1e-12), "production Sumw2 mismatch")
    return {"status": "SOURCE VERIFIED", "meaning": "same-input loader parity; no physics acceptance",
            "arithmetic_tolerance": {"relative": 1e-10, "absolute": 1e-12}}


def cross_check_evidence(document, identities):
    """Compare the captured Center parent identities without guessing JSON layout."""
    found = {}
    def walk(node):
        if isinstance(node, dict):
            if "pion_parent_id" in node and "t_edges" in node:
                key = node["pion_parent_id"]
                found[key] = node["t_edges"]
            for value in node.values():
                walk(value)
        elif isinstance(node, list):
            for value in node:
                walk(value)
    walk(document)
    require(found, "Center evidence has no parent identities")
    for identity in identities:
        key = identity["pion_parent_id"]
        require(key in found and list(found[key]) == identity["t_edges"], "Center evidence parent mismatch")
    return {"status": "SOURCE VERIFIED", "meaning": "parent ID/edge agreement only; artifact acceptance unverified"}


def accumulate_events(tree, config, full, mmcut, hole, cuts, subcuts, finder, mirror):
    """Inject source-owned selectors; numerical tests do not emulate ROOT cuts."""
    counters = dict.fromkeys(("seen", "invalid_weight", "nonfinite_weight",
                             "selected_invalid_mm", "rejected_cuts_or_hole"), 0)
    for event in tree:
        counters["seen"] += 1
        try:
            weight = float(getattr(event, "iter_weight", 1.))
        except (TypeError, ValueError):
            counters["invalid_weight"] += 1
            continue
        if not math.isfinite(weight):
            counters["nonfinite_weight"] += 1
            continue
        require(math.isfinite(weight * weight), "weight squared overflow; measurement invalid")
        allcuts = cuts(event, config["mm_min"], config["mm_max"])
        nommcuts = subcuts(event, config["mm_min"], config["mm_max"])
        outside = not hole.IsInside(event.phgcer_x_det, event.phgcer_y_det)
        if not outside or not nommcuts:
            counters["rejected_cuts_or_hole"] += 1
            continue
        mm, t = float(event.missmass), -float(event.t)
        if not math.isfinite(mm):
            counters["selected_invalid_mm"] += 1
            continue
        phi = float(event.phipq) * (180. / math.pi) if full.phi_edges is not None else None
        full.fill(mm, t, weight, phi, finder)
        mirror.Fill(mm, weight)
        if allcuts:
            mmcut.fill(mm, t, weight, phi, finder)
    return counters


def scan(config, settings, children=False, center_evidence=None):
    """Farm-only read-only adapter; never import main or invoke a fitting owner."""
    root = Path(__file__).resolve().parents[1]
    source_head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    validate_config(config, source_head)
    for part in ("utility", "cuts", "binning"):
        sys.path.insert(0, str(root / "src" / part))
    import ROOT
    import background_config as bg
    import find_bins
    from apply_cuts import set_val, apply_simc_cuts, apply_simc_sub_cuts
    from binning_helpers import find_2d_bin_indices
    from hgcer_hole import apply_HGCer_hole_cut, _hgcer_hole_path
    from pion_component_shapes import load_pion_simc_component_shape, _make_component_hist
    from pion_component_subtraction import build_t_bin_pion_parent_identity

    require(bg.resolve_particle_subtraction_mode(config) == "simc_shape_components", "wrong production mode")
    require(config.get("bg_active_profile") == "no_empirical_residual", "active profile missing/mismatch")
    require(config["analysis_runtime_config"] == bg.resolve_analysis_runtime_config(config["Q2"], config["W"]),
            "runtime configuration differs from current authoritative configuration")
    tcheck = find_bins._validate_authoritative_t_interval(config, "low")
    require(tcheck["valid"], "canonical t identity NOT VERIFIED: " + str(tcheck["validation_rejection_reasons"]))
    pcheck = find_bins._validate_authoritative_phi_interval(config, "low")
    pair = find_bins._validate_canonical_interval_pair(tcheck, pcheck)
    require(pcheck["valid"] and pair["valid"], "canonical pair identity NOT VERIFIED")
    t_axis, phi_axis = edges(tcheck["t_edges"]), edges(pcheck["phi_edges"])
    config = dict(config, canonical_t_binning=dict(tcheck["metadata"], phi_edges=phi_axis))
    axis_files = [file_identity(p) for p in find_bins._interval_paths(config).values()]
    n = int(config.get("bg_opt_mm_plot_nbins", bg.BG_OPT_MM_PLOT_NBINS))
    lo = float(config.get("bg_opt_mm_plot_min", bg.BG_OPT_MM_PLOT_MIN))
    hi = float(config.get("bg_opt_mm_plot_max", bg.BG_OPT_MM_PLOT_MAX))
    require(n > 0, "invalid MM bin count")
    mm_axis = edges([lo + (hi - lo) * i / n for i in range(n + 1)])
    report = {"schema": SCHEMA, "source_HEAD": source_head,
              "source_files": [file_identity(path) for path in
                               [Path(__file__), Path(bg.__file__), Path(find_bins.__file__),
                                root / "src/cuts/apply_cuts.py", root / "src/cuts/hgcer_hole.py",
                                root / "src/utility/pion_component_shapes.py",
                                root / "src/cuts/pion_component_subtraction.py",
                                root / "src/binning/binning_helpers.py"]],
              "worktree_inventory": subprocess.check_output(
                  ["git", "status", "--porcelain=v2", "--untracked-files=all"], cwd=root, text=True),
              "config_digest": hashlib.sha256(encoded(config)).hexdigest(),
              "canonical_inputs": axis_files, "canonical_pair": pair,
              "status": "NOT VERIFIED", "scientific_adequacy": "NOT VERIFIED", "settings": {},
              "limits": ["SIMC support only; no fit-data or subtraction validation",
                         "template systematics and scientific adequacy unresolved",
                         "config authenticity requires independent same-run provenance review"]}
    protected = list(axis_files)
    for setting in settings:
        set_val(config)
        identities = [dict(build_t_bin_pion_parent_identity(config, setting, i, t_axis[i:i + 2]),
                           t_edges=t_axis[i:i + 2]) for i in range(len(t_axis) - 1)]
        entry = {"status": "NOT VERIFIED", "parents": identities, "components": {}}
        report["settings"][setting] = entry
        if setting == "Center" and center_evidence is not None:
            entry["optional_evidence_check"] = cross_check_evidence(center_evidence, identities)
        files = bg.resolve_simc_pion_component_files(config, setting)
        tree_name = bg.resolve_simc_tree_name(config)
        hole_id = file_identity(_hgcer_hole_path(config["Q2"], config["W"], config["EPSSET"], setting))
        protected.append(hole_id)
        hole = apply_HGCer_hole_cut(config["Q2"], config["W"], config["EPSSET"], setting)
        # Raw SIMC coordinate windows. No data offset is applied to SIMC MM.
        regions = {}
        for target in ("pion_control", "kaon_nosub"):
            for name, windows in bg.resolve_particle_subtraction_component_fit_windows(
                    target, mm_offset_data=0., inp_dict=config, phi_setting=setting).items():
                regions[target + ":" + name] = windows
            for name, windows in bg.resolve_particle_subtraction_component_stage_amplitude_windows(
                    target, mm_offset_data=0., inp_dict=config, phi_setting=setting).items():
                regions[target + ":anchor:" + name] = windows
        regions["under_lambda:configured_analysis_mm_window"] = [(float(config["mm_min"]), float(config["mm_max"]))]
        for component in COMPONENTS:
            require(files.get(component), "missing configured SIMC input: " + component)
            identity = file_identity(files[component])
            protected.append(identity)
            source = ROOT.TFile.Open(identity["path"], "READ")
            require(source and not source.IsZombie(), "cannot open configured ROOT file")
            try:
                tree = source.Get(tree_name)
                require(tree and tree.InheritsFrom("TTree"), "missing/incompatible configured tree")
                branches = {str(b.GetName()) for b in tree.GetListOfBranches()}
                branches.update(str(b.GetName()) for b in tree.GetListOfLeaves())
                validate_branches(branches, phi=children)
                full = Population(mm_axis, t_axis, phi_axis if children else None, regions)
                mmcut = Population(mm_axis, t_axis, phi_axis if children else None, regions)
                # Fresh detached ROOT histogram reproduces TH1D filling exactly.
                mirror = _make_component_hist(component, setting, "detached_support", "mirror", n, lo, hi)
                counters = accumulate_events(tree, config, full, mmcut, hole, apply_simc_cuts,
                                             apply_simc_sub_cuts, find_2d_bin_indices, mirror)
                require(not counters["invalid_weight"] and not counters["selected_invalid_mm"],
                        "loader population cannot be reproduced: " + str(counters))
                actual = load_pion_simc_component_shape(identity["path"], config, setting, "kaon", component,
                                                        hgcer_cutg=hole, context="detached_support_parity")
                parity = production_parity(full.setting, actual["setting_shape_full"], actual["diagnostics"])
                require(counters["seen"] == actual["diagnostics"]["n_events_seen"]
                        and full.setting.report()["total"]["count"] == actual["diagnostics"]["n_events_passed"]
                        and mmcut.setting.report()["total"]["count"] == actual["diagnostics"]["n_events_passed_mm_window"],
                        "production population count mismatch")
                mirror.Scale(1. / mirror.Integral())
                for i in range(n + 2):
                    require(mirror.GetBinContent(i) == actual["setting_shape_full"].GetBinContent(i)
                            and mirror.GetBinError(i) == actual["setting_shape_full"].GetBinError(i),
                            "exact production ROOT histogram parity mismatch")
                entry["components"][component] = {"status": "SOURCE VERIFIED", "input": identity,
                    "tree": tree_name, "tree_entries": int(tree.GetEntries()), "branches": sorted(branches),
                    "iter_weight_default_used": "iter_weight" not in branches, "hole_input": hole_id,
                    "counts": counters, "region_metadata": regions, "full_template_population": full.report(),
                    "additional_mm_cut_population": mmcut.report(), "production_shape_parity": parity,
                    "scientific_adequacy": "NOT VERIFIED"}
            finally:
                source.Close()
    for identity in protected:
        require(file_identity(identity["path"]) == identity, "input changed during scan")
    for identity in report["source_files"]:
        require(file_identity(identity["path"]) == identity, "source changed during scan")
    require(subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip() == source_head,
            "HEAD changed during scan")
    report["measurement_status"] = "SOURCE VERIFIED"
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--config-sha256", required=True)
    parser.add_argument("--settings", nargs="+", choices=("Center", "Left"), required=True)
    parser.add_argument("--children", action="store_true")
    parser.add_argument("--center-evidence", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        require(not args.report.exists(), "report already exists; overwrite forbidden")
        identity = file_identity(args.config)
        require(identity["sha256"] == args.config_sha256, "config SHA256 mismatch")
        config = json.loads(args.config.read_text(encoding="utf-8"))
        evidence_id = file_identity(args.center_evidence) if args.center_evidence else None
        evidence = json.loads(args.center_evidence.read_text(encoding="utf-8")) if evidence_id else None
        report = scan(config, list(dict.fromkeys(args.settings)), args.children, evidence)
        require(file_identity(args.config) == identity, "config changed during scan")
        if evidence_id:
            require(file_identity(args.center_evidence) == evidence_id, "optional evidence changed during scan")
        report.update(config_input=identity, optional_center_evidence=evidence_id)
        with args.report.open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
        return 0
    except (ValueError, OSError, ImportError, KeyError, TypeError, RuntimeError) as exc:
        print("NOT VERIFIED: " + str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
