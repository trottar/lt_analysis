"""One read-only current-snapshot candidate diagnostic owner.

Locator: {"f1": {"Left-lowe": "path", ...}, "f3": "path", "f4": "path"}.
Relative paths resolve against the locator. Names and timestamps confer no
authority. --preflight-only writes nothing, including on failure. Full runs
require a passed preflight's exact --expected-snapshot-fingerprint and existing
external --outdir and --bundle-dir. No input generation or directory search.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import uuid
import zipfile

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src" / "cuts"))
import pion_hgcer_method_a_current_snapshot_contribution_debug as debug
import analyze_pion_hgcer_method_a_reweighted_contribution_debug as historical_io

MEMBERS = historical_io.MEMBERS
PAGE_TITLES = historical_io.PAGE_TITLES
BUNDLE_SCHEMA = "pion_hgcer_method_a_current_snapshot_contribution_debug_bundle/v1"
PREFLIGHT_SCHEMA = "pion_hgcer_method_a_current_snapshot_preflight/v1"
parse_json = historical_io.parse_json
json_bytes = historical_io.json_bytes
sha256 = historical_io.sha256
check_pdf = historical_io.check_pdf
create_zip = historical_io.create_zip
roots = historical_io.roots


def load_inputs(locator, observations):
    """Read exactly seven explicit regular files; collect observed digests first."""
    locator = Path(locator)
    debug.require(locator.is_file(), "input_manifest_missing")
    paths = parse_json(locator.read_bytes())
    debug.require(set(paths) == {"f1", "f3", "f4"} and isinstance(paths["f1"], dict) and
                  set(paths["f1"]) == set(debug.SETTINGS), "input_manifest_role_inventory_invalid")
    supplied = {**paths["f1"], "f3": paths["f3"], "f4": paths["f4"]}
    resolved = {}
    for role in debug.ROLES:
        name = supplied[role]
        debug.require(type(name) is str and bool(name.strip()), "input_manifest_path_invalid:" + role)
        p = Path(name)
        p = (p if p.is_absolute() else locator.resolve().parent/p).resolve()
        debug.require(p.is_file(), "snapshot_input_missing:" + role)
        debug.require(all(p != old and not p.samefile(old) for old in resolved.values()),
                      "input_manifest_duplicate_file:" + role)
        debug.require(not any(c in p.name for c in "/\\:"), "input_basename_invalid:" + role)
        resolved[role] = p
    raw = {role: p.read_bytes() for role, p in resolved.items()}
    hashes = {role: sha256(data) for role, data in raw.items()}
    for role in debug.ROLES:
        observations[role] = {"basename": resolved[role].name, "sha256": hashes[role],
                              "observed_fingerprint": None, "fingerprint_validated": False}
    payloads = {}
    for role in debug.ROLES:
        payload = parse_json(raw[role])
        payloads[role] = payload
        body = payload.get("contract" if role in debug.SETTINGS else
                           "acceptance_map" if role == "f3" else "correction", {})
        fp = body.get("fingerprint") if isinstance(body, dict) else None
        observations[role]["observed_fingerprint"] = fp if (isinstance(fp, str) and len(fp) == 64 and
                                                          all(c in "0123456789abcdef" for c in fp)) else None
        if role in debug.SETTINGS:
            setting = payload.get("setting", {})
            debug.require(setting.get("phi_setting", "") + "-" + setting.get("epsilon_filename_token", "") == role,
                          "snapshot_f1_role_assignment_invalid:" + role)
    return payloads, hashes, {role: p.name for role, p in resolved.items()}, resolved


def git_provenance(basenames):
    def value(args):
        return subprocess.check_output(["git", *args], cwd=REPO_ROOT, text=True).rstrip("\r\n")
    return {"git_head": value(["rev-parse", "HEAD"]), "git_branch": value(["branch", "--show-current"]),
            "worktree_dirty": bool(value(["status", "--porcelain", "--untracked-files=all"])),
            "input_basenames": basenames}


def preflight(locator):
    observations = {role: {"basename": None, "sha256": None, "observed_fingerprint": None,
                           "fingerprint_validated": False} for role in debug.ROLES}
    receipt = {"schema_version": PREFLIGHT_SCHEMA, "non_authoritative": True,
               "preflight": "FAIL", "internal_snapshot_consistency": "FAIL",
               "canonical_owner_equivalence": "BLOCKED", "canonical_owner_equivalence_evaluated": False,
               "observations": observations, "snapshot_fingerprint": None,
               "first_failure": None, "outputs_written": False}
    try:
        payloads, hashes, basenames, paths = load_inputs(locator, observations)
        snapshot, parents = debug.reproduce_snapshot([payloads[s] for s in debug.SETTINGS],
            payloads["f3"], payloads["f4"], hashes=hashes)
        # Retain no paths in the receipt. Raw file replacement during preflight also fails.
        for role, path in paths.items():
            debug.require(sha256(path.read_bytes()) == hashes[role], "snapshot_changed_during_preflight:" + role)
            observations[role]["fingerprint_validated"] = True
            observations[role]["scientific_fingerprints"] = snapshot["scientific_fingerprints"][role]
        receipt.update(preflight="PASS", internal_snapshot_consistency="PASS",
                       snapshot_fingerprint=snapshot["fingerprint"])
        return receipt, (snapshot, parents, basenames, paths)
    except (OSError, subprocess.SubprocessError) as exc:
        # OS errors can carry private absolute paths; report only their class.
        receipt["first_failure"] = "input_read_failed:" + type(exc).__name__
    except (ValueError, KeyError, TypeError, IndexError, OverflowError) as exc:
        receipt["first_failure"] = str(exc)
    return receipt, None


def verify_zip(path, expected=None):
    with zipfile.ZipFile(path) as archive:
        debug.require(archive.namelist() == list(MEMBERS) and archive.testzip() is None,
                      "zip_member_inventory_or_crc_invalid")
        members = {name: archive.read(name) for name in MEMBERS}
    if expected is not None:
        debug.require(members == expected, "zip_member_bytes_mismatch")
    manifest = parse_json(members["manifest.json"])
    debug.require(manifest["schema_version"] == BUNDLE_SCHEMA and manifest["status"] == "non_authoritative",
                  "bundle_manifest_schema_invalid")
    debug.require(manifest["fingerprint"] == debug.fingerprint({k: v for k, v in manifest.items()
                                                               if k != "fingerprint"}), "bundle_manifest_fingerprint_invalid")
    debug.require(set(manifest["members"]) == {"diagnostic.json", "diagnostic.pdf"}, "bundle_digest_inventory_invalid")
    for name, record in manifest["members"].items():
        debug.require(record == {"sha256": sha256(members[name]), "bytes": len(members[name])},
                      "bundle_member_digest_invalid:" + name)
    diagnostic = parse_json(members["diagnostic.json"])
    debug.validate_diagnostic(diagnostic)
    pages = check_pdf(members["diagnostic.pdf"])
    debug.require(manifest["diagnostic_fingerprint"] == diagnostic["fingerprint"] and
                  manifest["snapshot_fingerprint"] == diagnostic["snapshot"]["fingerprint"] and
                  manifest["pdf"] == {"pages": pages, "page_titles": list(PAGE_TITLES),
                                      "renderer_input": "persisted aggregates only"}, "bundle_receipt_invalid")
    debug.require(manifest["checks"] == {"persisted_json": True, "numerical_closure": True,
                  "pdf_completeness": True}, "bundle_checks_invalid")
    return manifest


def run(args):
    debug.require(args.expected_snapshot_fingerprint is not None, "passed_preflight_fingerprint_required")
    debug.digest(args.expected_snapshot_fingerprint)
    debug.require(args.outdir is not None and args.bundle_dir is not None, "external_output_roots_required")
    artifact_root, bundle_root = roots(args.outdir, args.bundle_dir)
    attempt = "kaonlt-current-snapshot-contribution-" + uuid.uuid4().hex
    final = bundle_root/(attempt+".zip")
    debug.require(not final.exists(), "output_already_exists")
    receipt, data = preflight(args.inputs_manifest)
    debug.require(data is not None, receipt["first_failure"])
    snapshot, parents, basenames, paths = data
    debug.require(snapshot["fingerprint"] == args.expected_snapshot_fingerprint,
                  "expected_snapshot_fingerprint_mismatch")
    diagnostic = debug.build_diagnostic(snapshot, parents, provenance=git_provenance(basenames))
    # Recheck the immutable bytes/bindings used above immediately before writing.
    for role, path in paths.items():
        debug.require(sha256(path.read_bytes()) == snapshot["input_file_hashes"][role],
                      "snapshot_changed_before_production:" + role)
    debug.validate_diagnostic(diagnostic)
    created = False
    with tempfile.TemporaryDirectory(dir=artifact_root, prefix="."+attempt+"-") as temporary:
        stage = Path(temporary)
        js, pdf = stage/"diagnostic.json", stage/"diagnostic.pdf"
        js.write_bytes(json_bytes(diagnostic))
        debug.validate_diagnostic(parse_json(js.read_bytes()))
        pdf_receipt = write_review_pdf(pdf, js)
        debug.require(pdf_receipt == {"pages": check_pdf(pdf.read_bytes()), "page_titles": list(PAGE_TITLES),
                                     "renderer_input": "persisted aggregates only"}, "renderer_receipt_mismatch")
        members = {"diagnostic.json": js.read_bytes(), "diagnostic.pdf": pdf.read_bytes()}
        manifest = debug.seal({"schema_version": BUNDLE_SCHEMA, "status": "non_authoritative",
            "diagnostic_fingerprint": diagnostic["fingerprint"], "snapshot_fingerprint": snapshot["fingerprint"],
            "members": {name: {"bytes": len(raw), "sha256": sha256(raw)} for name, raw in members.items()},
            "pdf": pdf_receipt, "checks": {"persisted_json": True, "numerical_closure": True, "pdf_completeness": True},
            "manifest_self_digest": "canonical fingerprint excluding fingerprint; raw self-hash cannot be embedded"})
        members["manifest.json"] = json_bytes(manifest)
        staged_zip = stage/"bundle.zip"
        create_zip(staged_zip, members)
        verify_zip(staged_zip, members)
        try:
            with final.open("xb") as destination:
                created = True
                with staged_zip.open("rb") as source:
                    shutil.copyfileobj(source, destination)
            verify_zip(final, members)
        except Exception:
            if created:
                final.unlink()
            raise
    return str(final)


def mm_display_geometry(mm_edges):
    """Coordinates from validated persisted edges; flow strips are categorical."""
    centers = [(a+b)/2 for a, b in zip(mm_edges, mm_edges[1:])]
    width = mm_edges[1]-mm_edges[0]
    return centers, (mm_edges[0]-width, mm_edges[-1]+width), (
        mm_edges[0]-width/2, mm_edges[-1]+width/2)


def label_mm_axis(axis, mm_edges):
    """Mark physical boundaries separately from the UF/OF display strips."""
    _, limits, flows = mm_display_geometry(mm_edges)
    axis.set_xlim(*limits)
    ticks = [mm_edges[0] + i*(mm_edges[-1]-mm_edges[0])/6 for i in range(7)]
    axis.set_xticks(ticks, [f"{x:.1f}" for x in ticks], fontsize=7)
    for boundary in (mm_edges[0], mm_edges[-1]):
        axis.axvline(boundary, color="black", linewidth=.7, linestyle="--")
    for x, label in zip(flows, ("UF", "OF")):
        axis.text(x, 1.02, label, transform=axis.get_xaxis_transform(),
                  ha="center", va="bottom", fontsize=7, rotation=90)
    axis.set_xlabel("Analysis missing mass [GeV]\nUF/OF: outside-range categories", fontsize=8)


def write_review_pdf(path, persisted_json):
    """Only read serialized aggregates; never accepts factors or F.1 objects."""
    value = parse_json(Path(persisted_json).read_bytes())
    debug.validate_diagnostic(value)
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import numpy as np

    parents = value["parents"]
    lineage = "current-snapshot candidate — detached diagnostic; Q4p4W2p74 / Left / lowe"
    units = "signed/absolute pion-background weight units (not final kaon yield)"
    with PdfPages(path, metadata={"Title": "Current-snapshot candidate contribution debug",
                                "Author": "KaonLT", "CreationDate": None, "ModDate": None}) as pdf:
        def save(fig, number):
            fig.suptitle(f"{number+1}. {PAGE_TITLES[number]}", fontsize=13)
            fig.text(.04, .025, lineage + "\nDetached, non-authoritative. Neff is descriptive; no yield uncertainty or absolute mis-ID calibration.", fontsize=7)
            fig.tight_layout(rect=(.02, .08, .98, .94))
            pdf.savefig(fig); plt.close(fig)

        fig = plt.figure(figsize=(11.7, 8.3))
        lines = ["Primary focus t0/t2; t1 internal reference. Exact shared F.4 reproduction; no child normalization.",
                 "B=sum(s*w0); B_rw=sum(s*w0*C); A=sum|b|; A_rw=sum|b_rw|; zeros retained.",
                 "Signed estimates and nonnegative support are distinct. Neither is an extracted kaon yield.",
                 "Neff=(sum nonnegative measure)^2/sum measure^2; not errors/covariance/Poisson counts.",
                 "Observed F.4 SHA: " + value["snapshot"]["input_file_hashes"]["f4"],
                 "Internally reproduced correction: " + value["snapshot"]["scientific_fingerprints"]["f4"]["body"],
                 "Observed diagnostic Git HEAD: " + value["provenance"]["git_head"],
                 "Snapshot seven input hashes and F.1 fingerprints are in diagnostic.json; paths omitted.",
                 "parent  count       B              B_rw           A              A_rw          |B|/A"]
        for p in parents:
            s = p["totals"]
            lines.append(f"t{p['canonical_t_index']}   {s['count']:6d}  " + "  ".join(f"{s[k]:.8g}" for k in ("B", "B_rw", "A", "A_rw")) + f"  {s['cancellation']['value']:.5g}")
            lines.append(f"     signed delta={s['signed_delta']:.8g}; sum b^2={s['sum_b2']:.8g}; sum b_rw^2={s['sum_b_rw2']:.8g}")
        lines.extend(["MM display: 40 fixed bins [0.8,1.4] GeV plus underflow/overflow; nine frozen phi children.",
                      "No clipping, weight cap, adaptive physics binning, alternate factor, production mutation or Method B.",
                      "Private canonical-five F.4 mismatch and absolute-SIMC units blockers remain unchanged."])
        fig.text(.04, .90, "\n\n".join(lines), va="top", family="monospace", fontsize=8)
        save(fig, 0)

        fig, axes = plt.subplots(3, 2, figsize=(11.7, 8.3))
        for t, p in enumerate(parents):
            for col, name in enumerate(("weights", "contributions")):
                h = p["histograms"][name]; edges = h["edges"]
                for field, label in (("baseline_counts", "baseline"), ("reweighted_counts", "reweighted")):
                    axes[t, col].stairs(h[field][1:-1], edges, label=label)
                axes[t, col].set_yscale("symlog", linthresh=1)
                axes[t, col].set_ylabel(f"t{t} records")
                axes[t, col].set_xlabel("w0 / w0*C" if col == 0 else "|b| / |b_rw| [background weight units]")
                axes[t, col].set_title("Common full-range bins; no clipping", fontsize=8)
                axes[t, col].legend(fontsize=7)
        save(fig, 1)

        fig, axes = plt.subplots(2, 2, figsize=(11.7, 8.3))
        for rank, style in (("C", "o-"), ("abs_b_rw", "s--"), ("abs_b", "^:")):
            for j, q in enumerate((.01, .05)):
                axes[0, j].plot(range(3), [p["tails"][rank][j]["reweighted_abs_share"]["value"] for p in parents], style, label="rank by " + rank)
                axes[0, j].set_title(f"Largest {q:.0%} records: share of total |b_rw|")
                axes[0, j].set_ylabel("nonnegative support share")
                axes[0, j].legend(fontsize=7)
        for axis, names in zip(axes[1], (("w0", "w_rw"), ("abs_b", "abs_b_rw"))):
            for name in names:
                axis.plot(range(3), [p["totals"]["neff"][name]["value"] for p in parents], "o-", label=name)
            axis.set_ylabel("descriptive nonnegative Neff"); axis.legend(fontsize=7)
        for axis in axes.flat:
            axis.set_xticks(range(3), ("t0", "t1", "t2")); axis.grid(alpha=.2)
        save(fig, 2)

        fig, axes = plt.subplots(3, 2, figsize=(11.7, 8.3))
        for t, p in enumerate(parents):
            for col, fields in enumerate((("B", "B_rw"), ("A", "A_rw"))):
                for j, field in enumerate(fields):
                    axes[t, col].bar(np.arange(4)+(j-.5)*.35, [p["sources"][s]["totals"][field] for s in debug.SOURCES], width=.35, label=field)
                axes[t, col].set_xticks(range(4), debug.SOURCES, fontsize=7)
                axes[t, col].set_ylabel(f"t{t} " + ("signed estimate" if col == 0 else "absolute support"))
                axes[t, col].axhline(0, color="black", linewidth=.5); axes[t, col].legend(fontsize=7)
        fig.text(.04, .065, units, fontsize=7)
        save(fig, 3)

        fig, axes = plt.subplots(3, 2, figsize=(11.7, 8.3))
        phi_edges = value["display_geometry"]["phi_edges"]
        mm_edges = value["display_geometry"]["mm_edges"]
        mm_centers, mm_limits, mm_flows = mm_display_geometry(mm_edges)
        labels = [f"{a:g}:{b:g}" for a, b in zip(phi_edges, phi_edges[1:])]
        for t, p in enumerate(parents):
            for col, marginal in enumerate((p["phi"], p["mm"])):
                for field, style in (("B", "-"), ("B_rw", "--"), ("signed_delta", ":")):
                    if col == 0:
                        axes[t, col].plot(range(len(marginal)), [s[field] for s in marginal], style, label=field)
                    else:
                        line, = axes[t, col].plot(mm_centers, [s[field] for s in marginal[1:-1]], style, label=field)
                        axes[t, col].plot(mm_flows, [marginal[0][field], marginal[-1][field]],
                                          linestyle="none", marker="o", markersize=3, color=line.get_color())
                axes[t, col].set_ylabel(f"t{t} signed background units"); axes[t, col].legend(fontsize=7)
                if col == 0:
                    axes[t, col].set_xticks(range(9), labels, rotation=35, fontsize=6)
                    axes[t, col].set_xlabel("canonical phi [degrees]")
                else:
                    label_mm_axis(axes[t, col], mm_edges)
                    for low, high in ((mm_limits[0], mm_edges[0]), (mm_edges[-1], mm_limits[1])):
                        axes[t, col].axvspan(low, high, color="grey", alpha=.12)
        save(fig, 4)

        for t, p in enumerate(parents):
            fig, axes = plt.subplots(2, 2, figsize=(11.7, 8.3))
            for axis, field in zip(axes.flat, ("signed_delta", "A_rw", "neff", "count")):
                matrix = [[s["neff"]["abs_b_rw"]["value"] if s["neff"]["abs_b_rw"]["available"] else np.nan for s in row] for row in p["joint"]] if field == "neff" else [[s[field] for s in row] for row in p["joint"]]
                array = np.asarray(matrix)
                options = {}
                if field == "signed_delta":
                    limit = max(float(np.max(np.abs(array))), 1e-30)
                    options = {"vmin": -limit, "vmax": limit, "cmap": "RdBu_r"}
                image = axis.imshow(array, aspect="auto", origin="lower", interpolation="nearest",
                                    extent=(*mm_limits, -.5, 8.5), **options)
                label = {"signed_delta": "Delta B [signed background units]", "A_rw": "sum |b_rw| [support units]", "neff": "Neff(|b_rw|); empty/zero cells unavailable", "count": "application records"}[field]
                axis.set_title(label, fontsize=8); fig.colorbar(image, ax=axis, shrink=.7)
                axis.set_yticks(range(9), labels, fontsize=6)
                axis.set_ylabel("phi [degrees]"); label_mm_axis(axis, mm_edges)
            save(fig, 5+t)
    return {"pages": len(PAGE_TITLES), "page_titles": list(PAGE_TITLES), "renderer_input": "persisted aggregates only"}


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--inputs-manifest", type=Path, required=True)
    result.add_argument("--preflight-only", action="store_true")
    result.add_argument("--expected-snapshot-fingerprint")
    result.add_argument("--outdir", type=Path)
    result.add_argument("--bundle-dir", type=Path)
    return result


def main(argv=None):
    args = parser().parse_args(argv)
    if args.preflight_only:
        receipt, _ = preflight(args.inputs_manifest)
        print(debug.canonical_json(receipt))
        return 0 if receipt["preflight"] == "PASS" else 1
    try:
        result = run(args)
    except (ValueError, OSError, KeyError, TypeError, IndexError, OverflowError,
            RuntimeError, subprocess.SubprocessError, zipfile.BadZipFile) as exc:
        print("BLOCKED: " + str(exc), file=sys.stderr)
        return 1
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
