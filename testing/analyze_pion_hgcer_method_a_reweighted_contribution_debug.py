"""One detached October 4 historical diagnostic invocation; no ROOT or analysis.

Locator JSON has exactly {"f1": {"Left-lowe": "path", ...}, "f3": "path",
"f4": "path"}. Relative paths resolve against the locator's directory. No
digests or authority overrides are accepted. --preflight-only validates the
same complete scientific chain without writing outputs. Paths are never saved.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import uuid
import zipfile

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src" / "cuts"))
import pion_hgcer_method_a_reweighted_contribution_debug as debug

MEMBERS = ("diagnostic.json", "diagnostic.pdf", "manifest.json")
PAGE_TITLES = (
    "Authority and descriptive parent support",
    "Actual baseline and reweighted distributions",
    "Tail concentration and nonnegative effective statistics",
    "Source-class signed estimates and absolute support",
    "Nine-phi and fixed missing-mass redistributions",
    "t0 joint locality", "t1 joint locality", "t2 joint locality",
)


def sha256(raw):
    return hashlib.sha256(raw).hexdigest()


def json_bytes(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+"\n").encode("utf-8")


def parse_json(raw):
    def pairs(items):
        result = {}
        for key, value in items:
            debug.require(key not in result, "duplicate_json_key:" + key)
            result[key] = value
        return result
    def invalid(value):
        raise debug.ContributionDebugError("nonfinite_json:" + value)
    value = json.loads(raw.decode("utf-8"), object_pairs_hook=pairs, parse_constant=invalid)
    debug.require(isinstance(value, dict), "json_object_required")
    return value


def load_inputs(locator):
    locator = Path(locator).resolve(strict=True)
    paths = parse_json(locator.read_bytes())
    debug.require(set(paths) == {"f1", "f3", "f4"} and isinstance(paths["f1"], dict)
                  and set(paths["f1"]) == set(debug.F1_HASHES), "input_manifest_role_inventory_invalid")
    supplied = {**paths["f1"], "f3": paths["f3"], "f4": paths["f4"]}
    resolved, seen = {}, set()
    for role, name in supplied.items():
        debug.require(isinstance(name, str) and bool(name.strip()), "input_manifest_path_invalid:" + role)
        p = Path(name)
        p = (p if p.is_absolute() else locator.parent/p).resolve()
        debug.require(p.is_file(), "historical_input_missing:" + role)
        debug.require(p not in seen, "input_manifest_duplicate_file:" + role)
        seen.add(p); resolved[role] = p
    # Verify ALL raw digests before parsing scientific data or computing anything.
    expected = {**debug.F1_HASHES, "f3": debug.F3_HASH, "f4": debug.F4_HASH}
    raw, hashes = {}, {}
    for role, p in resolved.items():
        raw[role] = p.read_bytes(); hashes[role] = sha256(raw[role])
        debug.require(hashes[role] == expected[role], "historical_raw_sha256_mismatch:" + role)
    payloads = {role: parse_json(data) for role, data in raw.items()}
    # Basenames can still contain private path separators on another platform.
    basenames = {role: p.name for role, p in resolved.items()}
    debug.require(all("/" not in n and "\\" not in n for n in basenames.values()), "input_basename_invalid")
    return ([payloads[s] for s in debug.F1_HASHES], payloads["f3"], payloads["f4"], hashes, basenames)


def git_provenance():
    def value(args):
        return subprocess.check_output(["git", *args], cwd=REPO_ROOT, text=True).rstrip("\r\n")
    return {"git_head": value(["rev-parse", "HEAD"]),
            "git_status_short": value(["status", "--short", "--untracked-files=all"]),
            "runtime_claim": "detached historical diagnostic only; no production promotion"}


def load_diagnostic(locator, reference_path=None):
    f1, f3, f4, hashes, basenames = load_inputs(locator)
    reference, digest = None, None
    if reference_path is not None:
        raw = Path(reference_path).read_bytes(); digest = sha256(raw)
        debug.require(digest == debug.REFERENCE_HASH, "historical_reference_raw_sha256_mismatch")
        reference = parse_json(raw)
    return debug.build_diagnostic(f1, f3, f4, hashes=hashes,
        provenance={**git_provenance(), "input_basenames": basenames},
        reference=reference, reference_sha=digest)


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
    lineage = "October 4 exact historical F.1/F.3/F.4; Q4p4W2p74 / Left / lowe"
    units = "signed/absolute pion-background weight units (not final kaon yield)"
    with PdfPages(path, metadata={"Title": "Historical Method-A contribution debug",
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
                 "Historical F.4 SHA: " + debug.F4_HASH,
                 "Historical correction: " + debug.F4_CORRECTION,
                 "Observed diagnostic Git HEAD: " + value["provenance"]["git_head"],
                 "Original seven input hashes and F.1 fingerprints are in diagnostic.json; paths omitted.",
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
        labels = [f"{a:g}:{b:g}" for a, b in zip(phi_edges, phi_edges[1:])]
        for t, p in enumerate(parents):
            for col, marginal in enumerate((p["phi"], p["mm"])):
                for field, style in (("B", "-"), ("B_rw", "--"), ("signed_delta", ":")):
                    axes[t, col].plot(range(len(marginal)), [s[field] for s in marginal], style, label=field)
                axes[t, col].set_ylabel(f"t{t} signed background units"); axes[t, col].legend(fontsize=7)
                if col == 0:
                    axes[t, col].set_xticks(range(9), labels, rotation=35, fontsize=6)
                    axes[t, col].set_xlabel("canonical phi [degrees]")
                else:
                    axes[t, col].set_xlabel("MM bin: 0=UF, 1..40=[0.8,1.4] GeV, 41=OF")
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
                image = axis.imshow(array, aspect="auto", origin="lower", **options)
                label = {"signed_delta": "Delta B [signed background units]", "A_rw": "sum |b_rw| [support units]", "neff": "Neff(|b_rw|); empty/zero cells unavailable", "count": "application records"}[field]
                axis.set_title(label, fontsize=8); fig.colorbar(image, ax=axis, shrink=.7)
                axis.set_yticks(range(9), labels, fontsize=6)
                axis.set_ylabel("phi [degrees]"); axis.set_xlabel("MM bin 0=UF, 1..40 normal, 41=OF")
            save(fig, 5+t)
    return {"pages": len(PAGE_TITLES), "page_titles": list(PAGE_TITLES), "renderer_input": "persisted aggregates only"}


def check_pdf(raw):
    # Matplotlib writes uncompressed page dictionaries and one page tree.
    # Check the independent object inventory AND declared tree count, not just
    # the renderer's reported save count. No external executable is launched.
    debug.require(raw.startswith(b"%PDF-") and raw.rstrip().endswith(b"%%EOF"), "pdf_container_invalid")
    count = len(re.findall(rb"/Type\s*/Page\b", raw))
    trees = re.findall(rb"/Type\s*/Pages\b[^>]*?/Count\s+(\d+)", raw, re.S)
    debug.require(count == len(PAGE_TITLES) and trees == [str(count).encode()], "pdf_page_completeness_failed")
    return count


def create_zip(path, members):
    """Deterministic ZIP metadata and order for a fixed set of member bytes."""
    with zipfile.ZipFile(path, "x", compression=zipfile.ZIP_DEFLATED) as archive:
        for name in MEMBERS:
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info, members[name])


def verify_zip(path, expected=None):
    with zipfile.ZipFile(path) as archive:
        debug.require(archive.namelist() == list(MEMBERS) and archive.testzip() is None, "zip_member_inventory_or_crc_invalid")
        members = {name: archive.read(name) for name in MEMBERS}
    if expected is not None:
        debug.require(members == expected, "zip_member_bytes_mismatch")
    manifest = parse_json(members["manifest.json"])
    debug.require(manifest["schema_version"] == "pion_hgcer_method_a_reweighted_contribution_debug_bundle/v1" and manifest["status"] == "non_authoritative", "bundle_manifest_schema_invalid")
    # A manifest cannot contain its own raw hash. It covers every other member;
    # its own canonical fingerprint and the ZIP CRC provide its self-check.
    core = dict(manifest); observed = core.pop("fingerprint", None)
    debug.require(debug.fingerprint(core) == observed, "bundle_manifest_fingerprint_invalid")
    debug.require(set(manifest["members"]) == {"diagnostic.json", "diagnostic.pdf"}, "bundle_digest_inventory_invalid")
    for name, record in manifest["members"].items():
        debug.require(record == {"sha256": sha256(members[name]), "bytes": len(members[name])}, "bundle_member_digest_invalid:" + name)
    diagnostic = parse_json(members["diagnostic.json"])
    debug.validate_diagnostic(diagnostic)
    pages = check_pdf(members["diagnostic.pdf"])
    debug.require(manifest["diagnostic_fingerprint"] == diagnostic["fingerprint"] and manifest["pdf"]["pages"] == pages and manifest["pdf"]["page_titles"] == list(PAGE_TITLES) and manifest["pdf"]["renderer_input"] == "persisted aggregates only", "bundle_receipt_invalid")
    return manifest


def roots(outdir, bundle_dir):
    outdir, bundle_dir = Path(outdir).resolve(), Path(bundle_dir).resolve()
    for root in (outdir, bundle_dir):
        debug.require(root.is_dir(), "existing_output_root_required")
        debug.require(not root.is_relative_to(REPO_ROOT), "checkout_output_forbidden")
    return outdir, bundle_dir


def run(args):
    artifact_root, bundle_root = roots(args.outdir, args.bundle_dir)
    attempt = "kaonlt-historical-contribution-" + uuid.uuid4().hex
    final = bundle_root/(attempt+".zip")
    debug.require(not final.exists(), "output_already_exists")
    diagnostic = load_diagnostic(args.inputs_manifest, args.reference_json)
    debug.validate_diagnostic(diagnostic)
    if args.preflight_only:
        return debug.canonical_json({"preflight": "PASS", "input_file_hashes": diagnostic["input_file_hashes"], "historical_f4_exact_reproduction": diagnostic["f4_exact_reproduction"], "outputs_written": False})
    created = False
    # Only this owner's directory is cleaned, even on gate failures.
    with tempfile.TemporaryDirectory(dir=artifact_root, prefix="."+attempt+"-") as temporary:
        stage = Path(temporary)
        js, pdf = stage/"diagnostic.json", stage/"diagnostic.pdf"
        js.write_bytes(json_bytes(diagnostic))
        debug.validate_diagnostic(parse_json(js.read_bytes()))
        receipt = write_review_pdf(pdf, js)
        debug.require(receipt["pages"] == check_pdf(pdf.read_bytes()), "renderer_receipt_mismatch")
        members = {"diagnostic.json": js.read_bytes(), "diagnostic.pdf": pdf.read_bytes()}
        manifest = debug.seal({"schema_version": "pion_hgcer_method_a_reweighted_contribution_debug_bundle/v1",
                               "status": "non_authoritative", "diagnostic_fingerprint": diagnostic["fingerprint"],
                               "members": {name: {"bytes": len(raw), "sha256": sha256(raw)} for name, raw in members.items()},
                               "pdf": receipt, "checks": {"persisted_json": True, "numerical_closure": True, "pdf_completeness": True},
                               "manifest_self_digest": "canonical fingerprint excluding fingerprint; raw self-hash cannot be embedded"})
        members["manifest.json"] = json_bytes(manifest)
        staged_zip = stage/"bundle.zip"
        create_zip(staged_zip, members); verify_zip(staged_zip, members)
        try:
            # Exclusive creation prevents races or overwrites; no accepted ZIP
            # remains when the final independent verification fails.
            with final.open("xb") as destination:
                created = True
                with staged_zip.open("rb") as source:
                    shutil.copyfileobj(source, destination)
            verify_zip(final, members)
            debug.require(final.is_file() and final.parent == bundle_root, "final_bundle_path_invalid")
        except Exception:
            if created:
                final.unlink()
            raise
    return str(final)


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--inputs-manifest", type=Path, required=True)
    result.add_argument("--outdir", type=Path, required=True)
    result.add_argument("--bundle-dir", type=Path, required=True)
    result.add_argument("--kinematic", choices=(debug.KINEMATIC,), required=True)
    result.add_argument("--phi", choices=("Left",), required=True)
    result.add_argument("--epsilon", choices=("lowe",), required=True)
    result.add_argument("--reference-json", type=Path)
    result.add_argument("--preflight-only", action="store_true")
    return result


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        result = run(args)
    except (ValueError, OSError, KeyError, TypeError, OverflowError, RuntimeError, subprocess.SubprocessError, zipfile.BadZipFile) as exc:
        print("BLOCKED: " + str(exc), file=sys.stderr)
        return 1
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
