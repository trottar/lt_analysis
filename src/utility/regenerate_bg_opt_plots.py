#! /usr/bin/python

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _infer_output_dir_for_csv(csv_path):
    parent = csv_path.parent
    if parent.name.lower() == "csv":
        return parent.parent / "plots"
    return None


def _resolve_search_root(path, pattern):
    if path.is_file():
        return [path]

    if not path.is_dir():
        return sorted(Path().glob(str(path)))

    lower_name = path.name.lower()
    if lower_name == "plots" and (path.parent / "csv").is_dir():
        return sorted((path.parent / "csv").rglob(pattern))
    if lower_name == "csv":
        return sorted(path.rglob(pattern))
    if (path / "csv").is_dir():
        return sorted((path / "csv").rglob(pattern))
    return sorted(path.rglob(pattern))


def _collect_csv_paths(inputs, mode, pattern):
    csv_paths = []
    seen = set()

    for raw_input in inputs:
        path = Path(raw_input).expanduser()
        candidates = _resolve_search_root(path, pattern)

        for candidate in candidates:
            if not candidate.is_file() or candidate.suffix.lower() != ".csv":
                continue
            name = candidate.name.lower()
            if mode != "all" and f"_bg_opt_{mode}" not in name:
                continue
            resolved = candidate.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            csv_paths.append(candidate)

    return sorted(csv_paths)


def _build_output_path(csv_path, output_dir=None, suffix=None):
    if output_dir is None:
        output_dir = _infer_output_dir_for_csv(csv_path)

    if output_dir is None:
        pdf_name = csv_path.with_suffix(".pdf").name
        if suffix:
            pdf_name = "{}_{}.pdf".format(csv_path.stem, suffix)
        return csv_path.with_name(pdf_name)

    output_dir = Path(output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    pdf_name = "{}.pdf".format(csv_path.stem if not suffix else "{}_{}".format(csv_path.stem, suffix))
    return output_dir / pdf_name


def _describe_input_behavior(raw_input):
    path = Path(raw_input).expanduser()
    if path.is_file():
        if path.parent.name.lower() == "csv":
            return "{} -> sibling plots/".format(path)
        return "{} -> same directory".format(path)
    if path.is_dir():
        lower_name = path.name.lower()
        if lower_name == "plots" and (path.parent / "csv").is_dir():
            return "{} -> uses sibling csv/ and writes back to plots/".format(path)
        if lower_name == "csv":
            return "{} -> writes to sibling plots/".format(path)
        if (path / "csv").is_dir():
            return "{} -> uses embedded csv/ and writes to embedded plots/".format(path)
        return "{} -> recursive search, writes next to each CSV".format(path)
    return "{} -> glob expansion".format(raw_input)


def _infer_run_dir_for_csv(csv_path):
    csv_path = Path(csv_path).expanduser()
    if csv_path.parent.name.lower() == "csv":
        return csv_path.parent.parent
    return csv_path.parent


def _infer_analysis_stem(csv_path):
    stem = Path(csv_path).stem
    for suffix in ("_bg_opt_low", "_bg_opt_high"):
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    return stem


def _artifact_paths_for_csv(csv_path):
    run_dir = _infer_run_dir_for_csv(csv_path)
    analysis_stem = _infer_analysis_stem(csv_path)
    return {
        "run_dir": run_dir,
        "analysis_stem": analysis_stem,
        "root": run_dir / "root" / "{}.root".format(analysis_stem),
        "json": run_dir / "json" / "{}.json".format(analysis_stem),
    }


def _clone_histogram(obj):
    if obj is None:
        return None
    cloned = obj.Clone()
    if hasattr(cloned, "SetDirectory"):
        try:
            cloned.SetDirectory(0)
        except Exception:
            pass
    return cloned


def _load_histogram(root_file, directory_name, histogram_name):
    current_dir = root_file.GetDirectory(directory_name)
    if not current_dir:
        return None
    hist = current_dir.Get(histogram_name)
    return _clone_histogram(hist)


def _load_archived_context(csv_path):
    artifacts = _artifact_paths_for_csv(csv_path)
    json_path = artifacts["json"]
    root_path = artifacts["root"]
    if not json_path.exists() or not root_path.exists():
        return None, None, artifacts

    try:
        with open(json_path, "r") as handle:
            payload = json.load(handle)
    except Exception as exc:
        print("WARNING: Failed to load archived JSON {}: {}".format(json_path, exc))
        return None, None, artifacts

    inpDict = payload.get("inpDict")
    histlist = payload.get("histlist", [])
    if not isinstance(histlist, list) or inpDict is None:
        return None, None, artifacts

    try:
        import ROOT
    except Exception as exc:
        print("WARNING: ROOT import failed while loading archived context for {}: {}".format(csv_path, exc))
        return None, None, artifacts

    root_file = ROOT.TFile.Open(str(root_path), "READ")
    if not root_file or root_file.IsZombie():
        print("WARNING: Failed to open archived ROOT file {}".format(root_path))
        return None, None, artifacts

    try:
        enriched_histlist = []
        for entry in histlist:
            hist_entry = dict(entry)
            phi_setting = str(hist_entry.get("phi_setting", "")).strip()
            if phi_setting:
                hist_entry["H_MM_DATA"] = _load_histogram(root_file, "{}/data".format(phi_setting), "H_MM_DATA")
                hist_entry["H_MM_pisub_DATA"] = _load_histogram(root_file, "{}/data".format(phi_setting), "H_MM_pisub_DATA")
                hist_entry["H_MM_fit1sub_DATA"] = _load_histogram(root_file, "{}/data".format(phi_setting), "H_MM_fit1sub_DATA")
                hist_entry["H_MM_fit2sub_DATA"] = _load_histogram(root_file, "{}/data".format(phi_setting), "H_MM_fit2sub_DATA")
                hist_entry["H_MM_full_SIMC"] = _load_histogram(root_file, "{}/simc".format(phi_setting), "H_MM_full_SIMC")
            enriched_histlist.append(hist_entry)
    finally:
        root_file.Close()

    return enriched_histlist, inpDict, artifacts


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Regenerate BG optimization diagnostics PDFs from saved *_bg_opt_*.csv files."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help=(
            "Path(s) to an archived run directory, a csv/ directory, a plots/ directory, "
            "a single *_bg_opt_*.csv file, or a glob."
        ),
    )
    parser.add_argument(
        "--mode",
        choices=("low", "high", "all"),
        default="all",
        help="When scanning directories/globs, restrict to low/high epsilon CSVs. Default: all.",
    )
    parser.add_argument(
        "--pattern",
        default="*_bg_opt_*.csv",
        help="Filename pattern to use when an input is a directory. Default: *_bg_opt_*.csv",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Optional directory to place regenerated PDFs. Default: infer plots/ from archive structure or use the CSV directory.",
    )
    parser.add_argument(
        "--suffix",
        default=None,
        help="Optional suffix to append to regenerated PDF filenames before .pdf.",
    )

    args = parser.parse_args(argv)

    repo_root = Path(__file__).resolve().parents[2]
    plotting_dir = repo_root / "src" / "plotting"
    if str(plotting_dir) not in sys.path:
        sys.path.append(str(plotting_dir))

    from bg_opt_diagnostics import plot_bg_optimization_diagnostics

    csv_paths = _collect_csv_paths(args.inputs, args.mode, args.pattern)
    if not csv_paths:
        print("No matching BG optimization CSV files found.")
        return 1

    print("Found {} BG optimization CSV file(s).".format(len(csv_paths)))
    for raw_input in args.inputs:
        print("  {}".format(_describe_input_behavior(raw_input)))
    print("The utility will use archived json/root files when present so MM pages can be rebuilt from saved histograms.")

    failures = []
    for csv_path in csv_paths:
        pdf_path = _build_output_path(csv_path, output_dir=args.output_dir, suffix=args.suffix)
        try:
            histlist, inpDict, artifacts = _load_archived_context(csv_path)
            if histlist is not None and inpDict is not None:
                print("Using archived context: {} and {}".format(artifacts["json"], artifacts["root"]))
            else:
                print("Archived context unavailable for {}; MM pages may be omitted.".format(csv_path))
            created = plot_bg_optimization_diagnostics(csv_path, pdf_path=pdf_path, histlist=histlist, inpDict=inpDict)
            print("Wrote {}".format(created))
        except Exception as exc:
            failures.append((csv_path, exc))
            print("FAILED {} -> {}".format(csv_path, exc))

    if failures:
        print("\n{} file(s) failed:".format(len(failures)))
        for csv_path, exc in failures:
            print("  {}: {}".format(csv_path, exc))
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
