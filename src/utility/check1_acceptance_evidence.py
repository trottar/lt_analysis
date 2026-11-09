"""Default-off, ROOT-independent Check 1 scalar evidence.

This module observes existing selections and fitted transfers. It never fits,
evaluates PID/proton models, traverses a tree, or changes a producer object.
All histogram arrays include ROOT underflow/overflow and contain bin integrals.
"""
from __future__ import annotations

from bisect import bisect_right
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import tempfile
import zipfile
from pathlib import PurePosixPath

import numpy as np

ENABLE_ENV = "LT_ANALYSIS_CHECK1_EVIDENCE"
AUTHORITY_ENV = "LT_ANALYSIS_CHECK1_AUTHORITY"
AUTHORITY_HASH_ENV = "LT_ANALYSIS_CHECK1_AUTHORITY_SHA256"
SCHEMA = "kaonlt_check1_evidence/v1"
REFERENCE_SHA256 = "02ff71998fd51d5e56f07484bb1c71185f78ed51d6ce409b0d34154bda7a6d7a"
SETTINGS = (("Left", "lowe"), ("Center", "lowe"),
            ("Left", "highe"), ("Center", "highe"), ("Right", "highe"))
SOURCES = ("prompt", "rand", "dummy", "dummy_rand")
CUTS = ("allcuts", "nommcuts")
EXECUTABLES = ("src/utility/check1_acceptance_evidence.py",
               "src/utility/bg_optimization.py",
               "src/cuts/pion_t_phi_template_comparison.py", "src/binning/calculate_yield.py",
               "testing/run_pion_t_template_comparison_gate.py",
               "testing/pion_check1_acceptance_evidence_profile_v1.json",
               "src/cuts/pion_background_chain_audit.py", "src/utility/prompt_trees.py",
               "src/main.py", "src/simc_ana/iter_weight.py", "src/simc_ana/compare_simc.py")
REFERENCE_SUMMARIES = {
    "lowe":"KaonLT_pion_t_phi_Q4p4W2p74_lowe_20261028T044038Z-run-summary.json",
    "highe":"KaonLT_HighEpsilon_Acquisition_20261030T022205-run-summary.json"}


def require(ok, reason):
    if not ok:
        raise ValueError("check1:" + reason)


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def digest(value):
    return hashlib.sha256(canonical_json(value)).hexdigest()


def file_hash(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def strict_json(raw):
    def pairs(items):
        result={}
        for k,v in items:
            require(k not in result,"duplicate_json_key:"+k);result[k]=v
        return result
    def invalid(value):
        raise ValueError("check1:nonfinite_json:"+value)
    return json.loads(raw,object_pairs_hook=pairs,parse_constant=invalid)


def safe_archive(archive):
    names=archive.namelist()
    require(len(names)==len(set(names)) and archive.testzip() is None,"archive_inventory_or_crc")
    for item in archive.infolist():
        path=PurePosixPath(item.filename)
        require(not path.is_absolute() and ".." not in path.parts and "\\" not in item.filename
                and not (item.external_attr >> 16) & 0o170000 == 0o120000,"unsafe_archive_member")
        require(item.file_size <= 160*1024*1024,"reference_member_size")


def load_reference(path):
    """Read the exact pinned seven-member packet, without extraction or writes."""
    path=Path(path)
    require(path.is_file() and file_hash(path)==REFERENCE_SHA256,"reference_digest")
    result=dict(reference_sha256=REFERENCE_SHA256,children={},parents={},summaries={},simc={})
    with zipfile.ZipFile(path) as archive:
        safe_archive(archive)
        expected={"reference_manifest.json","current_lowe.zip","current_highe.zip",
                  "SIMC_REFERENCE_lowe.npz","SIMC_REFERENCE_highe.npz",*REFERENCE_SUMMARIES.values()}
        require(set(archive.namelist())==expected,"reference_members")
        manifest=strict_json(archive.read("reference_manifest.json"))
        require(manifest["schema_version"]=="kaonlt_check1_reference_inputs/v1"
                and manifest["experimental_cohorts"]==["20261028T044038Z","20261030T022205"],"reference_schema_or_cohort")
        require(len(manifest["files"])==6 and {r["name"] for r in manifest["files"]}==expected-{"reference_manifest.json"},"reference_manifest_inventory")
        members={n:archive.read(n) for n in expected-{"reference_manifest.json"}}
        for r in manifest["files"]:
            raw=members[r["name"]]
            require(type(r["bytes"]) is int and len(raw)==r["bytes"]
                    and hashlib.sha256(raw).hexdigest()==r["sha256"],"reference_member_hash:"+r["name"])
        for eps in ("lowe","highe"):
            summary=strict_json(members[REFERENCE_SUMMARIES[eps]])
            require(summary["schema_version"]=="pion_t_template_comparison_run_summary/v1"
                    and summary["analysis_returncode"]==0,"reference_run_summary")
            require(summary["zip_identity"]["sha256"]==hashlib.sha256(members["current_"+eps+".zip"]).hexdigest()
                    and summary["zip_identity"]["bytes"]==len(members["current_"+eps+".zip"]),"reference_summary_zip_identity")
            result["summaries"][eps]=summary
            with zipfile.ZipFile(io.BytesIO(members["current_"+eps+".zip"])) as child:
                safe_archive(child)
                manifest_child=strict_json(child.read("manifest.json"))
                require(manifest_child["complete"] is True and manifest_child["errors"]==[],"reference_nested_manifest")
                declarations=manifest_child["global_artifacts"]
                basenames={"{}_kaon_Q4p4W2p74_{}_{}".format(s,e,suffix) for s,e in SETTINGS if e==eps
                           for suffix in ("pion-t-template-comparison.json","pion-t-template-comparison.pdf",
                                          "pion-t-phi-template-comparison.json")}
                require(len(declarations)==len(basenames) and set(child.namelist())==
                        {"global/","manifest.json","source_state.txt","source_checks.txt",
                         *(s+"_"+e+"/" for s,e in SETTINGS),*("global/"+b for b in basenames)},"reference_nested_inventory")
                for record in declarations.values():
                    name=record["archive_path"] if "archive_path" in record else "global/"+record["basename"]
                    raw=child.read(name)
                    require(hashlib.sha256(raw).hexdigest()==record["sha256"]
                            and len(raw)==record["byte_size"],"reference_nested_member_hash")
                for setting,e in SETTINGS:
                    if e!=eps:continue
                    stem="global/{}_kaon_Q4p4W2p74_{}_pion-t".format(setting,eps)
                    parent=strict_json(child.read(stem+"-template-comparison.json"))
                    doc=strict_json(child.read(stem+"-phi-template-comparison.json"))
                    require(parent["schema_version"]=="pion_t_template_comparison/v1"
                            and doc["schema_version"]=="pion_t_phi_template_comparison/v1"
                            and doc["phi_setting"]==setting and len(doc["cells"])==27
                            and len(parent["parents"])==3,"reference_setting_schema")
                    validate_coordinate_contract(doc)
                    result["children"][setting+":"+eps]=doc
                    result["parents"][setting+":"+eps]=parent
            with np.load(io.BytesIO(members["SIMC_REFERENCE_"+eps+".npz"]),allow_pickle=False) as loaded:
                require(len(loaded.files)==len(set(loaded.files)),"reference_npz_duplicates")
                result["simc"][eps]={k:loaded[k] for k in loaded.files}
            simc=result["simc"][eps]
            axis_keys={"t_bins","phi_bins","mm_edges","q2_edges","w_edges","theta_cm_edges",
                       "theta_cm_true_edges","q2_w_x_edges","q2_w_y_edges","t_vs_tmin_x_edges","t_vs_tmin_y_edges"}
            expected_keys=axis_keys|{"simc_"+var+"_"+s.lower()+"_"+field
                for s,e in SETTINGS if e==eps for var in ("mm","q2","w","theta_cm","theta_cm_true")
                for field in ("values","errors")}|{"simc_"+var+"_"+s.lower()+"_values"
                for s,e in SETTINGS if e==eps for var in ("q2_w","t_vs_tmin")}
            require(set(simc)==expected_keys,"reference_simulation_only_inventory")
            for k,v in simc.items():
                require(v.dtype==np.dtype("float64") and np.isfinite(v).all(),"reference_simc_dtype:"+k)
                if k in axis_keys:axis(v,k)
                elif k.endswith("_errors"):require(np.all(v>=0),"reference_negative_error")
            for setting,e in SETTINGS:
                if e!=eps:continue
                doc=result["children"][setting+":"+eps]
                require(simc["t_bins"].tolist()==doc["t_edges"] and simc["phi_bins"].tolist()==doc["phi_edges"]
                        and simc["simc_mm_"+setting.lower()+"_values"].shape==(3,9,100)
                        and simc["simc_mm_"+setting.lower()+"_errors"].shape==(3,9,100),"reference_native_axes_or_shape")
                for particle in ("kaon","pion"):
                    for role in ("Data","Dummy"):
                        require("{}_{}_Analysed_{}_Q4p4W2p74_{}.root".format(setting,particle,role,eps)
                                in summary["original_inputs"],"reference_input_role")
    return result


# These are execution/container lineage fields, not scientific quantities.
# The complete remaining child document is compared recursively, with no
# data-dependent choice of fields after a failed comparison.
# The existing coordinate fingerprint hashes plot_filename as well as the
# physical transform. Validate its complete self-identity first, retain both
# original/runtime fingerprints, and compare every transform/provenance value
# except the explicitly named execution paths below.
LINEAGE_KEYS = frozenset(("source_HEAD","source_files","source_identity_complete","runtime_inputs",
    "pion_parent_id","baseline_parent_fingerprint","canonical_interval_pair_id","canonical_interval_pair_hash",
    "path","plot_filename","resolved"))
COORDINATE_IDENTITY_KEYS = frozenset(("coordinate_fingerprint","kaon_data_coordinate_fingerprint",
                                     "pion_control_coordinate_fingerprint"))


def scientific_projection(value, coordinate=None):
    if isinstance(value,dict):
        candidate=(value if value.get("schema_version")=="kaon_data_coordinate/v1" else
                   value.get("coordinate") or value.get("kaon_data_coordinate"))
        if isinstance(candidate,dict) and candidate.get("schema_version")=="kaon_data_coordinate/v1":
            coordinate=candidate
    if isinstance(value,dict):
        result={}
        for k,v in value.items():
            if k in LINEAGE_KEYS:continue
            if k in COORDINATE_IDENTITY_KEYS and isinstance(v,str):
                require(coordinate is not None and v==coordinate["coordinate_fingerprint"],"coordinate_role_mismatch")
                # These exact aliases are assigned from the full coordinate
                # contract in pion_t_bin_parents. Keep the fields, validate
                # their original self-identity, and compare the entire physical
                # transform/provenance under its path-independent identity.
                physical={key:scientific_projection(item) for key,item in coordinate.items()
                          if key not in COORDINATE_IDENTITY_KEYS|{"fingerprint","plot_filename"}}
                result[k]=digest(physical)
            elif k=="root_filename" and isinstance(v,str):
                result[k]=PurePosixPath(v.replace("\\","/")).name
            elif k=="artifact_paths" and isinstance(v,list):
                result[k]=[PurePosixPath(p.replace("\\","/")).name if isinstance(p,str)
                           else scientific_projection(p,coordinate) for p in v]
            else:result[k]=scientific_projection(v,coordinate)
        return result
    if isinstance(value,list):return [scientific_projection(v,coordinate) for v in value]
    return value


def validate_coordinate_contract(child):
    coordinate=child["coordinate"]
    require(coordinate["schema_version"]=="kaon_data_coordinate/v1"
            and coordinate["source_particle"]=="kaon"
            and coordinate["coordinate_system"]=="kaon_derived_experimental_analysis"
            and coordinate["raw_mm_field"]=="MM" and coordinate["raw_t_expression"]=="-MandelT", "coordinate_contract")
    payload={k:v for k,v in coordinate.items() if k not in ("coordinate_fingerprint","fingerprint")}
    fingerprint=hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(",",":"),ensure_ascii=True,allow_nan=False).encode()).hexdigest()
    require(coordinate["coordinate_fingerprint"]==fingerprint,"coordinate_self_identity")
    def visit(value):
        if isinstance(value,dict):
            for k,v in value.items():
                if k in COORDINATE_IDENTITY_KEYS and isinstance(v,str):require(v==fingerprint,"coordinate_role_mismatch")
                visit(v)
        elif isinstance(value,list):
            for v in value:visit(v)
    visit(child)


def compare_science(actual,expected,path="reference"):
    if isinstance(expected,dict):
        require(isinstance(actual,dict) and set(actual)==set(expected),"scientific_keys:"+path)
        for k in expected:compare_science(actual[k],expected[k],path+"."+k)
    elif isinstance(expected,list):
        require(isinstance(actual,list) and len(actual)==len(expected),"scientific_length:"+path)
        for i,(a,b) in enumerate(zip(actual,expected)):compare_science(a,b,path+"["+str(i)+"]")
    elif type(expected) is float:
        require(type(actual) in (float,int) and math.isfinite(actual) and math.isclose(actual,expected,rel_tol=1e-8,abs_tol=1e-12),"scientific_value:"+path)
    else:
        require(type(actual) is type(expected) and actual==expected,"scientific_identity:"+path)


def diagnostic_value(value):
    """Bounded typed value; digest always covers the complete JSON value."""
    raw = canonical_json(value)
    result = dict(type=type(value).__name__, sha256=hashlib.sha256(raw).hexdigest())
    if len(raw) <= 512:
        result["value"] = value
    else:
        result.update(bytes=len(raw), summary="value omitted; complete canonical JSON digest retained")
        if isinstance(value, (dict, list, str)):
            result["length"] = len(value)
    return result


def first_science_difference(actual, expected, path):
    """Explain the unchanged comparator's first failure, in its own ordering."""
    try:
        compare_science(actual, expected, path)
    except ValueError as exc:
        if isinstance(expected, dict) and isinstance(actual, dict) and set(actual) == set(expected):
            for key in expected:
                difference = first_science_difference(actual[key], expected[key], path + "." + key)
                if difference is not None:
                    return difference
        elif isinstance(expected, list) and isinstance(actual, list) and len(actual) == len(expected):
            for index, (a, b) in enumerate(zip(actual, expected)):
                difference = first_science_difference(a, b, path + "[" + str(index) + "]")
                if difference is not None:
                    return difference
        return dict(path=path, category=str(exc).split(":", 2)[1],
                    actual=diagnostic_value(actual), expected=diagnostic_value(expected))
    return None


def first_exact_difference(actual, expected, path="parent_hash_operands"):
    """Hash preimages require exact JSON identity, including float spelling."""
    if canonical_json(actual) == canonical_json(expected):
        return None
    if isinstance(actual, dict) and isinstance(expected, dict) and set(actual) == set(expected):
        for key in sorted(expected):
            difference = first_exact_difference(actual[key], expected[key], path + "." + key)
            if difference is not None:
                return difference
    elif isinstance(actual, list) and isinstance(expected, list) and len(actual) == len(expected):
        for index, (a, b) in enumerate(zip(actual, expected)):
            difference = first_exact_difference(a, b, path + "[" + str(index) + "]")
            if difference is not None:
                return difference
    return dict(path=path, actual=diagnostic_value(actual), expected=diagnostic_value(expected))


def alignment_parent_preimage(alignment):
    """Mirror frozen _alignment_parent_hash for JSON-native saved records only.

    The producer's JSON-ready conversion is the identity on saved JSON values.
    Keep the float conversion and missing-component/operand semantics identical.
    This local operand selection alone does not establish setting-wide binding;
    setting_alignment_preimage owns that selection. Never repair stored output.
    """
    return dict(setting_key=alignment["parent_setting_key"],
                common_shift_gev=float(alignment["common_setting_shift_gev"]),
                resolved_configuration_hash=alignment["resolved_configuration_hash"],
                components={name: {key: ((alignment.get("components") or {}).get(name) or {}).get(key)
                    for key in ("residual_shift_gev", "window_expansion_gev", "effective_windows",
                                "immutable_source_hist_identifier", "immutable_source_hist_checksum")}
                    for name in ("pi_n", "pi_delta", "pi_sidis")})


_ALIGNMENT_COMPONENTS = ("pi_n", "pi_delta", "pi_sidis")
_PARENT_OPERANDS = ("residual_shift_gev", "window_expansion_gev", "effective_windows")


def setting_alignment_preimage(alignment):
    """Recover the setting-wide ancestor, never rehash a t-bin as its owner.

    The frozen producer retains setting-wide shifts/windows as parent_* on
    t-bin components. Identifiers/checksums and common shift remain inherited.
    Missing operands are an error, including on an unaccepted fallback record.
    """
    result = alignment_parent_preimage(alignment)
    require(set(alignment["components"]) == set(_ALIGNMENT_COMPONENTS), "scientific_lineage:components")
    for name in _ALIGNMENT_COMPONENTS:
        component = alignment["components"][name]
        for field in (*_PARENT_OPERANDS, "immutable_source_hist_identifier", "immutable_source_hist_checksum"):
            source = ("parent_" + field if field in _PARENT_OPERANDS
                      and alignment["analysis_scope"] != "setting-wide" else field)
            require(source in component and component[source] is not None,
                    "scientific_lineage:setting_preimage." + name + "." + source)
            result["components"][name][field] = deepcopy(component[source])
    return result


def _check1_hist_name(name, path):
    """Only the two source-proven uuid4().hex[:10] suffixes are erased.

    _build_hist_name owns the first; _clone_detached_hist owns the second.
    Preserve prefix, component, setting and caller context. Stable t/phi names
    are compared verbatim. Unknown generated syntax fails closed.
    """
    require(isinstance(name, str), "scientific_identity_syntax:" + path)
    if name.startswith("H_MM_component_shape"):
        match = re.fullmatch(
            r"(H_MM_component_shape_full_(?:pi_n|pi_delta|pi_sidis)_"
            r"(?:Left|Center|Right)_main_step3_pre_rand_sub)_[0-9a-f]{10}_[0-9a-f]{10}", name)
        require(match is not None, "scientific_identity_syntax:" + path)
        return match.group(1) + "_<generated>_<clone>"
    require(re.search(r"_[0-9a-f]{10}(?:_|$)", name) is None,
            "scientific_identity_syntax:" + path)
    return name


def _check1_identity_view(value, hashes, path):
    """Closed identity edits, on a copy, after raw per-side validation.

    Inventory: alignment hash plus components/applied_components/
    proposed_components name aliases, immutable_source_template_identifier,
    and creation_timestamp. No identifier object or scientific subtree is
    removed. All other keys, including checksums and axes, survive unchanged.
    """
    if isinstance(value, dict):
        result = {}
        for key, item in value.items():
            here = path + "." + key
            if key == "parent_alignment_hash":
                require(".pion_component_alignment." in here and isinstance(item, str) and item in hashes,
                        "scientific_lineage:" + here)
                result[key] = hashes[item]
            elif key in ("hist_name", "immutable_source_hist_name"):
                require(not isinstance(item, str) or not item.startswith("H_MM_component_shape")
                        or ".pion_component_alignment." in here or ".setting_wide_preimage." in here,
                        "scientific_identity_syntax:" + here)
                result[key] = _check1_hist_name(item, here)
            elif key == "creation_timestamp":
                # datetime.now(timezone.utc).isoformat(), in the frozen resolver.
                require(here.endswith(".pion_component_alignment.creation_timestamp") and isinstance(item, str),
                        "scientific_identity_syntax:" + here)
                try:
                    stamp = datetime.fromisoformat(item)
                except ValueError as error:
                    raise ValueError("check1:scientific_identity_syntax:" + here) from error
                require(stamp.tzinfo == timezone.utc and stamp.isoformat() == item,
                        "scientific_identity_syntax:" + here)
                result[key] = "<alignment creation UTC>"
            else:
                result[key] = _check1_identity_view(item, hashes, here)
        return result
    if isinstance(value, list):
        return [_check1_identity_view(v, hashes, path + "[" + str(i) + "]") for i, v in enumerate(value)]
    # Do not let another generated-name occurrence become a late comparison
    # failure or silently normalize it by spelling alone.
    require(not isinstance(value, str) or re.search(r"_[0-9a-f]{10}(?:_|$)", value) is None,
            "scientific_identity_syntax:" + path)
    return value


def _validate_check1_alignment(alignment, ancestor, path):
    require(isinstance(alignment, dict), "scientific_lineage:" + path)
    require(alignment.get("parent_alignment_hash") == ancestor, "scientific_lineage:" + path + ".parent_alignment_hash")
    require(digest(alignment["resolved_config"]) == alignment["resolved_configuration_hash"],
            "scientific_lineage:" + path + ".resolved_configuration_hash")
    for group in ("components", "applied_components", "proposed_components"):
        require(isinstance(alignment.get(group), dict), "scientific_lineage:" + path + "." + group)
        # A fallback intentionally has no proposed components.
        require(set(alignment[group]) == set(_ALIGNMENT_COMPONENTS) or
                (group == "proposed_components" and alignment[group] == {}),
                "scientific_lineage:" + path + "." + group)
        for name, component in alignment[group].items():
            here = path + "." + group + "." + name
            identifier = component.get("immutable_source_hist_identifier")
            require(isinstance(identifier, dict) and identifier.get("component") == name,
                    "scientific_lineage:" + here + ".immutable_source_hist_identifier")
            require(identifier == alignment["immutable_source_template_identifier"][name]
                    and component.get("immutable_source_hist_name") == identifier.get("hist_name"),
                    "scientific_lineage:" + here + ".immutable_source_hist_name")
            require(component.get("immutable_source_hist_checksum") == alignment["immutable_source_template_checksum"][name]
                    and re.fullmatch(r"[0-9a-f]{64}", component["immutable_source_hist_checksum"]) is not None,
                    "scientific_lineage:" + here + ".immutable_source_hist_checksum")
            # The resolver annotates only applied_components (also exposed as
            # components), after selecting the mixed map. Proposals have no
            # inherited link unless they were copied from a parent baseline.
            if not ((group == "proposed_components" and "parent_alignment_hash" not in component)
                    or component.get("parent_alignment_hash") == ancestor):
                error = ValueError("check1:scientific_lineage:" + here + ".parent_alignment_hash")
                error.check1_raw_operands = dict(actual=diagnostic_value(component.get("parent_alignment_hash")),
                    expected=diagnostic_value(ancestor), expectation="same-run setting-wide ancestor")
                raise error


def _validate_check1_candidate_alignment(parent, path):
    """Mirror the frozen producer's candidate_alignment_changed invariant."""
    baseline = parent["baseline_A"]["pion_component_alignment"]
    candidate = parent["candidate_B"]["pion_component_alignment"]
    if candidate != baseline:
        location = path + ".candidate_B.pion_component_alignment"
        error = ValueError("check1:scientific_lineage:" + location)
        error.check1_raw_operands = dict(expectation="same-run baseline_A.pion_component_alignment",
            difference=first_exact_difference(candidate, baseline, location))
        raise error


def check1_scientific_projection(document, key, coordinate=None):
    """Check 1 cross-run view; never used for producer/cache/runtime bindings.

    Validate each side's original ancestry before replacing its raw hash by
    the deterministic digest of its setting-wide scientific preimage. The
    original projection and compare_science semantics remain unchanged.
    """
    try:
        return _check1_scientific_projection(document, key, coordinate)
    except (KeyError, TypeError, IndexError) as error:
        raise ValueError("check1:scientific_lineage:" + key + ":missing_or_invalid_operand:" + str(error)) from error


def validate_check1_science_lineage(document, key):
    """Validate one raw side without constructing any comparison projection."""
    try:
        return _check1_science_lineage(document, key)
    except (KeyError, TypeError, IndexError) as error:
        raise ValueError("check1:scientific_lineage:" + key + ":missing_or_invalid_operand:" + str(error)) from error


def _check1_science_lineage(document, key):
    parents = document.get("parents")
    require(isinstance(parents, list) and len(parents) == 3, "scientific_lineage:" + key + ".parents")
    raw_ancestors = []
    preimages = []
    for index, parent in enumerate(parents):
        path = key + ".parents[" + str(index) + "]"
        alignment = parent["baseline_A"]["pion_component_alignment"]
        preimage = setting_alignment_preimage(alignment)
        ancestor = alignment.get("parent_alignment_hash")
        require(isinstance(ancestor, str) and digest(preimage) == ancestor,
                "scientific_lineage:" + path + ".setting_wide_preimage")
        preimages.append(preimage)
        raw_ancestors.append(ancestor)
        _validate_check1_alignment(alignment, ancestor, path + ".baseline_A")
        _validate_check1_candidate_alignment(parent, path)
        _validate_check1_alignment(parent["candidate_B"]["pion_component_alignment"], ancestor, path + ".candidate_B")
    require(all(p == preimages[0] for p in preimages), "scientific_lineage:" + key + ".setting_wide_ancestry")
    if "cells" in document:
        require(isinstance(document["cells"], list) and len(document["cells"]) == 27,
                "scientific_lineage:" + key + ".cells")
        for index, cell in enumerate(document["cells"]):
            alignment = cell.get("fit", {}).get("pion_component_alignment")
            require(alignment is not None or cell.get("fit", {}).get("status") == "unavailable",
                    "scientific_lineage:" + key + ".cells[" + str(index) + "].fit.pion_component_alignment")
            if alignment is not None:
                ti = cell.get("t_index")
                require(type(ti) is int and 0 <= ti < 3, "scientific_lineage:" + key + ".t_index")
                _validate_check1_alignment(alignment, raw_ancestors[ti], key + ".cells[" + str(index) + "].fit")
    return raw_ancestors, preimages


def _check1_scientific_projection(document, key, coordinate):
    raw_ancestors, preimages = validate_check1_science_lineage(document, key)
    hashes = {}
    # All raw parents/children must pass before any deterministic view exists.
    # Exact preimage inputs stay exact, as with the original hash strings.
    for index, (ancestor, preimage) in enumerate(zip(raw_ancestors, preimages)):
        hashes[ancestor] = digest(_check1_identity_view(preimage, {},
            key + ".parents[" + str(index) + "].setting_wide_preimage"))
    return _check1_identity_view(scientific_projection(document, coordinate), hashes, key)


def validate_check1_parent_binding(child, parent, key):
    """The separate artifact and embedded t parents are the same raw run."""
    missing = object()
    def operand(value):
        return dict(status="missing") if value is missing else dict(status="available", **diagnostic_value(value))
    def reject(location, artifact, index, branch, embedded, separate, reason, **context):
        error = ValueError("check1:scientific_lineage:" + location)
        error.check1_raw_operands = dict(artifact=artifact, parent_index=index, branch=branch, reason=reason,
            embedded=operand(embedded), separate=operand(separate), **context)
        if reason == "raw_alignment_mismatch":
            difference = first_exact_difference(embedded, separate,
                key + ".parents[" + str(index) + "]." + branch + ".pion_component_alignment")
            error.check1_raw_operands["difference"] = dict(path=difference["path"],
                embedded=difference["actual"], separate=difference["expected"])
        raise error
    documents = (child, parent)
    for artifact, document in zip(("embedded", "separate"), documents):
        if not isinstance(document, dict):
            reject(key, artifact, None, None, child, parent, "document_mapping_required")
    lists = (child.get("parents", missing), parent.get("parents", missing))
    for artifact, records in zip(("embedded", "separate"), lists):
        if not isinstance(records, list):
            reject(key + ".parents", artifact, None, None, *lists, "parent_list_required")
        if len(records) != 3:
            index = min(len(records), 3)
            values = [items[index] if index < len(items) else missing for items in lists
                      if isinstance(items, list)]
            if len(values) != 2:values = list(lists)
            reject(key + ".parents[" + str(index) + "]", artifact, index, None, *values,
                   "three_parents_required", embedded_parent_count=len(lists[0]),
                   separate_parent_count=len(lists[1]) if isinstance(lists[1], list) else None)
    for index in range(3):
        embedded, separate = lists[0][index], lists[1][index]
        location = key + ".parents[" + str(index) + "]"
        for artifact, record in (("embedded", embedded), ("separate", separate)):
            if not isinstance(record, dict):
                reject(location, artifact, index, None, embedded, separate, "parent_mapping_required")
        for branch in ("baseline_A", "candidate_B"):
            branches = (embedded.get(branch, missing), separate.get(branch, missing))
            for artifact, record in zip(("embedded", "separate"), branches):
                if not isinstance(record, dict):
                    reject(location + "." + branch, artifact, index, branch, *branches, "branch_mapping_required")
            alignments = tuple(record.get("pion_component_alignment", missing) for record in branches)
            for artifact, alignment in zip(("embedded", "separate"), alignments):
                if not isinstance(alignment, dict) or not alignment:
                    reject(location + "." + branch + ".pion_component_alignment", artifact, index, branch,
                           *alignments, "alignment_mapping_required")
            if alignments[0] != alignments[1]:
                reject(location + "." + branch + ".raw_artifact_binding", "embedded/separate", index, branch,
                       *alignments, "raw_alignment_mismatch")
        for record in (embedded, separate):
            _validate_check1_candidate_alignment(record, key + ".parents[" + str(index) + "]")


def alignment_identity_record(alignment, component):
    if not isinstance(alignment, dict):
        return dict(status="missing")
    selection = {key: diagnostic_value(alignment.get(key)) for key in (
        "accepted", "source", "analysis_scope", "persistence_status", "persistence_rejection_reasons")}
    components = alignment.get("applied_components")
    applied = components.get(component) if isinstance(components, dict) else None
    applied = applied if isinstance(applied, dict) else {}
    selection["applied_component"] = {key: diagnostic_value(applied.get(key)) for key in (
        "source", "component_applied_source", "component_fallback_used", "component_fallback_reason",
        "parent_alignment_hash")}
    return dict(status="available", parent_alignment_hash=diagnostic_value(alignment.get("parent_alignment_hash")),
                selection=selection)


_DIAGNOSTIC_UNSET = object()


def _check1_raw_operand(document, key, location):
    """Resolve only an admitted in-memory operand; never read an external path."""
    if not location.startswith(key + "."):
        return dict(status="unavailable", reason="invariant has no complete document location")
    path = location[len(key) + 1:]
    def value_at(value, operand_path):
        for token in re.findall(r"[^.\[\]]+|\[\d+\]", operand_path):
            value = value[int(token[1:-1])] if token.startswith("[") else value[token]
        return value
    try:
        if ".setting_wide_preimage" in path:
            head, tail = path.split(".setting_wide_preimage", 1)
            # The parent itself is intentionally not dumped into a receipt.
            index = int(re.fullmatch(r"parents\[(\d+)\]", head).group(1))
            document = setting_alignment_preimage(document["parents"][index]["baseline_A"]["pion_component_alignment"])
            path = tail.lstrip(".")
        try:
            value = value_at(document, path)
        except (KeyError, TypeError, IndexError):
            # Raw guard names omit this wrapper; syntax-gate paths already
            # describe the stored JSON and must be tried verbatim first.
            path = re.sub(r"\.(baseline_A|candidate_B|fit)\.(?!pion_component_alignment)",
                          r".\1.pion_component_alignment.", "." + path)[1:]
            value = value_at(document, path)
        return dict(status="available", path=location, operand=diagnostic_value(value))
    except (KeyError, TypeError, IndexError, ValueError, AttributeError):
        return dict(status="unavailable", reason="rejected operand missing or not retained", path=location)


def scientific_failure_diagnostic(actual, expected, source_sha, reference_sha, key, failure,
                                  actual_science=_DIAGNOSTIC_UNSET, expected_science=_DIAGNOSTIC_UNSET,
                                  raw_rejection=None, raw_documents=None):
    """Compact first-failure receipt using only admitted in-memory child records."""
    if raw_rejection is None:
        a = scientific_projection(actual) if actual_science is _DIAGNOSTIC_UNSET else actual_science
        b = scientific_projection(expected) if expected_science is _DIAGNOSTIC_UNSET else expected_science
        difference = first_science_difference(a, b, key)
        context_path = difference["path"] if difference else ""
    else:
        # No valid pair of comparison views exists. A raw-document mismatch
        # (for example a legitimate regenerated UUID) is not the rejected gate.
        difference = None
        context_path = failure.split(":", 2)[2]
    result = dict(schema_version="check1_scientific_failure/v1", source_sha=source_sha,
                  reference_sha256=reference_sha, setting=key.split(":")[0], epsilon=key.split(":")[1],
                  invariant=failure, first_difference=difference, diagnostic_only=True)
    if raw_rejection is not None:
        side = raw_rejection["side"]
        docs = raw_documents if raw_documents is not None else (actual, expected)
        result["projected_first_difference"] = dict(status="unavailable",
            reason="raw lineage or identity syntax rejected before comparison projections completed")
        result["raw_rejection"] = dict(invariant=failure, path=context_path, side=side)
        if "operands" in raw_rejection:
            # Guard-supplied operands describe the actual within-run binding,
            # including embedded/separate records, not the other run's value.
            result["raw_rejection"]["operands"] = raw_rejection["operands"]
        else:
            result["raw_rejection"]["rejected_operand"] = _check1_raw_operand(
                docs[0 if side == "actual" else 1], key, context_path)
    match = re.search(r"\.(cells|parents)\[(\d+)\]", context_path)
    if match is None:
        result["alignment_context"] = ("unavailable: rejected location does not identify a child cell or t parent"
            if raw_rejection is not None else "unavailable: first difference does not identify a child cell or t parent")
        return result
    index = int(match.group(2))
    component_match = re.search(r"\.(?:applied_components|components)\.(pi_n|pi_delta|pi_sidis)(?:\.|$)", context_path)
    component = component_match.group(1) if component_match else None
    result.update(component=component)
    if match.group(1) == "cells":result["cell_index"] = index
    preimages = {}
    for side, doc in (("actual", actual), ("expected", expected)):
        doc = doc if isinstance(doc, dict) else {}
        cells = doc.get("cells")
        cells = cells if isinstance(cells, list) else []
        cell = cells[index] if match.group(1) == "cells" and index < len(cells) else {}
        cell = cell if isinstance(cell, dict) else {}
        ti = cell.get("t_index") if match.group(1) == "cells" else index
        parents = doc.get("parents")
        parents = parents if isinstance(parents, list) else []
        parent = parents[ti] if type(ti) is int and 0 <= ti < len(parents) else {}
        parent = parent if isinstance(parent, dict) else {}
        baseline = parent.get("baseline_A")
        fit = cell.get("fit")
        pa = baseline.get("pion_component_alignment") if isinstance(baseline, dict) else None
        ca = fit.get("pion_component_alignment") if isinstance(fit, dict) else None
        pa = pa if isinstance(pa, dict) else None
        ca = ca if isinstance(ca, dict) else None
        applied = ca.get("applied_components") if ca is not None else None
        applied = applied.get(component) if isinstance(applied, dict) else None
        info = dict(t_index=ti, phi_index=cell.get("phi_index"),
                    parent_record="children[" + key + "].parents[" + str(ti) + "].baseline_A.pion_component_alignment",
                    parent=alignment_identity_record(pa, component), child=alignment_identity_record(ca, component))
        hashes = [pa.get("parent_alignment_hash") if isinstance(pa, dict) else None,
                  ca.get("parent_alignment_hash") if isinstance(ca, dict) else None,
                  applied.get("parent_alignment_hash") if isinstance(applied, dict) else None]
        info["lineage_consistency"] = ("missing" if any(h is None for h in hashes) else
                                        "consistent" if len(set(hashes)) == 1 else "inconsistent")
        info["parent_hash_consistency"] = dict(status="not_verified", reason="parent alignment unavailable")
        if isinstance(pa, dict):
            try:
                local_preimage = alignment_parent_preimage(pa)
                computed = digest(local_preimage)
                info["t_bin_local_rehash"] = dict(scope="t-bin local components; not the setting-wide owner",
                                                  recomputed_hash=computed)
                if isinstance(ca, dict):
                    predicate = deepcopy(local_preimage)
                    predicate["common_shift_gev"] = float(ca["common_setting_shift_gev"])
                    info["child_predicate_rehash"] = dict(scope="t-bin components with child invocation shift",
                                                           recomputed_hash=digest(predicate))
                try:
                    recovered = setting_alignment_preimage(pa)
                except (KeyError, TypeError, ValueError):
                    info["parent_hash_consistency"] = dict(status="not_verified",
                        reason="incomplete setting-wide parent-owned operands; local rehash is not an integrity test")
                    operands = local_preimage
                    info["parent_hash_operand_scope"] = "t-bin local only; setting-wide operands not verified"
                else:
                    preimages[side] = recovered
                    operands = recovered
                    info["parent_hash_operand_scope"] = "setting-wide ancestor"
                    computed = digest(recovered)
                    info["parent_hash_consistency"] = dict(status="matches" if computed == hashes[0] else "mismatch",
                        scope="setting-wide ancestor from parent-owned operands", recomputed_hash=computed)
                info["parent_hash_operands"] = {k: diagnostic_value(v) for k, v in operands.items() if k != "components"}
                info["parent_hash_operands"]["components"] = {name: {k: diagnostic_value(v) for k, v in values.items()}
                    for name, values in operands["components"].items()}
            except (KeyError, TypeError, ValueError):
                info["parent_hash_consistency"] = dict(status="not_verified", reason="incomplete or invalid parent preimage")
        result[side] = info
    result["first_parent_operand_difference"] = (dict(status="unavailable",
        reason="a cross-run operand difference is not the cause of a raw validation rejection")
        if raw_rejection is not None else first_exact_difference(preimages["actual"], preimages["expected"])
        if len(preimages) == 2 else dict(status="not_verified", reason="both parent preimages required"))
    # The owner bounds this value using its actual indented persistence format.
    # Compact JSON size cannot bound nested indentation or text-mode newlines.
    return result


def epsilon(inp):
    value = str(inp.get("EPSSET", "")).lower()
    require(value in ("low", "lowe", "high", "highe"), "epsilon")
    return "lowe" if value in ("low", "lowe") else "highe"


def requested(inp, setting, iteration=0, environ=None):
    env = os.environ if environ is None else environ
    return (env.get(ENABLE_ENV) == "1" and iteration == 0
            and not inp.get("bg_opt_prepass_mode")
            and str(inp.get("analysis_subtraction_context", "")).strip() != "bg_optimizer_pre_parent"
            and inp.get("ParticleType") == "kaon" and str(inp.get("Q2")) == "4p4"
            and str(inp.get("W")) == "2p74"
            and str(inp.get("EPSSET", "")).lower() in ("low", "lowe", "high", "highe")
            and (setting, epsilon(inp)) in SETTINGS)


def axis(values, label):
    a = np.asarray(values)
    require(a.ndim == 1 and a.dtype.kind in "fiu" and len(a) >= 2
            and np.isfinite(a).all() and np.all(a[1:] > a[:-1]), "axis:" + label)
    return np.asarray(a, dtype="<f8")


def member(value, edges, final_inclusive=False):
    i = bisect_right(edges, value) - 1
    if final_inclusive and value == edges[-1]:
        return len(edges) - 2
    return i if 0 <= i < len(edges) - 1 else -1


def close(actual, expected, label):
    a, b = np.asarray(actual), np.asarray(expected)
    require(a.shape == b.shape and a.dtype.kind in "fiu" and b.dtype.kind in "fiu"
            and np.isfinite(a).all() and np.isfinite(b).all()
            and np.allclose(a, b, rtol=1e-8, atol=1e-12), "closure:" + label)


def authority(inp, setting, outpath, environ=None):
    env = os.environ if environ is None else environ
    p = Path(env[AUTHORITY_ENV])
    require(p.is_file() and not p.is_symlink() and file_hash(p) == env[AUTHORITY_HASH_ENV], "authority_hash")
    doc = strict_json(p.read_text(encoding="utf-8"))
    require(doc["schema_version"]=="check1_owner_authority/v1","authority_schema")
    require(doc["reference_sha256"] == REFERENCE_SHA256
            and doc["settings"] == [list(s) for s in SETTINGS]
            and doc["private_output"] == str(Path(outpath).resolve())
            and (setting, epsilon(inp)) in SETTINGS, "authority_binding")
    root = Path(__file__).resolve().parents[2]
    for name in EXECUTABLES:
        require(file_hash(root / name) == doc["executables"][name]["sha256"], "executable_identity:" + name)
    return doc


def metadata(inp, setting, kind, axes, bindings, *, source_sha, attempt, executables,
             reference_sha256=REFERENCE_SHA256, lineage=None):
    require((setting, epsilon(inp)) in SETTINGS and kind in ("data", "simc"), "setting_or_kind")
    require(re.fullmatch("[0-9a-f]{40}", source_sha) is not None and attempt, "source_or_attempt")
    require(reference_sha256 == REFERENCE_SHA256 and set(executables) == set(EXECUTABLES), "reference_or_inventory")
    for name, identity in executables.items():
        require(re.fullmatch("[0-9a-f]{64}", identity["sha256"]) is not None
                and type(identity["bytes"]) is int and identity["bytes"] > 0, "executable:" + name)
    require(bindings and all(v.get("original") and v.get("runtime") and v.get("tree")
                            for v in bindings.values()), "source_bindings")
    expected_bindings={"simc:simc"} if kind=="simc" else {p+":"+s for p in ("kaon","pion") for s in SOURCES}
    require(set(bindings)==expected_bindings,"source_binding_inventory")
    for v in bindings.values():
        for role in ("original", "runtime"):
            require(re.fullmatch("[0-9a-f]{64}", v[role]["sha256"]) is not None
                    and type(v[role]["bytes"]) is int and v[role]["bytes"] > 0, "file_role_identity")
    if kind=="simc":
        validate_simc_binding(bindings["simc:simc"],setting,epsilon(inp))
    parsed = {k: axis(v, k).tolist() for k, v in axes.items()}
    require(set(parsed) == {"t", "mm", "native_mm", "canonical_t", "phi"}
            and len(parsed["t"]) == 101 and len(parsed["canonical_t"]) == 4
            and len(parsed["phi"]) == 10
            and parsed["canonical_t"][0] == parsed["t"][0] == float(inp["tmin"])
            and parsed["canonical_t"][-1] == parsed["t"][-1] == float(inp["tmax"]), "accepted_interval")
    require([float(inp["mm_min"]), float(inp["mm_max"])] == [1.10, 1.16], "signal_window")
    return dict(schema_version=SCHEMA, kind=kind, Q2="4p4", W="2p74", setting=setting,
                epsilon=epsilon(inp), iteration=0, branch="detached_Local_t", diagnostic_only=True,
                production_applied=False, scientifically_accepted=False, source_sha=source_sha,
                attempt=attempt, reference_sha256=reference_sha256, executables=deepcopy(executables),
                executable_inventory_sha256=digest(executables), source_bindings_sha256=digest(bindings),
                source_bindings=deepcopy(bindings), axes=parsed, lineage=deepcopy(lineage or {}),
                signal_window=[1.10, 1.16], units={"t":"GeV^2", "mm":"GeV", "phi":"degree",
                "histograms":"integrated bin weight; flows included"},
                axis_order=["source", "t_with_flows", "MM_with_flows"],
                covariance="conditional fixed-weight Sumw2; fitted/training/shared-systematic covariance unavailable",
                status="captured_pending_owner_verification")


ROW_DTYPES = {"source":"<i4", "entry":"<i8", "tree":"<U128", "t":"<f8", "mm":"<f8",
              "phi":"<f8", "ti":"<i4", "pi":"<i4", "allcuts":"?", "nommcuts":"?",
              "coefficient":"<f8", "weight":"<f8", "sumw2":"<f8"}
EXTRA_DTYPES = {"kaon":{"proton_factor":"<f8", "rf_accept":"?"},
                "pion":{"transfer":"<f8", "transfer_bin":"<i4", "transfer_fingerprint":"<U64"},
                "simc":{"iter_weight":"<f8", "normfac":"<f8"}}


def columns(rows, population):
    dtypes = dict(ROW_DTYPES, **EXTRA_DTYPES[population])
    ordered = sorted(rows, key=lambda r:(r["source"], r["tree"], r["entry"]))
    result = {}
    for key, dtype in dtypes.items():
        values = [r[key] for r in ordered]
        if np.dtype(dtype).kind in "iu":
            require(all(type(v) is int for v in values), "integer_column:" + key)
        elif np.dtype(dtype).kind == "b":
            require(all(type(v) is bool for v in values), "bool_column:" + key)
        elif np.dtype(dtype).kind == "U":
            require(all(isinstance(v,str) and len(v) <= np.dtype(dtype).itemsize // 4 for v in values), "text_column:"+key)
        result[key] = np.asarray(values, dtype=dtype)
    return result


def validate_rows(table, population, meta):
    dtypes = dict(ROW_DTYPES, **EXTRA_DTYPES[population])
    require(set(table) == set(dtypes), "columns:" + population)
    n = len(table["entry"])
    for k, dtype in dtypes.items():
        a = table[k]
        require(isinstance(a,np.ndarray) and a.ndim == 1 and len(a) == n
                and a.dtype == np.dtype(dtype), "dtype_or_length:" + k)
        if a.dtype.kind == "f":
            require(np.isfinite(a).all(), "nonfinite:" + k)
    keys = list(zip(table["source"].tolist(), table["tree"].tolist(), table["entry"].tolist()))
    require(keys == sorted(keys) and len(set(keys)) == n, "duplicate_or_unordered_identity")
    require(np.all(table["entry"] >= 0) and np.all(table["sumw2"] >= 0), "negative_count_or_variance")
    close(table["sumw2"], table["weight"]**2, "row_sumw2")
    require(np.all(~table["allcuts"] | table["nommcuts"]), "selection_flags")
    require(np.all(table["allcuts"] | table["nommcuts"]),"unselected_retained_row")
    for i in range(n):
        s = int(table["source"][i]); tree = str(table["tree"][i])
        key = population + ":" + ("simc" if population == "simc" else SOURCES[s] if 0 <= s < 4 else "invalid")
        require(key in meta["source_bindings"] and tree == meta["source_bindings"][key]["tree"], "row_source_binding")
        inclusive = population == "pion"
        require(int(table["ti"][i]) == member(float(table["t"][i]),meta["axes"]["canonical_t"],inclusive)
                and int(table["ti"][i]) >= 0, "t_membership")
        require(int(table["pi"][i]) == member(float(table["phi"][i]),meta["axes"]["phi"],inclusive), "phi_membership")
    if population == "kaon":
        close(table["weight"], table["coefficient"]*table["proton_factor"]*table["rf_accept"], "proton_RF")
    elif population == "pion":
        require(np.all(table["transfer_bin"] == np.searchsorted(meta["axes"]["mm"], table["mm"], side="right")), "transfer_bin")
        require(all(re.fullmatch("[0-9a-f]{64}",str(v)) for v in table["transfer_fingerprint"]), "transfer_fingerprint")
        close(table["weight"], table["coefficient"]*table["transfer"], "transfer_weight")
    else:
        require(np.all(table["source"] == 0), "simc_source")
        require(np.all(table["normfac"]>0) and (not n or np.all(table["normfac"]==table["normfac"][0])),"simc_normfac")
        close(table["coefficient"],table["normfac"], "simc_coefficient")
        close(table["weight"],table["iter_weight"]*table["normfac"], "normalization_once")


def projections(tables, meta):
    """Independently accumulate joint and marginal counts, sumw, Sumw2, absw."""
    arrays = {}
    ns = 1 if meta["kind"] == "simc" else 4
    te, me = meta["axes"]["t"],meta["axes"]["mm"]
    stages = (("simc", (("simc",1.),)),) if ns == 1 else (
        ("before", (("kaon",1.),)), ("background", (("pion",1.),)),
        ("after", (("kaon",1.),("pion",-1.))))
    for cut in CUTS:
        for stage, contributors in stages:
            base = cut + "/" + stage + "/"
            for moment in ("count","support_count","sumw","sumw2","absw"):
                dtype = "<i8" if moment in ("count","support_count") else "<f8"
                joint = np.zeros((ns,len(te)+1,len(me)+1),dtype=dtype)
                tp = np.zeros((ns,len(te)+1),dtype=dtype)
                mp = np.zeros((ns,len(me)+1),dtype=dtype)
                for pop, sign in contributors:
                    tab = tables[pop]
                    for i in np.flatnonzero(tab[cut]):
                        si = int(tab["source"][i]); ti = bisect_right(te,float(tab["t"][i])); mi = bisect_right(me,float(tab["mm"][i]))
                        w = float(tab["weight"][i])*sign
                        value = {"count":1, "support_count":int(w!=0.), "sumw":w, "sumw2":w*w, "absw":abs(w)}[moment]
                        joint[si,ti,mi] += value
                        tp[si,ti] += value
                        mp[si,mi] += value
                close(joint.sum(axis=2),tp,base+moment+":t")
                close(joint.sum(axis=1),mp,base+moment+":mm")
                arrays[base+moment+"_joint"] = joint
                arrays[base+moment+"_t"] = tp
                arrays[base+moment+"_mm"] = mp
    return arrays


def build_packet(meta, tables):
    require(set(tables) == ({"simc"} if meta["kind"] == "simc" else {"kaon","pion"}), "populations")
    for p,t in tables.items():
        validate_rows(t,p,meta)
    arrays = {"rows/"+p+"/"+k:a for p,t in tables.items() for k,a in t.items()}
    arrays.update(projections(tables,meta))
    slots = []
    for ti in range(3):
        for pi in range(9):
            support = {p:int(np.count_nonzero((t["ti"]==ti)&(t["pi"]==pi)&t["nommcuts"])) for p,t in tables.items()}
            nonzero = {p:int(np.count_nonzero((t["ti"]==ti)&(t["pi"]==pi)&t["nommcuts"]&(t["weight"]!=0.))) for p,t in tables.items()}
            slots.append(dict(t_index=ti,phi_index=pi,support=support,
                              nonzero_contributions=nonzero,
                              status="selected" if any(support.values()) else "empty",
                              reason=("observed_support" if any(nonzero.values()) else "observed_records_zero_weight")
                                     if any(support.values()) else "zero_observed_records"))
    doc = deepcopy(meta)
    doc["slots"] = slots
    doc["parent_slots"] = [dict(t_index=ti,
        support={p:int(np.count_nonzero((t["ti"]==ti)&t["nommcuts"])) for p,t in tables.items()},
        outside_phi={p:int(np.count_nonzero((t["ti"]==ti)&(t["pi"]<0)&t["nommcuts"])) for p,t in tables.items()},
        status="selected" if any(np.any((t["ti"]==ti)&t["nommcuts"]) for t in tables.values()) else "empty") for ti in range(3)]
    doc["rows"] = {p:len(t["entry"]) for p,t in tables.items()}
    doc["columns"] = {k:{"dtype":a.dtype.str,"shape":list(a.shape)} for k,a in sorted(arrays.items())}
    return doc, dict(sorted(arrays.items()))


def validate_packet(doc, arrays):
    require(doc.get("schema_version") == SCHEMA and doc.get("diagnostic_only") is True
            and doc.get("production_applied") is False and doc.get("scientifically_accepted") is False
            and doc.get("branch") == "detached_Local_t" and doc.get("iteration") == 0
            and doc.get("reference_sha256") == REFERENCE_SHA256, "metadata")
    checked = metadata(dict(EPSSET=doc["epsilon"],tmin=doc["axes"]["t"][0],
                            tmax=doc["axes"]["t"][-1],mm_min=1.10,mm_max=1.16),
                       doc["setting"],doc["kind"],doc["axes"],doc["source_bindings"],
                       source_sha=doc["source_sha"],attempt=doc["attempt"],executables=doc["executables"],
                       reference_sha256=doc["reference_sha256"],lineage=doc["lineage"])
    require(all(doc[k] == v for k,v in checked.items()), "metadata_contract")
    tables = {p:{k:arrays["rows/"+p+"/"+k] for k in dict(ROW_DTYPES,**EXTRA_DTYPES[p])} for p in doc["rows"]}
    rebuilt, expected = build_packet({k:v for k,v in doc.items() if k not in ("slots","rows","columns","npz")},tables)
    require(set(expected) == set(arrays) and rebuilt["columns"] == doc["columns"]
            and rebuilt["slots"] == doc["slots"] and rebuilt["parent_slots"]==doc["parent_slots"]
            and rebuilt["rows"] == doc["rows"], "packet_inventory")
    for k,a in expected.items():
        require(arrays[k].dtype == a.dtype and arrays[k].shape == a.shape, "array_type:"+k)
        if k.split("/")[-1].startswith(("count_","support_count_","sumw2_","absw_")):
            require(np.all(arrays[k]>=0),"negative_moment:"+k)
        if a.dtype.kind in "fiu":
            close(arrays[k],a,k)
        else:
            require(np.array_equal(arrays[k],a), "array_content:"+k)
    return tables


def basename(meta):
    return "{}_kaon_Q4p4W2p74_{}_check1_{}".format(meta["setting"],meta["epsilon"],meta["kind"])


def write_packet(outdir, doc, arrays):
    """Publish NPZ first, metadata last, without replacing a colliding name."""
    outdir = Path(outdir)
    stem = basename(doc)
    jp, npz = outdir/(stem+".json"),outdir/(stem+".npz")
    require(not os.path.lexists(jp) and not os.path.lexists(npz), "destination_exists")
    validate_packet(doc,arrays)
    published = []
    first_failure = None
    try:
        with tempfile.TemporaryDirectory(prefix=".check1-pending-",dir=outdir) as tmp:
            try:
                tn,tj = Path(tmp)/npz.name,Path(tmp)/jp.name
                with tn.open("xb") as stream:
                    np.savez_compressed(stream,**arrays)
                saved = deepcopy(doc)
                saved["npz"] = dict(basename=npz.name,bytes=tn.stat().st_size,sha256=file_hash(tn))
                tj.write_bytes(canonical_json(saved)+b"\n")
                read_packet(tj)
                for source,target in ((tn,npz),(tj,jp)):
                    info=source.stat()
                    os.link(source,target)  # atomic no-clobber; existing symlinks also fail
                    published.append((target,info.st_dev,info.st_ino))
            except BaseException as error:
                first_failure=error
                raise
    except BaseException as failure:
        original=first_failure if first_failure is not None else failure
        rollback_failures=[]
        for p,device,inode in reversed(published):
            try:
                info=p.lstat()
                if not p.is_symlink() and (info.st_dev,info.st_ino)==(device,inode):
                    p.unlink()
            except FileNotFoundError:
                pass
            except BaseException as error:
                rollback_failures.append((str(p),repr(error)))
        if rollback_failures:
            original.check1_rollback_failures=rollback_failures
        if original is failure:
            raise
        raise original from failure
    return saved


def read_packet(path):
    path = Path(path)
    require(path.is_file() and not path.is_symlink(), "metadata_file")
    doc = strict_json(path.read_text(encoding="utf-8"))
    record = doc["npz"]
    require(record["basename"] == basename(doc)+".npz", "npz_basename")
    npz = path.parent / record["basename"]
    require(npz.is_file() and not npz.is_symlink() and npz.stat().st_size == record["bytes"]
            and file_hash(npz) == record["sha256"], "npz_identity")
    with np.load(npz,allow_pickle=False) as loaded:
        require(len(loaded.files) == len(set(loaded.files)), "npz_duplicate_members")
        arrays = {k:loaded[k] for k in loaded.files}
    validate_packet(doc,arrays)
    return doc,arrays


def mm_projection(table,edges,ti,pi,cut,weight):
    values = np.zeros(len(edges)+1); variance = values.copy()
    count = 0
    for i in range(len(table["entry"])):
        if table[cut][i] and table["ti"][i] == ti and (pi is None or table["pi"][i] == pi):
            b = bisect_right(edges,float(table["mm"][i])); w = float(weight[i])
            values[b] += w; variance[b] += w*w; count += 1
    return dict(values=values.tolist(),variances=variance.tolist(),count=count)


def validate_data_closure(tables, observer):
    """Compare every old source/parent/child/outside-phi projection and B application."""
    for pop in ("kaon","pion"):
        tab = tables[pop]
        raw = tab["weight"] if pop == "kaon" else tab["coefficient"]
        for cut in CUTS:
            old = observer["observations"][pop][cut]; edges = old["edges"]
            for ti in range(3):
                checks = [(mm_projection(tab,edges,ti,None,cut,raw),old["parents"][ti])]
                for pi in range(9):
                    checks.append((mm_projection(tab,edges,ti,pi,cut,raw),old["children"][ti][pi]))
                checks.append((mm_projection(tab,edges,ti,-1,cut,raw),old["outside_phi"][ti]))
                for si,s in enumerate(SOURCES):
                    selected = dict(tab, **{cut:tab[cut] & (tab["source"]==si)})
                    checks.append((mm_projection(selected,edges,ti,None,cut,raw),old["sources"][ti][s]))
                for a,b in checks:
                    require(a["count"] == b["count"], "old_projection_count")
                    for k in ("values","variances"):
                        close(a[k],b[k],"old_"+pop+":"+cut+":"+k)
    pion = tables["pion"]
    for cell in observer["cells"]:
        ti,pi = cell["t_index"],cell["phi_index"]
        application = cell["applications"]["Local_t"]
        support = np.count_nonzero((pion["ti"]==ti)&(pion["pi"]==pi)&pion["nommcuts"])
        require(application["status"] == "available" or support == 0, "populated_transfer_unavailable")
        if application["status"] == "available":
            for cut in CUTS:
                projection = mm_projection(pion,observer["observations"]["pion"][cut]["edges"],ti,pi,cut,pion["weight"])
                for k in ("values","variances"):
                    close(projection[k],application["cuts"][cut][k],"Local_t:"+cut+":"+k)


def data_tables(hist, observer, old_rows, kaon_trees):
    """Join continuous cache scalars to the exact existing observer identities."""
    cache = hist["_authoritative_pion_control_source_cache"]
    coeff = observer["source_coefficients"]
    require(set(coeff) == set(SOURCES), "source_coefficients")
    nd,nu,nw = (float(hist[k]) for k in ("normfac_data","normfac_dummy","nWindows"))
    require(all(math.isfinite(v) and v > 0 for v in (nd,nu,nw))
            and coeff == dict(zip(SOURCES,(nd,-nd/nw,-nu,nu/nw))), "source_normalization")
    lookup = {"pion":{},"kaon":{}}
    for r in cache["records"]:
        key = (r["source_label"],int(r["entry_index"]))
        require(key not in lookup["pion"], "duplicate_pion_cache")
        require(r["coordinate_fingerprint"] == cache["coordinate_fingerprint"], "coordinate_fingerprint")
        lookup["pion"][key] = dict(t=float(r["adj_t"]),phi=float(r["phi_degrees"]),
            tree=r["source_tree_name"],cache_expected=dict(mm=float(r["adj_MM"]),ti=int(r["t_index"]),
                pi=-1 if r["phi_index"] is None else int(r["phi_index"]),allcuts=bool(r["allcuts"]),
                nommcuts=bool(r["nommcuts"]),weight=float(r["coefficient"])))
    for section in hist["_yield_data_event_cache"].values():
        n = len(section["adj_MM"])
        fields = ("adj_t","t_index","phi_index","allcuts","nommcuts","phi_shift_deg","source_label",
                  "source_entry_index","proton_cleaning_factor","proton_rf_accept")
        require(all(len(section[k]) == n for k in fields), "kaon_column_lengths")
        for i in range(n):
            s = str(section["source_label"][i]); s = "dummy" if s == "dummy_prompt" else s
            key=(s,int(section["source_entry_index"][i]))
            require(key not in lookup["kaon"], "duplicate_kaon_cache")
            lookup["kaon"][key]=dict(t=float(section["adj_t"][i]),phi=float(section["phi_shift_deg"][i]),
                tree=kaon_trees[s],proton_factor=float(section["proton_cleaning_factor"][i]),
                rf_accept=bool(section["proton_rf_accept"][i]),cache_expected=dict(mm=float(section["adj_MM"][i]),
                    ti=int(section["t_index"][i]),pi=int(section["phi_index"][i]),allcuts=bool(section["allcuts"][i]),
                    nommcuts=bool(section["nommcuts"][i]),weight=coeff[s]*float(section["proton_cleaning_factor"][i])
                    if bool(section["proton_rf_accept"][i]) else 0.))
    transfers = {}
    for ti in range(3):
        cells = [c for c in observer["cells"] if c["t_index"] == ti]
        require(len(cells) == 9 and {c["phi_index"] for c in cells} == set(range(9)), "cell_inventory")
        available = [c["applications"]["Local_t"]["weights"] for c in cells
                     if c["applications"]["Local_t"].get("status") == "available"]
        if available:
            require(all(v == available[0] for v in available), "parent_transfer_disagreement")
            transfers[ti] = available[0]
    tables = {}
    edges = observer["observations"]["pion"]["nommcuts"]["edges"]
    for pop in ("kaon","pion"):
        copied=[]; seen=set()
        for r in old_rows[pop]:
            key=(r["source"],r["entry"])
            require(key in lookup[pop] and key not in seen, "observer_cache_identity")
            seen.add(key)
            extra=dict(lookup[pop][key]);expected=extra.pop("cache_expected")
            require(all((-1 if r[k] is None else r[k])==v for k,v in expected.items()),"cache_scalar_or_flag_mismatch")
            base=dict(source=SOURCES.index(r["source"]),entry=int(r["entry"]),mm=float(r["mm"]),
                ti=int(r["ti"]),pi=-1 if r["pi"] is None else int(r["pi"]),allcuts=bool(r["allcuts"]),
                nommcuts=bool(r["nommcuts"]),coefficient=float(coeff[r["source"]]),**extra)
            if pop == "pion":
                require(r["ti"] in transfers, "populated_transfer_unavailable")
                weights=transfers[r["ti"]]; b=bisect_right(edges,r["mm"])
                require(len(weights) == len(edges)+1 and weights[0] == weights[-1] == 0
                        and all(math.isfinite(w) for w in weights), "transfer_axis")
                w=base["coefficient"]*weights[b]
                base.update(transfer=float(weights[b]),transfer_bin=b,transfer_fingerprint=digest(weights))
            else:
                w=float(r["weight"])
            base.update(weight=w,sumw2=w*w);copied.append(base)
        require(seen == set(lookup[pop]), "cache_rows_dropped")
        tables[pop]=columns(copied,pop)
    validate_data_closure(tables,observer)
    return tables


class SimcSink:
    def __init__(self):
        self.rows = []

    def observe(self, entry, mm, t, phi, iter_weight, allcuts, nommcuts, ti, pi):
        if not (allcuts or nommcuts):
            return
        self.rows.append(dict(entry=int(entry),mm=float(mm),t=float(t),phi=float(phi),
                              iter_weight=float(iter_weight),allcuts=bool(allcuts),nommcuts=bool(nommcuts),
                              ti=-1 if ti is None else int(ti),pi=-1 if pi is None else int(pi)))

    def table(self,normfac):
        return columns([dict(r,source=0,tree="h10",coefficient=float(normfac),normfac=float(normfac),
                             weight=r["iter_weight"]*normfac,sumw2=(r["iter_weight"]*normfac)**2)
                        for r in self.rows],"simc")


def native_simc_closure(table, axes, histograms):
    from pion_background_chain_audit import histogram
    for ti in range(3):
        for pi in range(9):
            h = histogram(histograms["mm"][ti][pi])
            require(h["edges"] == axes["native_mm"], "native_mm_axis")
            p = mm_projection(table,h["edges"],ti,pi,"allcuts",table["weight"])
            close(p["values"],h["values"],"native_simc_content")
            close(p["variances"],h["stored_sumw2"],"native_simc_sumw2")
            th = histogram(histograms["t"][ti][pi])
            require(th["edges"] == axes["t"], "fine_t_axis")
            values = np.zeros(len(axes["t"])+1); variances = values.copy()
            for i in np.flatnonzero(table["allcuts"] & (table["ti"]==ti)&(table["pi"]==pi)):
                b=bisect_right(axes["t"],float(table["t"][i])); w=float(table["weight"][i])
                values[b]+=w; variances[b]+=w*w
            close(values,th["values"],"fine_simc_content")
            close(variances,th["stored_sumw2"],"fine_simc_sumw2")


def source_bindings(auth,setting,eps,pop,trees,runtime=None):
    result={};runtime=runtime or {}
    if pop=="simc":
        require(trees=={"simc":"h10"},"simulation_tree_role")
        roles=simc_derivation(setting,eps)
        archive,hist,model,model_hist=(roles[k] for k in ("archive","archive_hist","model","model_hist"))
        require(auth.get("simc_derivations",{}).get(model)==roles,"simulation_derivation_authority")
        require(all(n in auth["original_inputs"] and n in auth["runtime_before"] for n in (archive,hist))
                and all(n in runtime for n in (model,model_hist)),"simulation_identity_missing")
        identity=lambda record:{k:record[k] for k in ("sha256","bytes")}
        binding=dict(tree="h10",basename=model,original_basename=archive,
            original=identity(auth["original_inputs"][archive]),runtime=identity(runtime[model]),
            archive_runtime=identity(auth["runtime_before"][archive]),
            archive_hist=dict(basename=hist,original=identity(auth["original_inputs"][hist]),
                              runtime=identity(auth["runtime_before"][hist])),
            runtime_hist=dict(basename=model_hist,runtime=identity(runtime[model_hist])),derivation=roles)
        validate_simc_binding(binding,setting,eps)
        return {"simc:simc":binding}
    for source,tree in trees.items():
        role="Dummy" if source in ("dummy","dummy_rand") else "Data"
        name="{}_{}_Analysed_{}_Q4p4W2p74_{}.root".format(setting,pop if pop=="pion" else "kaon",role,eps)
        result[pop+":"+source]=dict(tree=tree,basename=name,original=auth["original_inputs"][name],
                                   runtime=runtime.get(name,auth["runtime_before"][name]))
    return result


def simc_derivation(setting,eps):
    require((setting,eps) in SETTINGS,"simulation_setting")
    archive="Prod_Coin_Q4p4W2p74{}_{}.root".format(setting.lower(),eps)
    model="{}_kaon_Simc_Q4p4W2p74_{}.root".format(setting,eps)
    return dict(archive=archive,archive_hist=archive.replace(".root",".hist"),
                model=model,model_hist=model.replace(".root",".hist"),iteration=0,
                operation="copy archive ROOT and .hist; iter_weight on model ROOT",
                sources=["src/main.py","src/simc_ana/iter_weight.py","src/simc_ana/compare_simc.py"])


def validate_simc_binding(binding,setting,eps):
    roles=simc_derivation(setting,eps)
    require(binding.get("derivation")==roles and binding.get("basename")==roles["model"]
            and binding.get("original_basename")==roles["archive"] and binding.get("tree")=="h10"
            and binding.get("archive_hist",{}).get("basename")==roles["archive_hist"]
            and binding.get("runtime_hist",{}).get("basename")==roles["model_hist"],"simulation_derivation_binding")
    identities=[binding["original"],binding["runtime"],binding.get("archive_runtime",{}),
                binding["archive_hist"].get("original",{}),binding["archive_hist"].get("runtime",{}),
                binding["runtime_hist"].get("runtime",{})]
    for identity in identities:
        require(set(identity)=={"sha256","bytes"} and isinstance(identity["sha256"],str)
                and re.fullmatch("[0-9a-f]{64}",identity["sha256"]) is not None
                and type(identity["bytes"]) is int and identity["bytes"]>0,"simulation_file_role_identity")
    require(binding["original"]==binding["archive_runtime"]
            and binding["archive_hist"]["original"]==binding["archive_hist"]["runtime"]
            ==binding["runtime_hist"]["runtime"],"simulation_archive_copy_identity")


def export_data(hist,inp,observer,old_rows,outpath):
    from pion_background_chain_audit import histogram,snapshot
    from prompt_trees import get_prompt_tree_name,get_rand_tree_name
    observer=snapshot(observer)
    setting=hist["phi_setting"];eps=epsilon(inp);auth=authority(inp,setting,outpath)
    require(observer["status"]=="captured_pending_independent_validation"
            and observer["source_identity_complete"] and observer["production_parents_unchanged"]
            and observer["source_HEAD"]==auth["source_sha"],"observer_not_successful")
    coordinate=scientific_projection(observer["coordinate"])
    validate_coordinate_contract(observer)
    require(digest(coordinate)==auth["coordinates"][setting+":"+eps]
            and observer["selection_provenance"]["physical_pion_control_mask_fingerprint"]==auth["pid_fingerprints"][setting+":"+eps],"reference_coordinate_or_PID")
    active=bool((hist.get("_proton_contamination_cleaning_result_setting") or {}).get("accepted"))
    trees={s:(get_rand_tree_name if s in ("rand","dummy_rand") else get_prompt_tree_name)(
        "kaon",inp["EPSSET"],**({"rf_state":"noRF"} if active else {})) for s in SOURCES}
    tables=data_tables(hist,observer,old_rows,trees)
    runtime={}
    for role,identity in observer["runtime_inputs"].items():
        name=Path(identity["path"]).name
        require(name in auth["original_inputs"],"runtime_container_role")
        runtime[name]={k:identity[k] for k in ("sha256","bytes")}
    bindings=source_bindings(auth,setting,eps,"kaon",trees,runtime)
    pion_trees={s:hist["_authoritative_pion_control_source_cache"]["source_accounting"][s]["tree_name"] for s in SOURCES}
    bindings.update(source_bindings(auth,setting,eps,"pion",pion_trees))
    axes=dict(t=histogram(hist["H_t_DATA"])["edges"],mm=observer["observations"]["pion"]["nommcuts"]["edges"],
              native_mm=auth["native_mm"][eps],canonical_t=observer["t_edges"],phi=observer["phi_edges"])
    lineage=dict(coordinate=observer["coordinate"],
        reference_source_sha=auth["reference_source_sha"][setting+":"+eps],
        reference_coordinate_fingerprint=auth["reference_coordinate_fingerprints"][setting+":"+eps],
        selection_provenance=observer["selection_provenance"],
        source_coefficients=observer["source_coefficients"],selected_row_fingerprints=observer["selected_row_fingerprints"],
        membership={"kaon":"yield-loop half-open","pion":"canonical final endpoint inclusive"},
        original_child_json="{}_kaon_Q4p4W2p74_{}_pion-t-phi-template-comparison.json".format(setting,eps),
        parent_fingerprints=[digest(p) for p in observer["parents"]],
        cells=[dict(t_index=c["t_index"],phi_index=c["phi_index"],status=c["status"],reason=c.get("reason"),
                    transfer_status=c["applications"]["Local_t"]["status"],
                    transfer_fingerprint=digest(c["applications"]["Local_t"])) for c in observer["cells"]])
    child_path=Path(outpath)/lineage["original_child_json"]
    lineage["original_child_json_identity"]=dict(sha256=file_hash(child_path),bytes=child_path.stat().st_size)
    meta=metadata(inp,setting,"data",axes,bindings,source_sha=auth["source_sha"],attempt=auth["attempt"],
                  executables=auth["executables"],lineage=snapshot(lineage))
    return write_packet(outpath,*build_packet(meta,tables))


def export_simc(sink,tree_file,normfac,inp,setting,t_edges,phi_edges,histograms,outpath):
    from pion_background_chain_audit import histogram
    auth=authority(inp,setting,outpath);eps=epsilon(inp)
    axes=dict(t=histogram(histograms["t"][0][0])["edges"],mm=auth["broad_mm"][eps],
              native_mm=histogram(histograms["mm"][0][0])["edges"],canonical_t=list(t_edges),phi=list(phi_edges))
    table=sink.table(float(normfac));native_simc_closure(table,axes,histograms)
    roles=simc_derivation(setting,eps);path=Path(tree_file.GetName())
    require(path.name==roles["model"] and path.parent.resolve()==Path(outpath).resolve(),"simulation_runtime_role")
    runtime={}
    for name in (roles["model"],roles["model_hist"]):
        consumed=Path(outpath)/name
        require(consumed.is_file() and not consumed.is_symlink(),"simulation_runtime_file:"+name)
        runtime[name]=dict(sha256=file_hash(consumed),bytes=consumed.stat().st_size)
    bindings=source_bindings(auth,setting,eps,"simc",{"simc":"h10"},runtime)
    meta=metadata(inp,setting,"simc",axes,bindings,source_sha=auth["source_sha"],attempt=auth["attempt"],
                  executables=auth["executables"],lineage=dict(normfac=float(normfac),
                    reference_source_sha=auth["reference_source_sha"][setting+":"+eps],
                    coordinate="-evt.t; evt.missmass; evt.phipq*180/pi",membership="yield-loop half-open",
                    native_histogram_closure="PASS"))
    return write_packet(outpath,*build_packet(meta,{"simc":table}))
