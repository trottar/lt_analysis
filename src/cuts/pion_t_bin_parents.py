"""Particle-subtraction-stage ownership for authoritative per-t pion parents.

The numerical t-integrated input builder remains shared with averages so both
paths retain identical signed prompt/random/dummy selections.  This module is
the only normal-production owner of parent identity, freezing, diagnostics,
and PDF rendering.
"""

from __future__ import annotations

import os
import sys

from background_config import (
    resolve_particle_subtraction_mode,
    resolve_pion_subtraction_scope,
)
from pion_component_fits import (
    print_particle_subtraction_component_application_pages,
    print_particle_subtraction_component_fit_pages,
    print_particle_subtraction_kaon_lambda_comparison_page,
    record_particle_subtraction_page,
)
from pion_component_subtraction import (
    build_simc_shape_pion_control_weights,
    build_t_bin_pion_parent_identity,
    validate_authoritative_t_bin_pion_parent,
    validate_frozen_t_bin_pion_parent_collection,
)
from root_histogram_ownership import clone_root_histogram


def _is_component_t_bin_production(inp_dict):
    return (
        str((inp_dict or {}).get("ParticleType", "")).strip().lower() == "kaon"
        and resolve_particle_subtraction_mode(inp_dict) == "simc_shape_components"
        and resolve_pion_subtraction_scope(inp_dict) == "t_bin"
    )


def _shared_t_integrated_fit_records(*args, **kwargs):
    """Use the established signed-source histogram builder without ownership.

    Import lazily so ``rand_sub`` never imports ``ave_per_bin`` directly and
    so the generic averages module cannot become a particle-parent producer.
    """
    binning_path = os.path.normpath(
        os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "binning")
    )
    if binning_path not in sys.path:
        sys.path.append(binning_path)
    from ave_per_bin import process_hist_data

    return process_hist_data(*args, **kwargs)


def _diagnostic_application_from_fit_model(fit_result, *, t_index):
    """Build a detached parent-level model application for diagnostics only.

    The fit remains authoritative even if these optional display clones cannot
    be built.  The payload is deliberately never eligible for child use.
    """
    try:
        before_source = (
            fit_result.get("H_kaon_nosub_input")
            or fit_result.get("H_kaon_input")
        )
        pion_model = fit_result.get("H_kaon_pion_model")
        if before_source is None or pion_model is None:
            raise RuntimeError("missing parent fit-model histogram")

        before = clone_root_histogram(
            before_source,
            scope="pion_parent_t{}".format(int(t_index) + 1),
            role="diagnostic_before_pion",
            name="H_MM_parent_t{}_before_pion".format(int(t_index) + 1),
        )
        template = clone_root_histogram(
            pion_model,
            scope="pion_parent_t{}".format(int(t_index) + 1),
            role="diagnostic_pion_model",
            name="H_MM_parent_t{}_pion_model".format(int(t_index) + 1),
        )
        after = clone_root_histogram(
            before,
            scope="pion_parent_t{}".format(int(t_index) + 1),
            role="diagnostic_after_pion",
            name="H_MM_parent_t{}_after_pion".format(int(t_index) + 1),
        )
        after.Add(template, -1.0)

        weights = build_simc_shape_pion_control_weights(fit_result)
        return {
            "accepted": True,
            "fallback_used": False,
            "fallback_reason": "",
            "diagnostic_only": True,
            "application_authoritative": False,
            "production_applied": False,
            "application_context": "authoritative_parent_t_integrated_diagnostic",
            "analysis_scope": fit_result.get("analysis_scope"),
            "H_pion_control_model": weights.get("H_pion_control_model"),
            "H_kaon_pion_model": weights.get("H_kaon_pion_model") or template,
            "H_pion_weight_vs_MM": weights.get("H_pion_weight_vs_MM"),
            "H_pion_subtraction_template_MM": template,
            "H_pion_subtraction_template_MM_nosub": template,
            "H_MM_before_pion_subtraction": before,
            "H_MM_after_pion_subtraction": after,
            "H_MM_nosub_before_pion_subtraction": before,
            "H_MM_nosub_after_pion_subtraction": after,
            "diagnostics": {
                **dict(weights.get("diagnostics") or {}),
                "event_template_closure": {
                    "status": "parent_fit_model_diagnostic",
                    "child_event_templates_remain_downstream_only": True,
                },
            },
        }, {
            "status": "available",
            "mode": "detached_parent_fit_model",
            "event_template_owned_by": "downstream_child_application",
        }
    except Exception as exc:
        return None, {
            "status": "unavailable",
            "reason": "diagnostic_application_clone_failed",
            "detail": str(exc),
        }


def _parent_summary(parent):
    fit_result = parent.get("fit_result") or {}
    kaon_diagnostics = (fit_result.get("diagnostics") or {}).get("kaon") or {}
    return {
        "t_bin": int(parent["t_bin_index"]) + 1,
        "t_edges": list(parent["t_edges"]),
        "pion_parent_id": parent.get("pion_parent_id"),
        "A_n": fit_result.get("A_n"),
        "A_delta": fit_result.get("A_delta"),
        "A_sidis": fit_result.get("A_sidis"),
        "S_lambda": fit_result.get("S_lambda"),
        "chi2_ndf": fit_result.get("chi2_ndf_kaon"),
        "p_value": fit_result.get("fit_p_value_kaon"),
        "parent_fit_status": fit_result.get("fit_status_kaon"),
        "status": fit_result.get("fit_status_kaon"),
        "protected_fit_variant": kaon_diagnostics.get("fit_variant"),
        "protected_fit_attempted": kaon_diagnostics.get("protected_fit_attempted"),
        "protected_fit_succeeded": kaon_diagnostics.get("protected_fit_succeeded"),
        "parent_pi_delta_applied": kaon_diagnostics.get("pi_delta_applied"),
        "parent_failure_reason": kaon_diagnostics.get("failure_reason"),
        "diagnostic_application_status": dict(parent.get("diagnostic_application_status") or {}),
    }


def build_setting_t_bin_pion_parents(
    hist,
    inp_dict,
    *,
    tree_data,
    tree_dummy,
    n_windows,
    t_bins,
    phi_bins,
    kaon_signal_shape_payload,
    kaon_sigma0_shape_payload,
    proton_cleaning_result,
    parent_pion_alignment,
    diagnostic_application_builder=None,
):
    """Fit, validate, and freeze one RF-restored parent per canonical t bin."""
    if not _is_component_t_bin_production(inp_dict):
        return []
    if not isinstance(proton_cleaning_result, dict) or not bool(
        proton_cleaning_result.get("accepted")
    ):
        raise RuntimeError("authoritative_t_bin_pion_parents_require_committed_proton_cleaning")

    phi_setting = hist.get("phi_setting")
    manifest = hist.setdefault("pion_component_page_manifest", [])
    processed_dict = _shared_t_integrated_fit_records(
        tree_data,
        tree_dummy,
        t_bins,
        phi_bins,
        n_windows,
        phi_setting,
        inp_dict,
        particle_subtraction_scale_factor=None,
        kaon_signal_shape_payload=kaon_signal_shape_payload,
        kaon_sigma0_shape_payload=kaon_sigma0_shape_payload,
        proton_cleaning_result=proton_cleaning_result,
        parent_pion_alignment=parent_pion_alignment,
        t_integrated_fit_only=True,
    )

    parents = []
    for t_index in range(len(t_bins) - 1):
        entry = processed_dict.get("t_bin{}".format(t_index + 1)) or {}
        t_edges = [float(t_bins[t_index]), float(t_bins[t_index + 1])]
        fit_result = entry.get("particle_subtraction_component_fit")
        parent = {
            **build_t_bin_pion_parent_identity(inp_dict, phi_setting, t_index, t_edges),
            "t_bin_index": int(t_index),
            "t_edges": t_edges,
            "analysis_scope": "t_bin{}".format(t_index + 1),
            "input_selection": "no_rf_proton_cleaning_then_rf_restored",
            "source_target_state": "post_proton_post_rf",
            "fit_result": fit_result,
            "application_result": None,
            "fit_summary": entry.get("particle_subtraction_component_fit_summary"),
            "application_summary": None,
        }
        # Missing K-Lambda is discovered by the existing fitter and remains a
        # hard error.  Optional detached diagnostic products are independent.
        validate_authoritative_t_bin_pion_parent(
            parent, inp_dict, phi_setting, t_index, t_edges
        )
        if diagnostic_application_builder is None:
            diagnostic_payload, diagnostic_status = _diagnostic_application_from_fit_model(
                fit_result, t_index=t_index
            )
        else:
            try:
                diagnostic_payload = diagnostic_application_builder(
                    fit_result=fit_result,
                    processed_entry=entry,
                    t_index=t_index,
                    t_edges=t_edges,
                )
                diagnostic_status = {
                    "status": "available"
                    if isinstance(diagnostic_payload, dict)
                    else "unavailable",
                    "mode": "detached_parent_event_template",
                    "event_template_owned_by": "parent_diagnostic_application",
                }
                if isinstance(diagnostic_payload, dict) and not diagnostic_payload.get("accepted"):
                    diagnostic_status["status"] = "unavailable"
                    diagnostic_status["reason"] = diagnostic_payload.get("fallback_reason")
            except Exception as exc:
                diagnostic_payload = None
                diagnostic_status = {
                    "status": "unavailable",
                    "reason": "diagnostic_application_clone_failed",
                    "detail": str(exc),
                }
        parent["diagnostic_application_result"] = diagnostic_payload
        parent["diagnostic_application_status"] = diagnostic_status
        parents.append(parent)

    # A tuple makes the authoritative parent *collection* immutable to normal
    # consumers.  The ROOT state inside each validated parent remains owned by
    # that parent for the lifetime of the setting, but Step 6/yield/average
    # code cannot append, remove, or replace parent slots.
    frozen_parents = tuple(parents)
    hist["_pion_t_bin_parent_results"] = frozen_parents
    hist["_pion_t_bin_parent_source"] = "particle_subtraction_stage_rf_restored"
    hist["pion_t_parent_collection_frozen"] = True
    hist["pion_t_amplitude_table"] = [_parent_summary(parent) for parent in frozen_parents]
    setting_fit = hist.get("_particle_subtraction_component_fit_setting") or {}
    hist["pion_t_parent_diagnostics"] = {
        "phi_setting": phi_setting,
        "epsilon": str(inp_dict.get("EPSSET", "")).strip().lower(),
        "setting_wide": {
            "scope": "setting_wide_diagnostic",
            "A_n": setting_fit.get("A_n"),
            "A_delta": setting_fit.get("A_delta"),
            "A_sidis": setting_fit.get("A_sidis"),
            "S_lambda": setting_fit.get("S_lambda"),
            "chi2_ndf": setting_fit.get("chi2_ndf_kaon"),
            "p_value": setting_fit.get("fit_p_value_kaon"),
            "status": setting_fit.get("fit_status_kaon"),
            "diagnostic_only": bool(setting_fit.get("diagnostic_only", True)),
        },
        "parents": list(hist["pion_t_amplitude_table"]),
    }
    validate_frozen_t_bin_pion_parent_collection(hist, inp_dict, t_bins)
    print("[SIMC PION PARENTS] {} {} -- post-proton/RF-restored".format(
        phi_setting, inp_dict.get("EPSSET", "")
    ))
    print("  canonical t bins: {}".format(len(frozen_parents)))
    for row in hist["pion_t_amplitude_table"]:
        print(
            "  t{t_bin} {t_edges} parent={pion_parent_id} status={status}".format(
                **row
            )
        )
    print("  parents frozen before rand_sub PDF close")
    return frozen_parents


def _record_parent_page_if_absent(manifest, page_id, *, scope):
    if any(entry.get("page_id") == page_id for entry in manifest if isinstance(entry, dict)):
        return
    record_particle_subtraction_page(
        manifest, page_id, scope=scope, authoritative=True
    )


def _print_parent_overview_page(pdf_name, parent, title_prefix, manifest, page_prefix):
    """Emit a compact proton-analog parent overview without owning data hists."""
    try:
        import ROOT

        canvas = ROOT.TCanvas(
            "pion_parent_overview_{}".format(parent.get("pion_parent_id", "unknown")[-12:]),
            "Pion parent overview",
            900,
            650,
        )
        text = ROOT.TPaveText(0.10, 0.12, 0.90, 0.88, "NDC")
        text.SetFillStyle(0)
        text.SetBorderSize(0)
        text.AddText("{} AUTHORITATIVE PION PARENT".format(title_prefix))
        text.AddText("t{} [{:.3f}, {:.3f}]".format(
            int(parent["t_bin_index"]) + 1,
            float(parent["t_edges"][0]),
            float(parent["t_edges"][1]),
        ))
        text.AddText("input: post-proton, RF-restored")
        text.AddText("parent: {}".format(parent.get("pion_parent_id")))
        text.AddText("fit status: {}".format((parent.get("fit_result") or {}).get("fit_status_kaon")))
        text.AddText("diagnostic: {}".format((parent.get("diagnostic_application_status") or {}).get("status")))
        text.Draw()
        canvas.Print(pdf_name)
        canvas.Close()
        _record_parent_page_if_absent(
            manifest, "{}.overview".format(page_prefix), scope=parent["analysis_scope"]
        )
        return True
    except Exception as exc:
        parent["diagnostic_application_status"] = {
            **dict(parent.get("diagnostic_application_status") or {}),
            "overview_status": "unavailable",
            "overview_reason": str(exc),
        }
        return False


def render_setting_t_bin_pion_parent_pages(pdf_name, parents, inp_dict, *, title_prefix, page_manifest):
    """Render all authoritative parent sections before ``rand_sub`` returns."""
    for parent in parents:
        fit_result = parent["fit_result"]
        scope = parent["analysis_scope"]
        page_prefix = "pion.t_bin.{}".format(int(parent["t_bin_index"]))
        section_title = "{} t{} [{:.3f}, {:.3f}]".format(
            title_prefix,
            int(parent["t_bin_index"]) + 1,
            float(parent["t_edges"][0]),
            float(parent["t_edges"][1]),
        )
        print_particle_subtraction_component_fit_pages(
            pdf_name,
            fit_result,
            title_prefix=section_title,
            cut_window=(float(inp_dict["mm_min"]), float(inp_dict["mm_max"])),
            page_manifest=page_manifest,
            page_id_prefix=page_prefix,
            authoritative=True,
        )
        payload = parent.get("diagnostic_application_result")
        if isinstance(payload, dict):
            print_particle_subtraction_component_application_pages(
                pdf_name,
                payload,
                title_prefix=section_title,
                cut_window=(float(inp_dict["mm_min"]), float(inp_dict["mm_max"])),
                component_fit_result=fit_result,
                include_lambda_page=False,
                page_manifest=page_manifest,
                page_id_prefix=page_prefix,
                authoritative=True,
            )
        print_particle_subtraction_kaon_lambda_comparison_page(
            pdf_name,
            fit_result,
            payload,
            title_prefix=section_title,
            cut_window=(float(inp_dict["mm_min"]), float(inp_dict["mm_max"])),
            page_manifest=page_manifest,
            page_id_prefix=page_prefix,
            authoritative=True,
        )
        _print_parent_overview_page(
            pdf_name, parent, section_title, page_manifest, page_prefix
        )
        # Application pages retain a few legacy manifest names.  Add the
        # particle-stage aliases only when the corresponding page actually
        # printed, preserving the manifest's successful-output contract.
        existing_ids = {
            entry.get("page_id")
            for entry in page_manifest
            if isinstance(entry, dict)
        }
        aliases = {
            "pion_control_fit": ("pion_control_fit",),
            "protected_fit_or_status": ("protected_fit", "protected_status"),
            "weight": ("pion_weight",),
            "model_closure": ("model_closure",),
            "event_template_closure": ("event_template_closure",),
            "before_after": ("before_after",),
            "lambda_comparison": ("lambda_comparison",),
            "overview": ("overview",),
        }
        for stable_name, legacy_names in aliases.items():
            if any("{}.{}".format(page_prefix, name) in existing_ids for name in legacy_names):
                _record_parent_page_if_absent(
                    page_manifest,
                    "{}.{}".format(page_prefix, stable_name),
                    scope=scope,
                )


__all__ = (
    "build_setting_t_bin_pion_parents",
    "render_setting_t_bin_pion_parent_pages",
    "validate_frozen_t_bin_pion_parent_collection",
)
