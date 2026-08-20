"""Particle-subtraction-stage ownership for authoritative per-t pion parents.

The numerical t-integrated input builder remains shared with averages so both
paths retain identical signed prompt/random/dummy selections.  This module is
the only normal-production owner of parent identity, freezing, diagnostics,
and PDF rendering.
"""

from __future__ import annotations

import os
import sys
import textwrap

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
    evaluate_particle_subtraction_component_fit_result,
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


def _classify_diagnostic_exception(exc):
    detail = str(exc)
    normalized = detail.lower()
    if "missing" in normalized and "hist" in normalized:
        reason = "missing_source_histogram"
    elif "template" in normalized and "binning" in normalized:
        reason = "template_binning_mismatch"
    elif "template" in normalized:
        reason = "weight_build_failed"
    elif "tree" in normalized or "event" in normalized:
        reason = "event_template_failed"
    elif "clone" in normalized:
        reason = "histogram_clone_failed"
    else:
        reason = "proposal_build_failed"
    return {"status": "unavailable", "reason": reason, "detail": detail}


def _normalize_parent_diagnostic_outcome(outcome, fit_result, inp_dict):
    """Normalize legacy and proposal/final diagnostic callback results."""
    production_evaluation = evaluate_particle_subtraction_component_fit_result(
        fit_result,
        inp_dict,
    )
    if isinstance(outcome, dict) and "proposed_diagnostic_application_result" in outcome:
        proposed = outcome.get("proposed_diagnostic_application_result")
        proposed_status = dict(outcome.get("proposed_diagnostic_application_status") or {})
        final = outcome.get("final_diagnostic_application_result")
        final_status = dict(outcome.get("final_diagnostic_application_status") or {})
        production_evaluation = dict(outcome.get("production_evaluation") or production_evaluation)
    else:
        proposed = outcome if isinstance(outcome, dict) else None
        proposed_status = {
            "status": "available" if proposed is not None else "unavailable",
            "reason": None if proposed is not None else "proposal_unavailable",
            "detail": None,
        }
        final = proposed
        final_status = {
            "status": "available" if final is not None else "unavailable",
            "final_status": "applied_component" if final is not None else "unavailable",
            "final_reason": None if final is not None else "proposal_unavailable",
        }

    proposal_available = isinstance(proposed, dict) and bool(proposed.get("proposal_available", True))
    final_available = isinstance(final, dict)
    production_accepted = bool(production_evaluation.get("accepted"))
    rejection_reasons = [
        reason.strip()
        for reason in str(production_evaluation.get("reason") or "").split(";")
        if reason.strip()
    ]
    status = {
        "status": "available" if final_available else ("partial" if proposal_available else "unavailable"),
        "mode": "detached_parent_event_template",
        "event_template_owned_by": "parent_diagnostic_application",
        "proposal_status": proposed_status.get("status") or ("available" if proposal_available else "unavailable"),
        "proposal_reason": proposed_status.get("reason"),
        "proposal_detail": proposed_status.get("detail"),
        "production_evaluation": "accepted" if production_accepted else "rejected",
        "production_rejection_reasons": rejection_reasons,
        "fallback_mode": production_evaluation.get("fallback_mode"),
        "final_status": final_status.get("final_status") or ("applied_component" if final_available else "unavailable"),
        "final_reason": final_status.get("final_reason"),
        "detail": final_status.get("detail"),
    }
    return proposed, proposed_status, final, final_status, production_evaluation, status


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
        "proposal_status": (parent.get("proposed_diagnostic_application_status") or {}).get("status"),
        "production_gate_status": (parent.get("diagnostic_application_status") or {}).get("production_evaluation"),
        "production_gate_reason": "; ".join((parent.get("diagnostic_application_status") or {}).get("production_rejection_reasons") or []),
        "fallback_mode": (parent.get("diagnostic_application_status") or {}).get("fallback_mode"),
        "final_application_status": (parent.get("diagnostic_application_status") or {}).get("final_status"),
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
    hgcer_cutg=None,
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
        hgcer_cutg=hgcer_cutg,
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
            outcome = diagnostic_payload
        else:
            try:
                outcome = diagnostic_application_builder(
                    fit_result=fit_result,
                    processed_entry=entry,
                    t_index=t_index,
                    t_edges=t_edges,
                )
            except Exception as exc:
                outcome = None
                diagnostic_status = _classify_diagnostic_exception(exc)
                if (
                    bool(inp_dict.get("pion_parent_diagnostic_strict", False))
                    and diagnostic_status["reason"] == "proposal_build_failed"
                ):
                    raise
                production_evaluation = evaluate_particle_subtraction_component_fit_result(
                    fit_result,
                    inp_dict,
                )
                parent["proposed_diagnostic_application_result"] = None
                parent["proposed_diagnostic_application_status"] = dict(diagnostic_status)
                parent["final_diagnostic_application_result"] = None
                parent["final_diagnostic_application_status"] = {
                    "status": "unavailable",
                    "final_status": "unavailable",
                    "final_reason": diagnostic_status["reason"],
                    "detail": diagnostic_status["detail"],
                }
                parent["diagnostic_application_result"] = None
                parent["diagnostic_application_status"] = {
                    **diagnostic_status,
                    "proposal_status": "unavailable",
                    "production_evaluation": "accepted" if production_evaluation.get("accepted") else "rejected",
                    "production_rejection_reasons": [
                        reason.strip()
                        for reason in str(production_evaluation.get("reason") or "").split(";")
                        if reason.strip()
                    ],
                    "fallback_mode": production_evaluation.get("fallback_mode"),
                    "final_status": "unavailable",
                    "final_reason": diagnostic_status["reason"],
                }
                parents.append(parent)
                continue
        if diagnostic_application_builder is None:
            # Retain the legacy model-only fallback for callers that do not
            # provide the event-level parent builder.
            proposed = diagnostic_payload
            proposed_status = diagnostic_status
            final = diagnostic_payload
            final_status = {
                "status": diagnostic_status.get("status"),
                "final_status": "applied_component" if diagnostic_payload is not None else "unavailable",
                "final_reason": diagnostic_status.get("reason"),
                "detail": diagnostic_status.get("detail"),
            }
            production_evaluation = evaluate_particle_subtraction_component_fit_result(fit_result, inp_dict)
            diagnostic_status = {
                "status": diagnostic_status.get("status"),
                "mode": diagnostic_status.get("mode", "detached_parent_fit_model"),
                "proposal_status": diagnostic_status.get("status"),
                "proposal_reason": diagnostic_status.get("reason"),
                "proposal_detail": diagnostic_status.get("detail"),
                "production_evaluation": "accepted" if production_evaluation.get("accepted") else "rejected",
                "production_rejection_reasons": [
                    reason.strip()
                    for reason in str(production_evaluation.get("reason") or "").split(";")
                    if reason.strip()
                ],
                "fallback_mode": production_evaluation.get("fallback_mode"),
                "final_status": final_status["final_status"],
                "final_reason": final_status["final_reason"],
                "detail": final_status["detail"],
            }
        else:
            proposed, proposed_status, final, final_status, production_evaluation, diagnostic_status = (
                _normalize_parent_diagnostic_outcome(outcome, fit_result, inp_dict)
            )
        parent["proposed_diagnostic_application_result"] = proposed
        parent["proposed_diagnostic_application_status"] = proposed_status
        parent["final_diagnostic_application_result"] = final
        parent["final_diagnostic_application_status"] = final_status
        parent["production_application_policy"] = {
            "production_evaluation": "accepted" if production_evaluation.get("accepted") else "rejected",
            "fallback_mode": production_evaluation.get("fallback_mode"),
            "rejection_reasons": list(diagnostic_status.get("production_rejection_reasons") or []),
        }
        parent["diagnostic_application_result"] = final
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
            "[SIMC PION PARENT] {setting} {epsilon} t{t_bin} {t_edges} "
            "fit={status} protected={protected} proposal={proposal} "
            "production_gate={gate} final={final} gate_reason={reason} parent={parent}".format(
                setting=phi_setting,
                epsilon=inp_dict.get("EPSSET", ""),
                protected=row.get("protected_fit_variant"),
                proposal=row.get("proposal_status"),
                gate=row.get("production_gate_status"),
                final=row.get("final_application_status"),
                reason=(row.get("production_gate_reason") or "none")[:240],
                parent=str(row.get("pion_parent_id"))[-16:],
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
            850,
        )
        text = ROOT.TPaveText(0.07, 0.05, 0.93, 0.95, "NDC")
        text.SetFillStyle(0)
        text.SetBorderSize(0)
        text.SetTextAlign(12)
        text.SetTextSize(0.024)
        fit_result = parent.get("fit_result") or {}
        status = parent.get("diagnostic_application_status") or {}
        protected = ((fit_result.get("diagnostics") or {}).get("kaon") or {}).get(
            "pi_delta_signal_protected_fit"
        ) or {}
        lambda_retention = protected.get("lambda_retention") or {}
        text.AddText("{} AUTHORITATIVE PION PARENT".format(title_prefix))
        text.AddText("setting/epsilon: {} / {}".format(
            parent.get("phi_setting"), parent.get("epsilon")
        ))
        text.AddText("t{} [{:.3f}, {:.3f}]".format(
            int(parent["t_bin_index"]) + 1,
            float(parent["t_edges"][0]),
            float(parent["t_edges"][1]),
        ))
        text.AddText("input: post-proton, RF-restored")
        text.AddText("source target state: {}".format(parent.get("source_target_state")))
        text.AddText("parent: {}".format(parent.get("pion_parent_id")))
        text.AddText("fit status pion/kaon: {} / {}".format(
            fit_result.get("fit_status_pion"), fit_result.get("fit_status_kaon")
        ))
        text.AddText("production gate: {}  fallback: {}".format(
            status.get("production_evaluation"), status.get("fallback_mode")
        ))
        text.AddText("A_n={}  A_delta proposed={}  A_delta applied={}  A_sidis={}  S_lambda={}".format(
            fit_result.get("A_n"), fit_result.get("A_delta"),
            ((protected.get("protected_applied_amplitudes") or {}).get("pi_delta", fit_result.get("A_delta"))),
            fit_result.get("A_sidis"), fit_result.get("S_lambda"),
        ))
        text.AddText("protected: {}  chi2/ndf={}  p={}".format(
            protected.get("fit_variant") or protected.get("status") or "unavailable",
            protected.get("chi2_ndf"), protected.get("p_value"),
        ))
        text.AddText("Lambda retention: {}".format(
            lambda_retention.get("status") or lambda_retention.get("reason") or "unavailable"
        ))
        text.AddText("proposal: {}  final: {}".format(
            status.get("proposal_status"), status.get("final_status")
        ))
        text.AddText("diagnostic status: {}".format(status.get("status")))
        for reason in status.get("production_rejection_reasons") or []:
            for line in textwrap.wrap("gate reason: {}".format(reason), width=95):
                text.AddText(line)
        for line in textwrap.wrap("diagnostic reason: {}".format(
            status.get("final_reason") or status.get("proposal_reason") or "none"
        ), width=95):
            text.AddText(line)
        for line in textwrap.wrap("detail: {}".format(status.get("detail") or status.get("proposal_detail") or "none"), width=95):
            text.AddText(line)
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


def _print_parent_application_gate_page(pdf_name, parent, title_prefix, manifest, page_prefix):
    """Print the production-gate decision without conflating it with proposal availability."""
    try:
        import ROOT

        fit_result = parent.get("fit_result") or {}
        status = parent.get("diagnostic_application_status") or {}
        canvas = ROOT.TCanvas(
            "pion_parent_gate_{}".format(parent.get("pion_parent_id", "unknown")[-12:]),
            "Pion parent application gate",
            900,
            650,
        )
        text = ROOT.TPaveText(0.08, 0.08, 0.92, 0.92, "NDC")
        text.SetFillStyle(0)
        text.SetBorderSize(0)
        text.SetTextAlign(12)
        text.SetTextSize(0.026)
        text.AddText("{} application gate".format(title_prefix))
        text.AddText("component production evaluation: {}".format(status.get("production_evaluation")))
        text.AddText("fit status pion/kaon: {} / {}".format(
            fit_result.get("fit_status_pion"), fit_result.get("fit_status_kaon")
        ))
        text.AddText("proposal: {}  final: {}  fallback: {}".format(
            status.get("proposal_status"), status.get("final_status"), status.get("fallback_mode")
        ))
        text.AddText("pion validation: {}  kaon validation: {}".format(
            (fit_result.get("diagnostics") or {}).get("pion", {}).get("validation", {}).get("accepted"),
            (fit_result.get("diagnostics") or {}).get("kaon", {}).get("validation", {}).get("accepted"),
        ))
        reasons = status.get("production_rejection_reasons") or ["none"]
        for reason in reasons:
            for line in textwrap.wrap("reason: {}".format(reason), width=100):
                text.AddText(line)
        text.Draw()
        canvas.Print(pdf_name)
        canvas.Close()
        _record_parent_page_if_absent(
            manifest, "{}.application_gate".format(page_prefix), scope=parent["analysis_scope"]
        )
        return True
    except Exception as exc:
        parent["diagnostic_application_status"] = {
            **dict(parent.get("diagnostic_application_status") or {}),
            "application_gate_page_status": "unavailable",
            "application_gate_page_detail": str(exc),
        }
        return False


def _print_parent_proposal_final_pages(pdf_name, parent, title_prefix, manifest, page_prefix):
    """Additive proposed/final weight and before-after diagnostic pages."""
    try:
        import ROOT

        proposal = parent.get("proposed_diagnostic_application_result") or {}
        final = parent.get("final_diagnostic_application_result") or {}
        before = proposal.get("H_MM_nosub_before_pion_subtraction") or final.get("H_MM_nosub_before_pion_subtraction")
        proposed_after = proposal.get("H_MM_nosub_after_pion_subtraction")
        final_after = final.get("H_MM_nosub_after_pion_subtraction")
        proposed_weight = proposal.get("H_pion_weight_vs_MM")
        final_weight = final.get("H_pion_weight_vs_MM")
        if before is None and proposed_weight is None and final_weight is None:
            return False

        canvas = ROOT.TCanvas(
            "pion_parent_transition_{}".format(parent.get("pion_parent_id", "unknown")[-12:]),
            "Pion parent proposal/final",
            950,
            800,
        )
        canvas.Divide(1, 2)
        canvas.cd(1)
        legend_weight = ROOT.TLegend(0.56, 0.70, 0.90, 0.90)
        display_hists = []
        weight_drawn = False
        for histogram, label, color, style in (
            (proposed_weight, "proposed component weight", ROOT.kRed + 1, 2),
            (final_weight, "final applied weight", ROOT.kGreen + 2, 1),
        ):
            if histogram is None:
                continue
            display = clone_root_histogram(
                histogram,
                scope="{}".format(page_prefix),
                role="proposal_final_weight_display",
            )
            display.SetLineColor(color)
            display.SetLineStyle(style)
            display.SetTitle("{} proposed versus final pion weight;MM [GeV];weight".format(title_prefix))
            display.Draw("hist" if not weight_drawn else "hist same")
            display_hists.append(display)
            legend_weight.AddEntry(display, label, "l")
            weight_drawn = True
        if weight_drawn:
            legend_weight.Draw()

        canvas.cd(2)
        legend_spectrum = ROOT.TLegend(0.52, 0.62, 0.90, 0.90)
        spectrum_drawn = False
        for histogram, label, color, style in (
            (before, "post-proton/RF input", ROOT.kBlack, 1),
            (proposed_after, "proposed component-cleaned", ROOT.kRed + 1, 2),
            (final_after, "final applied/fallback-cleaned", ROOT.kGreen + 2, 1),
        ):
            if histogram is None:
                continue
            display = clone_root_histogram(
                histogram,
                scope="{}".format(page_prefix),
                role="proposal_final_spectrum_display",
            )
            display.SetLineColor(color)
            display.SetLineStyle(style)
            display.SetTitle("{} before / proposed / final;MM [GeV];yield".format(title_prefix))
            display.Draw("hist" if not spectrum_drawn else "hist same")
            display_hists.append(display)
            legend_spectrum.AddEntry(display, label, "l")
            spectrum_drawn = True
        if spectrum_drawn:
            legend_spectrum.Draw()
        canvas.Print(pdf_name)
        canvas.Close()
        _record_parent_page_if_absent(
            manifest, "{}.proposal_final_transition".format(page_prefix), scope=parent["analysis_scope"]
        )
        return True
    except Exception as exc:
        parent["diagnostic_application_status"] = {
            **dict(parent.get("diagnostic_application_status") or {}),
            "proposal_final_page_status": "unavailable",
            "proposal_final_page_detail": str(exc),
        }
        return False


def _print_parent_summary_page(pdf_name, parents, setting_wide_summary, manifest, title_prefix):
    """Add a compact setting-wide versus canonical-t parent table."""
    try:
        import ROOT

        canvas = ROOT.TCanvas("pion_parent_summary", "Pion parent summary", 1100, 700)
        text = ROOT.TPaveText(0.03, 0.05, 0.97, 0.95, "NDC")
        text.SetFillStyle(0)
        text.SetBorderSize(0)
        text.SetTextAlign(12)
        text.SetTextSize(0.020)
        text.AddText("{} setting-wide diagnostic versus authoritative t-bin pion parents".format(title_prefix))
        text.AddText("scope       A_n      A_delta    A_sidis    S_lambda   chi2/ndf   p       gate/final")
        rows = []
        if isinstance(setting_wide_summary, dict):
            rows.append(("setting", setting_wide_summary))
        rows.extend(("t{}".format(int(parent["t_bin_index"]) + 1), _parent_summary(parent)) for parent in parents)
        for label, row in rows:
            text.AddText("{:<10} {:>8} {:>10} {:>10} {:>10} {:>10} {:>7}  {}/{}".format(
                label,
                str(row.get("A_n")), str(row.get("A_delta")), str(row.get("A_sidis")),
                str(row.get("S_lambda")), str(row.get("chi2_ndf")), str(row.get("p_value")),
                row.get("production_gate_status", "diagnostic"), row.get("final_application_status", row.get("status")),
            ))
        text.Draw()
        canvas.Print(pdf_name)
        canvas.Close()
        _record_parent_page_if_absent(manifest, "pion.t_bin.summary", scope="setting")
        return True
    except Exception:
        return False


def _parent_plot_contract(page_manifest, parents, *, setting_wide_enabled):
    page_ids = [entry.get("page_id") for entry in page_manifest if isinstance(entry, dict)]
    duplicate_ids = sorted({page_id for page_id in page_ids if page_ids.count(page_id) > 1})
    required = []
    if setting_wide_enabled:
        required.append("pion.setting_wide.lambda_comparison")
    for parent in parents:
        prefix = "pion.t_bin.{}".format(int(parent["t_bin_index"]))
        required.extend((
            "{}.overview".format(prefix),
            "{}.application_gate".format(prefix),
            "{}.lambda_comparison".format(prefix),
        ))
    missing = [page_id for page_id in required if page_id not in page_ids]
    by_scope = {}
    for page_id in page_ids:
        if not isinstance(page_id, str):
            continue
        scope = "setting_wide" if page_id.startswith("pion.setting_wide") else ".".join(page_id.split(".")[:3])
        by_scope[scope] = by_scope.get(scope, 0) + 1
    report = {
        "page_counts": by_scope,
        "mandatory_lambda_pages": len([page_id for page_id in page_ids if str(page_id).endswith(".lambda_comparison")]),
        "mandatory_lambda_expected": len(parents) + (1 if setting_wide_enabled else 0),
        "duplicate_page_ids": duplicate_ids,
        "missing_mandatory_page_ids": missing,
    }
    print("[PLOT CONTRACT]")
    for scope, count in sorted(by_scope.items()):
        print("  {} pion pages: {}".format(scope, count))
    print("  mandatory Lambda pages: {} / expected {}".format(
        report["mandatory_lambda_pages"], report["mandatory_lambda_expected"]
    ))
    print("  duplicate page IDs: {}".format(len(duplicate_ids)))
    print("  missing mandatory page IDs: {}".format(missing))
    return report


def render_setting_t_bin_pion_parent_pages(
    pdf_name,
    parents,
    inp_dict,
    *,
    title_prefix,
    page_manifest,
    setting_wide_summary=None,
    setting_wide_enabled=True,
):
    """Render all authoritative parent sections before ``rand_sub`` returns."""
    _print_parent_summary_page(
        pdf_name, parents, setting_wide_summary, page_manifest, title_prefix
    )
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
        final_payload = parent.get("final_diagnostic_application_result")
        proposal_payload = parent.get("proposed_diagnostic_application_result")
        _print_parent_application_gate_page(
            pdf_name, parent, section_title, page_manifest, page_prefix
        )
        _print_parent_proposal_final_pages(
            pdf_name, parent, section_title, page_manifest, page_prefix
        )
        if isinstance(final_payload, dict):
            print_particle_subtraction_component_application_pages(
                pdf_name,
                final_payload,
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
            final_payload,
            title_prefix=section_title,
            cut_window=(float(inp_dict["mm_min"]), float(inp_dict["mm_max"])),
            page_manifest=page_manifest,
            page_id_prefix=page_prefix,
            authoritative=True,
            proposal_payload=proposal_payload,
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
    return _parent_plot_contract(
        page_manifest,
        parents,
        setting_wide_enabled=bool(setting_wide_enabled),
    )


__all__ = (
    "build_setting_t_bin_pion_parents",
    "render_setting_t_bin_pion_parent_pages",
    "validate_frozen_t_bin_pion_parent_collection",
)
