#! /usr/bin/python

from __future__ import annotations

import logging
import os

import matplotlib

matplotlib.use("Agg")
logging.getLogger("matplotlib").setLevel(logging.WARNING)
logging.getLogger("matplotlib.font_manager").setLevel(logging.WARNING)
logging.getLogger("matplotlib.backends.backend_pdf").setLevel(logging.WARNING)
import matplotlib.pyplot as plt


COMPONENT_STYLES = {
    "pi_n": {"color": "#1f77b4", "label": "pi-n"},
    "pi_delta": {"color": "#d62728", "label": "pi-delta"},
    "pi_sidis": {"color": "#2ca02c", "label": "pi-SIDIS"},
}
KAON_AUX_STYLES = {
    "k_lambda_signal": {"color": "#6a3d9a", "label": "K-Lambda"},
    "k_sigma0_signal": {"color": "#17becf", "label": "K-Sigma0"},
}


def _hist_to_stairs(hist):
    if hist is None:
        return None, None
    nbins = hist.GetNbinsX()
    edges = [hist.GetXaxis().GetBinLowEdge(i) for i in range(1, nbins + 2)]
    values = [hist.GetBinContent(i) for i in range(1, nbins + 1)]
    return edges, values


def _format_diag_line(component_name, diagnostics):
    if not diagnostics:
        return "{}: no diagnostics".format(component_name)
    return (
        "{}: seen={} passed={} mm_passed={} norm={} fallback={}".format(
            component_name,
            diagnostics.get("n_events_seen", 0),
            diagnostics.get("n_events_passed", 0),
            diagnostics.get("n_events_passed_mm_window", 0),
            diagnostics.get("normalized", False),
            diagnostics.get("fallback_used", False),
        )
    )


def plot_pion_component_background_payload(
    payload,
    phi_setting,
    inpDict,
    pdf_path,
    kaon_signal_payload=None,
    kaon_sigma0_payload=None,
):
    payload = payload or {}
    component_map = payload.get("components") or {}
    kaon_signal_payload = kaon_signal_payload or {}
    kaon_sigma0_payload = kaon_sigma0_payload or {}
    if not component_map and not kaon_signal_payload and not kaon_sigma0_payload:
        return None

    phi_setting = phi_setting or "Unknown"
    particle_type = inpDict.get("ParticleType", "particle")
    q2 = inpDict.get("Q2", "")
    w = inpDict.get("W", "")
    epsset = inpDict.get("EPSSET", "")
    mm_min = float(inpDict.get("mm_min", 0.0))
    mm_max = float(inpDict.get("mm_max", 0.0))

    fig, axes = plt.subplots(2, 1, figsize=(10.5, 8.5), constrained_layout=True)
    full_ax, cut_ax = axes

    for component_name in ("pi_n", "pi_delta", "pi_sidis"):
        component_payload = component_map.get(component_name) or {}
        diagnostics = component_payload.get("diagnostics", {})
        full_hist = component_payload.get("setting_shape_full")
        cut_hist = component_payload.get("setting_shape")
        style = COMPONENT_STYLES.get(
            component_name,
            {"color": "black", "label": component_name},
        )
        full_edges, full_values = _hist_to_stairs(full_hist)
        cut_edges, cut_values = _hist_to_stairs(cut_hist)
        legend_label = "{} | pass={} | fallback={}".format(
            style["label"],
            diagnostics.get("n_events_passed", 0),
            diagnostics.get("fallback_used", False),
        )

        if full_edges is not None:
            full_ax.stairs(
                full_values,
                full_edges,
                label=legend_label,
                color=style["color"],
                linewidth=1.8,
            )
        if cut_edges is not None:
            cut_ax.stairs(
                cut_values,
                cut_edges,
                label=style["label"],
                color=style["color"],
                linewidth=1.8,
            )

    aux_payloads = [
        ("k_lambda_signal", kaon_signal_payload),
        ("k_sigma0_signal", kaon_sigma0_payload),
    ]
    aux_present = []
    for aux_name, aux_payload in aux_payloads:
        aux_diag = aux_payload.get("diagnostics", {})
        aux_full_hist = aux_payload.get("setting_shape_full")
        aux_cut_hist = aux_payload.get("setting_shape")
        aux_full_edges, aux_full_values = _hist_to_stairs(aux_full_hist)
        aux_cut_edges, aux_cut_values = _hist_to_stairs(aux_cut_hist)
        aux_style = KAON_AUX_STYLES.get(aux_name, {"color": "black", "label": aux_name})

        if aux_full_edges is not None:
            full_ax.stairs(
                aux_full_values,
                aux_full_edges,
                label="{} | pass={} | fallback={}".format(
                    aux_style["label"],
                    aux_diag.get("n_events_passed", 0),
                    aux_diag.get("fallback_used", False),
                ),
                color=aux_style["color"],
                linewidth=1.8,
            )
            aux_present.append(aux_name)
        if aux_cut_edges is not None:
            cut_ax.stairs(
                aux_cut_values,
                aux_cut_edges,
                label=aux_style["label"],
                color=aux_style["color"],
                linewidth=1.8,
            )

    full_ax.axvline(mm_min, color="0.4", linestyle="--", linewidth=1.0)
    full_ax.axvline(mm_max, color="0.4", linestyle="--", linewidth=1.0)
    title_label = "pion-background SIMC components"
    if aux_present:
        title_label = "pion-background + kaon auxiliary SIMC components"
    full_ax.set_title(
        "{} {} {} | Q{} W{} {}".format(
            phi_setting,
            particle_type,
            title_label,
            q2,
            w,
            epsset,
        )
    )
    full_ax.set_xlabel("Missing Mass")
    full_ax.set_ylabel("Unit-normalized density")
    full_ax.grid(alpha=0.25)
    full_ax.legend(loc="upper right", fontsize=8)

    cut_ax.set_title("Same templates in the analysis MM-cut window")
    cut_ax.set_xlabel("Missing Mass")
    cut_ax.set_ylabel("Unit-normalized density")
    cut_ax.grid(alpha=0.25)
    cut_ax.legend(loc="upper right", fontsize=8)

    diagnostics = payload.get("diagnostics") or {}
    diag_lines = [
        "Mode: {}".format(payload.get("mode")),
        "Tree: {}".format(payload.get("tree_name")),
    ]
    for component_name in ("pi_n", "pi_delta", "pi_sidis"):
        diag_lines.append(
            _format_diag_line(component_name, diagnostics.get(component_name))
        )
    for aux_name, aux_payload in aux_payloads:
        aux_full_hist = aux_payload.get("setting_shape_full")
        aux_cut_hist = aux_payload.get("setting_shape")
        if aux_full_hist is None and aux_cut_hist is None:
            continue
        diag_lines.append(
            _format_diag_line(aux_name, aux_payload.get("diagnostics", {}))
        )
    fig.text(
        0.015,
        0.01,
        "\n".join(diag_lines),
        fontsize=8,
        family="monospace",
        va="bottom",
    )

    output_dir = os.path.dirname(pdf_path)
    if output_dir and not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    fig.savefig(pdf_path)
    plt.close(fig)
    return pdf_path


def plot_pion_component_backgrounds(hist, inpDict, pdf_path):
    return plot_pion_component_background_payload(
        hist.get("_simc_pion_component_payload"),
        hist.get("phi_setting", "Unknown"),
        inpDict,
        pdf_path,
        kaon_signal_payload=hist.get("_simc_kaon_signal_shape_payload"),
        kaon_sigma0_payload=hist.get("_simc_kaon_sigma0_shape_payload"),
    )
