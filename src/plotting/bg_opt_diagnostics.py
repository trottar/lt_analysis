#! /usr/bin/python

import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd


LOWER_IS_BETTER_METRICS = {
    "ratio_fail_count",
    "ratio_mean_dev",
    "ratio_rms",
    "kinematic_score",
    "composite_objective",
}


def _to_bool(value):
    if pd.isna(value):
        return False
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def _to_numeric_columns(df, columns):
    for column in columns:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")
    return df


def _format_bin_label(row):
    try:
        return "{}t x {}phi".format(int(row["requested_num_t_bins"]), int(row["requested_num_phi_bins"]))
    except Exception:
        return "unknown"


def _minmax_lower(series):
    series = pd.to_numeric(series, errors="coerce")
    valid = series.dropna()
    if valid.empty or valid.max() == valid.min():
        return pd.Series(0.0, index=series.index)
    return (series - valid.min()) / (valid.max() - valid.min())


def _minmax_higher(series):
    series = pd.to_numeric(series, errors="coerce")
    valid = series.dropna()
    if valid.empty or valid.max() == valid.min():
        return pd.Series(0.0, index=series.index)
    return (valid.max() - series) / (valid.max() - valid.min())


def _load_bg_opt_csv(csv_path):
    csv_path = Path(csv_path)
    df = pd.read_csv(csv_path)
    df.columns = [str(column).replace("selected_bi\nn_candidate", "selected_bin_candidate") for column in df.columns]

    bool_columns = ["valid", "fallback", "selected_bin_candidate", "selected_for_stage"]
    for column in bool_columns:
        if column not in df.columns:
            df[column] = False
        df["{}_bool".format(column)] = df[column].map(_to_bool)

    numeric_columns = [
        "bg_stat_scale2",
        "requested_num_t_bins",
        "requested_num_phi_bins",
        "actual_num_t_bins",
        "actual_num_phi_bins",
        "ratio_fail_count",
        "ratio_mean_dev",
        "ratio_rms",
        "kinematic_score",
        "valid_ratio_bins",
        "kinematic_points",
        "selection_score",
        "test_model_chi2",
    ]
    df = _to_numeric_columns(df, numeric_columns)
    df["bin_label"] = df.apply(_format_bin_label, axis=1)

    aggregate = df[df["row_kind"] == "bin_candidate"].copy()
    phi_selected = df[df["row_kind"] == "phi_selection"].copy()
    scale = df[
        (df["row_kind"] == "scale_candidate")
        & (df["stage"].isin(["coarse", "refined"]))
    ].copy()

    if not aggregate.empty and "selection_score" in aggregate.columns and aggregate["selection_score"].notna().any():
        aggregate["composite_objective"] = aggregate["selection_score"]
    elif not aggregate.empty:
        aggregate["composite_objective"] = (
            _minmax_lower(aggregate["ratio_fail_count"])
            + _minmax_lower(aggregate["ratio_rms"])
            + _minmax_lower(aggregate["ratio_mean_dev"])
            + _minmax_lower(aggregate["kinematic_score"])
            + _minmax_higher(aggregate["valid_ratio_bins"])
        )
    else:
        aggregate["composite_objective"] = pd.Series(dtype=float)

    return df, aggregate, phi_selected, scale


def _select_aggregate_row(aggregate):
    if aggregate.empty:
        return None

    selected = aggregate[aggregate["selected_bin_candidate_bool"]]
    if not selected.empty:
        return selected.sort_values(
            ["ratio_fail_count", "ratio_mean_dev", "ratio_rms", "kinematic_score"]
        ).iloc[0]

    return aggregate.sort_values(
        ["composite_objective", "ratio_fail_count", "ratio_rms", "ratio_mean_dev", "kinematic_score"]
    ).iloc[0]


def _selected_bin_tuple(selected_row):
    if selected_row is None:
        return None
    return (
        int(selected_row["requested_num_t_bins"]),
        int(selected_row["requested_num_phi_bins"]),
    )


def _format_metric(value, fmt="{:.3f}"):
    try:
        val = float(value)
    except Exception:
        return "nan"
    if not math.isfinite(val):
        return "inf"
    return fmt.format(val)


def _add_text_page(pdf, title, lines):
    fig, ax = plt.subplots(figsize=(8.5, 11.0))
    ax.axis("off")
    ax.set_title(title, fontsize=16, pad=18)
    fig.text(
        0.05,
        0.96,
        "\n".join(lines),
        va="top",
        ha="left",
        family="monospace",
        fontsize=10,
    )
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def _draw_heatmap(ax, frame, metric, title, selected_bin=None, fmt="{:.2f}", cmap="viridis"):
    pivot = frame.pivot_table(
        index="requested_num_t_bins",
        columns="requested_num_phi_bins",
        values=metric,
        aggfunc="mean",
    ).sort_index().sort_index(axis=1)

    if pivot.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(title)
        return None

    values = np.ma.masked_invalid(pivot.values.astype(float))
    im = ax.imshow(values, aspect="auto", cmap=cmap)
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels([str(int(val)) for val in pivot.columns])
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels([str(int(val)) for val in pivot.index])
    ax.set_xlabel("Requested phi bins")
    ax.set_ylabel("Requested t bins")
    ax.set_title(title)

    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            value = pivot.values[i, j]
            if pd.isna(value):
                label = "nan"
            else:
                label = fmt.format(value)
            ax.text(j, i, label, ha="center", va="center", fontsize=8)

    if selected_bin is not None:
        selected_t, selected_phi = selected_bin
        if selected_t in pivot.index and selected_phi in pivot.columns:
            row_idx = list(pivot.index).index(selected_t)
            col_idx = list(pivot.columns).index(selected_phi)
            ax.add_patch(
                Rectangle(
                    (col_idx - 0.5, row_idx - 0.5),
                    1.0,
                    1.0,
                    fill=False,
                    edgecolor="red",
                    linewidth=2.5,
                )
            )

    return im


def _add_heatmap_page(pdf, aggregate, metric, title, selected_bin=None, fmt="{:.2f}", cmap="viridis"):
    fig, ax = plt.subplots(figsize=(7.6, 5.8))
    im = _draw_heatmap(ax, aggregate, metric, title, selected_bin=selected_bin, fmt=fmt, cmap=cmap)
    if im is not None:
        fig.colorbar(im, ax=ax, label=metric)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def _add_top_candidates_page(pdf, aggregate, selected_bin):
    if aggregate.empty:
        return

    top = aggregate.sort_values(
        ["composite_objective", "ratio_fail_count", "ratio_rms", "ratio_mean_dev", "kinematic_score"]
    ).head(12)

    score_label = "selection score" if "selection_score" in aggregate.columns and aggregate["selection_score"].notna().any() else "composite objective"

    lines = [
        "Top aggregate bin candidates by {}".format(score_label),
        "",
        "rank  sel  bins        fail   mean_dev   rms      kin      valid   score",
        "-" * 80,
    ]
    for rank, (_, row) in enumerate(top.iterrows(), start=1):
        marker = "*"
        if selected_bin is not None:
            marker = "*" if (
                int(row["requested_num_t_bins"]) == selected_bin[0]
                and int(row["requested_num_phi_bins"]) == selected_bin[1]
            ) else " "
        lines.append(
                "{:>4}   {}   {:<10} {:>4}   {:>8}   {:>7}  {:>7}  {:>5}   {:>8}".format(
                    rank,
                    marker,
                    row["bin_label"],
                int(row["ratio_fail_count"]),
                _format_metric(row["ratio_mean_dev"], "{:.4f}"),
                _format_metric(row["ratio_rms"], "{:.4f}"),
                _format_metric(row["kinematic_score"], "{:.4f}"),
                int(row["valid_ratio_bins"]),
                    _format_metric(row["composite_objective"], "{:.4f}"),
                )
        )

    _add_text_page(pdf, "Aggregate Candidate Ranking", lines)


def _add_aggregate_tradeoff_page(pdf, aggregate, selected_bin):
    if aggregate.empty:
        return

    fig, ax = plt.subplots(figsize=(8.2, 6.0))
    scatter = ax.scatter(
        aggregate["ratio_fail_count"],
        aggregate["ratio_rms"],
        c=aggregate["kinematic_score"],
        s=40 + 10 * aggregate["valid_ratio_bins"].fillna(0.0),
        cmap="viridis_r",
        alpha=0.85,
    )

    for _, row in aggregate.iterrows():
        is_selected = selected_bin is not None and (
            int(row["requested_num_t_bins"]) == selected_bin[0]
            and int(row["requested_num_phi_bins"]) == selected_bin[1]
        )
        ax.annotate(
            row["bin_label"],
            (row["ratio_fail_count"], row["ratio_rms"]),
            textcoords="offset points",
            xytext=(4, 4),
            fontsize=8,
            fontweight="bold" if is_selected else "normal",
        )
        if is_selected:
            ax.scatter(
                [row["ratio_fail_count"]],
                [row["ratio_rms"]],
                s=220,
                facecolors="none",
                edgecolors="red",
                linewidths=2.0,
            )

    ax.set_xlabel("Ratio fail count")
    ax.set_ylabel("Ratio RMS")
    ax.set_title("Aggregate bin-candidate tradeoff")
    fig.colorbar(scatter, ax=ax, label="kinematic_score")
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def _add_phi_heatmap_page(pdf, phi_setting, phi_rows, selected_bin):
    if phi_rows.empty:
        return

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 8.0))
    axes = np.atleast_1d(axes).ravel()
    metrics = [
        ("bg_stat_scale2", "Selected BG_STAT_SCALE2", "{:.3f}", "magma"),
        ("ratio_fail_count", "Ratio fail count", "{:.0f}", "viridis_r"),
        ("ratio_rms", "Ratio RMS", "{:.3f}", "viridis_r"),
        ("kinematic_score", "Kinematic score", "{:.3f}", "viridis_r"),
    ]

    colorbars = []
    for ax, (metric, title, fmt, cmap) in zip(axes, metrics):
        im = _draw_heatmap(
            ax,
            phi_rows,
            metric,
            "{}: {}".format(phi_setting, title),
            selected_bin=selected_bin,
            fmt=fmt,
            cmap=cmap,
        )
        colorbars.append((ax, im, metric))

    for ax, im, metric in colorbars:
        if im is not None:
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=metric)

    fig.suptitle("Per-phi selected-scale behavior across binning: {}".format(phi_setting), fontsize=15)
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.97])
    pdf.savefig(fig)
    plt.close(fig)


def _plot_metric_lines(ax, frame, metric, title, selected_bin=None, selected_scale=None):
    if frame.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(title)
        return

    stage_styles = {
        "coarse": {"linestyle": "--", "alpha": 0.8},
        "refined": {"linestyle": "-", "alpha": 0.95},
    }
    colors = {"coarse": "#4c78a8", "refined": "#f58518"}

    for stage, group in frame.groupby("stage"):
        style = stage_styles.get(stage, {"linestyle": "-", "alpha": 0.9})
        group = group.sort_values("bg_stat_scale2")
        ax.plot(
            group["bg_stat_scale2"],
            group[metric],
            marker="o",
            linewidth=2.0,
            label=stage.capitalize(),
            color=colors.get(stage, None),
            **style
        )

        chosen = group[group["selected_for_stage_bool"]]
        if not chosen.empty:
            ax.scatter(
                chosen["bg_stat_scale2"],
                chosen[metric],
                s=80,
                marker="*",
                color="red",
                zorder=5,
            )

    if selected_scale is not None and math.isfinite(float(selected_scale)):
        ax.axvline(float(selected_scale), color="red", linewidth=1.5, alpha=0.7)

    ax.set_xlabel("BG_STAT_SCALE2")
    ax.set_ylabel(metric)
    ax.set_title(title)
    ax.grid(alpha=0.25)
    ax.legend()


def _add_selected_bin_scale_page(pdf, scale_rows, phi_setting, selected_bin, selected_scale):
    if selected_bin is None:
        return

    selected_t, selected_phi = selected_bin
    frame = scale_rows[
        (scale_rows["phi_setting"] == phi_setting)
        & (scale_rows["requested_num_t_bins"] == selected_t)
        & (scale_rows["requested_num_phi_bins"] == selected_phi)
    ].copy()
    if frame.empty:
        return

    frame = frame.drop_duplicates(["stage", "bg_stat_scale2", "ratio_fail_count", "ratio_mean_dev", "ratio_rms", "kinematic_score"])
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 8.0))
    axes = np.atleast_1d(axes).ravel()
    metrics = [
        ("ratio_fail_count", "Failure count vs BG_STAT_SCALE2"),
        ("ratio_mean_dev", "Mean ratio deviation vs BG_STAT_SCALE2"),
        ("ratio_rms", "Ratio RMS vs BG_STAT_SCALE2"),
        ("kinematic_score", "Kinematic score vs BG_STAT_SCALE2"),
    ]

    for ax, (metric, title) in zip(axes, metrics):
        _plot_metric_lines(ax, frame, metric, title, selected_bin=selected_bin, selected_scale=selected_scale)

    fig.suptitle(
        "Selected-binning scale scan: {} ({})".format(phi_setting, frame["bin_label"].iloc[0]),
        fontsize=15,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.97])
    pdf.savefig(fig)
    plt.close(fig)


def _plot_full_phase_space_metric(ax, frame, metric, title, selected_bin):
    if frame.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(title)
        return

    frame = frame.sort_values(["requested_num_t_bins", "requested_num_phi_bins", "bg_stat_scale2"])
    for bin_label, group in frame.groupby("bin_label"):
        group = group.sort_values("bg_stat_scale2")
        try:
            group_t = int(group["requested_num_t_bins"].iloc[0])
            group_phi = int(group["requested_num_phi_bins"].iloc[0])
        except Exception:
            group_t, group_phi = None, None

        is_selected = selected_bin is not None and group_t == selected_bin[0] and group_phi == selected_bin[1]
        ax.plot(
            group["bg_stat_scale2"],
            group[metric],
            marker="o",
            linewidth=2.4 if is_selected else 1.1,
            alpha=0.95 if is_selected else 0.45,
            label="{}{}".format(bin_label, " *" if is_selected else ""),
        )

        chosen = group[group["selected_for_stage_bool"]]
        if not chosen.empty:
            ax.scatter(
                chosen["bg_stat_scale2"],
                chosen[metric],
                s=70,
                marker="*",
                color="red" if is_selected else "black",
                zorder=5,
            )

    ax.set_xlabel("BG_STAT_SCALE2")
    ax.set_ylabel(metric)
    ax.set_title(title)
    ax.grid(alpha=0.25)


def _add_full_phase_space_page(pdf, scale_rows, phi_setting, stage, selected_bin):
    frame = scale_rows[
        (scale_rows["phi_setting"] == phi_setting)
        & (scale_rows["stage"] == stage)
    ].copy()
    if frame.empty:
        return

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.5))
    axes = np.atleast_1d(axes).ravel()
    metrics = [
        ("ratio_fail_count", "Failure count"),
        ("ratio_mean_dev", "Mean ratio deviation"),
        ("ratio_rms", "Ratio RMS"),
        ("kinematic_score", "Kinematic score"),
    ]
    for ax, (metric, title) in zip(axes, metrics):
        _plot_full_phase_space_metric(
            ax,
            frame,
            metric,
            "{} {}: {}".format(phi_setting, stage.capitalize(), title),
            selected_bin,
        )

    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=min(3, len(labels)), fontsize=8, frameon=False)

    fig.suptitle(
        "Full phase-space scale evolution: {} ({})".format(phi_setting, stage.capitalize()),
        fontsize=15,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.93])
    pdf.savefig(fig)
    plt.close(fig)


def _build_cover_lines(csv_path, df, aggregate, phi_selected, selected_row, selected_phi_rows, scale_rows):
    if df.empty:
        return ["No optimizer data found in {}".format(csv_path.name)]

    mode = str(df["mode"].dropna().iloc[0]) if "mode" in df and not df["mode"].dropna().empty else "unknown"
    selection_mode = str(df["selection_mode"].dropna().iloc[0]) if "selection_mode" in df and not df["selection_mode"].dropna().empty else "unknown"
    epsset = str(df["epsset"].dropna().iloc[0]) if "epsset" in df and not df["epsset"].dropna().empty else "unknown"
    particle = str(df["particle"].dropna().iloc[0]) if "particle" in df and not df["particle"].dropna().empty else "unknown"
    q2 = str(df["q2"].dropna().iloc[0]) if "q2" in df and not df["q2"].dropna().empty else "unknown"
    w = str(df["w"].dropna().iloc[0]) if "w" in df and not df["w"].dropna().empty else "unknown"
    outfilename = str(df["outfilename"].dropna().iloc[0]) if "outfilename" in df and not df["outfilename"].dropna().empty else csv_path.stem

    lines = [
        "Optimizer diagnostics for {}".format(csv_path.name),
        "",
        "particle:   {}".format(particle),
        "mode:       {}".format(mode),
        "selection:  {}".format(selection_mode),
        "epsset:     {}".format(epsset),
        "Q2 / W:     {} / {}".format(q2, w),
        "outfile:    {}".format(outfilename),
        "",
        "aggregate bin candidates scanned: {}".format(len(aggregate)),
        "phi-selection rows:               {}".format(len(phi_selected)),
        "scale-scan rows:                  {}".format(len(scale_rows)),
        "",
    ]

    if selected_row is not None:
        lines.extend([
            "Selected shared binning:",
            "  requested bins: {}t x {}phi".format(
                int(selected_row["requested_num_t_bins"]),
                int(selected_row["requested_num_phi_bins"]),
            ),
            "  aggregate fail count: {}".format(int(selected_row["ratio_fail_count"])),
            "  aggregate mean_dev:   {}".format(_format_metric(selected_row["ratio_mean_dev"], "{:.4f}")),
            "  aggregate rms:        {}".format(_format_metric(selected_row["ratio_rms"], "{:.4f}")),
            "  aggregate kin score:  {}".format(_format_metric(selected_row["kinematic_score"], "{:.4f}")),
            "  valid ratio bins:     {}".format(int(selected_row["valid_ratio_bins"])),
            "  selection score:      {}".format(_format_metric(selected_row.get("composite_objective", float("nan")), "{:.4f}")),
            "",
            "Selected BG_STAT_SCALE2 by phi:",
        ])
        for _, row in selected_phi_rows.sort_values("phi_setting").iterrows():
            lines.append(
                "  {:<8} scale={} fail={} mean_dev={} rms={} kin={} score={}".format(
                    row["phi_setting"],
                    _format_metric(row["bg_stat_scale2"], "{:.3f}"),
                    int(row["ratio_fail_count"]),
                    _format_metric(row["ratio_mean_dev"], "{:.4f}"),
                    _format_metric(row["ratio_rms"], "{:.4f}"),
                    _format_metric(row["kinematic_score"], "{:.4f}"),
                    _format_metric(row.get("selection_score", float("nan")), "{:.4f}"),
                )
            )

    return lines


def plot_bg_optimization_diagnostics(csv_path, pdf_path=None):
    csv_path = Path(csv_path)
    if pdf_path is None:
        pdf_path = csv_path.with_suffix(".pdf")
    pdf_path = Path(pdf_path)

    df, aggregate, phi_selected, scale_rows = _load_bg_opt_csv(csv_path)
    selected_row = _select_aggregate_row(aggregate)
    selected_bin = _selected_bin_tuple(selected_row)

    if selected_bin is not None:
        selected_phi_rows = phi_selected[
            (phi_selected["requested_num_t_bins"] == selected_bin[0])
            & (phi_selected["requested_num_phi_bins"] == selected_bin[1])
        ].copy()
    else:
        selected_phi_rows = phi_selected.copy()

    selected_phi_rows = selected_phi_rows.sort_values(["phi_setting", "ratio_fail_count", "ratio_rms"])

    with PdfPages(pdf_path) as pdf:
        _add_text_page(
            pdf,
            "BG Optimization Diagnostics",
            _build_cover_lines(csv_path, df, aggregate, phi_selected, selected_row, selected_phi_rows, scale_rows),
        )

        _add_top_candidates_page(pdf, aggregate, selected_bin)

        composite_title = "Aggregate bin scan: weighted selection score"
        if "selection_score" not in aggregate.columns or not aggregate["selection_score"].notna().any():
            composite_title = "Aggregate bin scan: composite objective"

        heatmap_specs = [
            ("ratio_fail_count", "Aggregate bin scan: ratio fail count", "{:.0f}", "viridis_r"),
            ("ratio_mean_dev", "Aggregate bin scan: mean ratio deviation", "{:.3f}", "viridis_r"),
            ("ratio_rms", "Aggregate bin scan: ratio RMS", "{:.3f}", "viridis_r"),
            ("kinematic_score", "Aggregate bin scan: kinematic score", "{:.3f}", "viridis_r"),
            ("valid_ratio_bins", "Aggregate bin scan: valid ratio bins", "{:.0f}", "viridis"),
            ("composite_objective", composite_title, "{:.3f}", "viridis_r"),
        ]
        for metric, title, fmt, cmap in heatmap_specs:
            _add_heatmap_page(
                pdf,
                aggregate,
                metric=metric,
                title=title,
                selected_bin=selected_bin,
                fmt=fmt,
                cmap=cmap,
            )

        _add_aggregate_tradeoff_page(pdf, aggregate, selected_bin)

        for phi_setting, phi_rows in phi_selected.groupby("phi_setting"):
            _add_phi_heatmap_page(pdf, phi_setting, phi_rows.copy(), selected_bin)

            selected_scale = float("nan")
            chosen = selected_phi_rows[selected_phi_rows["phi_setting"] == phi_setting]
            if not chosen.empty:
                selected_scale = float(chosen.iloc[0]["bg_stat_scale2"])

            _add_selected_bin_scale_page(pdf, scale_rows, phi_setting, selected_bin, selected_scale)
            for stage in ("coarse", "refined"):
                _add_full_phase_space_page(pdf, scale_rows, phi_setting, stage, selected_bin)

    return str(pdf_path)
