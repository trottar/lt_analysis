#! /usr/bin/python

import os

import numpy as np

from ltsep import Root

lt = Root(os.path.realpath(__file__), "Plot_LTSep")
OUTPATH = lt.OUTPATH


def _hist_to_arrays(hist):
    n_bins = hist.GetNbinsX()
    values = np.array([hist.GetBinContent(i) for i in range(1, n_bins + 1)], dtype=np.float64)
    errors = np.array([hist.GetBinError(i) for i in range(1, n_bins + 1)], dtype=np.float64)
    edges = np.array([hist.GetBinLowEdge(i) for i in range(1, n_bins + 2)], dtype=np.float64)
    return values, errors, edges


def _hist2d_to_arrays(hist):
    n_x = hist.GetNbinsX()
    n_y = hist.GetNbinsY()
    values = np.zeros((n_x, n_y), dtype=np.float64)
    for ix in range(1, n_x + 1):
        for iy in range(1, n_y + 1):
            values[ix - 1, iy - 1] = hist.GetBinContent(ix, iy)
    x_edges = np.array([hist.GetXaxis().GetBinLowEdge(i) for i in range(1, n_x + 2)], dtype=np.float64)
    y_edges = np.array([hist.GetYaxis().GetBinLowEdge(i) for i in range(1, n_y + 2)], dtype=np.float64)
    return values, x_edges, y_edges


def _matrix_to_arrays(hist_matrix):
    n_t = len(hist_matrix)
    n_phi = len(hist_matrix[0]) if n_t > 0 else 0
    sample_hist = None

    for row in hist_matrix:
        for hist in row:
            if hist is not None:
                sample_hist = hist
                break
        if sample_hist is not None:
            break

    if sample_hist is None:
        return None, None, None

    sample_values, sample_errors, edges = _hist_to_arrays(sample_hist)
    n_bins = len(sample_values)
    values = np.zeros((n_t, n_phi, n_bins), dtype=np.float64)
    errors = np.zeros((n_t, n_phi, n_bins), dtype=np.float64)

    for j in range(n_t):
        for k in range(n_phi):
            hist = hist_matrix[j][k]
            if hist is None:
                continue
            hist_values, hist_errors, _ = _hist_to_arrays(hist)
            values[j, k, :] = hist_values
            errors[j, k, :] = hist_errors

    return values, errors, edges


def _matrix2d_to_arrays(hist_matrix):
    n_t = len(hist_matrix)
    n_phi = len(hist_matrix[0]) if n_t > 0 else 0
    sample_hist = None

    for row in hist_matrix:
        for hist in row:
            if hist is not None:
                sample_hist = hist
                break
        if sample_hist is not None:
            break

    if sample_hist is None:
        return None, None, None

    sample_values, x_edges, y_edges = _hist2d_to_arrays(sample_hist)
    n_x, n_y = sample_values.shape
    values = np.zeros((n_t, n_phi, n_x, n_y), dtype=np.float64)

    for j in range(n_t):
        for k in range(n_phi):
            hist = hist_matrix[j][k]
            if hist is None:
                continue
            hist_values, _, _ = _hist2d_to_arrays(hist)
            values[j, k, :, :] = hist_values

    return values, x_edges, y_edges


def _cleanup_support(histlist):
    for hist in histlist:
        hist.pop("_xsect_support_data", None)
        hist.pop("_xsect_support_simc", None)


def _load_support_from_saved_histograms(hist, hist_prefix):
    t_bins = hist.get("t_bins", [])
    phi_bins = hist.get("phi_bins", [])
    n_t = max(len(t_bins) - 1, 0)
    n_phi = max(len(phi_bins) - 1, 0)
    if n_t == 0 or n_phi == 0:
        return None

    key_templates = {
        "mm": "H_MM_{}_{}_{}",
        "Q2": "H_Q2_{}_{}_{}",
        "W": "H_W_{}_{}_{}",
        "theta_cm": "H_theta_cm_{}_{}_{}",
        "t_vs_tmin": "H_t_vs_tmin_{}_{}_{}",
    }
    support_dict = {
        key: [[None for _ in range(n_phi)] for _ in range(n_t)]
        for key in key_templates
    }

    for support_key, template in key_templates.items():
        for j in range(n_t):
            for k in range(n_phi):
                support_hist = hist.get(template.format(hist_prefix, j, k))
                if support_hist is None:
                    return None
                support_dict[support_key][j][k] = support_hist

    return support_dict


def write_xsect_support(histlist, inpDict, output_file_lst=None):
    particle_type = inpDict["ParticleType"]
    q2 = inpDict["Q2"]
    w = inpDict["W"]
    eps_tag = "{}e".format(inpDict["EPSSET"])
    support_path = os.path.join(
        OUTPATH,
        "{}_xsect_support_Q{}W{}_{}.npz".format(particle_type, q2, w, eps_tag),
    )
    if os.path.exists(support_path):
        os.remove(support_path)

    payload = {}

    try:
        if not histlist:
            return None

        payload["t_bins"] = np.asarray(histlist[0]["t_bins"], dtype=np.float64)
        payload["phi_bins"] = np.asarray(histlist[0]["phi_bins"], dtype=np.float64)
        payload["settings"] = np.asarray([hist["phi_setting"] for hist in histlist], dtype="<U16")

        variable_map = (
            ("mm", "mm"),
            ("Q2", "q2"),
            ("W", "w"),
            ("theta_cm", "sin_theta_cm"),
        )
        map2d_variable_map = (
            ("t_vs_tmin", "t_vs_tmin"),
        )
        simc_only_variable_map = (
            ("theta_cm_true", "sin_theta_cm_true"),
        )

        for hist in histlist:
            setting_key = hist["phi_setting"].lower()
            data_support = hist.get("_xsect_support_data")
            simc_support = hist.get("_xsect_support_simc")

            if data_support is None:
                data_support = _load_support_from_saved_histograms(hist, "DATA")
                if data_support is not None:
                    hist["_xsect_support_data"] = data_support

            if data_support is None or simc_support is None:
                print("WARNING: Missing xsect support histograms for {} {}".format(inpDict["EPSSET"], hist["phi_setting"]))
                return None

            for support_key, file_key in variable_map:
                data_values, data_errors, edges = _matrix_to_arrays(data_support[support_key])
                simc_values, simc_errors, simc_edges = _matrix_to_arrays(simc_support[support_key])

                if data_values is None or simc_values is None:
                    print("WARNING: Incomplete xsect support histograms for {} {} {}".format(inpDict["EPSSET"], hist["phi_setting"], support_key))
                    return None

                payload["{}_edges".format(file_key)] = edges
                payload["data_{}_{}_values".format(file_key, setting_key)] = data_values
                payload["data_{}_{}_errors".format(file_key, setting_key)] = data_errors
                payload["simc_{}_{}_values".format(file_key, setting_key)] = simc_values
                payload["simc_{}_{}_errors".format(file_key, setting_key)] = simc_errors

                if not np.allclose(edges, simc_edges):
                    print("WARNING: Support histogram edges do not match for {} {}".format(hist["phi_setting"], support_key))
                    return None

            for support_key, file_key in map2d_variable_map:
                data_values, data_x_edges, data_y_edges = _matrix2d_to_arrays(data_support[support_key])
                simc_values, simc_x_edges, simc_y_edges = _matrix2d_to_arrays(simc_support[support_key])

                if data_values is None or simc_values is None:
                    print("WARNING: Incomplete 2D xsect support histograms for {} {} {}".format(inpDict["EPSSET"], hist["phi_setting"], support_key))
                    return None

                payload["{}_x_edges".format(file_key)] = data_x_edges
                payload["{}_y_edges".format(file_key)] = data_y_edges
                payload["data_{}_{}_values".format(file_key, setting_key)] = data_values
                payload["simc_{}_{}_values".format(file_key, setting_key)] = simc_values

                if not np.allclose(data_x_edges, simc_x_edges) or not np.allclose(data_y_edges, simc_y_edges):
                    print("WARNING: 2D support histogram edges do not match for {} {}".format(hist["phi_setting"], support_key))
                    return None

            for support_key, file_key in simc_only_variable_map:
                simc_values, simc_errors, simc_edges = _matrix_to_arrays(simc_support[support_key])

                if simc_values is None:
                    print("WARNING: Missing SIMC-only xsect support histograms for {} {} {}".format(inpDict["EPSSET"], hist["phi_setting"], support_key))
                    return None

                payload["{}_edges".format(file_key)] = simc_edges
                payload["simc_{}_{}_values".format(file_key, setting_key)] = simc_values
                payload["simc_{}_{}_errors".format(file_key, setting_key)] = simc_errors

        np.savez_compressed(support_path, **payload)
        if output_file_lst is not None:
            output_file_lst.append(support_path)
        return support_path
    finally:
        _cleanup_support(histlist)
