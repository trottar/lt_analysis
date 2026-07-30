"""Regression checks for pre-particle-subtraction binning inputs."""

from __future__ import annotations

import importlib.util
import math
import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
FIND_BINS_PATH = REPO_ROOT / "src" / "binning" / "find_bins.py"


class _FakeHistogram:
    def __init__(self, values):
        self.values = list(values)

    def Clone(self):
        return _FakeHistogram(self.values)

    def Scale(self, scale):
        self.values = [float(value) * float(scale) for value in self.values]


def _load_find_bins():
    fake_ltsep = types.ModuleType("ltsep")
    fake_ltsep.Misc = types.SimpleNamespace()
    fake_ltsep.Root = lambda *_args, **_kwargs: types.SimpleNamespace(LTANAPATH="")
    fake_utility = types.ModuleType("utility")
    fake_utility.flatten_hist = lambda histogram: list(histogram.values)
    fake_background_config = types.ModuleType("background_config")
    for name, value in {
        "MIN_PHI_BINS": 1,
        "MIN_T_BINS": 1,
        "PHI_BIN_MAX_DEG": 180.0,
        "PHI_BIN_MIN_DEG": -180.0,
        "PHI_BIN_MIN_EVENTS": 1,
        "T_BIN_ADJUST_MAX_ITERATIONS": 1,
        "T_BIN_ADJUST_TOLERANCE": 0.01,
        "T_BIN_EDGE_BIAS": 1.0,
        "T_BIN_MIN_EVENTS": 1,
    }.items():
        setattr(fake_background_config, name, value)
    fake_background_config.get_proton_contamination_cleaning_config = lambda **_kwargs: {
        "strict_mode": True,
        "t_binning": {
            "metadata_schema_version": 1,
            "require_metadata_sidecar": True,
            "allow_legacy_interval_without_metadata": False,
            "allow_authoritative_interval_file": True,
            "edge_tolerance": 1.0e-9,
            "canonical_bin_support_metric": "raw_event_count",
            "canonical_bin_support_threshold": 1,
            "allowed_support_metrics": (
                "raw_event_count", "absolute_weighted_support", "positive_weighted_support",
            ),
            "allowed_algorithm_identifiers": ("find_bins_adjust_t_bins",),
            "allowed_algorithm_versions": (1,),
        },
    }

    with mock.patch.dict(
        sys.modules,
        {
            "ltsep": fake_ltsep,
            "utility": fake_utility,
            "background_config": fake_background_config,
        },
    ):
        spec = importlib.util.spec_from_file_location("find_bins_test_module", FIND_BINS_PATH)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    return module


class PreParticleSubtractionBinningTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.find_bins = _load_find_bins()

    def test_pre_particle_subtraction_histograms_are_preferred(self):
        hist = {
            "phi_setting": "Center",
            "normfac_data": 1.0,
            "H_t_DATA": _FakeHistogram((10.0, 20.0, 30.0)),
            "H_ph_q_DATA": _FakeHistogram((1.0, 1.1, 1.2)),
            self.find_bins.PRE_PARTICLE_SUBTRACTION_T_HIST_KEY: _FakeHistogram((0.1, 0.2, 0.3)),
            self.find_bins.PRE_PARTICLE_SUBTRACTION_PHI_HIST_KEY: _FakeHistogram((0.1, 0.2, 0.3)),
        }

        t_values, phi_values = self.find_bins._collect_bin_samples(
            [hist],
            {"tmin": 0.0, "tmax": 0.4, "normfac_data": 1.0},
        )

        self.assertEqual(t_values.tolist(), [0.0, 0.2, 0.4])
        self.assertAlmostEqual(phi_values[1], 0.2 * (180.0 / math.pi))

    def test_final_histograms_remain_the_legacy_fallback(self):
        final_t = _FakeHistogram((0.1, 0.2, 0.3))
        final_phi = _FakeHistogram((0.1, 0.2, 0.3))
        selected_t, selected_phi = self.find_bins._get_binning_histograms(
            {"H_t_DATA": final_t, "H_ph_q_DATA": final_phi}
        )

        self.assertIs(selected_t, final_t)
        self.assertIs(selected_phi, final_phi)

    def _payload(self, weights=(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)):
        records = []
        for index, (t_value, weight) in enumerate(zip((0.0, 0.2, 0.4, 0.6, 0.8, 1.0), weights)):
            records.append({
                "source_label": "prompt",
                "entry_index": index,
                "adj_t": t_value,
                "phi_value": -1.0 + index * 0.4,
                "physical_coefficient": weight,
            })
        return [{"records": records, "source_stats": {"prompt": {"entries_selected": len(records)}}}]

    @staticmethod
    def _inp():
        return {
            "ParticleType": "kaon", "Q2": "3p0", "W": "3p14", "EPSSET": "low",
            "NumtBins": 2, "NumPhiBins": 2, "tmin": 0.0, "tmax": 1.0,
        }

    def test_metadata_backed_interval_round_trip_is_authoritative(self):
        with tempfile.TemporaryDirectory() as directory:
            self.find_bins.LTANAPATH = directory
            (Path(directory) / "src" / "kaon").mkdir(parents=True)
            inp = self._inp()
            computed = self.find_bins.resolve_canonical_analysis_bins_pre_subtraction(
                self._payload(), inp, quiet=True
            )
            self.assertEqual(computed["source"], "computed_from_pre_particle_subtraction_histograms")
            sidecar = Path(computed["metadata"]["metadata_file"])
            self.assertTrue(sidecar.exists())
            with sidecar.open(encoding="utf-8") as handle:
                metadata = json.load(handle)
            self.assertEqual(metadata["schema_version"], 1)
            self.assertEqual(metadata["t_edges"], computed["t_bins"].tolist())

            reused = self.find_bins.resolve_canonical_analysis_bins_pre_subtraction(
                self._payload(), self._inp(), quiet=True
            )
            self.assertEqual(reused["source"], "validated_authoritative_interval_file")
            self.assertEqual(reused["metadata"]["validation_rejection_reasons"], [])
            high_inp = self._inp()
            high_inp["EPSSET"] = "high"
            high = self.find_bins.resolve_canonical_analysis_bins_pre_subtraction(
                self._payload(), high_inp, quiet=True
            )
            self.assertEqual(high["t_bins"].tolist(), reused["t_bins"].tolist())
            self.assertEqual(high["phi_bins"].tolist(), reused["phi_bins"].tolist())
            self.assertEqual(high["metadata"]["consumer_epsilon"], "high")
            self.assertEqual(high["metadata"]["shared_from_epsilon"], "low")

    def test_phi_sidecar_mismatch_recomputes_low_pair_and_blocks_high(self):
        with tempfile.TemporaryDirectory() as directory:
            self.find_bins.LTANAPATH = directory
            (Path(directory) / "src" / "kaon").mkdir(parents=True)
            first = self.find_bins.resolve_canonical_analysis_bins_pre_subtraction(
                self._payload(), self._inp(), quiet=True
            )
            phi_metadata_path = Path(first["metadata"]["phi_metadata_file"])
            with phi_metadata_path.open(encoding="utf-8") as handle:
                phi_metadata = json.load(handle)
            phi_metadata["phi_edges"][1] += 0.5
            with phi_metadata_path.open("w", encoding="utf-8") as handle:
                json.dump(phi_metadata, handle)
            low = self.find_bins.resolve_canonical_analysis_bins_pre_subtraction(
                self._payload(), self._inp(), quiet=True
            )
            self.assertEqual(low["source"], "computed_from_pre_particle_subtraction_histograms")
            self.assertIn("phi_text_metadata_edges_mismatch", low["metadata"]["validation_rejection_reasons"])

            high_inp = self._inp()
            high_inp["EPSSET"] = "high"
            with phi_metadata_path.open("w", encoding="utf-8") as handle:
                json.dump({"schema_version": 999}, handle)
            with self.assertRaisesRegex(RuntimeError, "t/phi interval pair"):
                self.find_bins.resolve_canonical_analysis_bins_pre_subtraction(
                    self._payload(), high_inp, quiet=True
                )

    def test_legacy_writer_cannot_overwrite_canonical_sidecars(self):
        with tempfile.TemporaryDirectory() as directory:
            self.find_bins.LTANAPATH = directory
            (Path(directory) / "src" / "kaon").mkdir(parents=True)
            inp = self._inp()
            self.find_bins.resolve_canonical_analysis_bins_pre_subtraction(
                self._payload(), inp, quiet=True
            )
            with self.assertRaisesRegex(RuntimeError, "canonical_interval_overwrite_refused"):
                self.find_bins.write_bin_interval_files(
                    {**self._inp()}, np.array([0.0, 0.5, 1.0]), np.array([-180.0, 0.0, 180.0])
                )

    def test_missing_sidecar_is_rejected_then_prepass_replaces_it(self):
        with tempfile.TemporaryDirectory() as directory:
            self.find_bins.LTANAPATH = directory
            (Path(directory) / "src" / "kaon").mkdir(parents=True)
            first = self.find_bins.resolve_canonical_analysis_bins_pre_subtraction(
                self._payload(), self._inp(), quiet=True
            )
            Path(first["metadata"]["metadata_file"]).unlink()
            fallback = self.find_bins.resolve_canonical_analysis_bins_pre_subtraction(
                self._payload(), self._inp(), quiet=True
            )
            self.assertEqual(fallback["source"], "computed_from_pre_particle_subtraction_histograms")
            self.assertEqual(fallback["metadata"]["validation_status"], "computed_pre_particle_subtraction")
            self.assertIn("metadata_sidecar_missing", fallback["metadata"]["validation_rejection_reasons"])

    def test_text_metadata_edge_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            self.find_bins.LTANAPATH = directory
            (Path(directory) / "src" / "kaon").mkdir(parents=True)
            first = self.find_bins.resolve_canonical_analysis_bins_pre_subtraction(
                self._payload(), self._inp(), quiet=True
            )
            metadata_path = Path(first["metadata"]["metadata_file"])
            with metadata_path.open(encoding="utf-8") as handle:
                metadata = json.load(handle)
            metadata["t_edges"][1] += 0.01
            with metadata_path.open("w", encoding="utf-8") as handle:
                json.dump(metadata, handle)
            validated = self.find_bins._validate_authoritative_t_interval(self._inp(), "low")
            self.assertFalse(validated["valid"])
            self.assertIn("text_metadata_edges_mismatch", validated["validation_rejection_reasons"])

    def test_metadata_identity_and_configuration_mismatches_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            self.find_bins.LTANAPATH = directory
            (Path(directory) / "src" / "kaon").mkdir(parents=True)
            first = self.find_bins.resolve_canonical_analysis_bins_pre_subtraction(
                self._payload(), self._inp(), quiet=True
            )
            metadata_path = Path(first["metadata"]["metadata_file"])
            with metadata_path.open(encoding="utf-8") as handle:
                baseline = json.load(handle)
            mutations = (
                ("schema_version", 2, "metadata_schema_version_mismatch"),
                ("particle_type", "pion", "particle_type_mismatch"),
                ("Q2_token", "2p1", "q2_token_mismatch"),
                ("W_token", "2p95", "w_token_mismatch"),
                ("requested_num_t_bins", 99, "requested_bin_count_mismatch"),
                ("actual_num_t_bins", 99, "actual_bin_count_mismatch"),
                ("binning_config_hash", "wrong", "binning_config_hash_mismatch"),
            )
            for key, value, expected_reason in mutations:
                with self.subTest(key=key):
                    payload = dict(baseline)
                    payload[key] = value
                    with metadata_path.open("w", encoding="utf-8") as handle:
                        json.dump(payload, handle)
                    validated = self.find_bins._validate_authoritative_t_interval(self._inp(), "low")
                    self.assertIn(expected_reason, validated["validation_rejection_reasons"])

    def test_signed_cancellation_keeps_raw_support_meaningful(self):
        records = [{
            "source_label": "prompt", "entry_index": index, "adj_t": value,
            "phi_value": 0.1 + index * 0.1, "physical_coefficient": weight,
        } for index, (value, weight) in enumerate(((0.0, 1.0), (0.3, -1.0), (0.7, 1.0), (1.0, -1.0)))]
        support = self.find_bins._support_summaries(
            np.array([record["adj_t"] for record in records]),
            np.array([record["physical_coefficient"] for record in records]),
            np.array([0.0, 1.0]),
        )
        self.assertEqual(support["raw_event_count_by_t_bin"].tolist(), [4])
        self.assertEqual(support["signed_weighted_yield_by_t_bin"].tolist(), [0.0])
        self.assertEqual(support["absolute_weighted_support_by_t_bin"].tolist(), [4.0])


if __name__ == "__main__":
    unittest.main()
