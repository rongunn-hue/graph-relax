import json
import tempfile
import unittest
from pathlib import Path

import graph_first_experiment_global_linear as experiment


class GlobalLinearExperimentTests(unittest.TestCase):
    def run_experiment(self, output):
        base = Path("output/graph_first")
        return experiment.run(base / "iamp_graph.json", base / "iamp_nearness_rotation_direct_report.json", Path(output))

    def test_one_frozen_normalized_global_solve(self):
        with tempfile.TemporaryDirectory() as directory:
            report = self.run_experiment(directory)
            self.assertEqual(report["unknown_count"], 68)
            self.assertEqual(report["matrix"]["shape"][1], 68)
            self.assertEqual(report["matrix"]["shape"][0], sum(report["active_relation_counts"].values()) - len(report["skipped_rows"]))
            self.assertEqual(report["crowding"]["violation_count_before"], sum(report["crowding"]["counts_before"].values()))
            self.assertLess(abs(report["solution"]["mean_dx"]), 1e-10)
            self.assertLess(abs(report["solution"]["mean_dy"]), 1e-10)
            self.assertEqual(report["determinism"], {"runs": 2, "identical_serialized_calculation": True})

    def test_authoritative_reconnect_and_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            report = self.run_experiment(directory)
            validation = report["final_validation"]
            self.assertEqual((validation["vertices"], validation["edges"]), (34, 44))
            self.assertTrue(validation["endpoint_pairs_identical"])
            self.assertTrue(validation["reconnected_edges_valid"])
            for name in ("iamp_global_linear_start.svg", "iamp_global_linear_final.svg",
                         "iamp_global_linear_diagnostic.svg", "iamp_global_linear_report.json"):
                self.assertTrue((Path(directory) / name).is_file())
            loaded = json.loads((Path(directory) / "iamp_global_linear_report.json").read_text())
            self.assertEqual(loaded["final_coordinates"], report["final_coordinates"])


if __name__ == "__main__":
    unittest.main()
