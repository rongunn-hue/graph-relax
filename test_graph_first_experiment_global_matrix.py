import tempfile
import unittest
from pathlib import Path

import graph_first_experiment_global_matrix as experiment


class GlobalMatrixExperimentTests(unittest.TestCase):
    def test_matrix_invariants_and_serialized_determinism(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path("output/graph_first")
            report = experiment.run(base/"iamp_graph.json", base/"iamp_nearness_rotation_direct_report.json", Path(directory))
            self.assertEqual(report["matrices"]["C_dimensions"][1], 34)
            self.assertEqual(report["matrices"]["Q_dimensions"], [44, 34])
            self.assertEqual(report["matrices"]["M_dimensions"], [34, 34])
            self.assertEqual(report["matrices"]["M_nullity"], 1)
            self.assertEqual(report["matrices"]["KKT_rank"], 35)
            self.assertLess(report["centroid_constraint"]["maximum_absolute_error"], 1e-9)
            self.assertTrue(report["determinism"]["coordinates_identical_at_17_digit_serialization"])
            self.assertTrue(report["final_validation"]["endpoint_pairs_identical"])
            self.assertTrue(report["final_validation"]["reconnected_edges_valid"])


if __name__ == "__main__":
    unittest.main()
