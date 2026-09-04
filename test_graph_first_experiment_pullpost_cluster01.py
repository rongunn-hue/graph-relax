import unittest
from pathlib import Path
import json

import graph_first_experiment_crowding_analysis as crowding
import graph_first_experiment_pullpost_cluster01 as experiment


BASE = Path(__file__).parent/"output"/"graph_first"


class PullPostCluster01Tests(unittest.TestCase):
    def test_single_direct_move_and_reconnection(self):
        vertices, edges, coordinates, _, _ = crowding.load_inputs(
            BASE/"iamp_graph.json", BASE/"iamp_nearness_rotation_direct_report.json")
        points = json.loads((BASE/"iamp_crowding_points_report.json").read_text())
        final, movement = experiment.direct_pullpost_move(vertices, edges, coordinates, points)
        self.assertEqual("ACCEPTED", movement["status"])
        changed = [node for node in coordinates if coordinates[node] != final[node]]
        self.assertEqual(["component:TPB"], changed)
        self.assertGreaterEqual(movement["clearance_after"], 6.0-1e-9)
        self.assertEqual(0, movement["diagnostics"]["proper_unrelated_edge_crossing_count"])
        self.assertTrue(movement["reconnection_validation"]["valid"])


if __name__ == "__main__":
    unittest.main()
