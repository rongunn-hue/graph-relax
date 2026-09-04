import math
import json
import tempfile
import unittest
from pathlib import Path

import placement_geometry as geometry
import graph_first_experiment_endpoint_rotation as experiment


ROOT = Path(__file__).parent
BASE = ROOT / "output" / "graph_first"


class EndpointRotationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.graph, cls.start, cls.initial, _ = experiment.load_start(
            BASE / "iamp_graph.json", BASE / "iamp_nearness_mds.json",
            BASE / "iamp_nearness_mds_report.json")

    def test_start_satisfies_permanent_geometry_rules(self):
        self.assertEqual(0, self.initial["coincident_vertex_pair_count"])
        self.assertEqual(0, self.initial["vertex_on_unrelated_edge_interior_count"])
        self.assertGreaterEqual(self.initial["minimum_node_distance"] + 1e-9,
                                geometry.MIN_NODE_DISTANCE)

    def test_rotation_preserves_pivot_edge_length_and_moves_only_endpoint(self):
        moving, pivot = next(iter(self.graph.edges))
        candidate = dict(self.start)
        candidate[moving] = experiment.rotated_position(self.start[moving], self.start[pivot], 17)
        self.assertAlmostEqual(math.dist(self.start[moving], self.start[pivot]),
                               math.dist(candidate[moving], candidate[pivot]), places=9)
        self.assertTrue(all(candidate[node] == self.start[node] for node in self.graph if node != moving))

    def test_search_is_monotonic_and_every_move_preserves_length(self):
        report = json.loads((BASE / "iamp_nearness_rotation_report.json").read_text())
        moves = report["moves"]
        counts = [len(items) for items in report["crossing_history"]]
        self.assertTrue(all(after < before for before, after in zip(counts, counts[1:])))
        for move in moves:
            self.assertAlmostEqual(move["preserved_edge_length_before"],
                                   move["preserved_edge_length_after"], places=8)
            self.assertGreaterEqual(move["minimum_node_distance_after"] + 1e-9,
                                    geometry.MIN_NODE_DISTANCE)
        self.assertEqual(counts[-1], report["final_crossing_count"])

    def test_search_is_deterministic_on_crossed_generic_fixture(self):
        graph = __import__("networkx").Graph()
        graph.add_nodes_from((node, {"type": "COMPONENT", "name": node})
                             for node in "abcd")
        graph.add_edges_from((("a", "b"), ("c", "d")))
        positions = {"a": (-2.0, -2.0), "b": (2.0, 2.0),
                     "c": (-2.0, 2.0), "d": (2.0, -2.0)}
        first = experiment.run_search(graph, positions)
        second = experiment.run_search(graph, positions)
        self.assertEqual(first, second)


if __name__ == "__main__":
    unittest.main()
