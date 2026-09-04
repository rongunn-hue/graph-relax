import unittest
from pathlib import Path

import graph_first_experiment_crowding_analysis as experiment


BASE = Path(__file__).parent/"output"/"graph_first"


class CrowdingAnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.vertices, cls.edges, cls.coordinates, _, _ = experiment.load_inputs(
            BASE/"iamp_graph.json", BASE/"iamp_nearness_rotation_direct_report.json")

    def test_measurement_does_not_move_coordinates(self):
        before = dict(self.coordinates)
        experiment.analyze(self.vertices, self.edges, self.coordinates)
        self.assertEqual(before, self.coordinates)

    def test_every_item_has_clearance(self):
        result = experiment.analyze(self.vertices, self.edges, self.coordinates)
        self.assertEqual(34, len(result["vertex_clearance"]))
        self.assertEqual(44, len(result["edge_clearance"]))

    def test_relations_exclude_incident_geometry(self):
        result = experiment.analyze(self.vertices, self.edges, self.coordinates)
        for relation in result["relations"]:
            if relation["geometry_type"] == "vertex-edge":
                node = relation["item_a"][len("vertex:"):]
                self.assertNotIn(node, relation["item_b"][len("edge:"):].split("|"))
            elif relation["geometry_type"] == "edge-edge":
                first = set(relation["item_a"][len("edge:"):].split("|"))
                second = set(relation["item_b"][len("edge:"):].split("|"))
                self.assertFalse(first & second)


if __name__ == "__main__":
    unittest.main()
