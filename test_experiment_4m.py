import json
import unittest
from pathlib import Path

import experiment_4m as e4m
import graph_relax as gr

ROOT = Path(__file__).parent
INPUT = ROOT / "circuit_examples" / "iamp.circuit"
OUTPUT = ROOT / "output"


class Experiment4MTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.circuit = gr.parse_circuit(INPUT)
        cls.graph, cls.embedding, _, cls.positions = e4m.load_frozen(cls.circuit, OUTPUT)

    def test_anchor_and_fixed_processing_order(self):
        anchor, degree, order = e4m.anchor_and_order(self.graph)
        self.assertEqual(("component:U1", 8), (anchor, degree))
        self.assertEqual(33, len(order))
        self.assertEqual(order, sorted(order, key=lambda v: (-self.graph.degree(v), v)))

    def test_baseline_and_all_vertex_rotation(self):
        validation, _ = e4m.validate(self.circuit, self.graph, self.embedding, self.positions)
        self.assertTrue(validation["valid"])
        self.assertEqual("global_reflection", validation["all_vertex_cyclic_order_orientation"])
        self.assertEqual(44, validation["realized_incidence_count"])

    def test_candidate_is_on_anchor_ray_at_precision(self):
        start, anchor = (100.0, 50.0), (0.0, 0.0)
        point = e4m.candidate_position(start, anchor, 0.25)
        self.assertEqual((75.0, 37.5), point)
        self.assertEqual(1e-6, e4m.COORDINATE_PRECISION)

    def test_saved_result_is_valid_and_complete(self):
        metrics = json.loads((OUTPUT / "planar_anchor_compact_metrics.json").read_text())
        validation = json.loads((OUTPUT / "planar_anchor_compact_validation.json").read_text())
        self.assertLessEqual(metrics["pass_count"], 20)
        self.assertTrue(validation["anchor_fixed"])
        self.assertTrue(validation["validation"]["valid"])
        self.assertEqual(0, validation["validation"]["different_net_crossings"])
        self.assertEqual(44, validation["validation"]["realized_incidence_count"])
        self.assertEqual(sum(row["vertices_moved"] for row in metrics["passes"]),
                         metrics["accepted_move_count"])


if __name__ == "__main__":
    unittest.main()
