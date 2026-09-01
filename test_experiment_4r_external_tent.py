import json
import math
import tempfile
import unittest
from pathlib import Path

import experiment_4r_external_tent as ext
import graph_relax as gr


ROOT = Path(__file__).parent
INPUT = ROOT / "circuit_examples" / "iamp.circuit"
OUTPUT = ROOT / "output"


class AbstractExternalTent4RTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.circuit = gr.parse_circuit(INPUT)
        cls.graph, cls.embedding, cls.initial, cls.edges = ext.load_abstract_start(cls.circuit, OUTPUT)
        cls.model = ext.external_model(cls.initial)

    def test_h_zero_and_external_pole(self):
        center, farthest, radius, pole, _, _, _, _ = self.model
        zero = ext.positions_at_height(self.initial, self.model, 0.0)
        self.assertTrue(all(math.dist(zero[v], self.initial[v]) < 2e-6 for v in self.initial))
        self.assertAlmostEqual(math.dist(center, pole), 2 * radius, places=8)
        self.assertEqual(len(self.initial), 34)
        self.assertEqual(len(self.edges), 44)

    def test_rays_and_cable_equation(self):
        _, _, _, pole, cables, rays, hmax, _ = self.model
        low = ext.positions_at_height(self.initial, self.model, 0.2 * hmax)
        high = ext.positions_at_height(self.initial, self.model, 0.8 * hmax)
        for vertex in self.initial:
            self.assertLess(math.dist(high[vertex], pole), math.dist(low[vertex], pole))
            dx, dy = high[vertex][0] - pole[0], high[vertex][1] - pole[1]
            self.assertLess(abs(dx * rays[vertex][1] - dy * rays[vertex][0]), 2e-6)
            self.assertLess(abs(math.dist(high[vertex], pole) ** 2 + (0.8 * hmax) ** 2 - cables[vertex] ** 2), 0.02)

    def test_selected_saved_result_is_abstract_and_valid(self):
        document = json.loads((OUTPUT / "planar_external_tent_validation.json").read_text())
        self.assertFalse(document["physical_validation_participated"])
        self.assertTrue(document["selected_validation"]["abstract_graph_valid"])
        self.assertEqual(document["selected_validation"]["unrelated_edge_crossings"], 0)

    def test_deterministic_execution(self):
        def prepared_directory():
            directory = Path(tempfile.mkdtemp(prefix="abstract-external-tent-"))
            for name in ("node_planar_embedding.json", "planar_xy_compact_coordinates.json"):
                (directory / name).write_bytes((OUTPUT / name).read_bytes())
            return directory
        first, second = prepared_directory(), prepared_directory()
        self.assertEqual(ext.run(INPUT, first)["hashes"], ext.run(INPUT, second)["hashes"])


if __name__ == "__main__":
    unittest.main()
