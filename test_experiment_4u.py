import json
import math
import re
import unittest
from pathlib import Path

import experiment_4r_external_tent as abstract_graph
import experiment_4u as e4u
import graph_relax as gr


ROOT = Path(__file__).parent
INPUT = ROOT / "circuit_examples" / "iamp.circuit"
OUTPUT = ROOT / "output"


class TwoPole4UTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        circuit = gr.parse_circuit(INPUT)
        cls.graph, cls.embedding, cls.initial, cls.edges = abstract_graph.load_abstract_start(circuit, OUTPUT)
        cls.model = e4u.pair_model(cls.initial, "component:J_RA", "component:RMIN4")

    def test_zero_reproduces_original_and_poles_stay_fixed(self):
        self.assertEqual(e4u.positions_at_height(self.initial, self.model, 0.0), self.initial)
        positions = e4u.positions_at_height(self.initial, self.model, 25.0)
        for pole in (self.model["pole_a"], self.model["pole_b"]):
            self.assertEqual(positions[pole], self.initial[pole])

    def test_exact_two_circle_distances(self):
        height = 40.0
        positions = e4u.positions_at_height(self.initial, self.model, height)
        for vertex in self.initial:
            if vertex in (self.model["pole_a"], self.model["pole_b"]):
                continue
            for pole in (self.model["pole_a"], self.model["pole_b"]):
                original = math.dist(self.initial[vertex], self.initial[pole])
                expected = math.sqrt(original * original - height * height)
                self.assertAlmostEqual(math.dist(positions[vertex], self.initial[pole]), expected, places=7)

    def test_complete_pair_output_and_best_valid(self):
        document = json.loads((OUTPUT / "two_pole_all_pairs.json").read_text())
        self.assertEqual(document["pair_count"], 561)
        self.assertEqual(len(document["ranking"]), 561)
        self.assertTrue(document["ranking"][0]["abstract_graph_valid"])

    def test_all_visual_frames_share_viewport(self):
        paths = sorted(OUTPUT.glob("two_pole_*.svg"))
        paths = [path for path in paths if path.name != "two_pole_contraction.svg"]
        viewports = {re.search(r'viewBox="([^"]+)"', path.read_text()).group(1) for path in paths}
        self.assertEqual(len(viewports), 1)


if __name__ == "__main__":
    unittest.main()
