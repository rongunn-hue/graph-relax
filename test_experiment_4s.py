import json
import math
import tempfile
import unittest
from pathlib import Path

import experiment_4r_external_tent as abstract_graph
import experiment_4s as e4s
import graph_relax as gr


ROOT = Path(__file__).parent
INPUT = ROOT / "circuit_examples" / "iamp.circuit"
OUTPUT = ROOT / "output"


class SinglePole4STests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        circuit = gr.parse_circuit(INPUT)
        cls.graph, cls.embedding, cls.initial, cls.edges = abstract_graph.load_abstract_start(circuit, OUTPUT)

    def test_zero_reproduces_4k_and_every_vertex_stays_on_ray(self):
        zero = e4s.contracted_positions(self.initial, 0.0)
        self.assertEqual(zero, self.initial)
        for ratio in e4s.HEIGHT_RATIOS[1:]:
            positions = e4s.contracted_positions(self.initial, ratio)
            for vertex in self.initial:
                original = (self.initial[vertex][0] - e4s.POLE[0], self.initial[vertex][1] - e4s.POLE[1])
                moved = (positions[vertex][0] - e4s.POLE[0], positions[vertex][1] - e4s.POLE[1])
                self.assertAlmostEqual(original[0] * moved[1] - original[1] * moved[0], 0.0, places=7)

    def test_all_pairwise_distances_scale_exactly(self):
        baseline = e4s.pairwise_distances(self.initial)
        for ratio in e4s.HEIGHT_RATIOS:
            scale = e4s.contraction_factor(ratio)
            distances = e4s.pairwise_distances(e4s.contracted_positions(self.initial, ratio))
            self.assertLess(max(abs(distances[pair] / baseline[pair] - scale) for pair in baseline), 2e-10)

    def test_saved_result_is_abstract_and_all_valid(self):
        document = json.loads((OUTPUT / "planar_single_pole_validation.json").read_text())
        self.assertFalse(document["physical_validation_participated"])
        self.assertTrue(document["all_states_valid"])
        self.assertTrue(document["all_similarity_identities_verified"])

    def test_deterministic_serialization(self):
        def prepared():
            directory = Path(tempfile.mkdtemp(prefix="single-pole-"))
            for name in ("node_planar_embedding.json", "planar_xy_compact_coordinates.json"):
                (directory / name).write_bytes((OUTPUT / name).read_bytes())
            return directory
        first, second = prepared(), prepared()
        self.assertEqual(e4s.run(INPUT, first)["hashes"], e4s.run(INPUT, second)["hashes"])


if __name__ == "__main__":
    unittest.main()
