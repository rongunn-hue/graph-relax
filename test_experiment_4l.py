import json
import math
import tempfile
import unittest
from pathlib import Path

import experiment_4l as e4l
import graph_relax as gr

ROOT = Path(__file__).parent
INPUT = ROOT / "circuit_examples" / "iamp.circuit"
OUTPUT = ROOT / "output"


class Experiment4LTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.circuit = gr.parse_circuit(INPUT)
        cls.graph, cls.embedding, _, cls.old = e4l.load_frozen(cls.circuit, OUTPUT)
        cls.face_index, cls.outer, cls.faces = e4l.existing_outer_face(cls.embedding, cls.old)
        cls.points, cls.interior = e4l.barycentric_coordinates(cls.graph, cls.outer)

    def test_frozen_graph_embedding_and_existing_outer_face(self):
        self.assertEqual((34, 44), (self.graph.number_of_nodes(), self.graph.number_of_edges()))
        self.assertEqual(12, len(self.faces))
        self.assertEqual(4, self.face_index)
        self.assertEqual(14, len(self.outer))

    def test_direct_barycentric_equilibrium(self):
        for vertex in self.interior:
            neighbors = list(self.graph.neighbors(vertex))
            mean = (sum(self.points[n][0] for n in neighbors)/len(neighbors),
                    sum(self.points[n][1] for n in neighbors)/len(neighbors))
            self.assertLess(math.dist(self.points[vertex], mean), 1e-9)

    def test_degeneracy_is_scale_invariant(self):
        groups, pairs, zero_edges = e4l.diagnose(self.circuit, self.graph, self.points)
        self.assertEqual(2, len(pairs))
        self.assertEqual(7, len(zero_edges))
        for scale in (1, 72, 1_000_000):
            for first, second in zero_edges:
                self.assertEqual(tuple(value*scale for value in self.points[first]),
                                 tuple(value*scale for value in self.points[second]))

    def test_saved_failure_is_complete_and_deterministic(self):
        validation = json.loads((OUTPUT / "planar_barycentric_validation.json").read_text())
        self.assertFalse(validation["fully_valid"])
        self.assertFalse(validation["uniform_scaling_can_resolve"])
        self.assertIsNone(validation["selected_physical_scale"])
        with tempfile.TemporaryDirectory(prefix="graph-relax-4l-test-") as tmp:
            target = Path(tmp)
            for name in ("node_planar_embedding.json", "planar_node_initial_coordinates.json"):
                (target/name).write_bytes((OUTPUT/name).read_bytes())
            first = e4l.run(INPUT, target)["hashes"]
            second = e4l.run(INPUT, target)["hashes"]
            self.assertEqual(first, second)


if __name__ == "__main__":
    unittest.main()
