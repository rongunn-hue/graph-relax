import json
import unittest
from pathlib import Path

import experiment_4n as e4n
import graph_relax as gr

ROOT = Path(__file__).parent
INPUT = ROOT / "circuit_examples" / "iamp.circuit"
OUTPUT = ROOT / "output"


class Experiment4NTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.circuit = gr.parse_circuit(INPUT)
        cls.graph, cls.embedding, cls.positions, cls.order = e4n.load_start(cls.circuit, OUTPUT)

    def test_exact_start_steps_and_order(self):
        self.assertEqual((256., 128., 64., 32., 16., 8.), e4n.STEPS)
        self.assertEqual(34, len(self.order))
        self.assertEqual(self.order, sorted(self.graph, key=lambda v: (-self.graph.degree(v), v)))
        validation, _ = e4n.complete_validate(self.circuit, self.graph, self.embedding, self.positions)
        self.assertTrue(validation["valid"])

    def test_directions_are_deterministic_and_coarse(self):
        first = e4n.candidate_directions(self.graph, self.positions, self.order[0])
        second = e4n.candidate_directions(self.graph, self.positions, self.order[0])
        self.assertEqual(first, second)
        self.assertTrue(all("unit" in row for row in first))

    def test_saved_six_pass_result_is_valid(self):
        metrics = json.loads((OUTPUT / "planar_local_compact_metrics.json").read_text())
        validation = json.loads((OUTPUT / "planar_local_compact_validation.json").read_text())
        self.assertEqual(6, len(metrics["passes"]))
        self.assertEqual(list(e4n.STEPS), metrics["steps"])
        self.assertTrue(validation["exactly_six_passes"])
        self.assertTrue(validation["final_complete_validation"]["valid"])
        self.assertTrue(all(row["valid"] for row in validation["complete_pass_validations"]))


if __name__ == "__main__":
    unittest.main()
