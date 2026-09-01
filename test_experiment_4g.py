import json
import tempfile
import unittest
from pathlib import Path

import experiment_4g as e4g
import graph_relax as gr


ROOT = Path(__file__).parent
INPUT = ROOT / "test_data" / "iamp.circuit"
OUTPUT = ROOT / "output"


class Experiment4GTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.circuit = gr.parse_circuit(INPUT)
        cls.positions = e4g.load_frozen_state(cls.circuit, OUTPUT)
        cls.actual, _ = e4g.verify_frozen_state(cls.circuit, cls.positions, OUTPUT)
        cls.crossings, cls.interactions, cls.summary = e4g.diagnose(cls.circuit, cls.positions)

    def test_exact_frozen_baseline_and_enumeration(self):
        self.assertEqual((7, 0), (self.actual["net_crossings"], self.actual["component_overlaps"]))
        self.assertEqual([f"X{i:03d}" for i in range(1, 8)], [c["id"] for c in self.crossings])
        order = [(c["intersection"]["y"], c["intersection"]["x"]) for c in self.crossings]
        self.assertEqual(order, sorted(order))

    def test_relevant_components_and_identity_tracking(self):
        for crossing in self.crossings:
            expected = set(crossing["net_a_components"]) | set(crossing["net_b_components"])
            self.assertEqual(expected, set(crossing["relevant_components"]))
            for axis in ("x", "y"):
                candidate = crossing[f"best_{axis}_candidate"]
                if candidate is not None:
                    self.assertTrue(candidate["target_crossing_eliminated"])
                    self.assertTrue(set(candidate["move"]["components"]) <= expected)

    def test_axis_search_does_not_mutate_and_classifies(self):
        before = dict(self.positions)
        baseline = e4g.crossing_map(self.circuit, self.positions)
        cache = {}
        crossing = self.crossings[0]
        e4g.diagnostic_axis_search(self.circuit, self.positions, crossing, "X", baseline, cache)
        e4g.diagnostic_axis_search(self.circuit, self.positions, crossing, "Y", baseline, cache)
        self.assertEqual(before, self.positions)
        self.assertIn(crossing["classification"], {"X", "Y", "XY", "STRUCTURAL"})

    def test_interaction_matrix_and_report_determinism(self):
        ids = [c["id"] for c in self.crossings]
        self.assertEqual(ids, list(self.interactions))
        for row in self.interactions.values():
            self.assertEqual(ids, list(row))
            self.assertTrue(set(row.values()) <= {"E", "P", "R", "N/A"})
        first = e4g.text_report(self.actual, self.crossings, self.interactions, self.summary)
        second = e4g.text_report(self.actual, self.crossings, self.interactions, self.summary)
        self.assertEqual(first, second)

    def test_repeated_diagnosis_determinism(self):
        again = e4g.diagnose(self.circuit, dict(self.positions))
        self.assertEqual((self.crossings, self.interactions, self.summary), again)


if __name__ == "__main__":
    unittest.main()
