import json
import unittest
from pathlib import Path

import experiment_4n as e4n
import experiment_4o as e4o
import graph_relax as gr

ROOT = Path(__file__).parent
INPUT = ROOT/"circuit_examples"/"iamp.circuit"
OUTPUT = ROOT/"output"


class Experiment4OTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.circuit = gr.parse_circuit(INPUT)
        cls.graph, cls.embedding, cls.positions, cls.order = e4n.load_start(cls.circuit, OUTPUT)

    def test_fixed_threshold_schedule_and_angular_measure(self):
        self.assertEqual(.90, e4o.COMPARABLE_REDUCTION_FRACTION)
        self.assertEqual((256.,128.,64.,32.,16.,8.), e4n.STEPS)
        row = e4o.vertex_angular_measure(self.graph, self.positions, "component:U1")
        self.assertEqual(8, row["degree"])
        self.assertAlmostEqual(360.0, sum(row["gaps"]), places=9)

    def test_selection_prefers_angle_within_ninety_percent(self):
        candidates = [
            {"reduction":100.,"direction_index":0,"angular":{"angular_error":50.,"minimum_gap":10.}},
            {"reduction":90.,"direction_index":1,"angular":{"angular_error":20.,"minimum_gap":20.}},
            {"reduction":89.,"direction_index":2,"angular":{"angular_error":1.,"minimum_gap":100.}}]
        chosen=e4o.select_candidate(self.graph,"component:U1",candidates)
        self.assertEqual(1,chosen["direction_index"])

    def test_saved_result_has_six_valid_passes(self):
        metrics=json.loads((OUTPUT/"planar_circular_compact_metrics.json").read_text())
        validation=json.loads((OUTPUT/"planar_circular_compact_validation.json").read_text())
        self.assertEqual(6,len(metrics["passes"]))
        self.assertTrue(all(row["validation"]["valid"] for row in metrics["passes"]))
        self.assertTrue(validation["final_complete_validation"]["valid"])


if __name__=="__main__":unittest.main()
