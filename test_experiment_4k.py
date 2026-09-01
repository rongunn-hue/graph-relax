import json
import tempfile
import unittest
from pathlib import Path

import experiment_4j as e4j
import experiment_4k as e4k
import graph_relax as gr

ROOT = Path(__file__).parent
INPUT = ROOT / "circuit_examples" / "iamp.circuit"
OUTPUT = ROOT / "output"


class Experiment4KTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.circuit = gr.parse_circuit(INPUT)
        cls.loaded = e4j.load_stored(cls.circuit, OUTPUT)

    def test_exact_grid_and_independent_scaling(self):
        self.assertEqual(tuple(range(1, 53)), e4k.SCALE_VALUES)
        self.assertEqual(2704, len(e4k.SCALE_PAIRS))
        self.assertEqual(2704, len(set(e4k.SCALE_PAIRS)))
        a = e4k.geometry_at_scales(self.circuit, self.loaded[5], self.loaded[6], 52, 52)
        b = e4k.geometry_at_scales(self.circuit, self.loaded[5], self.loaded[6], 26, 52)
        ref = next(ref for ref in self.circuit.refs if a[0][ref][0] != 0)
        self.assertEqual(a[0][ref][0] / 2, b[0][ref][0])
        self.assertEqual(a[0][ref][1], b[0][ref][1])
        self.assertEqual((60.0, 40.0), (gr.WIDTH, gr.HEIGHT))

    def test_exact_4j_baseline(self):
        result = e4k.verify_baseline(self.circuit, self.loaded[3], self.loaded[4],
                                     self.loaded[5], self.loaded[6], OUTPUT)
        self.assertTrue(result[0]["fully_valid"])
        self.assertAlmostEqual(24262.045496, result[3]["total_connection_length"], places=5)

    def test_saved_complete_sweep_and_selection(self):
        doc = json.loads((OUTPUT / "planar_xy_scale_sweep.json").read_text())
        self.assertEqual(2704, doc["tested_pair_count"])
        self.assertEqual(2704, len(doc["records"]))
        valid = [r for r in doc["records"] if r["fully_valid"]]
        expected = min(valid, key=lambda r: (r["drawing_area"], r["total_connection_length"],
                                             max(r["sx"], r["sy"]), r["sx"], r["sy"]))
        self.assertEqual({"sx": expected["sx"], "sy": expected["sy"]}, doc["selected_pair"])
        self.assertEqual(len(valid), doc["fully_valid_pair_count"])
        self.assertEqual(min(r["sx"] for r in valid), doc["smallest_sx_in_any_valid_pair"])
        self.assertEqual(min(r["sy"] for r in valid), doc["smallest_sy_in_any_valid_pair"])

    def test_saved_validity_neighbors_and_cyclic_order(self):
        doc = json.loads((OUTPUT / "planar_xy_scale_sweep.json").read_text())
        chosen = next(r for r in doc["records"] if (r["sx"], r["sy"]) ==
                      (doc["selected_pair"]["sx"], doc["selected_pair"]["sy"]))
        self.assertTrue(chosen["fully_valid"])
        self.assertIn(chosen["cyclic_order_orientation"], ("same", "global_reflection"))
        self.assertEqual(44, chosen["realized_incidences"])
        self.assertIn("lower_x", doc["immediate_lower_neighbors"])
        self.assertIn("lower_y", doc["immediate_lower_neighbors"])

    def test_deterministic_serialization_from_saved_sweep(self):
        # Full ten-process determinism is performed by the experiment run; this
        # targeted test verifies stable JSON reserialization of its full grid.
        raw = (OUTPUT / "planar_xy_scale_sweep.json").read_bytes()
        with tempfile.TemporaryDirectory(prefix="graph-relax-4k-test-") as tmp:
            path = Path(tmp) / "copy.json"
            gr.write_json(path, json.loads(raw))
            self.assertEqual(raw, path.read_bytes())


if __name__ == "__main__":
    unittest.main()
