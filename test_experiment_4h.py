import unittest
from pathlib import Path

import experiment_4g as e4g
import experiment_4h as e4h
import graph_relax as gr


ROOT = Path(__file__).parent
INPUT = ROOT / "test_data" / "iamp.circuit"
OUTPUT = ROOT / "output"


class Experiment4HTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.circuit = gr.parse_circuit(INPUT)
        cls.positions = e4g.load_frozen_state(cls.circuit, OUTPUT)
        cls.actual, _ = e4g.verify_frozen_state(cls.circuit, cls.positions, OUTPUT)
        cls.initial = gr.initial_junctions(cls.circuit, cls.positions)
        cls.bounds = e4h.junction_bounds(cls.circuit, cls.positions)

    def test_baseline_movable_nets_and_centroids(self):
        self.assertEqual((7, 0), (self.actual["net_crossings"], self.actual["component_overlaps"]))
        expected = {net.name for net in self.circuit.nets if len(net.terminals) >= 3}
        self.assertEqual(expected, set(e4h.movable_nets(self.circuit)))
        self.assertEqual(self.initial, gr.initial_junctions(self.circuit, self.positions))

    def test_move_generation_bounds_and_body_exclusion(self):
        name = e4h.movable_nets(self.circuit)[0]
        moved = e4h.move_candidate(self.initial, name, -1, 0, e4h.MINIMUM_JUNCTION_STEP)
        self.assertEqual(self.initial[name][1], moved[name][1])
        self.assertNotEqual(self.initial[name][0], moved[name][0])
        self.assertTrue(e4h.in_bounds(self.initial[name], self.bounds[name]))
        unrelated = next(ref for ref in self.circuit.refs
                         if ref not in {t.rsplit('.', 1)[0] for n in self.circuit.nets
                                       if n.name == name for t in n.terminals})
        self.assertFalse(e4h.junction_outside_unrelated_bodies(
            self.circuit, self.positions, name, self.positions[unrelated]))

    def test_lexicographic_comparison_and_best_determinism(self):
        self.assertTrue(e4h.improves((6, 9999.0), (7, 1.0)))
        self.assertTrue(e4h.improves((7, 10.0), (7, 11.0)))
        self.assertFalse(e4h.improves((8, 0.0), (7, 100.0)))
        first = e4h.best_junction_move(self.circuit, self.positions, self.initial,
                                       self.bounds, e4h.INITIAL_JUNCTION_STEP)
        second = e4h.best_junction_move(self.circuit, self.positions, self.initial,
                                        self.bounds, e4h.INITIAL_JUNCTION_STEP)
        self.assertEqual(first, second)

    def test_search_component_immutability_logging_and_repeatability(self):
        before = dict(self.positions)
        first = e4h.junction_untangle(self.circuit, self.positions, self.initial, self.bounds,
                                      initial_step=e4h.INITIAL_JUNCTION_STEP,
                                      minimum_step=e4h.INITIAL_JUNCTION_STEP)
        second = e4h.junction_untangle(self.circuit, self.positions, self.initial, self.bounds,
                                       initial_step=e4h.INITIAL_JUNCTION_STEP,
                                       minimum_step=e4h.INITIAL_JUNCTION_STEP)
        self.assertEqual(first, second)
        self.assertEqual(before, self.positions)
        self.assertEqual(list(range(1, len(first[1]) + 1)), [m["move_number"] for m in first[1]])

    def test_crossing_identity_comparison(self):
        originals = __import__('json').loads((OUTPUT / "residual_crossings.json").read_text())["crossings"]
        rows, new_pairs = e4h.compare_original_crossings(
            originals, e4h.junction_crossing_map(self.circuit, self.positions, self.initial))
        self.assertEqual([f"X{i:03d}" for i in range(1, 8)], [row["id"] for row in rows])
        self.assertTrue(all(row["status"] in {"eliminated", "persists", "displaced/replaced"}
                            for row in rows))
        self.assertIsInstance(new_pairs, list)


if __name__ == "__main__":
    unittest.main()
