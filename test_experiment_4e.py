import tempfile
import unittest
from pathlib import Path

import experiment_4e as e4e
import graph_relax as gr


ROOT = Path(__file__).parent
INPUT = ROOT / "test_data" / "iamp.circuit"
OUTPUT = ROOT / "output"


class Experiment4ETests(unittest.TestCase):
    def setUp(self):
        self.circuit = gr.parse_circuit(INPUT)
        self.start = e4e.load_closeness_initial(self.circuit, OUTPUT)

    def test_exact_stored_start_reproduction(self):
        actual, expected = e4e.verify_source(self.circuit, self.start, OUTPUT)
        self.assertEqual(23, actual["net_crossings"])
        self.assertEqual(0, actual["component_overlaps"])
        self.assertLess(abs(actual["objective"] - expected["objective"]),
                        e4e.SOURCE_NUMERIC_TOLERANCE)

    def test_slide_is_x_only_and_overlap_is_rejected(self):
        ref = self.circuit.refs[0]
        candidate = e4e.slide_candidate(self.start, ref, -1, e4e.INITIAL_X_STEP)
        self.assertEqual(self.start[ref][1], candidate[ref][1])
        forced = dict(self.start)
        other = self.circuit.refs[1]
        forced[other] = forced[ref]
        self.assertTrue(e4e.has_overlap(self.circuit, forced))

    def test_legal_horizontal_swap_keeps_y(self):
        found = False
        refs = sorted(self.circuit.refs)
        for index, first in enumerate(refs):
            for second in refs[index + 1:]:
                if not e4e.vertical_spans_overlap(self.start[first], self.start[second]):
                    continue
                candidate = e4e.swap_candidate(self.start, first, second)
                if not e4e.has_overlap(self.circuit, candidate):
                    self.assertEqual(self.start[first][1], candidate[first][1])
                    self.assertEqual(self.start[second][1], candidate[second][1])
                    self.assertEqual(self.start[first][0], candidate[second][0])
                    found = True
                    break
            if found:
                break
        self.assertTrue(found)

    def test_lexicographic_comparison_and_best_move_determinism(self):
        self.assertTrue(e4e.improves((2, 100.0), (3, 1.0)))
        self.assertTrue(e4e.improves((3, 99.0), (3, 100.0)))
        self.assertFalse(e4e.improves((4, 1.0), (3, 100.0)))
        first = e4e.best_horizontal_move(self.circuit, self.start, e4e.INITIAL_X_STEP)
        second = e4e.best_horizontal_move(self.circuit, self.start, e4e.INITIAL_X_STEP)
        self.assertEqual(first, second)

    def test_move_logging_release_and_short_repeatability(self):
        first = e4e.horizontal_untangle(self.circuit, self.start,
                                        initial_step=e4e.INITIAL_X_STEP,
                                        minimum_step=e4e.INITIAL_X_STEP)
        second = e4e.horizontal_untangle(self.circuit, self.start,
                                         initial_step=e4e.INITIAL_X_STEP,
                                         minimum_step=e4e.INITIAL_X_STEP)
        self.assertEqual(first, second)
        self.assertTrue(all(first[0][ref][1] == self.start[ref][1] for ref in self.circuit.refs))
        self.assertEqual(list(range(1, len(first[1]) + 1)),
                         [move["move_number"] for move in first[1]])
        released = e4e.release_preprocessing(first[0])
        self.assertIsNot(released, first[0])
        released[self.circuit.refs[0]] = (0.0, 0.0)
        self.assertEqual((0.0, 0.0), released[self.circuit.refs[0]])


if __name__ == "__main__":
    unittest.main()
