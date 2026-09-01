import unittest
from pathlib import Path

import experiment_4e as e4e
import experiment_4f as e4f
import graph_relax as gr


ROOT = Path(__file__).parent
INPUT = ROOT / "test_data" / "iamp.circuit"
OUTPUT = ROOT / "output"


class Experiment4FTests(unittest.TestCase):
    def setUp(self):
        self.circuit = gr.parse_circuit(INPUT)
        self.start = e4e.load_closeness_initial(self.circuit, OUTPUT)

    def test_baseline_and_first_horizontal_equivalence(self):
        actual, _ = e4e.verify_source(self.circuit, self.start, OUTPUT)
        self.assertEqual((23, 0), (actual["net_crossings"], actual["component_overlaps"]))
        direct = e4e.horizontal_untangle(self.circuit, self.start)
        through_driver = e4f.run_axis_pass(self.circuit, self.start, "H")
        self.assertEqual(direct, through_driver)

    def test_vertical_candidates_fix_x_and_eligibility(self):
        ref = self.circuit.refs[0]
        slide = e4f.y_slide_candidate(self.start, ref, -1, e4f.INITIAL_Y_STEP)
        self.assertEqual(self.start[ref][0], slide[ref][0])
        found = False
        refs = sorted(self.circuit.refs)
        for index, first in enumerate(refs):
            for second in refs[index + 1:]:
                if e4f.horizontal_spans_overlap(self.start[first], self.start[second]):
                    swapped = e4f.vertical_swap_candidate(self.start, first, second)
                    self.assertEqual(self.start[first][0], swapped[first][0])
                    self.assertEqual(self.start[second][0], swapped[second][0])
                    found = True
                    break
            if found:
                break
        self.assertTrue(found)

    def test_vertical_overlap_rejection_and_best_determinism(self):
        forced = dict(self.start)
        forced[self.circuit.refs[1]] = forced[self.circuit.refs[0]]
        self.assertTrue(e4e.has_overlap(self.circuit, forced))
        first = e4f.best_vertical_move(self.circuit, self.start, e4f.INITIAL_Y_STEP)
        second = e4f.best_vertical_move(self.circuit, self.start, e4f.INITIAL_Y_STEP)
        self.assertEqual(first, second)

    def test_alternation_convergence_logging_and_release(self):
        circuit = gr.Circuit("tiny", ("A", "B"), {"A": ("1",), "B": ("1",)},
                             (gr.Net("N", ("A.1", "B.1")),))
        start = {"A": (-100.0, 0.0), "B": (100.0, 0.0)}
        result = e4f.alternating_untangle(circuit, start, maximum_passes=4)
        self.assertTrue(result[3]["converged"])
        self.assertEqual(["H", "V", "H"], [p["axis"] for p in result[1]])
        self.assertEqual(list(range(1, len(result[2]) + 1)),
                         [m["global_move_number"] for m in result[2]])
        released = e4f.release_preprocessing(result[0])
        self.assertIsNot(released, result[0])
        released["A"] = (0.0, 0.0)
        self.assertEqual((0.0, 0.0), released["A"])

    def test_short_repeated_run_determinism(self):
        first = e4f.vertical_untangle(self.circuit, self.start,
                                      initial_step=e4f.INITIAL_Y_STEP,
                                      minimum_step=e4f.INITIAL_Y_STEP)
        second = e4f.vertical_untangle(self.circuit, self.start,
                                       initial_step=e4f.INITIAL_Y_STEP,
                                       minimum_step=e4f.INITIAL_Y_STEP)
        self.assertEqual(first, second)
        self.assertTrue(all(first[0][ref][0] == self.start[ref][0] for ref in self.circuit.refs))


if __name__ == "__main__":
    unittest.main()
