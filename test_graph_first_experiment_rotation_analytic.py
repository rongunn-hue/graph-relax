import math
import unittest
from pathlib import Path

import graph_first_experiment_endpoint_rotation as sweep
import graph_first_experiment_rotation_analytic as analytic


ROOT = Path(__file__).parent
BASE = ROOT/"output"/"graph_first"


class AnalyticRotationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.graph, cls.start, _, _ = sweep.load_start(
            BASE/"iamp_graph.json", BASE/"iamp_nearness_mds.json",
            BASE/"iamp_nearness_mds_report.json")

    def test_circle_circle_events_satisfy_both_circles(self):
        angles, status = analytic.circle_circle_event_angles((0, 0), 5, (4, 0), 3)
        self.assertEqual("TWO_INTERSECTIONS", status)
        for angle in angles:
            point = (5*math.cos(angle), 5*math.sin(angle))
            self.assertAlmostEqual(5, math.dist((0, 0), point), places=9)
            self.assertAlmostEqual(3, math.dist((4, 0), point), places=9)

    def test_circle_segment_events_lie_on_circle_and_segment(self):
        angles = analytic.circle_segment_event_angles((0, 0), 2, (-3, 0), (3, 0))
        self.assertEqual(2, len(angles))
        self.assertEqual([-2.0, 2.0], sorted(round(2*math.cos(angle), 9) for angle in angles))

    def test_analytic_search_is_monotonic_and_preserves_radius(self):
        _, diagnostics, moves, _, history, stats = analytic.analytic_search(self.graph, self.start)
        counts = [len(items) for items in history]
        self.assertTrue(all(after < before for before, after in zip(counts, counts[1:])))
        self.assertGreater(stats["candidate_evaluations"], 0)
        for move in moves:
            self.assertAlmostEqual(move["preserved_radius"],
                                   move["preserved_edge_length_before"], places=12)
        self.assertEqual(counts[-1], diagnostics["proper_unrelated_edge_crossing_count"])


if __name__ == "__main__":
    unittest.main()
