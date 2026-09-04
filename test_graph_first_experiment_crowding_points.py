import unittest
from pathlib import Path

import graph_first_experiment_crowding_analysis as crowding
import graph_first_experiment_crowding_points as points


BASE = Path(__file__).parent/"output"/"graph_first"


class CrowdingPointTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.vertices, cls.edges, cls.coordinates, _, _ = crowding.load_inputs(
            BASE/"iamp_graph.json", BASE/"iamp_nearness_rotation_direct_report.json")
        cls.analysis = crowding.analyze(cls.vertices, cls.edges, cls.coordinates)

    def test_one_exact_midpoint_per_violation(self):
        events = points.localize_events(self.analysis)
        self.assertEqual(len(self.analysis["violations"]), len(events))
        for event in events:
            expected = [(event["pA"][0]+event["pB"][0])/2,
                        (event["pA"][1]+event["pB"][1])/2]
            self.assertEqual(expected, event["crowding_point"])

    def test_cluster_dimensions_use_event_points_only(self):
        events = points.localize_events(self.analysis)
        clusters, _, _ = points.cluster_events(events)
        by_id = {event["violation_id"]: event for event in events}
        for cluster in clusters:
            locations = [by_id[item]["crowding_point"]
                         for item in cluster["member_violation_ids"]]
            self.assertEqual(min(point[0] for point in locations),
                             cluster["crowding_point_bounding_box"]["min_x"])

    def test_measurement_preserves_coordinates(self):
        before = dict(self.coordinates)
        events = points.localize_events(self.analysis)
        points.cluster_events(events)
        self.assertEqual(before, self.coordinates)


if __name__ == "__main__":
    unittest.main()
