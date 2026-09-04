import ast
import math
import unittest
from pathlib import Path

import direct_clearance_geometry as clearance
import direct_crossing_repair as repair
import placement_geometry as geometry


class DirectClearanceGeometryTests(unittest.TestCase):
    def test_project_rendering_clearance_parameters(self):
        self.assertEqual(1.5, geometry.EDGE_STROKE_WIDTH)
        self.assertEqual(4.5, geometry.MIN_VISIBLE_EDGE_GAP)
        self.assertEqual(6.0, geometry.REQUIRED_CENTERLINE_CLEARANCE)

    def solve(self, pivot, moving, first, second, required=2.25):
        return clearance.nearest_clearance_escape(
            pivot, math.dist(pivot, moving), moving, first, second, required)

    def test_mathematical_escape_can_be_visually_too_close(self):
        pivot, moving = (0, 0), (10, 0)
        crossed = ((5, -3), (5, 3))
        angle = math.atan2(3, 5) + 1e-6
        almost = (10*math.cos(angle), 10*math.sin(angle))
        self.assertFalse(geometry.proper_segment_crossing(pivot, almost, *crossed))
        self.assertLess(geometry.segment_segment_distance(pivot, almost, *crossed), 2.25)

    def test_capsule_solution_reaches_requested_clearance(self):
        answer = self.solve((0, 0), (10, 0), (5, -3), (5, 3))
        self.assertGreaterEqual(answer["resulting_distance"]+geometry.GEOMETRY_TOLERANCE, 2.25)

    def test_interior_strip_boundary(self):
        answer = self.solve((0, 0), (10, 0), (5, -20), (5, 20))
        self.assertTrue(answer["boundary_feature"].startswith("STRIP"))

    def test_endpoint_a_cap_boundary(self):
        answer = self.solve((0, 0), (10, 0), (5, 0), (8, 4))
        self.assertEqual("CAP_A", answer["boundary_feature"])

    def test_endpoint_b_cap_boundary(self):
        answer = self.solve((0, 0), (10, 0), (2, 4), (5, 0))
        self.assertEqual("CAP_B", answer["boundary_feature"])

    def test_no_sampling_or_iterative_refinement(self):
        source = Path(clearance.__file__).read_text()
        tree = ast.parse(source)
        calls = {node.func.id for node in ast.walk(tree) if isinstance(node, ast.Call)
                 and isinstance(node.func, ast.Name)}
        for forbidden in ("range", "linspace", "bisect"):
            self.assertNotIn(forbidden, calls)
        self.assertNotIn("refine", source.lower())

    def test_opposite_endpoint_still_used_when_first_cannot_clear(self):
        outcomes = iter([
            {"valid": False, "reason": "TARGET_CLEARANCE_FAILED", "positions_evaluated": 1},
            {"valid": True, "resulting_total_crossings": 0,
             "signed_delta_radians": 0.2, "displacement": 1.0, "positions_evaluated": 1},
            {"valid": False, "positions_evaluated": 1},
            {"valid": False, "positions_evaluated": 1},
        ])
        crossing = {"edge_a": ["A", "B"], "edge_b": ["C", "D"]}
        results = repair.evaluate_crossing_orientations(
            crossing, lambda **unused: next(outcomes))
        selected = repair.select_reducing_orientation(results, 1)
        self.assertEqual(("B", "A"), (selected["moving"], selected["pivot"]))


if __name__ == "__main__":
    unittest.main()
