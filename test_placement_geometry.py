import math
import unittest

import placement_geometry as geometry


class PlacementGeometryTests(unittest.TestCase):
    def test_two_close_nodes_scale_to_exact_minimum(self):
        result, report = geometry.enforce_minimum_node_spacing(
            {"a": (0, 0), "b": (0.25, 0)}, 2.0)
        self.assertAlmostEqual(2.0, math.dist(result["a"], result["b"]), places=12)
        self.assertAlmostEqual(8.0, report["applied_uniform_scale_factor"], places=12)

    def test_already_spaced_graph_is_unchanged(self):
        original = {"a": (1.0, 2.0), "b": (4.0, 6.0)}
        result, report = geometry.enforce_minimum_node_spacing(original, 2.0)
        self.assertEqual(original, result)
        self.assertEqual(1.0, report["applied_uniform_scale_factor"])

    def test_scaling_preserves_ratios_angles_and_crossings(self):
        original = {"a": (0, 0), "b": (0.2, 0), "c": (0.2, 0.2), "d": (0, 0.2)}
        edges = [("a", "c"), ("b", "d")]
        before = geometry.graph_geometry_diagnostics(original, edges)
        result, _ = geometry.enforce_minimum_node_spacing(original, 1.0)
        after = geometry.graph_geometry_diagnostics(result, edges)
        self.assertAlmostEqual(math.dist(original["a"], original["c"]) /
                               math.dist(original["a"], original["b"]),
                               math.dist(result["a"], result["c"]) /
                               math.dist(result["a"], result["b"]), places=12)
        def angle(points):
            first = (points["a"][0] - points["b"][0], points["a"][1] - points["b"][1])
            second = (points["c"][0] - points["b"][0], points["c"][1] - points["b"][1])
            return math.acos(sum(x*y for x, y in zip(first, second)) /
                             (math.hypot(*first) * math.hypot(*second)))
        self.assertAlmostEqual(angle(original), angle(result), places=12)
        self.assertEqual(before["proper_unrelated_edge_crossings"],
                         after["proper_unrelated_edge_crossings"])

    def test_nonzero_anchor(self):
        result, report = geometry.enforce_minimum_node_spacing(
            {"a": (10, 10), "b": (10.5, 10)}, 1.0, anchor=(10, 10))
        self.assertEqual((10.0, 10.0), result["a"])
        self.assertEqual((11.0, 10.0), result["b"])
        self.assertEqual([10.0, 10.0], report["anchor"])

    def test_single_vertex_unchanged(self):
        original = {"only": (3, 7)}
        result, report = geometry.enforce_minimum_node_spacing(original)
        self.assertEqual(original, result)
        self.assertEqual(1.0, report["applied_uniform_scale_factor"])
        self.assertIsNone(report["final_minimum_node_distance"])

    def test_coincident_vertices_fail_without_division(self):
        original = {"a": (1, 1), "b": (1, 1)}
        result, report = geometry.enforce_minimum_node_spacing(original)
        self.assertEqual(original, result)
        self.assertEqual("FAILED_COINCIDENT_VERTICES", report["status"])
        self.assertIsNone(report["applied_uniform_scale_factor"])
        self.assertEqual([["a", "b"]], report["coincident_vertex_pairs"])

    def test_generic_coincidence_resolution_is_deterministic(self):
        original = {"b": (2, 3), "a": (2, 3), "other": (9, 9)}
        first, report = geometry.resolve_coincident_vertices(original, 2.0)
        second, _ = geometry.resolve_coincident_vertices(original, 2.0)
        self.assertEqual(first, second)
        self.assertEqual("RESOLVED", report["status"])
        self.assertAlmostEqual(2.0, math.dist(first["a"], first["b"]), places=12)
        self.assertEqual(0, report["remaining_coincident_vertex_pair_count"])

    def test_vertex_on_edge_is_independent_of_spacing(self):
        points = {"a": (0, 0), "b": (4, 0), "c": (2, 0)}
        diagnostics = geometry.graph_geometry_diagnostics(points, [("a", "b")])
        self.assertGreaterEqual(diagnostics["minimum_node_distance"], 1.0)
        self.assertEqual(0, diagnostics["coincident_vertex_pair_count"])
        self.assertEqual([{"vertex": "c", "edge": ["a", "b"]}],
                         diagnostics["vertices_on_unrelated_edge_interiors"])


if __name__ == "__main__":
    unittest.main()
