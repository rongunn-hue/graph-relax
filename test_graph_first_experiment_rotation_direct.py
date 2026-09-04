import ast
import math
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

import graph_first_experiment_endpoint_rotation as sweep
import graph_first_experiment_rotation_direct as direct
import placement_geometry as geometry
from direct_crossing_repair import endpoint_pivot_orientations as generic_orientations


ROOT = Path(__file__).parent
BASE = ROOT/"output"/"graph_first"


class DirectRotationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.graph, cls.start, _, _ = sweep.load_start(
            BASE/"iamp_graph.json", BASE/"iamp_nearness_mds.json",
            BASE/"iamp_nearness_mds_report.json")

    def test_segment_disk_clip_is_analytic_and_correct(self):
        clipped = direct.segment_disk_clip((-3, 0), (3, 0), (0, 0), 2)
        self.assertEqual([(-2.0, 0.0), (2.0, 0.0)], clipped)

    def test_direct_cardinality_and_monotonicity(self):
        _, diagnostics, moves, _, history, counters = direct.direct_search(self.graph, self.start)
        self.assertLessEqual(max((item["positions_evaluated"]
                                  for item in counters["positions_per_choice"]), default=0), 1)
        self.assertLessEqual(max(item["calculated"] for item in counters["solutions_per_crossing"]), 4)
        counts = [len(items) for items in history]
        self.assertTrue(all(after < before for before, after in zip(counts, counts[1:])))
        self.assertEqual(counts[-1], diagnostics["proper_unrelated_edge_crossing_count"])
        for move in moves:
            self.assertAlmostEqual(move["radius"], move["preserved_edge_length_before"], places=12)

    def test_only_opposite_endpoints_supply_pivots(self):
        crossing = {"edge_a": ["A", "B"], "edge_b": ["C", "D"]}
        self.assertEqual([
            {"moving": "A", "pivot": "B", "crossed_edge": ("C", "D")},
            {"moving": "B", "pivot": "A", "crossed_edge": ("C", "D")},
            {"moving": "C", "pivot": "D", "crossed_edge": ("A", "B")},
            {"moving": "D", "pivot": "C", "crossed_edge": ("A", "B")},
        ], generic_orientations(crossing))

    def test_active_uncrosser_has_no_150_spacing_branch(self):
        source = (ROOT/"graph_first_experiment_rotation_direct.py").read_text()
        self.assertNotIn("150", source)
        self.assertNotIn("spacing_forbidden", source)
        self.assertNotIn("minimum_spacing", source)

    def test_source_contains_no_angular_search_construct(self):
        source = (ROOT/"graph_first_experiment_rotation_direct.py").read_text()
        tree = ast.parse(source)
        names = {node.func.id for node in ast.walk(tree) if isinstance(node, ast.Call)
                 and isinstance(node.func, ast.Name)}
        self.assertNotIn("linspace", names)
        self.assertNotIn("range", names)
        self.assertNotIn("bisect", source.lower())
        self.assertNotIn("refine", source.lower())

    def test_direct_search_is_deterministic(self):
        first = direct.direct_search(self.graph, self.start)
        second = direct.direct_search(self.graph, self.start)
        self.assertEqual(first[0], second[0])
        self.assertEqual(first[2], second[2])
        self.assertEqual(first[4], second[4])

    def test_multiedge_vertex_is_planted_and_all_edges_reconnected(self):
        coordinates = {"M": (0.0, 0.0), "A": (2.0, 0.0),
                       "B": (0.0, 2.0), "X": (5.0, 5.0), "Y": (6.0, 5.0)}
        edges = [("M", "A"), ("M", "B"), ("X", "Y")]
        before = geometry.rebuild_straight_edge_segments(coordinates, edges)
        planted = dict(coordinates)
        planted["M"] = (1.0, 1.0)
        after = geometry.rebuild_straight_edge_segments(planted, edges)
        report = geometry.validate_rebuilt_edge_segments(planted, edges, after)
        self.assertTrue(report["valid"])
        before_by_edge = {tuple(item["edge"]): item for item in before}
        after_by_edge = {tuple(item["edge"]): item for item in after}
        for edge in [("A", "M"), ("B", "M")]:
            self.assertNotEqual(before_by_edge[edge], after_by_edge[edge])
            self.assertEqual((1.0, 1.0), after_by_edge[edge]["segment_end"])
        self.assertEqual(before_by_edge[("X", "Y")], after_by_edge[("X", "Y")])
        expected = sorted((min(first, second), max(first, second)) for first, second in edges)
        self.assertEqual(expected, [tuple(item["edge"]) for item in after])

    def test_motion_stages_render_post_move_reconnected_graph(self):
        final, _, moves, states, _, _ = direct.direct_search(self.graph, self.start)
        svg = direct.direct_motion_svg(self.graph, states, moves, sweep.viewbox_for(states))
        self.assertIn("POST-MOVE graph", svg)
        self.assertIn('data-active-vertex="net:GND"', svg)
        for move, state in zip(moves, states[1:]):
            self.assertEqual(tuple(move["end_position"]), state[move["moving_vertex"]])
            self.assertTrue(move["all_incident_edges_regenerated"])
        serialized_counts = []
        for state in states:
            serialized = {node: (float(direct.svg_number(point[0])),
                                 float(direct.svg_number(point[1])))
                          for node, point in state.items()}
            serialized_counts.append(
                geometry.graph_geometry_diagnostics(
                    serialized, self.graph.edges)["proper_unrelated_edge_crossing_count"])
        self.assertEqual([6, 3, 2, 1, 0], serialized_counts)
        self.assertEqual(final, states[-1])

    def test_motion_stops_on_clean_persistent_final_graph(self):
        final, _, moves, states, _, _ = direct.direct_search(self.graph, self.start)
        svg = direct.direct_motion_svg(self.graph, states, moves, sweep.viewbox_for(states))
        self.assertNotIn('repeatCount="indefinite"', svg)
        root = ET.fromstring(svg)
        ns = {"svg": "http://www.w3.org/2000/svg"}
        group = root.find("svg:g[@id='final-persistent-graph']", ns)
        self.assertIsNotNone(group)
        final_vertices = {item.attrib["data-final-vertex"]:
                          (float(item.attrib["cx"]), float(item.attrib["cy"]))
                          for item in group.findall("svg:circle", ns)}
        final_edges = {(item.attrib["data-final-edge-u"], item.attrib["data-final-edge-v"]):
                       ((float(item.attrib["x1"]), float(item.attrib["y1"])),
                        (float(item.attrib["x2"]), float(item.attrib["y2"])))
                       for item in group.findall("svg:line", ns)}
        self.assertEqual(34, len(final_vertices))
        self.assertEqual(44, len(final_edges))
        for node, point in final.items():
            self.assertAlmostEqual(point[0], final_vertices[node][0], places=6)
            self.assertAlmostEqual(point[1], final_vertices[node][1], places=6)
        for record in direct.rendered_edge_records(self.graph, final):
            key = tuple(record["edge"])
            self.assertAlmostEqual(record["segment_start"][0], final_edges[key][0][0], places=6)
            self.assertAlmostEqual(record["segment_start"][1], final_edges[key][0][1], places=6)
            self.assertAlmostEqual(record["segment_end"][0], final_edges[key][1][0], places=6)
            self.assertAlmostEqual(record["segment_end"][1], final_edges[key][1][1], places=6)


if __name__ == "__main__":
    unittest.main()
