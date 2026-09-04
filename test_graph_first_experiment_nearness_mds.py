import ast
import json
import math
import tempfile
import unittest
from itertools import combinations
from pathlib import Path
from unittest import mock

import networkx as nx
import numpy as np

import graph_first_experiment_nearness_mds as experiment


ROOT = Path(__file__).parent
GRAPH = ROOT / "output" / "graph_first" / "iamp_graph.json"


class NearnessMdsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.graph, _ = experiment.load_graph(GRAPH)
        cls.nodes = sorted(cls.graph)
        cls.distances = experiment.distance_matrix(cls.graph, cls.nodes)
        raw, cls.values, cls.gram, _, _ = experiment.classical_mds(cls.distances, cls.nodes)
        cls.center, cls.closeness = experiment.choose_center(cls.graph)
        cls.centered = experiment.center_coordinates(raw, cls.nodes, cls.center)
        cls.positions = {node: tuple(cls.centered[index]) for index, node in enumerate(cls.nodes)}

    def test_authoritative_graph_counts_and_connectivity(self):
        self.assertEqual(34, self.graph.number_of_nodes())
        self.assertEqual(44, self.graph.number_of_edges())
        self.assertTrue(nx.is_connected(self.graph))

    def test_distance_matrix_properties(self):
        self.assertTrue(np.array_equal(self.distances, self.distances.T))
        self.assertTrue(np.array_equal(np.diag(self.distances), np.zeros(len(self.nodes))))

    def test_center_is_derived_and_exactly_at_origin(self):
        degree = dict(self.graph.degree())
        expected = min(self.graph,
                       key=lambda node: (-self.closeness[node], -degree[node], node))
        self.assertEqual(expected, self.center)
        self.assertEqual((0.0, 0.0), self.positions[self.center])

    def test_coordinates_recompute_from_distance_matrix_only(self):
        count = len(self.nodes)
        j = np.eye(count) - np.ones((count, count)) / count
        expected_gram = -0.5 * j @ (self.distances ** 2) @ j
        self.assertTrue(np.allclose(expected_gram, self.gram, atol=1e-12))
        reconstructed, _, _, _, _ = experiment.classical_mds(self.distances.copy(), self.nodes)
        reconstructed = experiment.center_coordinates(reconstructed, self.nodes, self.center)
        self.assertTrue(np.allclose(reconstructed, self.centered, atol=1e-12))

    def test_source_has_no_forbidden_layout_calls(self):
        tree = ast.parse((ROOT / "graph_first_experiment_nearness_mds.py").read_text())
        called = {node.func.attr if isinstance(node.func, ast.Attribute) else node.func.id
                  for node in ast.walk(tree) if isinstance(node, ast.Call)
                  and isinstance(node.func, (ast.Attribute, ast.Name))}
        self.assertTrue({"random_layout", "spring_layout", "forceatlas2_layout",
                         "combinatorial_embedding_to_pos"}.isdisjoint(called))

    def test_crossings_coincidences_and_vertex_edge_events_independently(self):
        reported = experiment.validate(self.graph, self.positions)
        edges = sorted((min(a, b), max(a, b)) for a, b in self.graph.edges)
        crossings = []
        for (a, b), (c, d) in combinations(edges, 2):
            if {a, b} & {c, d}:
                continue
            values = (experiment.cross(self.positions[a], self.positions[b], self.positions[c]),
                      experiment.cross(self.positions[a], self.positions[b], self.positions[d]),
                      experiment.cross(self.positions[c], self.positions[d], self.positions[a]),
                      experiment.cross(self.positions[c], self.positions[d], self.positions[b]))
            epsilon = 1e-9 * max(1.0, math.dist(self.positions[a], self.positions[b]),
                                 math.dist(self.positions[c], self.positions[d]))
            first_straddles = ((values[0] > epsilon and values[1] < -epsilon) or
                               (values[0] < -epsilon and values[1] > epsilon))
            second_straddles = ((values[2] > epsilon and values[3] < -epsilon) or
                                (values[2] < -epsilon and values[3] > epsilon))
            if first_straddles and second_straddles:
                crossings.append({"edge_a": [a, b], "edge_b": [c, d]})
        coincidences = [list(pair) for pair in combinations(self.nodes, 2)
                        if math.dist(self.positions[pair[0]], self.positions[pair[1]]) <= 1e-9]
        on_edges = []
        for node in self.nodes:
            for a, b in edges:
                if node not in {a, b} and experiment.point_on_segment_interior(
                        self.positions[node], self.positions[a], self.positions[b]):
                    on_edges.append({"vertex": node, "edge": [a, b]})
        self.assertEqual(crossings, reported["proper_unrelated_edge_crossings"])
        self.assertEqual(coincidences, reported["coincident_distinct_vertices"])
        self.assertEqual(on_edges, reported["vertices_on_unrelated_edge_interiors"])

    def test_repeated_serialization_is_identical_and_declares_no_old_coordinates(self):
        with tempfile.TemporaryDirectory() as first_dir, tempfile.TemporaryDirectory() as second_dir:
            first, _ = experiment.run(GRAPH, Path(first_dir))
            second, _ = experiment.run(GRAPH, Path(second_dir))
            self.assertEqual([], first["old_coordinate_artifacts_read"])
            self.assertEqual([], second["old_coordinate_artifacts_read"])
            for name in ("iamp_nearness_mds.svg", "iamp_nearness_mds.json",
                         "iamp_nearness_mds_report.json"):
                self.assertEqual((Path(first_dir) / name).read_bytes(),
                                 (Path(second_dir) / name).read_bytes())

    def test_mds_pipeline_calls_shared_spacing_enforcement(self):
        with tempfile.TemporaryDirectory() as directory, mock.patch.object(
                experiment.geometry, "enforce_minimum_node_spacing",
                wraps=experiment.geometry.enforce_minimum_node_spacing) as shared:
            _, report = experiment.run(GRAPH, Path(directory))
        shared.assert_called_once()
        self.assertEqual(experiment.geometry.MIN_NODE_DISTANCE,
                         shared.call_args.args[1])
        self.assertEqual("SATISFIED",
                         report["minimum_spacing_enforcement"]["status"])
        self.assertEqual("RESOLVED", report["coincidence_resolution"]["status"])


if __name__ == "__main__":
    unittest.main()
