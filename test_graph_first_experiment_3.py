import json
import math
import os
import subprocess
import sys
import tempfile
import unittest
from itertools import combinations
from pathlib import Path

import graph_first_experiment_2 as experiment_2
import graph_first_experiment_3 as experiment


ROOT = Path(__file__).parent
BASE = ROOT / "output" / "graph_first"
GRAPH = BASE / "iamp_graph.json"
ANALYSIS = BASE / "iamp_graph_analysis.json"
BLUEPRINT = BASE / "iamp_topology_blueprint.json"


class GraphFirstExperiment3Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.graph, _, cls.analysis, cls.blueprint = experiment.load_inputs(
            GRAPH, ANALYSIS, BLUEPRINT)
        cls.embedding, cls.reference, _ = experiment.fresh_reference(
            cls.graph, cls.analysis, cls.blueprint)
        cls.positions, _ = experiment.radial_positions(cls.reference, cls.blueprint)

    def test_center_origin_and_every_radius_equals_bfs_layer(self):
        center = self.blueprint["structural_center_selection"]["chosen_vertex"]
        self.assertEqual((0.0, 0.0), self.positions[center])
        layers = {item["id"]: item["bfs_layer"] for item in self.blueprint["vertices"]}
        for node in self.graph:
            self.assertAlmostEqual(layers[node], math.hypot(*self.positions[node]), places=12)

    def test_graph_identity_and_incidence_are_unchanged(self):
        graph_doc = json.loads(GRAPH.read_text())
        self.assertEqual(34, self.graph.number_of_nodes())
        self.assertEqual(44, self.graph.number_of_edges())
        expected = {frozenset((edge["source"], edge["target"])) for edge in graph_doc["edges"]}
        self.assertEqual(expected, {frozenset(edge) for edge in self.graph.edges})

    def test_only_three_authoritative_noncoordinate_inputs_are_declared(self):
        with tempfile.TemporaryDirectory() as directory:
            coordinates, _ = experiment.run(GRAPH, ANALYSIS, BLUEPRINT, Path(directory))
        self.assertEqual([GRAPH.as_posix(), ANALYSIS.as_posix(), BLUEPRINT.as_posix()],
                         coordinates["authoritative_inputs"])
        self.assertEqual([], coordinates["old_coordinate_artifacts_read"])

    def test_independent_geometric_failures_are_recomputed(self):
        validation = experiment.validate_geometry(self.graph, self.positions, self.embedding)
        edges = sorted((min(a, b), max(a, b)) for a, b in self.graph.edges)

        def orient(a, b, c):
            return ((b[0] - a[0]) * (c[1] - a[1]) -
                    (b[1] - a[1]) * (c[0] - a[0]))

        independent_crossings = []
        for (a, b), (c, d) in combinations(edges, 2):
            if {a, b} & {c, d}:
                continue
            values = (orient(self.positions[a], self.positions[b], self.positions[c]),
                      orient(self.positions[a], self.positions[b], self.positions[d]),
                      orient(self.positions[c], self.positions[d], self.positions[a]),
                      orient(self.positions[c], self.positions[d], self.positions[b]))
            if values[0] * values[1] < -1e-18 and values[2] * values[3] < -1e-18:
                independent_crossings.append({"edge_a": [a, b], "edge_b": [c, d]})
        independent_coincidences = [list(pair) for pair in combinations(sorted(self.graph), 2)
                                    if math.dist(self.positions[pair[0]], self.positions[pair[1]]) <= 1e-9]
        independent_on_edge = []
        for node in sorted(self.graph):
            for first, second in edges:
                if node in {first, second}:
                    continue
                point, a, b = self.positions[node], self.positions[first], self.positions[second]
                length2 = math.dist(a, b) ** 2
                if length2 == 0:
                    continue
                t = ((point[0] - a[0]) * (b[0] - a[0]) +
                     (point[1] - a[1]) * (b[1] - a[1])) / length2
                if (1e-9 < t < 1 - 1e-9 and
                        abs(orient(a, b, point)) <= 1e-9 * max(1.0, math.sqrt(length2))):
                    independent_on_edge.append({"vertex": node, "edge": [first, second]})
        self.assertEqual(independent_crossings, validation["different_edge_crossings"])
        self.assertEqual(independent_coincidences, validation["coincident_distinct_vertices"])
        self.assertEqual(independent_on_edge, validation["vertices_on_unrelated_edge_interiors"])

    def test_repeated_outputs_are_byte_identical(self):
        with tempfile.TemporaryDirectory() as first_dir, tempfile.TemporaryDirectory() as second_dir:
            experiment.run(GRAPH, ANALYSIS, BLUEPRINT, Path(first_dir))
            experiment.run(GRAPH, ANALYSIS, BLUEPRINT, Path(second_dir))
            for name in ("iamp_planar_reference.svg", "iamp_topology_radial_initial.svg",
                         "iamp_topology_radial_initial.json", "iamp_topology_radial_report.json"):
                self.assertEqual((Path(first_dir) / name).read_bytes(),
                                 (Path(second_dir) / name).read_bytes())

    def test_clean_process_hash_seed_does_not_change_outputs(self):
        with tempfile.TemporaryDirectory() as first_dir, tempfile.TemporaryDirectory() as second_dir:
            for seed, directory in (("1", first_dir), ("987654", second_dir)):
                environment = dict(os.environ, PYTHONHASHSEED=seed)
                subprocess.run([
                    sys.executable, str(ROOT / "graph_first_experiment_3.py"),
                    "--graph", str(GRAPH), "--analysis", str(ANALYSIS),
                    "--blueprint", str(BLUEPRINT), "--output", directory,
                ], cwd=ROOT, env=environment, check=True, capture_output=True, text=True)
            for name in ("iamp_planar_reference.svg", "iamp_topology_radial_initial.svg",
                         "iamp_topology_radial_initial.json", "iamp_topology_radial_report.json"):
                self.assertEqual((Path(first_dir) / name).read_bytes(),
                                 (Path(second_dir) / name).read_bytes())


if __name__ == "__main__":
    unittest.main()
