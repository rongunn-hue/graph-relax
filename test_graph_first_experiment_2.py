import tempfile
import unittest
from pathlib import Path

import networkx as nx

import graph_first_experiment_2 as experiment


ROOT = Path(__file__).parent
GRAPH = ROOT / "output" / "graph_first" / "iamp_graph.json"
ANALYSIS = ROOT / "output" / "graph_first" / "iamp_graph_analysis.json"


class GraphFirstExperiment2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.graph, _, cls.analysis = experiment.load_authoritative_graph(GRAPH, ANALYSIS)
        cls.blueprint = experiment.build_blueprint(GRAPH, ANALYSIS)

    def test_structural_center_is_derived_by_declared_rule(self):
        degree = dict(self.graph.degree())
        closeness = nx.closeness_centrality(self.graph)
        expected = min(self.graph, key=lambda node: (-closeness[node], -degree[node], node))
        self.assertEqual(expected,
                         self.blueprint["structural_center_selection"]["chosen_vertex"])

    def test_bfs_layers_are_exact_and_cover_every_vertex(self):
        center = self.blueprint["structural_center_selection"]["chosen_vertex"]
        expected = nx.single_source_shortest_path_length(self.graph, center)
        actual = {node: int(layer) for layer, nodes in self.blueprint["bfs_layers"].items()
                  for node in nodes}
        self.assertEqual(expected, actual)
        self.assertEqual(set(self.graph), {item["id"] for item in self.blueprint["vertices"]})

    def test_block_decomposition_covers_edges_and_matches_cut_structure(self):
        blocks = self.blueprint["block_structure"]
        covered = {frozenset(edge) for block in blocks["blocks"] for edge in block["edges"]}
        self.assertEqual({frozenset(edge) for edge in self.graph.edges}, covered)
        self.assertEqual(sorted(nx.articulation_points(self.graph)), blocks["articulation_vertices"])
        self.assertEqual(sorted([sorted(edge) for edge in nx.bridges(self.graph)]), blocks["bridges"])

    def test_every_embedding_face_is_enumerated_once_and_uses_embedding_edges(self):
        embedding = experiment.restore_embedding(self.graph, self.analysis)
        faces = self.blueprint["planar_embedding"]["faces"]
        self.assertEqual(self.graph.number_of_edges() - self.graph.number_of_nodes() + 2,
                         len(faces))
        directed = []
        for face in faces:
            walk = face["boundary_walk"]
            directed.extend((walk[index], walk[(index + 1) % len(walk)])
                            for index in range(len(walk)))
        expected = {(node, neighbor) for node in embedding
                    for neighbor in embedding.neighbors_cw_order(node)}
        self.assertEqual(expected, set(directed))
        self.assertEqual(len(expected), len(directed))

    def test_serialization_is_byte_identical(self):
        with tempfile.TemporaryDirectory() as directory:
            first = Path(directory) / "first.json"
            second = Path(directory) / "second.json"
            experiment.run(GRAPH, ANALYSIS, first)
            experiment.run(GRAPH, ANALYSIS, second)
            self.assertEqual(first.read_bytes(), second.read_bytes())


if __name__ == "__main__":
    unittest.main()
