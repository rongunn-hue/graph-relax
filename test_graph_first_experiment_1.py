import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import networkx as nx

import graph_first_experiment_1 as experiment


ROOT = Path(__file__).parent
INPUT = ROOT / "circuit_examples" / "iamp.circuit"


class GraphFirstExperiment1Tests(unittest.TestCase):
    def setUp(self):
        self.circuit = experiment.parse_circuit(INPUT)
        self.graph = experiment.build_graph(self.circuit)

    def test_parser_is_deterministic_and_preserves_source_incidences(self):
        again = experiment.parse_circuit(INPUT)
        self.assertEqual(self.circuit, again)
        parsed = {(item.component, item.net, item.terminal) for item in self.circuit.incidences}
        graph = {(data["component"], data["net"], data["terminal"])
                 for _, _, data in self.graph.edges(data=True)}
        self.assertEqual(len(self.circuit.incidences), len(parsed))
        self.assertEqual(parsed, graph)

    def test_counts_are_derived_from_parsed_source(self):
        self.assertEqual(len(self.circuit.components) + len(self.circuit.nets),
                         self.graph.number_of_nodes())
        self.assertEqual(len(self.circuit.incidences), self.graph.number_of_edges())

    def test_graph_is_bipartite_with_declared_partitions(self):
        self.assertTrue(nx.is_bipartite(self.graph))
        for first, second in self.graph.edges:
            self.assertNotEqual(self.graph.nodes[first]["type"], self.graph.nodes[second]["type"])

    def test_planarity_is_reproducible_and_embedding_is_verified(self):
        first = experiment.analyze(self.circuit, self.graph)
        second = experiment.analyze(self.circuit, experiment.build_graph(self.circuit))
        self.assertEqual(first["planarity"], second["planarity"])
        self.assertTrue(first["planarity"]["embedding_validation"]["valid"])

    def test_serialization_is_byte_identical(self):
        with tempfile.TemporaryDirectory() as first_dir, tempfile.TemporaryDirectory() as second_dir:
            experiment.run(INPUT, Path(first_dir))
            experiment.run(INPUT, Path(second_dir))
            for name in ("iamp_graph.json", "iamp_graph_analysis.json"):
                first = (Path(first_dir) / name).read_bytes()
                second = (Path(second_dir) / name).read_bytes()
                self.assertEqual(first, second)
                self.assertEqual(hashlib.sha256(first).hexdigest(), hashlib.sha256(second).hexdigest())
                json.loads(first)


if __name__ == "__main__":
    unittest.main()
