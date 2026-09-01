import tempfile
import unittest
from pathlib import Path

import networkx as nx

import experiment_4p as e4p
import experiment_4q as e4q
import graph_relax as gr


ROOT = Path(__file__).parent
INPUT = ROOT / "test_data" / "iamp.circuit"
OUTPUT = ROOT / "output"


class Experiment4QTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.circuit = gr.parse_circuit(INPUT)
        (cls.graph, cls.candidates, cls.rims, cls.spokes,
         cls.document) = e4q.verify_baseline(cls.circuit, OUTPUT)
        cls.k, cls.solutions, cls.levels = e4q.exact_cardinality_search(cls.graph, cls.candidates)

    def test_exact_4p_reconstruction_and_candidate_set(self):
        self.assertEqual((78, 108), (self.graph.number_of_nodes(), self.graph.number_of_edges()))
        self.assertFalse(nx.check_planarity(self.graph)[0])
        self.assertEqual(44, len(self.candidates))
        expected = {(terminal, net.name) for net in self.circuit.nets for terminal in net.terminals}
        actual = {(candidate["terminal"], candidate["net"]) for candidate in self.candidates}
        self.assertEqual(expected, actual)

    def test_only_incidence_removable_and_constraints_immutable(self):
        edge_roles = {tuple(sorted(edge)): data["role"] for *edge, data in self.graph.edges(data=True)}
        self.assertTrue(all(edge_roles[tuple(sorted(e4q.edge_tuple(candidate)))] == "electrical_incidence"
                            for candidate in self.candidates))
        for solution in self.solutions:
            reduced = e4q.graph_without(self.graph, solution)
            self.assertTrue(all(reduced.has_edge(*edge) for edge in self.rims | self.spokes))

    def test_cardinality_order_exhaustive_minimum(self):
        self.assertEqual(2, self.k)
        self.assertEqual([0, 1, 2], [level["k"] for level in self.levels])
        self.assertTrue(all(level["tested_combinations"] == level["estimated_combinations"]
                            for level in self.levels))
        self.assertEqual([0, 0, 4], [level["planar_combinations"] for level in self.levels])

    def test_all_minimum_sets_verify_and_preserve_other_incidence(self):
        original = e4q.original_certificate_incidence_edges(OUTPUT)
        for index, solution in enumerate(self.solutions, 1):
            graph, result, checks = e4q.verify_solution(
                self.circuit, self.graph, self.candidates, solution, self.rims, self.spokes)
            self.assertTrue(result["verification"]["verified"])
            self.assertTrue(all(checks.values()))
            self.assertEqual(106, graph.number_of_edges())
            self.assertTrue({tuple(sorted(e4q.edge_tuple(candidate))) for candidate in solution} & original)

    def test_frequency_and_mandatory_results(self):
        frequency = e4q.frequency_analysis(self.candidates, self.solutions)
        self.assertEqual(["U1.8 -- +9V"], frequency["mandatory_incidences"])
        self.assertEqual(["U1"], frequency["mandatory_components"])
        self.assertEqual(["+9V"], frequency["mandatory_nets"])

    def test_repeated_run_determinism(self):
        first_dir = Path(tempfile.mkdtemp(prefix="graph-relax-4q-test-"))
        second_dir = Path(tempfile.mkdtemp(prefix="graph-relax-4q-test-"))
        for directory in (first_dir, second_dir):
            for name in ("planarity_report.json", "planarity_C_nonplanarity_certificate.json"):
                (directory / name).write_bytes((OUTPUT / name).read_bytes())
        self.assertEqual(e4q.run(INPUT, first_dir), e4q.run(INPUT, second_dir))


if __name__ == "__main__":
    unittest.main()
