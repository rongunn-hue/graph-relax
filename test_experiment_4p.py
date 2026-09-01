import tempfile
import unittest
from pathlib import Path

import networkx as nx

import experiment_4p as e4p
import graph_relax as gr


ROOT = Path(__file__).parent
INPUT = ROOT / "test_data" / "iamp.circuit"


class Experiment4PTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.circuit = gr.parse_circuit(INPUT)
        cls.a = e4p.incidence_graph(cls.circuit)
        cls.b = e4p.terminal_expanded_graph(cls.circuit)

    def test_exact_incidence_graph(self):
        self.assertEqual((34, 44), (self.a.number_of_nodes(), self.a.number_of_edges()))
        expected = {(e4p.component_vertex(t.rsplit('.', 1)[0]), e4p.net_vertex(n.name))
                    for n in self.circuit.nets for t in n.terminals}
        actual = {tuple(sorted(edge)) for edge in self.a.edges()}
        self.assertEqual({tuple(sorted(edge)) for edge in expected}, actual)

    def test_terminal_expansion_preserves_incidence(self):
        self.assertEqual((78, 88), (self.b.number_of_nodes(), self.b.number_of_edges()))
        self.assertTrue(e4p.incidence_preserved(self.circuit, self.b))
        self.assertEqual(44, sum(data["kind"] == "terminal" for _, data in self.b.nodes(data=True)))

    def test_ordinary_planarity_and_embedding_verification(self):
        for graph in (self.a, self.b):
            result = e4p.exact_planarity(graph)
            self.assertEqual("PLANAR", result["result"])
            self.assertTrue(result["verification"]["verified"])
            self.assertEqual(2, result["verification"]["euler_characteristic"])

    def test_fixed_orders_are_geometry_defined_and_complete(self):
        orders = e4p.cyclic_terminal_orders(self.circuit)
        self.assertEqual(set(e4p.component_vertex(ref) for ref in self.circuit.refs), set(orders))
        for ref in self.circuit.refs:
            self.assertEqual([e4p.terminal_vertex(f"{ref}.{pin}") for pin in self.circuit.pins[ref]],
                             orders[e4p.component_vertex(ref)])

    def test_wheel_relaxation_and_certificate_are_nonplanar(self):
        wheel = e4p.terminal_order_wheel_graph(self.circuit)
        self.assertEqual((78, 108), (wheel.number_of_nodes(), wheel.number_of_edges()))
        result = e4p.exact_planarity(wheel)
        self.assertEqual("NONPLANAR", result["result"])
        self.assertTrue(result["verification"]["certificate_independently_retested_nonplanar"])
        components, nets = e4p.certificate_participants(result["certificate"])
        self.assertTrue({"J_PWR", "R3", "U1"} <= set(components))
        self.assertTrue({"+9V", "-9V", "NODE_B", "NODE_C"} <= set(nets))

    def test_repeated_artifact_determinism(self):
        first_dir = Path(tempfile.mkdtemp(prefix="graph-relax-4p-test-"))
        second_dir = Path(tempfile.mkdtemp(prefix="graph-relax-4p-test-"))
        first = e4p.run(INPUT, first_dir)
        second = e4p.run(INPUT, second_dir)
        self.assertEqual(first, second)


if __name__ == "__main__":
    unittest.main()
