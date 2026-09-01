import unittest
from pathlib import Path

import graph_relax


INPUT = Path(__file__).parent / "test_data" / "iamp.circuit"


class JunctionAnnealingTests(unittest.TestCase):
    def setUp(self):
        self.circuit = graph_relax.parse_circuit(INPUT)
        self.positions = graph_relax.initial_positions(self.circuit)
        self.junctions = graph_relax.initial_junctions(self.circuit, self.positions)
        terminals = graph_relax.terminal_positions(self.circuit, self.positions)
        self.segments = graph_relax.junction_segments(self.circuit, terminals, self.junctions)

    def test_movable_junction_creation(self):
        expected = {n.name for n in self.circuit.nets if len(n.terminals) >= 3}
        self.assertEqual(expected, set(self.junctions))
        self.assertEqual(7, len(self.junctions))

    def test_two_terminal_nets_are_direct(self):
        two_terminal = {n.name for n in self.circuit.nets if len(n.terminals) == 2}
        counts = {name: sum(s[0] == name for s in self.segments) for name in two_terminal}
        self.assertTrue(all(count == 1 for count in counts.values()))

    def test_singletons_have_no_segments(self):
        singleton = {n.name for n in self.circuit.nets if len(n.terminals) == 1}
        self.assertTrue(singleton.isdisjoint(s[0] for s in self.segments))

    def test_junction_annealing_repeated_run_determinism(self):
        first = graph_relax.anneal_junction_state(
            self.circuit, self.positions, self.junctions, iterations=40)
        second = graph_relax.anneal_junction_state(
            self.circuit, self.positions, self.junctions, iterations=40)
        self.assertEqual(first, second)
        movable = graph_relax.junction_movable_objects(self.circuit, self.junctions)
        self.assertEqual(len(self.circuit.refs) + len(self.junctions), len(movable))


if __name__ == "__main__":
    unittest.main()
