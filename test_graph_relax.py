import tempfile
import unittest
from pathlib import Path

import graph_relax


INPUT = Path(__file__).parent / "test_data" / "iamp.circuit"


class GraphRelaxTests(unittest.TestCase):
    def test_iamp_counts(self):
        circuit = graph_relax.parse_circuit(INPUT)
        self.assertEqual(19, len(circuit.refs))
        self.assertEqual(15, len(circuit.nets))
        self.assertEqual(["R1_W", "R4_W", "RA_W"],
                         [n.name for n in circuit.nets if len(n.terminals) == 1])

    def test_centroid_star_segment_count(self):
        circuit = graph_relax.parse_circuit(INPUT)
        terminals = graph_relax.terminal_positions(circuit, graph_relax.initial_positions(circuit))
        expected = sum(len(n.terminals) for n in circuit.nets if len(n.terminals) > 1)
        self.assertEqual(expected, len(graph_relax.net_segments(circuit, terminals)))

    def test_relaxation_is_deterministic(self):
        circuit = graph_relax.parse_circuit(INPUT)
        start = graph_relax.initial_positions(circuit)
        a, sa, ra = graph_relax.relax(circuit, start)
        b, sb, rb = graph_relax.relax(circuit, start)
        self.assertEqual((a, sa, ra), (b, sb, rb))


if __name__ == "__main__":
    unittest.main()
