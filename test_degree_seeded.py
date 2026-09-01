import unittest
from pathlib import Path

import graph_relax


INPUT = Path(__file__).parent / "test_data" / "iamp.circuit"


class DegreeSeededTests(unittest.TestCase):
    def setUp(self):
        self.circuit = graph_relax.parse_circuit(INPUT)
        self.degrees = graph_relax.component_degrees(self.circuit)

    def test_degree_calculation_and_singletons_excluded(self):
        expected = {ref: set() for ref in self.circuit.refs}
        for net in self.circuit.nets:
            if len(net.terminals) > 1:
                for terminal in net.terminals:
                    expected[terminal.rsplit(".", 1)[0]].add(net.name)
        self.assertEqual({ref: len(nets) for ref, nets in expected.items()}, self.degrees)
        self.assertEqual(0, self.degrees["J_RA"] - 2)  # RA_W is a singleton.

    def test_deterministic_degree_sorting(self):
        ordered = sorted(self.circuit.refs, key=lambda ref: (-self.degrees[ref], ref))
        self.assertEqual(ordered, sorted(ordered, key=lambda ref: (-self.degrees[ref], ref)))

    def test_highest_degree_is_innermost_and_no_overlap(self):
        positions = graph_relax.degree_seeded_positions(self.circuit)
        radii = {ref: graph_relax.math.hypot(*point) for ref, point in positions.items()}
        for high in self.circuit.refs:
            for low in self.circuit.refs:
                if self.degrees[high] > self.degrees[low]:
                    self.assertLess(radii[high], radii[low])
        self.assertEqual(0, graph_relax.evaluate(self.circuit, positions)["component_overlaps"])

    def test_tied_maximum_uses_symmetric_ring(self):
        circuit = graph_relax.Circuit(
            "tie", ("A", "B"), {"A": ("1",), "B": ("1",)},
            (graph_relax.Net("N", ("A.1", "B.1")),))
        positions = graph_relax.degree_seeded_positions(circuit)
        self.assertNotEqual((0.0, 0.0), positions["A"])
        self.assertAlmostEqual(graph_relax.math.hypot(*positions["A"]),
                               graph_relax.math.hypot(*positions["B"]))
        self.assertAlmostEqual(positions["A"][0], -positions["B"][0], places=9)
        self.assertAlmostEqual(positions["A"][1], -positions["B"][1], places=9)

    def test_short_repeated_run_determinism(self):
        start = graph_relax.degree_seeded_positions(self.circuit)
        first = graph_relax.anneal(self.circuit, start, iterations=30)
        second = graph_relax.anneal(self.circuit, start, iterations=30)
        self.assertEqual(first, second)


if __name__ == "__main__":
    unittest.main()
