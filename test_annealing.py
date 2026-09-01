import unittest
from pathlib import Path

import graph_relax


INPUT = Path(__file__).parent / "test_data" / "iamp.circuit"


class AnnealingTests(unittest.TestCase):
    def test_move_accounting_and_best_energy(self):
        circuit = graph_relax.parse_circuit(INPUT)
        start = graph_relax.initial_positions(circuit)
        best, stats = graph_relax.anneal(circuit, start, iterations=20)
        self.assertEqual(20, stats["accepted_downhill_moves"] +
                         stats["accepted_uphill_moves"] + stats["rejected_moves"])
        self.assertLessEqual(graph_relax.evaluate(circuit, best)["objective"],
                             graph_relax.evaluate(circuit, start)["objective"])

    def test_annealing_is_deterministic(self):
        circuit = graph_relax.parse_circuit(INPUT)
        start = graph_relax.initial_positions(circuit)
        first = graph_relax.anneal(circuit, start, iterations=30)
        second = graph_relax.anneal(circuit, start, iterations=30)
        self.assertEqual(first, second)


if __name__ == "__main__":
    unittest.main()
