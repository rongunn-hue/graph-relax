import math
import unittest
from pathlib import Path

import experiment_4c as e4c
import graph_relax as gr


INPUT = Path(__file__).parent / "test_data" / "iamp.circuit"


class Experiment4CTests(unittest.TestCase):
    def setUp(self):
        self.circuit = gr.parse_circuit(INPUT)
        self.adjacency, self.distances, self.table, self.root, self.layers, self.connected = e4c.graph_analysis(self.circuit)

    def test_adjacency_multiterminal_and_singleton_exclusion(self):
        gnd_refs = {"C1", "C2", "J_EA", "J_EB", "J_OUT", "RMIN1", "J_PWR"}
        for first in gnd_refs:
            self.assertTrue(gnd_refs - {first} <= set(self.adjacency[first]))
        self.assertNotIn("J_R1", self.adjacency["J_R1"])
        self.assertEqual(tuple(sorted(set(self.adjacency["J_R1"]))), self.adjacency["J_R1"])

    def test_shortest_paths_closeness_and_root(self):
        self.assertTrue(self.connected)
        self.assertEqual(0, self.distances[self.root][self.root])
        for ref, row in self.table.items():
            expected = (len(self.circuit.refs) - 1) / row["distance_sum"]
            self.assertAlmostEqual(expected, row["closeness"])
        ranked = min(self.circuit.refs, key=lambda ref: (-self.table[ref]["closeness"],
                                                          self.table[ref]["distance_sum"],
                                                          -self.table[ref]["degree"], ref))
        self.assertEqual(ranked, self.root)

    def test_layers_and_deterministic_initialization(self):
        first, first_analysis = e4c.closeness_initial_positions(self.circuit)
        second, second_analysis = e4c.closeness_initial_positions(self.circuit)
        self.assertEqual((first, first_analysis), (second, second_analysis))
        self.assertEqual((0.0, 0.0), first[self.root])
        radii = {ref: math.hypot(*point) for ref, point in first.items()}
        for near in self.circuit.refs:
            for far in self.circuit.refs:
                if self.layers[near] < self.layers[far]:
                    self.assertLess(radii[near], radii[far])
        self.assertEqual(0, gr.evaluate(self.circuit, first)["component_overlaps"])

    def test_release_has_no_constraints(self):
        positions, _ = e4c.closeness_initial_positions(self.circuit)
        released = e4c.release_initialization(positions)
        self.assertEqual(positions, released)
        self.assertIsNot(positions, released)
        released[self.root] = (10.0, 20.0)
        self.assertEqual((10.0, 20.0), released[self.root])

    def test_short_repeated_run_determinism(self):
        positions, _ = e4c.closeness_initial_positions(self.circuit)
        first = gr.anneal(self.circuit, e4c.release_initialization(positions), iterations=30)
        second = gr.anneal(self.circuit, e4c.release_initialization(positions), iterations=30)
        self.assertEqual(first, second)


if __name__ == "__main__":
    unittest.main()
