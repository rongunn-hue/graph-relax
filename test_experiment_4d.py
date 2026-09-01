import unittest
from pathlib import Path

import experiment_4c as e4c
import experiment_4d as e4d
import graph_relax as gr


INPUT = Path(__file__).parent / "test_data" / "iamp.circuit"


class Experiment4DTests(unittest.TestCase):
    def setUp(self):
        self.circuit = gr.parse_circuit(INPUT)
        self.adjacency, _, _, self.root, _, _ = e4c.graph_analysis(self.circuit)
        self.affinities = e4d.affinity_values(self.adjacency, self.root)

    def test_direct_shared_root_exclusion_and_symmetry(self):
        for first, related in self.affinities.items():
            for second, value in related.items():
                direct = int(second in self.adjacency[first])
                shared = len((set(self.adjacency[first]) - {self.root}) &
                             (set(self.adjacency[second]) - {self.root}))
                self.assertEqual(direct + shared, value)
                self.assertEqual(value, self.affinities[second][first])
        a, b = "TPB", "TPO"
        with_root = len(set(self.adjacency[a]) & set(self.adjacency[b]))
        without_root = len((set(self.adjacency[a]) - {self.root}) &
                           (set(self.adjacency[b]) - {self.root}))
        self.assertEqual(with_root - 1, without_root)

    def test_weighted_barycenters_are_synchronous_and_deterministic(self):
        current, radius = e4d.provisional_positions(self.circuit, self.root)
        low, high = radius - gr.DEGREE_RING_SPACING / 2, radius + gr.DEGREE_RING_SPACING / 2
        first = e4d.affinity_proposals(current, self.affinities, self.root, low, high)
        second = e4d.affinity_proposals(dict(reversed(list(current.items()))),
                                        self.affinities, self.root, low, high)
        self.assertEqual(first, second)
        ref = sorted(self.affinities)[0]
        weights = self.affinities[ref]
        total = sum(weights.values())
        desired = (sum(w * current[n][0] for n, w in weights.items()) / total,
                   sum(w * current[n][1] for n, w in weights.items()) / total)
        self.assertEqual(e4d.clamp_radial(desired, current[ref], low, high), first[ref])

    def test_collision_resolution_and_completed_initialization(self):
        first = e4d.affinity_initial_positions(self.circuit)
        second = e4d.affinity_initial_positions(self.circuit)
        self.assertEqual(first, second)
        self.assertEqual(0, gr.evaluate(self.circuit, first[0])["component_overlaps"])

    def test_release_removes_constraints(self):
        positions, _, _ = e4d.affinity_initial_positions(self.circuit)
        released = e4d.release_initialization(positions)
        self.assertEqual(positions, released)
        self.assertIsNot(positions, released)
        released[self.root] = (1.0, 2.0)
        self.assertEqual((1.0, 2.0), released[self.root])

    def test_short_repeated_run_determinism(self):
        positions, _, _ = e4d.affinity_initial_positions(self.circuit)
        first = gr.anneal(self.circuit, e4d.release_initialization(positions), iterations=30)
        second = gr.anneal(self.circuit, e4d.release_initialization(positions), iterations=30)
        self.assertEqual(first, second)


if __name__ == "__main__":
    unittest.main()
