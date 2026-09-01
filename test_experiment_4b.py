import unittest
from pathlib import Path

import experiment_4b as e4b
import graph_relax as gr


INPUT = Path(__file__).parent / "test_data" / "iamp.circuit"


class Experiment4BTests(unittest.TestCase):
    def setUp(self):
        self.circuit = gr.parse_circuit(INPUT)
        self.positions, self.zmax = e4b.degree_seeded_positions_3d(self.circuit)

    def test_3d_distances_and_segment_separation(self):
        self.assertEqual(5.0, e4b.norm((3.0, 4.0, 0.0)))
        self.assertAlmostEqual(10.0, e4b.segment_distance_3d(
            (0, 0, 0), (10, 0, 0), (0, 0, 10), (10, 0, 10)))

    def test_z_bound_and_deterministic_movement(self):
        first = e4b.anneal_3d(self.circuit, self.positions, self.zmax, iterations=30)
        second = e4b.anneal_3d(self.circuit, self.positions, self.zmax, iterations=30)
        self.assertEqual(first, second)
        self.assertTrue(all(-self.zmax <= p[2] <= self.zmax for p in first[0].values()))

    def test_orientation_generation_is_deterministic(self):
        first, second = e4b.orientations(8, 3), e4b.orientations(8, 3)
        self.assertEqual(first, second)
        self.assertEqual(24, len(first))

    def test_rigid_projection_and_2d_reconstruction(self):
        orientation = e4b.orientations(2, 1)[0]
        projected = e4b.project_centers(self.positions, orientation)
        a, b = self.circuit.refs[:2]
        original_delta = e4b.subtract(self.positions[a], self.positions[b])
        projected_distance = gr.math.dist(projected[a], projected[b])
        self.assertLessEqual(projected_distance, e4b.norm(original_delta) + 1e-9)
        terminals = gr.terminal_positions(self.circuit, projected)
        self.assertEqual(sum(len(pins) for pins in self.circuit.pins.values()), len(terminals))
        self.assertIn("objective", gr.evaluate(self.circuit, projected))

    def test_projection_selection_and_short_repeatability(self):
        subset = e4b.orientations(4, 2)
        first = e4b.projection_search(self.circuit, self.positions, subset)
        second = e4b.projection_search(self.circuit, self.positions, subset)
        self.assertEqual(first, second)
        objectives = [candidate["objective"] for candidate in first[4]]
        self.assertEqual(min(objectives), gr.rounded(first[1]["objective"]))


if __name__ == "__main__":
    unittest.main()
