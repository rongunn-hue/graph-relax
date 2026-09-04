import hashlib
import tempfile
import unittest
from pathlib import Path

import graph_first_experiment_corridor_half_available as experiment


class HalfAvailableCorridorTests(unittest.TestCase):
    def run_once(self,directory):
        b=Path('output/graph_first');return experiment.run(b/'iamp_graph.json',b/'iamp_corridor_separation_report.json',Path(directory))

    def test_finite_side_a_and_unbounded_side_b_stop_without_motion(self):
        with tempfile.TemporaryDirectory() as d:
            r=self.run_once(d)
            self.assertAlmostEqual(r['side_a_available_space']['D_A'],307.98510481975853,10)
            self.assertAlmostEqual(r['side_a_available_space']['delta_A'],r['side_a_available_space']['D_A']/2,12)
            self.assertIsNone(r['side_b_available_space']['D_B']);self.assertIsNone(r['side_b_available_space']['delta_B'])
            self.assertFalse(r['operation_applied']);self.assertTrue(r['coordinates_unchanged']);self.assertTrue(r['graph_was_not_modified'])

    def test_baseline_and_connectivity_preserved(self):
        with tempfile.TemporaryDirectory() as d:
            r=self.run_once(d);self.assertEqual((r['vertex_count'],r['electrical_incidence_count']),(34,44));self.assertEqual(r['baseline_metrics']['proper_crossings'],0);self.assertEqual(r['baseline_metrics']['clearance_violations'],6);self.assertTrue(r['edge_reconnection']['valid'])

    def test_deterministic_outputs(self):
        with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
            self.run_once(a);self.run_once(b)
            for name in ('iamp_corridor_half_available.svg','iamp_corridor_half_available_report.json'):
                self.assertEqual(hashlib.sha256((Path(a)/name).read_bytes()).digest(),hashlib.sha256((Path(b)/name).read_bytes()).digest())


if __name__=='__main__':unittest.main()
