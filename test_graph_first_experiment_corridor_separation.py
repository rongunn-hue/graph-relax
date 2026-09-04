import hashlib
import json
import math
import tempfile
import unittest
from pathlib import Path

import graph_first_experiment_corridor_separation as experiment
import placement_geometry as geometry


class CorridorSeparationTests(unittest.TestCase):
    def run_once(self, directory):
        base=Path('output/graph_first')
        return experiment.run(base/'iamp_graph.json',base/'iamp_sequential_local_face_pulls_report.json',Path(directory))

    def test_topological_boundaries_outward_motion_and_reconnection(self):
        with tempfile.TemporaryDirectory() as d:
            r=self.run_once(d)
            self.assertEqual((r['vertex_count'],r['electrical_incidence_count']),(34,44))
            self.assertEqual(r['before_metrics']['proper_crossings'],0)
            self.assertEqual(r['before_metrics']['clearance_violations'],6)
            self.assertEqual(r['corridor']['rmin1_side_boundary']['ordered_vertices'],['component:U1','net:EA','component:J_EA','net:GND'])
            self.assertEqual(r['corridor']['j_out_side_boundary']['ordered_vertices'],['component:U1','net:EOUT','component:J_OUT','net:GND'])
            self.assertFalse(set(r['rmin1_side_movable_vertices']) & set(r['j_out_side_movable_vertices']))
            self.assertTrue(all(x['absolute_signed_distance_increased'] for x in r['movement']))
            self.assertTrue(r['edge_reconnection']['valid'])
            self.assertEqual(r['edge_reconnection']['graph_edge_count'],44)
            self.assertEqual(r['after_metrics']['proper_crossings'],0)

    def test_fixed_vertices_and_uniform_displacement(self):
        with tempfile.TemporaryDirectory() as d:
            r=self.run_once(d);moved={x['vertex'] for x in r['movement']}
            for v,p in r['coordinates_before'].items():
                if v not in moved:self.assertEqual(p,r['coordinates_after'][v])
            for x in r['movement']:
                self.assertAlmostEqual(math.hypot(*x['displacement']),experiment.DELTA,12)

    def test_byte_deterministic(self):
        with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
            self.run_once(a);self.run_once(b)
            for name in ('iamp_corridor_separation.json','iamp_corridor_separation.svg'):
                # Report has the required longer filename.
                actual='iamp_corridor_separation_report.json' if name.endswith('.json') else name
                self.assertEqual(hashlib.sha256((Path(a)/actual).read_bytes()).digest(),hashlib.sha256((Path(b)/actual).read_bytes()).digest())


if __name__=='__main__':unittest.main()
