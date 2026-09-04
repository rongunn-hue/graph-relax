import hashlib,math,tempfile,unittest
from pathlib import Path
import graph_first_experiment_vertical_halfspace_opening as experiment

class VerticalOpeningTests(unittest.TestCase):
 def run_once(self,d):
  b=Path('output/graph_first');return experiment.run(b/'iamp_graph.json',b/'iamp_corridor_separation_report.json',Path(d))
 def test_baseline_vertical_only_and_reconnection(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.run_once(d);self.assertEqual((r['vertex_count'],r['electrical_incidence_count']),(34,44));self.assertEqual(r['before_metrics']['proper_crossings'],0);self.assertEqual(r['before_metrics']['clearance_violations'],6);self.assertTrue(r['edge_reconnection']['valid'])
   for v,p in r['coordinates_before'].items():self.assertEqual(p[0],r['coordinates_after'][v][0])
   for v in r['upper_movable_vertices']:self.assertLess(r['coordinates_after'][v][1],r['coordinates_before'][v][1])
   for v in r['lower_movable_vertices']:self.assertGreater(r['coordinates_after'][v][1],r['coordinates_before'][v][1])
 def test_face_depth_and_half_distance(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.run_once(d);f=r['face_region_depth'];self.assertAlmostEqual(f['MOVE_UP'],f['H_UPPER']/2,12);self.assertAlmostEqual(f['MOVE_DOWN'],f['H_LOWER']/2,12);self.assertGreater(f['MOVE_UP'],0);self.assertGreater(f['MOVE_DOWN'],0)
 def test_deterministic(self):
  with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
   self.run_once(a);self.run_once(b)
   for n in ('iamp_vertical_halfspace_opening.svg','iamp_vertical_halfspace_opening_report.json'):self.assertEqual(hashlib.sha256((Path(a)/n).read_bytes()).digest(),hashlib.sha256((Path(b)/n).read_bytes()).digest())
if __name__=='__main__':unittest.main()
