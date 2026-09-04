import hashlib,tempfile,unittest
from pathlib import Path
import graph_first_experiment_polygon_area_expansion as ex
class Tests(unittest.TestCase):
 def go(self,d):
  b=Path('output/graph_first');return ex.run(b/'iamp_graph.json',b/'iamp_nearness_rotation_direct_report.json',Path(d))
 def test_expansion_and_validity(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.go(d);self.assertEqual(r['final_validation']['proper_crossings'],0);self.assertEqual((r['final_validation']['vertices'],r['final_validation']['edges']),(34,44));self.assertTrue(r['final_validation']['connectivity_unchanged'])
   for x in r['operations']:
    self.assertTrue(x['fixed_coordinates_unchanged'] if 'fixed_coordinates_unchanged' in x else True)
    if x['status']=='EXPANDED':self.assertGreater(x['area_after'],x['area_before']);self.assertGreater(x['metrics_after']['minimum_nearest_boundary_node_distance'],x['metrics_before']['minimum_nearest_boundary_node_distance']);self.assertGreater(x['metrics_after']['average_nearest_boundary_node_distance'],x['metrics_before']['average_nearest_boundary_node_distance'])
 def test_determinism(self):
  with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
   x=self.go(a);y=self.go(b);self.assertEqual(x,y)
   for n in ('iamp_polygon_area_expansion_start.svg','iamp_polygon_area_expansion_final.svg','iamp_polygon_area_expansion_report.json'):self.assertEqual(hashlib.sha256((Path(a)/n).read_bytes()).digest(),hashlib.sha256((Path(b)/n).read_bytes()).digest())
if __name__=='__main__':unittest.main()
