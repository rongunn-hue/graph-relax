import hashlib,tempfile,unittest
from pathlib import Path
import graph_first_experiment_three_case_decrowd as ex
class Tests(unittest.TestCase):
 def go(self,d):
  b=Path('output/graph_first');return ex.run(b/'iamp_graph.json',b/'iamp_nearness_rotation_direct_report.json',Path(d))
 def test_three_stages_and_validity(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.go(d);self.assertIn('component:TPB',[x['vertex'] for x in r['initial_crowding_map']['non_polygon_nodes']]);self.assertEqual((r['final_validation']['vertices'],r['final_validation']['edges'],r['final_validation']['proper_crossings']),(34,44,0));self.assertEqual(r['final_validation']['coincidences'],0);self.assertEqual(r['final_validation']['vertex_on_unrelated_edge'],0)
 def test_determinism(self):
  with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
   x=self.go(a);y=self.go(b);self.assertEqual(x,y)
   for n in ('iamp_three_case_decrowd_start.svg','iamp_three_case_decrowd_polygon.svg','iamp_three_case_decrowd_nodes.svg','iamp_three_case_decrowd_final.svg','iamp_three_case_decrowd_report.json'):self.assertEqual(hashlib.sha256((Path(a)/n).read_bytes()).digest(),hashlib.sha256((Path(b)/n).read_bytes()).digest())
if __name__=='__main__':unittest.main()
