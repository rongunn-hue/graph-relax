import hashlib,tempfile,unittest
from pathlib import Path
import graph_first_experiment_metric_protected as ex
class Tests(unittest.TestCase):
 def go(self,d):
  b=Path('output/graph_first');return ex.run(b/'iamp_graph.json',b/'iamp_nearness_rotation_direct_report.json',Path(d))
 def test_protections_and_hard_validity(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.go(d);p=r['protection_report'];self.assertTrue(p['all_satisfied_at_final']);self.assertTrue(all(x['final_d']>=x['D_lock']-1e-7 for x in p['creation_order']));self.assertEqual((r['final_validation']['vertices'],r['final_validation']['edges'],r['final_validation']['proper_crossings'],r['final_validation']['coincidences'],r['final_validation']['vertex_on_unrelated_edge']),(34,44,0,0,0));self.assertTrue(r['final_validation']['endpoint_pairs_identical'])
 def test_determinism(self):
  with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
   x=self.go(a);y=self.go(b);self.assertEqual(x,y)
   for n in ['iamp_metric_protected_start.svg','iamp_metric_protected_polygon.svg','iamp_metric_protected_nodes.svg','iamp_metric_protected_final.svg','iamp_metric_protected_report.json']:self.assertEqual(hashlib.sha256((Path(a)/n).read_bytes()).digest(),hashlib.sha256((Path(b)/n).read_bytes()).digest())
if __name__=='__main__':unittest.main()
