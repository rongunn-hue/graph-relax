import hashlib,tempfile,unittest
from pathlib import Path
import graph_first_experiment_common_metric_guard as ex
class Tests(unittest.TestCase):
 def go(self,d):
  b=Path('output/graph_first');return ex.run(b/'iamp_graph.json',b/'iamp_nearness_rotation_direct_report.json',Path(d))
 def test_guard_and_validity(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.go(d);self.assertEqual((r['final_validation']['vertices'],r['final_validation']['edges'],r['final_validation']['proper_crossings'],r['final_validation']['coincidences'],r['final_validation']['vertex_on_unrelated_edge']),(34,44,0,0,0));self.assertTrue(r['final_validation']['endpoint_pairs_identical'])
   for stage in ('node_stage','edge_stage'):
    for o in r[stage]['operations']:
     for c in o.get('candidate_results',[o]):
      if c.get('decision')=='ACCEPTED':self.assertFalse(c['gate']['new_crowding']);self.assertFalse(c['gate']['existing_crowding_worsened'])
 def test_determinism(self):
  with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
   x=self.go(a);y=self.go(b);self.assertEqual(x,y)
   for n in ['iamp_common_metric_guard_start.svg','iamp_common_metric_guard_polygon.svg','iamp_common_metric_guard_nodes.svg','iamp_common_metric_guard_final.svg','iamp_common_metric_guard_report.json']:self.assertEqual(hashlib.sha256((Path(a)/n).read_bytes()).digest(),hashlib.sha256((Path(b)/n).read_bytes()).digest())
if __name__=='__main__':unittest.main()
