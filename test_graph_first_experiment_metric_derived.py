import hashlib,tempfile,unittest
from pathlib import Path
import graph_first_experiment_metric_derived as ex
class Tests(unittest.TestCase):
 def go(self,d):
  b=Path('output/graph_first');return ex.run(b/'iamp_graph.json',b/'iamp_nearness_rotation_direct_report.json',Path(d))
 def test_analytic_scale_aggregate_and_validity(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.go(d);self.assertTrue(all(o['candidate_count']<=1 for o in r['polygon_stage']['operations']));self.assertTrue(all(x['s_pair'] is None or x['s_pair']>=1 for o in r['polygon_stage']['operations'] for x in o['deficient_NN_relationships']));self.assertEqual((r['final_validation']['vertices'],r['final_validation']['edges'],r['final_validation']['proper_crossings'],r['final_validation']['coincidences'],r['final_validation']['vertex_on_unrelated_edge']),(34,44,0,0,0));self.assertTrue(r['final_validation']['endpoint_pairs_identical'])
   for stage in ('node_stage','edge_stage'):
    for o in r[stage]['operations']:
     for c in o.get('candidate_results',[o]):
      if c.get('decision')=='ACCEPTED':self.assertLess(c['gate']['P_after'],c['gate']['P_before'])
 def test_determinism(self):
  with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
   x=self.go(a);y=self.go(b);self.assertEqual(x,y)
   for n in ['iamp_metric_derived_start.svg','iamp_metric_derived_polygon.svg','iamp_metric_derived_nodes.svg','iamp_metric_derived_final.svg','iamp_metric_derived_report.json']:self.assertEqual(hashlib.sha256((Path(a)/n).read_bytes()).digest(),hashlib.sha256((Path(b)/n).read_bytes()).digest())
if __name__=='__main__':unittest.main()
