import hashlib,tempfile,unittest
from pathlib import Path
import graph_first_experiment_corrected_crowding as ex
class Tests(unittest.TestCase):
 def go(self,d):
  b=Path('output/graph_first');return ex.run(b/'iamp_graph.json',b/'iamp_nearness_rotation_direct_report.json',Path(d))
 def test_detector_and_validity(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.go(d);self.assertFalse(r['explicit_audits']['TPO']['flagged']);self.assertIn('SHARED_ENDPOINT_EDGE_OVERLAY',r['explicit_audits']['TPB']['categories']);self.assertTrue(all(x['vertex'] not in x['edge'] for x in r['initial_maps']['node_edge']));self.assertEqual((r['final_validation']['vertices'],r['final_validation']['edges'],r['final_validation']['crossings']),(34,44,0));self.assertEqual(r['final_validation']['coincidences'],0);self.assertEqual(r['final_validation']['vertex_on_unrelated_edge'],0)
 def test_determinism(self):
  with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
   x=self.go(a);y=self.go(b);self.assertEqual(x,y)
   for n in ('iamp_corrected_crowding_start.svg','iamp_corrected_crowding_polygon.svg','iamp_corrected_crowding_nodes.svg','iamp_corrected_crowding_final.svg','iamp_corrected_crowding_report.json'):self.assertEqual(hashlib.sha256((Path(a)/n).read_bytes()).digest(),hashlib.sha256((Path(b)/n).read_bytes()).digest())
if __name__=='__main__':unittest.main()
