import hashlib,tempfile,unittest
from pathlib import Path
import graph_first_experiment_node_separation_150 as experiment
class Tests(unittest.TestCase):
 def run1(self,d):
  b=Path("output/graph_first");return experiment.run(b/"iamp_graph.json",b/"iamp_nearness_mds.json",b/"iamp_nearness_mds_report.json",Path(d))
 def test_all_pairs_and_connectivity(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.run1(d);self.assertTrue(r["success"]);self.assertEqual(r["distinct_pairs_audited"],561);self.assertEqual(r["final_pairs_below_150"],0);self.assertGreaterEqual(r["final_minimum_pairwise_distance"],150);self.assertEqual((r["vertex_count"],r["electrical_edge_count"],r["rebuilt_edge_count"]),(34,44,44))
 def test_determinism(self):
  with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
   x=self.run1(a);y=self.run1(b);self.assertEqual(x["corrections"],y["corrections"]);self.assertEqual(x["final_coordinates"],y["final_coordinates"])
   for n in ("iamp_node_separation_150.svg","iamp_node_separation_150_report.json"):self.assertEqual(hashlib.sha256((Path(a)/n).read_bytes()).digest(),hashlib.sha256((Path(b)/n).read_bytes()).digest())
if __name__=="__main__":unittest.main()
