import hashlib,tempfile,unittest
from pathlib import Path
import graph_first_experiment_uncross_from_sep150 as experiment
class Tests(unittest.TestCase):
 def run1(self,d):
  b=Path("output/graph_first");return experiment.run(b/"iamp_graph.json",b/"iamp_node_separation_150_report.json",Path(d))
 def test_literal_spacing_and_directness(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.run1(d);self.assertEqual(r["starting_audit"]["pairs_audited"],561);self.assertGreaterEqual(r["starting_audit"]["minimum_pairwise_distance"],150);self.assertGreaterEqual(r["final_audit"]["minimum_pairwise_distance"],150);self.assertTrue(r["full_150_spacing_preserved"]);self.assertLessEqual(r["maximum_positions_per_orientation"],1);self.assertLessEqual(r["maximum_orientations_per_crossing"],4);self.assertEqual((r["final_audit"]["vertex_count"],r["final_audit"]["edge_count"]),(34,44))
   self.assertEqual(r["relaxed_reducing_candidate_count"],0)
   self.assertEqual(len(r["relaxed_spacing_fallback_audits"]),3)
   self.assertEqual(sum(len(x["solutions"]) for x in r["relaxed_spacing_fallback_audits"]),12)
   self.assertFalse(r["spacing_policy"]["fallback"]["parameter_sweep"])
 def test_determinism(self):
  with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
   x=self.run1(a);y=self.run1(b);self.assertEqual(x["crossing_history"],y["crossing_history"]);self.assertEqual(x["orientation_audits"],y["orientation_audits"]);self.assertEqual(x["final_coordinates"],y["final_coordinates"])
   for n in ("iamp_uncross_from_sep150.svg","iamp_uncross_from_sep150_report.json"):self.assertEqual(hashlib.sha256((Path(a)/n).read_bytes()).digest(),hashlib.sha256((Path(b)/n).read_bytes()).digest())
if __name__=="__main__":unittest.main()
