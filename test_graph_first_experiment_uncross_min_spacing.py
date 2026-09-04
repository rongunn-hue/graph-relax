import hashlib,tempfile,unittest
from pathlib import Path
import graph_first_experiment_uncross_min_spacing as experiment
class Tests(unittest.TestCase):
 def run1(self,d):
  b=Path("output/graph_first");return experiment.run(b/"iamp_graph.json",b/"iamp_nearness_mds.json",b/"iamp_nearness_mds_report.json",b/"iamp_nearness_rotation_direct_report.json",Path(d))
 def test_spacing_and_direct_cardinality(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.run1(d);self.assertGreaterEqual(r["sanitation"]["post_sanitation_minimum_node_distance"]+1e-9,r["S_MIN"]);self.assertTrue(r["uncrossing"]["S_MIN_preserved_after_every_accepted_move"]);self.assertLessEqual(r["uncrossing"]["maximum_positions_per_orientation"],1);self.assertLessEqual(r["uncrossing"]["maximum_orientations_per_crossing"],4);self.assertEqual((r["vertex_count"],r["electrical_edge_count"]),(34,44))
   self.assertEqual(r["sanitation"]["all_distinct_pairs_audited"],561);self.assertEqual(r["final"]["all_distinct_pair_count_audited"],561);self.assertTrue(r["uncrossing"]["moving_pivot_pairs_globally_valid_after_every_move"])
 def test_determinism(self):
  with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
   x=self.run1(a);y=self.run1(b);self.assertEqual(x["sanitation"],y["sanitation"]);self.assertEqual(x["uncrossing"]["moves"],y["uncrossing"]["moves"]);self.assertEqual(x["final"]["coordinates"],y["final"]["coordinates"])
   for n in ("iamp_uncross_with_min_spacing.svg","iamp_uncross_with_min_spacing_report.json"):self.assertEqual(hashlib.sha256((Path(a)/n).read_bytes()).digest(),hashlib.sha256((Path(b)/n).read_bytes()).digest())
if __name__=="__main__":unittest.main()
