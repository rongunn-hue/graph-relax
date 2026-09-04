import hashlib,tempfile,unittest
from pathlib import Path
import graph_first_experiment_voronoi_centroid_move as experiment

class CentroidMoveTests(unittest.TestCase):
 def run1(self,d):
  b=Path("output/graph_first"); return experiment.run(b/"iamp_graph.json",b/"iamp_global_node_density_report.json",Path(d))
 def test_one_exact_attempt_and_final_invariants(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.run1(d); self.assertEqual(r["attempt_count"],34); self.assertEqual(r["accepted_moves"]+r["rejected_moves"],34)
   self.assertTrue(r["frozen_initial_partition"]); self.assertEqual({x["complete_edge_regeneration_count"] for x in r["attempts"]},{44})
   self.assertEqual((r["vertex_count"],r["electrical_incidence_count"]),(34,44)); self.assertEqual(r["final_validation"]["proper_crossings"],0)
   self.assertAlmostEqual(r["final_partition_area_error"],0,7)
 def test_determinism(self):
  with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
   ra=self.run1(a); rb=self.run1(b); self.assertEqual(ra["attempts"],rb["attempts"]); self.assertEqual(ra["final_coordinates"],rb["final_coordinates"])
   for n in ("iamp_voronoi_centroid_move.svg","iamp_voronoi_centroid_after.svg","iamp_voronoi_centroid_move_report.json"):
    self.assertEqual(hashlib.sha256((Path(a)/n).read_bytes()).digest(),hashlib.sha256((Path(b)/n).read_bytes()).digest())
if __name__=="__main__": unittest.main()
