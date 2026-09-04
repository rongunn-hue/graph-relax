import hashlib,tempfile,unittest
from pathlib import Path
import graph_first_experiment_voronoi_centroid_pass2 as experiment

class Pass2Tests(unittest.TestCase):
 def run1(self,d):
  b=Path("output/graph_first"); return experiment.run(b/"iamp_graph.json",b/"iamp_voronoi_centroid_move_report.json",Path(d))
 def test_exact_once_and_invariants(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.run1(d); self.assertEqual(len(r["pass2_attempts"]),34); self.assertEqual(r["pass2_accepted_moves"]+r["pass2_rejected_moves"],34)
   self.assertTrue(r["frozen_pass2_partition"]); self.assertEqual({x["complete_edge_regeneration_count"] for x in r["pass2_attempts"]},{44})
   self.assertEqual((r["vertex_count"],r["electrical_incidence_count"]),(34,44)); self.assertFalse(any(r["pass2_final_validation"][k] for k in ("proper_crossings","coincident_vertices","vertex_on_unrelated_edge_events")))
   self.assertAlmostEqual(r["pass2_final_partition_area_error"],0,7)
 def test_determinism(self):
  with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
   ra=self.run1(a); rb=self.run1(b); self.assertEqual(ra["pass2_attempts"],rb["pass2_attempts"]); self.assertEqual(ra["pass2_final_coordinates"],rb["pass2_final_coordinates"])
   for n in ("iamp_voronoi_centroid_pass2.svg","iamp_voronoi_centroid_pass2_report.json"):
    self.assertEqual(hashlib.sha256((Path(a)/n).read_bytes()).digest(),hashlib.sha256((Path(b)/n).read_bytes()).digest())
if __name__=="__main__": unittest.main()
