import hashlib,tempfile,unittest
from pathlib import Path
import graph_first_experiment_voronoi_threshold80_monotonic as experiment
class Tests(unittest.TestCase):
 def run1(self,d):
  b=Path("output/graph_first");return experiment.run(b/"iamp_graph.json",b/"iamp_voronoi_centroid_move_report.json",b/"iamp_voronoi_threshold80_report.json",b/"iamp_voronoi_threshold80_deficiency_report.json",Path(d))
 def test_monotonic_rule(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.run1(d);self.assertEqual(len(r["attempts"]),r["active_count"]);self.assertTrue(all(a["graph_valid"] and a["total_deficiency_decreased"] and a["all_currently_deficient_non_worsening"] for a in r["attempts"] if a["decision"]=="ACCEPTED"));self.assertEqual({a["complete_edge_regeneration_count"] for a in r["attempts"]},{44});self.assertLess(r["final_total_deficiency"],r["starting_total_deficiency"]);self.assertFalse(any(r["final_validation"][k] for k in ("proper_crossings","coincident_vertices","vertex_on_unrelated_edge_events")))
 def test_determinism(self):
  with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
   x=self.run1(a);y=self.run1(b);self.assertEqual(x["attempts"],y["attempts"]);self.assertEqual(x["final_coordinates"],y["final_coordinates"])
   for n in ("iamp_voronoi_threshold80_monotonic.svg","iamp_voronoi_threshold80_monotonic_report.json"):self.assertEqual(hashlib.sha256((Path(a)/n).read_bytes()).digest(),hashlib.sha256((Path(b)/n).read_bytes()).digest())
if __name__=="__main__":unittest.main()
