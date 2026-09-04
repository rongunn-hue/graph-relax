import hashlib,tempfile,unittest
from pathlib import Path
import graph_first_experiment_voronoi_threshold80_deficiency as experiment
class Tests(unittest.TestCase):
 def run1(self,d):
  b=Path("output/graph_first");return experiment.run(b/"iamp_graph.json",b/"iamp_voronoi_centroid_move_report.json",b/"iamp_voronoi_threshold80_report.json",Path(d))
 def test_acceptance_and_validity(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.run1(d);self.assertEqual(len(r["attempts"]),r["active_count"]);self.assertTrue(all(x["graph_valid"] and x["deficiency_improving"] for x in r["attempts"] if x["decision"]=="ACCEPTED"));self.assertTrue(all(x["complete_edge_regeneration_count"]==44 for x in r["attempts"]));self.assertLess(r["final_total_deficiency"],r["starting_total_deficiency"]);self.assertFalse(any(r["final_validation"][k] for k in ("proper_crossings","coincident_vertices","vertex_on_unrelated_edge_events")))
 def test_determinism(self):
  with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
   ra=self.run1(a);rb=self.run1(b);self.assertEqual(ra["active_vertices"],rb["active_vertices"]);self.assertEqual(ra["attempts"],rb["attempts"]);self.assertEqual(ra["final_coordinates"],rb["final_coordinates"])
   for n in ("iamp_voronoi_threshold80_deficiency.svg","iamp_voronoi_threshold80_deficiency_report.json"):self.assertEqual(hashlib.sha256((Path(a)/n).read_bytes()).digest(),hashlib.sha256((Path(b)/n).read_bytes()).digest())
if __name__=="__main__":unittest.main()
