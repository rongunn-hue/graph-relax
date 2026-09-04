import hashlib,tempfile,unittest
from pathlib import Path
import graph_first_experiment_voronoi_threshold80 as experiment
class ThresholdTests(unittest.TestCase):
 def run1(self,d):
  b=Path("output/graph_first");return experiment.run(b/"iamp_graph.json",b/"iamp_voronoi_centroid_move_report.json",Path(d))
 def test_threshold_protocol(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.run1(d);self.assertEqual(r["active_count"]+r["frozen_count"],34);self.assertEqual(len(r["attempts"]),r["active_count"]);self.assertEqual({x["complete_edge_regeneration_count"] for x in r["attempts"]},{44});self.assertTrue(all(x["pass1_starting_area"]>=r["minimum_acceptable_area"] and not x["centroid_attempted"] for x in r["per_vertex"] if x["classification"]=="FROZEN"));self.assertFalse(any(r["final_validation"][k] for k in ("proper_crossings","coincident_vertices","vertex_on_unrelated_edge_events")))
 def test_determinism(self):
  with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
   ra=self.run1(a);rb=self.run1(b);self.assertEqual(ra["active_vertices"],rb["active_vertices"]);self.assertEqual(ra["attempts"],rb["attempts"]);self.assertEqual(ra["final_coordinates"],rb["final_coordinates"])
   for n in ("iamp_voronoi_threshold80.svg","iamp_voronoi_threshold80_report.json"):self.assertEqual(hashlib.sha256((Path(a)/n).read_bytes()).digest(),hashlib.sha256((Path(b)/n).read_bytes()).digest())
if __name__=="__main__":unittest.main()
