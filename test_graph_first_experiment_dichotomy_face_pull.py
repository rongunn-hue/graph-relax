import hashlib,tempfile,unittest
from pathlib import Path
import graph_first_experiment_dichotomy_face_pull as experiment
class DichotomyTests(unittest.TestCase):
 def run1(self,d):
  b=Path("output/graph_first");return experiment.run(b/"iamp_graph.json",b/"iamp_corridor_separation_report.json",Path(d))
 def test_existing_poles_safe_half_and_validity(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.run1(d);self.assertEqual((r["vertex_count"],r["electrical_incidence_count"]),(34,44))
   for x in r["trace"]:
    if x["status"]=="ACCEPTED":
     self.assertIn(x["pull_pole_vertex_id"],x["boundary_vertex_ids"]);self.assertTrue(x["edge_reconnection"]["valid"]);self.assertEqual(x["metrics_after"]["proper_crossings"],0)
     for m in x["movements"]:self.assertNotEqual(m["original_coordinates"],m["new_coordinates"])
 def test_no_duplicate_movement_and_determinism(self):
  with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
   r=self.run1(a);self.run1(b);vs=[m["vertex"] for x in r["trace"] if x["status"]=="ACCEPTED" for m in x["movements"]];self.assertEqual(len(vs),len(set(vs)))
   for n in ("iamp_dichotomy_face_pull.svg","iamp_dichotomy_face_pull_clean.svg","iamp_dichotomy_face_pull_report.json"):self.assertEqual(hashlib.sha256((Path(a)/n).read_bytes()).digest(),hashlib.sha256((Path(b)/n).read_bytes()).digest())
if __name__=="__main__":unittest.main()
