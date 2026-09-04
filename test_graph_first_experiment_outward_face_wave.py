import hashlib,tempfile,unittest
from pathlib import Path
import graph_first_experiment_outward_face_wave as experiment

class FaceWaveTests(unittest.TestCase):
 def run1(self,d):
  b=Path("output/graph_first");return experiment.run(b/"iamp_graph.json",b/"iamp_corridor_separation_report.json",Path(d))
 def test_alternates_and_reaches_c_vertices(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.run1(d);moves=[x for x in r["trace"] if x["status"] in ("MOVED_UP","MOVED_DOWN")]
   self.assertEqual([x["side"] for x in moves[:4]],["TOP","BOTTOM","TOP","BOTTOM"]);self.assertTrue(r["second_top_reaches_c1_face"]);self.assertTrue(r["second_bottom_reaches_c2_face"])
   for x in moves:self.assertIn(x["next_face_id"],x["all_adjacent_faces"])
 def test_geometry_and_movement_invariants(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.run1(d);seen=set()
   for i,s in enumerate(r["stages"][1:],1):
    op=next(x for x in r["trace"] if x.get("side")==s["side"] and x.get("face_id")==s["face_id"] and x.get("status") in ("MOVED_UP","MOVED_DOWN"));prior=r["stages"][i-1]["coordinates"]
    self.assertLessEqual(op["actual_move"],op["half_gap"]+1e-9);self.assertAlmostEqual(op["actual_move"],min(op["max_valid_move"],op["half_gap"]),9);self.assertTrue(s["edge_reconnection"]["valid"]);self.assertEqual(s["metrics"]["proper_crossings"],0);self.assertEqual(s["metrics"]["coincident_vertices"],0);self.assertEqual(s["metrics"]["vertex_on_unrelated_edge_events"],0)
    for v in op["moved_vertices"]:
     self.assertNotIn(v,seen);seen.add(v);self.assertEqual(prior[v][0],s["coordinates"][v][0]);self.assertLess(s["coordinates"][v][1],prior[v][1]) if s["side"]=="TOP" else self.assertGreater(s["coordinates"][v][1],prior[v][1])
 def test_deterministic(self):
  with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
   r=self.run1(a);self.run1(b);self.assertEqual((r["vertex_count"],r["electrical_incidence_count"]),(34,44))
   for n in ("iamp_outward_face_wave.svg","iamp_outward_face_wave_report.json"):self.assertEqual(hashlib.sha256((Path(a)/n).read_bytes()).digest(),hashlib.sha256((Path(b)/n).read_bytes()).digest())
if __name__=="__main__":unittest.main()
