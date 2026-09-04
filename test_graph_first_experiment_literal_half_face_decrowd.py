import hashlib,tempfile,unittest
from pathlib import Path
import graph_first_experiment_literal_half_face_decrowd as experiment

class LiteralHalfTests(unittest.TestCase):
 def run1(self,d):
  b=Path("output/graph_first");return experiment.run(b/"iamp_graph.json",b/"iamp_corridor_separation_report.json",Path(d))
 def test_literal_half_and_hard_only_reduction(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.run1(d);moves=[x for x in r["trace"] if x.get("status") in ("MOVED_UP","MOVED_DOWN")]
   for x in moves:
    self.assertEqual(x["target_formula"],"TARGET = full_gap / 2");self.assertFalse(x["target_uses_readability_clearance"]);self.assertAlmostEqual(x["half_gap"],x["full_gap"]/2,12)
    if x["exact_half_hard_valid"]:self.assertAlmostEqual(x["actual_move"],x["half_gap"],9)
    else:self.assertLess(x["actual_move"],x["half_gap"]);self.assertTrue(x["exact_half_hard_failures"])
 def test_alternation_c1_c2_and_valid_stages(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.run1(d);a=r["literal_half_assertions"];self.assertEqual(a["alternating_execution"],["TOP","BOTTOM","TOP","BOTTOM"]);self.assertTrue(a["second_top_contains_C1"]);self.assertTrue(a["second_bottom_contains_C2"])
   for s in r["stages"]:self.assertEqual(s["metrics"]["proper_crossings"],0);self.assertEqual(s["metrics"]["coincident_vertices"],0);self.assertEqual(s["metrics"]["vertex_on_unrelated_edge_events"],0);self.assertEqual(len(s["coordinates"]),34)
 def test_deterministic_outputs(self):
  with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
   self.run1(a);self.run1(b)
   for n in ("iamp_literal_half_face_decrowd.svg","iamp_literal_half_face_decrowd_report.json"):self.assertEqual(hashlib.sha256((Path(a)/n).read_bytes()).digest(),hashlib.sha256((Path(b)/n).read_bytes()).digest())
if __name__=="__main__":unittest.main()
