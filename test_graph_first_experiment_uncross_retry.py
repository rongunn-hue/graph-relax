import hashlib,tempfile,unittest
from pathlib import Path
import graph_first_experiment_uncross_retry as experiment

class RetryTests(unittest.TestCase):
 def run1(self,d):
  b=Path("output/graph_first")
  return experiment.run(b/"iamp_graph.json",b/"iamp_uncross_from_sep150_report.json",Path(d))
 def test_authoritative_input_and_validity(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.run1(d)
   self.assertTrue(r["starting_coordinates_are_source_final_coordinates"])
   self.assertFalse(r["spacing_is_hard_minimum"]); self.assertFalse(r["spacing_sweep_used"])
   self.assertEqual((r["final"]["vertex_count"],r["final"]["edge_count"]),(34,44))
   self.assertEqual(r["final"]["coincident_vertex_count"],0)
   self.assertEqual(r["final"]["vertex_on_unrelated_edge_count"],0)
   self.assertTrue(all(m["crossings_after"]<m["crossings_before"] for m in r["accepted_moves"]))
 def test_determinism(self):
  with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
   x=self.run1(a); y=self.run1(b); self.assertEqual(x,y)
   for n in ("iamp_uncross_retry_from_current.svg","iamp_uncross_retry_from_current_report.json"):
    self.assertEqual(hashlib.sha256((Path(a)/n).read_bytes()).digest(),hashlib.sha256((Path(b)/n).read_bytes()).digest())
if __name__=="__main__": unittest.main()
