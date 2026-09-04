import hashlib,tempfile,unittest
from pathlib import Path
import graph_first_experiment_global_node_density as experiment
class DensityTests(unittest.TestCase):
 def run1(self,d):
  b=Path("output/graph_first");return experiment.run(b/"iamp_graph.json",b/"iamp_corridor_separation_report.json",Path(d))
 def test_partition_and_invariants(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.run1(d);self.assertEqual((r["vertex_count"],r["electrical_incidence_count"],r["proper_crossings"]),(34,44,0));self.assertTrue(r["coordinates_unchanged"]);self.assertEqual(len(r["ranked_vertices"]),34);self.assertEqual(r["negative_cell_count"],0);self.assertAlmostEqual(r["sum_clipped_cell_areas"],r["domain"]["area"],7)
 def test_ratios_and_determinism(self):
  with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
   r=self.run1(a);self.run1(b)
   for x in r["ranked_vertices"]:self.assertAlmostEqual(x["crowding_ratio"],r["target_area"]/x["voronoi_area"],12)
   for n in ("iamp_global_node_density.svg","iamp_global_node_density_report.json"):self.assertEqual(hashlib.sha256((Path(a)/n).read_bytes()).digest(),hashlib.sha256((Path(b)/n).read_bytes()).digest())
if __name__=="__main__":unittest.main()
