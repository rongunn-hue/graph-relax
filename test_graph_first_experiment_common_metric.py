import hashlib, tempfile, unittest
from pathlib import Path
import graph_first_experiment_common_metric as ex

class Tests(unittest.TestCase):
 def go(self,d):
  b=Path("output/graph_first");return ex.run(b/"iamp_graph.json",b/"iamp_nearness_rotation_direct_report.json",Path(d))
 def test_metric_detector_and_validity(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.go(d);self.assertEqual((r["starting_validation"]["vertices"],r["starting_validation"]["edges"],r["starting_validation"]["proper_crossings"]),(34,44,0));self.assertTrue(all(x["affected"]["vertex"] not in x["affected"]["edge"] for x in r["initial_map"]["node_edge"]));self.assertTrue(all(x["d"]<x["D"] and abs(x["ratio"]-x["d"]/x["D"])<1e-12 for rows in r["initial_map"].values() for x in rows));self.assertTrue(r["audits"]["TPO"]["own_incident_edges_excluded_from_node_edge"]);self.assertTrue(any("component:TPB" in str(x["affected"]) for x in r["initial_map"]["shared_endpoint_overlay"]));self.assertEqual((r["final_validation"]["vertices"],r["final_validation"]["edges"],r["final_validation"]["proper_crossings"],r["final_validation"]["coincidences"],r["final_validation"]["vertex_on_unrelated_edge"]),(34,44,0,0,0))
 def test_annotations_and_determinism(self):
  with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
   x=self.go(a);y=self.go(b);self.assertEqual(x,y)
   names=["iamp_common_metric_direct_start.svg","iamp_common_metric_direct_polygon.svg","iamp_common_metric_direct_nodes.svg","iamp_common_metric_direct_final.svg","iamp_common_metric_direct_report.json"]
   for n in names:self.assertEqual(hashlib.sha256((Path(a)/n).read_bytes()).digest(),hashlib.sha256((Path(b)/n).read_bytes()).digest())
   for n,key in zip(names[:4],["initial_map",None,None,"final_map"]):self.assertIn("EVENT_ANNOTATION_COUNT=",(Path(a)/n).read_text())
if __name__=="__main__":unittest.main()
