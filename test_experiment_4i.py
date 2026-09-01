import json,tempfile,unittest
from pathlib import Path
import experiment_4i as e4i, experiment_4p as e4p, graph_relax as gr

ROOT=Path(__file__).parent; INPUT=ROOT/"circuit_examples"/"iamp.circuit"; OUTPUT=ROOT/"output"

class Experiment4ITests(unittest.TestCase):
 def setUp(self): self.circuit=gr.parse_circuit(INPUT)
 def test_graph_embedding_and_stable_relabeling(self):
  graph,source,integer,mapping,reverse,positions=e4i.stable_embedding(self.circuit)
  self.assertEqual((34,44),(len(graph),graph.number_of_edges()));self.assertEqual(list(range(34)),sorted(reverse));integer.check_structure()
  self.assertEqual((0,64,0,21),(min(x for x,y in positions.values()),max(x for x,y in positions.values()),min(y for x,y in positions.values()),max(y for x,y in positions.values())))
 def test_ray_rectangle_attachment(self):
  self.assertEqual((30.0,0.0),e4i.ray_rectangle_attachment((0,0),(100,0)))
  self.assertEqual((0.0,20.0),e4i.ray_rectangle_attachment((0,0),(0,100)))
 def test_scale_72_and_free_attachments_validate(self):
  graph,source,integer,mapping,reverse,positions=e4i.stable_embedding(self.circuit)
  realization=e4i.realize(self.circuit,reverse,positions,72); validation=e4i.validate_geometry(self.circuit,graph,source,*realization)
  self.assertTrue(validation["valid"]);self.assertEqual((0,0,0,0),(validation["component_overlaps"],validation["clearance_violations"],validation["different_net_crossings"],validation["unrelated_component_body_intersections"]))
  self.assertEqual("global_reflection",validation["cyclic_order_orientation"])
  self.assertEqual(44,len(realization[3]));self.assertTrue(all(e4i.perimeter_parameter(realization[2][t],realization[0][t.rsplit('.',1)[0]])>=0 for t in realization[2]))
 def test_rejects_inconsistent_local_reversal(self):
  self.assertTrue(e4i.cyclic_equal([2,1,3],[1,3,2]));self.assertFalse(e4i.cyclic_equal([1,3,2],[1,2,3]))
 def test_serialization_and_repeated_run(self):
  def prep():
   d=Path(tempfile.mkdtemp(prefix="graph-relax-4i-test-"))
   for n in ("closeness_metrics.json","closeness_hv_annealed_metrics.json","junction_untangled_metrics.json"):(d/n).write_bytes((OUTPUT/n).read_bytes())
   return d
  a,b=prep(),prep();self.assertEqual(e4i.run(INPUT,a)["hashes"],e4i.run(INPUT,b)["hashes"])
  self.assertTrue(json.loads((a/"planar_node_initial_validation.json").read_text())["valid"])

if __name__=="__main__":unittest.main()
