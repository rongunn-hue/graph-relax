import json,tempfile,unittest
from pathlib import Path
import experiment_4j as e4j, graph_relax as gr

ROOT=Path(__file__).parent;INPUT=ROOT/"circuit_examples"/"iamp.circuit";OUTPUT=ROOT/"output"

class Experiment4JTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.circuit=gr.parse_circuit(INPUT);cls.loaded=e4j.load_stored(cls.circuit,OUTPUT)
 def test_exact_stored_baseline_and_sequence(self):
  self.assertEqual(tuple(range(72,0,-1)),e4j.SCALES)
  result=e4j.verify_baseline(self.circuit,self.loaded[2],self.loaded[3],self.loaded[4],self.loaded[5],self.loaded[6],OUTPUT)
  self.assertTrue(result[1]["valid"]);self.assertAlmostEqual(34084.492436,result[2]["total_connection_length"],places=5)
 def test_uniform_scaling_and_dimensions(self):
  a=e4j.geometry_at_scale(self.circuit,self.loaded[5],self.loaded[6],72);b=e4j.geometry_at_scale(self.circuit,self.loaded[5],self.loaded[6],36)
  ref=self.circuit.refs[-1];self.assertEqual((a[0][ref][0]/2,a[0][ref][1]/2),b[0][ref]);self.assertEqual((gr.WIDTH,gr.HEIGHT),(60.0,40.0))
 def test_complete_sweep_selection_and_failure(self):
  records,geometries,smallest,failed=e4j.sweep(self.circuit,self.loaded[3],self.loaded[4],self.loaded[5],self.loaded[6])
  self.assertEqual(72,len(records));self.assertEqual(52,smallest);self.assertEqual(51,failed["scale"]);self.assertEqual(5,failed["clearance_violations"]);self.assertFalse(failed["fully_valid"])
  self.assertTrue(all(r["cyclic_order_orientation"]=="global_reflection" for r in records if r["fully_valid"]))
 def test_failure_classification(self):
  record={"component_overlaps":0,"clearance_violations":5,"different_net_crossings":0,"unrelated_component_body_intersections":0,"electrical_incidence_exact":True,"cyclic_order_orientation":"global_reflection"}
  self.assertEqual(["clearance violations: 5"],e4j.failure_reasons(record))
 def test_serialization_repeatability(self):
  def prep():
   d=Path(tempfile.mkdtemp(prefix="graph-relax-4j-test-"))
   for n in ("node_planar_embedding.json","planar_node_initial_coordinates.json","planar_node_initial_validation.json","planar_node_initial_metrics.json"):(d/n).write_bytes((OUTPUT/n).read_bytes())
   return d
  a,b=prep(),prep();self.assertEqual(e4j.run(INPUT,a)["hashes"],e4j.run(INPUT,b)["hashes"]);self.assertEqual(72,len(json.loads((a/"planar_scale_sweep.json").read_text())["records"]))

if __name__=="__main__":unittest.main()
