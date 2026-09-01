import json,math,tempfile,unittest
from pathlib import Path
import experiment_4n as e4n
import experiment_4q_tent as tent
import graph_relax as gr
ROOT=Path(__file__).parent;INPUT=ROOT/"circuit_examples"/"iamp.circuit";OUTPUT=ROOT/"output"
class Tent4QTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.circuit=gr.parse_circuit(INPUT);cls.graph,cls.embedding,cls.initial,cls.order=e4n.load_start(cls.circuit,OUTPUT);cls.model=tent.tent_model(cls.initial)
 def test_zero_reproduces_and_rays_are_preserved(self):
  center,cables,directions,centered,minimum,owner=self.model
  zero=tent.positions_at_height(self.initial,center,cables,directions,centered,0.)
  for v in self.initial:
   self.assertLess(math.dist(zero[v],self.initial[v]),2e-6)
  moved=tent.positions_at_height(self.initial,center,cables,directions,centered,.5*minimum)
  for v in directions:
   dx,dy=moved[v][0]-center[0],moved[v][1]-center[1];ux,uy=directions[v]
   self.assertLess(abs(dx*uy-dy*ux),2e-6)
 def test_radius_decreases_and_direct_computation_has_no_drift(self):
  center,cables,directions,centered,minimum,_=self.model
  a=tent.positions_at_height(self.initial,center,cables,directions,centered,.2*minimum)
  b=tent.positions_at_height(self.initial,center,cables,directions,centered,.7*minimum)
  self.assertTrue(all(math.dist(b[v],center)<=math.dist(a[v],center)+2e-6 for v in directions))
  self.assertEqual(b,tent.positions_at_height(self.initial,center,cables,directions,centered,.7*minimum))
 def test_saved_selected_geometry_is_valid(self):
  v=json.loads((OUTPUT/"planar_tent_compact_validation.json").read_text());self.assertTrue(v["selected_complete_validation"]["valid"])
 def test_repeated_execution_is_deterministic(self):
  def prep():
   d=Path(tempfile.mkdtemp(prefix="graph-relax-tent-"))
   for n in ("node_planar_embedding.json","planar_xy_compact_coordinates.json"):(d/n).write_bytes((OUTPUT/n).read_bytes())
   return d
  a,b=prep(),prep();self.assertEqual(tent.run(INPUT,a)["hashes"],tent.run(INPUT,b)["hashes"])
if __name__=="__main__":unittest.main()
