import json
import unittest
from pathlib import Path

import experiment_4n as e4n
import experiment_4p_active_circular as active
import graph_relax as gr

ROOT=Path(__file__).parent;INPUT=ROOT/"circuit_examples"/"iamp.circuit";OUTPUT=ROOT/"output"


class ActiveCircular4PTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.circuit=gr.parse_circuit(INPUT);cls.graph,cls.embedding,cls.positions,cls.order=e4n.load_start(cls.circuit,OUTPUT)
 def test_fixed_parameters_and_direction_basis(self):
  self.assertEqual((1.10,2.0,2),(active.LENGTH_ALLOWANCE,active.ANGULAR_TOLERANCE_DEGREES,active.PASSES_PER_STEP))
  _,directions=active.movement_directions(self.graph,self.positions,self.order[0])
  self.assertEqual(8,len(directions));self.assertEqual(list(active.DIRECTION_LABELS[1:]),[x[0] for x in directions])
  self.assertTrue(all(abs(sum(v*v for v in unit)-1)<1e-12 for _,unit in directions))
 def test_selection_allows_ten_percent_and_uses_two_degree_band(self):
  vertex="component:U1";current=100.
  rows=[{"local_length":100.,"direction_index":0,"angular":{"angular_error":20.,"minimum_gap":10.}},
        {"local_length":109.,"direction_index":1,"angular":{"angular_error":10.,"minimum_gap":20.}},
        {"local_length":80.,"direction_index":2,"angular":{"angular_error":11.5,"minimum_gap":15.}},
        {"local_length":111.,"direction_index":3,"angular":{"angular_error":0.,"minimum_gap":50.}}]
  self.assertEqual(2,active.select_candidate(self.graph,vertex,rows,current)["direction_index"])
 def test_saved_result_has_twelve_valid_passes(self):
  metrics=json.loads((OUTPUT/"planar_active_circular_metrics.json").read_text())
  validation=json.loads((OUTPUT/"planar_active_circular_validation.json").read_text())
  self.assertEqual(12,len(metrics["passes"]));self.assertTrue(validation["exactly_twelve_passes"])
  self.assertTrue(all(row["validation"]["valid"] for row in metrics["passes"]))
  self.assertTrue(validation["final_complete_validation"]["valid"])


if __name__=="__main__":unittest.main()
