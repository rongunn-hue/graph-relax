import json,random,unittest
from pathlib import Path
import experiment_4n as e4n
import experiment_4p_monte_carlo as mc
import graph_relax as gr
ROOT=Path(__file__).parent;INPUT=ROOT/"circuit_examples"/"iamp.circuit";OUTPUT=ROOT/"output"
class MonteCarlo4PTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.circuit=gr.parse_circuit(INPUT);cls.graph,cls.embedding,cls.positions,cls.order=e4n.load_start(cls.circuit,OUTPUT)
 def test_fixed_sampling_parameters_and_rng(self):
  self.assertEqual((314159,128,.10),(mc.RNG_SEED,mc.SAMPLES_PER_VERTEX,mc.ELITE_FRACTION))
  self.assertEqual(random.Random(314159).random(),random.Random(mc.RNG_SEED).random())
  self.assertEqual(26112,len(self.order)*mc.SAMPLES_PER_VERTEX*len(e4n.STEPS))
 def test_elite_mean_and_small_population(self):
  current=(0.,0.);small=[{"point":(1.,0.)} for _ in range(9)]
  elite,mean=mc.elite_proposal(small,current);self.assertEqual(1,len(elite));self.assertIsNone(mean)
  large=[{"point":(float(i),0.)} for i in range(20)]
  elite,mean=mc.elite_proposal(large,current);self.assertEqual(2,len(elite));self.assertEqual((.5,0.),mean)
 def test_saved_result_has_six_valid_passes(self):
  m=json.loads((OUTPUT/"planar_monte_carlo_compact_metrics.json").read_text());v=json.loads((OUTPUT/"planar_monte_carlo_compact_validation.json").read_text())
  self.assertEqual(6,len(m["passes"]));self.assertEqual(26112,m["total_sampled_candidates"])
  self.assertTrue(all(p["validation"]["valid"] for p in m["passes"]));self.assertTrue(v["final_complete_validation"]["valid"])
if __name__=="__main__":unittest.main()
