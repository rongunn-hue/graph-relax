import hashlib,tempfile,unittest
from pathlib import Path
import graph_first_experiment_outward_boundary_wave as experiment

class OutwardBoundaryWaveTests(unittest.TestCase):
 def run_once(self,d):
  b=Path('output/graph_first');return experiment.run(b/'iamp_graph.json',b/'iamp_corridor_separation_report.json',Path(d))
 def test_frozen_sequences_and_exterior_exclusion(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.run_once(d);self.assertTrue(r['boundary_sequences']['constructed_from_original_embedding']);self.assertTrue(r['boundary_sequences']['frozen_before_movement']);self.assertFalse(r['sequence_assertions']['exterior_processed']);self.assertFalse(r['sequence_assertions']['duplicate_boundary_processing'])
 def test_half_gap_anchor_x_and_one_move(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.run_once(d);before=r['stages'][0]['coordinates'];seen=set()
   for stage in r['operations']:
    self.assertAlmostEqual(stage['half_gap'],stage['full_gap']/2,12)
    self.assertAlmostEqual(stage['actual_move'],min(stage['max_valid_move'],stage['half_gap']),12)
    state=next(x for x in r['stages'] if x.get('operation')==stage)
    for anchor in stage['anchors']:self.assertEqual(before[anchor],state['coordinates'][anchor])
    for v in stage['eligible_unmoved_vertices']:
     self.assertNotIn(v,seen);seen.add(v);self.assertEqual(before[v][0],state['coordinates'][v][0]);self.assertLess(state['coordinates'][v][1],before[v][1]) if stage['side']=='top' else self.assertGreater(state['coordinates'][v][1],before[v][1])
 def test_connectivity_and_determinism(self):
  with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
   x=self.run_once(a);self.run_once(b);self.assertEqual((x['vertex_count'],x['electrical_incidence_count']),(34,44));self.assertTrue(x['sequence_assertions']['final_edge_reconnection_valid'])
   for n in ('iamp_outward_boundary_wave.svg','iamp_outward_boundary_wave_report.json'):self.assertEqual(hashlib.sha256((Path(a)/n).read_bytes()).digest(),hashlib.sha256((Path(b)/n).read_bytes()).digest())
 def test_limited_bottom_does_not_stop_or_undo_top(self):
  with tempfile.TemporaryDirectory() as d:
   r=self.run_once(d);top,bottom=r['operations'];self.assertTrue(top['limiting_validity_search']['cap_is_valid']);self.assertFalse(bottom['limiting_validity_search']['cap_is_valid']);self.assertLess(bottom['actual_move'],bottom['half_gap'])
   self.assertEqual(r['side_progress']['top'],{'available_transitions':1,'processed_transitions':1,'accepted_transitions':1,'rejected_transitions':0,'exhausted':True})
   self.assertEqual(r['side_progress']['bottom'],{'available_transitions':1,'processed_transitions':1,'accepted_transitions':1,'rejected_transitions':0,'exhausted':True})
   self.assertNotEqual(r['stages'][0]['coordinates']['net:EA'],r['stages'][-1]['coordinates']['net:EA'])
   self.assertEqual(r['stages'][1]['coordinates']['net:EA'],r['stages'][-1]['coordinates']['net:EA'])
if __name__=='__main__':unittest.main()
