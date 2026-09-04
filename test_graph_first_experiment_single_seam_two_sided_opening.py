import json,tempfile,unittest
from pathlib import Path
import graph_first_experiment_single_seam_two_sided_opening as experiment

class SingleSeamAuditTests(unittest.TestCase):
    def test_no_ambiguous_or_one_face_event_is_moved(self):
        b=Path('output/graph_first');seq=b/'iamp_sequential_local_face_pulls_report.json';before=next(s for s in json.loads(seq.read_text())['stages'] if s['stage']=='AFTER face:003')['coordinates']
        with tempfile.TemporaryDirectory() as d:r=experiment.run(b/'iamp_graph.json',b/'iamp_nearness_rotation_direct_report.json',seq,Path(d))
        self.assertFalse(r['primitive_tested']);self.assertIsNone(r['selected_event']);self.assertTrue(r['coordinates_unchanged']);self.assertEqual(6,r['clearance_violations']);self.assertEqual(0,r['proper_crossings']);self.assertTrue(all(e['eligibility']=='INELIGIBLE' for e in r['event_eligibility_audit']))

if __name__=='__main__':unittest.main()
