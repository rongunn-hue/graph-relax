import tempfile,unittest,math
from pathlib import Path
import graph_first_experiment_face004_directed_opening as experiment

class DirectedOpeningTests(unittest.TestCase):
    def test_only_direction_changes_and_edges_reconnect(self):
        b=Path('output/graph_first')
        with tempfile.TemporaryDirectory() as d:r=experiment.run(b/'iamp_graph.json',b/'iamp_nearness_rotation_direct_report.json',b/'iamp_sequential_local_face_pulls_report.json',b/'iamp_all_faces_report.json',b/'iamp_face004_direction_test_report.json',b/'iamp_crowding_points_report.json',Path(d))
        self.assertEqual(6,r['pre_metrics']['clearance_violations']);self.assertEqual(34,r['pre_metrics']['vertex_count']);self.assertEqual(44,r['pre_metrics']['electrical_incidence_count'])
        for m in r['movements']:self.assertAlmostEqual(m['old_radial_displacement_magnitude'],m['new_directed_displacement_magnitude'],places=12)
        for segment in r['complete_rebuilt_edge_segments']:
            a,b=segment['edge'];self.assertEqual(r['final_coordinates'][a],segment['start']);self.assertEqual(r['final_coordinates'][b],segment['end'])

if __name__=='__main__':unittest.main()
