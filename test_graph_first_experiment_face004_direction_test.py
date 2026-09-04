import json,tempfile,unittest,math
from pathlib import Path
import graph_first_experiment_face004_direction_test as experiment

class Face004DirectionTests(unittest.TestCase):
    def test_exact_reversal_same_source_same_magnitude_and_reconnection(self):
        b=Path('output/graph_first')
        with tempfile.TemporaryDirectory() as d:r=experiment.run(b/'iamp_graph.json',b/'iamp_nearness_rotation_direct_report.json',b/'iamp_sequential_local_face_pulls_report.json',b/'iamp_crowding_points_report.json',Path(d))
        self.assertTrue(r['case_a_reproduces_sequential_result']);self.assertTrue(r['displacement_magnitudes_equal'])
        pre=r['pre_coordinates'];pole=r['pole_coordinate']
        for v in r['movable_vertices']:
            da=[r['case_a']['coordinates'][v][i]-pre[v][i] for i in (0,1)];db=[r['case_b']['coordinates'][v][i]-pre[v][i] for i in (0,1)]
            self.assertAlmostEqual(da[0],-db[0],places=12);self.assertAlmostEqual(da[1],-db[1],places=12)
        for case in ('case_a','case_b'):
            m=r[case]['metrics'];self.assertEqual(34,m['vertex_count']);self.assertEqual(44,m['electrical_incidence_count']);self.assertEqual(0,m['stale_edge_endpoint_count'])

if __name__=='__main__':unittest.main()
