import json,tempfile,unittest
from pathlib import Path
import graph_first_experiment_face007_differential_owned as experiment

class DifferentialOwnedTests(unittest.TestCase):
    def test_deterministic_exclusive_simultaneous_pull(self):
        b=Path('output/graph_first');args=(b/'iamp_graph.json',b/'iamp_nearness_rotation_direct_report.json',b/'iamp_all_faces_report.json',b/'iamp_face_poles_report.json',b/'iamp_crowding_points_report.json',b/'iamp_face_pull_rmin1_free_report.json')
        with tempfile.TemporaryDirectory() as d:
            first=experiment.run(*args,Path(d));second=experiment.run(*args,Path(d))
        self.assertEqual(first['owner_face_by_vertex'],second['owner_face_by_vertex']);self.assertEqual(first['final_coordinates'],second['final_coordinates'])
        members=[v for values in first['pull_cell_vertices'].values() for v in values];self.assertEqual(len(members),len(set(members)))
        self.assertEqual(0,sum(first['ownership_validation'].values()));self.assertEqual(34,first['vertex_count']);self.assertEqual(44,first['electrical_incidence_count'])
        for segment in first['complete_rebuilt_edge_segments']:
            a,b=segment['edge'];self.assertEqual(first['final_coordinates'][a],segment['start']);self.assertEqual(first['final_coordinates'][b],segment['end'])

if __name__=='__main__':unittest.main()
