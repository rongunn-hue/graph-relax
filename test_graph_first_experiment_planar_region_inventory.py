import json,tempfile,unittest
from pathlib import Path
import graph_first_experiment_planar_region_inventory as experiment

class RegionInventoryTests(unittest.TestCase):
    def test_darts_faces_adjacency_and_no_mutation(self):
        base=Path('output/graph_first');direct=base/'iamp_nearness_rotation_direct_report.json';before=json.loads(direct.read_text())['final_coordinates']
        with tempfile.TemporaryDirectory() as directory:r=experiment.run(base/'iamp_graph.json',direct,base/'iamp_face_poles_report.json',base/'iamp_crowding_points_report.json',Path(directory))
        self.assertEqual(34,r['vertex_count']);self.assertEqual(44,r['electrical_incidence_count']);self.assertEqual(0,r['proper_crossing_count']);self.assertEqual(12,r['total_faces']);self.assertEqual(11,r['bounded_faces'])
        self.assertEqual(1,sum(f['face_type']=='EXTERIOR' for f in r['faces']));self.assertEqual(88,sum(len(f['ordered_boundary_darts']) for f in r['faces']))
        self.assertEqual(44,len(r['graph_edge_face_sides']));self.assertEqual(3,len(r['known_region_matches']))
        for face in r['faces']:
            self.assertEqual(set(face['distinct_boundary_vertices']),set(face['ordered_boundary_walk']))
        self.assertEqual(before,json.loads(direct.read_text())['final_coordinates'])

if __name__=='__main__':unittest.main()
