import json,tempfile,unittest,xml.etree.ElementTree as ET
from pathlib import Path
import graph_first_experiment_sequential_local_face_pulls as experiment

class SequentialLocalFacePullTests(unittest.TestCase):
    def test_persistent_ownership_one_move_and_reconnected_final_svg(self):
        b=Path('output/graph_first');direct=b/'iamp_nearness_rotation_direct_report.json';before=json.loads(direct.read_text())['final_coordinates']
        with tempfile.TemporaryDirectory() as d:
            out=Path(d);r=experiment.run(b/'iamp_graph.json',direct,b/'iamp_all_faces_report.json',b/'iamp_crowding_points_report.json',b/'iamp_face_pull_rmin1_free_report.json',out);root=ET.parse(out/'iamp_sequential_local_face_pulls.svg').getroot()
        owners=r['persistent_owner_face_by_vertex'];members=[v for values in r['persistent_pull_cell_vertices'].values() for v in values]
        self.assertEqual(len(members),len(set(members)));self.assertEqual(set(owners),set(r['movement_stage_by_vertex']));self.assertEqual(34,r['vertex_count']);self.assertEqual(44,r['electrical_incidence_count'])
        final=next(x for x in root if x.get('id')==f"stage-{len(r['stages'])-1}");self.assertEqual(44,len([x for x in final if x.get('data-edge')]))
        self.assertEqual(before,json.loads(direct.read_text())['final_coordinates'])

if __name__=='__main__':unittest.main()
