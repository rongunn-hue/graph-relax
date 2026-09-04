import json,tempfile,unittest,xml.etree.ElementTree as ET
from pathlib import Path
import graph_first_experiment_face_pull_three as experiment

class ThreeFacePullTests(unittest.TestCase):
    def test_topology_only_selection_reconnection_and_cumulative_output(self):
        base=Path('output/graph_first');direct=base/'iamp_nearness_rotation_direct_report.json';original=json.loads(direct.read_text())['final_coordinates']
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory);report=experiment.run(base/'iamp_graph.json',direct,base/'iamp_face_poles_report.json',output)
            root=ET.parse(output/'iamp_face_pull_three_cumulative.svg').getroot();final_group=next(x for x in root if x.get('id')=='stage-3')
        self.assertEqual(34,report['vertex_count']);self.assertEqual(44,report['electrical_incidence_count'])
        for data in report['topological_classification'].values():
            for record in data['classification']:
                self.assertEqual(record['classification']=='FREE',record['degree']==2 and len(record['face_boundary_edges'])==2 and not record['additional_incident_edges'])
        for stage in report['independent_tests']+report['cumulative_stages']:
            self.assertEqual(0,stage['stale_edge_endpoint_count'])
        self.assertEqual(original,json.loads(direct.read_text())['final_coordinates'])
        self.assertEqual(44,len([x for x in final_group if x.get('data-edge')]))

if __name__=='__main__':unittest.main()
