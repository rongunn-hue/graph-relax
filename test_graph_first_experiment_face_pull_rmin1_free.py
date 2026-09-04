import json
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path
import graph_first_experiment_face_pull_rmin1_free as experiment

class RestrictedFacePullTests(unittest.TestCase):
    def test_only_topologically_free_vertices_move(self):
        base=Path('output/graph_first');direct=base/'iamp_nearness_rotation_direct_report.json';original=json.loads(direct.read_text())['final_coordinates']
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)
            report=experiment.run(base/'iamp_graph.json',direct,base/'iamp_face_poles_report.json',base/'iamp_face_pull_rmin1_report.json',output)
            root=ET.parse(output/'iamp_face_pull_rmin1_free.svg').getroot()
        self.assertEqual(34,report['vertex_count']);self.assertEqual(44,report['electrical_incidence_count'])
        free=set(report['free_vertices']);anchors=set(report['anchor_vertices']);self.assertTrue(free);self.assertTrue(anchors)
        for record in report['boundary_classification']:
            self.assertEqual(record['classification']=='FREE',record['degree']==2 and not record['additional_incident_edges'])
        for stage in report['stages']:
            self.assertEqual(0,stage['stale_edge_endpoint_count'])
            for node in set(original)-free:self.assertEqual(original[node],stage['coordinates'][node])
        final=next(element for element in root if element.get('id')=='stage-5')
        coordinates=report['stages'][-1]['coordinates']
        circles={element.get('data-vertex'):element for element in final if element.get('data-vertex')}
        lines={element.get('data-edge'):element for element in final if element.get('data-edge')}
        self.assertEqual(set(coordinates),set(circles));self.assertEqual(44,len(lines))
        for node,point in coordinates.items():
            self.assertEqual(float(circles[node].get('cx')),point[0]);self.assertEqual(float(circles[node].get('cy')),point[1])
        for edge,line in lines.items():
            first,second=edge.split('|');self.assertEqual([float(line.get('x1')),float(line.get('y1'))],coordinates[first]);self.assertEqual([float(line.get('x2')),float(line.get('y2'))],coordinates[second])
        self.assertEqual(original,json.loads(direct.read_text())['final_coordinates'])

if __name__=='__main__':unittest.main()
