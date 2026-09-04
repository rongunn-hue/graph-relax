import json
import tempfile
import unittest
from pathlib import Path

import graph_first_experiment_face_pull_rmin1 as experiment

class FacePullTests(unittest.TestCase):
    def test_direct_snapshots_move_only_boundary_and_scale_area(self):
        base=Path('output/graph_first');direct=base/'iamp_nearness_rotation_direct_report.json'
        original=json.loads(direct.read_text())['final_coordinates']
        with tempfile.TemporaryDirectory() as directory:
            report=experiment.run(base/'iamp_graph.json',direct,base/'iamp_face_poles_report.json',Path(directory))
        self.assertEqual(34,report['vertex_count']);self.assertEqual(44,report['electrical_incidence_count'])
        boundary=set(report['face_boundary_vertices'])
        for stage in report['stages']:
            self.assertAlmostEqual(stage['s']**2,stage['face_area_ratio'],places=12)
            for node in set(original)-boundary:self.assertEqual(original[node],stage['coordinates'][node])
        self.assertEqual(original,json.loads(direct.read_text())['final_coordinates'])

if __name__=='__main__':unittest.main()
