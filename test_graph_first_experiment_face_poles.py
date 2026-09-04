import json
import tempfile
import unittest
from pathlib import Path

import graph_first_experiment_face_poles as experiment


class FacePoleTests(unittest.TestCase):
    def test_faces_and_poles_are_measurement_only(self):
        base=Path('output/graph_first'); direct=base/'iamp_nearness_rotation_direct_report.json'
        before=json.loads(direct.read_text())['final_coordinates']
        with tempfile.TemporaryDirectory() as directory:
            report=experiment.run(base/'iamp_graph.json',direct,Path(directory))
        after=json.loads(direct.read_text())['final_coordinates']
        self.assertEqual(before,after)
        self.assertEqual(report['vertex_count'],34); self.assertEqual(report['electrical_incidence_count'],44)
        self.assertEqual(report['proper_crossing_count'],0); self.assertEqual(report['face_count'],12)
        self.assertEqual(len(report['selected_faces']),3)
        for face in report['selected_faces']:
            self.assertGreater(face['pole_to_boundary_minimum_clearance'],0)
            self.assertFalse(face['strictly_interior_graph_vertices'])
            self.assertFalse(face['strictly_interior_graph_edges'])

if __name__=='__main__': unittest.main()
