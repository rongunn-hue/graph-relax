import json
import unittest
from pathlib import Path

import graph_first_experiment_cluster02_section_analysis as experiment
import graph_first_experiment_crowding_analysis as crowding

BASE=Path(__file__).parent/"output"/"graph_first"

class Cluster02SectionAnalysisTests(unittest.TestCase):
    def test_analysis_is_measurement_only_and_covers_eight_events(self):
        _,_,coordinates,_,_=crowding.load_inputs(BASE/"iamp_graph.json",BASE/"iamp_nearness_rotation_direct_report.json")
        before=dict(coordinates);points=json.loads((BASE/"iamp_crowding_points_report.json").read_text())
        result=experiment.analyze(points,coordinates)
        self.assertEqual(8,result["violation_count"])
        self.assertEqual(before,coordinates)
        self.assertEqual(set(result["side_groups"]),{"SIDE_A","SIDE_B","ON_AXIS"})

if __name__=="__main__":unittest.main()
