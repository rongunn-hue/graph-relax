import tempfile
import unittest
from pathlib import Path

import experiment_4r as e4r
import graph_relax as gr


ROOT = Path(__file__).parent
INPUT = ROOT / "test_data" / "iamp.circuit"
OUTPUT = ROOT / "output"


class Experiment4RTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.circuit = gr.parse_circuit(INPUT)

    def test_loads_all_exact_4q_embeddings(self):
        for set_id in e4r.SET_IDS:
            document, graph, embedding = e4r.load_and_verify_set(self.circuit, OUTPUT, set_id)
            self.assertEqual(78, graph.number_of_nodes())
            self.assertEqual(106, graph.number_of_edges())
            self.assertEqual(set(graph), set(embedding))

    def test_rotation_preservation_exposes_mixed_reversals(self):
        expected = {"SET001": {"J_PWR", "U1"}, "SET002": {"J_PWR", "U1"},
                    "SET003": {"J_PWR", "U1"}, "SET004": {"J_PWR", "J_R4", "U1"}}
        for set_id in e4r.SET_IDS:
            document, _, _ = e4r.load_and_verify_set(self.circuit, OUTPUT, set_id)
            rows = e4r.component_rotation_analysis(self.circuit, document)
            reversed_refs = {row["component"] for row in rows if row["state"] == "reversed_order"}
            self.assertEqual(expected[set_id], reversed_refs)
            validation = e4r.geometry_validation(self.circuit, document, rows)
            self.assertTrue(validation["mixed_independent_reversals"])
            self.assertFalse(validation["primary_realization_valid"])

    def test_gadget_collapse_is_valid_planar_embedding(self):
        for set_id in e4r.SET_IDS:
            document, _, _ = e4r.load_and_verify_set(self.circuit, OUTPUT, set_id)
            collapsed = e4r.collapse_embedding(self.circuit, document)
            collapsed.check_structure()
            self.assertEqual(34, len(collapsed))
            self.assertEqual(84, collapsed.number_of_edges())

    def test_failure_geometry_has_no_false_connectivity_or_layers(self):
        with tempfile.TemporaryDirectory(prefix="graph-relax-4r-test-") as directory:
            out = Path(directory)
            for set_id in e4r.SET_IDS:
                source = OUTPUT / f"constrained_planarization_embedding_{set_id}.json"
                (out / source.name).write_bytes(source.read_bytes())
            result = e4r.run(INPUT, out)
            for set_id in e4r.SET_IDS:
                geometry = __import__('json').loads(
                    (out / f"planar_realization_{set_id}_geometry.json").read_text())
                self.assertEqual([], geometry["layer0_polylines"])
                self.assertEqual([], geometry["layer1_polylines"])
                self.assertEqual("FAILED", result["results"][set_id])

    def test_comparison_ranking_and_deterministic_serialization(self):
        def prepare():
            out = Path(tempfile.mkdtemp(prefix="graph-relax-4r-repeat-"))
            for set_id in e4r.SET_IDS:
                source = OUTPUT / f"constrained_planarization_embedding_{set_id}.json"
                (out / source.name).write_bytes(source.read_bytes())
            return out
        first, second = prepare(), prepare()
        result_a, result_b = e4r.run(INPUT, first), e4r.run(INPUT, second)
        self.assertEqual(result_a, result_b)
        comparison = __import__('json').loads((first / "planar_realization_comparison.json").read_text())
        self.assertEqual(list(e4r.SET_IDS), comparison["ranking"])


if __name__ == "__main__":
    unittest.main()
