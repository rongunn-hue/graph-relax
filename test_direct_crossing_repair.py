import unittest

import direct_crossing_repair as repair


CROSSING = {"edge_a": ["A", "B"], "edge_b": ["C", "D"]}


def result(valid, crossings, rotation, displacement=1.0, positions=1):
    return {"valid": valid, "resulting_total_crossings": crossings,
            "signed_delta_radians": rotation, "displacement": displacement,
            "positions_evaluated": positions}


class GenericDirectCrossingRepairTests(unittest.TestCase):
    def evaluate(self, outcomes):
        iterator = iter(outcomes)
        values = repair.evaluate_crossing_orientations(
            CROSSING, lambda **unused: next(iterator))
        return values, repair.select_reducing_orientation(values, 3)

    def test_first_endpoint_can_fail_while_opposite_succeeds(self):
        values, selected = self.evaluate([
            result(False, None, 0.1), result(True, 2, 0.2),
            result(False, None, 0.3), result(False, None, 0.4)])
        self.assertEqual(("B", "A"), (selected["moving"], selected["pivot"]))

    def test_first_can_create_crossing_while_opposite_succeeds(self):
        values, selected = self.evaluate([
            result(True, 3, 0.1), result(True, 2, 0.2),
            result(False, None, 0.3), result(False, None, 0.4)])
        self.assertEqual("REJECTED_TOTAL_CROSSINGS_NOT_REDUCED", values[0]["selection_status"])
        self.assertEqual("B", selected["moving"])

    def test_both_orientations_can_fail(self):
        _, selected = self.evaluate([
            result(False, None, 0.1), result(False, None, 0.2),
            result(False, None, 0.3), result(False, None, 0.4)])
        self.assertIsNone(selected)

    def test_both_orientations_can_succeed(self):
        values, selected = self.evaluate([
            result(True, 2, 0.3), result(True, 1, 0.4),
            result(False, None, 0.5), result(False, None, 0.6)])
        self.assertEqual("B", selected["moving"])
        self.assertEqual(2, sum(value["resulting_total_crossings"] < 3
                                for value in values if value["valid"]))

    def test_deterministic_selection_for_multiple_successes(self):
        _, selected = self.evaluate([
            result(True, 1, 0.3), result(True, 1, -0.2),
            result(True, 1, 0.2, 2.0), result(False, None, 0.1)])
        self.assertEqual(("B", "A"), (selected["moving"], selected["pivot"]))

    def test_no_orientation_evaluates_more_than_one_position(self):
        values, _ = self.evaluate([result(True, 2, 0.1) for _ in range(4)])
        self.assertEqual(4, len(values))
        self.assertEqual(1, max(value["positions_evaluated"] for value in values))
        with self.assertRaises(AssertionError):
            repair.evaluate_crossing_orientations(
                CROSSING, lambda **unused: result(True, 2, 0.1, positions=2))


if __name__ == "__main__":
    unittest.main()
