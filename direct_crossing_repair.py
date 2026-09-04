"""Generic deterministic orchestration for direct endpoint crossing repair.

This module knows only graph vertices, graph edges, current coordinates, and a
single-position geometric orientation solver.  Every proper crossing is always
tested in all four endpoint/pivot orientations; no preferred endpoint exists.
"""

from __future__ import annotations


def crossing_key(crossing):
    return (tuple(crossing["edge_a"]), tuple(crossing["edge_b"]))


def endpoint_pivot_orientations(crossing):
    """Return A/B, B/A, C/D, D/C exactly once and in stable order."""
    first = tuple(crossing["edge_a"])
    second = tuple(crossing["edge_b"])
    return [
        {"moving": first[0], "pivot": first[1], "crossed_edge": second},
        {"moving": first[1], "pivot": first[0], "crossed_edge": second},
        {"moving": second[0], "pivot": second[1], "crossed_edge": first},
        {"moving": second[1], "pivot": second[0], "crossed_edge": first},
    ]


def evaluate_crossing_orientations(crossing, solve_orientation):
    """Calculate at most one position for each of the four orientations."""
    results = []
    for index, orientation in enumerate(endpoint_pivot_orientations(crossing)):
        result = solve_orientation(**orientation)
        positions_tested = int(result.get("positions_evaluated", 0))
        if positions_tested > 1:
            raise AssertionError("an endpoint/pivot orientation evaluated multiple positions")
        result = dict(result)
        result.update({"orientation_index": index, **orientation,
                       "positions_evaluated": positions_tested})
        results.append(result)
    if len(results) > 4:
        raise AssertionError("more than four direct solutions evaluated for one crossing")
    return results


def select_reducing_orientation(results, current_crossing_count):
    """Select deterministically by total crossings, never target removal alone."""
    candidates = []
    for result in results:
        if not result.get("valid"):
            result["selection_status"] = "REJECTED_GEOMETRICALLY_INVALID"
            continue
        after = result.get("resulting_total_crossings")
        if after is None or after >= current_crossing_count:
            result["selection_status"] = "REJECTED_TOTAL_CROSSINGS_NOT_REDUCED"
            continue
        result["selection_status"] = "ELIGIBLE_STRICT_REDUCTION"
        key = (after, abs(result["signed_delta_radians"]), result["displacement"],
               result["moving"], result["pivot"])
        candidates.append((key, result))
    if not candidates:
        return None
    selected = min(candidates, key=lambda item: item[0])[1]
    selected["selection_status"] = "SELECTED"
    return selected


def best_move_across_crossings(crossings, current_crossing_count, solver_factory):
    """Evaluate all four orientations of every crossing and choose globally."""
    all_results = []
    eligible = []
    audits = []
    for crossing in crossings:
        key = crossing_key(crossing)
        results = evaluate_crossing_orientations(crossing, solver_factory(key))
        selected = select_reducing_orientation(results, current_crossing_count)
        audits.append({"crossing": {"edge_a": list(key[0]), "edge_b": list(key[1])},
                       "solutions": results})
        all_results.extend(results)
        if selected is not None:
            selected["target"] = key
            eligible.append(selected)
    if not eligible:
        return None, audits
    selected = min(eligible, key=lambda item: (
        item["resulting_total_crossings"], abs(item["signed_delta_radians"]),
        item["displacement"], item["moving"], item["pivot"]))
    for audit in audits:
        for result in audit["solutions"]:
            if result.get("selection_status") == "SELECTED":
                result["selection_status"] = (
                    "ACCEPTED_GLOBAL_MOVE" if result is selected
                    else "REJECTED_BETTER_GLOBAL_MOVE_AVAILABLE")
            elif result.get("selection_status") == "ELIGIBLE_STRICT_REDUCTION":
                result["selection_status"] = "REJECTED_DETERMINISTIC_SELECTION"
    return selected, audits
