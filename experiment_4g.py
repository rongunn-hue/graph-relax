#!/usr/bin/env python3
"""Experiment 4G: deterministic diagnosis of frozen 4F crossings."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import experiment_4e as e4e
import experiment_4f as e4f
import graph_relax as gr


SOURCE_TOLERANCE = 2e-5
POINT_TOLERANCE = 1e-7


def load_frozen_state(circuit, output_dir):
    document = json.loads((output_dir / "closeness_hv_annealed_coordinates.json").read_text())
    components = document["components"]
    if set(components) != set(circuit.refs):
        raise ValueError("stored 4F component set does not match input")
    return {ref: (float(components[ref]["x"]), float(components[ref]["y"]))
            for ref in circuit.refs}


def verify_frozen_state(circuit, positions, output_dir):
    actual = gr.evaluate(circuit, positions)
    expected = json.loads((output_dir / "closeness_hv_annealed_metrics.json").read_text())[
        "states"]["experiment_4f_annealed_and_quenched"]
    if actual["net_crossings"] != 7 or actual["net_crossings"] != expected["net_crossings"]:
        raise ValueError(f"frozen 4F crossing count is {actual['net_crossings']}, expected 7")
    if actual["component_overlaps"] != 0 or actual["component_overlaps"] != expected["component_overlaps"]:
        raise ValueError("frozen 4F state does not reproduce zero overlaps")
    scalar_paths = (("total_connection_length",), ("objective",),
                    ("energy", "length"), ("energy", "overlap"),
                    ("energy", "clearance"), ("energy", "crossing"),
                    ("energy", "crowding"))
    for path in scalar_paths:
        left, right = actual, expected
        for part in path:
            left, right = left[part], right[part]
        if abs(left - right) > SOURCE_TOLERANCE:
            raise ValueError(f"frozen 4F metric {'.'.join(path)} exceeds serialization tolerance")
    return actual, expected


def intersection_point(a, b, c, d):
    dx1, dy1 = b[0] - a[0], b[1] - a[1]
    dx2, dy2 = d[0] - c[0], d[1] - c[1]
    denominator = dx1 * dy2 - dy1 * dx2
    if abs(denominator) < 1e-15:
        raise ValueError("proper crossing unexpectedly has parallel segments")
    t = ((c[0] - a[0]) * dy2 - (c[1] - a[1]) * dx2) / denominator
    return a[0] + t * dx1, a[1] + t * dy1


def segment_identity(net, terminal):
    return net, terminal


def crossing_identity(first, second):
    return tuple(sorted((segment_identity(first[0], first[1]),
                         segment_identity(second[0], second[1]))))


def crossing_map(circuit, positions):
    segments = gr.net_segments(circuit, gr.terminal_positions(circuit, positions))
    result = {}
    for index, first in enumerate(segments):
        for second in segments[index + 1:]:
            if first[0] == second[0] or not gr.proper_intersection(first[2], first[3], second[2], second[3]):
                continue
            key = crossing_identity(first, second)
            ordered = sorted((first, second), key=lambda item: (item[0], item[1]))
            result[key] = {"segments": ordered,
                           "intersection": intersection_point(first[2], first[3], second[2], second[3])}
    return result


def net_components(circuit):
    return {net.name: sorted({terminal.rsplit(".", 1)[0] for terminal in net.terminals})
            for net in circuit.nets}


def enumerate_crossings(circuit, positions):
    memberships = net_components(circuit)
    records = []
    for key, item in crossing_map(circuit, positions).items():
        first, second = item["segments"]
        point = item["intersection"]
        immediate = sorted({first[1].rsplit(".", 1)[0], second[1].rsplit(".", 1)[0]})
        relevant = sorted(set(memberships[first[0]]) | set(memberships[second[0]]))
        record = {
            "identity": [list(part) for part in key],
            "net_a": first[0], "net_b": second[0],
            "intersection": {"x": gr.rounded(point[0]), "y": gr.rounded(point[1])},
            "segments": [segment_document(first), segment_document(second)],
            "immediate_segment_anchor_components": immediate,
            "net_a_components": memberships[first[0]], "net_b_components": memberships[second[0]],
            "relevant_components": relevant,
        }
        records.append((point, first[0], second[0], first[1], second[1], record))
    records.sort(key=lambda row: (row[0][1], row[0][0], *row[1:5]))
    result = []
    for index, row in enumerate(records, 1):
        record = row[-1]
        record["id"] = f"X{index:03d}"
        result.append(record)
    return result


def segment_document(segment):
    net, terminal, start, end = segment
    return {"net": net, "segment_identity": [net, terminal],
            "type": "terminal_to_centroid", "terminal": terminal,
            "terminal_component": terminal.rsplit(".", 1)[0],
            "endpoints": [
                {"id": terminal, "kind": "component_terminal", "component": terminal.rsplit(".", 1)[0],
                 "x": gr.rounded(start[0]), "y": gr.rounded(start[1])},
                {"id": f"centroid:{net}", "kind": "derived_net_centroid", "component": None,
                 "x": gr.rounded(end[0]), "y": gr.rounded(end[1])}],
            "length": gr.rounded(math.dist(start, end))}


def step_values(initial, minimum):
    values, step = [], initial
    while True:
        values.append(step)
        if step <= minimum:
            return values
        step = max(minimum, step / 2.0)


def move_tie(move):
    if move["move_type"] == "slide":
        return 0, move["components"][0], "", 0 if move["direction"] < 0 else 1
    return 1, move["components"][0], move["components"][1], 0


def candidate_positions(baseline, axis, move_type, components, amount=None):
    result = dict(baseline)
    coordinate = 0 if axis == "X" else 1
    if move_type == "slide":
        ref = components[0]
        point = list(result[ref])
        point[coordinate] += amount
        result[ref] = tuple(point)
    else:
        first, second = components
        a, b = list(result[first]), list(result[second])
        a[coordinate], b[coordinate] = b[coordinate], a[coordinate]
        result[first], result[second] = tuple(a), tuple(b)
    return result


def candidate_cache_record(circuit, baseline, axis, move_type, components, amount=None):
    positions = candidate_positions(baseline, axis, move_type, components, amount)
    if e4e.has_overlap(circuit, positions):
        return None
    metrics = gr.evaluate(circuit, positions)
    crossings = crossing_map(circuit, positions)
    coordinate = 0 if axis == "X" else 1
    old = [baseline[ref][coordinate] for ref in components]
    new = [positions[ref][coordinate] for ref in components]
    move = {"axis": axis, "move_type": move_type, "components": list(components),
            "old_coordinates": old, "new_coordinates": new,
            "direction": 0 if amount is None else (-1 if amount < 0 else 1),
            "total_displacement": sum(abs(a - b) for a, b in zip(old, new))}
    return {"move": move, "metrics": metrics, "crossings": crossings}


def diagnostic_axis_search(circuit, baseline, crossing, axis, baseline_crossings, cache):
    relevant = crossing["relevant_components"]
    initial = e4e.INITIAL_X_STEP if axis == "X" else e4f.INITIAL_Y_STEP
    minimum = e4e.MINIMUM_X_STEP if axis == "X" else e4f.MINIMUM_Y_STEP
    target = tuple(tuple(part) for part in crossing["identity"])
    candidates = []
    for ref in relevant:
        for step in step_values(initial, minimum):
            for direction in (-1, 1):
                cache_key = (axis, "slide", ref, direction * step)
                if cache_key not in cache:
                    cache[cache_key] = candidate_cache_record(
                        circuit, baseline, axis, "slide", (ref,), direction * step)
                if cache[cache_key] is not None and target not in cache[cache_key]["crossings"]:
                    candidates.append(cache[cache_key])
    for index, first in enumerate(relevant):
        for second in relevant[index + 1:]:
            eligible = (e4e.vertical_spans_overlap(baseline[first], baseline[second]) if axis == "X"
                        else e4f.horizontal_spans_overlap(baseline[first], baseline[second]))
            if not eligible:
                continue
            cache_key = (axis, "swap", first, second)
            if cache_key not in cache:
                cache[cache_key] = candidate_cache_record(circuit, baseline, axis, "swap", (first, second))
            if cache[cache_key] is not None and target not in cache[cache_key]["crossings"]:
                candidates.append(cache[cache_key])
    if not candidates:
        return None
    before_keys = set(baseline_crossings)
    def rank(item):
        new = set(item["crossings"]) - before_keys
        return (item["metrics"]["net_crossings"], len(new),
                item["metrics"]["total_connection_length"], item["move"]["total_displacement"],
                move_tie(item["move"]))
    best = min(candidates, key=rank)
    after_keys = set(best["crossings"])
    new_keys = sorted(after_keys - before_keys)
    eliminated = sorted((before_keys - after_keys) - {target})
    return {"move": best["move"], "target_crossing_eliminated": True,
            "qualifies_without_crossing_increase": best["metrics"]["net_crossings"] <= len(before_keys),
            "total_crossings_before": len(before_keys),
            "total_crossings_after": best["metrics"]["net_crossings"],
            "new_crossings": identities_document(new_keys),
            "eliminated_non_target_crossings": identities_document(eliminated),
            "connection_length_before": gr.rounded(gr.evaluate(circuit, baseline)["total_connection_length"]),
            "connection_length_after": gr.rounded(best["metrics"]["total_connection_length"]),
            "component_overlaps": best["metrics"]["component_overlaps"],
            "clearance_violations": best["metrics"]["clearance_violations"],
            "_crossings": best["crossings"], "_rank": rank(best)}


def identities_document(keys):
    return [[list(first), list(second)] for first, second in keys]


def classify(x_result, y_result):
    xq = x_result is not None and x_result["qualifies_without_crossing_increase"]
    yq = y_result is not None and y_result["qualifies_without_crossing_increase"]
    if xq or yq:
        if xq and (not yq or x_result["_rank"] <= y_result["_rank"]):
            return "X", "a horizontal move is the best non-crossing-increasing target elimination", x_result
        return "Y", "a vertical move is the best non-crossing-increasing target elimination", y_result
    available = [result for result in (x_result, y_result) if result is not None]
    if available:
        chosen = min(available, key=lambda result: (result["_rank"], 0 if result is x_result else 1))
        return "XY", "single-axis elimination exists only with an increase in total crossings", chosen
    return "STRUCTURAL", "no tested single-axis slide or swap eliminates the target", None


def interaction_codes(original, candidate):
    if candidate is None:
        return {record["id"]: "N/A" for record in original}
    after = candidate["_crossings"]
    after_net_pairs = {tuple(sorted((key[0][0], key[1][0]))) for key in after}
    result = {}
    for record in original:
        key = tuple(tuple(part) for part in record["identity"])
        if key in after:
            old = record["intersection"]
            new = after[key]["intersection"]
            result[record["id"]] = ("P" if math.dist((old["x"], old["y"]), new) <= POINT_TOLERANCE else "R")
        elif tuple(sorted((record["net_a"], record["net_b"]))) in after_net_pairs:
            result[record["id"]] = "R"
        else:
            result[record["id"]] = "E"
    return result


def add_locality(crossings, positions):
    all_bounds = [gr.rect_bounds(positions[ref]) for ref in positions]
    total = (min(b[0] for b in all_bounds), min(b[1] for b in all_bounds),
             max(b[2] for b in all_bounds), max(b[3] for b in all_bounds))
    total_diagonal = math.hypot(total[2] - total[0], total[3] - total[1])
    for record in crossings:
        refs = record["relevant_components"]
        bounds = [gr.rect_bounds(positions[ref]) for ref in refs]
        local = (min(b[0] for b in bounds), min(b[1] for b in bounds),
                 max(b[2] for b in bounds), max(b[3] for b in bounds))
        diagonal = math.hypot(local[2] - local[0], local[3] - local[1])
        point = (record["intersection"]["x"], record["intersection"]["y"])
        nearest = min(math.dist(point, positions[ref]) for ref in refs)
        record["locality"] = {"nearest_involved_component_center_distance": gr.rounded(nearest),
                              "relevant_component_body_bounding_box": {
                                  "min_x": gr.rounded(local[0]), "min_y": gr.rounded(local[1]),
                                  "max_x": gr.rounded(local[2]), "max_y": gr.rounded(local[3])},
                              "local_box_diagonal": gr.rounded(diagonal),
                              "total_layout_box_diagonal": gr.rounded(total_diagonal),
                              "fraction_of_total_layout_diagonal": gr.rounded(diagonal / total_diagonal)}


def component_summary(crossings):
    ids = {}
    for record in crossings:
        for ref in record["relevant_components"]:
            ids.setdefault(ref, []).append(record["id"])
    rows = [{"component": ref, "residual_crossings_influenced": len(values), "crossing_ids": values}
            for ref, values in sorted(ids.items())]
    remaining = {record["id"]: set(record["relevant_components"]) for record in crossings}
    groups = []
    while remaining:
        seed = min(remaining)
        group, members = {seed}, set(remaining.pop(seed))
        changed = True
        while changed:
            changed = False
            for identifier, refs in list(remaining.items()):
                if members & refs:
                    group.add(identifier); members |= refs; remaining.pop(identifier); changed = True
        groups.append({"crossing_ids": sorted(group), "components": sorted(members)})
    return rows, groups


def public_result(result):
    if result is None:
        return None
    return {key: value for key, value in result.items() if not key.startswith("_")}


def diagnose(circuit, positions):
    frozen = dict(positions)
    baseline_crossings = crossing_map(circuit, frozen)
    crossings = enumerate_crossings(circuit, frozen)
    if len(crossings) != 7:
        raise ValueError(f"diagnosis requires seven crossings, found {len(crossings)}")
    add_locality(crossings, frozen)
    cache, interactions = {}, {}
    classifications = {name: 0 for name in ("X", "Y", "XY", "STRUCTURAL")}
    for record in crossings:
        x_result = diagnostic_axis_search(circuit, frozen, record, "X", baseline_crossings, cache)
        y_result = diagnostic_axis_search(circuit, frozen, record, "Y", baseline_crossings, cache)
        classification, reason, chosen = classify(x_result, y_result)
        classifications[classification] += 1
        record["classification"] = classification
        record["classification_reason"] = reason
        record["best_x_candidate"] = public_result(x_result)
        record["best_y_candidate"] = public_result(y_result)
        record["selected_diagnostic_axis"] = None if chosen is None else chosen["move"]["axis"]
        interactions[record["id"]] = interaction_codes(crossings, chosen)
    if frozen != positions:
        raise AssertionError("diagnosis mutated frozen baseline coordinates")
    rows, groups = component_summary(crossings)
    nets = sorted({record[name] for record in crossings for name in ("net_a", "net_b")})
    removable = classifications["X"] + classifications["Y"]
    interacting = sum(any(code in ("E", "R") and column != row
                          for column, code in interactions[row].items()) for row in interactions)
    summary = {"classification_counts": classifications,
               "removable_without_increasing_total_crossings": removable,
               "crossings_whose_selected_move_interacts_with_another": interacting,
               "distinct_nets": len(nets), "net_names": nets,
               "distinct_components": len(rows), "component_involvement": rows,
               "component_membership_groups": groups,
               "region_concentration": ("one component-membership-connected group" if len(groups) == 1
                                        else f"{len(groups)} disconnected component-membership groups")}
    return crossings, interactions, summary


def write_annotated_svg(path, source_path, crossings):
    source = source_path.read_text(encoding="utf-8")
    marks = []
    for record in crossings:
        point = record["intersection"]
        marks.append(f'<circle cx="{point["x"]}" cy="{point["y"]}" r="6" fill="none" stroke="#d00" stroke-width="2"/>')
        marks.append(f'<text x="{point["x"] + 8}" y="{point["y"] - 8}" font-size="12" fill="#d00">{record["id"]}</text>')
    path.write_text(source.replace("</svg>", "\n" + "\n".join(marks) + "\n</svg>"), encoding="utf-8")


def text_report(actual, crossings, interactions, summary):
    lines = ["Experiment 4G residual crossing diagnosis", "",
             f"4F final crossings: {actual['net_crossings']}",
             f"4F final overlaps: {actual['component_overlaps']}", ""]
    for record in crossings:
        x, y = record["intersection"]["x"], record["intersection"]["y"]
        lines.append(f"{record['id']}  {record['net_a']} x {record['net_b']}  ({x}, {y})  {record['classification']}")
        lines.append(f"  reason: {record['classification_reason']}")
        lines.append(f"  components: {', '.join(record['relevant_components'])}")
        for axis in ("x", "y"):
            candidate = record[f"best_{axis}_candidate"]
            lines.append(f"  best {axis.upper()}: none" if candidate is None else
                         f"  best {axis.upper()}: {candidate['move']} -> {candidate['total_crossings_after']} crossings, "
                         f"{len(candidate['new_crossings'])} new")
    lines.extend(["", "Classification counts: " + json.dumps(summary["classification_counts"], sort_keys=True),
                  f"Removable without increasing crossings: {summary['removable_without_increasing_total_crossings']}",
                  f"Selected moves interacting with another crossing: {summary['crossings_whose_selected_move_interacts_with_another']}",
                  f"Distinct nets/components: {summary['distinct_nets']}/{summary['distinct_components']}",
                  f"Membership concentration: {summary['region_concentration']}", "",
                  "Interaction matrix (rows target, columns original crossing):",
                  "ID " + " ".join(record["id"] for record in crossings)])
    lines.extend(row + " " + " ".join(interactions[row][record["id"]] for record in crossings)
                 for row in interactions)
    return "\n".join(lines) + "\n"


def run(input_path, output_dir):
    circuit = gr.parse_circuit(input_path)
    positions = load_frozen_state(circuit, output_dir)
    actual, expected = verify_frozen_state(circuit, positions, output_dir)
    crossings, interactions, summary = diagnose(circuit, positions)
    crossing_path = output_dir / "residual_crossings.json"
    interaction_path = output_dir / "residual_crossing_interactions.json"
    gr.write_json(crossing_path, {"source": "closeness_hv_annealed_coordinates.json",
                                  "source_metrics_reproduced": gr.rounded_metrics(actual),
                                  "source_metrics_reported": expected,
                                  "crossing_predicate": "graph_relax.proper_intersection (strict proper 2-D crossing, epsilon 1e-9)",
                                  "segment_model": "terminal-to-derived-arithmetic-centroid for every non-singleton net",
                                  "crossings": crossings, "summary": summary})
    gr.write_json(interaction_path, {"codes": {"E": "eliminated", "P": "persists at same intersection",
                                                 "R": "repositioned or replaced between same nets",
                                                 "N/A": "no target-eliminating move"},
                                         "column_order": [record["id"] for record in crossings],
                                         "matrix": interactions})
    (output_dir / "residual_crossing_report.txt").write_text(
        text_report(actual, crossings, interactions, summary), encoding="utf-8")
    write_annotated_svg(output_dir / "residual_crossings.svg",
                        output_dir / "closeness_hv_annealed.svg", crossings)
    return {"crossings": crossings, "interactions": interactions, "summary": summary,
            "hashes": {"residual_crossings.json": hashlib.sha256(crossing_path.read_bytes()).hexdigest(),
                       "residual_crossing_interactions.json": hashlib.sha256(interaction_path.read_bytes()).hexdigest()}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.input, args.output_dir), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
