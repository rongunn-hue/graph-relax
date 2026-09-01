#!/usr/bin/env python3
"""Experiment 4H: junction-only untangling from frozen 4F components."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import experiment_4g as e4g
import graph_relax as gr


INITIAL_JUNCTION_STEP = max(gr.WIDTH, gr.HEIGHT) + gr.BODY_CLEARANCE
MINIMUM_JUNCTION_STEP = INITIAL_JUNCTION_STEP / 256.0
BOUND_MARGIN = INITIAL_JUNCTION_STEP
DIRECTIONS = ((-1.0, 0.0, "-X"), (1.0, 0.0, "+X"),
              (0.0, -1.0, "-Y"), (0.0, 1.0, "+Y"))


def movable_nets(circuit):
    return tuple(sorted(net.name for net in circuit.nets if len(net.terminals) >= 3))


def junction_bounds(circuit, positions, margin=BOUND_MARGIN):
    terminals = gr.terminal_positions(circuit, positions)
    result = {}
    for net in circuit.nets:
        if len(net.terminals) < 3:
            continue
        points = [terminals[name] for name in net.terminals]
        result[net.name] = (min(p[0] for p in points) - margin,
                            min(p[1] for p in points) - margin,
                            max(p[0] for p in points) + margin,
                            max(p[1] for p in points) + margin)
    return result


def in_bounds(point, bounds):
    return bounds[0] <= point[0] <= bounds[2] and bounds[1] <= point[1] <= bounds[3]


def junction_outside_unrelated_bodies(circuit, positions, net_name, point):
    net = next(net for net in circuit.nets if net.name == net_name)
    related = {terminal.rsplit(".", 1)[0] for terminal in net.terminals}
    for ref in circuit.refs:
        if ref in related:
            continue
        left, top, right, bottom = gr.rect_bounds(positions[ref])
        if left < point[0] < right and top < point[1] < bottom:
            return False
    return True


def junction_crossing_map(circuit, positions, junctions):
    terminals = gr.terminal_positions(circuit, positions)
    segments = gr.junction_segments(circuit, terminals, junctions)
    result = {}
    for index, first in enumerate(segments):
        for second in segments[index + 1:]:
            if first[0] == second[0] or not gr.proper_intersection(first[2], first[3], second[2], second[3]):
                continue
            identities = tuple(sorted(((first[0], first[1]), (second[0], second[1]))))
            result[identities] = {
                "net_pair": tuple(sorted((first[0], second[0]))),
                "intersection": e4g.intersection_point(first[2], first[3], second[2], second[3]),
                "segments": (first, second)}
    return result


def preprocessing_key(metrics):
    return metrics["net_crossings"], metrics["total_connection_length"]


def improves(candidate, current):
    return candidate[0] < current[0] or (candidate[0] == current[0] and candidate[1] < current[1] - 1e-12)


def move_candidate(junctions, net_name, dx, dy, step):
    candidate = dict(junctions)
    x, y = candidate[net_name]
    candidate[net_name] = (x + dx * step, y + dy * step)
    return candidate


def best_junction_move(circuit, positions, junctions, bounds, step):
    current_metrics = gr.evaluate_junction_state(circuit, positions, junctions)
    current_key = preprocessing_key(current_metrics)
    best = None
    for net_name in movable_nets(circuit):
        for direction_index, (dx, dy, label) in enumerate(DIRECTIONS):
            candidate = move_candidate(junctions, net_name, dx, dy, step)
            point = candidate[net_name]
            if not in_bounds(point, bounds[net_name]):
                continue
            if not junction_outside_unrelated_bodies(circuit, positions, net_name, point):
                continue
            metrics = gr.evaluate_junction_state(circuit, positions, candidate)
            key = preprocessing_key(metrics)
            if not improves(key, current_key):
                continue
            record = (key + (net_name, direction_index), candidate, metrics,
                      {"net": net_name, "direction": label})
            if best is None or record[0] < best[0]:
                best = record
    return best, current_metrics


def junction_untangle(circuit, positions, initial, bounds, *,
                      initial_step=INITIAL_JUNCTION_STEP,
                      minimum_step=MINIMUM_JUNCTION_STEP):
    frozen_components = dict(positions)
    junctions = dict(initial)
    moves = []
    step = initial_step
    halvings = sweeps = crossing_moves = 0
    while True:
        sweeps += 1
        best, before = best_junction_move(circuit, frozen_components, junctions, bounds, step)
        if best is None:
            if step > minimum_step:
                step = max(minimum_step, step / 2.0)
                halvings += 1
                continue
            reason = "no lexicographically improving legal junction move at minimum step"
            break
        _, candidate, after, description = best
        net_name = description["net"]
        reduced = after["net_crossings"] < before["net_crossings"]
        crossing_moves += int(reduced)
        moves.append({"move_number": len(moves) + 1, **description,
                      "step": gr.rounded(step),
                      "old_xy": [gr.rounded(value) for value in junctions[net_name]],
                      "new_xy": [gr.rounded(value) for value in candidate[net_name]],
                      "crossings_before": before["net_crossings"],
                      "crossings_after": after["net_crossings"],
                      "connection_length_before": gr.rounded(before["total_connection_length"]),
                      "connection_length_after": gr.rounded(after["total_connection_length"]),
                      "reduced_crossings": reduced})
        junctions = candidate
    if frozen_components != positions:
        raise AssertionError("junction search mutated frozen component coordinates")
    return junctions, moves, {"accepted_junction_moves": len(moves),
                              "crossing_reducing_junction_moves": crossing_moves,
                              "initial_step": initial_step, "minimum_step": minimum_step,
                              "step_halvings": halvings, "search_sweeps": sweeps,
                              "stopping_reason": reason}


def point_document(point):
    return {"x": gr.rounded(point[0]), "y": gr.rounded(point[1])}


def bounds_document(bounds):
    return {name: {"min_x": gr.rounded(value[0]), "min_y": gr.rounded(value[1]),
                   "max_x": gr.rounded(value[2]), "max_y": gr.rounded(value[3])}
            for name, value in sorted(bounds.items())}


def coordinate_document(circuit, positions, initial, final, bounds):
    return {"circuit": circuit.name,
            "components_frozen_from": "closeness_hv_annealed_coordinates.json",
            "components": {ref: point_document(positions[ref]) for ref in circuit.refs},
            "initial_junctions": {name: point_document(initial[name]) for name in sorted(initial)},
            "final_junctions": {name: point_document(final[name]) for name in sorted(final)},
            "junction_bounds": bounds_document(bounds)}


def compare_original_crossings(original_records, final_crossings):
    final_by_pair = {}
    for item in final_crossings.values():
        final_by_pair.setdefault(item["net_pair"], []).append(item["intersection"])
    rows = []
    original_pairs = set()
    for original in original_records:
        pair = tuple(sorted((original["net_a"], original["net_b"])))
        original_pairs.add(pair)
        candidates = final_by_pair.get(pair, [])
        old = (original["intersection"]["x"], original["intersection"]["y"])
        if not candidates:
            status = "eliminated"
            final_points = []
        elif any(math.dist(old, point) <= e4g.POINT_TOLERANCE for point in candidates):
            status = "persists"
            final_points = candidates
        else:
            status = "displaced/replaced"
            final_points = candidates
        rows.append({"id": original["id"], "net_a": original["net_a"], "net_b": original["net_b"],
                     "status": status,
                     "original_intersection": original["intersection"],
                     "final_same_net_pair_intersections": [point_document(p) for p in final_points]})
    new_pairs = sorted(set(final_by_pair) - original_pairs)
    return rows, [{"net_a": pair[0], "net_b": pair[1],
                   "intersections": [point_document(p) for p in final_by_pair[pair]]}
                  for pair in new_pairs]


def write_diagnostic_svg(path, circuit, positions, junctions, crossing_status):
    gr.write_junction_svg(path, circuit, positions, junctions)
    source = path.read_text(encoding="utf-8")
    marks = []
    final_crossings = junction_crossing_map(circuit, positions, junctions)
    for index, item in enumerate(sorted(final_crossings.values(),
                                        key=lambda value: (value["intersection"][1], value["intersection"][0],
                                                           value["net_pair"])), 1):
        x, y = item["intersection"]
        label = f"H{index:03d}"
        marks.append(f'<circle cx="{x:.6f}" cy="{y:.6f}" r="6" fill="none" stroke="#d00" stroke-width="2"/>')
        marks.append(f'<text x="{x+8:.6f}" y="{y-8:.6f}" font-size="12" fill="#d00">{label}</text>')
    path.write_text(source.replace("</svg>", "\n" + "\n".join(marks) + "\n</svg>"), encoding="utf-8")


def run(input_path, output_dir):
    circuit = gr.parse_circuit(input_path)
    positions = e4g.load_frozen_state(circuit, output_dir)
    frozen_metrics, frozen_report = e4g.verify_frozen_state(circuit, positions, output_dir)
    initial = gr.initial_junctions(circuit, positions)
    bounds = junction_bounds(circuit, positions)
    if set(initial) != set(movable_nets(circuit)):
        raise AssertionError("movable junction set differs from 3+-terminal nets")
    initial_metrics = gr.evaluate_junction_state(circuit, positions, initial)
    if initial_metrics["net_crossings"] != 7:
        raise ValueError(f"initial junction model has {initial_metrics['net_crossings']} crossings, expected 7")
    final, moves, search = junction_untangle(circuit, positions, initial, bounds)
    final_metrics = gr.evaluate_junction_state(circuit, positions, final)
    original = json.loads((output_dir / "residual_crossings.json").read_text())["crossings"]
    status, new_pairs = compare_original_crossings(original, junction_crossing_map(circuit, positions, final))
    coordinate_path = output_dir / "junction_untangled_coordinates.json"
    move_path = output_dir / "junction_untangled_moves.json"
    gr.write_json(coordinate_path, coordinate_document(circuit, positions, initial, final, bounds))
    gr.write_json(move_path, {"direction_order": [item[2] for item in DIRECTIONS],
                              "preprocessing_objective": ["net_crossings", "total_connection_length"],
                              "moves": moves})
    report = {"frozen_4f_baseline": gr.rounded_metrics(frozen_metrics),
              "frozen_4f_reported": frozen_report,
              "junction_model_initial": gr.rounded_metrics(initial_metrics),
              "junction_only_final": gr.rounded_metrics(final_metrics),
              "movable_junction_count": len(initial),
              "movable_nets": [{"net": name, "initial": point_document(initial[name]),
                                  "bounds": bounds_document({name: bounds[name]})[name]}
                                 for name in sorted(initial)],
              "search": search, "original_4g_crossing_status": status,
              "new_net_pair_crossings": new_pairs,
              "component_coordinates_unchanged": True,
              "coordinate_sha256": hashlib.sha256(coordinate_path.read_bytes()).hexdigest(),
              "move_log_sha256": hashlib.sha256(move_path.read_bytes()).hexdigest()}
    gr.write_json(output_dir / "junction_untangled_metrics.json", report)
    write_diagnostic_svg(output_dir / "junction_untangled.svg", circuit, positions, final, status)
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.input, args.output_dir), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
