#!/usr/bin/env python3
"""Experiment 4E: deterministic horizontal untangling of stored 4C state."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import graph_relax as gr


INITIAL_X_STEP = gr.WIDTH + gr.BODY_CLEARANCE
MINIMUM_X_STEP = INITIAL_X_STEP / 256.0
SOURCE_NUMERIC_TOLERANCE = 2e-5


def load_closeness_initial(circuit, output_dir):
    graph = json.loads((output_dir / "closeness_graph.json").read_text(encoding="utf-8"))
    stored = graph["initialization_coordinates"]
    if set(stored) != set(circuit.refs):
        raise ValueError("stored 4C initialization component set does not match input")
    return {ref: (float(stored[ref]["x"]), float(stored[ref]["y"]))
            for ref in circuit.refs}


def verify_source(circuit, positions, output_dir):
    actual = gr.evaluate(circuit, positions)
    report = json.loads((output_dir / "closeness_metrics.json").read_text(encoding="utf-8"))
    expected = report["closeness_initial"]
    if actual["net_crossings"] != expected["net_crossings"]:
        raise ValueError(f"stored 4C crossings differ: {actual['net_crossings']} != {expected['net_crossings']}")
    if actual["component_overlaps"] != expected["component_overlaps"]:
        raise ValueError(f"stored 4C overlaps differ: {actual['component_overlaps']} != {expected['component_overlaps']}")
    if abs(actual["total_connection_length"] - expected["total_connection_length"]) > SOURCE_NUMERIC_TOLERANCE:
        raise ValueError("stored 4C connection length exceeds serialization tolerance")
    if abs(actual["objective"] - expected["objective"]) > SOURCE_NUMERIC_TOLERANCE:
        raise ValueError("stored 4C objective exceeds serialization tolerance")
    return actual, expected


def has_overlap(circuit, positions):
    return any(gr.overlap_area(positions[first], positions[second]) > 0.0
               for index, first in enumerate(circuit.refs)
               for second in circuit.refs[index + 1:])


def vertical_spans_overlap(first, second):
    return abs(first[1] - second[1]) < gr.HEIGHT


def slide_candidate(positions, ref, direction, step):
    candidate = dict(positions)
    x, y = candidate[ref]
    candidate[ref] = (x + direction * step, y)
    return candidate


def swap_candidate(positions, first, second):
    candidate = dict(positions)
    x1, y1 = candidate[first]
    x2, y2 = candidate[second]
    candidate[first], candidate[second] = (x2, y1), (x1, y2)
    return candidate


def preprocessing_key(metrics):
    return metrics["net_crossings"], metrics["total_connection_length"]


def improves(candidate_key, current_key):
    return candidate_key[0] < current_key[0] or \
           (candidate_key[0] == current_key[0] and candidate_key[1] < current_key[1] - 1e-12)


def best_horizontal_move(circuit, positions, step):
    current_metrics = gr.evaluate(circuit, positions)
    current_key = preprocessing_key(current_metrics)
    best = None
    for ref in sorted(circuit.refs):
        for direction in (-1, 1):
            candidate_positions = slide_candidate(positions, ref, direction, step)
            if has_overlap(circuit, candidate_positions):
                continue
            metrics = gr.evaluate(circuit, candidate_positions)
            key = preprocessing_key(metrics)
            if not improves(key, current_key):
                continue
            tie = (0, ref, "", 0 if direction < 0 else 1)
            record = (key + tie, candidate_positions, metrics,
                      {"move_type": "slide", "component": ref, "direction": direction})
            if best is None or record[0] < best[0]:
                best = record
    refs = sorted(circuit.refs)
    for index, first in enumerate(refs):
        for second in refs[index + 1:]:
            if not vertical_spans_overlap(positions[first], positions[second]):
                continue
            candidate_positions = swap_candidate(positions, first, second)
            if has_overlap(circuit, candidate_positions):
                continue
            metrics = gr.evaluate(circuit, candidate_positions)
            key = preprocessing_key(metrics)
            if not improves(key, current_key):
                continue
            tie = (1, first, second, 0)
            record = (key + tie, candidate_positions, metrics,
                      {"move_type": "swap", "components": (first, second)})
            if best is None or record[0] < best[0]:
                best = record
    return best, current_metrics


def horizontal_untangle(circuit, start, *, initial_step=INITIAL_X_STEP,
                        minimum_step=MINIMUM_X_STEP):
    positions = dict(start)
    original_y = {ref: point[1] for ref, point in positions.items()}
    step = initial_step
    halvings = slides = swaps = sweeps = 0
    moves = []
    while True:
        sweeps += 1
        best, before = best_horizontal_move(circuit, positions, step)
        if best is None:
            if step > minimum_step:
                step = max(minimum_step, step / 2.0)
                halvings += 1
                continue
            reason = "no lexicographically improving legal move at minimum X step"
            break
        _, new_positions, after, description = best
        move_number = len(moves) + 1
        if description["move_type"] == "slide":
            ref = description["component"]
            coordinate_change = {"component": ref, "old_x": positions[ref][0],
                                 "new_x": new_positions[ref][0]}
            slides += 1
        else:
            first, second = description["components"]
            coordinate_change = {"components": [first, second],
                                 "old_x": [positions[first][0], positions[second][0]],
                                 "new_x": [new_positions[first][0], new_positions[second][0]]}
            swaps += 1
        moves.append({"move_number": move_number, "move_type": description["move_type"],
                      **coordinate_change,
                      "crossings_before": before["net_crossings"],
                      "crossings_after": after["net_crossings"],
                      "connection_length_before": gr.rounded(before["total_connection_length"]),
                      "connection_length_after": gr.rounded(after["total_connection_length"])})
        positions = new_positions
    if any(positions[ref][1] != original_y[ref] for ref in circuit.refs):
        raise AssertionError("horizontal preprocessing changed a Y coordinate")
    return positions, moves, {"accepted_slides": slides, "accepted_swaps": swaps,
                              "initial_x_step": initial_step, "minimum_x_step": minimum_step,
                              "step_halvings": halvings, "search_sweeps": sweeps,
                              "stopping_reason": reason}


def release_preprocessing(positions):
    return dict(positions)


def change(final, baseline):
    absolute = final - baseline
    return {"absolute": gr.rounded(absolute) if isinstance(absolute, float) else absolute,
            "percent": None if baseline == 0 else gr.rounded(100.0 * absolute / baseline)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    circuit = gr.parse_circuit(args.input)
    start = load_closeness_initial(circuit, args.output_dir)
    source_actual, source_report = verify_source(circuit, start, args.output_dir)
    horizontal, moves, search = horizontal_untangle(circuit, start)
    horizontal_metrics = gr.evaluate(circuit, horizontal)
    gr.write_svg(args.output_dir / "closeness_horizontal.svg", circuit, horizontal)
    horizontal_coordinates = args.output_dir / "closeness_horizontal_coordinates.json"
    gr.write_json(horizontal_coordinates, gr.coordinate_document(circuit, horizontal))
    gr.write_json(args.output_dir / "closeness_horizontal_moves.json", {
        "preprocessing_objective": ["net_crossings", "total_connection_length"],
        "moves": moves})
    horizontal_report = {
        "source_verification": {"serialization_tolerance": SOURCE_NUMERIC_TOLERANCE,
                                "stored_state_metrics": gr.rounded_metrics(source_actual),
                                "frozen_4c_initial_metrics": source_report},
        "search": search, "before": gr.rounded_metrics(source_actual),
        "after": gr.rounded_metrics(horizontal_metrics),
        "coordinate_sha256": hashlib.sha256(horizontal_coordinates.read_bytes()).hexdigest()}
    gr.write_json(args.output_dir / "closeness_horizontal_metrics.json", horizontal_report)
    best, anneal_stats = gr.anneal(circuit, release_preprocessing(horizontal))
    final, quench_sweeps, reason = gr.relax(circuit, best)
    final_metrics = gr.evaluate(circuit, final)
    gr.write_svg(args.output_dir / "closeness_horizontal_annealed.svg", circuit, final)
    final_coordinates = args.output_dir / "closeness_horizontal_annealed_coordinates.json"
    gr.write_json(final_coordinates, gr.coordinate_document(circuit, final))
    closeness_report = json.loads((args.output_dir / "closeness_metrics.json").read_text())
    four_c_final = closeness_report["states"]["experiment_4c_closeness_annealed_and_quenched"]
    changes = {
        "total_connection_length": change(final_metrics["total_connection_length"], four_c_final["total_connection_length"]),
        "component_overlaps": change(final_metrics["component_overlaps"], four_c_final["component_overlaps"]),
        "clearance_violations": change(final_metrics["clearance_violations"], four_c_final["clearance_violations"]),
        "net_crossings": change(final_metrics["net_crossings"], four_c_final["net_crossings"]),
        "crowding": change(final_metrics["energy"]["crowding"], four_c_final["energy"]["crowding"]),
        "objective": change(final_metrics["objective"], four_c_final["objective"]),
    }
    final_report = {
        "states": {"experiment_4c_initial": source_report,
                   "experiment_4e_horizontal_preprocessed": gr.rounded_metrics(horizontal_metrics),
                   "experiment_4e_annealed_and_quenched": gr.rounded_metrics(final_metrics)},
        "changes_from_experiment_4c_final": changes,
        "preprocessing": search,
        "annealing": {"seed": gr.ANNEAL_SEED, "initial_temperature": gr.INITIAL_TEMPERATURE,
                      "final_temperature": gr.FINAL_TEMPERATURE, "cooling_schedule": "geometric",
                      "initial_move_distance": gr.INITIAL_MOVE_DISTANCE,
                      "final_move_distance": gr.FINAL_MOVE_DISTANCE,
                      "move_distance_schedule": "geometric", "iterations": gr.ANNEAL_ITERATIONS,
                      "accepted_downhill_moves": anneal_stats["accepted_downhill_moves"],
                      "accepted_uphill_moves": anneal_stats["accepted_uphill_moves"],
                      "rejected_moves": anneal_stats["rejected_moves"],
                      "best_energy_encountered": gr.rounded(anneal_stats["best_energy"])},
        "quench": {"final_energy": gr.rounded(final_metrics["objective"]),
                   "sweeps": quench_sweeps, "stopping_reason": reason},
        "horizontal_coordinate_sha256": hashlib.sha256(horizontal_coordinates.read_bytes()).hexdigest(),
        "final_coordinate_sha256": hashlib.sha256(final_coordinates.read_bytes()).hexdigest()}
    gr.write_json(args.output_dir / "closeness_horizontal_annealed_metrics.json", final_report)
    print(json.dumps({"horizontal": horizontal_report, "final": final_report},
                     sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
