#!/usr/bin/env python3
"""Experiment 4F: alternating horizontal and vertical untangling."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import experiment_4e as e4e
import graph_relax as gr


INITIAL_Y_STEP = gr.HEIGHT + gr.BODY_CLEARANCE
MINIMUM_Y_STEP = INITIAL_Y_STEP / 256.0
MAXIMUM_AXIS_PASSES = 20


def horizontal_spans_overlap(first, second):
    return abs(first[0] - second[0]) < gr.WIDTH


def y_slide_candidate(positions, ref, direction, step):
    candidate = dict(positions)
    x, y = candidate[ref]
    candidate[ref] = (x, y + direction * step)
    return candidate


def vertical_swap_candidate(positions, first, second):
    candidate = dict(positions)
    x1, y1 = candidate[first]
    x2, y2 = candidate[second]
    candidate[first], candidate[second] = (x1, y2), (x2, y1)
    return candidate


def best_vertical_move(circuit, positions, step):
    current_metrics = gr.evaluate(circuit, positions)
    current_key = e4e.preprocessing_key(current_metrics)
    best = None
    for ref in sorted(circuit.refs):
        for direction in (-1, 1):
            candidate_positions = y_slide_candidate(positions, ref, direction, step)
            if e4e.has_overlap(circuit, candidate_positions):
                continue
            metrics = gr.evaluate(circuit, candidate_positions)
            key = e4e.preprocessing_key(metrics)
            if not e4e.improves(key, current_key):
                continue
            tie = (0, ref, "", 0 if direction < 0 else 1)
            record = (key + tie, candidate_positions, metrics,
                      {"move_type": "slide", "component": ref, "direction": direction})
            if best is None or record[0] < best[0]:
                best = record
    refs = sorted(circuit.refs)
    for index, first in enumerate(refs):
        for second in refs[index + 1:]:
            if not horizontal_spans_overlap(positions[first], positions[second]):
                continue
            candidate_positions = vertical_swap_candidate(positions, first, second)
            if e4e.has_overlap(circuit, candidate_positions):
                continue
            metrics = gr.evaluate(circuit, candidate_positions)
            key = e4e.preprocessing_key(metrics)
            if not e4e.improves(key, current_key):
                continue
            tie = (1, first, second, 0)
            record = (key + tie, candidate_positions, metrics,
                      {"move_type": "swap", "components": (first, second)})
            if best is None or record[0] < best[0]:
                best = record
    return best, current_metrics


def vertical_untangle(circuit, start, *, initial_step=INITIAL_Y_STEP,
                      minimum_step=MINIMUM_Y_STEP):
    positions = dict(start)
    original_x = {ref: point[0] for ref, point in positions.items()}
    step = initial_step
    halvings = slides = swaps = sweeps = 0
    moves = []
    while True:
        sweeps += 1
        best, before = best_vertical_move(circuit, positions, step)
        if best is None:
            if step > minimum_step:
                step = max(minimum_step, step / 2.0)
                halvings += 1
                continue
            reason = "no lexicographically improving legal move at minimum Y step"
            break
        _, new_positions, after, description = best
        if description["move_type"] == "slide":
            ref = description["component"]
            coordinate_change = {"component": ref, "old_y": positions[ref][1],
                                 "new_y": new_positions[ref][1]}
            slides += 1
        else:
            first, second = description["components"]
            coordinate_change = {"components": [first, second],
                                 "old_y": [positions[first][1], positions[second][1]],
                                 "new_y": [new_positions[first][1], new_positions[second][1]]}
            swaps += 1
        moves.append({"move_number": len(moves) + 1, "move_type": description["move_type"],
                      **coordinate_change, "crossings_before": before["net_crossings"],
                      "crossings_after": after["net_crossings"],
                      "connection_length_before": gr.rounded(before["total_connection_length"]),
                      "connection_length_after": gr.rounded(after["total_connection_length"])})
        positions = new_positions
    if any(positions[ref][0] != original_x[ref] for ref in circuit.refs):
        raise AssertionError("vertical preprocessing changed an X coordinate")
    return positions, moves, {"accepted_slides": slides, "accepted_swaps": swaps,
                              "initial_y_step": initial_step, "minimum_y_step": minimum_step,
                              "step_halvings": halvings, "search_sweeps": sweeps,
                              "stopping_reason": reason}


def run_axis_pass(circuit, positions, axis):
    if axis == "H":
        return e4e.horizontal_untangle(circuit, positions)
    if axis == "V":
        return vertical_untangle(circuit, positions)
    raise ValueError(f"unknown axis {axis}")


def coordinate_snapshot(positions):
    return {ref: {"x": gr.rounded(point[0]), "y": gr.rounded(point[1])}
            for ref, point in positions.items()}


def alternating_untangle(circuit, start, *, maximum_passes=MAXIMUM_AXIS_PASSES):
    positions = dict(start)
    passes, moves = [], []
    latest_moves = {"H": None, "V": None}
    for pass_number in range(1, maximum_passes + 1):
        axis = "H" if pass_number % 2 else "V"
        before = gr.evaluate(circuit, positions)
        new_positions, local_moves, stats = run_axis_pass(circuit, positions, axis)
        after = gr.evaluate(circuit, new_positions)
        for local in local_moves:
            record = dict(local)
            record["global_move_number"] = len(moves) + 1
            record["pass_number"] = pass_number
            record["axis"] = axis
            record["reduced_crossings"] = record["crossings_after"] < record["crossings_before"]
            record.pop("move_number", None)
            moves.append(record)
        accepted = len(local_moves)
        latest_moves[axis] = accepted
        passes.append({"pass_number": pass_number, "axis": axis,
                       "accepted_move_count": accepted,
                       "accepted_slides": stats["accepted_slides"],
                       "accepted_swaps": stats["accepted_swaps"],
                       "crossings_before": before["net_crossings"],
                       "crossings_after": after["net_crossings"],
                       "connection_length_before": gr.rounded(before["total_connection_length"]),
                       "connection_length_after": gr.rounded(after["total_connection_length"]),
                       "ordinary_metrics_before": gr.rounded_metrics(before),
                       "ordinary_metrics_after": gr.rounded_metrics(after),
                       "coordinates": coordinate_snapshot(new_positions)})
        positions = new_positions
        if latest_moves["H"] == 0 and latest_moves["V"] == 0:
            reason = "latest horizontal and vertical passes both accepted zero moves"
            converged = True
            break
    else:
        reason = "20-pass safety limit reached while moves were still occurring"
        converged = False
    return positions, passes, moves, {"stopping_reason": reason, "converged": converged,
                                     "maximum_axis_passes": maximum_passes}


def release_preprocessing(positions):
    return dict(positions)


def change(final, baseline):
    absolute = final - baseline
    return {"absolute": gr.rounded(absolute) if isinstance(absolute, float) else absolute,
            "percent": None if baseline == 0 else gr.rounded(100.0 * absolute / baseline)}


def comparison_changes(final_metrics, baseline):
    return {
        "total_connection_length": change(final_metrics["total_connection_length"], baseline["total_connection_length"]),
        "component_overlaps": change(final_metrics["component_overlaps"], baseline["component_overlaps"]),
        "clearance_violations": change(final_metrics["clearance_violations"], baseline["clearance_violations"]),
        "net_crossings": change(final_metrics["net_crossings"], baseline["net_crossings"]),
        "crowding": change(final_metrics["energy"]["crowding"], baseline["energy"]["crowding"]),
        "objective": change(final_metrics["objective"], baseline["objective"])}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    circuit = gr.parse_circuit(args.input)
    start = e4e.load_closeness_initial(circuit, args.output_dir)
    source_actual, source_report = e4e.verify_source(circuit, start, args.output_dir)
    preprocessed, passes, moves, stopping = alternating_untangle(circuit, start)
    preprocessed_metrics = gr.evaluate(circuit, preprocessed)
    gr.write_json(args.output_dir / "closeness_hv_passes.json", {
        "start_metrics": gr.rounded_metrics(source_actual), "passes": passes, **stopping})
    gr.write_json(args.output_dir / "closeness_hv_moves.json", {
        "preprocessing_objective": ["net_crossings", "total_connection_length"],
        "crossing_reducing_move_numbers": [m["global_move_number"] for m in moves if m["reduced_crossings"]],
        "moves": moves})
    gr.write_svg(args.output_dir / "closeness_hv.svg", circuit, preprocessed)
    preprocessed_coordinates = args.output_dir / "closeness_hv_coordinates.json"
    gr.write_json(preprocessed_coordinates, gr.coordinate_document(circuit, preprocessed))
    h_passes = [p for p in passes if p["axis"] == "H"]
    v_passes = [p for p in passes if p["axis"] == "V"]
    summary = {"horizontal_passes": len(h_passes), "vertical_passes": len(v_passes),
               "horizontal_slides": sum(p["accepted_slides"] for p in h_passes),
               "horizontal_swaps": sum(p["accepted_swaps"] for p in h_passes),
               "vertical_slides": sum(p["accepted_slides"] for p in v_passes),
               "vertical_swaps": sum(p["accepted_swaps"] for p in v_passes),
               "total_preprocessing_moves": len(moves),
               "crossing_reducing_moves": sum(m["reduced_crossings"] for m in moves), **stopping}
    horizontal_metrics_report = json.loads(
        (args.output_dir / "closeness_horizontal_annealed_metrics.json").read_text())
    four_e_horizontal = horizontal_metrics_report["states"]["experiment_4e_horizontal_preprocessed"]
    metrics_report = {"source_verification": {"stored": gr.rounded_metrics(source_actual),
                                               "frozen_4c": source_report},
                      "states": {"experiment_4c_initial": source_report,
                                 "experiment_4e_horizontal_preprocessed": four_e_horizontal,
                                 "experiment_4f_hv_preprocessed": gr.rounded_metrics(preprocessed_metrics)},
                      "summary": summary,
                      "coordinate_sha256": hashlib.sha256(preprocessed_coordinates.read_bytes()).hexdigest()}
    gr.write_json(args.output_dir / "closeness_hv_metrics.json", metrics_report)
    best, anneal_stats = gr.anneal(circuit, release_preprocessing(preprocessed))
    final, quench_sweeps, quench_reason = gr.relax(circuit, best)
    final_metrics = gr.evaluate(circuit, final)
    gr.write_svg(args.output_dir / "closeness_hv_annealed.svg", circuit, final)
    final_coordinates = args.output_dir / "closeness_hv_annealed_coordinates.json"
    gr.write_json(final_coordinates, gr.coordinate_document(circuit, final))
    closeness_report = json.loads((args.output_dir / "closeness_metrics.json").read_text())
    four_c_final = closeness_report["states"]["experiment_4c_closeness_annealed_and_quenched"]
    four_e_final = horizontal_metrics_report["states"]["experiment_4e_annealed_and_quenched"]
    final_report = {
        "states": {"experiment_4c_initial": source_report,
                   "experiment_4e_horizontal_preprocessed": four_e_horizontal,
                   "experiment_4f_hv_preprocessed": gr.rounded_metrics(preprocessed_metrics),
                   "experiment_4f_annealed_and_quenched": gr.rounded_metrics(final_metrics)},
        "changes_from_experiment_4c_final": comparison_changes(final_metrics, four_c_final),
        "changes_from_experiment_4e_final": comparison_changes(final_metrics, four_e_final),
        "preprocessing_summary": summary,
        "crossing_sequence": [{"state": "4C start", "crossings": source_actual["net_crossings"]}] +
                             [{"state": f"{p['axis']}{(p['pass_number'] + 1) // 2}",
                               "crossings": p["crossings_after"]} for p in passes],
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
                   "sweeps": quench_sweeps, "stopping_reason": quench_reason},
        "preprocessed_coordinate_sha256": hashlib.sha256(preprocessed_coordinates.read_bytes()).hexdigest(),
        "final_coordinate_sha256": hashlib.sha256(final_coordinates.read_bytes()).hexdigest()}
    gr.write_json(args.output_dir / "closeness_hv_annealed_metrics.json", final_report)
    print(json.dumps({"preprocessing": metrics_report, "final": final_report},
                     sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
