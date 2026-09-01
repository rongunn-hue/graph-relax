#!/usr/bin/env python3
"""Experiment 4O: coarse local compaction with angular-distribution preference."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

import experiment_4i as e4i
import experiment_4j as e4j
import experiment_4n as e4n
import experiment_4p as e4p
import graph_relax as gr


COMPARABLE_REDUCTION_FRACTION = 0.90


def vertex_angular_measure(graph, positions, vertex):
    neighbors = sorted(graph.neighbors(vertex))
    degree = len(neighbors)
    if degree < 2:
        return None
    center = positions[vertex]
    angles = sorted((math.degrees(math.atan2(positions[neighbor][1]-center[1],
                                             positions[neighbor][0]-center[0])) % 360.0)
                    for neighbor in neighbors)
    gaps = [angles[index+1]-angles[index] for index in range(degree-1)]
    gaps.append(angles[0]+360.0-angles[-1])
    ideal = 360.0/degree
    return {"degree": degree, "minimum_gap": min(gaps),
            "ideal_gap": ideal,
            "angular_error": sum(abs(gap-ideal) for gap in gaps),
            "angles": angles, "gaps": gaps}


def angular_summary(graph, positions):
    rows = {vertex: vertex_angular_measure(graph, positions, vertex)
            for vertex in sorted(graph) if graph.degree(vertex) >= 2}
    return {"eligible_vertex_count": len(rows),
            "global_minimum_incident_angular_gap": gr.rounded(
                min(row["minimum_gap"] for row in rows.values())),
            "mean_per_vertex_minimum_angular_gap": gr.rounded(
                sum(row["minimum_gap"] for row in rows.values())/len(rows)),
            "mean_angular_error": gr.rounded(
                sum(row["angular_error"] for row in rows.values())/len(rows)),
            "vertices": {vertex: {key: gr.rounded(value) if key not in ("degree", "angles", "gaps") else
                                   ([gr.rounded(item) for item in value] if isinstance(value, list) else value)
                                  for key, value in row.items()}
                         for vertex, row in rows.items()}}


def eligible_candidates(circuit, graph, embedding, positions, realization, vertex, step):
    original_length = e4n.local_length(vertex, realization[3])
    candidates, trials = [], 0
    for direction_index, direction in enumerate(e4n.candidate_directions(graph, positions, vertex)):
        if (direction["target_distance"] is not None and
                step > direction["target_distance"]+e4n.EPSILON):
            continue
        point = e4n.candidate_point(positions[vertex], direction["unit"], step)
        if point == positions[vertex]:
            continue
        trial_positions = dict(positions); trial_positions[vertex] = point
        try:
            trial_geometry = e4n.geometry(circuit, trial_positions)
        except ValueError:
            continue
        trials += 1
        if not e4n.locally_valid(circuit, embedding, trial_positions, vertex, trial_geometry):
            continue
        new_length = e4n.local_length(vertex, trial_geometry[3])
        reduction = original_length-new_length
        if reduction <= 1e-9:
            continue
        angular = vertex_angular_measure(graph, trial_positions, vertex)
        candidates.append({"reduction": reduction, "direction_index": direction_index,
                           "direction": direction["label"], "point": point,
                           "positions": trial_positions, "geometry": trial_geometry,
                           "length_before": original_length, "length_after": new_length,
                           "angular": angular})
    return candidates, trials


def select_candidate(graph, vertex, candidates):
    if not candidates:
        return None
    if graph.degree(vertex) == 1:
        return min(candidates, key=lambda row: (-row["reduction"], row["direction_index"]))
    best_reduction = max(row["reduction"] for row in candidates)
    comparable = [row for row in candidates
                  if row["reduction"]+1e-9 >= COMPARABLE_REDUCTION_FRACTION*best_reduction]
    return min(comparable, key=lambda row: (row["angular"]["angular_error"],
                                            -row["angular"]["minimum_gap"],
                                            -row["reduction"], row["direction_index"]))


def positions_from_coordinates(circuit, path):
    document = json.loads(path.read_text())
    positions = {e4p.component_vertex(ref):
                 (float(document["components"][ref]["x"]), float(document["components"][ref]["y"]))
                 for ref in circuit.refs}
    positions.update({e4p.net_vertex(net.name):
                      (float(document["net_junctions"][net.name]["x"]),
                       float(document["net_junctions"][net.name]["y"]))
                      for net in circuit.nets})
    return positions


def percent_change(new, old):
    return gr.rounded(100*(new-old)/old)


def run(input_path, output_dir):
    started = time.perf_counter()
    circuit = gr.parse_circuit(input_path)
    graph, embedding, positions, order = e4n.load_start(circuit, output_dir)
    baseline_validation, realization = e4n.complete_validate(circuit, graph, embedding, positions)
    if not baseline_validation["valid"]:
        raise ValueError("stored Experiment 4K starting geometry failed reproduction")
    initial_positions = dict(positions)
    moves, passes = [], []
    for pass_number, step in enumerate(e4n.STEPS, 1):
        moved = trials = 0
        movement = 0.0
        for vertex in order:
            candidates, vertex_trials = eligible_candidates(
                circuit, graph, embedding, positions, realization, vertex, step)
            trials += vertex_trials
            chosen = select_candidate(graph, vertex, candidates)
            if chosen is None:
                continue
            old_point = positions[vertex]
            positions, realization = chosen["positions"], chosen["geometry"]
            actual_movement = math.dist(old_point, chosen["point"])
            moved += 1; movement += actual_movement
            moves.append({"move_number": len(moves)+1, "pass": pass_number, "step": step,
                          "vertex": vertex, "degree": graph.degree(vertex),
                          "direction": chosen["direction"], "old": e4i.point_doc(old_point),
                          "new": e4i.point_doc(chosen["point"]),
                          "movement_distance": gr.rounded(actual_movement),
                          "local_length_before": gr.rounded(chosen["length_before"]),
                          "local_length_after": gr.rounded(chosen["length_after"]),
                          "local_length_reduction": gr.rounded(chosen["reduction"]),
                          "resulting_minimum_gap": None if chosen["angular"] is None else gr.rounded(chosen["angular"]["minimum_gap"]),
                          "resulting_angular_error": None if chosen["angular"] is None else gr.rounded(chosen["angular"]["angular_error"])})
        validation, checked_geometry = e4n.complete_validate(circuit, graph, embedding, positions)
        if not validation["valid"]:
            raise RuntimeError(f"complete validation failed after step {step}")
        realization = checked_geometry
        measured = e4i.measure(circuit, realization[0], realization[3])
        extents = e4j.geometry_extents(realization[0], realization[1])
        passes.append({"pass": pass_number, "step": step,
                       "vertices_moved": moved, "accepted_moves": moved,
                       "candidate_trials": trials,
                       "total_movement_distance": gr.rounded(movement),
                       "total_connection_length": gr.rounded(measured["total_connection_length"]),
                       "angular": angular_summary(graph, positions),
                       "width": extents["width"], "height": extents["height"],
                       "area": gr.rounded(extents["width"]*extents["height"]),
                       "validation": validation})
    final_validation, final_geometry = e4n.complete_validate(circuit, graph, embedding, positions)
    if not final_validation["valid"]:
        raise RuntimeError("final complete validation failed")
    final_metrics = e4i.measure(circuit, final_geometry[0], final_geometry[3])
    final_extents = e4j.geometry_extents(final_geometry[0], final_geometry[1])
    final_angular = angular_summary(graph, positions)
    four_k_angular = angular_summary(graph, initial_positions)
    four_n_positions = positions_from_coordinates(
        circuit, output_dir/"planar_local_compact_coordinates.json")
    four_n_angular = angular_summary(graph, four_n_positions)
    four_n_metrics = json.loads((output_dir/"planar_local_compact_metrics.json").read_text())
    comparisons = {
        "experiment_4k": {"width": 3166, "height": 1112, "area": 3520592,
                          "total_connection_length": 23095.744883, "angular": four_k_angular,
                          "runtime_seconds": None},
        "experiment_4n": {"width": 2401.838680, "height": 984.746882,
                          "area": 2365203.151197, "total_connection_length": 13979.580061,
                          "angular": four_n_angular,
                          "runtime_seconds": four_n_metrics["runtime_seconds"]}}
    elapsed = time.perf_counter()-started
    coordinate_path = output_dir/"planar_circular_compact_coordinates.json"
    components, nets, attachments, edges = final_geometry
    gr.write_json(coordinate_path, {"circuit": circuit.name,
        "source": "exact stored Experiment 4K geometry", "processing_order": order,
        "steps": list(e4n.STEPS), "comparable_reduction_fraction": COMPARABLE_REDUCTION_FRACTION,
        "component_dimensions": {"width": gr.WIDTH, "height": gr.HEIGHT},
        "graph_vertex_positions": {vertex: e4i.point_doc(positions[vertex]) for vertex in sorted(positions)},
        "components": {ref: e4i.point_doc(components[ref]) for ref in circuit.refs},
        "net_junctions": {name: e4i.point_doc(nets[name]) for name in sorted(nets)},
        "attachments": {terminal: e4i.point_doc(attachments[terminal]) for terminal in sorted(attachments)},
        "edges": [{"terminal": edge["terminal"], "component": edge["component"], "net": edge["net"],
                   "points": [e4i.point_doc(point) for point in edge["points"]]} for edge in edges]})
    metrics_path = output_dir/"planar_circular_compact_metrics.json"
    metrics_document = {"processing_order": order, "steps": list(e4n.STEPS),
        "comparable_reduction_fraction": COMPARABLE_REDUCTION_FRACTION,
        "passes": passes, "moves": moves, "total_accepted_moves": len(moves),
        "geometry_extents": final_extents,
        "drawing_area": gr.rounded(final_extents["width"]*final_extents["height"]),
        "metrics": gr.rounded_metrics(final_metrics), "angular": final_angular,
        "runtime_seconds": gr.rounded(elapsed), "comparison": comparisons}
    gr.write_json(metrics_path, metrics_document)
    validation_path = output_dir/"planar_circular_compact_validation.json"
    gr.write_json(validation_path, {"baseline_4k_reproduced": True,
        "embedding_unchanged": True, "exactly_six_passes": len(passes)==6,
        "complete_pass_validations": [row["validation"] for row in passes],
        "final_complete_validation": final_validation})
    svg_path = output_dir/"planar_circular_compact.svg"
    e4i.write_svg(svg_path, circuit, *final_geometry)
    paths = (coordinate_path, metrics_path, validation_path, svg_path)
    return {"metrics": metrics_document, "validation": final_validation,
            "hashes": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.input, args.output_dir), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
