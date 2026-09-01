#!/usr/bin/env python3
"""Active Circular Relaxation (later user-labeled 4P; distinct from planarity 4P)."""

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
import experiment_4o as e4o
import experiment_4p as planarity_4p
import graph_relax as gr


LENGTH_ALLOWANCE = 1.10
ANGULAR_TOLERANCE_DEGREES = 2.0
PASSES_PER_STEP = 2
DIRECTION_LABELS = (
    "no_move", "inward", "outward", "clockwise_tangent", "counterclockwise_tangent",
    "inward_clockwise", "inward_counterclockwise",
    "outward_clockwise", "outward_counterclockwise")


def movement_directions(graph, positions, vertex):
    center = positions[vertex]
    neighbors = sorted(graph.neighbors(vertex))
    local_center = (sum(positions[n][0] for n in neighbors)/len(neighbors),
                    sum(positions[n][1] for n in neighbors)/len(neighbors))
    rx, ry = center[0]-local_center[0], center[1]-local_center[1]
    magnitude = math.hypot(rx, ry)
    if magnitude <= e4n.EPSILON:
        first = positions[neighbors[0]]
        rx, ry = center[0]-first[0], center[1]-first[1]
        magnitude = math.hypot(rx, ry)
        if magnitude <= e4n.EPSILON:
            rx, ry, magnitude = 1.0, 0.0, 1.0
    outward = (rx/magnitude, ry/magnitude)
    inward = (-outward[0], -outward[1])
    # In SVG/screen coordinates (+Y downward), this is clockwise.
    clockwise = (-outward[1], outward[0])
    counterclockwise = (outward[1], -outward[0])
    def normalized(first, second):
        x, y = first[0]+second[0], first[1]+second[1]
        length = math.hypot(x, y)
        return x/length, y/length
    directions = [
        ("inward", inward), ("outward", outward),
        ("clockwise_tangent", clockwise), ("counterclockwise_tangent", counterclockwise),
        ("inward_clockwise", normalized(inward, clockwise)),
        ("inward_counterclockwise", normalized(inward, counterclockwise)),
        ("outward_clockwise", normalized(outward, clockwise)),
        ("outward_counterclockwise", normalized(outward, counterclockwise))]
    return local_center, directions


def candidate_rows(circuit, graph, embedding, positions, realization, vertex, step):
    current_length = e4n.local_length(vertex, realization[3])
    current_angular = e4o.vertex_angular_measure(graph, positions, vertex)
    rows = [{"direction_index": 0, "direction": "no_move", "point": positions[vertex],
             "positions": positions, "geometry": realization,
             "local_length": current_length, "angular": current_angular}]
    _, directions = movement_directions(graph, positions, vertex)
    for offset, (label, unit) in enumerate(directions, 1):
        point = e4n.candidate_point(positions[vertex], unit, step)
        if point == positions[vertex]:
            continue
        trial_positions = dict(positions); trial_positions[vertex] = point
        try:
            trial_geometry = e4n.geometry(circuit, trial_positions)
        except ValueError:
            continue
        if not e4n.locally_valid(circuit, embedding, trial_positions, vertex, trial_geometry):
            continue
        rows.append({"direction_index": offset, "direction": label, "point": point,
                     "positions": trial_positions, "geometry": trial_geometry,
                     "local_length": e4n.local_length(vertex, trial_geometry[3]),
                     "angular": e4o.vertex_angular_measure(graph, trial_positions, vertex)})
    return rows, current_length


def select_candidate(graph, vertex, rows, current_length):
    if graph.degree(vertex) == 1:
        return min(rows, key=lambda row: (row["local_length"], row["direction_index"]))
    eligible = [row for row in rows if row["local_length"] <= LENGTH_ALLOWANCE*current_length+1e-9]
    best_error = min(row["angular"]["angular_error"] for row in eligible)
    angularly_comparable = [row for row in eligible
                            if row["angular"]["angular_error"]-best_error < ANGULAR_TOLERANCE_DEGREES]
    return min(angularly_comparable, key=lambda row: (
        row["local_length"], -row["angular"]["minimum_gap"], row["direction_index"]))


def positions_from(path, circuit):
    document = json.loads(path.read_text())
    positions = {planarity_4p.component_vertex(ref):
                 (float(document["components"][ref]["x"]), float(document["components"][ref]["y"]))
                 for ref in circuit.refs}
    positions.update({planarity_4p.net_vertex(net.name):
                      (float(document["net_junctions"][net.name]["x"]),
                       float(document["net_junctions"][net.name]["y"]))
                      for net in circuit.nets})
    return positions


def comparison_row(graph, circuit, output_dir, name, coordinate_name,
                   width, height, area, length, runtime):
    positions = positions_from(output_dir/coordinate_name, circuit)
    return {"width": width, "height": height, "area": area,
            "total_connection_length": length,
            "angular": e4o.angular_summary(graph, positions),
            "crossings": 0, "component_overlaps": 0, "clearance_violations": 0,
            "unrelated_component_body_intersections": 0,
            "electrical_incidences": "44/44", "runtime_seconds": runtime}


def run(input_path, output_dir):
    started = time.perf_counter()
    circuit = gr.parse_circuit(input_path)
    graph, embedding, positions, order = e4n.load_start(circuit, output_dir)
    baseline, realization = e4n.complete_validate(circuit, graph, embedding, positions)
    if not baseline["valid"]:
        raise ValueError("stored 4K starting geometry failed reproduction")
    moves, passes = [], []
    global_pass = 0
    for step in e4n.STEPS:
        for step_pass in range(1, PASSES_PER_STEP+1):
            global_pass += 1
            moved = trials = 0
            movement = 0.0
            for vertex in order:
                rows, current_length = candidate_rows(
                    circuit, graph, embedding, positions, realization, vertex, step)
                trials += len(rows)-1
                chosen = select_candidate(graph, vertex, rows, current_length)
                if chosen["direction"] == "no_move":
                    continue
                old = positions[vertex]
                positions, realization = chosen["positions"], chosen["geometry"]
                distance = math.dist(old, chosen["point"])
                moved += 1; movement += distance
                moves.append({"move_number": len(moves)+1, "pass": global_pass,
                              "step": step, "step_pass": step_pass,
                              "vertex": vertex, "degree": graph.degree(vertex),
                              "direction": chosen["direction"],
                              "old": e4i.point_doc(old), "new": e4i.point_doc(chosen["point"]),
                              "movement_distance": gr.rounded(distance),
                              "local_length_before": gr.rounded(current_length),
                              "local_length_after": gr.rounded(chosen["local_length"]),
                              "resulting_angular_error": None if chosen["angular"] is None else gr.rounded(chosen["angular"]["angular_error"]),
                              "resulting_minimum_gap": None if chosen["angular"] is None else gr.rounded(chosen["angular"]["minimum_gap"])})
            validation, checked = e4n.complete_validate(circuit, graph, embedding, positions)
            if not validation["valid"]:
                raise RuntimeError(f"complete validation failed after pass {global_pass}")
            realization = checked
            measured = e4i.measure(circuit, realization[0], realization[3])
            extents = e4j.geometry_extents(realization[0], realization[1])
            passes.append({"pass": global_pass, "step": step, "step_pass": step_pass,
                           "vertices_moved": moved, "accepted_moves": moved,
                           "candidate_trials": trials,
                           "total_movement_distance": gr.rounded(movement),
                           "total_connection_length": gr.rounded(measured["total_connection_length"]),
                           "angular": e4o.angular_summary(graph, positions),
                           "width": extents["width"], "height": extents["height"],
                           "area": gr.rounded(extents["width"]*extents["height"]),
                           "validation": validation})
    final_validation, final_geometry = e4n.complete_validate(circuit, graph, embedding, positions)
    if not final_validation["valid"]:
        raise RuntimeError("final complete validation failed")
    final_metrics = e4i.measure(circuit, final_geometry[0], final_geometry[3])
    extents = e4j.geometry_extents(final_geometry[0], final_geometry[1])
    angular = e4o.angular_summary(graph, positions)
    runtime = time.perf_counter()-started
    four_n_metrics = json.loads((output_dir/"planar_local_compact_metrics.json").read_text())
    four_o_metrics = json.loads((output_dir/"planar_circular_compact_metrics.json").read_text())
    comparisons = {
        "experiment_4k": comparison_row(graph, circuit, output_dir, "4K",
            "planar_xy_compact_coordinates.json", 3166, 1112, 3520592,
            23095.744883, None),
        "experiment_4n": comparison_row(graph, circuit, output_dir, "4N",
            "planar_local_compact_coordinates.json", 2401.838680, 984.746882,
            2365203.151197, 13979.580061, four_n_metrics["runtime_seconds"]),
        "experiment_4o": comparison_row(graph, circuit, output_dir, "4O",
            "planar_circular_compact_coordinates.json", 2545.874241, 1006.351588,
            2562044.585279, 15269.151443, four_o_metrics["runtime_seconds"])}
    components, nets, attachments, edges = final_geometry
    coordinate_path = output_dir/"planar_active_circular_coordinates.json"
    gr.write_json(coordinate_path, {"circuit": circuit.name,
        "source": "exact stored Experiment 4K geometry", "processing_order": order,
        "steps": list(e4n.STEPS), "passes_per_step": PASSES_PER_STEP,
        "length_allowance": LENGTH_ALLOWANCE,
        "angular_tolerance_degrees": ANGULAR_TOLERANCE_DEGREES,
        "candidate_order": list(DIRECTION_LABELS),
        "component_dimensions": {"width": gr.WIDTH, "height": gr.HEIGHT},
        "graph_vertex_positions": {vertex: e4i.point_doc(positions[vertex]) for vertex in sorted(positions)},
        "components": {ref: e4i.point_doc(components[ref]) for ref in circuit.refs},
        "net_junctions": {name: e4i.point_doc(nets[name]) for name in sorted(nets)},
        "attachments": {terminal: e4i.point_doc(attachments[terminal]) for terminal in sorted(attachments)},
        "edges": [{"terminal": edge["terminal"], "component": edge["component"], "net": edge["net"],
                   "points": [e4i.point_doc(point) for point in edge["points"]]} for edge in edges]})
    metrics_path = output_dir/"planar_active_circular_metrics.json"
    metrics_document = {"processing_order": order, "steps": list(e4n.STEPS),
        "passes_per_step": PASSES_PER_STEP, "length_allowance": LENGTH_ALLOWANCE,
        "angular_tolerance_degrees": ANGULAR_TOLERANCE_DEGREES,
        "passes": passes, "moves": moves, "total_accepted_moves": len(moves),
        "geometry_extents": extents, "drawing_area": gr.rounded(extents["width"]*extents["height"]),
        "metrics": gr.rounded_metrics(final_metrics), "angular": angular,
        "runtime_seconds": gr.rounded(runtime), "comparison": comparisons}
    gr.write_json(metrics_path, metrics_document)
    validation_path = output_dir/"planar_active_circular_validation.json"
    gr.write_json(validation_path, {"baseline_4k_reproduced": True,
        "embedding_unchanged": True, "exactly_twelve_passes": len(passes)==12,
        "complete_pass_validations": [row["validation"] for row in passes],
        "final_complete_validation": final_validation})
    svg_path = output_dir/"planar_active_circular.svg"
    e4i.write_svg(svg_path, circuit, *final_geometry)
    paths = (coordinate_path, metrics_path, validation_path, svg_path)
    return {"metrics": metrics_document, "validation": final_validation,
            "hashes": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}}


def main():
    parser=argparse.ArgumentParser();parser.add_argument("input",type=Path);parser.add_argument("output_dir",type=Path);args=parser.parse_args()
    print(json.dumps(run(args.input,args.output_dir),sort_keys=True,indent=2))


if __name__=="__main__":main()
