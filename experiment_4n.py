#!/usr/bin/env python3
"""Experiment 4N: six coarse local incident-length compaction passes."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

import networkx as nx

import experiment_4i as e4i
import experiment_4j as e4j
import experiment_4p as e4p
import graph_relax as gr


STEPS = (256.0, 128.0, 64.0, 32.0, 16.0, 8.0)
EPSILON = 1e-7


def load_start(circuit, output_dir):
    embedding_document = json.loads((output_dir / "node_planar_embedding.json").read_text())
    coordinates = json.loads((output_dir / "planar_xy_compact_coordinates.json").read_text())
    if (coordinates["x_scale"], coordinates["y_scale"]) != (49, 52):
        raise ValueError("Experiment 4N requires the exact stored 4K (49,52) geometry")
    graph = e4p.incidence_graph(circuit)
    if e4p.graph_document(graph) != embedding_document["graph"]:
        raise ValueError("stored embedding graph differs from source graph")
    embedding = nx.PlanarEmbedding()
    embedding.set_data(embedding_document["cyclic_neighbor_order_clockwise"])
    embedding.check_structure()
    positions = {e4p.component_vertex(ref):
                 (float(coordinates["components"][ref]["x"]),
                  float(coordinates["components"][ref]["y"]))
                 for ref in circuit.refs}
    positions.update({e4p.net_vertex(net.name):
                      (float(coordinates["net_junctions"][net.name]["x"]),
                       float(coordinates["net_junctions"][net.name]["y"]))
                      for net in circuit.nets})
    order = sorted(graph, key=lambda vertex: (-graph.degree(vertex), vertex))
    return graph, embedding, positions, order


def geometry(circuit, positions):
    components = {ref: positions[e4p.component_vertex(ref)] for ref in circuit.refs}
    nets = {net.name: positions[e4p.net_vertex(net.name)] for net in circuit.nets}
    attachments, edges = {}, []
    for net in circuit.nets:
        for terminal in net.terminals:
            ref = terminal.rsplit(".", 1)[0]
            attachment = e4i.ray_rectangle_attachment(components[ref], nets[net.name])
            attachments[terminal] = attachment
            edges.append({"terminal": terminal, "component": ref, "net": net.name,
                          "points": [attachment, nets[net.name]]})
    return components, nets, attachments, edges


def cyclic_orientation(embedding, positions):
    states = []
    for vertex in sorted(embedding):
        required = list(embedding.neighbors_cw_order(vertex))
        if len(required) < 3:
            continue
        center = positions[vertex]
        geometric = sorted(required, key=lambda neighbor: (
            math.atan2(positions[neighbor][1]-center[1], positions[neighbor][0]-center[0]), neighbor))
        if e4i.cyclic_equal(geometric, required):
            states.append("same")
        elif e4i.cyclic_equal(geometric, list(reversed(required))):
            states.append("global_reflection")
        else:
            return "inconsistent_local_reversal"
    return states[0] if states and len(set(states)) == 1 else "inconsistent_local_reversal"


def complete_validate(circuit, graph, embedding, positions):
    try:
        realization = geometry(circuit, positions)
    except ValueError as error:
        return {"valid": False, "construction_error": str(error)}, None
    result = e4i.validate_geometry(circuit, graph, embedding, *realization)
    orientation = cyclic_orientation(embedding, positions)
    result["all_vertex_cyclic_order_orientation"] = orientation
    result["missing_incidences"] = 0 if result["electrical_incidence_exact"] else None
    result["added_incidences"] = 0 if result["electrical_incidence_exact"] else None
    result["valid"] = bool(result["valid"] and orientation in ("same", "global_reflection") and
                           result["expected_incidence_count"] == 44 and
                           result["realized_incidence_count"] == 44)
    return result, realization


def affected_edge_indexes(edges, vertex):
    kind, name = vertex.split(":", 1)
    if kind == "component":
        return {index for index, edge in enumerate(edges) if edge["component"] == name}
    return {index for index, edge in enumerate(edges) if edge["net"] == name}


def locally_valid(circuit, embedding, positions, moved_vertex, realization):
    components, _, _, edges = realization
    kind, name = moved_vertex.split(":", 1)
    if kind == "component":
        for other in circuit.refs:
            if other == name:
                continue
            if gr.overlap_area(components[name], components[other]) > 0:
                return False
            if gr.rect_gap(components[name], components[other]) < gr.BODY_CLEARANCE:
                return False
    affected = affected_edge_indexes(edges, moved_vertex)
    # Only pairs containing an affected edge can have changed crossing status.
    tested_pairs = set()
    for first_index in sorted(affected):
        first = edges[first_index]
        for second_index, second in enumerate(edges):
            if first_index == second_index or first["net"] == second["net"]:
                continue
            pair = tuple(sorted((first_index, second_index)))
            if pair in tested_pairs:
                continue
            tested_pairs.add(pair)
            if gr.proper_intersection(*first["points"], *second["points"]):
                return False
    # Changed edges must avoid every unrelated body.
    for index in sorted(affected):
        edge = edges[index]
        for ref in circuit.refs:
            if ref != edge["component"] and gr.segment_rect_distance(
                    *edge["points"], gr.rect_bounds(components[ref])) <= EPSILON:
                return False
    # A moved component body must also avoid all unchanged edges.
    if kind == "component":
        bounds = gr.rect_bounds(components[name])
        for index, edge in enumerate(edges):
            if index not in affected and edge["component"] != name and gr.segment_rect_distance(
                    *edge["points"], bounds) <= EPSILON:
                return False
    return cyclic_orientation(embedding, positions) in ("same", "global_reflection")


def local_length(vertex, edges):
    indexes = affected_edge_indexes(edges, vertex)
    return sum(math.dist(*edges[index]["points"]) for index in indexes)


def candidate_directions(graph, positions, vertex):
    center = positions[vertex]
    raw = []
    for neighbor in sorted(graph.neighbors(vertex)):
        raw.append((f"neighbor:{neighbor}", positions[neighbor], True))
    neighbors = sorted(graph.neighbors(vertex))
    centroid = (sum(positions[n][0] for n in neighbors)/len(neighbors),
                sum(positions[n][1] for n in neighbors)/len(neighbors))
    raw.append(("neighbor_centroid", centroid, True))
    raw.extend((("left", (center[0]-1.0, center[1]), False),
                ("right", (center[0]+1.0, center[1]), False),
                ("up", (center[0], center[1]-1.0), False),
                ("down", (center[0], center[1]+1.0), False)))
    result, seen = [], set()
    for label, target, bounded in raw:
        dx, dy = target[0]-center[0], target[1]-center[1]
        distance = math.hypot(dx, dy)
        if distance <= EPSILON:
            continue
        unit = (dx/distance, dy/distance)
        key = (round(unit[0], 12), round(unit[1], 12))
        if key in seen:
            continue
        seen.add(key)
        result.append({"label": label, "unit": unit,
                       "target_distance": distance if bounded else None})
    return result


def candidate_point(center, unit, step):
    return gr.rounded(center[0]+unit[0]*step), gr.rounded(center[1]+unit[1]*step)


def best_move(circuit, graph, embedding, positions, realization, vertex, step):
    original_length = local_length(vertex, realization[3])
    best = None
    trials = 0
    for direction_index, direction in enumerate(candidate_directions(graph, positions, vertex)):
        if direction["target_distance"] is not None and step > direction["target_distance"] + EPSILON:
            continue
        point = candidate_point(positions[vertex], direction["unit"], step)
        if point == positions[vertex]:
            continue
        trial_positions = dict(positions); trial_positions[vertex] = point
        try:
            trial_geometry = geometry(circuit, trial_positions)
        except ValueError:
            continue
        trials += 1
        if not locally_valid(circuit, embedding, trial_positions, vertex, trial_geometry):
            continue
        new_length = local_length(vertex, trial_geometry[3])
        reduction = original_length-new_length
        if reduction <= 1e-9:
            continue
        if best is None or reduction > best[0]+1e-9:
            best = (reduction, direction_index, direction["label"], point,
                    trial_positions, trial_geometry, original_length, new_length)
    return best, trials


def run(input_path, output_dir):
    started = time.perf_counter()
    circuit = gr.parse_circuit(input_path)
    graph, embedding, positions, order = load_start(circuit, output_dir)
    baseline, realization = complete_validate(circuit, graph, embedding, positions)
    if not baseline["valid"]:
        raise ValueError("stored 4K starting geometry failed reproduction")
    moves, passes = [], []
    for pass_number, step in enumerate(STEPS, 1):
        moved, movement, trials = 0, 0.0, 0
        for vertex in order:
            old_point = positions[vertex]
            best, vertex_trials = best_move(circuit, graph, embedding, positions,
                                            realization, vertex, step)
            trials += vertex_trials
            if best is None:
                continue
            reduction, _, label, point, positions, realization, before, after = best
            actual_movement = math.dist(old_point, point)
            moved += 1; movement += actual_movement
            moves.append({"move_number": len(moves)+1, "pass": pass_number, "step": step,
                          "vertex": vertex, "direction": label,
                          "old": e4i.point_doc(old_point),
                          "new": e4i.point_doc(point), "local_length_before": gr.rounded(before),
                          "movement_distance": gr.rounded(actual_movement),
                          "local_length_after": gr.rounded(after), "local_length_reduction": gr.rounded(reduction)})
        validation, validated_geometry = complete_validate(circuit, graph, embedding, positions)
        if not validation["valid"]:
            raise RuntimeError(f"complete validation failed after step {step}")
        realization = validated_geometry
        metrics = e4i.measure(circuit, realization[0], realization[3])
        extents = e4j.geometry_extents(realization[0], realization[1])
        passes.append({"pass": pass_number, "step": step, "vertices_moved": moved,
                       "accepted_moves": moved, "candidate_trials": trials,
                       "total_movement_distance": gr.rounded(movement),
                       "total_connection_length": gr.rounded(metrics["total_connection_length"]),
                       "width": extents["width"], "height": extents["height"],
                       "area": gr.rounded(extents["width"]*extents["height"]),
                       "validation": validation})
    final_validation, final_geometry = complete_validate(circuit, graph, embedding, positions)
    if not final_validation["valid"]:
        raise RuntimeError("final complete validation failed")
    final_metrics = e4i.measure(circuit, final_geometry[0], final_geometry[3])
    extents = e4j.geometry_extents(final_geometry[0], final_geometry[1])
    elapsed = time.perf_counter()-started
    coordinate_path = output_dir / "planar_local_compact_coordinates.json"
    components, nets, attachments, edges = final_geometry
    gr.write_json(coordinate_path, {"circuit": circuit.name,
        "source": "exact stored Experiment 4K geometry", "processing_order": order,
        "steps": list(STEPS), "component_dimensions": {"width": gr.WIDTH, "height": gr.HEIGHT},
        "graph_vertex_positions": {vertex: e4i.point_doc(positions[vertex]) for vertex in sorted(positions)},
        "components": {ref: e4i.point_doc(components[ref]) for ref in circuit.refs},
        "net_junctions": {name: e4i.point_doc(nets[name]) for name in sorted(nets)},
        "attachments": {terminal: e4i.point_doc(attachments[terminal]) for terminal in sorted(attachments)},
        "edges": [{"terminal": edge["terminal"], "component": edge["component"], "net": edge["net"],
                   "points": [e4i.point_doc(point) for point in edge["points"]]} for edge in edges]})
    metrics_path = output_dir / "planar_local_compact_metrics.json"
    metrics_document = {"processing_order": order, "steps": list(STEPS), "passes": passes,
                        "moves": moves, "total_accepted_moves": len(moves),
                        "geometry_extents": extents,
                        "drawing_area": gr.rounded(extents["width"]*extents["height"]),
                        "metrics": gr.rounded_metrics(final_metrics),
                        "runtime_seconds": gr.rounded(elapsed),
                        "comparison_4k": {"width": 3166, "height": 1112, "area": 3520592,
                                          "total_connection_length": 23095.744883,
                                          "crossings": 0, "overlaps": 0,
                                          "clearance_violations": 0, "body_intersections": 0}}
    gr.write_json(metrics_path, metrics_document)
    validation_path = output_dir / "planar_local_compact_validation.json"
    gr.write_json(validation_path, {"baseline_4k_reproduced": True,
                                    "embedding_unchanged": True,
                                    "exactly_six_passes": len(passes) == 6,
                                    "complete_pass_validations": [row["validation"] for row in passes],
                                    "final_complete_validation": final_validation})
    svg_path = output_dir / "planar_local_compact.svg"
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
