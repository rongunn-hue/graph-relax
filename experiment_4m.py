#!/usr/bin/env python3
"""Experiment 4M: greedily pull frozen 4I vertices toward maximum-degree U1."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import networkx as nx

import experiment_4i as e4i
import experiment_4j as e4j
import experiment_4p as e4p
import graph_relax as gr


MAX_PASSES = 20
COORDINATE_PRECISION = 10 ** (-gr.ROUND_DIGITS)


def load_frozen(circuit, output_dir):
    embedding_document = json.loads((output_dir / "node_planar_embedding.json").read_text())
    coordinates = json.loads((output_dir / "planar_node_initial_coordinates.json").read_text())
    if coordinates["scale"] != 72:
        raise ValueError("Experiment 4M requires the exact scale-72 4I realization")
    graph = e4p.incidence_graph(circuit)
    if e4p.graph_document(graph) != embedding_document["graph"]:
        raise ValueError("stored 4I graph differs from source graph")
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
    return graph, embedding, embedding_document, positions


def anchor_and_order(graph):
    ordered = sorted(graph, key=lambda vertex: (-graph.degree(vertex), vertex))
    anchor = ordered[0]
    return anchor, graph.degree(anchor), [vertex for vertex in ordered if vertex != anchor]


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


def cyclic_equal(first, second):
    return e4i.cyclic_equal(first, second)


def all_vertex_cyclic_orientation(embedding, positions):
    rows, states = [], []
    for vertex in sorted(embedding):
        required = list(embedding.neighbors_cw_order(vertex))
        if len(required) < 3:
            state = "ambiguous_degree_below_3"
            geometric = sorted(required)
        else:
            center = positions[vertex]
            geometric = sorted(required, key=lambda neighbor: (
                math.atan2(positions[neighbor][1]-center[1], positions[neighbor][0]-center[0]),
                neighbor))
            if cyclic_equal(geometric, required):
                state = "same"
            elif cyclic_equal(geometric, list(reversed(required))):
                state = "global_reflection"
            else:
                state = "inconsistent"
            states.append(state)
        rows.append({"vertex": vertex, "required_clockwise": required,
                     "geometric_angular_order": geometric, "state": state})
    unique = set(states)
    orientation = (next(iter(unique)) if len(unique) == 1 and
                   unique <= {"same", "global_reflection"} else
                   "inconsistent_local_reversal")
    return orientation, rows


def validate(circuit, graph, embedding, positions):
    try:
        realization = geometry(circuit, positions)
    except ValueError as error:
        return {"valid": False, "construction_error": str(error)}, None
    base = e4i.validate_geometry(circuit, graph, embedding, *realization)
    orientation, rows = all_vertex_cyclic_orientation(embedding, positions)
    base["all_vertex_cyclic_order_orientation"] = orientation
    base["all_vertex_cyclic_order_rows"] = rows
    base["missing_incidences"] = 0 if base["electrical_incidence_exact"] else None
    base["added_incidences"] = 0 if base["electrical_incidence_exact"] else None
    base["valid"] = bool(base["valid"] and orientation in ("same", "global_reflection") and
                         base["expected_incidence_count"] == 44 and
                         base["realized_incidence_count"] == 44)
    return base, realization


def rounded_point(point):
    return (gr.rounded(point[0]), gr.rounded(point[1]))


def candidate_position(start, anchor, fraction):
    return rounded_point((start[0] + fraction * (anchor[0]-start[0]),
                          start[1] + fraction * (anchor[1]-start[1])))


def furthest_continuously_valid(circuit, graph, embedding, positions, vertex, anchor):
    start, target = positions[vertex], positions[anchor]
    distance = math.dist(start, target)
    if distance <= COORDINATE_PRECISION:
        return start, 0.0, 0
    evaluations = 0
    def check(fraction):
        nonlocal evaluations
        evaluations += 1
        trial = dict(positions)
        trial[vertex] = candidate_position(start, target, fraction)
        return validate(circuit, graph, embedding, trial)[0]["valid"], trial[vertex]
    endpoint_valid, endpoint = check(1.0)
    if endpoint_valid:
        return endpoint, math.dist(start, endpoint), evaluations
    low, high, best = 0.0, 1.0, start
    # Find the end of the connected valid interval containing the current point.
    while distance * (high-low) > COORDINATE_PRECISION:
        middle = (low+high)/2.0
        valid, point = check(middle)
        if valid:
            low, best = middle, point
        else:
            high = middle
    # Rounding can repeat candidates; independently certify the returned state.
    trial = dict(positions); trial[vertex] = best
    if not validate(circuit, graph, embedding, trial)[0]["valid"]:
        raise RuntimeError("selected path point failed final independent validation")
    movement = math.dist(start, best)
    return best, movement if movement >= COORDINATE_PRECISION else 0.0, evaluations


def compact(circuit, graph, embedding, initial_positions, anchor, order):
    positions = dict(initial_positions)
    moves, pass_rows = [], []
    for pass_number in range(1, MAX_PASSES+1):
        pass_moves = 0
        for vertex in order:
            old = positions[vertex]
            point, distance, evaluations = furthest_continuously_valid(
                circuit, graph, embedding, positions, vertex, anchor)
            if distance > 0:
                positions[vertex] = point
                pass_moves += 1
                moves.append({"move_number": len(moves)+1, "pass": pass_number,
                              "vertex": vertex, "degree": graph.degree(vertex),
                              "old": e4i.point_doc(old), "new": e4i.point_doc(point),
                              "movement_distance": gr.rounded(distance),
                              "path_evaluations": evaluations})
        pass_rows.append({"pass": pass_number, "vertices_moved": pass_moves})
        if pass_moves == 0:
            return positions, moves, pass_rows, "complete pass produced no movement"
    return positions, moves, pass_rows, "20-pass safety limit reached"


def coordinate_document(circuit, positions, realization, validation, anchor, order):
    components, nets, attachments, edges = realization
    return {"circuit": circuit.name, "source": "exact stored Experiment 4I scale-72 realization",
            "anchor": anchor, "processing_order": order,
            "component_dimensions": {"width": gr.WIDTH, "height": gr.HEIGHT},
            "graph_vertex_positions": {vertex: e4i.point_doc(positions[vertex]) for vertex in sorted(positions)},
            "components": {ref: e4i.point_doc(components[ref]) for ref in circuit.refs},
            "net_junctions": {name: e4i.point_doc(nets[name]) for name in sorted(nets)},
            "attachments": {terminal: e4i.point_doc(attachments[terminal]) for terminal in sorted(attachments)},
            "edges": [{"terminal": edge["terminal"], "component": edge["component"], "net": edge["net"],
                       "points": [e4i.point_doc(point) for point in edge["points"]]} for edge in edges],
            "cyclic_order_orientation": validation["all_vertex_cyclic_order_orientation"]}


def percent_change(new, old):
    return gr.rounded(100*(new-old)/old)


def run(input_path, output_dir):
    circuit = gr.parse_circuit(input_path)
    graph, embedding, _, initial_positions = load_frozen(circuit, output_dir)
    baseline_validation, baseline_geometry = validate(circuit, graph, embedding, initial_positions)
    if not baseline_validation["valid"]:
        raise ValueError("stored Experiment 4I geometry does not reproduce as valid")
    anchor, anchor_degree, order = anchor_and_order(graph)
    positions, moves, passes, stopping_reason = compact(
        circuit, graph, embedding, initial_positions, anchor, order)
    validation, realization = validate(circuit, graph, embedding, positions)
    if not validation["valid"]:
        raise RuntimeError("final anchor-compacted realization is invalid")
    metrics = e4i.measure(circuit, realization[0], realization[3])
    extents = e4j.geometry_extents(realization[0], realization[1])
    area = extents["width"]*extents["height"]
    coordinate_path = output_dir / "planar_anchor_compact_coordinates.json"
    gr.write_json(coordinate_path, coordinate_document(
        circuit, positions, realization, validation, anchor, order))
    metrics_path = output_dir / "planar_anchor_compact_metrics.json"
    metrics_document = {"anchor": anchor, "anchor_degree": anchor_degree,
                        "processing_order": order, "passes": passes,
                        "pass_count": len(passes), "stopping_reason": stopping_reason,
                        "moved_vertex_count": len(set(move["vertex"] for move in moves)),
                        "accepted_move_count": len(moves), "movements": moves,
                        "geometry_extents": extents, "drawing_area": gr.rounded(area),
                        "metrics": gr.rounded_metrics(metrics),
                        "comparison_4k": {"width": 3166, "height": 1112, "area": 3520592,
                                          "total_connection_length": 23095.744883},
                        "change_from_4k_percent": {
                            "width": percent_change(extents["width"], 3166),
                            "height": percent_change(extents["height"], 1112),
                            "area": percent_change(area, 3520592),
                            "connection_length": percent_change(metrics["total_connection_length"], 23095.744883)}}
    gr.write_json(metrics_path, metrics_document)
    validation_path = output_dir / "planar_anchor_compact_validation.json"
    gr.write_json(validation_path, {"baseline_4i_reproduced": True,
                                    "embedding_unchanged": True,
                                    "anchor_fixed": positions[anchor] == initial_positions[anchor],
                                    "coordinate_precision": COORDINATE_PRECISION,
                                    "maximum_passes": MAX_PASSES,
                                    "validation": validation})
    svg_path = output_dir / "planar_anchor_compact.svg"
    e4i.write_svg(svg_path, circuit, *realization)
    paths = (coordinate_path, metrics_path, validation_path, svg_path)
    return {"anchor": anchor, "anchor_degree": anchor_degree,
            "pass_count": len(passes), "stopping_reason": stopping_reason,
            "moved_vertex_count": metrics_document["moved_vertex_count"],
            "metrics": metrics_document,
            "hashes": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.input, args.output_dir), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
