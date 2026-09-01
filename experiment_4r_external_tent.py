#!/usr/bin/env python3
"""Experiment 4R: external-pole contraction of a pure geometric graph."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import time
from pathlib import Path

import networkx as nx

import experiment_4p as graph_authority
import graph_relax as gr


SWEEP_FRACTIONS = (0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60,
                   0.70, 0.80, 0.90, 0.95, 0.98)
SVG_FRACTIONS = (0.00, 0.10, 0.25, 0.50, 0.75, 0.90, 0.98)
TOLERANCE = 1e-7


def point_doc(point):
    return {"x": gr.rounded(point[0]), "y": gr.rounded(point[1])}


def load_abstract_start(circuit, output_dir):
    embedding_doc = json.loads((output_dir / "node_planar_embedding.json").read_text())
    coordinate_doc = json.loads((output_dir / "planar_xy_compact_coordinates.json").read_text())
    if (coordinate_doc.get("x_scale"), coordinate_doc.get("y_scale")) != (49, 52):
        raise ValueError("stored Experiment 4K scale pair is not (49, 52)")
    graph = graph_authority.incidence_graph(circuit)
    expected_vertices = set(graph)
    expected_edges = {tuple(sorted(edge)) for edge in graph.edges()}
    stored_vertices = {vertex["id"] for vertex in embedding_doc["graph"]["vertices"]}
    stored_edges = {tuple(sorted((edge["u"], edge["v"]))) for edge in embedding_doc["graph"]["edges"]}
    if stored_vertices != expected_vertices or stored_edges != expected_edges:
        raise ValueError("stored 4I embedding does not match authoritative incidence graph")
    embedding = nx.PlanarEmbedding()
    embedding.set_data({vertex: embedding_doc["cyclic_neighbor_order_clockwise"][vertex]
                        for vertex in sorted(graph)})
    embedding.check_structure()
    positions = {}
    for ref, point in coordinate_doc["components"].items():
        positions[graph_authority.component_vertex(ref)] = (point["x"], point["y"])
    for net, point in coordinate_doc["net_junctions"].items():
        positions[graph_authority.net_vertex(net)] = (point["x"], point["y"])
    if set(positions) != expected_vertices or graph.number_of_nodes() != 34 or graph.number_of_edges() != 44:
        raise ValueError("stored 4K point geometry does not reproduce the 34-vertex/44-edge graph")
    return graph, embedding, positions, tuple(sorted(expected_edges))


def external_model(initial):
    vertices = sorted(initial)
    center = (sum(initial[v][0] for v in vertices) / len(vertices),
              sum(initial[v][1] for v in vertices) / len(vertices))
    distances = {v: math.dist(initial[v], center) for v in vertices}
    farthest = sorted(vertices, key=lambda v: (-distances[v], v))[0]
    radius = distances[farthest]
    direction = ((initial[farthest][0] - center[0]) / radius,
                 (initial[farthest][1] - center[1]) / radius)
    pole = (center[0] + 2 * radius * direction[0],
            center[1] + 2 * radius * direction[1])
    cables = {v: math.dist(initial[v], pole) for v in vertices}
    rays = {v: ((initial[v][0] - pole[0]) / cables[v],
                (initial[v][1] - pole[1]) / cables[v]) for v in vertices}
    hmax, owner = min((length, vertex) for vertex, length in cables.items())
    return center, farthest, radius, pole, cables, rays, hmax, owner


def positions_at_height(initial, model, height):
    _, _, _, pole, cables, rays, hmax, _ = model
    if not 0 <= height < hmax:
        raise ValueError("height outside strict fixed-cable limit")
    result = {}
    for vertex in sorted(initial):
        planar_radius = math.sqrt(cables[vertex] ** 2 - height ** 2)
        result[vertex] = (
            gr.rounded(pole[0] + rays[vertex][0] * planar_radius),
            gr.rounded(pole[1] + rays[vertex][1] * planar_radius),
        )
    return result


def orientation(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def point_on_segment_interior(point, a, b, tolerance=TOLERANCE):
    if abs(orientation(a, b, point)) > tolerance:
        return False
    dot = (point[0] - a[0]) * (point[0] - b[0]) + (point[1] - a[1]) * (point[1] - b[1])
    return dot < -tolerance


def segments_cross_or_overlap(a, b, c, d, tolerance=TOLERANCE):
    o1, o2, o3, o4 = orientation(a, b, c), orientation(a, b, d), orientation(c, d, a), orientation(c, d, b)
    if ((o1 > tolerance and o2 < -tolerance) or (o1 < -tolerance and o2 > tolerance)) and \
       ((o3 > tolerance and o4 < -tolerance) or (o3 < -tolerance and o4 > tolerance)):
        return True
    if all(abs(value) <= tolerance for value in (o1, o2, o3, o4)):
        axis = 0 if abs(a[0] - b[0]) >= abs(a[1] - b[1]) else 1
        overlap = min(max(a[axis], b[axis]), max(c[axis], d[axis])) - max(min(a[axis], b[axis]), min(c[axis], d[axis]))
        return overlap > tolerance
    return False


def cyclic_orientation(graph, embedding, positions):
    results = {}
    violations = []
    for vertex in sorted(graph):
        if graph.degree(vertex) < 3:
            continue
        geometric = sorted(graph.neighbors(vertex), key=lambda n: math.atan2(
            positions[n][1] - positions[vertex][1], positions[n][0] - positions[vertex][0]))
        required = list(embedding.neighbors_cw_order(vertex))
        def rotations(sequence):
            return [sequence[i:] + sequence[:i] for i in range(len(sequence))]
        if geometric in rotations(required):
            results[vertex] = "direct"
        elif geometric in rotations(list(reversed(required))):
            results[vertex] = "reflection"
        else:
            results[vertex] = "violation"
            violations.append(vertex)
    orientations = {value for value in results.values() if value != "violation"}
    global_orientation = next(iter(orientations)) if len(orientations) == 1 and not violations else None
    if len(orientations) > 1:
        violations.extend(sorted(results))
    return global_orientation, sorted(set(violations)), results


def validate_abstract(graph, embedding, positions, expected_vertices, expected_edges):
    vertices = sorted(positions)
    coincident = []
    for index, first in enumerate(vertices):
        for second in vertices[index + 1:]:
            if math.dist(positions[first], positions[second]) <= TOLERANCE:
                coincident.append([first, second])
    crossings = []
    edge_list = sorted(expected_edges)
    for index, edge_a in enumerate(edge_list):
        for edge_b in edge_list[index + 1:]:
            if set(edge_a) & set(edge_b):
                continue
            if segments_cross_or_overlap(positions[edge_a[0]], positions[edge_a[1]],
                                         positions[edge_b[0]], positions[edge_b[1]]):
                crossings.append([list(edge_a), list(edge_b)])
    vertex_on_edge = []
    for vertex in vertices:
        for edge in edge_list:
            if vertex in edge:
                continue
            if point_on_segment_interior(positions[vertex], positions[edge[0]], positions[edge[1]]):
                vertex_on_edge.append({"vertex": vertex, "edge": list(edge)})
    global_orientation, cyclic_violations, per_vertex = cyclic_orientation(graph, embedding, positions)
    exact_vertices = set(positions) == set(expected_vertices)
    exact_edges = {tuple(sorted(edge)) for edge in graph.edges()} == set(expected_edges)
    valid = exact_vertices and exact_edges and not crossings and not coincident and not vertex_on_edge and not cyclic_violations
    return {
        "all_original_vertices_exist": exact_vertices,
        "all_original_edges_preserved": exact_edges,
        "unrelated_edge_crossings": len(crossings),
        "crossing_details": crossings,
        "coincident_vertex_pairs": len(coincident),
        "coincident_details": coincident,
        "vertices_on_unrelated_edges": len(vertex_on_edge),
        "vertex_on_edge_details": vertex_on_edge,
        "cyclic_order_violations": len(cyclic_violations),
        "cyclic_violation_vertices": cyclic_violations,
        "cyclic_order_orientation": global_orientation,
        "per_vertex_cyclic_orientation": per_vertex,
        "abstract_graph_valid": valid,
    }


def evaluate(graph, embedding, initial, edges, model, fraction):
    height = fraction * model[6]
    positions = positions_at_height(initial, model, height)
    validation = validate_abstract(graph, embedding, positions, initial, edges)
    xs = [point[0] for point in positions.values()]
    ys = [point[1] for point in positions.values()]
    width, drawing_height = max(xs) - min(xs), max(ys) - min(ys)
    centroid = (sum(xs) / len(xs), sum(ys) / len(ys))
    edge_length = sum(math.dist(positions[a], positions[b]) for a, b in edges)
    return {
        "fraction": gr.rounded(fraction),
        "h": gr.rounded(height),
        "width": gr.rounded(width),
        "height": gr.rounded(drawing_height),
        "area": gr.rounded(width * drawing_height),
        "total_edge_length": gr.rounded(edge_length),
        "centroid_x": gr.rounded(centroid[0]),
        "centroid_y": gr.rounded(centroid[1]),
        "centroid_displacement": gr.rounded(math.dist(model[0], centroid)),
        "unrelated_edge_crossings": validation["unrelated_edge_crossings"],
        "coincident_vertex_pairs": validation["coincident_vertex_pairs"],
        "vertices_on_unrelated_edges": validation["vertices_on_unrelated_edges"],
        "cyclic_order_violations": validation["cyclic_order_violations"],
        "cyclic_order_orientation": validation["cyclic_order_orientation"],
        "abstract_graph_valid": validation["abstract_graph_valid"],
    }, positions, validation


def write_svg(path, positions, edges, title, guides=None):
    all_points = list(positions.values()) + ([] if guides is None else [guides["center"], guides["pole"]])
    left, right = min(p[0] for p in all_points) - 30, max(p[0] for p in all_points) + 30
    top, bottom = min(p[1] for p in all_points) - 30, max(p[1] for p in all_points) + 30
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.6f} {top:.6f} {right-left:.6f} {bottom-top:.6f}">',
             f'<rect x="{left:.6f}" y="{top:.6f}" width="{right-left:.6f}" height="{bottom-top:.6f}" fill="white"/>',
             f'<title>{title}</title>']
    if guides is not None:
        lines.append('<g stroke="#789" stroke-width="0.8" opacity="0.22">')
        for vertex in sorted(positions):
            lines.append(f'<line x1="{guides["pole"][0]:.6f}" y1="{guides["pole"][1]:.6f}" x2="{positions[vertex][0]:.6f}" y2="{positions[vertex][1]:.6f}"/>')
        lines.append('</g>')
    lines.append('<g stroke="#37474f" stroke-width="1.2" fill="none">')
    for first, second in edges:
        lines.append(f'<line x1="{positions[first][0]:.6f}" y1="{positions[first][1]:.6f}" x2="{positions[second][0]:.6f}" y2="{positions[second][1]:.6f}"/>')
    lines.append('</g><g fill="#1565c0">')
    for vertex in sorted(positions):
        lines.append(f'<circle cx="{positions[vertex][0]:.6f}" cy="{positions[vertex][1]:.6f}" r="3"><title>{vertex}</title></circle>')
    lines.append('</g>')
    if guides is not None:
        for label, point, color in (("C", guides["center"], "#d32f2f"), ("P", guides["pole"], "#7b1fa2")):
            lines.append(f'<circle cx="{point[0]:.6f}" cy="{point[1]:.6f}" r="5" fill="{color}"/><text x="{point[0]+7:.6f}" y="{point[1]-7:.6f}" font-family="monospace" font-size="12" fill="{color}">{label}</text>')
    lines.append('</svg>')
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_csv(path, records):
    fields = list(records[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def run(input_path, output_dir):
    started = time.perf_counter()
    output_dir.mkdir(parents=True, exist_ok=True)
    circuit = gr.parse_circuit(input_path)
    graph, embedding, initial, edges = load_abstract_start(circuit, output_dir)
    model = external_model(initial)
    center, farthest, radius, pole, cables, _, hmax, owner = model
    sweep = []
    payloads = {}
    for fraction in SWEEP_FRACTIONS:
        payload = evaluate(graph, embedding, initial, edges, model, fraction)
        sweep.append(payload[0])
        payloads[fraction] = payload
    valid = [record for record in sweep if record["abstract_graph_valid"]]
    if not valid:
        raise RuntimeError("even h=0 failed abstract graph validation")
    selected_record = max(valid, key=lambda record: record["h"])
    selected_fraction = selected_record["fraction"]
    selected_record, selected_positions, selected_validation = payloads[selected_fraction]
    baseline = sweep[0]
    ratios = {key: gr.rounded(selected_record[key] / baseline[key]) for key in ("width", "height", "area", "total_edge_length")}

    svg_paths = []
    diagnostic_samples = {}
    for fraction in SVG_FRACTIONS:
        record, positions, validation = payloads.get(fraction, evaluate(graph, embedding, initial, edges, model, fraction))
        name = f'planar_external_tent_abstract_{fraction:.2f}'.replace('.', '_') + '.svg'
        path = output_dir / name
        write_svg(path, positions, edges, f'Abstract external tent h/Hmax={fraction:.2f}')
        svg_paths.append(path)
        diagnostic_samples[f'{fraction:.2f}'] = record

    compact_path = output_dir / "planar_external_tent_compact.svg"
    diagnostic_path = output_dir / "planar_external_tent_diagnostic.svg"
    write_svg(compact_path, selected_positions, edges, "Selected abstract external-tent graph")
    write_svg(diagnostic_path, selected_positions, edges, "Selected abstract external-tent graph with pole rays",
              {"center": center, "pole": pole})

    coordinate_path = output_dir / "planar_external_tent_coordinates.json"
    gr.write_json(coordinate_path, {
        "model": "pure abstract geometric graph: zero-size vertices and straight center-to-center edges",
        "source": "stored Experiment 4K point coordinates",
        "center": point_doc(center), "farthest_vertex": farthest,
        "farthest_coordinates": point_doc(initial[farthest]), "radius": gr.rounded(radius),
        "pole": point_doc(pole), "minimum_cable_length": gr.rounded(hmax),
        "minimum_cable_vertex": owner, "selected_fraction": selected_fraction,
        "selected_height": selected_record["h"],
        "vertices": {v: point_doc(selected_positions[v]) for v in sorted(selected_positions)},
        "edges": [list(edge) for edge in edges],
    })
    metrics_path = output_dir / "planar_external_tent_metrics.json"
    gr.write_json(metrics_path, {
        "model": "abstract graph only", "physical_geometry_participated": False,
        "center": point_doc(center), "farthest_vertex": farthest, "radius": gr.rounded(radius),
        "pole": point_doc(pole), "minimum_cable_length": gr.rounded(hmax), "minimum_cable_vertex": owner,
        "sweep": sweep, "svg_diagnostic_samples": diagnostic_samples,
        "selected": selected_record, "selected_ratios_to_h0": ratios,
    })
    validation_path = output_dir / "planar_external_tent_validation.json"
    gr.write_json(validation_path, {
        "validation_model": "abstract graph only", "physical_validation_participated": False,
        "vertex_count": len(initial), "edge_count": len(edges),
        "selected_fraction": selected_fraction, "selected_validation": selected_validation,
    })
    csv_path = output_dir / "planar_external_tent_sweep.csv"
    write_csv(csv_path, sweep)
    paths = [coordinate_path, metrics_path, validation_path, compact_path, diagnostic_path, csv_path] + svg_paths
    return {
        "runtime_seconds": gr.rounded(time.perf_counter() - started),
        "selected": selected_record, "ratios": ratios, "sweep": sweep,
        "hashes": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.input, args.output_dir), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
