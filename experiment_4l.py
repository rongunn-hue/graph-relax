#!/usr/bin/env python3
"""Experiment 4L: one direct barycentric realization of the frozen 4I embedding."""

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


EPSILON = 1e-10


def load_frozen(circuit, output_dir):
    document = json.loads((output_dir / "node_planar_embedding.json").read_text())
    coordinates = json.loads((output_dir / "planar_node_initial_coordinates.json").read_text())
    graph = e4p.incidence_graph(circuit)
    if (graph.number_of_nodes(), graph.number_of_edges()) != (34, 44):
        raise ValueError("4L source graph is not the frozen 34-vertex/44-edge graph")
    if e4p.graph_document(graph) != document["graph"]:
        raise ValueError("4L source graph differs from stored 4I graph")
    embedding = nx.PlanarEmbedding()
    embedding.set_data(document["cyclic_neighbor_order_clockwise"])
    embedding.check_structure()
    verified, _ = e4p.verify_embedding(graph, embedding)
    if not verified:
        raise ValueError("stored 4I embedding failed verification")
    old_points = {e4p.component_vertex(ref):
                  (coordinates["components"][ref]["x"] / 72.0,
                   coordinates["components"][ref]["y"] / 72.0)
                  for ref in circuit.refs}
    old_points.update({e4p.net_vertex(net.name):
                       (coordinates["net_junctions"][net.name]["x"] / 72.0,
                        coordinates["net_junctions"][net.name]["y"] / 72.0)
                       for net in circuit.nets})
    return graph, embedding, document, old_points


def enumerate_faces(embedding):
    marked, faces = set(), []
    for first in sorted(embedding):
        for second in embedding.neighbors_cw_order(first):
            if (first, second) not in marked:
                faces.append(embedding.traverse_face(first, second, marked))
    return faces


def signed_area(face, points):
    return 0.5 * sum(points[face[index]][0] * points[face[(index + 1) % len(face)]][1] -
                     points[face[(index + 1) % len(face)]][0] * points[face[index]][1]
                     for index in range(len(face)))


def existing_outer_face(embedding, stored_points):
    """Recover the boundary face of the exact stored 4I point drawing."""
    faces = enumerate_faces(embedding)
    indexed = [(abs(signed_area(face, stored_points)), index, face)
               for index, face in enumerate(faces)]
    _, index, face = max(indexed, key=lambda row: (row[0], -row[1]))
    return index, list(face), faces


def solve_linear(matrix, vector):
    """Deterministic Gaussian elimination with partial pivoting."""
    size = len(vector)
    rows = [list(map(float, matrix[index])) + [float(vector[index])]
            for index in range(size)]
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: (abs(rows[row][column]), -row))
        if abs(rows[pivot][column]) <= EPSILON:
            raise ValueError("singular barycentric equilibrium system")
        rows[column], rows[pivot] = rows[pivot], rows[column]
        divisor = rows[column][column]
        rows[column] = [value / divisor for value in rows[column]]
        for row in range(size):
            if row == column:
                continue
            factor = rows[row][column]
            if factor:
                rows[row] = [rows[row][offset] - factor * rows[column][offset]
                             for offset in range(size + 1)]
    return [rows[index][-1] for index in range(size)]


def barycentric_coordinates(graph, outer_face):
    count = len(outer_face)
    boundary = {vertex: (math.cos(-math.pi / 2 + 2 * math.pi * index / count),
                         math.sin(-math.pi / 2 + 2 * math.pi * index / count))
                for index, vertex in enumerate(outer_face)}
    interior = sorted(set(graph) - set(outer_face))
    lookup = {vertex: index for index, vertex in enumerate(interior)}
    matrix = [[0.0] * len(interior) for _ in interior]
    x_rhs = [0.0] * len(interior)
    y_rhs = [0.0] * len(interior)
    for vertex in interior:
        row = lookup[vertex]
        neighbors = sorted(graph.neighbors(vertex))
        matrix[row][row] = float(len(neighbors))
        for neighbor in neighbors:
            if neighbor in lookup:
                matrix[row][lookup[neighbor]] -= 1.0
            else:
                x_rhs[row] += boundary[neighbor][0]
                y_rhs[row] += boundary[neighbor][1]
    xs = solve_linear(matrix, x_rhs)
    ys = solve_linear(matrix, y_rhs)
    points = dict(boundary)
    points.update({vertex: (xs[lookup[vertex]], ys[lookup[vertex]]) for vertex in interior})
    return points, interior


def coincident_groups(points):
    groups = []
    remaining = set(points)
    while remaining:
        first = min(remaining)
        group = sorted(vertex for vertex in remaining
                       if math.dist(points[first], points[vertex]) <= EPSILON)
        remaining.difference_update(group)
        if len(group) > 1:
            groups.append({"vertices": group, "coordinate": e4i.point_doc(points[first])})
    return groups


def diagnose(circuit, graph, points):
    groups = coincident_groups(points)
    component_pairs = []
    zero_incidence_edges = []
    for group in groups:
        components = sorted(vertex for vertex in group["vertices"] if vertex.startswith("component:"))
        for index, first in enumerate(components):
            for second in components[index + 1:]:
                component_pairs.append([first, second])
    for first, second in sorted(graph.edges()):
        if math.dist(points[first], points[second]) <= EPSILON:
            zero_incidence_edges.append([first, second])
    return groups, component_pairs, zero_incidence_edges


def normalized_extents(points):
    xs = [point[0] for point in points.values()]
    ys = [point[1] for point in points.values()]
    return {"min_x": gr.rounded(min(xs)), "max_x": gr.rounded(max(xs)),
            "min_y": gr.rounded(min(ys)), "max_y": gr.rounded(max(ys)),
            "width": gr.rounded(max(xs)-min(xs)), "height": gr.rounded(max(ys)-min(ys))}


def write_failure_svg(path, graph, points, groups):
    scale, margin = 300.0, 45.0
    transformed = {vertex: (point[0] * scale, point[1] * scale) for vertex, point in points.items()}
    xs = [point[0] for point in transformed.values()]
    ys = [point[1] for point in transformed.values()]
    left, top = min(xs)-margin, min(ys)-margin
    width, height = max(xs)-min(xs)+2*margin, max(ys)-min(ys)+2*margin
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left:.3f} {top:.3f} {width:.3f} {height:.3f}">',
             f'<rect x="{left:.3f}" y="{top:.3f}" width="{width:.3f}" height="{height:.3f}" fill="white"/>',
             '<g stroke="#999" stroke-width="1" fill="none">']
    for first, second in sorted(graph.edges()):
        a, b = transformed[first], transformed[second]
        lines.append(f'<line x1="{a[0]:.3f}" y1="{a[1]:.3f}" x2="{b[0]:.3f}" y2="{b[1]:.3f}"/>')
    lines.append('</g><g font-family="monospace" font-size="8">')
    for vertex in sorted(points):
        x, y = transformed[vertex]
        color = "#b00020" if any(vertex in group["vertices"] for group in groups) else ("#1769aa" if vertex.startswith("net:") else "#222")
        lines.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="3" fill="{color}"/>')
        lines.append(f'<text x="{x+5:.3f}" y="{y-5:.3f}" fill="{color}">{vertex}</text>')
    lines.append('</g><text x="{:.3f}" y="{:.3f}" font-family="monospace" font-size="10" fill="#b00020">FAILED: coincident barycentric vertices prevent rectangle expansion</text></svg>'.format(left+5, top+14))
    path.write_text("\n".join(lines)+"\n", encoding="utf-8")


def run(input_path, output_dir):
    circuit = gr.parse_circuit(input_path)
    graph, embedding, embedding_document, stored_points = load_frozen(circuit, output_dir)
    face_index, outer_face, faces = existing_outer_face(embedding, stored_points)
    points, interior = barycentric_coordinates(graph, outer_face)
    groups, component_pairs, zero_edges = diagnose(circuit, graph, points)
    # Equal normalized coordinates remain equal under every finite uniform scale.
    impossible = bool(component_pairs or zero_edges)
    if not impossible:
        raise RuntimeError("unexpectedly nondegenerate realization; scale search is not implemented for this outcome")
    coordinate_path = output_dir / "planar_barycentric_coordinates.json"
    coordinate_document = {
        "status": "NO_VALID_PHYSICAL_SCALE", "method": "direct barycentric equilibrium",
        "outer_boundary": "unit regular convex polygon", "outer_face_index": face_index,
        "outer_face_vertex_order": outer_face,
        "normalized_coordinates": {vertex: e4i.point_doc(points[vertex]) for vertex in sorted(points)},
        "normalized_extents": normalized_extents(points), "physical_scale": None,
        "component_dimensions": {"width": gr.WIDTH, "height": gr.HEIGHT}}
    gr.write_json(coordinate_path, coordinate_document)
    validation_path = output_dir / "planar_barycentric_validation.json"
    validation = {
        "status": "FAILED", "graph_vertex_count": 34, "graph_edge_count": 44,
        "frozen_embedding_preserved": True, "face_count": len(faces),
        "selected_outer_face_index": face_index, "selected_outer_face": outer_face,
        "linear_system_solved": True, "interior_vertex_count": len(interior),
        "normalized_coordinate_extents": normalized_extents(points),
        "coincident_vertex_groups": groups,
        "coincident_component_pairs": component_pairs,
        "zero_length_incidence_edges": zero_edges,
        "uniform_scaling_can_resolve": False,
        "failure_reason": "Equal normalized coordinates remain equal under every uniform scale; component rectangles overlap and ray/rectangle attachment is undefined for zero-length incidence rays.",
        "selected_physical_scale": None,
        "electrical_source_incidence_count": graph.number_of_edges(),
        "physical_electrical_incidences_realized": 44-len(zero_edges),
        "missing_or_unrealizable_incidences": len(zero_edges),
        "added_incidences": 0,
        "component_overlaps": len(component_pairs),
        "clearance_violations": None, "different_net_crossings": None,
        "unrelated_component_body_intersections": None,
        "cyclic_order_orientation": "NOT_EVALUABLE_FOR_DEGENERATE_REALIZATION",
        "fully_valid": False}
    gr.write_json(validation_path, validation)
    metrics_path = output_dir / "planar_barycentric_metrics.json"
    metrics = {
        "status": "NO_VALID_PHYSICAL_REALIZATION", "selected_physical_scale": None,
        "drawing_width": None, "drawing_height": None, "drawing_area": None,
        "total_connection_length": None, "crowding": None, "objective": None,
        "crossings": None, "component_overlaps": len(component_pairs),
        "clearance_violations": None, "unrelated_component_body_intersections": None,
        "comparison_4k": {"width": 3166, "height": 1112, "area": 3520592,
                          "total_connection_length": 23095.744883,
                          "crossings": 0, "component_overlaps": 0,
                          "clearance_violations": 0, "unrelated_component_body_intersections": 0},
        "comparison_possible": False}
    gr.write_json(metrics_path, metrics)
    svg_path = output_dir / "planar_barycentric.svg"
    write_failure_svg(svg_path, graph, points, groups)
    paths = (coordinate_path, metrics_path, validation_path, svg_path)
    return {"status": "FAILED", "outer_face": outer_face,
            "normalized_extents": normalized_extents(points),
            "coincident_component_pairs": component_pairs,
            "zero_length_incidence_edges": zero_edges,
            "hashes": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.input, args.output_dir), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
