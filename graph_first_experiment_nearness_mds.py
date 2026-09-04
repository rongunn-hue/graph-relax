#!/usr/bin/env python3
"""One-shot classical MDS placement from graph shortest-path nearness."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import networkx as nx
import numpy as np

import graph_first_experiment_2 as experiment_2
import placement_geometry as geometry


TOLERANCE = geometry.GEOMETRY_TOLERANCE


def load_graph(path: Path) -> tuple[nx.Graph, dict]:
    document = json.loads(path.read_text(encoding="utf-8"))
    graph = nx.Graph()
    for vertex in document["vertices"]:
        graph.add_node(vertex["id"], type=vertex["type"], name=vertex["name"])
    for edge in document["edges"]:
        graph.add_edge(edge["source"], edge["target"])
    if graph.number_of_nodes() != len(document["vertices"]):
        raise ValueError("vertex reconstruction mismatch")
    if graph.number_of_edges() != len(document["edges"]):
        raise ValueError("incidence reconstruction mismatch")
    if not nx.is_connected(graph):
        raise ValueError("classical graph-distance MDS requires a connected graph")
    return graph, document


def distance_matrix(graph: nx.Graph, nodes: list[str]) -> np.ndarray:
    lengths = dict(nx.all_pairs_shortest_path_length(graph))
    return np.array([[float(lengths[first][second]) for second in nodes]
                     for first in nodes], dtype=float)


def deterministic_degenerate_basis(values: np.ndarray, vectors: np.ndarray):
    """Canonicalize eigenspaces using projections of lexical standard axes."""
    result = vectors.copy()
    groups = []
    start = 0
    while start < len(values):
        end = start + 1
        scale = max(1.0, abs(values[start]))
        while end < len(values) and abs(values[end] - values[start]) <= 1e-11 * scale:
            end += 1
        group = list(range(start, end))
        groups.append(group)
        if len(group) > 1:
            subspace = vectors[:, group]
            projector = subspace @ subspace.T
            basis = []
            for axis in range(len(values)):
                candidate = projector[:, axis].copy()
                for prior in basis:
                    candidate -= np.dot(candidate, prior) * prior
                norm = np.linalg.norm(candidate)
                if norm > 1e-10:
                    basis.append(candidate / norm)
                if len(basis) == len(group):
                    break
            if len(basis) != len(group):
                raise RuntimeError("could not canonicalize degenerate eigenspace")
            result[:, group] = np.column_stack(basis)
        start = end
    return result, groups


def classical_mds(distances: np.ndarray, nodes: list[str]):
    count = len(nodes)
    centering = np.eye(count) - np.ones((count, count)) / count
    gram = -0.5 * centering @ (distances ** 2) @ centering
    values, vectors = np.linalg.eigh(gram)
    order = np.argsort(values)[::-1]
    values = values[order]
    vectors = vectors[:, order]
    vectors, groups = deterministic_degenerate_basis(values, vectors)
    positive = [index for index, value in enumerate(values) if value > TOLERANCE]
    if len(positive) < 2:
        raise ValueError("distance matrix has fewer than two positive MDS eigenvalues")
    chosen = positive[:2]
    coordinates = vectors[:, chosen] * np.sqrt(values[chosen])

    # Each nondegenerate axis still has sign ambiguity. Use the lexically first
    # vertex among equal greatest magnitudes and make that coordinate positive.
    sign_anchors = []
    for axis in range(2):
        magnitudes = np.abs(coordinates[:, axis])
        maximum = float(np.max(magnitudes))
        candidates = [index for index, value in enumerate(magnitudes)
                      if math.isclose(float(value), maximum, rel_tol=0.0, abs_tol=1e-12)]
        anchor = min(candidates, key=lambda index: nodes[index])
        if coordinates[anchor, axis] < 0:
            coordinates[:, axis] *= -1
        sign_anchors.append(nodes[anchor])
    return coordinates, values, gram, groups, sign_anchors


def choose_center(graph: nx.Graph) -> tuple[str, dict[str, float]]:
    closeness = nx.closeness_centrality(graph)
    degree = dict(graph.degree())
    center = min(graph, key=lambda node: (-closeness[node], -degree[node], node))
    return center, closeness


def center_coordinates(coordinates: np.ndarray, nodes: list[str], center: str) -> np.ndarray:
    centered = coordinates - coordinates[nodes.index(center)]
    centered[nodes.index(center)] = (0.0, 0.0)
    return centered


cross = geometry.cross
point_on_segment_interior = geometry.point_on_segment_interior


def validate(graph: nx.Graph, positions: dict[str, tuple[float, float]]) -> dict:
    diagnostics = geometry.graph_geometry_diagnostics(positions, graph.edges)
    return {
        "minimum_node_distance": diagnostics["minimum_node_distance"],
        "proper_unrelated_edge_crossing_count": diagnostics["proper_unrelated_edge_crossing_count"],
        "proper_unrelated_edge_crossings": diagnostics["proper_unrelated_edge_crossings"],
        "coincident_distinct_vertex_pair_count": diagnostics["coincident_vertex_pair_count"],
        "coincident_distinct_vertices": diagnostics["coincident_vertex_pairs"],
        "vertex_on_unrelated_edge_interior_count": diagnostics["vertex_on_unrelated_edge_interior_count"],
        "vertices_on_unrelated_edge_interiors": diagnostics["vertices_on_unrelated_edge_interiors"],
        "electrical_incidences_preserved": graph.number_of_edges(),
    }


def average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda index: (values[index], index))
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and math.isclose(
                values[order[end]], values[order[start]], rel_tol=0.0, abs_tol=1e-12):
            end += 1
        rank = (start + 1 + end) / 2.0
        for offset in range(start, end):
            ranks[order[offset]] = rank
        start = end
    return ranks


def pearson(first: list[float], second: list[float]) -> float:
    a, b = np.asarray(first), np.asarray(second)
    return float(np.corrcoef(a, b)[0, 1])


def measurements(graph: nx.Graph, nodes: list[str], distances: np.ndarray,
                 positions: dict[str, tuple[float, float]], center: str,
                 gram: np.ndarray, eigenvalues: np.ndarray) -> dict:
    xs = [positions[node][0] for node in nodes]
    ys = [positions[node][1] for node in nodes]
    graph_pairs, euclidean_pairs = [], []
    minimum_distance = math.inf
    for first in range(len(nodes)):
        for second in range(first + 1, len(nodes)):
            graph_pairs.append(float(distances[first, second]))
            geometric = math.dist(positions[nodes[first]], positions[nodes[second]])
            euclidean_pairs.append(geometric)
            minimum_distance = min(minimum_distance, geometric)
    coordinate_matrix = np.array([positions[node] for node in nodes])
    reconstructed_gram = coordinate_matrix @ coordinate_matrix.T
    squared_error = sum((a - b) ** 2 for a, b in zip(graph_pairs, euclidean_pairs))
    return {
        "width": max(xs) - min(xs),
        "height": max(ys) - min(ys),
        "aspect_ratio_width_over_height": ((max(xs) - min(xs)) / (max(ys) - min(ys))),
        "total_straight_edge_length": sum(math.dist(positions[a], positions[b])
                                          for a, b in graph.edges),
        "minimum_vertex_to_vertex_distance": minimum_distance,
        "maximum_vertex_radius_from_center": max(math.dist((0.0, 0.0), positions[node])
                                                  for node in nodes),
        "nearness_preservation": {
            "unordered_pair_count": len(graph_pairs),
            "pearson_correlation": pearson(graph_pairs, euclidean_pairs),
            "spearman_rank_correlation": pearson(average_ranks(graph_pairs),
                                                  average_ranks(euclidean_pairs)),
            "kruskal_stress_1": math.sqrt(squared_error /
                                          sum(value * value for value in graph_pairs)),
            "root_mean_squared_distance_error": math.sqrt(squared_error / len(graph_pairs)),
            "classical_mds_strain_frobenius": float(np.linalg.norm(gram - reconstructed_gram)),
            "positive_eigenvalue_count": sum(float(value) > TOLERANCE for value in eigenvalues),
            "leading_eigenvalues": [float(value) for value in eigenvalues[:6]],
        },
    }


def svg(graph: nx.Graph, positions, center: str, spacing_status: str) -> str:
    scale = 135.0
    margin = 125.0
    scaled = {node: (positions[node][0] * scale, positions[node][1] * scale) for node in graph}
    xs, ys = [p[0] for p in scaled.values()], [p[1] for p in scaled.values()]
    min_x, max_x = min(xs) - margin, max(xs) + margin
    min_y, max_y = min(ys) - margin, max(ys) + margin
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{min_x:.6f} {min_y:.6f} '
        f'{max_x-min_x:.6f} {max_y-min_y:.6f}" width="1400" height="1000">',
        '<rect x="-10000" y="-10000" width="20000" height="20000" fill="white"/>',
        f'<text x="{min_x+14:.6f}" y="{min_y+25:.6f}" font-family="sans-serif" '
        f'font-size="17" fill="{("#b00020" if spacing_status != "SATISFIED" else "#222")}">'
        f'IAMP graph nearness — classical MDS'
        f'{" — REJECTED: COINCIDENT VERTICES" if spacing_status != "SATISFIED" else ""}</text>',
    ]
    for first, second in sorted((min(a, b), max(a, b)) for a, b in graph.edges):
        x1, y1 = scaled[first]
        x2, y2 = scaled[second]
        parts.append(f'<line x1="{x1:.6f}" y1="{y1:.6f}" x2="{x2:.6f}" y2="{y2:.6f}" '
                     'stroke="#777" stroke-width="1.3"/>')
    for node in sorted(graph):
        x, y = scaled[node]
        component = graph.nodes[node]["type"] == "COMPONENT"
        if node == center:
            parts.append(f'<circle cx="{x:.6f}" cy="{y:.6f}" r="9" fill="#ffd600" '
                         'stroke="#111" stroke-width="2.2"/>')
        elif component:
            parts.append(f'<circle cx="{x:.6f}" cy="{y:.6f}" r="5" fill="#1565c0" '
                         'stroke="white" stroke-width="1"/>')
        else:
            parts.append(f'<rect x="{x-4.5:.6f}" y="{y-4.5:.6f}" width="9" height="9" '
                         'fill="#d84315" stroke="white" stroke-width="1"/>')
        prefix = "C:" if component else "N:"
        parts.append(f'<text x="{x+7:.6f}" y="{y-7:.6f}" font-family="sans-serif" '
                     f'font-size="10" fill="#111">{prefix}{graph.nodes[node]["name"]}</text>')
    parts.append('</svg>')
    return "\n".join(parts) + "\n"


def write_json(path: Path, document: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def run(graph_path: Path, output_dir: Path):
    graph, graph_doc = load_graph(graph_path)
    nodes = sorted(graph)
    distances = distance_matrix(graph, nodes)
    raw, eigenvalues, gram, degenerate_groups, anchors = classical_mds(distances, nodes)
    center, closeness = choose_center(graph)
    centered = center_coordinates(raw, nodes, center)
    raw_positions = {node: (float(centered[index, 0]), float(centered[index, 1]))
                     for index, node in enumerate(nodes)}
    raw_validation = validate(graph, raw_positions)
    resolved_positions, coincidence_resolution = geometry.resolve_coincident_vertices(
        raw_positions, geometry.MIN_NODE_DISTANCE)
    positions, spacing = geometry.enforce_minimum_node_spacing(
        resolved_positions, geometry.MIN_NODE_DISTANCE, anchor=(0.0, 0.0))
    validation = validate(graph, positions)
    metric = measurements(graph, nodes, distances, positions, center, gram, eigenvalues)
    center_distances = nx.single_source_shortest_path_length(graph, center)
    vertices = [{
        "id": node,
        "type": graph.nodes[node]["type"],
        "graph_distance_from_center": center_distances[node],
        "closeness_centrality": closeness[node],
        "x": positions[node][0],
        "y": positions[node][1],
        "euclidean_distance_from_center": math.dist((0.0, 0.0), positions[node]),
    } for node in nodes]
    coordinate_doc = {
        "schema": "graph-relax.nearness-classical-mds.v1",
        "authoritative_input": graph_path.as_posix(),
        "old_coordinate_artifacts_read": [],
        "method": "one-shot classical MDS of all-pairs unweighted graph distances",
        "placement_status": ("ACCEPTED" if spacing["status"] == "SATISFIED"
                             else "REJECTED_MINIMUM_SPACING"),
        "placement_pipeline": {
            "coordinate_generation": "classical MDS",
            "coincidence_resolution": coincidence_resolution,
            "minimum_spacing_enforcement": spacing,
            "geometry_validation": validation,
        },
        "center_vertex": center,
        "canonicalization": {
            "translation": "maximum-closeness vertex translated exactly to origin",
            "axis_sign_rule": "greatest absolute coordinate; lexical vertex ID tie-break; positive",
            "axis_sign_anchor_vertices": anchors,
            "degenerate_eigenvalue_groups": degenerate_groups,
        },
        "vertex_order": nodes,
        "distance_matrix": distances.astype(int).tolist(),
        "vertices": vertices,
        "edges": [{"source": min(a, b), "target": max(a, b)}
                  for a, b in sorted((min(a, b), max(a, b)) for a, b in graph.edges)],
    }
    report = {
        "schema": "graph-relax.nearness-classical-mds-report.v1",
        "placement_status": ("ACCEPTED" if spacing["status"] == "SATISFIED"
                             else "REJECTED_MINIMUM_SPACING"),
        "graph_validation": {
            "component_vertices": sum(graph.nodes[node]["type"] == "COMPONENT" for node in graph),
            "net_vertices": sum(graph.nodes[node]["type"] == "NET" for node in graph),
            "total_vertices": graph.number_of_nodes(),
            "electrical_incidence_edges": graph.number_of_edges(),
            "connected": nx.is_connected(graph),
            "device_id_used_in_graph": False,
        },
        "maximum_closeness": {"vertex": center, "value": closeness[center]},
        "eigendecomposition": {
            "eigenvalues_descending": [float(value) for value in eigenvalues],
            "selected_eigenvalues": [float(eigenvalues[0]), float(eigenvalues[1])],
            "degenerate_eigenvalue_groups": degenerate_groups,
            "selected_eigenvalue_degeneracy": len(degenerate_groups[0]) > 1 or len(degenerate_groups[1]) > 1,
        },
        "measurements": metric,
        "minimum_spacing_enforcement": spacing,
        "coincidence_resolution": coincidence_resolution,
        "geometry_validation": validation,
        "raw_geometry_validation": raw_validation,
        "vertices": vertices,
        "interpretation_constraint": (
            "The abstract graph is planar; crossings here are placement artifacts."
        ),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(output_dir / "iamp_nearness_mds.json", coordinate_doc)
    write_json(output_dir / "iamp_nearness_mds_report.json", report)
    (output_dir / "iamp_nearness_mds.svg").write_text(
        svg(graph, positions, center, spacing["status"]), encoding="utf-8")
    (output_dir / "iamp_nearness_mds_raw.svg").write_text(
        svg(graph, raw_positions, center, "RAW_UNSANITIZED"), encoding="utf-8")
    return coordinate_doc, report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--graph", type=Path,
                        default=Path("output/graph_first/iamp_graph.json"))
    parser.add_argument("--output", type=Path, default=Path("output/graph_first"))
    args = parser.parse_args()
    _, report = run(args.graph, args.output)
    print(json.dumps({
        "maximum_closeness": report["maximum_closeness"],
        "measurements": report["measurements"],
        "geometry_validation": report["geometry_validation"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
