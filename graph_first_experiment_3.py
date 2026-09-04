#!/usr/bin/env python3
"""Graph-first Experiment 3: first topology-driven radial placement."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import networkx as nx
from networkx.algorithms.planar_drawing import (
    combinatorial_embedding_to_pos,
    triangulate_embedding,
)

import graph_first_experiment_2 as experiment_2
import placement_geometry as geometry


TOLERANCE = geometry.GEOMETRY_TOLERANCE
AUTHORITATIVE_FILENAMES = (
    "iamp_graph.json",
    "iamp_graph_analysis.json",
    "iamp_topology_blueprint.json",
)


def cyclic_equivalent(actual: list[str], expected: list[str]) -> bool:
    if len(actual) != len(expected):
        return False
    return any(actual == expected[index:] + expected[:index] for index in range(len(expected)))


def cyclic_relation(actual: list[str], expected: list[str]) -> str:
    if cyclic_equivalent(actual, expected):
        return "DIRECT"
    reversed_expected = list(reversed(expected))
    if cyclic_equivalent(actual, reversed_expected):
        return "REFLECTED"
    return "MISMATCH"


def load_inputs(graph_path: Path, analysis_path: Path, blueprint_path: Path):
    graph, graph_doc, analysis_doc = experiment_2.load_authoritative_graph(
        graph_path, analysis_path
    )
    blueprint = json.loads(blueprint_path.read_text(encoding="utf-8"))
    if blueprint["graph_counts"] != {
        "vertices": graph.number_of_nodes(), "edges": graph.number_of_edges()
    }:
        raise ValueError("blueprint graph count mismatch")
    stored_authority = Path(blueprint["authoritative_input"])
    if not stored_authority.is_absolute():
        stored_authority = Path.cwd() / stored_authority
    if stored_authority.resolve() != graph_path.resolve():
        raise ValueError("blueprint authority mismatch")
    return graph, graph_doc, analysis_doc, blueprint


def same_cyclic_face(first: list[str], second: list[str]) -> bool:
    return cyclic_equivalent(first, second) or cyclic_equivalent(first, list(reversed(second)))


def fresh_reference(graph: nx.Graph, analysis: dict, blueprint: dict):
    embedding = experiment_2.restore_embedding(graph, analysis)
    # The drawing implementation uses sets internally. Stable integer labels
    # prevent Python's randomized string hashes from changing traversal order.
    ordered_nodes = sorted(graph)
    to_integer = {node: index for index, node in enumerate(ordered_nodes)}
    from_integer = {index: node for node, index in to_integer.items()}
    integer_embedding = nx.PlanarEmbedding()
    integer_embedding.set_data({
        to_integer[node]: [to_integer[neighbor] for neighbor in embedding.neighbors_cw_order(node)]
        for node in ordered_nodes
    })
    integer_embedding.check_structure()
    _, integer_outer_face = triangulate_embedding(integer_embedding, fully_triangulate=False)
    algorithm_outer_face = [from_integer[node] for node in integer_outer_face]
    selected_id = blueprint["planar_embedding"]["best_outer_face_candidates"][0]
    selected_face = next(face["boundary_walk"] for face in blueprint["planar_embedding"]["faces"]
                         if face["id"] == selected_id)
    integer_positions = combinatorial_embedding_to_pos(integer_embedding, fully_triangulate=False)
    positions = {from_integer[node]: point for node, point in integer_positions.items()}
    center = blueprint["structural_center_selection"]["chosen_vertex"]
    cx, cy = positions[center]
    translated = {
        node: (float(positions[node][0] - cx), float(positions[node][1] - cy))
        for node in sorted(graph)
    }
    return embedding, translated, {
        "algorithm": "NetworkX combinatorial_embedding_to_pos; Chrobak-Payne grid drawing",
        "deterministic_relabeling": {
            "method": "lexically sorted graph vertex IDs mapped to consecutive integers",
            "mapping": {node: to_integer[node] for node in ordered_nodes},
        },
        "outer_face_api_support": False,
        "outer_face_note": (
            "The routine cannot accept a requested face ID. With fully_triangulate=False it "
            "mechanically leaves a longest encountered face untriangulated."
        ),
        "selected_topological_face_id": selected_id,
        "selected_topological_face_boundary": selected_face,
        "algorithm_outer_face_boundary": list(algorithm_outer_face),
        "algorithm_outer_face_matches_selected": same_cyclic_face(
            list(algorithm_outer_face), selected_face
        ),
    }


def radial_positions(reference: dict[str, tuple[float, float]], blueprint: dict):
    layers = {item["id"]: item["bfs_layer"] for item in blueprint["vertices"]}
    center = blueprint["structural_center_selection"]["chosen_vertex"]
    angles = {}
    positions = {}
    for node in sorted(reference):
        if node == center:
            angles[node] = None
            positions[node] = (0.0, 0.0)
            continue
        x, y = reference[node]
        if math.hypot(x, y) <= TOLERANCE:
            raise ValueError(f"reference vertex {node} coincides with structural center")
        angle = math.atan2(y, x)
        radius = layers[node]
        angles[node] = angle
        positions[node] = (radius * math.cos(angle), radius * math.sin(angle))
    return positions, angles


cross = geometry.cross
point_on_segment_interior = geometry.point_on_segment_interior
proper_crossing = geometry.proper_segment_crossing


def geometric_cyclic_orders(graph: nx.Graph, positions, embedding: nx.PlanarEmbedding) -> dict:
    results = {}
    determinate_relations = []
    for node in sorted(graph):
        if graph.degree(node) < 3:
            continue
        x, y = positions[node]
        entries = sorted(
            ((math.atan2(positions[neighbor][1] - y, positions[neighbor][0] - x), neighbor)
             for neighbor in graph.neighbors(node)),
            key=lambda item: (item[0], item[1]),
        )
        ambiguous = any(abs(entries[index][0] - entries[index - 1][0]) <= TOLERANCE
                        for index in range(1, len(entries)))
        order = [neighbor for _, neighbor in entries]
        expected = list(embedding.neighbors_cw_order(node))
        relation = "UNDEFINED_COLLINEAR_RAYS" if ambiguous else cyclic_relation(order, expected)
        if relation in {"DIRECT", "REFLECTED"}:
            determinate_relations.append(relation)
        results[node] = {
            "geometric_counterclockwise_order": order,
            "embedding_stored_order": expected,
            "relation": relation,
        }
    consistent = (not any(item["relation"].startswith("MISMATCH") for item in results.values())
                  and len(set(determinate_relations)) <= 1)
    orientation = (determinate_relations[0] if determinate_relations and consistent
                   else "INCONSISTENT_OR_UNDEFINED")
    return {"vertices": results, "globally_consistent": consistent, "orientation": orientation}


def validate_geometry(graph: nx.Graph, positions, embedding: nx.PlanarEmbedding) -> dict:
    diagnostics = geometry.graph_geometry_diagnostics(positions, graph.edges)
    crossings = diagnostics["proper_unrelated_edge_crossings"]
    coincidences = diagnostics["coincident_vertex_pairs"]
    vertex_on_edge = diagnostics["vertices_on_unrelated_edge_interiors"]
    cyclic = geometric_cyclic_orders(graph, positions, embedding)
    valid = not crossings and not coincidences and not vertex_on_edge
    return {
        "different_edge_crossings": crossings,
        "different_edge_crossing_count": len(crossings),
        "minimum_node_distance": diagnostics["minimum_node_distance"],
        "coincident_distinct_vertices": coincidences,
        "coincident_distinct_vertex_pair_count": len(coincidences),
        "vertices_on_unrelated_edge_interiors": vertex_on_edge,
        "vertex_on_unrelated_edge_interior_count": len(vertex_on_edge),
        "electrical_incidences_expected": graph.number_of_edges(),
        "electrical_incidences_preserved": graph.number_of_edges(),
        "electrical_incidence_validation": True,
        "embedding_order_validation": cyclic,
        "abstract_graph_valid": valid,
    }


def measurements(graph: nx.Graph, positions, blueprint: dict) -> dict:
    xs = [point[0] for point in positions.values()]
    ys = [point[1] for point in positions.values()]
    total_length = sum(math.dist(positions[a], positions[b]) for a, b in graph.edges)
    by_layer = {}
    for layer, members in sorted(blueprint["bfs_layers"].items(), key=lambda item: int(item[0])):
        radii = [math.hypot(*positions[node]) for node in members]
        angles = sorted(math.atan2(positions[node][1], positions[node][0]) % (2 * math.pi)
                        for node in members if int(layer) > 0)
        if len(angles) >= 2:
            gaps = [angles[(index + 1) % len(angles)] - angles[index]
                    for index in range(len(angles))]
            gaps[-1] += 2 * math.pi
            minimum_angle = min(gaps)
        else:
            minimum_angle = None
        by_layer[layer] = {
            "vertex_count": len(members),
            "mean_radius": sum(radii) / len(radii),
            "minimum_radius": min(radii),
            "maximum_radius": max(radii),
            "minimum_angular_separation_radians": minimum_angle,
            "minimum_angular_separation_degrees": (
                math.degrees(minimum_angle) if minimum_angle is not None else None
            ),
        }
    return {
        "width": max(xs) - min(xs),
        "height": max(ys) - min(ys),
        "total_straight_edge_length": total_length,
        "radial_behavior_by_bfs_layer": by_layer,
    }


def svg_document(graph: nx.Graph, positions, title: str, rejected: bool) -> str:
    scale = 115.0
    margin = 90.0
    xs = [point[0] * scale for point in positions.values()]
    ys = [point[1] * scale for point in positions.values()]
    min_x, max_x = min(xs) - margin, max(xs) + margin
    min_y, max_y = min(ys) - margin, max(ys) + margin
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{min_x:.6f} {min_y:.6f} '
        f'{max_x-min_x:.6f} {max_y-min_y:.6f}">',
        '<rect x="-10000" y="-10000" width="20000" height="20000" fill="white"/>',
        f'<text x="{min_x+12:.6f}" y="{min_y+24:.6f}" font-family="sans-serif" '
        f'font-size="16" fill="{("#b00020" if rejected else "#222")}">'
        f'{title}{" — REJECTED" if rejected else ""}</text>',
    ]
    for first, second in sorted((min(a, b), max(a, b)) for a, b in graph.edges):
        x1, y1 = positions[first][0] * scale, positions[first][1] * scale
        x2, y2 = positions[second][0] * scale, positions[second][1] * scale
        parts.append(f'<line x1="{x1:.6f}" y1="{y1:.6f}" x2="{x2:.6f}" y2="{y2:.6f}" '
                     'stroke="#666" stroke-width="1.4"/>')
    for node in sorted(graph):
        x, y = positions[node][0] * scale, positions[node][1] * scale
        component = graph.nodes[node]["type"] == "COMPONENT"
        color = "#1565c0" if component else "#d84315"
        shape = (f'<circle cx="{x:.6f}" cy="{y:.6f}" r="5"' if component else
                 f'<rect x="{x-4.5:.6f}" y="{y-4.5:.6f}" width="9" height="9"')
        parts.append(shape + f' fill="{color}" stroke="white" stroke-width="1"/>')
        parts.append(f'<text x="{x+7:.6f}" y="{y-7:.6f}" font-family="sans-serif" '
                     f'font-size="10" fill="#111">{graph.nodes[node]["name"]}</text>')
    parts.append('</svg>')
    return "\n".join(parts) + "\n"


def write_json(path: Path, document: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def run(graph_path: Path, analysis_path: Path, blueprint_path: Path, output_dir: Path):
    graph, graph_doc, analysis, blueprint = load_inputs(graph_path, analysis_path, blueprint_path)
    embedding, raw_reference, reference_info = fresh_reference(graph, analysis, blueprint)
    reference, reference_spacing = geometry.enforce_minimum_node_spacing(
        raw_reference, geometry.MIN_NODE_DISTANCE, anchor=(0.0, 0.0))
    if reference_spacing["status"] != "SATISFIED":
        raise RuntimeError("fresh planar reference failed the permanent node-spacing invariant")
    reference_validation = validate_geometry(graph, reference, embedding)
    if not reference_validation["abstract_graph_valid"]:
        raise RuntimeError("fresh planar reference failed geometric planarity validation")
    raw_radial, angles = radial_positions(reference, blueprint)
    radial, radial_spacing = geometry.enforce_minimum_node_spacing(
        raw_radial, geometry.MIN_NODE_DISTANCE, anchor=(0.0, 0.0))
    validation = validate_geometry(graph, radial, embedding)
    metrics = measurements(graph, radial, blueprint)
    center = blueprint["structural_center_selection"]["chosen_vertex"]
    center_expected = list(embedding.neighbors_cw_order(center))
    center_geometric = validation["embedding_order_validation"]["vertices"][center]
    radius_checks = {
        node: {
            "bfs_layer": next(item["bfs_layer"] for item in blueprint["vertices"] if item["id"] == node),
            "actual_radius": math.hypot(*radial[node]),
        }
        for node in sorted(graph)
    }
    radius_valid = all(math.isclose(item["actual_radius"], item["bfs_layer"], abs_tol=TOLERANCE)
                       for item in radius_checks.values())
    coordinates_doc = {
        "schema": "graph-relax.topology-radial-initial.v1",
        "authoritative_inputs": [graph_path.as_posix(), analysis_path.as_posix(), blueprint_path.as_posix()],
        "old_coordinate_artifacts_read": [],
        "construction": "embedding-derived angle; BFS-layer radius",
        "structural_center": center,
        "vertices": [{
            "id": node,
            "type": graph.nodes[node]["type"],
            "bfs_layer": radius_checks[node]["bfs_layer"],
            "reference_angle_radians": angles[node],
            "x": radial[node][0], "y": radial[node][1],
        } for node in sorted(graph)],
        "edges": [{"source": min(a, b), "target": max(a, b)}
                  for a, b in sorted((min(a, b), max(a, b)) for a, b in graph.edges)],
    }
    report = {
        "schema": "graph-relax.topology-radial-report.v1",
        "candidate_status": "ACCEPTED" if validation["abstract_graph_valid"] else "REJECTED",
        "graph_counts": {"vertices": graph.number_of_nodes(), "edges": graph.number_of_edges()},
        "reference_realization": reference_info,
        "reference_minimum_spacing_enforcement": reference_spacing,
        "reference_validation": reference_validation,
        "radial_minimum_spacing_enforcement": radial_spacing,
        "radial_validation": validation,
        "radius_equals_bfs_layer": radius_valid,
        "radius_checks": radius_checks,
        "center_neighbor_order": {
            "embedding_order": center_expected,
            "geometric_order": center_geometric["geometric_counterclockwise_order"],
            "relation": center_geometric["relation"],
        },
        "measurements": metrics,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(output_dir / "iamp_topology_radial_initial.json", coordinates_doc)
    write_json(output_dir / "iamp_topology_radial_report.json", report)
    (output_dir / "iamp_planar_reference.svg").write_text(
        svg_document(graph, reference, "Fresh planar reference", False), encoding="utf-8")
    (output_dir / "iamp_topology_radial_initial.svg").write_text(
        svg_document(graph, radial, "Topology radial initial", not validation["abstract_graph_valid"]),
        encoding="utf-8")
    return coordinates_doc, report


def main() -> None:
    parser = argparse.ArgumentParser()
    base = Path("output/graph_first")
    parser.add_argument("--graph", type=Path, default=base / "iamp_graph.json")
    parser.add_argument("--analysis", type=Path, default=base / "iamp_graph_analysis.json")
    parser.add_argument("--blueprint", type=Path, default=base / "iamp_topology_blueprint.json")
    parser.add_argument("--output", type=Path, default=base)
    args = parser.parse_args()
    _, report = run(args.graph, args.analysis, args.blueprint, args.output)
    print(json.dumps({
        "candidate_status": report["candidate_status"],
        "reference_realization": report["reference_realization"],
        "validation": report["radial_validation"],
        "measurements": report["measurements"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
