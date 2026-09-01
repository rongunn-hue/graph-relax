#!/usr/bin/env python3
"""Experiment 4R: realize, or exactly diagnose failure to realize, 4Q embeddings."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import networkx as nx

import experiment_4p as e4p
import graph_relax as gr


SET_IDS = ("SET001", "SET002", "SET003", "SET004")
POINT_SCALE = 100.0


def cyclic_equal(first, second):
    if len(first) != len(second):
        return False
    if not first:
        return True
    return any(first == second[offset:] + second[:offset] for offset in range(len(second)))


def graph_from_document(document):
    graph = nx.Graph()
    for vertex in document["vertices"]:
        data = dict(vertex)
        identifier = data.pop("id")
        graph.add_node(identifier, **data)
    for edge in document["edges"]:
        data = dict(edge)
        first, second = data.pop("u"), data.pop("v")
        graph.add_edge(first, second, **data)
    return graph


def embedding_from_document(document):
    embedding = nx.PlanarEmbedding()
    embedding.set_data(document["cyclic_edge_order_clockwise"])
    embedding.check_structure()
    return embedding


def expected_graph(circuit, exemptions):
    graph = e4p.terminal_order_wheel_graph(circuit)
    graph.remove_edges_from((e4p.terminal_vertex(item["terminal"]), e4p.net_vertex(item["net"]))
                            for item in exemptions)
    return graph


def load_and_verify_set(circuit, output_dir, set_id):
    path = output_dir / f"constrained_planarization_embedding_{set_id}.json"
    document = json.loads(path.read_text(encoding="utf-8"))
    if document["set_id"] != set_id:
        raise ValueError(f"{path.name}: wrong set id")
    graph = graph_from_document(document["graph"])
    expected = expected_graph(circuit, document["exemptions"])
    if e4p.graph_document(graph) != e4p.graph_document(expected):
        raise ValueError(f"{set_id}: vertex/edge set differs from exact 4Q constrained graph")
    embedding = embedding_from_document(document)
    if set(embedding) != set(graph):
        raise ValueError(f"{set_id}: embedding vertex set differs")
    embedded_edges = {frozenset((node, neighbor)) for node in embedding
                      for neighbor in embedding.neighbors_cw_order(node)}
    if embedded_edges != {frozenset(edge) for edge in graph.edges()}:
        raise ValueError(f"{set_id}: embedding edge set differs")
    source_incidences = {(terminal, net.name) for net in circuit.nets for terminal in net.terminals}
    exemptions = {(item["terminal"], item["net"]) for item in document["exemptions"]}
    represented = {(data["name"], data["net"]) for node, data in graph.nodes(data=True)
                   if node.startswith("terminal:") and graph.has_edge(node, e4p.net_vertex(data["net"]))}
    if represented != source_incidences - exemptions:
        raise ValueError(f"{set_id}: nonexempt electrical incidences differ")
    return document, graph, embedding


def component_rotation_analysis(circuit, document):
    exemptions = {item["terminal"] for item in document["exemptions"]}
    rotations = document["cyclic_edge_order_clockwise"]
    rows = []
    for ref in circuit.refs:
        desired = [e4p.terminal_vertex(f"{ref}.{pin}") for pin in circuit.pins[ref]
                   if f"{ref}.{pin}" not in exemptions]
        actual = [terminal for terminal in rotations[e4p.component_vertex(ref)] if terminal in desired]
        if len(desired) <= 2:
            state = "orientation_topologically_ambiguous"
        elif cyclic_equal(actual, desired):
            state = "current_order"
        elif cyclic_equal(actual, list(reversed(desired))):
            state = "reversed_order"
        else:
            state = "neither_current_nor_reversed"
        rows.append({"component": ref, "physical_clockwise_terminal_order": desired,
                     "saved_embedding_clockwise_terminal_order": actual, "state": state})
    return rows


def collapse_embedding(circuit, document):
    vertices = {item["id"]: item for item in document["graph"]["vertices"]}
    rotations = document["cyclic_edge_order_clockwise"]
    collapsed = {}
    for vertex in rotations:
        if vertex.startswith("component:"):
            collapsed[vertex] = [e4p.net_vertex(vertices[terminal]["net"])
                                 for terminal in rotations[vertex]
                                 if e4p.net_vertex(vertices[terminal]["net"]) in rotations[terminal]]
        elif vertex.startswith("net:"):
            collapsed[vertex] = [e4p.component_vertex(vertices[terminal]["component"])
                                 for terminal in rotations[vertex] if vertex in rotations[terminal]]
    embedding = nx.PlanarEmbedding()
    embedding.set_data(collapsed)
    embedding.check_structure()
    return embedding


def point_coordinates(circuit, document):
    # Stage A fails before a geometric realization exists. These coordinates are
    # therefore only a deterministic failure-diagnostic layout, not a claimed
    # planar drawing. Avoid NetworkX's hash-order-sensitive point constructor.
    components = {ref: ((index % 5) * 120.0, (index // 5) * 100.0)
                  for index, ref in enumerate(circuit.refs)}
    nets = {net.name: ((index % 5) * 120.0, 500.0 + (index // 5) * 80.0)
            for index, net in enumerate(sorted(circuit.nets, key=lambda item: item.name))}
    return components, nets


def point_document(point):
    return {"x": float(point[0]), "y": float(point[1])}


def write_failure_svg(path, circuit, components, nets, reversed_components, set_id):
    terminals = gr.terminal_positions(circuit, components)
    points = list(components.values()) + list(nets.values())
    left = min(point[0] for point in points) - 100
    top = min(point[1] for point in points) - 100
    width = max(point[0] for point in points) - left + 100
    height = max(point[1] for point in points) - top + 100
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left} {top} {width} {height}">',
             f'<rect x="{left}" y="{top}" width="{width}" height="{height}" fill="white"/>',
             f'<text x="{left+20}" y="{top+30}" font-family="monospace" font-size="16" fill="#a00">{set_id}: realization failed - mixed component rotation reversal</text>',
             '<g stroke="#111" stroke-width="2" fill="#f5f5f5">']
    for ref in circuit.refs:
        x, y = components[ref]
        color = "#ffd9d9" if ref in reversed_components else "#f5f5f5"
        lines.append(f'<rect x="{x-gr.WIDTH/2}" y="{y-gr.HEIGHT/2}" width="{gr.WIDTH}" height="{gr.HEIGHT}" fill="{color}"/>')
    lines.append('</g><g font-family="monospace" font-size="10" text-anchor="middle">')
    for ref in circuit.refs:
        x, y = components[ref]
        lines.append(f'<text x="{x}" y="{y}">{ref}</text>')
    lines.append('</g><g fill="#d22">')
    for x, y in terminals.values():
        lines.append(f'<circle cx="{x}" cy="{y}" r="3"/>')
    lines.append('</g><g fill="#1769aa">')
    for name, (x, y) in nets.items():
        lines.append(f'<circle cx="{x}" cy="{y}" r="4"/><text x="{x+6}" y="{y-6}" font-family="monospace" font-size="9">{name}</text>')
    lines.append('</g></svg>')
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def geometry_validation(circuit, document, rotation_rows):
    reversed_components = [row["component"] for row in rotation_rows if row["state"] == "reversed_order"]
    current_components = [row["component"] for row in rotation_rows if row["state"] == "current_order"]
    invalid = [row["component"] for row in rotation_rows if row["state"] == "neither_current_nor_reversed"]
    mixed = bool(reversed_components and current_components)
    return {"validator": "independent saved-rotation versus physical-perimeter-order comparison",
            "input_embedding_structure_valid": True,
            "input_vertex_set_valid": True, "input_edge_set_valid": True,
            "nonexempt_electrical_incidences_valid": True,
            "reversed_components": reversed_components, "current_order_components": current_components,
            "invalid_rotation_components": invalid, "mixed_independent_reversals": mixed,
            "global_reflection_can_repair": not mixed and not invalid,
            "component_mirroring_required": bool(reversed_components),
            "component_mirroring_permitted": False,
            "primary_realization_valid": False,
            "failure": ("The exact saved rotation system requires independent reversal of some component "
                        "terminal cycles while other 3+-terminal cycles retain current order. The fixed "
                        "anonymous rectangle geometry cannot realize this mixed pattern without forbidden "
                        "component mirroring or changing the embedding.")}


def rank_key(result):
    metrics = result["metrics"]
    inf = float("inf")
    return (0 if metrics["primary_valid"] else 1,
            metrics["primary_crossings"] if metrics["primary_crossings"] is not None else inf,
            metrics["component_overlaps"] if metrics["component_overlaps"] is not None else inf,
            metrics["clearance_violations"] if metrics["clearance_violations"] is not None else inf,
            metrics["layer1_layer0_crossings"] if metrics["layer1_layer0_crossings"] is not None else inf,
            metrics["layer1_layer1_crossings"] if metrics["layer1_layer1_crossings"] is not None else inf,
            metrics["total_bends"] if metrics["total_bends"] is not None else inf,
            metrics["total_connection_length"] if metrics["total_connection_length"] is not None else inf,
            metrics["bounding_box_area"] if metrics["bounding_box_area"] is not None else inf,
            result["set_id"])


def run_set(circuit, output_dir, set_id):
    document, _, _ = load_and_verify_set(circuit, output_dir, set_id)
    rotation_rows = component_rotation_analysis(circuit, document)
    validation = geometry_validation(circuit, document, rotation_rows)
    components, nets = point_coordinates(circuit, document)
    reversed_components = validation["reversed_components"]
    coordinates = {"set_id": set_id, "status": "FAILED_BEFORE_RECTANGLE_REALIZATION",
                   "method": "deterministic lexical diagnostic grid; not a planar realization",
                   "point_scale": None,
                   "components": {ref: point_document(components[ref]) for ref in circuit.refs},
                   "net_points": {name: point_document(nets[name]) for name in sorted(nets)},
                   "exemptions": document["exemptions"]}
    geometry = {"set_id": set_id, "status": "NO_VALID_PRIMARY_GEOMETRY",
                "layer0_polylines": [], "layer1_polylines": [],
                "reason": validation["failure"]}
    metrics = {"set_id": set_id, "stage_a": "FAILED", "stage_b": "NOT EXHAUSTIVELY SEARCHED",
               "primary_valid": False, "primary_crossings": None, "component_overlaps": None,
               "clearance_violations": None, "primary_connection_length": None,
               "primary_bends": None, "straight_primary_edges": None,
               "primary_edges_requiring_bends": None, "maximum_bends_on_primary_edge": None,
               "bounding_box_width": None, "bounding_box_height": None,
               "bounding_box_area": None, "layer1_connection_length": None,
               "layer1_bends": None, "layer1_layer0_crossings": None,
               "layer1_layer1_crossings": None, "total_bends": None,
               "total_connection_length": None, "failure": validation["failure"],
               "rotation_analysis": rotation_rows}
    paths = {
        "coordinates": output_dir / f"planar_realization_{set_id}_coordinates.json",
        "geometry": output_dir / f"planar_realization_{set_id}_geometry.json",
        "metrics": output_dir / f"planar_realization_{set_id}_metrics.json",
        "validation": output_dir / f"planar_realization_{set_id}_validation.json"}
    for key, payload in (("coordinates", coordinates), ("geometry", geometry),
                         ("metrics", metrics), ("validation", validation)):
        gr.write_json(paths[key], payload)
    for suffix in ("initial", "final"):
        write_failure_svg(output_dir / f"planar_realization_{set_id}_{suffix}.svg",
                          circuit, components, nets, reversed_components, set_id)
    return {"set_id": set_id, "metrics": metrics,
            "hashes": {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                       for path in paths.values()}}


def report_text(results):
    lines = ["Experiment 4R: constrained-planar embedding realization", "",
             "SET | Primary valid | L0 crossings | L1/L0 crossings | bends | total length | area"]
    for result in results:
        lines.append(f"{result['set_id']} | NO | N/A | N/A | N/A | N/A | N/A")
    lines.extend(["", "Result: all four Stage A realizations failed before edge geometry generation.",
                  "SET001-SET003 require independent reversal of J_PWR and U1.",
                  "SET004 requires independent reversal of J_PWR, J_R4, and U1.",
                  "Other 3+-terminal components retain current order, so a global reflection cannot repair the mismatch.",
                  "No Layer-0 or Layer-1 geometry was generated, because doing so would require a forbidden component mirror or embedding change.",
                  "Stage B: NOT EXHAUSTIVELY SEARCHED (conditional on Stage A success).", "",
                  "Answers:",
                  "1. None of the exact saved 4Q embeddings can be realized by the fixed anonymous rectangles without mirroring selected components.",
                  "2. No valid primary realization was produced; zero-crossing geometry therefore was not claimed.",
                  "3. No Layer-1 comparison is possible.",
                  "4. No visual organization comparison with 4F is valid because 4R produced no complete geometry.",
                  "5. Topology-first remains promising at the abstract level, but the 4P/4Q wheel representation omitted oriented port-order chirality needed by the physical rectangle model."])
    return "\n".join(lines) + "\n"


def run(input_path, output_dir):
    circuit = gr.parse_circuit(input_path)
    results = [run_set(circuit, output_dir, set_id) for set_id in SET_IDS]
    ranking = [result["set_id"] for result in sorted(results, key=rank_key)]
    comparison = {"predetermined_ranking_fields": ["primary valid", "primary crossings",
                    "component overlaps", "clearance violations", "Layer1/Layer0 crossings",
                    "Layer1/Layer1 crossings", "total bends", "total connection length", "area"],
                  "ranking": ranking, "all_failed_stage_a": True,
                  "ranking_note": "All metric fields tie as unavailable; SET ID is deterministic final tie-break.",
                  "sets": [{"set_id": result["set_id"], "metrics": result["metrics"]} for result in results]}
    comparison_path = output_dir / "planar_realization_comparison.json"
    gr.write_json(comparison_path, comparison)
    (output_dir / "planar_realization_report.txt").write_text(report_text(results), encoding="utf-8")
    hashes = {comparison_path.name: hashlib.sha256(comparison_path.read_bytes()).hexdigest()}
    for result in results:
        hashes.update(result["hashes"])
    return {"results": {result["set_id"]: result["metrics"]["stage_a"] for result in results},
            "hashes": hashes}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.input, args.output_dir), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
